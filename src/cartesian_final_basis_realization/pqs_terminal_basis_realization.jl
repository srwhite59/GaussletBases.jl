struct CartesianTerminalBasisBlock
    unit_key::Symbol
    support_indices::Vector{Int}
    support_states::Vector{NTuple{3,Int}}
    coefficients::Union{Nothing,Matrix{Float64}}
    column_range::UnitRange{Int}
end
struct CartesianTerminalBasisRealization
    blocks::Vector{CartesianTerminalBasisBlock}
    final_dimension::Int
    max_cross_overlap::Float64
end
const _TERMINAL_WORKSPACE_BYTES = 64 * 1024^2
function _support_cross(left, right, overlaps)
    sx, sy, sz = overlaps
    matrix = Matrix{Float64}(undef, length(left), length(right))
    @inbounds for j in eachindex(right), i in eachindex(left)
        ix, iy, iz = left[i]
        jx, jy, jz = right[j]
        matrix[i, j] = sx[ix, jx] * sy[iy, jy] * sz[iz, jz]
    end
    return matrix
end
function _support_action(left, right, coefficients, overlaps)
    result = zeros(Float64, length(left), size(coefficients, 2))
    maxcols = max(1, _TERMINAL_WORKSPACE_BYTES ÷ (8 * max(length(left), 1)))
    for firstcol in 1:maxcols:length(right)
        cols = firstcol:min(firstcol + maxcols - 1, length(right))
        result .+= _support_cross(left, @view(right[cols]), overlaps) *
                   @view(coefficients[cols, :])
    end
    return result
end
function _block_action(left, right_states, right_coefficients, overlaps)
    action = _support_action(
        left.support_states, right_states, right_coefficients, overlaps)
    return isnothing(left.coefficients) ? action : transpose(left.coefficients) * action
end
function _block_pair_matrix(left, right, overlaps)
    if isnothing(right.coefficients)
        eye = Matrix{Float64}(I, length(right.support_states), length(right.support_states))
        return _block_action(left, right.support_states, eye, overlaps)
    end
    return _block_action(left, right.support_states, right.coefficients, overlaps)
end
function _support_weights(states, bundles)
    wx = _nested_axis_pgdg(bundles, :x).weights
    wy = _nested_axis_pgdg(bundles, :y).weights
    wz = _nested_axis_pgdg(bundles, :z).weights
    return [wx[s[1]] * wy[s[2]] * wz[s[3]] for s in states]
end
function _push_block!(blocks, key, indices, states, coefficients, nextcol)
    count = isnothing(coefficients) ? length(indices) : size(coefficients, 2)
    push!(blocks, CartesianTerminalBasisBlock(
        key, Int.(indices), NTuple{3,Int}.(states), coefficients,
        nextcol:(nextcol + count - 1)))
    return nextcol + count
end
function _append_direct!(blocks, record, bundles, overlaps, nextcol, identity_atol)
    support = record.support_record
    error = norm(_support_cross(support.support_states, support.support_states, overlaps) -
                 Matrix{Float64}(I, support.support_count, support.support_count), Inf)
    error <= identity_atol ||
        throw(ArgumentError("terminal direct sector overlap is not identity"))
    all(weight -> isfinite(weight) && weight > 0,
        _support_weights(support.support_states, bundles)) ||
        throw(ArgumentError("terminal direct sector weights must be finite and positive"))
    return _push_block!(
        blocks, support.unit_key, support.support_indices, support.support_states,
        nothing, nextcol)
end
function _shell_seed(record, contract, bundles)
    support = record.support_record
    source_shape = Tuple(Int(value) for value in support.source_mode_shape)
    dims = _nested_axis_lengths(bundles)
    full_sides =
        _nested_projected_q_shell_full_sides(bundles, support.outer_box, source_shape)
    full_coefficients = _nested_product_coefficients(full_sides..., dims)
    modes = _nested_projected_q_shell_boundary_comx_product_modes(source_shape)
    indices = _nested_box_support_indices(support.outer_box..., dims)
    states = NTuple{3,Int}[_cartesian_unflat_index(index, dims) for index in indices]
    coefficients = Matrix{Float64}(full_coefficients[indices, modes.column_indices])
    rule = get(contract.metadata, :raw_product_source_retained_rule, nothing)
    !isnothing(rule) && Int(rule.retained_count) == size(coefficients, 2) ||
        throw(ArgumentError("terminal PQS retained rule does not match seed columns"))
    return (; indices, states, coefficients)
end
function _subtract_previous(state, block, residual)
    contribution = isnothing(block.coefficients) ? residual : block.coefficients * residual
    row = Dict{Int,Int}()
    indices = Int[]
    states = NTuple{3,Int}[]
    for (index, point) in Iterators.flatten((
        zip(state.indices, state.states),
        zip(block.support_indices, block.support_states),
    ))
        haskey(row, index) && continue
        row[index] = length(indices) + 1
        push!(indices, index)
        push!(states, point)
    end
    coefficients = zeros(Float64, length(indices), size(state.coefficients, 2))
    for (i, index) in enumerate(state.indices)
        coefficients[row[index], :] .+= @view state.coefficients[i, :]
    end
    for (i, index) in enumerate(block.support_indices)
        coefficients[row[index], :] .-= @view contribution[i, :]
    end
    return (; indices, states, coefficients)
end
function _canonicalize!(coefficients, states, bundles, weight_atol)
    weights = vec(transpose(coefficients) * _support_weights(states, bundles))
    for col in eachindex(weights)
        isfinite(weights[col]) && abs(weights[col]) > weight_atol ||
            throw(ArgumentError("terminal PQS final weight is not finite positive gaugeable"))
        weights[col] < 0 && (coefficients[:, col] .*= -1)
    end
    return coefficients
end
function _realize_shell(record, contract, blocks, bundles, overlaps,
    projection_atol, rank_atol, identity_atol, weight_atol)
    state = _shell_seed(record, contract, bundles)
    for pass in 1:2
        changed = false
        for block in blocks
            residual = _block_action(block, state.states, state.coefficients, overlaps)
            norm(residual, Inf) <= projection_atol && continue
            pass == 2 && throw(ArgumentError("terminal PQS projection residual remains above tolerance"))
            state = _subtract_previous(state, block, residual)
            changed = true
        end
        changed || break
    end
    gram = transpose(state.coefficients) *
           _support_action(state.states, state.states, state.coefficients, overlaps)
    overlap = Symmetric((gram + transpose(gram)) ./ 2)
    rank = count(>(rank_atol), eigvals(overlap))
    rank == record.retained_count ||
        throw(ArgumentError("terminal PQS projected Gram rank does not match retained count"))
    coefficients = _canonicalize!(
        state.coefficients * inv(sqrt(overlap)), state.states, bundles, weight_atol)
    error = norm(transpose(coefficients) *
                 _support_action(state.states, state.states, coefficients, overlaps) -
                 Matrix{Float64}(I, size(coefficients, 2), size(coefficients, 2)), Inf)
    error <= identity_atol ||
        throw(ArgumentError("terminal PQS realized shell overlap is not identity"))
    return state.indices, state.states, coefficients
end
function pqs_terminal_basis_realization(
    support_records, retained_records, transform_contracts, bundles;
    identity_atol::Real = 1.0e-8,
    cross_atol::Real = 1.0e-8,
    projection_atol::Real = 1.0e-12,
    rank_atol::Real = 1.0e-10,
    weight_atol::Real = 1.0e-14,
)
    length(support_records) == length(retained_records) ||
        throw(ArgumentError("terminal support and retained records must match"))
    contracts = Dict(contract.unit_key => contract for contract in transform_contracts)
    overlaps = (_nested_axis_pgdg(bundles, :x).overlap,
        _nested_axis_pgdg(bundles, :y).overlap,
        _nested_axis_pgdg(bundles, :z).overlap)
    blocks = CartesianTerminalBasisBlock[]
    nextcol = 1
    for (support, record) in zip(support_records, retained_records)
        support.unit_key == record.support_record.unit_key ||
            throw(ArgumentError("terminal support and retained record order mismatch"))
        if record.transform_kind === :direct_identity_transform_contract
            nextcol = _append_direct!(
                blocks, record, bundles, overlaps, nextcol, identity_atol)
        elseif record.transform_kind === :pqs_source_modes_boundary_selection_shell_realization_contract
            contract = get(contracts, record.transform_contract_unit_key, nothing)
            isnothing(contract) && throw(ArgumentError("missing terminal transform contract"))
            indices, states, coefficients = _realize_shell(
                record, contract, blocks, bundles, overlaps,
                projection_atol, rank_atol, identity_atol, weight_atol)
            nextcol = _push_block!(
                blocks, record.support_record.unit_key, indices, states,
                coefficients, nextcol)
        else
            throw(ArgumentError("unsupported terminal transform kind for PQS basis realization"))
        end
    end
    max_cross = maximum(
        (norm(_block_pair_matrix(blocks[left], blocks[right], overlaps), Inf)
         for right in 2:length(blocks) for left in 1:(right - 1));
        init = 0.0,
    )
    max_cross <= cross_atol ||
        throw(ArgumentError("terminal cross-block overlap exceeds tolerance"))
    return CartesianTerminalBasisRealization(blocks, nextcol - 1, max_cross)
end

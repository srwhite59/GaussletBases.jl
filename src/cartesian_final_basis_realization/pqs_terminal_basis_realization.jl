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
struct CartesianParentResidualFunctionBlock
    source_unit_key::Symbol
    support_indices::Vector{Int}
    support_states::Vector{NTuple{3,Int}}
    coefficients::Matrix{Float64}
    diagnostics::NamedTuple
end
const _TERMINAL_WORKSPACE_BYTES = 64 * 1024^2
function _support_cross!(matrix, left, right, overlaps)
    sx, sy, sz = overlaps
    @inbounds for j in eachindex(right), i in eachindex(left)
        ix, iy, iz = left[i]
        jx, jy, jz = right[j]
        matrix[i, j] = sx[ix, jx] * sy[iy, jy] * sz[iz, jz]
    end
    return matrix
end
function _support_identity_error(states, overlaps)
    sx, sy, sz = overlaps
    error = 0.0
    @inbounds for j in eachindex(states), i in eachindex(states)
        ix, iy, iz = states[i]
        jx, jy, jz = states[j]
        target = i == j ? 1.0 : 0.0
        value = sx[ix, jx] * sy[iy, jy] * sz[iz, jz]
        error = max(error, abs(value - target))
    end
    return error
end
function _matrix_identity_error(matrix)
    error = 0.0
    @inbounds for j in axes(matrix, 2), i in axes(matrix, 1)
        target = i == j ? 1.0 : 0.0
        error = max(error, abs(matrix[i, j] - target))
    end
    return error
end
function _buffer_view!(buffer, nrow, ncol)
    matrix = buffer[]
    if size(matrix, 1) < nrow || size(matrix, 2) < ncol
        matrix = Matrix{Float64}(undef, nrow, ncol)
        buffer[] = matrix
    end
    return @view matrix[1:nrow, 1:ncol]
end
function _support_action(left, right, coefficients, overlaps)
    result = zeros(Float64, length(left), size(coefficients, 2))
    scratch = Ref(Matrix{Float64}(undef, 0, 0))
    return _support_action!(result, left, right, coefficients, overlaps, scratch)
end
function _support_action!(result, left, right, coefficients, overlaps, scratch)
    fill!(result, 0.0)
    maxcols = max(1, _TERMINAL_WORKSPACE_BYTES ÷ (8 * max(length(left), 1)))
    cross_buffer = _buffer_view!(scratch, length(left), min(length(right), maxcols))
    for firstcol in 1:maxcols:length(right)
        cols = firstcol:min(firstcol + maxcols - 1, length(right))
        ncols = length(cols)
        cross = @view cross_buffer[:, 1:ncols]
        _support_cross!(cross, left, @view(right[cols]), overlaps)
        mul!(result, cross, @view(coefficients[cols, :]), 1.0, 1.0)
    end
    return result
end
function _support_weights(states, bundles)
    wx = _nested_axis_pgdg(bundles, :x).weights
    wy = _nested_axis_pgdg(bundles, :y).weights
    wz = _nested_axis_pgdg(bundles, :z).weights
    return [wx[s[1]] * wy[s[2]] * wz[s[3]] for s in states]
end

function _prf_terminal_overlap(basis, states, coefficients, overlaps)
    result = zeros(Float64, basis.final_dimension, size(coefficients, 2))
    for block in basis.blocks
        action = _support_action(
            block.support_states, states, coefficients, overlaps)
        result[block.column_range, :] .= isnothing(block.coefficients) ?
            action : transpose(block.coefficients) * action
    end
    return result
end

function _prf_moment_summary(states, coefficients, pgdg)
    S = Tuple(axis.overlap for axis in pgdg)
    position = Tuple(axis.position for axis in pgdg)
    x2 = Tuple(axis.x2 for axis in pgdg)
    function expectation(axis, factors)
        matrix = zeros(Float64, length(states), length(states))
        _support_cross!(matrix, states, states, ntuple(
            index -> index == axis ? factors[index] : S[index], 3))
        return transpose(coefficients) * matrix * coefficients
    end
    position_matrices = ntuple(axis -> expectation(axis, position), 3)
    x2_matrices = ntuple(axis -> expectation(axis, x2), 3)
    centers = NTuple{3,Float64}[
        ntuple(axis -> position_matrices[axis][column, column], 3)
        for column in axes(coefficients, 2)
    ]
    spreads = NTuple{3,Float64}[
        ntuple(axis -> sqrt(max(0.0,
            x2_matrices[axis][column, column] - centers[column][axis]^2)), 3)
        for column in axes(coefficients, 2)
    ]
    return (; centers, spreads,
        total_spreads = Float64[sqrt(sum(abs2, spread)) for spread in spreads])
end

function build_parent_residual_function_block(
    basis::CartesianTerminalBasisRealization,
    source_block::CartesianTerminalBasisBlock,
    bundles,
    target_parent_indices,
    target_columns;
    support_atol::Real = 1.0e-12,
    orthogonality_atol::Real = 1.0e-10,
    identity_atol::Real = 1.0e-8,
)
    any(block -> block === source_block, basis.blocks) || throw(ArgumentError(
        "parent residual source block must be an exact block of the terminal basis"))
    support_tolerance = Float64(support_atol)
    orthogonality_tolerance = Float64(orthogonality_atol)
    identity_tolerance = Float64(identity_atol)
    all(value -> isfinite(value) && value >= 0.0,
        (support_tolerance, orthogonality_tolerance, identity_tolerance)) ||
        throw(ArgumentError("parent residual tolerances must be finite and nonnegative"))
    indices = Int.(target_parent_indices)
    length(unique(indices)) == length(indices) || throw(ArgumentError(
        "parent residual target support contains duplicate parent rows"))
    targets = Matrix{Float64}(target_columns)
    size(targets, 1) == length(indices) || throw(DimensionMismatch(
        "parent residual target rows must match target parent indices"))
    size(targets, 2) > 0 || throw(ArgumentError(
        "parent residual construction requires at least one selected target"))
    all(isfinite, targets) || throw(ArgumentError(
        "parent residual target columns must be finite"))

    dims = _nested_axis_lengths(bundles)
    all(index -> 1 <= index <= prod(dims), indices) || throw(ArgumentError(
        "parent residual target support lies outside the resolved parent"))
    support_rows = Dict(index => row for (row, index) in pairs(source_block.support_indices))
    local_targets = zeros(Float64, length(source_block.support_indices), size(targets, 2))
    leakage = 0.0
    for (row, index) in pairs(indices)
        local_row = get(support_rows, index, 0)
        if iszero(local_row)
            leakage = max(leakage, maximum(abs, @view targets[row, :]))
        else
            local_targets[local_row, :] .= @view targets[row, :]
        end
    end
    leakage <= support_tolerance || throw(ArgumentError(
        "parent residual target columns leak outside the exact source support"))

    pgdg = Tuple(_nested_axis_pgdg(bundles, axis) for axis in (:x, :y, :z))
    overlaps = Tuple(axis.overlap for axis in pgdg)
    selected_overlap = _prf_terminal_overlap(
        basis, source_block.support_states, local_targets, overlaps)
    source_overlap = @view selected_overlap[source_block.column_range, :]
    cross_contamination = 0.0
    for block in basis.blocks
        block === source_block && continue
        cross_contamination = max(cross_contamination,
            isempty(block.column_range) ? 0.0 : maximum(abs, @view selected_overlap[block.column_range, :]))
    end
    cross_contamination <= orthogonality_tolerance || throw(ArgumentError(
        "support-local targets have material overlap with another terminal block"))

    projected = local_targets - (isnothing(source_block.coefficients) ?
        Matrix(source_overlap) : source_block.coefficients * source_overlap)
    support_metric = zeros(Float64,
        length(source_block.support_states), length(source_block.support_states))
    _support_cross!(support_metric, source_block.support_states,
        source_block.support_states, overlaps)
    raw_metric = transpose(projected) * support_metric * projected
    metric = Symmetric((raw_metric + transpose(raw_metric)) ./ 2)
    eigenvalues = eigvals(metric)
    all(isfinite, eigenvalues) || throw(ArgumentError(
        "parent residual selected-target metric is not finite"))
    metric_scale = max(maximum(abs, eigenvalues), 1.0)
    metric_tolerance = max(1.0e-12, eps(Float64) * metric_scale * size(metric, 1))
    minimum(eigenvalues) >= -metric_tolerance || throw(ArgumentError(
        "parent residual selected-target metric is materially indefinite"))
    minimum(eigenvalues) > metric_tolerance || throw(ArgumentError(
        "parent residual selected targets lose numerical rank after terminal projection"))
    coefficients = projected * inv(sqrt(metric))
    pivots = Int[]
    for column in axes(coefficients, 2)
        pivot = argmax(abs.(@view coefficients[:, column]))
        coefficients[pivot, column] < 0.0 && (coefficients[:, column] .*= -1)
        push!(pivots, source_block.support_indices[pivot])
    end

    terminal_overlap = _prf_terminal_overlap(
        basis, source_block.support_states, coefficients, overlaps)
    prf_metric = transpose(coefficients) * support_metric * coefficients
    terminal_error = norm(terminal_overlap, Inf)
    identity_error = norm(prf_metric - I, Inf)
    terminal_error <= orthogonality_tolerance || throw(ArgumentError(
        "parent residual block is not orthogonal to the complete terminal basis"))
    identity_error <= identity_tolerance || throw(ArgumentError(
        "parent residual block is not parent-metric orthonormal"))
    charge_sums = vec(sum(abs2, coefficients; dims = 1))
    moments = _prf_moment_summary(source_block.support_states, coefficients, pgdg)
    support_bounds = ntuple(axis -> extrema(state[axis] for state in source_block.support_states), 3)
    diagnostics = (;
        support = (; unit_key = source_block.unit_key, support_bounds,
            support_count = length(source_block.support_indices), maximum_leakage = leakage),
        dimensions = (; parent_count = prod(dims), terminal_count = basis.final_dimension,
            target_count = size(targets, 2), prf_count = size(coefficients, 2)),
        projection = (; selected_target_metric_eigenvalues = eigenvalues,
            projection_loss = norm(local_targets - projected),
            cross_block_contamination = cross_contamination,
            metric_tolerance),
        validation = (; terminal_orthogonality_error = terminal_error,
            parent_metric_identity_error = identity_error),
        phase_pivots = pivots,
        parent_charge_sums = charge_sums,
        moments,
    )
    return CartesianParentResidualFunctionBlock(
        source_block.unit_key, copy(source_block.support_indices),
        copy(source_block.support_states), coefficients, diagnostics)
end

function _validate_parent_residual_collection(
    basis::CartesianTerminalBasisRealization,
    bundles,
    prfs::AbstractVector{<:CartesianParentResidualFunctionBlock},
    orthogonality_atol::Real,
    identity_atol::Real,
)
    !isempty(prfs) || throw(ArgumentError(
        "parent residual collection requires at least one block"))
    orthogonality_tolerance = Float64(orthogonality_atol)
    identity_tolerance = Float64(identity_atol)
    all(value -> isfinite(value) && value >= 0.0,
        (orthogonality_tolerance, identity_tolerance)) || throw(ArgumentError(
        "parent residual collection tolerances must be finite and nonnegative"))

    dims = _nested_axis_lengths(bundles)
    overlaps = Tuple(_nested_axis_pgdg(bundles, axis).overlap for axis in (:x, :y, :z))
    ranges = UnitRange{Int}[]
    next_column = 1
    terminal_row_sums = zeros(Float64, basis.final_dimension)
    for prf in prfs
        length(prf.support_indices) == length(prf.support_states) ==
            size(prf.coefficients, 1) || throw(DimensionMismatch(
            "parent residual support and coefficient rows must match"))
        size(prf.coefficients, 2) > 0 || throw(ArgumentError(
            "parent residual block must contain at least one column"))
        all(isfinite, prf.coefficients) || throw(ArgumentError(
            "parent residual coefficients must be finite"))
        matches = findall(block -> block.unit_key === prf.source_unit_key, basis.blocks)
        length(matches) == 1 || throw(ArgumentError(
            "parent residual source block identity is not unique in the terminal basis"))
        source = basis.blocks[only(matches)]
        prf.support_indices == source.support_indices &&
            prf.support_states == source.support_states || throw(ArgumentError(
            "parent residual support no longer matches its source terminal block"))
        for (index, state) in zip(prf.support_indices, prf.support_states)
            _cartesian_unflat_index(index, dims) == state || throw(ArgumentError(
                "parent residual index/state mapping is inconsistent"))
        end
        count = size(prf.coefficients, 2)
        push!(ranges, next_column:(next_column + count - 1))
        next_column += count
        terminal_row_sums .+= vec(sum(abs, _prf_terminal_overlap(
            basis, prf.support_states, prf.coefficients, overlaps); dims = 2))
    end

    terminal_error = maximum(terminal_row_sums)
    gram = zeros(Float64, next_column - 1, next_column - 1)
    for right in eachindex(prfs), left in 1:right
        support_action = _support_action(
            prfs[left].support_states, prfs[right].support_states,
            prfs[right].coefficients, overlaps)
        block = transpose(prfs[left].coefficients) * support_action
        gram[ranges[left], ranges[right]] .= block
        left != right && (gram[ranges[right], ranges[left]] .= transpose(block))
    end
    identity_error = norm(gram - I, Inf)
    terminal_error <= orthogonality_tolerance || throw(ArgumentError(
        "combined parent residual blocks are not terminal-orthogonal"))
    identity_error <= identity_tolerance || throw(ArgumentError(
        "combined parent residual blocks are not parent-metric orthonormal"))
    return ranges, (; terminal_orthogonality_error = terminal_error,
        parent_metric_identity_error = identity_error,
        block_count = length(prfs), prf_count = next_column - 1)
end
function _push_block!(blocks, key, indices, states, coefficients, nextcol)
    count = isnothing(coefficients) ? length(indices) : size(coefficients, 2)
    push!(blocks, CartesianTerminalBasisBlock(
        key, Int.(indices), NTuple{3,Int}.(states), coefficients,
        nextcol:(nextcol + count - 1)))
    return nextcol + count
end
function _validate_block_support!(block, support, seen)
    length(block.support_indices) == length(support.support_indices) &&
        length(block.support_states) == length(support.support_states) ||
        throw(ArgumentError("terminal block support length does not match owned support"))
    local_seen = Set{Int}()
    for i in eachindex(block.support_indices)
        index = block.support_indices[i]
        index == support.support_indices[i] &&
            block.support_states[i] == support.support_states[i] ||
            throw(ArgumentError("terminal block support does not match owned support"))
        index in local_seen &&
            throw(ArgumentError("terminal block support contains duplicate parent rows"))
        push!(local_seen, index)
        index in seen &&
            throw(ArgumentError("terminal block supports are not disjoint"))
        push!(seen, index)
    end
    return nothing
end
function _append_direct!(blocks, record, bundles, overlaps, nextcol, identity_atol)
    support = record.support_record
    error = _support_identity_error(support.support_states, overlaps)
    error <= identity_atol ||
        throw(ArgumentError("terminal direct sector overlap is not identity"))
    all(weight -> isfinite(weight) && weight > 0,
        _support_weights(support.support_states, bundles)) ||
        throw(ArgumentError("terminal direct sector weights must be finite and positive"))
    return _push_block!(
        blocks, support.unit_key, support.support_indices, support.support_states,
        nothing, nextcol)
end
function _validate_shell_seed_contract(support, rule, raw_plan, modes, source_shape)
    isnothing(rule) &&
        throw(ArgumentError("terminal PQS seed requires retained source-mode rule"))
    rule.retained_rule_kind === :boundary_comx_product_mode_selection ||
        throw(ArgumentError("terminal PQS seed requires boundary retained rule"))
    rule.transform_kind === :source_mode_column_selector ||
        throw(ArgumentError("terminal PQS seed requires source-mode column selector"))
    rule.source_mode_dims == source_shape ||
        throw(ArgumentError("terminal PQS seed retained-rule source-mode mismatch"))
    rule.source_mode_ordering === :x_major_y_major_z_fast ||
        throw(ArgumentError("terminal PQS seed requires x-major/y-major/z-fast ordering"))
    rule.retained_column_indices == modes.column_indices ||
        throw(ArgumentError("terminal PQS seed retained columns do not match generated boundary columns"))
    rule.retained_mode_indices == modes.mode_indices ||
        throw(ArgumentError("terminal PQS seed retained modes do not match generated boundary modes"))
    rule.retained_count == length(modes.column_indices) ||
        throw(ArgumentError("terminal PQS seed retained count does not match boundary columns"))
    isnothing(raw_plan) && return nothing
    raw_plan.source_intervals == support.outer_box ||
        throw(ArgumentError("terminal PQS seed raw source intervals do not match support outer box"))
    raw_plan.source_shape == length.(support.outer_box) ||
        throw(ArgumentError("terminal PQS seed raw source shape does not match support outer box"))
    raw_plan.source_mode_dims == source_shape ||
        throw(ArgumentError("terminal PQS seed raw source-mode dimensions mismatch"))
    raw_plan.source_mode_ordering == rule.source_mode_ordering ||
        throw(ArgumentError("terminal PQS seed raw source-mode ordering mismatch"))
    raw_plan.source_mode_count == prod(source_shape) ||
        throw(ArgumentError("terminal PQS seed raw source-mode count mismatch"))
    return nothing
end
function _shell_seed_axis_facts(contract)
    facts = get(contract.metadata, :raw_product_source_axis_transform_facts, nothing)
    isnothing(facts) && return nothing
    facts isa Tuple && length(facts) == 3 ||
        throw(ArgumentError("terminal PQS seed requires exactly three source-axis transform facts"))
    all(fact -> fact isa CRPS.AxisSourceTransformFact &&
        fact.coefficient_status === :not_materialized, facts) && return nothing
    return facts
end
function _validate_shell_seed_axis_fact(fact, axis, interval, source_dim)
    fact isa CRPS.AxisSourceTransformFact ||
        throw(ArgumentError("terminal PQS seed source-axis transform fact type mismatch"))
    fact.axis == axis ||
        throw(ArgumentError("terminal PQS seed source-axis transform fact axis mismatch"))
    fact.coefficient_status === :materialized ||
        throw(ArgumentError("terminal PQS seed source-axis transform fact must be materialized"))
    fact.source_interval == interval ||
        throw(ArgumentError("terminal PQS seed source-axis transform fact interval mismatch"))
    fact.source_mode_dim == source_dim ||
        throw(ArgumentError("terminal PQS seed source-axis transform fact source-mode dimension mismatch"))
    size(fact.coefficient_matrix) == (length(interval), source_dim) ||
        throw(ArgumentError("terminal PQS seed source-axis transform fact coefficient shape mismatch"))
    return Matrix{Float64}(fact.coefficient_matrix)
end
function _shell_seed_full_coefficients_from_axis_facts(facts, support, source_shape, dims)
    matrices = ntuple(
        axis -> _validate_shell_seed_axis_fact(
            facts[axis],
            axis,
            support.outer_box[axis],
            source_shape[axis],
        ),
        3,
    )
    ncols = prod(source_shape)
    row_indices = Int[]
    col_indices = Int[]
    values = Float64[]
    column = 0
    for ixcol in 1:source_shape[1], iycol in 1:source_shape[2], izcol in 1:source_shape[3]
        column += 1
        for (ix_local, ix) in enumerate(support.outer_box[1])
            vx = matrices[1][ix_local, ixcol]
            iszero(vx) && continue
            for (iy_local, iy) in enumerate(support.outer_box[2])
                vy = matrices[2][iy_local, iycol]
                iszero(vy) && continue
                for (iz_local, iz) in enumerate(support.outer_box[3])
                    vz = matrices[3][iz_local, izcol]
                    iszero(vz) && continue
                    push!(row_indices, _cartesian_flat_index(ix, iy, iz, dims))
                    push!(col_indices, column)
                    push!(values, vx * vy * vz)
                end
            end
        end
    end
    return SparseArrays.sparse(row_indices, col_indices, values, prod(dims), ncols)
end
function _shell_seed_full_coefficients(record, contract, bundles, source_shape, dims)
    support = record.support_record
    facts = _shell_seed_axis_facts(contract)
    if !isnothing(facts)
        return _shell_seed_full_coefficients_from_axis_facts(facts, support, source_shape, dims)
    end
    full_sides = _nested_projected_q_shell_full_sides(bundles, support.outer_box, source_shape)
    return _nested_product_coefficients(full_sides..., dims)
end
function _shell_seed(record, contract, bundles)
    support = record.support_record
    source_shape = Tuple(Int(value) for value in support.source_mode_shape)
    dims = _nested_axis_lengths(bundles)
    full_coefficients = _shell_seed_full_coefficients(record, contract, bundles, source_shape, dims)
    modes = _nested_projected_q_shell_boundary_comx_product_modes(source_shape)
    indices = support.support_indices
    states = support.support_states
    coefficients = Matrix{Float64}(full_coefficients[indices, modes.column_indices])
    rule = get(contract.metadata, :raw_product_source_retained_rule, nothing)
    raw_plan = get(contract.metadata, :raw_product_source_plan, nothing)
    _validate_shell_seed_contract(support, rule, raw_plan, modes, source_shape)
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
function _realize_shell(record, contract, bundles, overlaps, identity_atol, weight_atol)
    state = _shell_seed(record, contract, bundles)
    gram = transpose(state.coefficients) *
           _support_action(state.states, state.states, state.coefficients, overlaps)
    overlap = Symmetric((gram + transpose(gram)) ./ 2)
    coefficients = _canonicalize!(
        state.coefficients * inv(sqrt(overlap)), state.states, bundles, weight_atol)
    error = _matrix_identity_error(transpose(coefficients) *
                 _support_action(state.states, state.states, coefficients, overlaps))
    error <= identity_atol ||
        throw(ArgumentError("terminal PQS realized shell overlap is not identity"))
    return state.indices, state.states, coefficients
end
function pqs_terminal_basis_realization(
    support_records, retained_records, transform_contracts, bundles;
    identity_atol::Real = 1.0e-8,
    weight_atol::Real = 1.0e-14,
)
    length(support_records) == length(retained_records) ||
        throw(ArgumentError("terminal support and retained records must match"))
    contracts = Dict(contract.unit_key => contract for contract in transform_contracts)
    overlaps = (_nested_axis_pgdg(bundles, :x).overlap,
        _nested_axis_pgdg(bundles, :y).overlap,
        _nested_axis_pgdg(bundles, :z).overlap)
    blocks = CartesianTerminalBasisBlock[]
    seen_support = Set{Int}()
    nextcol = 1
    for (support, record) in zip(support_records, retained_records)
        support.unit_key == record.support_record.unit_key ||
            throw(ArgumentError("terminal support and retained record order mismatch"))
        if record.transform_kind === :direct_identity_transform_contract
            nextcol = _append_direct!(
                blocks, record, bundles, overlaps, nextcol, identity_atol)
            _validate_block_support!(last(blocks), support, seen_support)
        elseif record.transform_kind === :pqs_source_modes_boundary_selection_shell_realization_contract
            contract = get(contracts, record.transform_contract_unit_key, nothing)
            isnothing(contract) && throw(ArgumentError("missing terminal transform contract"))
            indices, states, coefficients = _realize_shell(
                record, contract, bundles, overlaps, identity_atol, weight_atol)
            nextcol = _push_block!(
                blocks, record.support_record.unit_key, indices, states,
                coefficients, nextcol)
            _validate_block_support!(last(blocks), support, seen_support)
        elseif record.transform_kind === :compact_thin_slab_face_product_contract
            contract = get(contracts, record.transform_contract_unit_key, nothing)
            isnothing(contract) && throw(ArgumentError("missing terminal transform contract"))
            indices, states, coefficients =
                _terminal_compact_thin_slab_block(contract.source_cpbs, contract.metadata, bundles)
            error = _matrix_identity_error(transpose(coefficients) *
                    _support_action(states, states, coefficients, overlaps))
            error <= identity_atol ||
                throw(ArgumentError("terminal compact thin slab overlap is not identity"))
            nextcol = _push_block!(
                blocks, record.support_record.unit_key, indices, states,
                coefficients, nextcol)
            _validate_block_support!(last(blocks), support, seen_support)
        else
            throw(ArgumentError("unsupported terminal transform kind for PQS basis realization"))
        end
    end
    return CartesianTerminalBasisRealization(blocks, nextcol - 1, 0.0)
end

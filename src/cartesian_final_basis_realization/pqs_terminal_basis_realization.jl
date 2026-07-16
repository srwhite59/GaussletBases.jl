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
struct CartesianParentBackedInjectedComposition
    terminal_basis::CartesianTerminalBasisRealization
    parent_residual_blocks::Vector{CartesianParentResidualFunctionBlock}
    parent_residual_column_ranges::Vector{UnitRange{Int}}
    parent_backed_dimension::Int
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

function _validate_parent_backed_injected_composition(
    composition::CartesianParentBackedInjectedComposition,
    bundles;
    orthogonality_atol::Real = 1.0e-10,
    identity_atol::Real = 1.0e-8,
)
    basis = composition.terminal_basis
    next_column = 1
    for block in basis.blocks
        count = length(block.column_range)
        block.column_range == next_column:(next_column + count - 1) ||
            throw(ArgumentError(
                "parent-backed terminal column ranges are not contiguous and native ordered"))
        length(block.support_indices) == length(block.support_states) ||
            throw(DimensionMismatch(
                "parent-backed terminal support indices and states differ"))
        if isnothing(block.coefficients)
            length(block.support_indices) == count || throw(DimensionMismatch(
                "parent-backed terminal identity block dimensions differ"))
        else
            size(block.coefficients) == (length(block.support_indices), count) ||
                throw(DimensionMismatch(
                    "parent-backed terminal coefficient dimensions differ"))
            all(isfinite, block.coefficients) || throw(ArgumentError(
                "parent-backed terminal coefficients must be finite"))
        end
        next_column += count
    end
    next_column - 1 == basis.final_dimension || throw(DimensionMismatch(
        "parent-backed terminal dimension does not match native column ranges"))

    prfs = composition.parent_residual_blocks
    source_keys = Symbol[prf.source_unit_key for prf in prfs]
    length(unique(source_keys)) == length(source_keys) || throw(ArgumentError(
        "parent-backed PRF source identities must be unique"))
    support_identities = [(prf.support_indices, prf.support_states) for prf in prfs]
    length(unique(support_identities)) == length(support_identities) ||
        throw(ArgumentError("parent-backed PRF support identities must be unique"))
    ranges, geometry = _validate_parent_residual_collection(
        basis, bundles, prfs, orthogonality_atol, identity_atol)
    composition.parent_residual_column_ranges == ranges || throw(ArgumentError(
        "parent-backed PRF column ranges are stale or reordered"))
    prf_dimension = sum(size(prf.coefficients, 2) for prf in prfs)
    composition.parent_backed_dimension == basis.final_dimension + prf_dimension ||
        throw(DimensionMismatch(
            "parent-backed dimension does not match terminal and PRF dimensions"))

    weights = Tuple(_nested_axis_pgdg(bundles, axis).weights for axis in (:x, :y, :z))
    final_weights = reduce(vcat,
        (_terminal_ida_final_weights(block, weights, 1.0e-12) for block in basis.blocks))
    return (; terminal_dimension = basis.final_dimension, prf_dimension,
        parent_backed_dimension = composition.parent_backed_dimension,
        prf_column_ranges = ranges, geometry,
        minimum_final_ida_weight = minimum(final_weights),
        maximum_final_ida_weight = maximum(final_weights))
end

function _terminal_block_parent_coefficients(block::CartesianTerminalBasisBlock)
    count = length(block.column_range)
    if isnothing(block.coefficients)
        length(block.support_indices) == count || throw(DimensionMismatch(
            "terminal identity block support and column counts differ"))
        return Matrix{Float64}(I, count, count)
    end
    size(block.coefficients) == (length(block.support_indices), count) ||
        throw(DimensionMismatch("terminal block coefficient dimensions are inconsistent"))
    return copy(block.coefficients)
end

function _parent_backed_injection_lowdin(columns, support_metric, label)
    raw_gram = transpose(columns) * support_metric * columns
    gram = Symmetric((raw_gram + transpose(raw_gram)) ./ 2)
    values = eigvals(gram)
    all(isfinite, values) || throw(ArgumentError("$label metric is not finite"))
    tolerance = max(1.0e-12, eps(Float64) * max(maximum(abs, values), 1.0) * length(values))
    minimum(values) > tolerance || throw(ArgumentError("$label loses numerical rank"))
    return columns * inv(sqrt(gram)), values
end

function build_parent_backed_injected_composition(
    basis::CartesianTerminalBasisRealization,
    bundles,
    requests,
)
    request_values = collect(requests)
    !isempty(request_values) || throw(ArgumentError(
        "parent-backed injection requires at least one request"))
    source_indices = Int[]
    old_prfs = CartesianParentResidualFunctionBlock[]
    request_prfs = Vector{Vector{CartesianParentResidualFunctionBlock}}()
    for request in request_values
        hasproperty(request, :source_block) && hasproperty(request, :prfs) &&
            hasproperty(request, :target_coordinates) || throw(ArgumentError(
            "parent-backed injection requests require source_block, prfs, and target_coordinates"))
        source = request.source_block
        index = findfirst(block -> block === source, basis.blocks)
        isnothing(index) && throw(ArgumentError(
            "parent-backed injection source must be an exact terminal block"))
        index in source_indices && throw(ArgumentError(
            "parent-backed injection source blocks must be unique"))
        push!(source_indices, index)
        prfs = CartesianParentResidualFunctionBlock[request.prfs...]
        !isempty(prfs) || throw(ArgumentError(
            "parent-backed injection request requires parent residual blocks"))
        all(prf -> prf.source_unit_key === source.unit_key &&
            prf.support_indices == source.support_indices &&
            prf.support_states == source.support_states, prfs) ||
            throw(ArgumentError(
                "parent-backed injection PRFs must match the exact source block"))
        append!(old_prfs, prfs)
        push!(request_prfs, prfs)
    end
    for right in eachindex(old_prfs), left in firstindex(old_prfs):(right - 1)
        old_prfs[left] === old_prfs[right] && throw(ArgumentError(
            "parent-backed injection PRF blocks must be unique"))
    end
    _validate_parent_residual_collection(basis, bundles, old_prfs, 1.0e-10, 1.0e-8)

    pgdg = Tuple(_nested_axis_pgdg(bundles, axis) for axis in (:x, :y, :z))
    overlaps = Tuple(axis.overlap for axis in pgdg)
    weights = Tuple(axis.weights for axis in pgdg)
    replacement_blocks = Dict{Int,CartesianTerminalBasisBlock}()
    raw_complements = Matrix{Float64}[]
    request_contexts = NamedTuple[]
    for (request_index, (request, prfs, block_index)) in enumerate(
            zip(request_values, request_prfs, source_indices))
        source = basis.blocks[block_index]
        G = _terminal_block_parent_coefficients(source)
        R = hcat((prf.coefficients for prf in prfs)...)
        old_span = hcat(G, R)
        support_metric = zeros(Float64, length(source.support_states),
            length(source.support_states))
        _support_cross!(support_metric, source.support_states,
            source.support_states, overlaps)
        old_gram = transpose(old_span) * support_metric * old_span
        norm(old_gram - I, Inf) <= 1.0e-8 || throw(ArgumentError(
            "parent-backed injection source span is not metric-orthonormal"))

        target_coordinates = Matrix{Float64}(request.target_coordinates)
        size(target_coordinates, 1) == size(old_span, 2) ||
            throw(DimensionMismatch(
                "parent-backed injection target coordinate rows must match [G,R]"))
        0 < size(target_coordinates, 2) <= size(G, 2) || throw(ArgumentError(
            "parent-backed injection target count must be positive and no larger than G"))
        all(isfinite, target_coordinates) || throw(ArgumentError(
            "parent-backed injection target coordinates must be finite"))
        norm(transpose(target_coordinates) * target_coordinates - I, Inf) <= 1.0e-10 ||
            throw(ArgumentError(
                "parent-backed injection targets must be orthonormal in [G,R] coordinates"))
        Y = old_span * target_coordinates
        C = transpose(G) * support_metric * Y
        singulars = svdvals(C)
        minimum(singulars) > 1.0e-10 || throw(ArgumentError(
            "parent-backed injection targets do not have full rank in the terminal block"))
        Qperp = nullspace(transpose(C); atol = 1.0e-12, rtol = 1.0e-10)
        size(Qperp, 2) == size(G, 2) - size(Y, 2) || throw(ArgumentError(
            "parent-backed injection terminal complement has the wrong dimension"))
        F = hcat(Y, G * Qperp)
        norm(transpose(F) * support_metric * F - I, Inf) <= 1.0e-8 ||
            throw(ArgumentError("parent-backed injection scaffold is not orthonormal"))
        projected_old = F * (transpose(F) * support_metric * G)
        Ginj, projected_values = _parent_backed_injection_lowdin(
            projected_old, support_metric, "parent-backed projected terminal seeds")

        coordinates = transpose(old_span) * support_metric * Ginj
        complement_coordinates = nullspace(transpose(coordinates);
            atol = 1.0e-12, rtol = 1.0e-10)
        old_prf_count = size(R, 2)
        size(complement_coordinates, 2) == old_prf_count || throw(ArgumentError(
            "parent-backed injection changed the parent residual count"))
        raw_complement, _ = _parent_backed_injection_lowdin(
            old_span * complement_coordinates, support_metric,
            "parent-backed exact complement")

        replacement = CartesianTerminalBasisBlock(source.unit_key,
            copy(source.support_indices), copy(source.support_states), Ginj,
            source.column_range)
        final_weights = _terminal_ida_final_weights(replacement, weights, 1.0e-12)
        replacement_blocks[block_index] = replacement
        push!(raw_complements, raw_complement)
        push!(request_contexts, (; request_index,
            source_unit_key = source.unit_key,
            terminal_column_range = source.column_range,
            terminal_count = size(G, 2), prf_count = old_prf_count,
            target_count = size(Y, 2), target_projection_singular_values = singulars,
            projected_seed_gram_eigenvalues = projected_values,
            old_span, Y, Ginj, support_metric,
            minimum_final_ida_weight = minimum(final_weights),
            maximum_final_ida_weight = maximum(final_weights)))
    end

    blocks = copy(basis.blocks)
    for index in source_indices
        blocks[index] = replacement_blocks[index]
    end
    injected_basis = CartesianTerminalBasisRealization(
        blocks, basis.final_dimension, basis.max_cross_overlap)
    new_prfs = CartesianParentResidualFunctionBlock[]
    for (block_index, raw_complement) in zip(source_indices, raw_complements)
        source = injected_basis.blocks[block_index]
        push!(new_prfs, build_parent_residual_function_block(
            injected_basis, source, bundles, source.support_indices, raw_complement))
    end
    all(index -> blocks[index] === basis.blocks[index],
        setdiff(eachindex(blocks), source_indices)) || throw(ArgumentError(
        "parent-backed injection changed an unaffected terminal block"))
    parent_backed_dimension = basis.final_dimension + sum(
        size(prf.coefficients, 2) for prf in new_prfs)
    ranges, geometry = _validate_parent_residual_collection(
        injected_basis, bundles, new_prfs, 1.0e-10, 1.0e-8)
    provisional = CartesianParentBackedInjectedComposition(injected_basis, new_prfs,
        ranges, parent_backed_dimension, (;))
    validation = _validate_parent_backed_injected_composition(provisional, bundles)
    request_diagnostics = NamedTuple[]
    for (context, prf) in zip(request_contexts, new_prfs)
        final_complement = prf.coefficients
        new_span = hcat(context.Ginj, final_complement)
        target_recovery = svdvals(
            transpose(context.Y) * context.support_metric * context.Ginj)
        span_singulars = svdvals(
            transpose(context.old_span) * context.support_metric * new_span)
        g_identity_error = norm(
            transpose(context.Ginj) * context.support_metric * context.Ginj - I, Inf)
        r_identity_error = norm(
            transpose(final_complement) * context.support_metric * final_complement - I, Inf)
        cross_error = norm(
            transpose(context.Ginj) * context.support_metric * final_complement, Inf)
        maximum(abs.(1.0 .- target_recovery)) <= 1.0e-8 || throw(ArgumentError(
            "parent-backed injection target recovery failed after final PRF rebuild"))
        maximum(abs.(1.0 .- span_singulars)) <= 1.0e-8 || throw(ArgumentError(
            "parent-backed injection final columns did not preserve the combined span"))
        max(g_identity_error, r_identity_error) <= 1.0e-8 &&
            cross_error <= 1.0e-10 || throw(ArgumentError(
            "parent-backed final terminal/complement metric validation failed"))
        push!(request_diagnostics, (; request_index = context.request_index,
            source_unit_key = context.source_unit_key,
            terminal_column_range = context.terminal_column_range,
            terminal_count = context.terminal_count, prf_count = context.prf_count,
            target_count = context.target_count,
            target_projection_singular_values = context.target_projection_singular_values,
            projected_seed_gram_eigenvalues = context.projected_seed_gram_eigenvalues,
            target_recovery_singular_values = target_recovery,
            old_new_span_singular_values = span_singulars,
            terminal_identity_error = g_identity_error,
            complement_identity_error = r_identity_error,
            terminal_complement_cross_error = cross_error,
            minimum_final_ida_weight = context.minimum_final_ida_weight,
            maximum_final_ida_weight = context.maximum_final_ida_weight))
    end
    diagnostics = (; dimensions = (;
            terminal_dimension = basis.final_dimension,
            parent_residual_dimension = parent_backed_dimension - basis.final_dimension,
            parent_backed_dimension,
            request_count = length(request_values)),
        source_block_indices = source_indices,
        requests = request_diagnostics,
        geometry = validation.geometry)
    return CartesianParentBackedInjectedComposition(injected_basis, new_prfs,
        ranges, parent_backed_dimension, diagnostics)
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

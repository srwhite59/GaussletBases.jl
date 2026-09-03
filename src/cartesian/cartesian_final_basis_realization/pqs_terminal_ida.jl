function _check_terminal_ida_weights(weights, n, name)
    length(weights) >= n ||
        throw(DimensionMismatch("terminal IDA $name weights are too short"))
    all(isfinite, weights) ||
        throw(ArgumentError("terminal IDA $name weights are not finite"))
    return weights
end

function _terminal_ida_support_weights(block, weights)
    wx, wy, wz = weights
    return Float64[
        Float64(wx[state[1]]) * Float64(wy[state[2]]) * Float64(wz[state[3]])
        for state in block.support_states
    ]
end

function _terminal_ida_final_weights(block, weights, weight_atol)
    support_weights = _terminal_ida_support_weights(block, weights)
    final_weights = isnothing(block.coefficients) ?
        support_weights :
        vec(transpose(block.coefficients) * support_weights)
    all(weight -> isfinite(weight) && weight > weight_atol, final_weights) ||
        throw(ArgumentError("terminal IDA final weights must be finite and positive"))
    return final_weights
end

function _normalize_terminal_ida_block!(block, left_weights, right_weights)
    size(block) == (length(left_weights), length(right_weights)) ||
        throw(DimensionMismatch("terminal IDA block and weight dimensions differ"))
    @inbounds for col in axes(block, 2), row in axes(block, 1)
        block[row, col] /= left_weights[row] * right_weights[col]
    end
    return block
end

function assemble_terminal_ida_interaction!(
    destination,
    basis::CartesianTerminalBasisRealization,
    coefficients,
    raw_pair_terms_x,
    raw_pair_terms_y,
    raw_pair_terms_z,
    weights_x,
    weights_y,
    weights_z;
    weight_atol = 1.0e-12,
    symmetry_atol = 1.0e-10,
)
    size(destination) == (basis.final_dimension, basis.final_dimension) ||
        throw(DimensionMismatch("terminal IDA destination size is wrong"))
    coeffs = Float64.(coefficients)
    !isempty(coeffs) && all(isfinite, coeffs) ||
        throw(ArgumentError("terminal IDA coefficients must be finite"))
    max_state = (0, 0, 0)
    for block in basis.blocks, state in block.support_states
        max_state = ntuple(axis -> max(max_state[axis], state[axis]), 3)
    end
    nterms = length(coeffs)
    raw_terms = (
        _check_terminal_factor_terms(_terminal_factor_terms(raw_pair_terms_x), nterms, max_state[1], "raw pair x"),
        _check_terminal_factor_terms(_terminal_factor_terms(raw_pair_terms_y), nterms, max_state[2], "raw pair y"),
        _check_terminal_factor_terms(_terminal_factor_terms(raw_pair_terms_z), nterms, max_state[3], "raw pair z"),
    )
    weights = (
        _check_terminal_ida_weights(weights_x, max_state[1], "x"),
        _check_terminal_ida_weights(weights_y, max_state[2], "y"),
        _check_terminal_ida_weights(weights_z, max_state[3], "z"),
    )
    block_weights = [
        _terminal_ida_final_weights(block, weights, Float64(weight_atol))
        for block in basis.blocks
    ]
    for right in eachindex(basis.blocks), left in 1:right
        left_block = basis.blocks[left]
        right_block = basis.blocks[right]
        rows = left_block.column_range
        cols = right_block.column_range
        block = _terminal_gaussian_sum_action(left_block, right_block, coeffs, raw_terms)
        _normalize_terminal_ida_block!(block, block_weights[left], block_weights[right])
        if left == right
            norm(block - transpose(block), Inf) <= symmetry_atol ||
                throw(ArgumentError("terminal IDA diagonal block is not symmetric"))
        end
        destination[rows, cols] .+= block
        left != right && (destination[cols, rows] .+= transpose(block))
    end
    norm(destination - transpose(destination), Inf) <= symmetry_atol ||
        throw(ArgumentError("terminal IDA destination is not symmetric"))
    return destination
end

function parent_gaussian_direct_resource(bundles, expansion::CoulombGaussianExpansion)
    _r3_validate_pgdg_expansion(bundles, expansion)
    pgdg = Tuple(_nested_axis_pgdg(bundles, axis) for axis in (:x, :y, :z))
    dims = _nested_axis_lengths(bundles)
    nterms = length(expansion)
    raw_terms = ntuple(axis -> _check_terminal_factor_terms(
        pgdg[axis].pair_factor_terms_raw, nterms, dims[axis], "parent raw pair $(axis)"), 3)
    weights = ntuple(axis -> _check_terminal_ida_weights(
        pgdg[axis].weights, dims[axis], "parent $(axis)"), 3)
    axis_centers = ntuple(axis -> Float64.(pgdg[axis].centers), 3)
    all(axis -> length(axis_centers[axis]) == dims[axis], 1:3) ||
        throw(DimensionMismatch("parent Gaussian-direct center dimensions are inconsistent"))
    centers = Vector{NTuple{3,Float64}}(undef, prod(dims))
    onsite = Vector{Float64}(undef, prod(dims))
    for parent_index in eachindex(centers)
        state = _cartesian_unflat_index(parent_index, dims)
        centers[parent_index] = ntuple(axis -> axis_centers[axis][state[axis]], 3)
        weight = prod(weights[axis][state[axis]] for axis in 1:3)
        value = 0.0
        @inbounds for term in 1:nterms
            value += expansion.coefficients[term] * prod(
                raw_terms[axis][term, state[axis], state[axis]] for axis in 1:3)
        end
        onsite[parent_index] = value / abs2(weight)
    end
    return _parent_gaussian_direct_resource(centers, onsite, expansion)
end

function _parent_backed_charge_matrix(coefficients, charge_atol, name)
    values = Matrix{Float64}(coefficients)
    all(isfinite, values) || throw(ArgumentError("$name coefficients must be finite"))
    charges = abs2.(values)
    sums = vec(sum(charges; dims = 1))
    all(value -> isfinite(value) && abs(value - 1.0) <= charge_atol, sums) ||
        throw(ArgumentError("$name parent charges must sum to one without renormalization"))
    return charges, sums
end

function _parent_direct_charge_block(
    left_indices, left_charges, right_indices, right_charges, kernel;
    tile_size::Int = 128,
)
    nleft = isnothing(left_charges) ? length(left_indices) : size(left_charges, 2)
    result = zeros(Float64, nleft, size(right_charges, 2))
    for left_first in 1:tile_size:length(left_indices)
        left_rows = left_first:min(left_first + tile_size - 1, length(left_indices))
        for right_first in 1:tile_size:length(right_indices)
            right_rows = right_first:min(right_first + tile_size - 1, length(right_indices))
            K = Matrix{Float64}(undef, length(left_rows), length(right_rows))
            @inbounds for (jlocal, j) in pairs(right_rows), (ilocal, i) in pairs(left_rows)
                K[ilocal, jlocal] = kernel(left_indices[i], right_indices[j])
            end
            action = K * @view(right_charges[right_rows, :])
            if isnothing(left_charges)
                result[left_rows, :] .+= action
            else
                result .+= transpose(@view(left_charges[left_rows, :])) * action
            end
        end
    end
    return result
end

function _parent_backed_direct_blocks(
    basis, prfs, bundles, kernel;
    charge_atol::Real = 5.0e-8,
    orthogonality_atol::Real = 1.0e-10,
    identity_atol::Real = 1.0e-8,
    symmetry_atol::Real = 1.0e-10,
)
    tolerance = Float64(charge_atol)
    symmetry_tolerance = Float64(symmetry_atol)
    all(value -> isfinite(value) && value >= 0.0,
        (tolerance, symmetry_tolerance)) || throw(ArgumentError(
        "parent-backed direct tolerances must be finite and nonnegative"))
    ranges, geometry = _validate_parent_residual_collection(
        basis, bundles, prfs, Float64(orthogonality_atol), Float64(identity_atol))
    prf_charges = Matrix{Float64}[]
    prf_charge_sums = Vector{Float64}[]
    for (index, prf) in pairs(prfs)
        charges, sums = _parent_backed_charge_matrix(
            prf.coefficients, tolerance, "parent residual block $index")
        push!(prf_charges, charges)
        push!(prf_charge_sums, sums)
    end
    terminal_charge_sums = zeros(Float64, basis.final_dimension)
    G_R = zeros(Float64, basis.final_dimension, last(ranges[end]))
    for block in basis.blocks
        block_charges = isnothing(block.coefficients) ? nothing :
            first(_parent_backed_charge_matrix(
                block.coefficients, tolerance, "terminal block $(block.unit_key)"))
        terminal_charge_sums[block.column_range] .= isnothing(block_charges) ?
            1.0 : vec(sum(block_charges; dims = 1))
        for index in eachindex(prfs)
            G_R[block.column_range, ranges[index]] .= _parent_direct_charge_block(
                block.support_indices, block_charges,
                prfs[index].support_indices, prf_charges[index], kernel)
        end
    end
    R_R = zeros(Float64, size(G_R, 2), size(G_R, 2))
    for right in eachindex(prfs), left in 1:right
        block = _parent_direct_charge_block(
            prfs[left].support_indices, prf_charges[left],
            prfs[right].support_indices, prf_charges[right], kernel)
        if left == right
            norm(block - transpose(block), Inf) <= symmetry_tolerance ||
                throw(ArgumentError("parent-backed direct diagonal block is not symmetric"))
            block = (block + transpose(block)) ./ 2
        end
        R_R[ranges[left], ranges[right]] .= block
        left != right && (R_R[ranges[right], ranges[left]] .= transpose(block))
    end
    all(isfinite, G_R) && all(isfinite, R_R) || throw(ArgumentError(
        "parent-backed direct blocks must be finite"))
    return (; G_R, R_R, prf_column_ranges = ranges,
        terminal_charge_sums, prf_charge_sums,
        diagnostics = (; geometry, charge_atol = tolerance))
end

function parent_gaussian_direct_blocks(
    basis::CartesianTerminalBasisRealization,
    prfs::AbstractVector{<:CartesianParentResidualFunctionBlock},
    bundles,
    resource::_ParentGaussianDirectResource,
    expansion::CoulombGaussianExpansion,
)
    _r3_validate_pgdg_expansion(bundles, expansion)
    dims = _nested_axis_lengths(bundles)
    prod(dims) == length(resource.centers) ||
        throw(DimensionMismatch("parent Gaussian-direct resource has the wrong parent size"))
    all(value -> isfinite(value) && value > 0.0, resource.onsite_values) &&
        resource.gaussian_exponents == (pi / 2) .* abs2.(resource.onsite_values) ||
        throw(ArgumentError("parent Gaussian-direct resource onsite calibration is invalid"))
    resource.validation.term_count == length(expansion) &&
        resource.coulomb_expansion_fingerprint ==
            _coulomb_expansion_fingerprint(expansion) || throw(ArgumentError(
        "parent Gaussian-direct resource Coulomb expansion identity is inconsistent"))
    pgdg = Tuple(_nested_axis_pgdg(bundles, axis) for axis in (:x, :y, :z))
    for index in eachindex(resource.centers)
        state = _cartesian_unflat_index(index, dims)
        expected = ntuple(axis -> Float64(pgdg[axis].centers[state[axis]]), 3)
        resource.centers[index] == expected || throw(ArgumentError(
            "parent Gaussian-direct resource centers do not match the resolved parent"))
    end
    for block in basis.blocks, (index, state) in zip(block.support_indices, block.support_states)
        _cartesian_unflat_index(index, dims) == state ||
            throw(ArgumentError("terminal block parent index/state mapping is inconsistent"))
    end
    for prf in prfs, (index, state) in zip(prf.support_indices, prf.support_states)
        _cartesian_unflat_index(index, dims) == state ||
            throw(ArgumentError("parent residual index/state mapping is inconsistent"))
    end
    kernel = (left, right) -> _parent_gaussian_direct_value(resource, left, right)
    result = _parent_backed_direct_blocks(basis, prfs, bundles, kernel)
    return merge(result, (; resource_validation = resource.validation,
        coulomb_expansion_fingerprint = resource.coulomb_expansion_fingerprint))
end

function parent_backed_injected_interaction_base_blocks(
    composition::CartesianParentBackedInjectedComposition,
    bundles;
    expansion,
)
    validation = _validate_parent_backed_injected_composition(composition, bundles)
    _r3_validate_pgdg_expansion(bundles, expansion)
    basis = composition.terminal_basis
    pgdg = Tuple(_nested_axis_pgdg(bundles, axis) for axis in (:x, :y, :z))
    terminal = zeros(Float64, basis.final_dimension, basis.final_dimension)
    assemble_terminal_ida_interaction!(terminal, basis, expansion.coefficients,
        pgdg[1].pair_factor_terms_raw, pgdg[2].pair_factor_terms_raw,
        pgdg[3].pair_factor_terms_raw,
        pgdg[1].weights, pgdg[2].weights, pgdg[3].weights)

    resource = parent_gaussian_direct_resource(bundles, expansion)
    direct = parent_gaussian_direct_blocks(basis,
        composition.parent_residual_blocks, bundles, resource, expansion)
    direct.prf_column_ranges == composition.parent_residual_column_ranges ||
        throw(ArgumentError("parent-backed direct PRF ranges differ from the composition"))
    nG = basis.final_dimension
    nR = composition.parent_backed_dimension - nG
    size(direct.G_R) == (nG, nR) && size(direct.R_R) == (nR, nR) ||
        throw(DimensionMismatch("parent-backed direct interaction block dimensions differ"))
    all(isfinite, terminal) && all(isfinite, direct.G_R) && all(isfinite, direct.R_R) ||
        throw(ArgumentError("parent-backed base interaction blocks must be finite"))
    norm(terminal - transpose(terminal), Inf) <= 1.0e-10 &&
        norm(direct.R_R - transpose(direct.R_R), Inf) <= 1.0e-10 ||
        throw(ArgumentError("parent-backed base interaction diagonal blocks must be symmetric"))
    return (; terminal, terminal_parent_residual = direct.G_R,
        parent_residual = direct.R_R,
        diagnostics = (; validation, direct = direct.diagnostics,
            direct_resource = resource.validation,
            coulomb_expansion_fingerprint = direct.coulomb_expansion_fingerprint))
end

function parent_ida_direct_blocks(
    basis::CartesianTerminalBasisRealization,
    prfs::AbstractVector{<:CartesianParentResidualFunctionBlock},
    bundles,
    expansion::CoulombGaussianExpansion,
)
    _r3_validate_pgdg_expansion(bundles, expansion)
    pgdg = Tuple(_nested_axis_pgdg(bundles, axis) for axis in (:x, :y, :z))
    dims = _nested_axis_lengths(bundles)
    terms = Tuple(axis.pair_factor_terms for axis in pgdg)
    all(all(isfinite, term) for term in terms) || throw(ArgumentError(
        "parent IDA comparator factors must be finite"))
    residual_support = sum(length(prf.support_indices) for prf in prfs)
    terminal_support = sum(length(block.support_indices) for block in basis.blocks)
    pair_count = terminal_support * residual_support + residual_support^2
    max_pair_count = 100_000_000
    pair_count <= max_pair_count || throw(ArgumentError(
        "bounded parent IDA comparator pair count exceeds its validation limit"))
    states = NTuple{3,Int}[_cartesian_unflat_index(index, dims) for index in 1:prod(dims)]
    kernel = function(left, right)
        a, b = states[left], states[right]
        value = 0.0
        @inbounds for term in eachindex(expansion.coefficients)
            value += expansion.coefficients[term] * prod(
                terms[axis][term, a[axis], b[axis]] for axis in 1:3)
        end
        return value
    end
    result = _parent_backed_direct_blocks(basis, prfs, bundles, kernel)
    return merge(result, (; diagnostics = merge(
        result.diagnostics, (; pair_count, max_pair_count))))
end

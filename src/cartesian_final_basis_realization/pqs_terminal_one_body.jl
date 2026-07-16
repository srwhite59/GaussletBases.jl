function _check_terminal_axis_factor(axis, n, name)
    size(axis, 1) >= n && size(axis, 2) >= n ||
        throw(DimensionMismatch("terminal one-body $name factor is too small"))
    all(isfinite, axis) || throw(ArgumentError("terminal one-body $name factor is not finite"))
    norm(axis - transpose(axis), Inf) <= 1.0e-10 ||
        throw(ArgumentError("terminal one-body $name factor is not symmetric"))
    return axis
end

function _add_terminal_product_block!(
    destination, left, right, axes, factor, mirror, action_buffer, tile_buffer, block_buffer)
    ncols = isnothing(right.coefficients) ?
        length(right.support_states) :
        size(right.coefficients, 2)
    action = _buffer_view!(action_buffer, length(left.support_states), ncols)
    fill!(action, 0.0)
    maxcols = max(1, _TERMINAL_WORKSPACE_BYTES ÷ (8 * max(length(left.support_states), 1)))
    tile_buffer_view = _buffer_view!(
        tile_buffer, length(left.support_states), min(length(right.support_states), maxcols))
    for firstcol in 1:maxcols:length(right.support_states)
        cols = firstcol:min(firstcol + maxcols - 1, length(right.support_states))
        tile = @view tile_buffer_view[:, 1:length(cols)]
        _support_cross!(tile, left.support_states, @view(right.support_states[cols]), axes)
        if isnothing(right.coefficients)
            @view(action[:, cols]) .+= tile
        else
            mul!(action, tile, @view(right.coefficients[cols, :]), 1.0, 1.0)
        end
    end
    rows = left.column_range
    cols = right.column_range
    if isnothing(left.coefficients)
        destination[rows, cols] .+= factor .* action
        mirror && (destination[cols, rows] .+= factor .* transpose(action))
        return destination
    end
    block = _buffer_view!(block_buffer, size(left.coefficients, 2), ncols)
    mul!(block, transpose(left.coefficients), action)
    destination[rows, cols] .+= factor .* block
    mirror && (destination[cols, rows] .+= factor .* transpose(block))
    return destination
end

function _terminal_factor_terms(factors)
    factors isa AbstractArray{<:Real,3} && return factors
    isempty(factors) && throw(ArgumentError("terminal Gaussian factor packet is empty"))
    n = length(factors)
    dims = size(first(factors))
    terms = Array{Float64}(undef, n, dims...)
    for term in 1:n
        size(factors[term]) == dims ||
            throw(DimensionMismatch("terminal Gaussian factor matrices must match"))
        terms[term, :, :] .= factors[term]
    end
    return terms
end

function _check_terminal_factor_terms(terms, nterms, n, name)
    size(terms, 1) == nterms ||
        throw(DimensionMismatch("terminal Gaussian $name term count mismatch"))
    size(terms, 2) >= n && size(terms, 3) >= n ||
        throw(DimensionMismatch("terminal Gaussian $name factor is too small"))
    all(isfinite, terms) ||
        throw(ArgumentError("terminal Gaussian $name factor is not finite"))
    for term in 1:nterms
        matrix = @view terms[term, :, :]
        norm(matrix - transpose(matrix), Inf) <= 1.0e-10 ||
            throw(ArgumentError("terminal Gaussian $name factor term is not symmetric"))
    end
    return terms
end

function _fill_terminal_gaussian_sum_action!(
    action, left, right, coefficients, factor_terms, tile_buffer)
    fx, fy, fz = factor_terms
    nterms = length(coefficients)
    fill!(action, 0.0)
    maxcols = max(1, _TERMINAL_WORKSPACE_BYTES ÷ (8 * max(length(left.support_states), 1)))
    tile_buffer_view = _buffer_view!(
        tile_buffer, length(left.support_states), min(length(right.support_states), maxcols))
    for firstcol in 1:maxcols:length(right.support_states)
        cols = firstcol:min(firstcol + maxcols - 1, length(right.support_states))
        tile = @view tile_buffer_view[:, 1:length(cols)]
        @inbounds for (jj, j) in enumerate(cols), i in eachindex(left.support_states)
            ix, iy, iz = left.support_states[i]
            jx, jy, jz = right.support_states[j]
            value = 0.0
            for term in 1:nterms
                value += coefficients[term] * fx[term, ix, jx] *
                    fy[term, iy, jy] * fz[term, iz, jz]
            end
            tile[i, jj] = value
        end
        if isnothing(right.coefficients)
            @views action[:, cols] .+= tile
        else
            mul!(action, tile, @view(right.coefficients[cols, :]), 1.0, 1.0)
        end
    end
    return action
end

function _terminal_gaussian_sum_action(left, right, coefficients, factor_terms)
    ncols = isnothing(right.coefficients) ?
        length(right.support_states) :
        size(right.coefficients, 2)
    action = zeros(Float64, length(left.support_states), ncols)
    tile_buffer = Ref(Matrix{Float64}(undef, 0, 0))
    _fill_terminal_gaussian_sum_action!(
        action, left, right, coefficients, factor_terms, tile_buffer)
    if isnothing(left.coefficients)
        return action
    end
    block = Matrix{Float64}(undef, size(left.coefficients, 2), ncols)
    mul!(block, transpose(left.coefficients), action)
    return block
end

function _add_terminal_gaussian_sum_block!(
    destination, left, right, coefficients, factor_terms, factor, mirror,
    action_buffer, tile_buffer, block_buffer)
    ncols = isnothing(right.coefficients) ?
        length(right.support_states) :
        size(right.coefficients, 2)
    action = _buffer_view!(action_buffer, length(left.support_states), ncols)
    _fill_terminal_gaussian_sum_action!(
        action, left, right, coefficients, factor_terms, tile_buffer)
    rows = left.column_range
    cols = right.column_range
    if isnothing(left.coefficients)
        @views destination[rows, cols] .+= factor .* action
        mirror && @views destination[cols, rows] .+= factor .* transpose(action)
        return destination
    end
    block = _buffer_view!(block_buffer, size(left.coefficients, 2), ncols)
    mul!(block, transpose(left.coefficients), action)
    @views destination[rows, cols] .+= factor .* block
    mirror && @views destination[cols, rows] .+= factor .* transpose(block)
    return destination
end

function _accumulate_terminal_gaussian_sum!(
    destination,
    basis::CartesianTerminalBasisRealization,
    coefficients,
    factors_x,
    factors_y,
    factors_z;
    scale = -1.0,
)
    action_buffer = Ref(Matrix{Float64}(undef, 0, 0))
    tile_buffer = Ref(Matrix{Float64}(undef, 0, 0))
    block_buffer = Ref(Matrix{Float64}(undef, 0, 0))
    return _accumulate_terminal_gaussian_sum!(
        destination, basis, coefficients, factors_x, factors_y, factors_z,
        action_buffer, tile_buffer, block_buffer; scale = scale)
end

function _accumulate_terminal_gaussian_sum!(
    destination,
    basis::CartesianTerminalBasisRealization,
    coefficients,
    factors_x,
    factors_y,
    factors_z,
    action_buffer,
    tile_buffer,
    block_buffer;
    scale = -1.0,
)
    size(destination) == (basis.final_dimension, basis.final_dimension) ||
        throw(DimensionMismatch("terminal Gaussian-sum destination size is wrong"))
    max_state = (0, 0, 0)
    for block in basis.blocks, state in block.support_states
        max_state = ntuple(axis -> max(max_state[axis], state[axis]), 3)
    end
    coeffs = Float64.(coefficients)
    all(isfinite, coeffs) ||
        throw(ArgumentError("terminal Gaussian-sum coefficients must be finite"))
    nterms = length(coeffs)
    factor_terms = (
        _check_terminal_factor_terms(_terminal_factor_terms(factors_x), nterms, max_state[1], "x"),
        _check_terminal_factor_terms(_terminal_factor_terms(factors_y), nterms, max_state[2], "y"),
        _check_terminal_factor_terms(_terminal_factor_terms(factors_z), nterms, max_state[3], "z"),
    )
    factor = Float64(scale)
    for right in eachindex(basis.blocks), left in 1:right
        left_block = basis.blocks[left]
        right_block = basis.blocks[right]
        _add_terminal_gaussian_sum_block!(
            destination, left_block, right_block, coeffs, factor_terms, factor, left != right,
            action_buffer, tile_buffer, block_buffer)
    end
    return destination
end

function assemble_terminal_product_operator!(
    destination,
    basis::CartesianTerminalBasisRealization,
    axis_x,
    axis_y,
    axis_z;
    scale = 1.0,
)
    action_buffer = Ref(Matrix{Float64}(undef, 0, 0))
    tile_buffer = Ref(Matrix{Float64}(undef, 0, 0))
    block_buffer = Ref(Matrix{Float64}(undef, 0, 0))
    return _assemble_terminal_product_operator!(
        destination, basis, axis_x, axis_y, axis_z,
        action_buffer, tile_buffer, block_buffer; scale = scale)
end

function _assemble_terminal_product_operator!(
    destination,
    basis::CartesianTerminalBasisRealization,
    axis_x,
    axis_y,
    axis_z,
    action_buffer,
    tile_buffer,
    block_buffer;
    scale = 1.0,
)
    size(destination) == (basis.final_dimension, basis.final_dimension) ||
        throw(DimensionMismatch("terminal product destination size is wrong"))
    max_state = (0, 0, 0)
    for block in basis.blocks, state in block.support_states
        max_state = ntuple(axis -> max(max_state[axis], state[axis]), 3)
    end
    axes = (
        _check_terminal_axis_factor(axis_x, max_state[1], "x"),
        _check_terminal_axis_factor(axis_y, max_state[2], "y"),
        _check_terminal_axis_factor(axis_z, max_state[3], "z"),
    )
    factor = Float64(scale)
    for right in eachindex(basis.blocks), left in 1:right
        left_block = basis.blocks[left]
        right_block = basis.blocks[right]
        _add_terminal_product_block!(
            destination, left_block, right_block, axes, factor, left != right,
            action_buffer, tile_buffer, block_buffer)
    end
    return destination
end

_parent_residual_operator_block(prf::CartesianParentResidualFunctionBlock) = (;
    support_states = prf.support_states,
    coefficients = prf.coefficients,
    column_range = 1:size(prf.coefficients, 2))

function _parent_residual_product_blocks(
    basis::CartesianTerminalBasisRealization,
    prfs::AbstractVector{<:CartesianParentResidualFunctionBlock},
    ranges,
    axis_x,
    axis_y,
    axis_z;
    scale::Real = 1.0,
    symmetry_atol::Real = 1.0e-10,
)
    max_state = (0, 0, 0)
    for block in basis.blocks, state in block.support_states
        max_state = ntuple(axis -> max(max_state[axis], state[axis]), 3)
    end
    for prf in prfs, state in prf.support_states
        max_state = ntuple(axis -> max(max_state[axis], state[axis]), 3)
    end
    axes = (
        _check_terminal_axis_factor(axis_x, max_state[1], "x"),
        _check_terminal_axis_factor(axis_y, max_state[2], "y"),
        _check_terminal_axis_factor(axis_z, max_state[3], "z"),
    )
    operator_blocks = [_parent_residual_operator_block(prf) for prf in prfs]
    G_R = zeros(Float64, basis.final_dimension, last(ranges[end]))
    for (right_index, right) in pairs(operator_blocks), block in basis.blocks
        action = _support_action(
            block.support_states, right.support_states, right.coefficients, axes)
        G_R[block.column_range, ranges[right_index]] .= isnothing(block.coefficients) ?
            action : transpose(block.coefficients) * action
    end
    R_R = zeros(Float64, size(G_R, 2), size(G_R, 2))
    for right_index in eachindex(operator_blocks), left_index in 1:right_index
        left = operator_blocks[left_index]
        right = operator_blocks[right_index]
        block = transpose(left.coefficients) * _support_action(
            left.support_states, right.support_states, right.coefficients, axes)
        if left_index == right_index
            norm(block - transpose(block), Inf) <= symmetry_atol || throw(ArgumentError(
                "parent residual product diagonal block is not symmetric"))
            block = (block + transpose(block)) ./ 2
        end
        R_R[ranges[left_index], ranges[right_index]] .= block
        left_index != right_index &&
            (R_R[ranges[right_index], ranges[left_index]] .= transpose(block))
    end
    factor = Float64(scale)
    G_R .*= factor
    R_R .*= factor
    return (; G_R, R_R)
end

function _parent_residual_gaussian_sum_blocks(
    basis::CartesianTerminalBasisRealization,
    prfs::AbstractVector{<:CartesianParentResidualFunctionBlock},
    ranges,
    coefficients,
    factors_x,
    factors_y,
    factors_z;
    scale::Real = -1.0,
    symmetry_atol::Real = 1.0e-10,
)
    coeffs = Float64.(coefficients)
    !isempty(coeffs) && all(isfinite, coeffs) || throw(ArgumentError(
        "parent residual Gaussian-sum coefficients must be finite"))
    max_state = (0, 0, 0)
    for block in basis.blocks, state in block.support_states
        max_state = ntuple(axis -> max(max_state[axis], state[axis]), 3)
    end
    for prf in prfs, state in prf.support_states
        max_state = ntuple(axis -> max(max_state[axis], state[axis]), 3)
    end
    nterms = length(coeffs)
    factor_terms = (
        _check_terminal_factor_terms(_terminal_factor_terms(factors_x), nterms, max_state[1], "PRF x"),
        _check_terminal_factor_terms(_terminal_factor_terms(factors_y), nterms, max_state[2], "PRF y"),
        _check_terminal_factor_terms(_terminal_factor_terms(factors_z), nterms, max_state[3], "PRF z"),
    )
    operator_blocks = [_parent_residual_operator_block(prf) for prf in prfs]
    G_R = zeros(Float64, basis.final_dimension, last(ranges[end]))
    for (right_index, right) in pairs(operator_blocks), block in basis.blocks
        G_R[block.column_range, ranges[right_index]] .=
            _terminal_gaussian_sum_action(block, right, coeffs, factor_terms)
    end
    R_R = zeros(Float64, size(G_R, 2), size(G_R, 2))
    for right_index in eachindex(operator_blocks), left_index in 1:right_index
        block = _terminal_gaussian_sum_action(
            operator_blocks[left_index], operator_blocks[right_index], coeffs, factor_terms)
        if left_index == right_index
            norm(block - transpose(block), Inf) <= symmetry_atol || throw(ArgumentError(
                "parent residual Gaussian-sum diagonal block is not symmetric"))
            block = (block + transpose(block)) ./ 2
        end
        R_R[ranges[left_index], ranges[right_index]] .= block
        left_index != right_index &&
            (R_R[ranges[right_index], ranges[left_index]] .= transpose(block))
    end
    factor = Float64(scale)
    G_R .*= factor
    R_R .*= factor
    return (; G_R, R_R)
end

function parent_residual_one_body_blocks(
    basis::CartesianTerminalBasisRealization,
    bundles,
    prfs::AbstractVector{<:CartesianParentResidualFunctionBlock},
    atom_locations,
    nuclear_charges;
    expansion,
    orthogonality_atol::Real = 1.0e-10,
    identity_atol::Real = 1.0e-8,
)
    length(atom_locations) == length(nuclear_charges) || throw(DimensionMismatch(
        "parent residual atom locations and nuclear charges must match"))
    locations = NTuple{3,Float64}[Tuple(Float64.(location)) for location in atom_locations]
    charges = Float64.(nuclear_charges)
    all(location -> all(isfinite, location), locations) || throw(ArgumentError(
        "parent residual nuclear locations must be finite"))
    all(isfinite, charges) || throw(ArgumentError(
        "parent residual nuclear charges must be finite"))
    _r3_validate_pgdg_expansion(bundles, expansion)
    ranges, geometry = _validate_parent_residual_collection(
        basis, bundles, prfs, orthogonality_atol, identity_atol)
    pgdg = Tuple(_nested_axis_pgdg(bundles, axis) for axis in (:x, :y, :z))
    S = Tuple(axis.overlap for axis in pgdg)
    function add_blocks(left, right)
        left.G_R .+= right.G_R
        left.R_R .+= right.R_R
        return left
    end
    kinetic = _parent_residual_product_blocks(
        basis, prfs, ranges, pgdg[1].kinetic, S[2], S[3])
    kinetic = add_blocks(kinetic, _parent_residual_product_blocks(
        basis, prfs, ranges, S[1], pgdg[2].kinetic, S[3]))
    kinetic = add_blocks(kinetic, _parent_residual_product_blocks(
        basis, prfs, ranges, S[1], S[2], pgdg[3].kinetic))
    unit_nuclear = [begin
        factors = ntuple(axis -> _r3a_centered_factor_terms(
            pgdg[axis], expansion, location[axis]), 3)
        _parent_residual_gaussian_sum_blocks(
            basis, prfs, ranges, expansion.coefficients, factors...)
    end for location in locations]
    H1 = (; G_R = copy(kinetic.G_R), R_R = copy(kinetic.R_R))
    length(unit_nuclear) == length(charges) || throw(DimensionMismatch(
        "parent residual unit-nuclear center count differs from charges"))
    for index in eachindex(charges)
        H1.G_R .+= charges[index] .* unit_nuclear[index].G_R
        H1.R_R .+= charges[index] .* unit_nuclear[index].R_R
    end
    position = (;
        x = _parent_residual_product_blocks(basis, prfs, ranges, pgdg[1].position, S[2], S[3]),
        y = _parent_residual_product_blocks(basis, prfs, ranges, S[1], pgdg[2].position, S[3]),
        z = _parent_residual_product_blocks(basis, prfs, ranges, S[1], S[2], pgdg[3].position))
    x2 = (;
        x = _parent_residual_product_blocks(basis, prfs, ranges, pgdg[1].x2, S[2], S[3]),
        y = _parent_residual_product_blocks(basis, prfs, ranges, S[1], pgdg[2].x2, S[3]),
        z = _parent_residual_product_blocks(basis, prfs, ranges, S[1], S[2], pgdg[3].x2))
    blocks = Any[kinetic, H1, unit_nuclear...,
        position.x, position.y, position.z, x2.x, x2.y, x2.z]
    all(block -> all(isfinite, block.G_R) && all(isfinite, block.R_R), blocks) ||
        throw(ArgumentError("parent residual one-body blocks must be finite"))
    return (; kinetic, nuclear_attraction_unit_by_center = unit_nuclear, H1,
        position, x2, prf_column_ranges = ranges,
        diagnostics = (; terminal_dimension = basis.final_dimension,
            prf_dimension = last(ranges[end]), center_count = length(locations), geometry))
end

function _parent_backed_operator_matrix(GG, blocks, composition)
    nG = composition.terminal_basis.final_dimension
    nB = composition.parent_backed_dimension
    size(GG) == (nG, nG) || throw(DimensionMismatch(
        "parent-backed terminal operator dimension mismatch"))
    size(blocks.G_R) == (nG, nB - nG) &&
        size(blocks.R_R) == (nB - nG, nB - nG) || throw(DimensionMismatch(
        "parent-backed residual operator dimensions mismatch"))
    matrix = zeros(Float64, nB, nB)
    residual_range = (nG + 1):nB
    matrix[1:nG, 1:nG] .= GG
    matrix[1:nG, residual_range] .= blocks.G_R
    matrix[residual_range, 1:nG] .= transpose(blocks.G_R)
    matrix[residual_range, residual_range] .= blocks.R_R
    all(isfinite, matrix) || throw(ArgumentError(
        "parent-backed operator matrix must be finite"))
    norm(matrix - transpose(matrix), Inf) <= 1.0e-10 || throw(ArgumentError(
        "parent-backed operator matrix must be symmetric"))
    return Matrix{Float64}((matrix + transpose(matrix)) ./ 2)
end

function parent_backed_injected_one_body_operators(
    composition::CartesianParentBackedInjectedComposition,
    bundles,
    atom_locations,
    nuclear_charges;
    expansion,
)
    _validate_parent_backed_injected_composition(composition, bundles)
    basis = composition.terminal_basis
    charges = Float64.(nuclear_charges)
    length(atom_locations) == length(charges) || throw(DimensionMismatch(
        "parent-backed atom location and nuclear charge counts differ"))
    all(isfinite, charges) || throw(ArgumentError(
        "parent-backed nuclear charges must be finite"))
    pgdg = _r3_validate_pgdg_expansion(bundles, expansion)
    S = Tuple(axis.overlap for axis in pgdg)
    function terminal_product(axes...)
        matrix = zeros(Float64, basis.final_dimension, basis.final_dimension)
        assemble_terminal_product_operator!(matrix, basis, axes...)
        return matrix
    end
    terminal_kinetic = terminal_product(pgdg[1].kinetic, S[2], S[3])
    terminal_kinetic .+= terminal_product(S[1], pgdg[2].kinetic, S[3])
    terminal_kinetic .+= terminal_product(S[1], S[2], pgdg[3].kinetic)
    terminal_position = (;
        x = terminal_product(pgdg[1].position, S[2], S[3]),
        y = terminal_product(S[1], pgdg[2].position, S[3]),
        z = terminal_product(S[1], S[2], pgdg[3].position))
    terminal_x2 = (;
        x = terminal_product(pgdg[1].x2, S[2], S[3]),
        y = terminal_product(S[1], pgdg[2].x2, S[3]),
        z = terminal_product(S[1], S[2], pgdg[3].x2))
    locations = NTuple{3,Float64}[Tuple(Float64.(location)) for location in atom_locations]
    all(location -> all(isfinite, location), locations) || throw(ArgumentError(
        "parent-backed nuclear locations must be finite"))
    terminal_units = [begin
        factors = ntuple(axis -> _r3a_centered_factor_terms(
            pgdg[axis], expansion, location[axis]), 3)
        matrix = zeros(Float64, basis.final_dimension, basis.final_dimension)
        _accumulate_terminal_gaussian_sum!(matrix, basis,
            expansion.coefficients, factors...)
    end for location in locations]
    residual = parent_residual_one_body_blocks(basis, bundles,
        composition.parent_residual_blocks, locations, charges; expansion)
    length(terminal_units) == length(residual.nuclear_attraction_unit_by_center) ==
        length(locations) || throw(DimensionMismatch(
        "parent-backed terminal and PRF unit-nuclear center counts differ"))
    kinetic = _parent_backed_operator_matrix(
        terminal_kinetic, residual.kinetic, composition)
    unit_nuclear = Matrix{Float64}[_parent_backed_operator_matrix(
        terminal_units[index], residual.nuclear_attraction_unit_by_center[index],
        composition) for index in eachindex(locations)]
    position = (;
        x = _parent_backed_operator_matrix(terminal_position.x, residual.position.x, composition),
        y = _parent_backed_operator_matrix(terminal_position.y, residual.position.y, composition),
        z = _parent_backed_operator_matrix(terminal_position.z, residual.position.z, composition))
    x2 = (;
        x = _parent_backed_operator_matrix(terminal_x2.x, residual.x2.x, composition),
        y = _parent_backed_operator_matrix(terminal_x2.y, residual.x2.y, composition),
        z = _parent_backed_operator_matrix(terminal_x2.z, residual.x2.z, composition))
    H1 = copy(kinetic)
    for index in eachindex(charges)
        H1 .+= charges[index] .* unit_nuclear[index]
    end
    all(isfinite, H1) && norm(H1 - transpose(H1), Inf) <= 1.0e-10 ||
        throw(ArgumentError("parent-backed physical H1 must be finite and symmetric"))
    return (; kinetic, nuclear_attraction_unit_by_center = unit_nuclear,
        one_body_hamiltonian = H1, position, x2,
        diagnostics = (; terminal_dimension = basis.final_dimension,
            parent_residual_dimension = composition.parent_backed_dimension -
                basis.final_dimension,
            parent_backed_dimension = composition.parent_backed_dimension,
            center_count = length(locations)))
end

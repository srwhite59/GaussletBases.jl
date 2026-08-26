
function _pqs_source_box_route_driver_centered_factor_terms(pgdg, expansion, center)
    center == pgdg.center && return pgdg.gaussian_factor_terms
    ops = mapped_ordinary_one_body_operators(
        pgdg.basis; exponents = expansion.exponents, center, backend = pgdg.backend)
    return ops.gaussian_factors
end

function _pqs_source_box_route_driver_validate_pgdg_expansion(pgdg, expansion)
    length(pgdg) == 3 || throw(DimensionMismatch("PQS IDA requires three PGDG axes"))
    for (axis_name, axis) in zip((:x, :y, :z), pgdg)
        length(axis.exponents) == length(expansion) &&
            axis.exponents == expansion.exponents ||
            throw(ArgumentError("PQS $(axis_name)-axis Coulomb exponent sequence mismatch"))
    end
    return nothing
end

function _pqs_source_box_route_driver_terminal_products(terminal_basis_realization, pgdg)
    S = Tuple(axis.overlap for axis in pgdg)
    T = Tuple(axis.kinetic for axis in pgdg)
    n = terminal_basis_realization.final_dimension
    K = zeros(Float64, n, n)
    C = CartesianFinalBasisRealization
    buffers = C._terminal_operator_buffers(terminal_basis_realization)
    C._assemble_terminal_product_operator!(K, terminal_basis_realization, T[1], S[2], S[3], buffers...)
    C._assemble_terminal_product_operator!(K, terminal_basis_realization, S[1], T[2], S[3], buffers...)
    C._assemble_terminal_product_operator!(K, terminal_basis_realization, S[1], S[2], T[3], buffers...)
    return (; kinetic = K)
end

function _pqs_source_box_route_driver_terminal_unit_nuclear(terminal_basis_realization, expansion, pgdg, atom_locations::Vector{NTuple{3,Float64}})
    _pqs_source_box_route_driver_validate_pgdg_expansion(pgdg, expansion)
    n = terminal_basis_realization.final_dimension
    C = CartesianFinalBasisRealization
    U = Matrix{Float64}[]
    buffers = C._terminal_operator_buffers(terminal_basis_realization)
    for location in atom_locations
        matrix = zeros(Float64, n, n)
        factors = ntuple(axis -> _pqs_source_box_route_driver_centered_factor_terms(
            pgdg[axis], expansion, location[axis]), 3)
        C._accumulate_terminal_gaussian_sum!(
            matrix, terminal_basis_realization, expansion.coefficients,
            factors[1], factors[2], factors[3], buffers...)
        push!(U, matrix)
    end
    return U
end

function _pqs_source_box_route_driver_terminal_vee(terminal_basis_realization, expansion, pgdg)
    _pqs_source_box_route_driver_validate_pgdg_expansion(pgdg, expansion)
    V = zeros(Float64, terminal_basis_realization.final_dimension,
        terminal_basis_realization.final_dimension)
    CartesianFinalBasisRealization.assemble_terminal_ida_interaction!(
        V, terminal_basis_realization, expansion.coefficients,
        pgdg[1].pair_factor_terms_raw, pgdg[2].pair_factor_terms_raw,
        pgdg[3].pair_factor_terms_raw,
        pgdg[1].weights, pgdg[2].weights, pgdg[3].weights)
    return V
end

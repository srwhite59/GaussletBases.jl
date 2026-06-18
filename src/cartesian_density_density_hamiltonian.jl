struct _CartesianDensityDensityHamiltonian{T}
    one_body::Matrix{T}
    density_interaction::Matrix{T}
    orbital_to_density::Matrix{T}
    nup::Int
    ndn::Int
    constant_energy::T
    orbital_basis::Tuple
    density_basis::Tuple
    nuclear_charges::Vector{T}
    nuclear_positions::Matrix{T}
end

function _cartesian_nuclear_position_matrix(positions)
    rows = collect(positions)
    matrix = zeros(Float64, length(rows), 3)
    for (row, position) in pairs(rows)
        length(position) == 3 ||
            throw(DimensionMismatch("nuclear positions must have three coordinates"))
        matrix[row, :] .= Float64.(position)
    end
    return matrix
end

_cartesian_dense_float_matrix(matrix::Matrix{Float64}) = matrix
_cartesian_dense_float_matrix(matrix) = Matrix{Float64}(matrix)

_cartesian_float_vector(values::Vector{Float64}) = values

function _cartesian_float_vector(values)
    return Float64[Float64(value) for value in values]
end

function _cartesian_nuclear_position_matrix(positions::Matrix{Float64})
    size(positions, 2) == 3 ||
        throw(DimensionMismatch("nuclear position matrix must have three columns"))
    return positions
end

function _cartesian_nuclear_position_matrix(positions::AbstractMatrix{<:Real})
    size(positions, 2) == 3 ||
        throw(DimensionMismatch("nuclear position matrix must have three columns"))
    return _cartesian_dense_float_matrix(positions)
end

function _cartesian_check_symmetric_finite_matrix(
    name::AbstractString,
    matrix::AbstractMatrix{<:Real},
)
    dense = _cartesian_dense_float_matrix(matrix)
    size(dense, 1) == size(dense, 2) ||
        throw(DimensionMismatch("$(name) must be square"))
    all(isfinite, dense) ||
        throw(ArgumentError("$(name) contains non-finite entries"))
    norm(dense - transpose(dense), Inf) <= 1.0e-8 ||
        throw(ArgumentError("$(name) must be symmetric"))
    return dense
end

function _CartesianDensityDensityHamiltonian(
    one_body,
    density_interaction,
    orbital_to_density,
    nup::Integer,
    ndn::Integer,
    constant_energy::Real;
    orbital_basis,
    density_basis,
    nuclear_charges,
    nuclear_positions,
)
    one_body_matrix =
        _cartesian_check_symmetric_finite_matrix("one-body Hamiltonian", one_body)
    density_matrix =
        _cartesian_check_symmetric_finite_matrix(
            "density interaction",
            density_interaction,
        )
    transform = _cartesian_dense_float_matrix(orbital_to_density)
    size(transform) == (size(density_matrix, 1), size(one_body_matrix, 1)) ||
        throw(DimensionMismatch("orbital-to-density transform shape mismatch"))
    all(isfinite, transform) ||
        throw(ArgumentError("orbital-to-density transform contains non-finite entries"))
    nup_value = Int(nup)
    ndn_value = Int(ndn)
    nup_value >= 0 && ndn_value >= 0 && nup_value + ndn_value > 0 ||
        throw(ArgumentError("spin-sector electron counts must be non-negative and nonzero"))
    norb = size(one_body_matrix, 1)
    nup_value <= norb && ndn_value <= norb ||
        throw(ArgumentError("spin-sector electron counts must not exceed orbital dimension"))
    constant = Float64(constant_energy)
    isfinite(constant) ||
        throw(ArgumentError("constant energy must be finite"))
    charges = _cartesian_float_vector(nuclear_charges)
    all(isfinite, charges) ||
        throw(ArgumentError("nuclear charges contain non-finite entries"))
    positions = _cartesian_nuclear_position_matrix(nuclear_positions)
    size(positions, 1) == length(charges) ||
        throw(DimensionMismatch("nuclear position row count must match charges"))
    all(isfinite, positions) ||
        throw(ArgumentError("nuclear positions contain non-finite entries"))
    return _CartesianDensityDensityHamiltonian{Float64}(
        one_body_matrix,
        density_matrix,
        transform,
        nup_value,
        ndn_value,
        constant,
        Tuple(orbital_basis),
        Tuple(density_basis),
        charges,
        positions,
    )
end

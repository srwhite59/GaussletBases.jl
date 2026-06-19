struct _CartesianIDAHamiltonian{T}
    kinetic::Matrix{T}
    nuclear_attraction_unit_by_center::Vector{Matrix{T}}
    electron_electron_ida::Matrix{T}
    nup::Int
    ndn::Int
    nuclear_charges::Vector{T}
    nuclear_positions::Matrix{T}
    nuclear_repulsion::T
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

function _CartesianIDAHamiltonian(
    kinetic,
    nuclear_attraction_unit_by_center,
    electron_electron_ida,
    nup::Integer,
    ndn::Integer,
    nuclear_repulsion::Real;
    nuclear_charges,
    nuclear_positions,
)
    kinetic_matrix =
        _cartesian_check_symmetric_finite_matrix("kinetic energy", kinetic)
    electron_electron_matrix =
        _cartesian_check_symmetric_finite_matrix(
            "IDA electron-electron interaction",
            electron_electron_ida,
        )
    size(electron_electron_matrix) == size(kinetic_matrix) ||
        throw(DimensionMismatch("IDA one-body and electron-electron dimensions differ"))
    center_matrices = Matrix{Float64}[
        _cartesian_check_symmetric_finite_matrix(
            "unit nuclear attraction",
            matrix,
        ) for matrix in nuclear_attraction_unit_by_center
    ]
    all(size(matrix) == size(kinetic_matrix) for matrix in center_matrices) ||
        throw(DimensionMismatch("unit nuclear attraction dimensions differ"))
    nup_value = Int(nup)
    ndn_value = Int(ndn)
    nup_value >= 0 && ndn_value >= 0 && nup_value + ndn_value > 0 ||
        throw(ArgumentError("spin-sector electron counts must be non-negative and nonzero"))
    norb = size(kinetic_matrix, 1)
    nup_value <= norb && ndn_value <= norb ||
        throw(ArgumentError("spin-sector electron counts must not exceed orbital dimension"))
    charges = _cartesian_float_vector(nuclear_charges)
    all(isfinite, charges) ||
        throw(ArgumentError("nuclear charges contain non-finite entries"))
    length(center_matrices) == length(charges) ||
        throw(DimensionMismatch("unit nuclear attraction count must match charges"))
    positions = _cartesian_nuclear_position_matrix(nuclear_positions)
    size(positions, 1) == length(charges) ||
        throw(DimensionMismatch("nuclear position row count must match charges"))
    all(isfinite, positions) ||
        throw(ArgumentError("nuclear positions contain non-finite entries"))
    repulsion = Float64(nuclear_repulsion)
    isfinite(repulsion) ||
        throw(ArgumentError("nuclear repulsion must be finite"))
    return _CartesianIDAHamiltonian{Float64}(
        kinetic_matrix,
        center_matrices,
        electron_electron_matrix,
        nup_value,
        ndn_value,
        charges,
        positions,
        repulsion,
    )
end

function _cartesian_ida_one_body(
    ham::_CartesianIDAHamiltonian;
    charge_multipliers = ham.nuclear_charges,
)
    multipliers = _cartesian_float_vector(charge_multipliers)
    length(multipliers) == length(ham.nuclear_attraction_unit_by_center) ||
        throw(DimensionMismatch("charge multiplier count must match centers"))
    all(isfinite, multipliers) ||
        throw(ArgumentError("charge multipliers contain non-finite entries"))
    matrix = copy(ham.kinetic)
    for (charge, nuclear) in zip(multipliers, ham.nuclear_attraction_unit_by_center)
        matrix .+= charge .* nuclear
    end
    return matrix
end

function _cartesian_ida_nuclear_repulsion(
    ham::_CartesianIDAHamiltonian;
    charge_multipliers = ham.nuclear_charges,
)
    multipliers = _cartesian_float_vector(charge_multipliers)
    length(multipliers) == length(ham.nuclear_charges) ||
        throw(DimensionMismatch("charge multiplier count must match centers"))
    all(isfinite, multipliers) ||
        throw(ArgumentError("charge multipliers contain non-finite entries"))
    repulsion = 0.0
    for right in 2:length(multipliers), left in 1:(right - 1)
        distance = norm(ham.nuclear_positions[right, :] .- ham.nuclear_positions[left, :])
        distance > 0.0 ||
            throw(ArgumentError("nuclear repulsion requires distinct centers"))
        repulsion += multipliers[left] * multipliers[right] / distance
    end
    return repulsion
end

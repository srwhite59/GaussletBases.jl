"""
    CartesianIDAHamiltonian(
        kinetic,
        nuclear_attraction_unit_by_center,
        electron_electron_ida,
        nup,
        ndn;
        nuclear_charges,
        nuclear_positions,
    )

One-basis Cartesian IDA Hamiltonian data. The basis is fixed and localized:
`kinetic`, each uncharged unit-nuclear attraction matrix, and
`electron_electron_ida` all share the same `n x n` basis. Dense operator
matrices are treated as owned/read-only by the Hamiltonian object. Nuclear
positions are stored as `ncenter x 3`, and nuclear repulsion is derived from the
stored physical charges and positions.
"""
struct CartesianIDAHamiltonian{T}
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

_cartesian_float_vector(values::Vector{Float64}) = copy(values)

function _cartesian_float_vector(values)
    return Float64[Float64(value) for value in values]
end

function _cartesian_nuclear_position_matrix(positions::Matrix{Float64})
    size(positions, 2) == 3 ||
        throw(DimensionMismatch("nuclear position matrix must have three columns"))
    return copy(positions)
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

function CartesianIDAHamiltonian(
    kinetic,
    nuclear_attraction_unit_by_center,
    electron_electron_ida,
    nup::Integer,
    ndn::Integer;
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
    center_matrices =
        _cartesian_ida_center_matrices(nuclear_attraction_unit_by_center)
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
    repulsion = _cartesian_ida_nuclear_repulsion(charges, positions)
    return CartesianIDAHamiltonian{Float64}(
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

"""
    one_body_hamiltonian(ham::CartesianIDAHamiltonian; center_weights = ones(ncenter))

Assemble `K + sum_A center_weights[A] * Z_A * U_A` in the Hamiltonian's
localized IDA basis.
"""
function one_body_hamiltonian(
    ham::CartesianIDAHamiltonian;
    center_weights = ones(length(ham.nuclear_charges)),
)
    weights = _cartesian_ida_center_weights(center_weights, ham.nuclear_charges)
    matrix = copy(ham.kinetic)
    for (weight, charge, nuclear) in
        zip(weights, ham.nuclear_charges, ham.nuclear_attraction_unit_by_center)
        matrix .+= weight * charge .* nuclear
    end
    return matrix
end

"""
    nuclear_repulsion(ham::CartesianIDAHamiltonian; center_weights = ones(ncenter))

Return `sum_{A<B} (w_A Z_A)(w_B Z_B)/|R_A - R_B|` using stored physical
charges and positions.
"""
function nuclear_repulsion(
    ham::CartesianIDAHamiltonian;
    center_weights = ones(length(ham.nuclear_charges)),
)
    weights = _cartesian_ida_center_weights(center_weights, ham.nuclear_charges)
    return _cartesian_ida_nuclear_repulsion(
        weights .* ham.nuclear_charges,
        ham.nuclear_positions,
    )
end

function _cartesian_ida_density_matrix(
    ham::CartesianIDAHamiltonian,
    density,
    name::AbstractString,
)
    matrix = _cartesian_dense_float_matrix(density)
    size(matrix) == size(ham.kinetic) ||
        throw(DimensionMismatch("$(name) must have size $(size(ham.kinetic))"))
    all(isfinite, matrix) ||
        throw(ArgumentError("$(name) contains non-finite entries"))
    norm(matrix - transpose(matrix), Inf) <= 1.0e-8 ||
        throw(ArgumentError("$(name) must be symmetric"))
    return 0.5 .* (matrix .+ transpose(matrix))
end

function _cartesian_ida_spin_densities(
    ham::CartesianIDAHamiltonian,
    density_alpha,
    density_beta,
)
    alpha = _cartesian_ida_density_matrix(ham, density_alpha, "density_alpha")
    beta = _cartesian_ida_density_matrix(ham, density_beta, "density_beta")
    return alpha, beta
end

function _cartesian_ida_approximate_interaction_energy(
    ham::CartesianIDAHamiltonian,
    density_alpha,
    density_beta,
)
    alpha, beta = _cartesian_ida_spin_densities(ham, density_alpha, density_beta)
    interaction = ham.electron_electron_ida
    total_row_density = diag(alpha) + diag(beta)
    direct = 0.5 * dot(total_row_density, interaction * total_row_density)
    exchange_alpha = 0.5 * sum(interaction .* alpha .* transpose(alpha))
    exchange_beta = 0.5 * sum(interaction .* beta .* transpose(beta))
    return Float64(direct - exchange_alpha - exchange_beta)
end

function _cartesian_ida_approximate_interaction_fock(
    ham::CartesianIDAHamiltonian,
    same_spin_density::AbstractMatrix{Float64},
    total_row_density::AbstractVector{Float64},
)
    interaction = ham.electron_electron_ida
    fock = Matrix(Diagonal(interaction * total_row_density)) .-
        interaction .* same_spin_density
    return 0.5 .* (fock .+ transpose(fock))
end

function _cartesian_ida_approximate_interaction_fock_alpha(
    ham::CartesianIDAHamiltonian,
    density_alpha,
    density_beta,
)
    alpha, beta = _cartesian_ida_spin_densities(ham, density_alpha, density_beta)
    return _cartesian_ida_approximate_interaction_fock(
        ham,
        alpha,
        diag(alpha) + diag(beta),
    )
end

function _cartesian_ida_approximate_interaction_fock_beta(
    ham::CartesianIDAHamiltonian,
    density_alpha,
    density_beta,
)
    alpha, beta = _cartesian_ida_spin_densities(ham, density_alpha, density_beta)
    return _cartesian_ida_approximate_interaction_fock(
        ham,
        beta,
        diag(alpha) + diag(beta),
    )
end

function _cartesian_ida_approximate_electronic_energy(
    ham::CartesianIDAHamiltonian,
    density_alpha,
    density_beta,
)
    alpha, beta = _cartesian_ida_spin_densities(ham, density_alpha, density_beta)
    one_body = one_body_hamiltonian(ham)
    return Float64(
        sum(one_body .* (alpha + beta)) +
        _cartesian_ida_approximate_interaction_energy(ham, alpha, beta)
    )
end

function _cartesian_ida_approximate_electronic_fock_alpha(
    ham::CartesianIDAHamiltonian,
    density_alpha,
    density_beta,
)
    return one_body_hamiltonian(ham) +
        _cartesian_ida_approximate_interaction_fock_alpha(
            ham,
            density_alpha,
            density_beta,
        )
end

function _cartesian_ida_approximate_electronic_fock_beta(
    ham::CartesianIDAHamiltonian,
    density_alpha,
    density_beta,
)
    return one_body_hamiltonian(ham) +
        _cartesian_ida_approximate_interaction_fock_beta(
            ham,
            density_alpha,
            density_beta,
        )
end

function _cartesian_ida_center_weights(center_weights, charges)
    weights = _cartesian_float_vector(center_weights)
    length(weights) == length(charges) ||
        throw(DimensionMismatch("center weight count must match centers"))
    all(isfinite, weights) ||
        throw(ArgumentError("center weights contain non-finite entries"))
    return weights
end

function _cartesian_ida_nuclear_repulsion(charges, positions)
    repulsion = 0.0
    for right in 2:length(charges), left in 1:(right - 1)
        distance = norm(positions[right, :] .- positions[left, :])
        distance > 0.0 ||
            throw(ArgumentError("nuclear repulsion requires distinct centers"))
        repulsion += charges[left] * charges[right] / distance
    end
    return repulsion
end

function _cartesian_ida_center_matrices(matrices)
    return Matrix{Float64}[
        _cartesian_check_symmetric_finite_matrix(
            "unit nuclear attraction",
            matrix,
        ) for matrix in matrices
    ]
end

function _cartesian_ida_center_matrices(tensor::AbstractArray{<:Real,3})
    return Matrix{Float64}[
        _cartesian_check_symmetric_finite_matrix(
            "unit nuclear attraction",
            view(tensor, :, :, center),
        ) for center in axes(tensor, 3)
    ]
end

function _cartesian_ida_center_tensor(ham::CartesianIDAHamiltonian)
    norb = size(ham.kinetic, 1)
    center_count = length(ham.nuclear_attraction_unit_by_center)
    tensor = Array{Float64}(undef, norb, norb, center_count)
    for center in 1:center_count
        tensor[:, :, center] .= ham.nuclear_attraction_unit_by_center[center]
    end
    return tensor
end

"""
    write_cartesian_ida_hamiltonian(path, ham::CartesianIDAHamiltonian)

Write the minimal versioned Cartesian IDA Hamiltonian JLD2 artifact.
"""
function write_cartesian_ida_hamiltonian(path, ham::CartesianIDAHamiltonian)
    jldopen(String(path), "w") do file
        file["artifact_kind"] = :cartesian_ida_hamiltonian
        file["format_version"] = 1
        file["kinetic"] = ham.kinetic
        file["nuclear_attraction_unit_by_center"] =
            _cartesian_ida_center_tensor(ham)
        file["electron_electron_ida"] = ham.electron_electron_ida
        file["nuclear_charges"] = ham.nuclear_charges
        file["nuclear_positions"] = ham.nuclear_positions
        file["nup"] = ham.nup
        file["ndn"] = ham.ndn
    end
    return path
end

"""
    read_cartesian_ida_hamiltonian(path)

Read a minimal Cartesian IDA Hamiltonian JLD2 artifact.
"""
function read_cartesian_ida_hamiltonian(path)
    return jldopen(String(path), "r") do file
        file["artifact_kind"] === :cartesian_ida_hamiltonian ||
            throw(ArgumentError("artifact_kind must be :cartesian_ida_hamiltonian"))
        Int(file["format_version"]) == 1 ||
            throw(ArgumentError("unsupported Cartesian IDA Hamiltonian format version"))
        CartesianIDAHamiltonian(
            file["kinetic"],
            file["nuclear_attraction_unit_by_center"],
            file["electron_electron_ida"],
            Int(file["nup"]),
            Int(file["ndn"]);
            nuclear_charges = file["nuclear_charges"],
            nuclear_positions = file["nuclear_positions"],
        )
    end
end

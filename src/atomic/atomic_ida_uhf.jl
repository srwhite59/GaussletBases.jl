"""
    density_matrix(Cocc)

Build the spatial one-particle density matrix from occupied orbital
coefficients in the current orthonormal atomic IDA orbital basis.

If `Cocc` is a vector, it is interpreted as one occupied orbital. If it is a
matrix, its columns are interpreted as occupied orbitals.
"""
function density_matrix(coefficients::AbstractVector{<:Real})
    orbitals = Vector{Float64}(coefficients)
    density = orbitals * transpose(orbitals)
    return 0.5 .* (density .+ transpose(density))
end

function density_matrix(coefficients::AbstractMatrix{<:Real})
    occupied = Matrix{Float64}(coefficients)
    density = occupied * transpose(occupied)
    return 0.5 .* (density .+ transpose(density))
end

function _validate_spin_density(ops::AtomicIDAOperators, density::AbstractMatrix{<:Real}, name::AbstractString)
    size(density) == size(ops.one_body.hamiltonian) ||
        throw(DimensionMismatch("$(name) must have size $(size(ops.one_body.hamiltonian))"))
    return Matrix{Float64}(density)
end

function _validate_occupation_count(norbitals::Int, nocc::Integer, name::AbstractString)
    nocc >= 0 || throw(ArgumentError("$(name) must be nonnegative"))
    nocc <= norbitals || throw(ArgumentError("$(name) cannot exceed the orbital dimension"))
    return Int(nocc)
end

function _occupied_columns(coefficients::AbstractMatrix{<:Real}, nocc::Int)
    return nocc == 0 ? coefficients[:, 1:0] : coefficients[:, 1:nocc]
end

function _diagonalize_fock(fock::AbstractMatrix{<:Real}, nocc::Int)
    eig = eigen(Hermitian(Matrix{Float64}(fock)))
    occupied = _occupied_columns(eig.vectors, nocc)
    density = density_matrix(occupied)
    return (
        fock = Matrix{Float64}(fock),
        orbital_energies = Vector{Float64}(real(eig.values)),
        coefficients = Matrix{Float64}(eig.vectors),
        occupied_coefficients = Matrix{Float64}(occupied),
        density = density,
    )
end

"""
    uhf_energy(
        ops::AtomicIDAOperators,
        density_alpha::AbstractMatrix,
        density_beta::AbstractMatrix,
    )

Build the UHF total energy for the present atomic IDA approximation.

The densities are spatial one-particle density matrices in the channel-major
orbital basis. The energy is evaluated as

```text
E = tr[h (rho_alpha + rho_beta)]
    + 1/2 tr[J[rho_alpha + rho_beta] (rho_alpha + rho_beta)]
    - 1/2 tr[K[rho_alpha] rho_alpha]
    - 1/2 tr[K[rho_beta] rho_beta]
```

This is still the current atomic IDA mean-field model. It is not a fully
general four-index Coulomb energy.
"""
function uhf_energy(
    ops::AtomicIDAOperators,
    density_alpha::AbstractMatrix{<:Real},
    density_beta::AbstractMatrix{<:Real},
)
    density_alpha_matrix = _validate_spin_density(ops, density_alpha, "density_alpha")
    density_beta_matrix = _validate_spin_density(ops, density_beta, "density_beta")

    total_density = density_alpha_matrix + density_beta_matrix
    one_body = ops.one_body.hamiltonian
    direct = direct_matrix(ops, total_density)
    exchange_alpha = exchange_matrix(ops, density_alpha_matrix)
    exchange_beta = exchange_matrix(ops, density_beta_matrix)

    return sum(one_body .* total_density) +
           0.5 * sum(direct .* total_density) -
           0.5 * sum(exchange_alpha .* density_alpha_matrix) -
           0.5 * sum(exchange_beta .* density_beta_matrix)
end

"""
    uhf_step(
        ops::AtomicIDAOperators,
        density_alpha::AbstractMatrix,
        density_beta::AbstractMatrix;
        nalpha::Int,
        nbeta::Int,
    )

Take one undamped UHF fixed-point update in the current atomic IDA model.

The input densities are interpreted as `rho_alpha` and `rho_beta` in the
current orthonormal channel-major orbital basis. The function builds
`F_alpha`, `F_beta`, diagonalizes them, occupies the lowest `nalpha` and
`nbeta` orbitals, and returns the updated densities together with the Fock
data and a simple density residual.
"""
function uhf_step(
    ops::AtomicIDAOperators,
    density_alpha::AbstractMatrix{<:Real},
    density_beta::AbstractMatrix{<:Real};
    nalpha::Int,
    nbeta::Int,
)
    density_alpha_matrix = _validate_spin_density(ops, density_alpha, "density_alpha")
    density_beta_matrix = _validate_spin_density(ops, density_beta, "density_beta")
    norbitals = size(ops.one_body.hamiltonian, 1)
    nalpha = _validate_occupation_count(norbitals, nalpha, "nalpha")
    nbeta = _validate_occupation_count(norbitals, nbeta, "nbeta")

    fock_alpha = fock_matrix_alpha(ops, density_alpha_matrix, density_beta_matrix)
    fock_beta = fock_matrix_beta(ops, density_alpha_matrix, density_beta_matrix)
    alpha_solution = _diagonalize_fock(fock_alpha, nalpha)
    beta_solution = _diagonalize_fock(fock_beta, nbeta)

    residual_alpha = norm(alpha_solution.density - density_alpha_matrix, Inf)
    residual_beta = norm(beta_solution.density - density_beta_matrix, Inf)
    residual = max(residual_alpha, residual_beta)
    energy = uhf_energy(ops, alpha_solution.density, beta_solution.density)

    return (
        density_alpha = alpha_solution.density,
        density_beta = beta_solution.density,
        fock_alpha = alpha_solution.fock,
        fock_beta = beta_solution.fock,
        coefficients_alpha = alpha_solution.coefficients,
        coefficients_beta = beta_solution.coefficients,
        occupied_coefficients_alpha = alpha_solution.occupied_coefficients,
        occupied_coefficients_beta = beta_solution.occupied_coefficients,
        orbital_energies_alpha = alpha_solution.orbital_energies,
        orbital_energies_beta = beta_solution.orbital_energies,
        energy = energy,
        residual_alpha = residual_alpha,
        residual_beta = residual_beta,
        residual = residual,
    )
end

function _initial_uhf_density(ops::AtomicIDAOperators, nocc::Int)
    norbitals = size(ops.one_body.hamiltonian, 1)
    nocc = _validate_occupation_count(norbitals, nocc, "occupation count")
    eig = eigen(Hermitian(ops.one_body.hamiltonian))
    occupied = _occupied_columns(eig.vectors, nocc)
    return density_matrix(occupied)
end

"""
    uhf_scf(
        ops::AtomicIDAOperators;
        nalpha::Int,
        nbeta::Int,
        maxiter::Int = 50,
        damping::Real = 0.25,
        tol::Real = 1e-8,
    )

Run a minimal damped UHF fixed-point iteration for the current atomic IDA
model.

This is intentionally a small kernel:

- initialize from the one-body atomic Hamiltonian
- build `F_alpha` and `F_beta`
- occupy the lowest orbitals
- rebuild the densities
- mix with simple damping

No DIIS, no large SCF framework, and no occupation management beyond fixed
`nalpha` and `nbeta`.

The damping convention is:

```text
rho_next = (1 - d) rho_new + d rho_old
```

so `damping = 0` means a full fixed-point update.
"""
function uhf_scf(
    ops::AtomicIDAOperators;
    nalpha::Int,
    nbeta::Int,
    maxiter::Int = 50,
    damping::Real = 0.25,
    tol::Real = 1.0e-8,
)
    norbitals = size(ops.one_body.hamiltonian, 1)
    nalpha = _validate_occupation_count(norbitals, nalpha, "nalpha")
    nbeta = _validate_occupation_count(norbitals, nbeta, "nbeta")
    maxiter >= 1 || throw(ArgumentError("uhf_scf requires maxiter >= 1"))
    0 <= damping < 1 || throw(ArgumentError("uhf_scf requires 0 <= damping < 1"))
    tol > 0 || throw(ArgumentError("uhf_scf requires tol > 0"))

    density_alpha = _initial_uhf_density(ops, nalpha)
    density_beta = _initial_uhf_density(ops, nbeta)
    energies = Float64[]
    residuals = Float64[]
    energy_change = Inf
    converged = false
    iterations = 0

    for iteration in 1:maxiter
        iterations = iteration
        step = uhf_step(ops, density_alpha, density_beta; nalpha = nalpha, nbeta = nbeta)
        mixed_alpha = (1 - damping) .* step.density_alpha .+ damping .* density_alpha
        mixed_beta = (1 - damping) .* step.density_beta .+ damping .* density_beta
        mixed_alpha = 0.5 .* (mixed_alpha .+ transpose(mixed_alpha))
        mixed_beta = 0.5 .* (mixed_beta .+ transpose(mixed_beta))

        residual = max(
            norm(mixed_alpha - density_alpha, Inf),
            norm(mixed_beta - density_beta, Inf),
        )
        energy = uhf_energy(ops, mixed_alpha, mixed_beta)
        energy_change = isempty(energies) ? Inf : abs(energy - energies[end])
        push!(energies, energy)
        push!(residuals, residual)

        density_alpha = mixed_alpha
        density_beta = mixed_beta

        if residual <= tol && energy_change <= tol
            converged = true
            break
        end
    end

    final_alpha = fock_matrix_alpha(ops, density_alpha, density_beta)
    final_beta = fock_matrix_beta(ops, density_alpha, density_beta)
    alpha_solution = _diagonalize_fock(final_alpha, nalpha)
    beta_solution = _diagonalize_fock(final_beta, nbeta)
    final_energy = uhf_energy(ops, density_alpha, density_beta)

    return (
        density_alpha = density_alpha,
        density_beta = density_beta,
        fock_alpha = final_alpha,
        fock_beta = final_beta,
        coefficients_alpha = alpha_solution.coefficients,
        coefficients_beta = beta_solution.coefficients,
        occupied_coefficients_alpha = alpha_solution.occupied_coefficients,
        occupied_coefficients_beta = beta_solution.occupied_coefficients,
        orbital_energies_alpha = alpha_solution.orbital_energies,
        orbital_energies_beta = beta_solution.orbital_energies,
        energy = final_energy,
        energies = energies,
        residuals = residuals,
        energy_change = energy_change,
        iterations = iterations,
        converged = converged,
    )
end

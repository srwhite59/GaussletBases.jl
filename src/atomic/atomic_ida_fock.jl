"""
    fock_matrix(ops::AtomicIDAOperators, density::AbstractMatrix)

Build the narrow spinless / same-spin Fock-style effective one-body matrix for
the current atomic IDA approximation:

```math
F = h + J - K
```

where:

- `h` is the existing atomic one-body Hamiltonian
- `J` is `direct_matrix(ops, density)`
- `K` is `exchange_matrix(ops, density)`

This helper treats the supplied density as both the Hartree density and the
same-spin exchange density. It is therefore useful as a narrow algebraic or
spinless-model helper, but it does **not** by itself fix a physically complete
HF spin convention. A spin-aware mean-field step should instead use
`fock_matrix_alpha` and `fock_matrix_beta`.
"""
function fock_matrix(ops::AtomicIDAOperators, density::AbstractMatrix{<:Real})
    size(density) == size(ops.one_body.hamiltonian) ||
        throw(DimensionMismatch("density matrix must have size $(size(ops.one_body.hamiltonian))"))

    one_body = ops.one_body.hamiltonian
    direct = direct_matrix(ops, density)
    exchange = exchange_matrix(ops, density)
    fock = one_body + direct - exchange
    return 0.5 .* (fock .+ transpose(fock))
end

"""
    fock_matrix_alpha(ops::AtomicIDAOperators, density_alpha::AbstractMatrix, density_beta::AbstractMatrix)
    fock_matrix_beta(ops::AtomicIDAOperators, density_alpha::AbstractMatrix, density_beta::AbstractMatrix)

Build the UHF-style spin-aware Fock matrices for the current atomic IDA
approximation.

For the present channel-major spatial-orbital basis, the first spin-aware
convention is:

```math
F^\alpha = h + J[\rho^\alpha + \rho^\beta] - K[\rho^\alpha]
```

```math
F^\beta = h + J[\rho^\alpha + \rho^\beta] - K[\rho^\beta]
```

So:

- the direct/Hartree term depends on the total density
- the exchange term depends only on the same-spin density

This is still a small assembly layer. It does not manage occupations, mixing,
or SCF iterations.
"""
function fock_matrix_alpha(
    ops::AtomicIDAOperators,
    density_alpha::AbstractMatrix{<:Real},
    density_beta::AbstractMatrix{<:Real},
)
    size(density_alpha) == size(ops.one_body.hamiltonian) ||
        throw(DimensionMismatch("density_alpha must have size $(size(ops.one_body.hamiltonian))"))
    size(density_beta) == size(ops.one_body.hamiltonian) ||
        throw(DimensionMismatch("density_beta must have size $(size(ops.one_body.hamiltonian))"))

    total_density = Matrix{Float64}(density_alpha) + Matrix{Float64}(density_beta)
    fock = ops.one_body.hamiltonian + direct_matrix(ops, total_density) - exchange_matrix(ops, density_alpha)
    return 0.5 .* (fock .+ transpose(fock))
end

"""
    fock_matrix_beta(ops::AtomicIDAOperators, density_alpha::AbstractMatrix, density_beta::AbstractMatrix)

Build the beta-spin UHF-style Fock matrix for the current atomic IDA
approximation:

```math
F^\beta = h + J[\rho^\alpha + \rho^\beta] - K[\rho^\beta]
```

Use this together with `fock_matrix_alpha` when assembling a spin-aware UHF
step in the current density-density / IDA model.
"""
function fock_matrix_beta(
    ops::AtomicIDAOperators,
    density_alpha::AbstractMatrix{<:Real},
    density_beta::AbstractMatrix{<:Real},
)
    size(density_alpha) == size(ops.one_body.hamiltonian) ||
        throw(DimensionMismatch("density_alpha must have size $(size(ops.one_body.hamiltonian))"))
    size(density_beta) == size(ops.one_body.hamiltonian) ||
        throw(DimensionMismatch("density_beta must have size $(size(ops.one_body.hamiltonian))"))

    total_density = Matrix{Float64}(density_alpha) + Matrix{Float64}(density_beta)
    fock = ops.one_body.hamiltonian + direct_matrix(ops, total_density) - exchange_matrix(ops, density_beta)
    return 0.5 .* (fock .+ transpose(fock))
end

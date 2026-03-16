"""
    fock_matrix(ops::AtomicIDAOperators, density::AbstractMatrix)

Build the narrow Fock-style effective one-body matrix for the current atomic
IDA approximation:

```math
F = h + J - K
```

where:

- `h` is the existing atomic one-body Hamiltonian
- `J` is `direct_matrix(ops, density)`
- `K` is `exchange_matrix(ops, density)`

This is a small helper layer, not a full SCF framework. It does not choose
occupations, update orbitals, or manage iterations.
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

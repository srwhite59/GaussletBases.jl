# Atomic and ordinary workflows

For the workflow pages that explain when to use these entry points, start with:

- [Current atomic branch](../explanations/current_atomic_branch.md)
- [Current ordinary branch](../explanations/current_ordinary_branch.md)
- [Example guide](../howto/example_guide.md)

## Small doctest

```jldoctest atomic_and_ordinary
julia> using GaussletBases

julia> rb = build_basis(RadialBasisSpec(:G10; rmax = 8.0, mapping = AsinhMapping(c = 0.1, s = 0.2), xgaussian_count = 0));

julia> grid = radial_quadrature(rb; accuracy = :high);

julia> radial_ops = atomic_operators(rb, grid; Z = 1.0, lmax = 0);

julia> ida = atomic_ida_operators(radial_ops; lmax = 0);

julia> length(orbitals(ida)) == length(rb)
true
```

## Atomic one-body and IDA layers

```@docs
AtomicOneBodyOperators
atomic_one_body_operators
AtomicIDAOperators
atomic_ida_operators
direct_matrix
exchange_matrix
fock_matrix
fock_matrix_alpha
fock_matrix_beta
density_matrix
uhf_energy
uhf_step
uhf_scf
```

## Ordinary mapped and hybrid line

```@docs
MappedOrdinaryOneBody1D
LegacyAtomicGaussianShell
LegacyAtomicGaussianSupplement
LegacySGaussianData
legacy_atomic_gaussian_supplement
legacy_s_gaussian_data
mapped_ordinary_one_body_operators
mapped_cartesian_hydrogen_energy
ordinary_sho_hamiltonian
ordinary_sho_spectrum
CartesianProductOrbital3D
OrdinaryCartesianIDAOperators
ordinary_cartesian_ida_operators
OrdinaryCartesianOrbital3D
OrdinaryCartesianOperators3D
ordinary_cartesian_product_operators
nested_cartesian_operators
ordinary_cartesian_qiu_white_operators
gto_overlap_matrix
gto_occupancy_matrix
GaussletBases.ordinary_cartesian_1s2_check
ordinary_cartesian_vee_expectation
```

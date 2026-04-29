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
legacy_bond_aligned_diatomic_gaussian_supplement
legacy_bond_aligned_heteronuclear_gaussian_supplement
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

For molecular Gaussian supplements, `max_width` is a primitive-level
core/locality cutoff: primitives wider than the cutoff are removed before
placing the supplement on each nucleus. A contracted shell disappears only when
no primitive remains after this filtering.

## Hydrogenic core corrections

`apply_ordinary_cartesian_corrections` is a post-assembly correction for
`OrdinaryCartesianOperators3D`. The public path uses the projector one-body
reference correction to shift the lowest hydrogenic core eigenvalue to
`-Z^2/2`. ESOI is explicit opt-in through
`HydrogenicCoreProjectorCorrectionSpec(; include_esoi=true)` and calibrates the
same `1s^2` Coulomb scalar to `5Z/8`.

Because these are post-assembly matrix corrections, the returned operators keep
the total matrices authoritative and drop invalid decomposition sidecars:
`kinetic_one_body === nothing`, `nuclear_one_body_by_center === nothing`, and
`nuclear_term_storage == :total_only`. Internal diagnostic modes such as
`:local_exact` are not part of the public correction surface.

For branch or counterpoise calculations, use
`ordinary_cartesian_corrected_branch`. It assembles one branch one-body matrix
with caller-provided `nuclear_charges`, applies one or more
`HydrogenicCoreBranchCorrectionSpec`s sequentially, and returns corrected
matrices rather than a new operator payload. The branch surface default
`orbital_selector = :localized_lowest` chooses a center-local subspace by
nearest carried nucleus and degenerates to the full space for one-center
payloads without carried nuclei. `orbital_selector = :global_lowest` remains an
explicit debug/reference selector for single-correction calls. Multi-correction
calls require localized selectors and reject duplicate localized center indices.
Branch diagnostics report `branch_nuclear_charges`, `correction_count`,
`overlap_error`, optional `corrected_center_indices`, and a tuple of
per-correction diagnostics in `corrections`.

The intended counterpoise pattern is to build one operator payload that retains
per-center nuclear terms, then reuse it for the full branch and each fragment
branch:

```julia
spec_a = HydrogenicCoreBranchCorrectionSpec(; Z = Z, nucleus = nuclei[1])
spec_b = HydrogenicCoreBranchCorrectionSpec(; Z = Z, nucleus = nuclei[2])

full = ordinary_cartesian_corrected_branch(
    operators;
    nuclear_charges = [Z, Z],
    corrections = [spec_a, spec_b],
)
atom_a = ordinary_cartesian_corrected_branch(
    operators;
    nuclear_charges = [Z, 0.0],
    corrections = spec_a,
)
atom_b = ordinary_cartesian_corrected_branch(
    operators;
    nuclear_charges = [0.0, Z],
    corrections = spec_b,
)
```

```@docs
HydrogenicCoreProjectorCorrectionSpec
HydrogenicCoreBranchCorrectionSpec
OrdinaryCartesianCorrectionResult
OrdinaryCartesianBranchCorrectionResult
apply_ordinary_cartesian_corrections
ordinary_cartesian_corrected_branch
```

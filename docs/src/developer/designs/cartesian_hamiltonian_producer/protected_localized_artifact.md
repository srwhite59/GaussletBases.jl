# Protected-Localized Artifact Contract

Status: implemented, internal, and opt-in under
`HP-RG-PROTECT-ART-FN-01` / `HP-RG-PROTECT-ART-TEST-01` and
`HP-RG-PROTECT-ARTLOC-FN-01` / `HP-RG-PROTECT-ARTLOC-TEST-01`.

This page is the canonical persistence and row-locality contract for the
protected-localized Hamiltonian. The numerical construction of the protected
fixed and localized bases, exact `H1_L`, and inherited-site `Vee_L` belongs to
the [protected-localized basis convention](protected_localized_basis.md). This
page records those completed objects; it does not redefine their construction.

## Identity And Ownership

`src/cartesian_ida_hamiltonian.jl` owns the artifact writer and reader:

- `write_protected_localized_ida_hamiltonian`;
- `read_protected_localized_ida_hamiltonian`.

The on-disk identity is fixed:

```text
artifact_kind      = :protected_localized_inherited_site_ida_hamiltonian
format_version     = 1
convention_id      = :protected_localized_inherited_site_ida_v1
convention_version = 1
```

`src/cartesian_residual_gaussians/augmented_operators.jl` owns
`protected_localized_row_locality`, which constructs optional locality
metadata before the writer validates and stores it.

Protected ladder members may consume this artifact contract. Their manifests,
cross-basis transfers, restart sidecars, and summaries are governed separately
by [protected-localized ladder bundles](protected_localized_ladder.md) and are
not part of this schema.

## Native Artifact Contract

The canonical matrices are finite, symmetric, equal-size `H1_L` and `Vee_L`
in native protected-localized order. The artifact also carries:

- `nup`, `ndn`, and `final_dimension`;
- nuclear charges and positions, with derived nuclear repulsion;
- sector counts and basis-specific contiguous ranges;
- diagnostics and source/current provenance;
- basis controls and geometry inputs;
- a validated Coulomb-expansion summary;
- localized-ordering metadata;
- optional validated `row_locality`.

The native sector law is:

```text
final_dimension = base + compact_R
M               = base + compact_R
Z               = protected_Z + broad_Z
Qperp           = final_dimension - Z
```

The stored ranges are basis-specific and remain in their original,
unpermuted coordinate order:

```text
M / inherited site order: base_G, compact_R
F / fixed construction order: protected_Z, broad_Z, Qperp
L / native localized order: 1:final_dimension
```

`H1_L`, `Vee_L`, and row-locality vectors use native `L` row numbering. The
`base_G` / `compact_R` locality labels follow the inherited `M` site order.
The `protected_Z` / `broad_Z` / `Qperp` ranges describe the pre-localization
fixed `F` construction; after multiplication by `W`, they are not contiguous
localized-`L` row sectors and must not be used to classify `L` rows. A
consumer may make a permuted working copy, but no permutation changes or
reinterprets the stored matrices or basis-specific ranges.

Unless a writer supplies reviewed overrides, localized ordering is:

```text
order          = :inherited_M_site_order
transform      = :angular_style_projected_localization
parent_M_order = :base_G_then_compact_R
fixed_F_order  = :protected_Z_then_broad_Z_then_Qperp
interaction    = :inherited_pre_injection_site_order
exact_one_body = :protected_localized_dense_transform
```

## Row Locality

For the protected-localized transform `loc`, the actual native `L`-basis
coefficient matrix is

```text
ML = hcat(loc.C, loc.Qp) * loc.W
```

Given inherited native `M`-space position matrices `X_M`, `Y_M`, and `Z_M`,
the centers are their diagonal expectations in `L`:

```text
center_x = diag(ML' * X_M * ML)
center_y = diag(ML' * Y_M * ML)
center_z = diag(ML' * Z_M * ML)
```

These numerical expectations, not labels or manifest geometry, are center
authority. When matching second moments are supplied, spreads are:

```text
spread_axis = sqrt(max(<axis^2>_L - center_axis^2, 0))
```

The required locality vectors are native-order `center_x`, `center_y`, and
`center_z`, plus `sector_label` and `native_sector_index`. Sector metadata
labels rows `1:base` as `base_G` and the remaining native rows as `compact_R`,
with sector-local indices starting at one.

The permutations have one precise direction:

- `z_order_to_native[k]` is the native row at z-sorted position `k`;
- `native_to_z_order[i]` is the z-sorted position of native row `i`.

`z_order_to_native` sorts by `(center_z, native_index)`, and
`native_to_z_order` is its validated inverse. These vectors are metadata only.
The existing convention never stores z-sorted `H1_L`, `Vee_L`, or sector
ranges as canonical data.

## Validation And Compatibility

The writer validates before writing:

- finite, symmetric, same-size `H1_L` and `Vee_L`;
- nonnegative electron counts that do not exceed `final_dimension`;
- matching nuclear-charge and nuclear-position counts;
- the native sector laws above;
- required provenance, basis-control, and geometry-input fields;
- the Coulomb-expansion summary, plus complete localized-ordering metadata;
- row-locality shape and values when locality is supplied.

The protected reader rejects missing or unrecognized artifact kind, format
version, convention ID, or convention version. It also rejects inconsistent
matrix dimensions, electron or sector counts, required metadata, diagnostics,
missing or non-two-endpoint sector-range records, and malformed locality data.

New writes require a complete validated Coulomb-expansion summary. Legacy
protected-localized artifacts from before that summary was introduced may
omit it; readback then returns `coulomb_expansion = nothing`, and consumers
must not infer `:high` accuracy from the absence.

Legacy protected-localized artifacts may omit `row_locality` entirely. The
`row_locality/center_x` key is the compatibility presence marker: when it is
absent, readback returns `row_locality = nothing`. When it is present, required
vectors must have `final_dimension` entries, centers must be finite,
permutations must cover `1:final_dimension` and be mutual inverses, z order
must be monotone with native-index tie breaking, and sector metadata must
agree with native `base_G` / `compact_R` rows. Spreads are optional, but
`spread_x`, `spread_y`, and `spread_z` must be present together, finite, and
nonnegative; partial spreads are invalid.

The ordinary `read_cartesian_ida_hamiltonian` accepts only
`:cartesian_ida_hamiltonian`. It rejects the protected-localized artifact kind;
consumers must call the protected reader explicitly.

## Implementation And Evidence

Source ownership:

- `src/cartesian_ida_hamiltonian.jl`;
- `src/cartesian_residual_gaussians/augmented_operators.jl`.

Shared committed reader regression:

- `test/ida/cartesian_ida_hamiltonian_runtests.jl`.

There is no dedicated committed protected-artifact test file. The bounded
protected write/readback, rejection, native-center, permutation, sector, and
legacy-no-locality smokes used for implementation acceptance are recorded by:

- commit `fd105b751` and manager running-log Pass 299;
- commit `3fe2af697` and manager running-log Pass 301.

Both artifact persistence and row-locality metadata are implemented. Their
validation IDs describe the accepted gates; they are not future or
not-yet-implemented lifecycle entries.

## Exclusions

This contract does not own or authorize:

- protected-basis numerical construction or alternative interaction rules;
- ladder manifests, transfers, restart sidecars, or bundle summaries, which
  belong to [protected-localized ladder bundles](protected_localized_ladder.md);
- EGOI targets or corrections, which remain governed by
  [retained-GTO local-product EGOI](retained_gto_egoi.md);
- occupied-first or additive-reference construction;
- rho0, screened-Hartree, or corrected artifact variants;
- solver behavior or public-driver/default workflow;
- matrix reordering under the existing convention;
- changes to ordinary Cartesian IDA, PQS, WL, or RG artifact semantics;
- Cr2-specific production claims.

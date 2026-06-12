Pass 006 response: PQS electron-nuclear source-factor audit

Status: completed as a read-only audit. I made no production source or test
changes.

Inspected surfaces:

- `src/cartesian_pair_block_materialization/white_lindsey_electron_nuclear.jl`
- `src/cartesian_pair_block_materialization/pqs_source_safe_terms.jl`
- `src/cartesian_pair_block_materialization/pqs_source_one_body.jl`
- `src/cartesian_nested_faces.jl`
- `src/ordinary_qw_raw_blocks.jl`
- `src/ordinary_qw_operator_assembly.jl`
- `src/CartesianContractedParentMetrics.jl`
- `src/CartesianCPBBlockProviders.jl`
- `src/CartesianParentAxisFactors.jl`
- `src/ordinary_coulomb.jl`
- `docs/src/developer/pqs_source_box_operator_framework.md`
- `docs/src/developer/numerical_contracts.md`

Findings:

1. Existing Gaussian/Coulomb authority

`CoulombGaussianExpansion` in `src/ordinary_coulomb.jl` stores positive
`coefficients` and positive `exponents` for the approximation
`1/r ~= sum_t c_t exp(-zeta_t r^2)`. The object itself carries no nuclear
charge and no attraction sign.

The active WL electron-nuclear path in
`white_lindsey_electron_nuclear.jl` builds an axis context from:

- a parent axis bundle;
- `coulomb_expansion.coefficients`;
- `coulomb_expansion.exponents`;
- a by-center record containing charge and location.

It calls `_mapped_ordinary_gausslet_1d_bundle(...; exponents, center, backend)`
for each axis and reads
`pgdg_intermediate.gaussian_factor_terms`. The local support block contracts

```text
sum_t (-c_t) * Gx_t * Gy_t * Gz_t
```

using `Float64[-value for value in expansion.coefficients]`.

The center charge is recorded in metadata, but the local/global by-center
matrix is uncharged:

- `nuclear_charge_recorded = true`
- `nuclear_charge_applied = false`
- `centers_summed = false`
- `uncharged_by_center_convention = true`

The Hamiltonian assembly convention then applies `+Z * V_center`, where
`V_center` is already the unit-charge negative attraction matrix.

The old nested/QW code agrees with this split in its by-center storage path.
`_qwrg_contracted_nuclear_axis_term_tables(...)` builds centered per-axis
Gaussian term tables once per center/axis, and
`_qwrg_bond_aligned_direct_contracted_nuclear_one_body_by_center(...)` /
`_qwrg_bond_aligned_staged_by_center_nuclear_one_body_by_center(...)` use
`Float64[-value for value in expansion.coefficients]`. Final reassembly in
`_assemble_one_body_hamiltonian(...)` multiplies each stored per-center nuclear
matrix by the requested nuclear charge.

2. PQS source-space status

The active PQS source-mode path in
`pqs_source_safe_terms.jl` currently supports only caller-supplied 1D factor
matrices for:

- overlap;
- position;
- x2;
- kinetic.

It forms raw product-source blocks through `_pqs_source_pair_product_result`.
Retained source-mode contraction is already available through
`pqs_source_pair_retained_one_body_block(source_result, left_rule, right_rule)`,
which applies:

```text
O_retained = O_source[left_retained_columns, right_retained_columns]
```

The selector/batch layer in `pqs_source_one_body.jl` intentionally restricts
retained terms to `:overlap` and `:kinetic` today.

I did not find an active CPBM/PQS helper that materializes source-space
electron-nuclear Gaussian factor matrices for the new retained source-mode
route. There are old private CCPM helpers:

- `_pqs_pqs_source_box_local_gaussian_sum_block(...)`
- `_pqs_pqs_source_box_centered_local_gaussian_sum_block(...)`
- `_pqs_pqs_source_box_nuclear_attraction_by_center(...)`
- related product/PQS variants.

Those are retirement-era private source-box/shadow helpers, not the new CPBM
route authority. They are useful as convention/oracle references only.

3. CCPM convention risk

The CCPM local-Gaussian lane builds a positive Gaussian sum first:

```text
sum_t c_t * Gx_t * Gy_t * Gz_t
```

Then `_source_box_nuclear_attraction_by_center(...)` applies `-Z` and also
returns a derived total block. That convention is documented in
`pqs_source_box_operator_framework.md`.

For the new PQS retained source-mode route, I recommend following the active
post-CPB/WL by-center convention instead:

- build one unit-charge by-center electron-nuclear matrix using `-c_t`;
- record charge and center identity;
- do not multiply by `Z`;
- do not sum centers;
- apply charges only in a later one-electron diagnostic/Hamiltonian assembly
  helper.

This avoids reintroducing the old CCPM physical wrapper as a production
provider contract.

Recommended implementation seam:

Add the new PQS source-space electron-nuclear code in the CPBM PQS source
layer, adjacent to the current safe one-body helpers:

- likely new source file or small section near
  `src/cartesian_pair_block_materialization/pqs_source_safe_terms.jl`;
- exported only if later needed by the module contract, matching the existing
  CPBM export style;
- term name should be aligned with the active WL/CPB lane, e.g.
  `:source_electron_nuclear_by_center` and
  `:retained_source_electron_nuclear_by_center`.

Smallest helper shape:

```text
pqs_source_pair_electron_nuclear_by_center_block(
    record;
    coulomb_expansion,
    center_record or center,
    gaussian_axis_terms_1d = nothing,
    parent/source axis object = ...
)
```

The implementation should:

1. Reuse `_is_ready_pqs_source_pair_record`, `_pqs_source_pair_dims`,
   `_pqs_source_mode_ordering`, `CRPS.source_mode_indices`, and
   `_fill_source_mode_product_block!` style where possible.
2. Build or accept per-axis, term-first Gaussian factor arrays with shape
   `(nterms, left_axis_source_count, right_axis_source_count)`.
3. Use `-expansion.coefficients` for the unit-charge by-center attraction
   matrix.
4. Preserve by-center separation and charge metadata:
   `nuclear_charge_recorded = true`, `nuclear_charge_applied = false`,
   `centers_summed = false`.
5. Return a raw product-source-mode `PairBlockMaterializationResult`.
6. Reuse the existing retained contraction:
   `pqs_source_pair_retained_one_body_block(source_result)`.

Minimal missing helper:

The gap is not the 3D contraction. That is already patterned by
`_pqs_source_pair_product_result` and the old CCPM local-Gaussian helper. The
gap is a module-owned way to produce PQS source-mode 1D Gaussian factor arrays
from the current parent/source axis facts without using CCPM as production
authority.

Concretely, the next implementation needs a narrow 1D helper that, for one
axis and one center coordinate, returns term-first matrices compatible with the
PQS source-mode dimensions/order:

```text
(nterms, n_left_source_axis, n_right_source_axis)
```

It should reuse the same analytic Gaussian factor convention as
`_mapped_ordinary_gausslet_1d_bundle(...).pgdg_intermediate.gaussian_factor_terms`
or the existing `gaussian_factor_matrices(...)` machinery, depending on which
source-axis object the PQS route currently owns. It should not route through
CCPM `_pqs_pqs_source_box_*` helpers except as a test oracle.

Blockers/convention risks:

- Need to identify the actual source-axis object available to the new PQS
  source route. If the route only has dimensions and retained rules, it cannot
  generate centered Gaussian factors yet.
- Do not use CCPM `_source_box_nuclear_attraction_by_center(...)` as the new
  provider surface, because it applies `-Z` and reports a total block.
- Do not claim PGDG provenance unless the helper really consumes
  `pgdg_intermediate.gaussian_factor_terms`.
- Do not introduce source/retained weights or IDA language for this one-body
  term.
- Keep sign/charge language explicit: the recommended new block is negative
  unit-charge by-center data, not a positive `gaussian_sum` and not a
  charge-applied Hamiltonian contribution.

First physics diagnostic after implementation:

Use a one-center H or He+ H1 diagnostic, not RHF:

- build retained source-mode overlap;
- build retained source-mode kinetic;
- build one retained source-mode electron-nuclear by-center matrix;
- assemble `H1 = K + Z * V_center` only in the diagnostic assembly layer;
- report the lowest H1 eigenvalue and variational distance from the
  hydrogenic reference.

For H, the reference is `-0.5` Hartree. For He+, the one-electron reference is
`-2.0` Hartree. This catches the sign, charge stage, factor convention, and
retained contraction before any IDA/RHF interpretation.

Nonclaims to preserve:

- no shell realization;
- no Lowdin cleanup;
- no IDA or density-density;
- no Hamiltonian/driver adoption;
- no full-parent CPB or direct Cartesian fallback;
- no QW/fixed-block route adoption;
- no public/default route or artifacts;
- no center summation inside the provider;
- no nuclear charge application inside the by-center provider block.

Deletion/shrinkage report:

- Nothing was deleted or simplified in this pass because the blurb required a
  read-only audit and the existing code surfaces are still needed as either
  active WL/CPB code or historical CCPM references.
- No new test or probe artifact was added.
- Once the CPBM PQS retained source-mode nuclear helper exists, the likely
  retirement candidates are the CCPM `_pqs_pqs_source_box_*local_gaussian*`
  and `_pqs_pqs_source_box_*nuclear_attraction*` test pressure/summaries that
  still preserve the old positive-Gaussian-plus-physical-wrapper vocabulary.
  They should remain oracle-only until the new CPBM helper has a direct
  equivalence check.

Validation:

- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.

-- repo-doer@macmini

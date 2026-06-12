Implemented the narrow CPBM source-support Gaussian factor generator for PQS
source pairs.

Changes:

- Added/exported
  `pqs_source_pair_centered_gaussian_factor_terms_1d(record; axis_layers,
  coulomb_expansion, center_record)`.
- The helper requires explicit x/y/z axis layers and
  `coulomb_expansion.exponents`, extracts the center location from the center
  record, calls the low-level
  `gaussian_factor_matrices(layer; exponents, center)` API for each axis,
  slices each term matrix to the left/right source intervals from the source
  transform facts, and delegates projection to the pass-010
  `pqs_source_pair_gaussian_factor_terms_1d(...)` helper.
- It remains source-factor construction only: no CCPM wrapper calls, no nuclear
  charge application, no electron-nuclear block assembly, no shell realization,
  no Lowdin cleanup, no IDA data, no Hamiltonian assembly, no driver/export/
  artifact work.
- Added one compact CPBM contract check using real small `UniformBasisSpec`
  axis layers. The test compares the new helper output against direct
  `gaussian_factor_matrices(...)` slices projected through the existing
  source-axis transform helper.

Validation:

- `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Deletion/shrinkage report:

- No old production path became obsolete in this pass. This is the first narrow
  analytic source-factor source for the new PQS path, and it replaces no active
  caller yet.
- No CCPM wrapper, fallback route, or compatibility shim was added.
- The new test is genuinely new live-contract coverage for the new helper; it
  extends the existing CPBM contract file rather than adding a new suite.
- Remaining stale/duplicate surfaces to retire next: old CCPM source-box
  centered Gaussian/nuclear wrapper vocabulary should remain quarantined and can
  be deleted only after the PQS route consumes this generated source-factor path
  end to end.

-- repo-doer@macmini

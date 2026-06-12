Implemented the repo-owned PQS source-axis transform fact builder in CPBM.

Changes:

- Added/exported
  `pqs_source_axis_transform_facts_from_pgdg_axes(axis_sources;
  source_intervals, source_mode_dims, enforce_symmetric_odd = false)`.
- Added `src/cartesian_pair_block_materialization/pqs_source_axis_transforms.jl`.
- The helper accepts x/y/z PGDG intermediates or mapped ordinary gausslet
  bundles, calls the low-level `_nested_doside_1d(...)` kernel directly, and
  returns materialized `CartesianRawProductSources.AxisSourceTransformFact`s.
- It records compact diagnostics:
  requested/resolved source dims, coefficient shapes, per-axis overlap errors,
  max overlap error, source mode adjustment flag, transform source, and explicit
  no shell/Lowdin/IDA/Hamiltonian/driver/artifact claims.
- It does not call `_cartesian_raw_product_box_plan(...)`, old fixed-block
  matrices, CCPM wrappers, shell realization, Lowdin cleanup, IDA,
  density-density, RHF, global drivers, exports, or artifacts.

Focused test:

- Extended the existing CPBM contract file with one compact testset.
- The test uses the pass-015 13-point / source interval `2:12` / source dims
  `(5, 5, 5)` fixture.
- It compares the new helper's per-axis coefficient matrices against the old
  `_cartesian_source_box_axis_transform_plan(...)` result as an oracle.
- It verifies the facts are materialized, shaped `(11, 5)`, match the oracle
  local coefficients, have small overlap error, and carry the intended nonclaim
  flags.

Validation:

- `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Deletion/shrinkage report:

- No old code was deleted in this pass because the old nested transform helper
  is still used as an oracle in the new focused contract check.
- The new helper removes the need for future PQS source-box probes to call the
  old `_cartesian_raw_product_box_plan(...)` wrapper just to obtain source-axis
  transform facts.
- The new test is genuinely new live-contract coverage for the repo-owned
  transform-builder seam and replaces the pass-015 probe's ad hoc dependence on
  old raw-plan output for future validation.
- Remaining stale/duplicate surfaces to retire next: once downstream PQS H1
  probes consume `pqs_source_axis_transform_facts_from_pgdg_axes(...)`, the old
  nested source-box transform wrappers should be quarantined as oracle-only and
  any tests that present them as the route source should be shrunk or deleted.

-- repo-doer@macmini

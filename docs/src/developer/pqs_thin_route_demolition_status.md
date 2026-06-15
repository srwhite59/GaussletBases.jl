# PQS Thin Route Demolition Status

This note tracks the non-adabatic thinning branch
`demolition/pqs-thin-route`. The branch is allowed to break package load and
tests while report/status/test scaffolding is removed. Do not treat this as a
merge-ready state until the breakages below are reviewed and repaired or
accepted.

## Checkpoint 1 - Report Scaffold Tests

Commit:

- `5596b15c` - `Demolish Cartesian/PQS report scaffold tests`

Deleted:

- exact CRC/terminal print-line tests;
- flat terminal route glue tests;
- staged low-order policy tests that mainly asserted old status vocabulary;
- driver module-boundary/report tests that protected private field clouds;
- direct helper code for CRC print-line and terminal-shellification report
  aliases.

Validation before commit:

- `git diff --check` passed.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.

Line-count impact:

- 4,832 deleted lines.

## Checkpoint 2 - Route-Global/Metadata Demolition Cut

Status:

- uncommitted demolition cut for review.

Deleted:

- route-core sidecar adapter include and its sidecar tests;
- terminal retained-unit and pair-stage fingerprint tests;
- policy-metadata/report smoke tests;
- legacy contracted-parent source-box fixture include and file;
- old route-global one-body adapter files;
- White-Lindsey adapter-summary file;
- direct tests for the removed route-global/adapter-summary layer.

Intended breakage:

- old route-global WL adapter APIs are no longer defined by
  `CartesianPairBlockMaterialization`;
- old raw-box/fixture producer calls from diatomic diagnostic scaffolding may
  now point at deleted legacy fixture helpers;
- tests that assert route-global status names, adapter summaries, source-backed
  fixture metadata, or print/report mirrors are expected to fail.

Scaffolding class:

- report/status aliases;
- route-core sidecar mirrors;
- flat route-global one-body adapter vocabulary;
- legacy source-box fixture/oracle helpers;
- test-only fingerprint and metadata inventories.

Validation:

- `git diff --check` passed.
- Caller grep still finds references to removed scaffold names in docs, old
  tests, and several source seams that should be reviewed before repair:
  `src/pqs_source_box_route_driver_helpers.jl` still references route-core
  sidecar inventory helpers, `src/pqs_source_box_diatomic_complete_core_shell.jl`
  still references `_pqs_pqs_product_raw_box_route_producer(...)`, and
  `white_lindsey_unit_coefficients.jl` / old tests still reference the removed
  White-Lindsey adapter-summary helpers.
- Package load was tried once and intentionally not repaired. First error:
  `_white_lindsey_unit_metadata_value` is undefined during the
  `precompile_workloads.jl` decomposed WL one-body context. This is expected
  fallout from removing `white_lindsey_adapter_summary.jl`, not evidence that a
  PQS numerical kernel was edited.

Line-count impact:

- 12,561 deleted lines in the uncommitted second cut.

Current repair strategy:

- Do not rebuild the removed scaffold.
- First classify package-load failures as either undefined old adapter names or
  real numerical-kernel loss.
- If a live route still needs one removed numerical primitive, restore only that
  primitive behind a compact route-owned object, not the old report/status
  field cloud.

## Checkpoint 3 - Precompile and Route-Global Caller Cut

Status:

- uncommitted demolition cut for review.

Deleted:

- `precompile_workloads.jl` and its include from `GaussletBases.jl`;
- White-Lindsey boundary-stratum coefficient adapter files;
- decomposed White-Lindsey route-global one-body and density-density files;
- route-global combined GTO layout, mixed-block, matrix assembly,
  final-basis, density-density, and atom-GTO route wrapper files;
- direct tests for those old route-global, one-body selector, GTO readiness,
  and WL acceptance surfaces.

Simplified:

- Removed White-Lindsey selector entries from the mixed one-body dispatch
  surface after deleting the implementation they pointed to.
- Removed route-core sidecar calls from the driver helper instead of replacing
  the deleted sidecar adapters.

Validation:

- `git diff --check` passed.
- Caller grep still finds only docs/old logs for deleted route-global/GTO names,
  plus one active source diagnostic seam:
  `_pqs_source_box_route_driver_diatomic_raw_box_route_payload(...)` still calls
  `_pqs_pqs_product_raw_box_route_producer(...)`.
- Package load was attempted once after the cut and now passes:
  `julia --project=. -e 'using GaussletBases; println("load ok")'`.

Line-count impact:

- 14,568 deleted lines in the uncommitted third cut.

Remaining live-looking breakage:

- The diatomic raw-box route payload is still a live source caller of the
  deleted legacy contracted-parent fixture producer. It did not block package
  load, but it is the next demolition/repair decision point.
- Lower-level WL boundary-stratum route vocabulary still exists in retained
  units, pair-operator planning, and terminal lowering. This pass removed the
  old route-global implementation, not the lower-level route vocabulary.

Recommended next cut:

- Audit and either delete or replace
  `_pqs_source_box_route_driver_diatomic_raw_box_route_payload(...)` and its
  downstream status/report fields. Do not restore
  `_pqs_pqs_product_raw_box_route_producer(...)`.
- After that, re-audit remaining WL boundary-stratum vocabulary to separate
  lower-level route semantics from deleted old route-global implementation
  scaffolding.

## Checkpoint 4 - Driver Scaffold Test Cut

Status:

- uncommitted demolition cut for review.
- Important correction: an initial over-aggressive attempt deleted
  `bin/cartesian_ham_builder.jl`, driver inputs, and mixed route-construction
  source files. That was rolled back before this checkpoint. The driver entry
  point, input-driven workflow, and route-construction files are protected for
  now. The demolition target is the bulky private report/status scaffolding
  around them, not the existence of a driver.

Deleted:

- direct driver/report tests for fake-PQS, He PQS, one-center config smoke,
  atom-growth materialized checkpoints, parent-contract field clouds, and
  standard source-box route setup;
- runner includes that anchored two deleted scaffold tests.

Simplified:

- None of the protected driver/source route files remain deleted in the current
  diff.
- The remaining cut removes old tests that primarily assert private driver
  field clouds, old setup surfaces, or report-driven endpoint fixtures.

Validation:

- `git diff --check` passed after restoring the protected driver/source files.
- Package load was attempted once and passes:
  `julia --project=. -e 'using GaussletBases; println("load ok")'`.

Line-count impact:

- Current corrected fourth-cut diff is 1,347 deleted test/runner lines plus this
  status-note update.
- Current total branch pressure relative to `main` after correcting the
  overreach: 214 insertions and 33,308 deletions.

Current breakage assessment:

- The driver and route source files are restored and should still be treated as
  live until their construction/math can be separated from report/status
  bureaucracy.
- The diatomic raw-box diagnostic seam remains the next source-level decision:
  cut report/status callers without deleting the core driver or identifiable
  construction path.

Recommended next cut:

- Split decisions inside mixed files:
  keep or extract construction/math, and cut report/status/alias/export
  scaffolding.
- Start with `_pqs_source_box_route_driver_diatomic_raw_box_route_payload(...)`
  and downstream fields, but do not delete `bin/cartesian_ham_builder.jl`, all
  driver inputs, or the independent H2 PQS construction path without explicit
  approval.

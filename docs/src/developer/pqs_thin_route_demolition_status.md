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

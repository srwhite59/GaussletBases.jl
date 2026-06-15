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

## Checkpoint 5 - Route-Status Test and Artifact-Key Cut

Status:

- uncommitted demolition cut for review.
- The protected driver entry point and driver inputs remain present.

Deleted:

- synthetic route-core, pair-operator-plan, pair-block-materialization, and
  terminal-lowering contract tests that mostly asserted status/materialized
  sidecars and private route-core summaries;
- old shellification module/plan integration tests that preserved
  metadata-only route summaries, materialization flags, and private
  low-order planning vocabulary;
- old atom-growth, high-order recipe audit, fixed-block timing, QW adapter, and
  residual/compression policy tests;
- the endpoint manifest mirror.

Simplified:

- Removed the flattened He/H2 artifact-key writer cloud from
  `pqs_source_box_route_driver_reporting.jl`.
- Kept the basic driver artifact path: `cartesian_save` still writes durable
  `report` and optional `materialization` payloads to JLD2.
- Shrunk printed/TSV materialization fields to a compact status/dimension
  summary instead of the old low-order/ham/artifact status field list.

Validation:

- `git diff --check` passed.
- Package load was attempted once and passes:
  `julia --project=. -e 'using GaussletBases; println("load ok")'`.
- Caller grep for the deleted artifact writers and deleted test names finds no
  active source/test/bin callers. Remaining hits are historical docs/logs.
- `cartesian_selected_terminal_lowering_contract_inventory` still has active
  source/driver callers, so only its stale test was deleted in this cut.

Line-count impact:

- 9,148 deleted lines in the uncommitted fifth cut.
- Current total branch pressure relative to `main`: 218 insertions and 42,456
  deletions.

Current breakage assessment:

- Package load remains green.
- Old helper/status tests are intentionally gone.
- No driver input file or driver entry point was deleted.
- The flattened artifact contract is intentionally broken; future survival
  checks should read compact semantic facts from the durable report or a smaller
  driver-owned artifact, not exact legacy key lists.

Recommended next cut:

- Audit the remaining source-side selected-terminal-lowering inventory and
  pair/materialization summary fields for active driver construction use versus
  report/status-only use.
- Keep the driver alive, but continue reducing `pqs_source_box_route_driver_*`
  report/status field clouds.

## Checkpoint 6 - Low-Order Materializer Collapse

Status:

- uncommitted demolition cut for review.
- The protected driver entry point and driver inputs remain present.

Deleted/simplified:

- Replaced the 2,655-line route-configured low-order/WL materialization
  implementation with a compact blocked materialization stage.
- Removed the old route-configured one-center/diatomic materializer probes,
  basis adapters, Ham adapters, WL preflight helpers, and artifact export
  decision tree.
- Simplified driver printing so it no longer expects the removed low-order,
  basis-artifact, and Ham-artifact status field cloud.
- Replaced deleted CCPM standard setup and parent-axis readiness/count helpers
  with small local driver records. This keeps the driver workflow alive without
  restoring the old CCPM setup scaffolding.
- Bypassed low-order terminal-route unit summary bookkeeping for source-box
  routes.

Validation:

- `git diff --check` passed.
- Package load passed:
  `julia --project=. -e 'using GaussletBases; println("load ok")'`.
- A protected driver readiness smoke completed with saving disabled:
  `h2_pqs_q5_independent_source_box_r4.jl save_artifact=false save_tsv=false`.

Line-count impact:

- 2,713 deleted lines in the uncommitted sixth cut.
- Current total branch pressure relative to `main`: 438 insertions and 45,169
  deletions.

Current breakage assessment:

- Package load is green.
- The basic input-driven driver workflow still runs for the independent H2 PQS
  readiness input.
- Old low-order/WL route-configured materialization, basis bundle export, and
  Ham bundle export through this private driver are intentionally blocked by
  `:route_configured_low_order_materializer_removed`.
- This cut does not validate final-basis/H1/H1-J/RHF driver inputs.

Recommended next cut:

- Continue separating driver construction from report/status vocabulary in
  `pqs_source_box_route_driver_helpers.jl`.
- Re-run compact driver smokes for final-basis/H1 stages only after the next
  source-side thinning pass, not as broad test gates.

## Checkpoint 7 - Detailed Driver Print Cut

Status:

- uncommitted demolition cut for review.

Deleted/simplified:

- Removed the large detailed print helper that dumped system, parent, unit,
  pair, and diagnostics field clouds.
- Kept the compact driver summary and the compact materialization print block.

Validation:

- `git diff --check` passed.
- Package load passed.
- The protected H2 independent PQS readiness driver smoke completed with
  saving disabled.

Line-count impact:

- 150 deleted lines in the uncommitted seventh cut.
- Current total branch pressure relative to `main`: 438 insertions and 45,319
  deletions.

Current breakage assessment:

- Package load and the readiness driver smoke remain green.
- Detailed text observability is intentionally reduced; future checks should
  use compact driver facts rather than print-line content.

Recommended next cut:

- Continue source-side thinning in `pqs_source_box_route_driver_helpers.jl`,
  focusing on report/status-only fields that are not needed for driver
  construction.

## Checkpoint 8 - Low-Order Report Field Cloud Cut

Status:

- uncommitted demolition cut for review.

Deleted/simplified:

- Removed the flat `low_order_*`, `route_core_*`, `pqs_prototype_*`, and
  `lw_complete_shell_*` report-field expansion from the driver report.
- Replaced the report-stage low-order summary with a compact status snapshot:
  status, materialization status/blocker, and broad authority booleans.
- Removed the report-stage PQS prototype and LW complete-shell report summary
  helpers that existed only to preserve detailed driver/report observability.
- Removed the retained-unit transform test block that directly asserted the
  private driver terminal-route wrapper. The remaining file still tests the
  retained-unit transform module contract directly.

Validation:

- `git diff --check` passed.
- Caller grep for the deleted report-stage helper/field names found no active
  source/bin/test/docs callers.
- Package load passed.
- The protected H2 independent PQS readiness driver smoke completed with
  saving disabled.

Line-count impact:

- 1,126 deleted lines in the uncommitted eighth cut after restoring compact
  driver survival summaries.

Current breakage assessment:

- This cut intentionally breaks consumers that key on deleted low-order flat
  report fields. That is report/status scaffolding, not the H2 survival route.
- The protected driver entry point and driver inputs remain present.
- Package load and the protected readiness driver smoke remain green.

Recommended next cut:

- If load and the readiness driver smoke remain green, continue auditing the
  remaining terminal-route state bookkeeping for construction use versus
  report-only use.

## Checkpoint 9 - Terminal Shellification Private Test Cut

Status:

- uncommitted demolition cut for review.

Deleted/simplified:

- Deleted `test/nested/cartesian_terminal_shellification_geometry_runtests.jl`.
  The test directly asserted private terminal-shellification helper names,
  scaffold object kinds, metadata/status vocabulary, and error strings.
- Removed that file from the slow integration runner.
- Kept `src/cartesian_terminal_shellification_geometry.jl` and the
  `CartesianShellification.raw_terminal_geometry` implementation because the
  driver still uses those helpers during construction.

Validation:

- `git diff --check` passed.
- Caller grep for the deleted test include found no active source/bin/test
  callers. The only remaining mention is this demolition note.
- Package load passed.
- The protected H2 independent PQS readiness driver smoke completed with
  saving disabled.

Line-count impact:

- About 858 test lines deleted before validation.

Current breakage assessment:

- This cut removes private helper/status test pressure only. It does not remove
  source geometry construction or driver inputs.
- Package load and the protected readiness driver smoke remain green.

Recommended next cut:

- Continue deleting private test pressure around old route/report scaffolding
  before attempting to simplify the terminal-route state source family.

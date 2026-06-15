# PQS Thin Route Demolition Status

## Promoted Thin Route and Temporary Line Ladders

The thin Cartesian/PQS route has been promoted to `main`. The demolition branch
name was retired after the fast-forward promotion.

Cartesian validation during the remaining migration should use temporary
driver line ladders, not the deleted helper/schema/status tests. A line ladder
answers one migration question: can this surviving Cartesian sub-line still run
the driver far enough to prove it exists?

Current line ladders are run with:

```text
julia --project=. tools/run_cartesian_line_ladder.jl --line=wl_atomic
julia --project=. tools/run_cartesian_line_ladder.jl --line=wl_diatomic
julia --project=. tools/run_cartesian_line_ladder.jl --line=pqs_atomic
julia --project=. tools/run_cartesian_line_ladder.jl --line=pqs_diatomic
```

The full 2x2x2 matrix remains:

```text
julia --project=. tools/run_cartesian_driver_ladder.jl
```

These ladders are migration scaffolding, not permanent architecture. When a
line is merged into the common driver route, its temporary ladder should be
deleted or folded into the main matrix. Do not reintroduce standalone
Cartesian nested helper tests, exact print-string tests, route payload schema
tests, or `status`/`available`/`blocker`/`readiness` field-cloud tests.

The old executable Cartesian test groups `nested`, `ordinary`, and `diatomic`
were removed from the default test dispatcher. Their remaining validation role
is covered by the driver ladders while the Cartesian lines are being merged.

## Demolition History

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

## Checkpoint 10 - Old Cartesian Export/Metadata Test Cut

Status:

- uncommitted demolition cut for review.

Deleted/simplified:

- Deleted the old Cartesian basis-bundle export test and overlap-projector
  bundle test. These preserved JLD2 key manifests and legacy artifact/export
  contracts rather than the protected driver mechanism.
- Deleted the unreferenced metadata-only parent Coulomb source summary test.
- Removed the deleted export tests from the slow integration runner.
- Kept compact module-contract tests such as the contracted-parent scaffold
  test because they exercise active constructors/audits rather than old driver
  report fields.

Validation:

- `git diff --check` passed.
- Caller grep for the deleted export/metadata test names found no active
  source/bin/test/docs callers.
- Package load passed.
- The protected H2 independent PQS readiness driver smoke completed with
  saving disabled.

Line-count impact:

- About 776 test lines deleted before validation.

Current breakage assessment:

- This cut removes old artifact/export/metadata test pressure. It does not
  remove driver inputs, the driver entry point, or numerical kernels.
- Package load and the protected readiness driver smoke remain green.

Recommended next cut:

- Continue using the driver as the survival test mechanism; avoid restoring old
  JLD2 manifest and export-key tests unless a public artifact contract is
  deliberately rebuilt.

## Checkpoint 11 - Legacy WL Seed and Atom-Growth Helper Cut

Status:

- uncommitted demolition cut for review.

Deleted/simplified:

- Deleted `src/white_lindsey_materialized_seed.jl`.
- Deleted `src/cartesian_atom_growth_route_driver_helpers.jl`.
- Removed both files from `src/GaussletBases.jl`.

Expected old-path fallout:

- Remaining source callers after deletion:
  - `src/cartesian_shellization_route.jl` still references
    `_white_lindsey_low_order_materialized_seed_inventory`.
  - `src/pqs_source_box_route_driver_helpers.jl` still references three
    `_pqs_source_box_route_driver_atom_growth_*` helper families from old
    low-order branches.
- These are old WL/atom-growth low-order scaffolding callers, not the protected
  independent H2 PQS readiness driver path. They should be cut or quarantined
  in a later repair/delete pass rather than restored.

Validation:

- `git diff --check` passed.
- Caller grep found the expected old-path source fallout listed above.
- Package load was attempted and failed first in `cartesian_bundle_export.jl`
  because the old export layer typed against the deleted
  `_WhiteLindseyLowOrderHamBundleAdapter`.

Line-count impact:

- 2,452 source lines deleted before validation.

Current breakage assessment:

- This cut intentionally breaks old WL seed and atom-growth low-order paths if
  those branches are invoked.
- The protected driver entry point and driver inputs remain present.
- The first load failure is old Cartesian bundle export scaffolding, not the H2
  survival route.

Recommended next cut:

- Remove or block the four remaining old-path callers listed above instead of
  restoring the deleted helper files.

## Checkpoint 12 - Cartesian Bundle Export Layer Cut

Status:

- uncommitted demolition cut for review.

Deleted/simplified:

- Deleted `src/cartesian_bundle_export.jl` and `src/cartesian_bundle_io.jl`.
- Removed the corresponding public exports, empty generic declarations, and
  includes from `src/GaussletBases.jl`.
- Deleted the nested hybrid orbital-transfer test that depended on writing and
  reading Cartesian basis bundles.
- Removed that test from the slow integration runner.

Expected old-path fallout:

- The ordinary-test `cartesian_basis_bundle_payload` call was removed in
  checkpoint 13 instead of restoring the old export API.

Validation:

- `git diff --check` passed.
- Caller grep found only the expected old export/test fallout and the four
  old low-order source caller fallouts.
- Package load passed after this cut.
- The protected H2 independent PQS readiness driver smoke completed with
  saving disabled.

Line-count impact:

- 1,815 source/test lines deleted before validation.

Current breakage assessment:

- This cut intentionally removes the old Cartesian basis/Hamiltonian bundle
  export API from the demolition branch.
- The protected driver entry point and driver inputs remain present.
- Package load and the protected readiness driver smoke are green.

Recommended next cut:

- If package load returns green, continue removing old-path callers such as the
  ordinary-test bundle assertion and the four low-order source caller fallouts.

## Checkpoint 13 - Deleted-Helper Caller Block Cut

Status:

- uncommitted demolition cut for review.

Deleted/simplified:

- Removed the one-center route materializer call to the deleted WL seed
  inventory helper by leaving the obsolete inventory slot empty.
- Replaced the old atom-growth low-order unit/transform/pair helper calls with
  explicit `nothing` or `:blocked_atom_growth_route_removed` summaries.
- Removed the ordinary-test bundle payload call and bundle-key assertions.

Validation:

- `git diff --check` passed.
- Caller grep for deleted WL seed, atom-growth, and bundle export symbols is
  clean outside this demolition note.
- Package load passed.
- The protected H2 independent PQS readiness driver smoke completed with
  saving disabled.

Line-count impact:

- Small source/test cut: about 18 net deleted lines before validation.

Current breakage assessment:

- Caller grep for the deleted WL seed, atom-growth, and bundle export symbols is
  now clean outside this demolition note.
- The old atom-growth/WL branches are now explicitly blocked rather than
  depending on deleted helper files.
- Package load and the protected readiness driver smoke are green.

Recommended next cut:

- Continue deleting old export/report/status tests; avoid reintroducing bundle
  or atom-growth compatibility just to satisfy old callers.

## Checkpoint 14 - Shellization Route Scaffold Cut

Status:

- uncommitted demolition cut for review.

Deleted/simplified:

- Deleted `src/cartesian_shellization_route.jl`, which had become old
  route-configured low-order/WL materialization scaffolding.
- Removed the include from `src/GaussletBases.jl`.
- Replaced its one surviving driver utility use with a tiny local
  `_cartesian_parent_location_tuple` in `pqs_source_box_route_driver_helpers.jl`.

Validation:

- `git diff --check` passed.
- Caller grep found no old route-scaffold code references. Remaining mentions
  are this demolition note plus the new local parent-location helper.
- Package load passed.
- The protected H2 independent PQS readiness driver smoke completed with
  saving disabled.

Line-count impact:

- About 1,014 net source lines deleted before validation.

Current breakage assessment:

- This removes old shellization-route materializer/report scaffolding. The
  module-owned shellification and terminal-geometry code remains.
- The protected driver entry point and driver inputs remain present.
- Package load and the protected readiness driver smoke are green.

Recommended next cut:

- Audit remaining low-order shellification plan materializers in
  `cartesian_shellification_plan.jl`; cut only old WL/atom-growth materializer
  surfaces, not shared shellification primitives.

## Checkpoint 15 - Old Low-Order Shellification Plan Cut

Status:

- uncommitted demolition cut for review.

Deleted/simplified:

- Deleted `src/cartesian_shellification_plan.jl`, the old private low-order
  shellification plan/materializer layer.
- Removed its include from `src/GaussletBases.jl`.
- Collapsed the driver helper's atom-growth shell-stage planner to an explicit
  `:blocked_atom_growth_route_removed` tombstone instead of carrying the old
  WL/atom-growth shellification scaffold.

Validation:

- `git diff --check` passed.
- Caller grep for the deleted shellification-plan helpers is clean outside this
  demolition note.
- Package load passed.
- The protected H2 independent PQS readiness driver smoke completed with saving
  disabled.

Line-count impact:

- About 2,396 net source lines deleted before validation.

Current breakage assessment:

- This removes old low-order/WL/atom-growth shellification materializer
  scaffolding.
- The protected driver entry point, driver inputs, and independent H2 readiness
  route remain present and green.
- The module-owned `CartesianShellification` and terminal-lowering route spine
  remain intact.

Recommended next cut:

- Continue auditing old driver/report status surfaces. Avoid deleting
  construction-bearing independent H2 PQS source/final/H1/H1-J/RHF code; target
  report/status wrappers and obsolete old-flat tests.

## Checkpoint 16 - Old Nested Reporting and Experimental Export Cut

Status:

- uncommitted demolition cut for review.

Deleted/simplified:

- Deleted `src/cartesian_nested_reporting.jl`, including the doside/COMX trace
  diagnostic writer and old nested fixed-block timing report implementation.
- Deleted `src/bond_aligned_diatomic_geometry_export.jl`.
- Deleted `src/experimental_chain_export.jl`.
- Removed their public exports, empty generic declarations, and includes from
  `src/GaussletBases.jl`.
- Deleted ordinary/diatomic/docs tests and public-doc references that mainly
  protected exact experimental export metadata, text output, and doside trace
  diagnostics.

Validation:

- `git diff --check` passed.
- Caller grep for deleted export/trace APIs is clean in source/tests/bin/public
  docs.
- Package load passed.
- The protected H2 independent PQS readiness driver smoke completed with saving
  disabled.

Line-count impact:

- About 1,942 net source/test/doc lines deleted before validation.

Current breakage assessment:

- `cartesian_nested_atomic.jl` still has optional timing-call sites for
  `_nested_capture_timeg_report`. Package load does not require that path, but
  optional nested fixed-block timing is now a live-looking fallout if invoked.
- The deleted export/trace surfaces were old diagnostic/producer-side
  bureaucracy, not part of the protected Cartesian/PQS driver survival path.
- The protected driver entry point and driver inputs remain present and green.

Recommended next cut:

- Either leave optional nested timing broken until the smaller system is
  repaired, or restore only a tiny timing helper without restoring the deleted
  doside/export reporting layer. Continue targeting report/status wrappers
  rather than numerical kernels.

## Checkpoint 17 - Compact WL Capability Repair

Status:

- uncommitted repair checkpoint for review.
- Correction: WL is first-class construction capability, not old bloat. The
  demolition branch should delete old WL/QW report/export/scaffold surfaces, but
  it must preserve or compactly rebuild WL physics-producing driver capability.

Classification of deleted WL surfaces:

- Real numerical/construction kernels: the deleted
  `white_lindsey_overlap/kinetic/position/x2/electron_nuclear/one_body` block
  files, density-density construction, and unit-coefficient adapters contained
  real WL operator construction logic.
- Route/final-basis construction: `white_lindsey_materialized_seed.jl` and
  decomposed unit-pair inventory carried useful construction material mixed with
  old route/report inventory.
- Scaffold/report glue: adapter summaries, seed-oracle summaries, route-global
  status wrappers, precompile workloads, and broad nested WL tests mostly
  protected old flat vocabulary and were not restored.
- Tests protecting physics output: old He/H WL acceptance tests protected real
  scalar capability, but were too broad and scaffold-heavy to restore wholesale.

Restored capability:

- Added a compact atomic pure-gausslet WL materialization path through
  `cartesian_materialization(...)`.
- The path uses the surviving `one_center_atomic_full_parent_fixed_block`
  numerical builder, computes the finite/symmetric one-body Hamiltonian, and
  reports a small scalar summary: retained dimension, overlap identity error,
  H1 lowest eigenvalue, and symmetry/finite checks.
- Added compact JLD2 writing only when `save_basis_artifact` or
  `save_ham_artifact` is requested. These files contain basis/operator arrays
  and a short provenance summary, not the old report/status field cloud.
- Added an opt-in driver input and survival smoke:
  `test/driver_inputs/he_wl_q5_pure_gausslet_h1.jl` and
  `test/nested/cartesian_ham_builder_he_wl_q5_pure_gausslet_h1_runtests.jl`.

Validation:

- `git diff --check` passed.
- Package load passed.
- WL driver survival smoke passed, including compact Ham JLD2 write/read:
  retained dimension `419`, H1 lowest `-1.991344469963436`, H1 finite,
  H1 symmetry error `0.0`, overlap identity error about `5.22e-15`.
- Protected H2 independent PQS readiness driver smoke completed with saving
  disabled and still reported support counts
  `(atom_contact_core = 275, shared_shell_1 = 578, shared_shell_2 = 362)`.

Line-count impact:

- Repair adds about 231 lines before review: 185 tracked source lines plus a
  small 11-line input and 35-line opt-in smoke.
- This is an intentional capability-regression repair, not a general exception
  for rebuilding deleted report scaffolding.

Current breakage assessment:

- WL atomic pure-gausslet H1 is physics-producing again through the driver.
- WL density/H1-J/RHF and diatomic WL materialization are not restored in this
  checkpoint.
- Compact basis/Ham JLD2 output exists for this WL atomic seam; broader
  consumer artifact contracts still need deliberate design.

Recommended next repair/delete step:

- Keep repairing WL by compact construction capability, not old adapter
  summaries. Next WL work should decide the minimal compact path for
  density/H1-J/RHF or diatomic pure-gausslet WL, while keeping report/status
  alias surfaces deleted.

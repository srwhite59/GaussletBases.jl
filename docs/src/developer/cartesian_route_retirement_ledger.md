# Cartesian Route Retirement Ledger

This developer ledger tracks old Cartesian, White--Lindsey, and PQS route
surfaces while the implementation moves toward the module spine:

```text
CartesianCPB
-> CartesianShellification
-> CartesianTerminalLowering
-> CartesianRetainedUnits
-> CartesianRetainedUnitTransformContracts
-> CartesianUnitPairs
-> CartesianPairOperatorPlans
-> CartesianPairBlockMaterialization
```

Use this file to answer: "Which old authority surface is this pass moving
toward adapter, oracle, or retired?" The ledger is not a public API contract.
It is a cleanup map for preventing old route glue, report aliases, and oracle
paths from becoming new route authority.

| Old surface / file | Current role | Replacement module/stage | Equivalence test | Status | Deletion condition | Notes |
| --- | --- | --- | --- | --- | --- | --- |
| `src/pqs_source_box_route_driver_helpers.jl` | Legacy route-driver glue and report compatibility aliases. | Module-owned route state/summaries and future route adapter over the module spine. | Focused driver/report compatibility checks comparing compact structured summaries to the still-required report aliases. | compatibility | Final reports and downstream consumers no longer depend on flat scalar aliases; route code consumes module-owned objects directly. | No new route concepts should be introduced here. |
| `src/cartesian_pair_block_materialization/white_lindsey_*` and old White--Lindsey nested kernels such as `_nested_doside_1d`, `_nested_face_product`, `_nested_edge_product`, `_nested_corner_piece` | Old-kernel reuse behind the White--Lindsey boundary-stratum adapter. | `CartesianPairBlockMaterialization` White--Lindsey boundary-stratum adapter. | Focused local adapter equivalence tests against selected old-kernel/fixed-block slices for retained-unit coefficients and low-order one-body blocks. | adapter/oracle | New spine can produce and validate all low-order one-body blocks and global overlap/one-body placement without treating old kernels as route authority. | Old kernels may feed adapters or oracle comparisons; they should not decide route structure. |
| Old White--Lindsey materialized seed / fixed-block safe one-body matrices | Oracle/reference for selected safe one-body global matrix slices. | `CartesianPairBlockMaterialization` route-global safe one-body path: route-local one-body block collection -> term-specific placement plan -> term-specific dense global retained matrix. | `test/nested/cartesian_pair_block_global_overlap_oracle_runtests.jl` covers selected old-oracle overlap slices through `White-Lindsey local adapter overlap blocks -> local one-body block collection -> overlap placement plan -> global overlap matrix`: face/face `162:170, 162:170`, face/edge `162:170, 210:212`, and transpose-filled edge/face `210:212, 162:170`. Observed max absolute errors were `1.7636953211319766e-15`, `2.6343504968563524e-46`, and `2.6343504968563524e-46`, respectively. `test/nested/cartesian_pair_block_global_kinetic_oracle_runtests.jl` covers selected White--Lindsey global kinetic equivalence through `White-Lindsey local adapter kinetic blocks -> local one-body block collection -> kinetic placement plan -> global kinetic matrix` on the same slices. Observed kinetic max absolute errors were `3.4638958368304884e-14`, `3.1907967999947086e-30`, and `3.1907967999947086e-30`, respectively. `test/nested/cartesian_pair_block_global_position_oracle_runtests.jl` covers selected global `position_x`, `position_y`, and `position_z` equivalence through `White-Lindsey local adapter position blocks -> local one-body block collection -> position placement plan -> global position matrix` on the same slices. Observed position_x max absolute errors were `1.7435358081325719e-15`, `9.857551802004026e-46`, and `9.857551802004026e-46`; position_y errors were `8.459873058838228e-16`, `8.758115402030107e-46`, and `8.758115402030107e-46`; position_z errors were `8.474501317837309e-16`, `4.926439913641935e-47`, and `4.926439913641935e-47`. `test/nested/cartesian_pair_block_global_x2_oracle_runtests.jl` covers selected global `x2_x`, `x2_y`, and `x2_z` equivalence through `White-Lindsey local adapter x2 blocks -> local one-body block collection -> x2 placement plan -> global x2 matrix` on the same slices. Observed x2_x max absolute errors were `1.844473375912524e-15`, `6.0243961376664915e-34`, and `6.0243961376664915e-34`; x2_y errors were `4.541938536395539e-16`, `1.7964170116877891e-32`, and `1.7964170116877891e-32`; x2_z errors were `4.577863961878137e-16`, `2.7369110631344083e-47`, and `2.7369110631344083e-47`. Synthetic global retained matrix pilots now exist for `:overlap`, `:kinetic`, `:position_x`, `:position_y`, `:position_z`, `:x2_x`, `:x2_y`, and `:x2_z`. | adapter/oracle-only; partially replaced for selected face/face and face/edge safe one-body slices | Broader retained-unit inventory equivalence, broader local/global block equivalence, and route-driver adoption must land before deleting old authority. | This is selected-slice overlap, kinetic, position, and x2 equivalence plus synthetic safe one-body global assembly coverage, not full-route replacement. Keep the old path as validation authority. No full White--Lindsey route assembly, driver integration, Hamiltonian assembly, term summing, Coulomb, IDA/MWG, PQS projection/Lowdin, exports, or artifacts are claimed here. |
| Raw PQS support-row contraction / shell-row fallback paths | Debug/oracle only. | Raw product source-space blocks plus shell-realization bridge. | PQS source-space safe-term blocks plus explicit projection/Lowdin realization match the support-row or shell-row oracle cases. | oracle | PQS source-space blocks plus explicit projection/Lowdin realization match oracle cases. | Must not become the production PQS algorithm. Shell/support rows are compatibility and debug surfaces. |
| Flat `terminal_shellification_*` / `route_core_typed_pair_operator_*` report aliases | Final report compatibility. | Compact summaries from `terminal_route_state`, `pair_operator_summary`, and module summaries. | Report compatibility checks comparing alias values to compact module summaries, without copying large staged metadata objects. | compatibility | Reports, tests, and users consume structured summaries rather than flat aliases. | Do not add new flat driver fields. Derive required aliases at the report boundary only. |
| Local one-body block collection | New local result organizer, not old code. | Future global placement/assembly layer. | Focused tests that local entries, skipped records, summaries, and materialization flags match the mixed one-body block-set consumer. | active | Not a deletion target; becomes typed if it crosses a stable module boundary. | It organizes local pair-block results and skips only. It does not place global matrices, sum terms, build Hamiltonians, export artifacts, or change IDA/MWG semantics. |

## Rules For Future Passes

- New work should move one ledger row toward adapter, oracle, or retired.
- Do not add new flat driver fields.
- Do not add metadata-only readiness layers unless they replace a duplicated
  surface, feed overlap/global placement, or retire an old route surface.
- Old code may be used as adapter or oracle, not as new route authority.
- Deletion requires an explicit equivalence test and a stated deletion
  condition.

## Current Replacement Pressure

- Selected global overlap, kinetic, position_x, position_y, position_z, x2_x,
  x2_y, and x2_z equivalence are now covered for White--Lindsey face/face and
  face/edge slices.
- Synthetic global retained matrix pilots now exist for the safe one-body
  terms `:overlap`, `:kinetic`, `:position_x`, `:position_y`, `:position_z`,
  `:x2_x`, `:x2_y`, and `:x2_z`.
- A private route-shaped global one-body adapter now composes
  `PairBlockMaterializationPlan` inputs through the local block collection,
  term-specific placement plans, and term-specific dense global retained matrix
  pilots for the individual safe one-body terms `:overlap`, `:kinetic`,
  `:position_x`, `:position_y`, `:position_z`, `:x2_x`, `:x2_y`, and `:x2_z`.
  This is module-level adapter coverage, not route-driver wiring, term summing,
  one-body Hamiltonian construction, or Hamiltonian assembly.
- A private route-shaped safe one-body matrix-set adapter now wraps individual
  `route_global_one_body_matrix(...; term = ...)` calls and returns
  term-separated result objects plus a compact summary. It does not copy dense
  matrices out of term results, allows partial success with blocked per-term
  results for unsupported terms, and is still not term summing, a one-body
  Hamiltonian object, or Hamiltonian assembly. PQS source-space records remain
  blocked inside term results until shell projection/Lowdin realization exists.
- A private overlap-only route-state adapter now accepts structured state that
  already carries a `PairBlockMaterializationPlan` and delegates to the existing
  route-global overlap matrix path. This is replacement pressure only: it is
  not production route-driver wiring and does not retire the old route by
  itself.
- A private overlap-only driver-facing hook,
  `driver_global_overlap_result`, delegates structured state carrying
  `pair_block_materialization_plan` to the route-state overlap adapter. It
  blocks with `:missing_pair_block_materialization_plan` when no structured
  pair-block materialization plan is present. It is not production driver
  wiring and does not retire the old route by itself.
- A private opt-in overlap driver option,
  `private_global_overlap_requested`, carries
  `private_global_overlap_result` plus `private_global_overlap_summary` from the
  materialization stage. It returns the structured overlap result when enough
  structured state, global dimension, and overlap inputs are present, or a
  compact blocked summary otherwise. This remains nonproduction and off by
  default; it is not old-route retirement.
- Selected White--Lindsey old-oracle equivalence currently covers overlap,
  kinetic, position_x, position_y, position_z, x2_x, x2_y, and x2_z.
- `src/cartesian_pair_block_materialization/one_body_global_matrix_helpers.jl`
  owns only behavior-neutral symmetric placement validation/insertion. The
  term-specific global matrix files still own result statuses, blockers,
  object kinds, and materialization flags.
- Route-global adapter validation should stay focused:
  `cartesian_pair_block_route_global_one_body_adapter_runtests.jl` is the
  broader individual-term adapter contract at about 40--45 seconds, and
  `cartesian_pair_block_route_global_matrix_set_smoke_runtests.jl` is the
  focused matrix-set smoke at about 40 seconds after slimming. Do not add these
  casually to broader default runners; old-oracle tests remain gate/oracle
  tests, not routine per-pass checks.
- The next replacement slice should decide whether to broaden selected
  safe-term oracle coverage or move toward carefully scoped driver wiring.
- The current global one-body pilots do not assemble Hamiltonians, build
  one-body Hamiltonian objects, sum terms, build Coulomb, touch IDA/MWG, build
  density-density, nuclear attraction, or Gaussian local terms, perform PQS
  projection/Lowdin, wire route drivers, export artifacts, or assemble a full
  White--Lindsey route.

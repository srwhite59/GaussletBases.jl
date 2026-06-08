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
| Old fixed-block / materialized seed overlap path | Oracle for selected White--Lindsey local/global overlap checks. | Local one-body block collection plus overlap placement plan and global overlap matrix. | `test/nested/cartesian_pair_block_global_overlap_oracle_runtests.jl` covers selected old-oracle overlap slices through `White-Lindsey local adapter overlap blocks -> local one-body block collection -> overlap placement plan -> global overlap matrix`: face/face `162:170, 162:170`, face/edge `162:170, 210:212`, and transpose-filled edge/face `210:212, 162:170`. Observed max absolute errors were `1.7636953211319766e-15`, `2.6343504968563524e-46`, and `2.6343504968563524e-46`, respectively. | oracle; selected global-overlap equivalence covered | Full retirement still requires broader focused equivalence tests covering retained-unit inventory, local blocks, placement, and global overlap/one-body assembly. | This is selected-slice global overlap equivalence only, not full-route overlap replacement. Keep as validation authority. No kinetic, position, x2, driver integration, Hamiltonian, Coulomb, IDA/MWG, or PQS realization coverage is claimed here. |
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

- Selected global overlap equivalence is now covered for White--Lindsey
  face/face and face/edge slices.
- The next replacement slice should extend the same placement/global assembly
  pattern to one additional safe one-body term, likely kinetic, before driver
  wiring.

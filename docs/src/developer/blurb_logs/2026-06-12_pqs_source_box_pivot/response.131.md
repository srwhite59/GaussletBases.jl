Pass 131 response - private complete-core/shell Ham payload seam

Commit:
- `08c851ba Add private PQS complete core shell Ham payload`

Files changed:
- `src/pqs_source_box_route_driver_helpers.jl`
- `test/nested/pqs_source_box_route_driver_complete_core_shell_ham_payload_runtests.jl`

Implementation summary:
- Added private helper `_pqs_source_box_route_driver_complete_core_shell_ham_payload(...)`.
- Added private blocked-payload and missing-input blocker helpers for the Ham seam.
- Added compact slot `complete_core_shell_ham_payload` to `_PQSCompleteCoreShellDiagnosticRoutePayload`.
- Wired `_pqs_source_box_route_driver_complete_core_shell_diagnostic_route_payload(...)` to build the Ham payload from route-owned objects:
  - `source_payload`
  - `final_basis`
  - materialized H1 payload and `final_hamiltonian`
  - `density_inputs`
  - H1/J diagnostic route payload
  - nested `density_interaction`
- Kept report-facing fields unchanged. No public API, report aliases, export/artifact writing, RHF route wiring, SCF/DIIS work, full four-index Coulomb tensor, or IDA/MWG promotion.

Final payload/slot:
- Slot: `complete_core_shell_ham_payload`
- Object kind: `:pqs_source_box_complete_core_shell_ham_payload`
- Success status: `:materialized_pqs_source_box_complete_core_shell_ham_payload`
- Electron-electron representation: `:pre_final_density_interaction`
- Required labels preserved:
  - `density_gauge = :pre_final_localized_positive_weight`
  - `raw_pair_factor_convention = :raw_numerator`
  - `support_row_order = :core_then_shell`
  - `signed_final_weight_division_used = false`
  - `raw_no_division_used = false`
  - `density_normalized_pair_terms_used_as_authority = false`
  - `public_api = false`
  - `exports_materialized = false`
  - `artifacts_materialized = false`
  - `rhf_product_surface = false`
  - `serious_hf_claim = false`

Validation:
```text
julia --project=. test/nested/pqs_source_box_route_driver_complete_core_shell_ham_payload_runtests.jl
```
- Passed: 21/21.
- Reported test time: `1m52.3s`.
- Cold command included package precompilation: `GaussletBases` precompiled in about 56.5s.

```text
julia --project=. -e 'using GaussletBases; println("load ok")'
```
- Passed: `load ok`.

```text
git diff --check
```
- Passed.

Git status after commit:
```text
## main...origin/main [ahead 1]
```

Deletion/shrinkage report:
- deleted: none.
- simplified: complete core/shell downstream Ham handoff now has one private object boundary instead of requiring consumers to coordinate H1/J report aliases and RHF validator state.
- quarantined: RHF remains validator-only; report/export/artifact paths remain outside this private payload; pre-final density interaction is explicitly labeled and not promoted to IDA/MWG or full Coulomb tensor semantics.
- not deleted because: existing H1/J report aliases and materialization/export fields are compatibility surfaces and this pass was scoped to adding the private seam only.
- exact remaining caller/blocker: the new payload is available through `assembly.complete_core_shell_diagnostic_route_payload.complete_core_shell_ham_payload`; no downstream CR2/HFDMRG/export caller has been added yet.

-- repo-doer@macmini

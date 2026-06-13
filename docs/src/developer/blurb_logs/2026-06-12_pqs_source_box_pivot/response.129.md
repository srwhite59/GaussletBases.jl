Pass 129 response - no-edit route RHF adoption audit

Scope:
- No source/test/doc edits.
- No Julia commands run; read-only inspection was sufficient.
- Inspected the requested route-driver, private RHF helper, framework, and fixture-policy surfaces.

Current source inventory:
- `src/pqs_source_box_route_driver_helpers.jl`
  - `_pqs_source_box_route_driver_complete_core_shell_source_plan_payload(...)` builds the private source-plan payload: region plan, source plan, Coulomb expansion, status/blocker/missing-input summary.
  - `_pqs_source_box_route_driver_complete_core_shell_h1_payload(...)` builds the complete core/shell final basis and H1 payload.
  - `_pqs_source_box_route_driver_complete_core_shell_density_inputs(...)` derives axis weights and raw pair-factor terms from the source-plan bundles and Coulomb expansion; its summary explicitly marks private diagnostic use and rejects retained diagnostic weights as authority.
  - `_pqs_source_box_route_driver_complete_core_shell_h1_j_diagnostic_payload(...)` consumes region/source plan, final basis, H1 payload, density inputs, and Coulomb expansion to materialize or block the H1/J diagnostic payload.
  - `_PQSCompleteCoreShellDiagnosticRoutePayload` is already the compact route-owned authority for the complete core/shell diagnostic path. It carries `status`, `blocker`, `route_family`, `source_payload`, `final_basis`, `h1_payload`, `density_inputs`, `h1_j_payload`, `missing_inputs`, `summary`, and `metadata`.
  - `cartesian_assembly(...)` builds `complete_core_shell_diagnostic_route_payload`, then keeps compatibility aliases for `complete_core_shell_h1_j_diagnostic_payload`, summary, status, and blocker.
  - `_pqs_source_box_route_driver_complete_core_shell_h1_j_report_fields(...)` exposes only summary/status/scalar aliases: final dimension, H1 energy, self-Coulomb, density gauge, and driver-route materialization.
- `src/pqs_multilayer_complete_core_shell_rhf.jl`
  - `_pqs_multilayer_complete_core_shell_rhf_input_contract(...)` requires real route-owned objects: source plan, final basis, H1 payload, density inputs, Coulomb expansion, explicit electron count, and explicit fixture role. No scalar report aliases are sufficient.
  - `_pqs_multilayer_complete_core_shell_rhf_scf_control_payload(...)` is the private SCF control object. It defaults to `mixing_kind = :fixed_point`, but when called with `mixing_kind = :fock_diis` and omitted `max_history`, it resolves `max_history = 8`.
  - `_pqs_multilayer_complete_core_shell_rhf_scf_payload(...)` consumes the input contract plus H1/J density interaction, optional initial density, and SCF controls. Materialized summaries remain private diagnostic/non-public and explicitly keep `driver_route_materialized = false`, `route_report_materialized = false`, `public_api = false`, `exports_materialized = false`, and `artifacts_materialized = false`.

Recommended request object / keyword shape:
- Add one private route request object or compact keyword group, not broad recipe/report plumbing:
  - `private_complete_core_shell_rhf_requested::Bool`
  - `electron_count`
  - `fixture_role`
  - optional `scf_control_payload`
  - optional SCF control overrides only when no control payload is supplied: `mixing_kind = :fock_diis`, `max_iterations`, `density_atol`, `energy_atol`, `residual_atol`, `trace_atol`, `idempotency_atol`, `max_history = nothing`, `diis_start_iteration`, `diis_regularization`, `diis_coefficient_max_abs`
  - `metadata`
- I would name the private request helper/object `_pqs_source_box_route_driver_complete_core_shell_rhf_diagnostic_request(...)` or `_PQSCompleteCoreShellRHFDiagnosticRequest` if a struct is preferred. Keep it private and local to `src/pqs_source_box_route_driver_helpers.jl`.
- Required request discipline:
  - RHF must be opt-in via the private request flag.
  - `electron_count` must be explicit; do not infer it from nuclear charge or center records.
  - `fixture_role` must be explicit; use `:route_smoke` for the current one-center driver dry-run role.
  - Default route-driver RHF request should pass `mixing_kind = :fock_diis` and omit `max_history` so the validated private SCF default resolves to 8.

Recommended consumption point:
- Consume the request immediately after the existing complete core/shell diagnostic route payload has `source_payload`, `final_basis`, `h1_payload`, `density_inputs`, and `h1_j_payload`.
- The cleanest seam is inside `_pqs_source_box_route_driver_complete_core_shell_diagnostic_route_payload(...)`, after `h1_j_payload` is built. That function already has the complete route-owned object bundle needed by `_pqs_multilayer_complete_core_shell_rhf_input_contract(...)`.
- Avoid building RHF from `cartesian_report(...)` aliases. The report only has summaries/scalars and is the wrong authority boundary.

Recommended private payload / slot:
- Helper name:
  - `_pqs_source_box_route_driver_complete_core_shell_rhf_diagnostic_payload(...)`
- Slot name:
  - `complete_core_shell_rhf_diagnostic_payload`
- Best placement:
  - Add `rhf_diagnostic_payload` or `complete_core_shell_rhf_diagnostic_payload` as a new field on `_PQSCompleteCoreShellDiagnosticRoutePayload`.
  - Then expose an assembly-level compatibility alias only if needed by the focused route-driver test.
- Payload shape should stay compact:
  - `object_kind`
  - `status`
  - `blocker`
  - `route_family`
  - `requested`
  - `input_contract`
  - `scf_control_payload`
  - `scf_payload`
  - `missing_inputs`
  - `summary`
  - `metadata`
- Do not add report fields in the first implementation pass. If later needed, expose a summary-only group: status, blocker, requested, converged, iteration count, residual metric/value, final energy, and nonclaim flags. Do not copy the SCF field cloud into route reports.

Required blockers:
- `:complete_core_shell_rhf_not_requested`
- `:not_applicable_complete_core_shell_rhf_non_pqs_source_box_route`
- `:missing_complete_core_shell_h1_j_diagnostic_payload`
- `:missing_complete_core_shell_h1_j_route_inputs`
- `:missing_electron_count`
- `:missing_fixture_role`
- `:invalid_electron_count`
- `:unsupported_open_shell_rhf_input`
- `:unsupported_fixture_role`
- `:invalid_scf_controls`
- `:missing_rhf_input_contract`
- `:missing_rhf_scf_inputs`
- `:scf_not_converged`

Required nonclaims:
- No route production materialization.
- No public API.
- No exports.
- No artifacts.
- No GTO behavior.
- No IDA/MWG behavior.
- No density-density production route behavior.
- No physics endpoint claim.
- No serious-HF claim.
- H1/J remains private diagnostic support for route smoke only.
- `rhf_materialized` may be true inside the private SCF payload when converged, but the route/report/public flags must remain false: `driver_route_materialized = false`, `route_report_materialized = false`, `public_api = false`, `exports_materialized = false`, `artifacts_materialized = false`.

Smallest next implementation pass:
- File to touch:
  - `src/pqs_source_box_route_driver_helpers.jl`
- Focused test to touch or add:
  - Prefer a narrow route-driver/report test that exercises `cartesian_assembly(...)` for the one-center source-box driver dry-run and requests private complete core/shell RHF diagnostic materialization.
  - If an existing focused route-driver H1/J test already covers this route payload, extend that test minimally; otherwise add one narrow nested test. Do not broaden to integration suites.
- Implementation steps:
  1. Add private RHF diagnostic request helper/object with explicit `electron_count`, explicit `fixture_role`, and default `mixing_kind = :fock_diis`.
  2. Add `_pqs_source_box_route_driver_complete_core_shell_rhf_diagnostic_payload(...)` that blocks when not requested or when route/H1/J/RHF inputs are unavailable.
  3. Build `_pqs_multilayer_complete_core_shell_rhf_input_contract(...)` from `source_payload.source_plan`, `final_basis`, `h1_payload`, `density_inputs`, `source_payload.coulomb_expansion`, `electron_count`, and `fixture_role`.
  4. Build `_pqs_multilayer_complete_core_shell_rhf_scf_payload(...)` from the input contract plus the existing `h1_j_payload`/density interaction and request control payload.
  5. Add the resulting private RHF payload as a compact slot on `_PQSCompleteCoreShellDiagnosticRoutePayload`.
  6. Keep report-facing fields unchanged in the first pass unless the focused test needs a private assembly alias.
- This seam is clear. It does not require broad route option plumbing or a new report-field cloud. The existing H1/J payload structure is sufficient for the RHF input contract.

Git status:
```text
## main...origin/main
```

Deletion/shrinkage forecast:
- deleted: none in this no-edit audit.
- simplified: after the payload slot exists, route assembly can stop coordinating RHF readiness through adjacent scalar/status groups.
- quarantined: RHF stays private diagnostic route-smoke only, behind explicit request/electron-count/fixture-role gates.
- not deleted because: current H1/J report aliases are compatibility surface and the task explicitly forbids report field changes.
- exact remaining caller/blocker: no route caller exists yet; first implementation must add the private request path and one focused assembly-level exercise before any report or downstream cleanup.

-- repo-doer@macmini

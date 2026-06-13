Pass 110 response - RHF route-integration seam audit

Scope:
- No code changes.
- No Julia commands.
- Read-only inspection only.

1. Current driver payload sufficiency:
- The current `_PQSCompleteCoreShellDiagnosticRoutePayload` already carries the route-owned data needed for private RHF SCF except the explicit RHF diagnostic inputs:
  - `source_payload` with region/source plan and Coulomb expansion;
  - `final_basis`;
  - materialized `h1_payload`;
  - route-owned `density_inputs`;
  - materialized H1/J diagnostic payload when density inputs are available.
- That is enough structured route data for the private RHF input contract, one-step payload, initial-density payload, and SCF payload once explicit `electron_count` and `fixture_role` are supplied.
- It does not currently carry `electron_count`, `fixture_role`, or SCF controls, and it should not infer them from center charges or nuclei.

2. Smallest acceptable explicit-input seam:
- Introduce one private structured request object, not report-field scalars.
- Recommended name:
  `_PQSCompleteCoreShellRHFDiagnosticRequest`
  or a NamedTuple payload from `_pqs_source_box_route_driver_complete_core_shell_rhf_diagnostic_request(...)`.
- Minimal fields:
  - `status`
  - `blocker`
  - `requested::Bool`
  - `electron_count`
  - `fixture_role`
  - `scf_controls` as one compact object/NamedTuple
  - `missing_inputs`
  - `metadata`
- `electron_count` and `fixture_role` should come only from this explicit private diagnostic request. They should not be derived from route center records, nuclear charge, system metadata, or compact fixture geometry.

3. Consumption point:
- Consume it inside `_pqs_source_box_route_driver_complete_core_shell_diagnostic_route_payload(...)`, adjacent to H1/J payload construction.
- The diagnostic route payload is the right owner because it already coordinates source plan, final basis, H1, density inputs, and H1/J.
- Add structured fields to `_PQSCompleteCoreShellDiagnosticRoutePayload` only when implementation is approved, for example:
  - `rhf_request`
  - `rhf_input_contract`
  - `rhf_scf_payload`
  - `rhf_summary`
- Do not consume this in `cartesian_report(...)` in the first implementation pass.

4. Public options/report fields:
- Keep route-driver public options unchanged in the first implementation pass.
- Keep report-facing fields unchanged in the first implementation pass.
- If status visibility is needed, keep it inside the structured assembly/diagnostic route payload first, not as a scalar report-field cloud.

5. Validation recommendation:
- First implementation pass should use a synthetic/private helper test.
- A cold route dry-run is not required just to introduce the request object and private diagnostic payload seam.
- A focused route dry-run becomes appropriate only when a real route-driver path is made to pass an explicit RHF request into `cartesian_assembly` and materialize RHF on the compact 223-dimensional fixture.

Recommended next implementation pass:
- Add the private request object/helper and a private route-owned RHF diagnostic payload helper in `src/pqs_source_box_route_driver_helpers.jl`.
- Wire only `_pqs_source_box_route_driver_complete_core_shell_diagnostic_route_payload(...)` to carry a blocked/not-requested RHF diagnostic state by default and to materialize private RHF SCF only when an explicit request object is supplied.
- Add one synthetic/private test that builds a compact diagnostic route payload from synthetic pieces or a narrow helper call and verifies:
  - missing request gives not-requested or blocked-with-missing-explicit-input status;
  - supplied request carries `electron_count` and `fixture_role`;
  - public report fields are unchanged.

Exact code surfaces to touch next if approved:
- `src/pqs_source_box_route_driver_helpers.jl`
  - `_PQSCompleteCoreShellDiagnosticRoutePayload`
  - `_pqs_source_box_route_driver_complete_core_shell_diagnostic_route_payload(...)`
  - new private RHF request/helper adjacent to the H1/J diagnostic helper family
- likely one new focused test under `test/nested/`
- no report helper change in the first pass unless manager explicitly asks for status exposure.

Blockers and risks:
- Blocker: no existing explicit `electron_count`/`fixture_role` carrier in the route-driver diagnostic payload.
- Blocker: using center charges would violate the explicit-input contract and fixture policy.
- Risk: adding public recipe options now would imply route behavior adoption before fixture policy and physics endpoint decisions are made.
- Risk: adding scalar report aliases would recreate the field-cloud problem the compact payloads were meant to avoid.
- Risk: materialized private SCF on the compact route smoke could be mistaken for physics endpoint validation unless the request carries `fixture_role = :route_smoke` and nonclaims remain explicit.
- Risk: first real route dry-run may mostly measure compilation; use focused validation only when actual route consumption is introduced.

Validation status:

```text
## main...origin/main
```

Deletion/shrinkage forecast:
- deleted: none in the next request-object pass.
- simplified: later, report aliases can stay unchanged while internal diagnostic state moves into a compact object.
- quarantined: RHF SCF should remain private diagnostic-only and request-gated.
- not deleted because: H1/J report fields are still existing compatibility surface and should not be churned during the first RHF request pass.
- exact remaining caller/blocker: no route-driver caller can supply explicit `electron_count` and `fixture_role` yet; the missing explicit diagnostic request object is the next blocker.

-- repo-doer@macmini

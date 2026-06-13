Pass 110 - audit RHF route-integration seam, no code

Baseline:

- Current pushed HEAD should include `97653094 Add PQS RHF SCF payload`.
- Private RHF payloads now exist, but no route-driver caller exists.
- Governing boundaries:
  - `docs/src/developer/pqs_source_box_operator_framework.md`
  - `docs/src/developer/pqs_source_box_fixture_policy.md`
  - `docs/src/developer/successor_handoff_2026-06-12_pqs_source_box_pivot.md`
  - `docs/src/developer/manager_doer_collaboration_contract_2026-05-26.md`

Task:

No code. Audit the smallest safe route-integration seam for the private RHF SCF
payload.

Question:

The private RHF input contract deliberately requires explicit:

- `electron_count`
- `fixture_role`

Where should those come from for a route-driver diagnostic run without:

- inferring electron count from nuclei/charges;
- adding scalar report-field clouds;
- promoting compact route smokes to physics endpoints;
- adding public API or production route behavior?

Inspect only. Relevant surfaces likely include:

- `src/pqs_source_box_route_driver_helpers.jl`
- `_PQSCompleteCoreShellDiagnosticRoutePayload`
- `_pqs_source_box_route_driver_complete_core_shell_diagnostic_route_payload(...)`
- current complete-core/shell H1/H1J report-field helpers
- current route-driver option/config surfaces if any
- `docs/src/developer/pqs_source_box_fixture_policy.md`

Deliverables:

1. State whether the current driver payload already carries enough structured
   data for private RHF SCF except for explicit diagnostic inputs.
2. Identify the smallest acceptable object/seam for explicit diagnostic inputs,
   e.g. an internal `rhf_diagnostic_request` or `complete_core_shell_rhf_request`
   object, if appropriate.
3. State where that object should be consumed:
   - inside the diagnostic route payload;
   - adjacent to H1/J payload;
   - or nowhere yet.
4. State whether route-driver public options/report fields should remain
   unchanged in the first implementation pass.
5. State whether a cold route dry-run validation would be required for the first
   wiring pass, or whether a synthetic/private helper test is enough initially.
6. List exact blockers and risks.

Decision rules:

- Do not edit files.
- Do not run a heavy route dry-run for this audit.
- If a command is likely to exceed 60 seconds, do not run it for this pass.
- Do not request interactive escalation/approval; if approval is truly needed,
  write `ATTENTION.md` and stop.

Report back:

- Recommended next implementation pass, or a clear block.
- Exact code surfaces to touch next if implementation is recommended.
- Validation status (`git status --short --branch` is enough).
- Deletion/shrinkage forecast:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini

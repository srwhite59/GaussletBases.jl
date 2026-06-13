Pass 129 - no-edit route RHF adoption audit

Purpose:

Audit the smallest private route-driver boundary needed to expose the existing
private complete-core/shell RHF diagnostic through the one-center PQS source-box
driver path. Do not implement it yet.

Why now:

The private SCF helper now has a validated default Fock-DIIS control path for
the compact route-smoke fixture. Before wiring anything into the driver, we need
to pin down:

- the exact private driver slot;
- the compact request object and required fields;
- how to keep route/report/public behavior unchanged unless explicitly
  requested;
- blocker/nonclaim vocabulary.

This remains route-smoke/control diagnostics. Serious HF belongs to
`codexhome/work/hfdmrg`, and downstream CR2 validation belongs to the CR2 agent
after this line is ready.

Read these surfaces:

- `src/pqs_source_box_route_driver_helpers.jl`
  - `_pqs_source_box_route_driver_complete_core_shell_density_inputs`
  - `_pqs_source_box_route_driver_complete_core_shell_h1_j_diagnostic_payload`
  - `_PQSCompleteCoreShellDiagnosticRoutePayload`
  - the `cartesian_assembly` section that consumes the diagnostic route payload
  - existing H1/J report aliases and nonclaim flags
- `src/pqs_multilayer_complete_core_shell_rhf.jl`
  - `_pqs_multilayer_complete_core_shell_rhf_input_contract`
  - `_pqs_multilayer_complete_core_shell_rhf_scf_control_payload`
  - `_pqs_multilayer_complete_core_shell_rhf_scf_payload`
- `docs/src/developer/pqs_source_box_fixture_policy.md`
- `docs/src/developer/pqs_source_box_operator_framework.md`
- If helpful for boundary context only:
  - `/Users/srw/Library/CloudStorage/Dropbox/codexhome/work/hfdmrg/README.md`
  - `/Users/srw/Library/CloudStorage/Dropbox/codexhome/work/cr2/AGENTS.md`

Audit questions:

1. What is the smallest private driver request object or keyword group?
   Recommended shape to evaluate:
   - `private_complete_core_shell_rhf_requested`
   - explicit `electron_count`
   - explicit `fixture_role`
   - optional SCF controls or `scf_control_payload`
   - default `mixing_kind = :fock_diis`

2. Where should the request be consumed?
   Candidate: immediately after the existing complete-core/shell diagnostic
   route payload has `source_payload`, `final_basis`, `h1_payload`,
   `density_inputs`, and `h1_j_payload`.

3. What compact payload/slot should be added if implemented later?
   Candidate names to assess:
   - `_pqs_source_box_route_driver_complete_core_shell_rhf_diagnostic_payload`
   - `complete_core_shell_rhf_diagnostic_payload`
   - summary-only report aliases, if any, should stay minimal.

4. What blockers should exist before materialization?
   Consider:
   - not requested;
   - non-PQS route;
   - missing H1/J diagnostic payload;
   - missing explicit electron count;
   - missing fixture role;
   - unsupported/open-shell electron count;
   - unsupported fixture role;
   - SCF not converged.

5. What fields must remain nonclaims?
   Confirm no route production materialization, no public API, no exports,
   no artifacts, no GTO/IDA/MWG behavior, no physics endpoint, and no serious-HF
   claim.

Trust boundary:

- No file edits.
- No Julia test/probe required unless a tiny load/read-only check is needed.
- No route wiring.
- No report field additions.
- No public API.
- No fixture promotion.
- No `hfdmrg` execution and no CR2 execution.

Decision rule:

- If the seam is clear, recommend the smallest implementation pass and the
  exact files/tests it should touch.
- If the seam would require broad route option plumbing or new report-field
  clouds, recommend a smaller request-object design pass instead.
- If H1/J payload structure is insufficient for RHF input contract, report the
  exact missing object rather than proposing scalar aliases.
- Do not request interactive escalation/approval; if approval is truly needed,
  write `ATTENTION.md` and stop.

Report back:

- Current source inventory with file/function names.
- Recommended request object/keyword shape.
- Recommended private payload/slot name.
- Required blockers/nonclaims.
- Smallest next implementation pass, or why not.
- Git status.
- Deletion/shrinkage forecast:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini

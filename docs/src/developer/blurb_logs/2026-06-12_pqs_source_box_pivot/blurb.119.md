Pass 119 - local compact PQS probe with residual-aware SCF payload

Baseline:

- Current pushed HEAD should include `5957c694 Add PQS RHF SCF residual diagnostics`.
- Private SCF payload now reports compact residual diagnostics in its summary.

Task:

Run a local ignored compact route-owned PQS RHF probe using the updated
residual-aware private SCF payload.

Purpose:

Confirm that the real compact route-smoke fixture now reports residual
diagnostics directly through `_pqs_multilayer_complete_core_shell_rhf_scf_payload(...)`,
without relying on a separate local residual script.

Scope:

- Local ignored `tmp/work` probe only.
- No tracked source/test/doc edits.
- No damping/DIIS/acceleration implementation.
- No route-driver wiring, report fields, public API, exports, artifacts, GTO,
  IDA/MWG, or fixture promotion.

Use:

- explicit `electron_count = 4`
- explicit `fixture_role = :route_smoke`
- current private RHF input contract, initial density, and SCF payload
- `max_iterations = 50` unless you have a good reason to use 25 or 100

Report:

- elapsed time;
- SCF status/blocker;
- iteration count;
- last/final total energy;
- density change and energy change;
- residual diagnostics from `scf.summary.residual_diagnostics`:
  - residual status/blocker;
  - density trace error;
  - closed-shell idempotency error;
  - commutator residual;
  - spatial commutator residual;
  - labels/rules;
- nonclaim flags.

Decision rules:

- If the residual-aware payload does not expose the expected diagnostics on the
  nonconverged real fixture, report the blocker rather than patching.
- If it does expose them and they match the local pass-117 scale, recommend the
  next pass as a design-only acceleration/mixing contract for the private SCF
  helper.
- Do not request interactive escalation/approval; if approval is truly needed,
  write `ATTENTION.md` and stop.

Validation/status:

- Run only the local probe and `git status --short --branch`.

Report back:

- Artifact paths.
- Probe results.
- Recommended next pass.
- Deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini

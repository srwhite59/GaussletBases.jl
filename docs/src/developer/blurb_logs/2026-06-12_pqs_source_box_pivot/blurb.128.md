Pass 128 - confirm compact route-smoke uses default Fock-DIIS history 8, no code

Purpose:

Confirm that after `460d428e`, the compact route-smoke RHF probe converges when
the caller requests `mixing_kind = :fock_diis` but omits `max_history`, so the
new private default history 8 is actually used end to end.

Why now:

Passes 125 and 126 showed explicit `max_history = 8` converges the compact
route-smoke fixture under strict gates. Pass 127 changed the private default
from 6 to 8. Before considering any route-facing RHF slot, confirm the default
path itself produces the same behavior.

Scope:

- Local ignored `tmp/work` probe only.
- No tracked source/test/doc edits.
- No tolerance changes.
- No route wiring/report fields/public API/exports/artifacts/GTO/IDA/MWG.
- No fixture promotion.
- No `hfdmrg` or CR2 comparison; serious HF and downstream Cr2 validation stay
  outside this private PQS route-smoke control probe.

Use:

- explicit `electron_count = 4`
- explicit `fixture_role = :route_smoke`
- tracked `_pqs_multilayer_complete_core_shell_rhf_scf_payload(...)`
- tracked `_pqs_multilayer_complete_core_shell_rhf_scf_control_payload(...)`
- `mixing_kind = :fock_diis`
- omit `max_history`
- assert/report that the control payload resolves `max_history == 8`
- keep `diis_regularization = 1.0e-12`
- keep `diis_coefficient_max_abs = 25.0`
- `max_iterations = 100`
- strict gates: `density_atol = residual_atol = 1.0e-8`,
  `energy_atol = 1.0e-10`

Record:

- artifact paths for the local probe script, summary, and table;
- elapsed time;
- resolved control `max_history`;
- converged?;
- iteration count;
- final/last total energy;
- density change;
- update density change;
- energy change;
- commutator residual;
- spatial commutator residual;
- trace error;
- idempotency error;
- DIIS used/fallback/solve/pathology counts;
- fallback iteration numbers if available without broad instrumentation;
- whether convergence was blocked by final recomputed diagnostics.

Decision rules:

- If the omitted-history default run converges like the explicit history-8
  runs, recommend a no-edit route-RHF adoption audit next: identify the exact
  private driver slot and request object needed, without implementation.
- If it does not converge or resolves the wrong history, stop and report the
  mismatch; do not patch in this pass.
- Do not request interactive escalation/approval; if approval is truly needed,
  write `ATTENTION.md` and stop.

Validation/status:

- Run only the local ignored probe and `git status --short --branch`.

Report back:

- Probe table/summary.
- Recommended next pass.
- Git status.
- Deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini

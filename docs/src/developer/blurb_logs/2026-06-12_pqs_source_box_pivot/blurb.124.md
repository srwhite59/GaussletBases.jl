Pass 124 - validate tracked Fock-DIIS on compact route-smoke fixture

Baseline:

- Current pushed HEAD should include `718071ba Add PQS RHF Fock DIIS controls`.
- Private SCF helper now has tracked bounded Fock-DIIS controls.

Task:

Run a local ignored compact route-owned PQS RHF probe using the tracked private
Fock-DIIS path.

Scope:

- Local ignored `tmp/work` probe only.
- No tracked source/test/doc edits.
- No route wiring/report fields/public API/exports/artifacts/GTO/IDA/MWG.
- No fixture promotion.

Use:

- explicit `electron_count = 4`
- explicit `fixture_role = :route_smoke`
- `_pqs_multilayer_complete_core_shell_rhf_scf_control_payload(...)` with:
  - `mixing_kind = :fock_diis`
  - `max_history = 6`
  - `diis_start_iteration = 2`
  - `diis_regularization = 1.0e-12`
  - `diis_coefficient_max_abs = 25.0`
  - `max_iterations = 100`
  - `density_atol = 1.0e-8`
  - `energy_atol = 1.0e-10`
  - `residual_atol = 1.0e-8`

Report:

- elapsed time;
- SCF status/blocker;
- iteration count;
- final total energy if materialized, otherwise last total energy;
- density change / update density change;
- energy change;
- residual diagnostics:
  - commutator residual;
  - spatial commutator residual;
  - trace error;
  - idempotency error;
- DIIS counters:
  - used count;
  - fallback count;
  - solve failure count;
  - coefficient pathology count;
- whether final convergence was blocked by final recomputed diagnostics.

Decision rules:

- If tracked Fock-DIIS converges on the compact route-smoke fixture, report it
  clearly but do not promote to physics acceptance.
- If it reduces residual but still misses tolerance, report the blocker and
  recommend whether to adjust tolerances, max iterations, history, or design.
- Do not patch production code in this pass.
- Do not request interactive escalation/approval; if approval is truly needed,
  write `ATTENTION.md` and stop.

Validation/status:

- Run only the local ignored probe and `git status --short --branch`.

Report back:

- Artifact paths.
- Probe result table/summary.
- Recommended next pass.
- Deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini

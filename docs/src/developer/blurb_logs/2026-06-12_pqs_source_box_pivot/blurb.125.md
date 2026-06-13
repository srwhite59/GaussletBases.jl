Pass 125 - local tracked-Fock-DIIS asymptote probe, no code

Baseline:

- Current pushed HEAD should include `775b0096 Record PQS tracked DIIS probe`.
- Pass 124 validated the tracked private Fock-DIIS path but it remained blocked
  at strict `1e-8` density/residual gates after 100 iterations.

Task:

Run a local ignored asymptote/parameter probe using the tracked private
Fock-DIIS path on the compact route-owned PQS route-smoke fixture.

Scope:

- Local ignored `tmp/work` probe only.
- No tracked source/test/doc edits.
- No route wiring/report fields/public API/exports/artifacts/GTO/IDA/MWG.
- No fixture promotion.
- No production tolerance/default changes.

Use:

- explicit `electron_count = 4`
- explicit `fixture_role = :route_smoke`
- tracked `_pqs_multilayer_complete_core_shell_rhf_scf_payload(...)`
- tracked `_pqs_multilayer_complete_core_shell_rhf_scf_control_payload(...)`

Probe cases:

1. Main asymptote:
   - `mixing_kind = :fock_diis`
   - `max_history = 6`
   - `diis_regularization = 1.0e-12`
   - `diis_coefficient_max_abs = 25.0`
   - `max_iterations = 200`
   - strict gates: `density_atol = residual_atol = 1.0e-8`,
     `energy_atol = 1.0e-10`

2. Optional comparison if cheap after fixture build:
   - same controls but `max_history = 8`
   - same strict gates

Report for each case:

- converged?
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
- whether convergence was blocked by final recomputed diagnostics.

Decision rules:

- If history 6 converges by 200 iterations, recommend a small tracked control
  default/max-iteration follow-up.
- If history 6 plateaus above `1e-8` but history 8 converges cleanly, recommend
  a local confirmatory run before changing defaults.
- If both plateau near `2e-8` to `3e-8`, recommend a design/tolerance policy
  pass, not route wiring.
- Do not patch code in this pass.
- Do not request interactive escalation/approval; if approval is truly needed,
  write `ATTENTION.md` and stop.

Validation/status:

- Run only the local ignored probe and `git status --short --branch`.

Report back:

- Artifact paths.
- Probe table/summary.
- Recommended next pass.
- Deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini

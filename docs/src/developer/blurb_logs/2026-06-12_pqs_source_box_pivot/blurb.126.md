Pass 126 - confirm tracked Fock-DIIS history 8, no code

Purpose:

Confirm that the tracked private Fock-DIIS path converges the compact
route-smoke fixture with `max_history = 8` under the strict gates before we
consider changing any private default.

Why now:

Pass 125 showed:

- `max_history = 6` plateaus by 200 iterations with density delta around
  `1.22e-7` and commutator residual around `3.12e-8`.
- `max_history = 8` converges by iteration 34 with density delta
  `8.33e-9` and commutator residual `2.95e-9`.
- The history-8 run had 3 coefficient-pathology fallbacks but no solve
  failures.

This is promising but should be confirmed once before changing tracked
defaults. The user also reminded us that serious HF belongs to `hfdmrg`; this
PQS RHF seam remains private route-smoke/control diagnostics only.

Scope:

- Local ignored `tmp/work` probe only.
- No tracked source/test/doc edits.
- No private default changes.
- No tolerance changes.
- No route wiring/report fields/public API/exports/artifacts/GTO/IDA/MWG.
- No fixture promotion and no serious-HF comparison against `hfdmrg`.

Use:

- explicit `electron_count = 4`
- explicit `fixture_role = :route_smoke`
- tracked `_pqs_multilayer_complete_core_shell_rhf_scf_payload(...)`
- tracked `_pqs_multilayer_complete_core_shell_rhf_scf_control_payload(...)`
- `mixing_kind = :fock_diis`
- `max_history = 8`
- `diis_regularization = 1.0e-12`
- `diis_coefficient_max_abs = 25.0`
- `max_iterations = 100`
- strict gates: `density_atol = residual_atol = 1.0e-8`,
  `energy_atol = 1.0e-10`

Record:

- artifact paths for the local probe script, summary, and table;
- elapsed time;
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

- If history 8 again converges under the strict gates with bounded fallbacks,
  recommend a small tracked follow-up changing only the private Fock-DIIS
  default `max_history` from 6 to 8, with focused synthetic validation.
- If history 8 does not converge or shows unstable fallback/pathology behavior,
  recommend a control-policy design pass rather than changing defaults.
- Do not patch code in this pass.
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

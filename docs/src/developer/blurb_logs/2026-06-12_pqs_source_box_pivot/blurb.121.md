Pass 121 - local compact PQS Fock-DIIS experiment, no production code

Baseline:

- Current pushed HEAD should include `06c8391c Record PQS RHF control design`.
- Pass 120 recommends local-only Fock DIIS before adding any private
  SCF-control implementation.

Task:

Run a local ignored Fock-DIIS experiment on the compact route-owned PQS RHF
route-smoke fixture.

Scope:

- Local ignored `tmp/work` scripts/summaries only.
- No tracked source/test/doc edits.
- No production/private DIIS code.
- No route wiring, report fields, public API, exports, artifacts, GTO, IDA/MWG,
  or fixture promotion.

Experiment:

- Use explicit `electron_count = 4`.
- Use explicit `fixture_role = :route_smoke`.
- Reuse the current private RHF input contract, one-step payload, and residual
  conventions.
- Implement Fock DIIS in the local script:
  - build ordinary final-basis Fock from the current density;
  - compute error vector `vec(F * D - D * F)` where `D = P / occupancy`;
  - keep a bounded history of symmetrized Fock matrices and error vectors;
  - solve the standard constrained DIIS least-squares system with a small
    regularization if needed;
  - diagonalize the DIIS-mixed Fock to rebuild idempotent closed-shell `P`;
  - do not mix density matrices directly.

Suggested sweep:

- `max_history = 4` and `6`;
- `diis_start_iteration = 2`;
- `diis_regularization = 1.0e-12` and, only if needed, `1.0e-10`;
- `max_iterations = 100`;
- convergence tolerances matching the current private payload:
  - `density_atol = 1.0e-8`
  - `energy_atol = 1.0e-10`
  - use/report commutator residual against a comparable threshold, e.g.
    `residual_atol = 1.0e-8`.

Report for each run:

- max history and regularization;
- converged?
- iteration count;
- final/last total energy;
- fixed-point density delta;
- energy delta;
- commutator residual;
- trace error;
- idempotency error;
- any DIIS solve failures or coefficient pathologies.

Decision rules:

- If Fock DIIS reduces the commutator residual substantially and/or converges,
  recommend a bounded private SCF-control implementation pass.
- If Fock DIIS fails or is unstable, do not implement it; recommend the next
  audit/experiment.
- If the local script hits an algebraic issue such as singular DIIS matrix,
  report it with the regularization/history that failed.
- Do not request interactive escalation/approval; if approval is truly needed,
  write `ATTENTION.md` and stop.

Validation/status:

- Run only the local ignored experiment and `git status --short --branch`.

Report back:

- Artifact paths.
- Sweep table.
- Best run and rationale, if any.
- Recommended next pass.
- Deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini

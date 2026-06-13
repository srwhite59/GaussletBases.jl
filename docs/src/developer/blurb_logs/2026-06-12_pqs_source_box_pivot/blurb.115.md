Pass 115 - local compact PQS RHF damping sweep, no production code

Baseline:

- Current pushed HEAD should include `69d533ec Record PQS RHF iteration trace`.
- Pass 114 showed monotone energy decrease but density residual plateau/drift
  for the undamped private SCF loop.

Task:

Run a local ignored damping/mixing sweep for the same compact route-owned PQS
RHF route-smoke fixture.

This is a probe only. Do not change tracked production code/tests/docs.

Goal:

Determine whether simple density damping fixes the compact route-smoke RHF SCF
nonconvergence before we design any production/private SCF-control seam.

Allowed:

- Create/update ignored `tmp/work` scripts and summaries.
- Reuse the existing compact route fixture and private one-step payloads.
- Implement damping inside the local probe script only.

Suggested local mixing rule:

```text
P_candidate = 2 * C_occ * C_occ'
P_next = (1 - alpha) * P_current + alpha * P_candidate
```

where `alpha` is the damping/mixing factor. Keep density symmetrized. Track
trace drift; if trace drifts materially, report it rather than silently
renormalizing unless the probe explicitly labels that choice.

Suggested sweep:

- `alpha = 1.0` as undamped baseline if cheap/reused;
- `alpha = 0.75`;
- `alpha = 0.5`;
- `alpha = 0.25`;
- optionally `alpha = 0.1` if the first three still plateau.

Use the same explicit diagnostic inputs:

- `electron_count = 4`
- `fixture_role = :route_smoke`

Recommended limits:

- `max_iterations = 100`
- `density_atol = 1.0e-8`
- `energy_atol = 1.0e-10`

Report for each alpha:

- converged?
- iteration count;
- final/last total energy;
- final/last density change;
- final/last energy change;
- density trace;
- monotone energy yes/no;
- any obvious oscillation/two-cycle.

Decision rules:

- Do not add damping/mixing to production code in this pass.
- If no alpha converges, recommend a residual/Fock convention audit.
- If one or more alpha values converge, recommend a small private SCF-control
  implementation pass with explicit damping controls and tests.
- If the cold route build dominates runtime, report elapsed time. Do not use
  `/usr/bin/time`.
- Do not request interactive escalation/approval; if approval is truly needed,
  write `ATTENTION.md` and stop.

Validation/status:

- Run only the local ignored probe and `git status --short --branch`.

Report back:

- Sweep table.
- Best alpha and rationale, if any.
- Recommended next pass.
- Local artifact paths.
- Deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini

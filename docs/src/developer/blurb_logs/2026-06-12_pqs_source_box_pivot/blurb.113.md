Pass 113 - local ignored compact PQS RHF SCF probe, no tracked code

Baseline:

- Current pushed HEAD should include `f614cd5f Align PQS RHF SCF final payload`.
- The private SCF payload now returns a final one-step payload recomputed on the
  returned final density.

Task:

Run the first local real compact PQS RHF SCF probe using the existing compact
driver/source-plan/H1/J diagnostic ingredients.

This is a probe only:

- local ignored `tmp/work` script is allowed;
- no tracked source/test/doc edits;
- no route-driver wiring;
- no report fields;
- no public API;
- no fixture promotion.

Goal:

Determine whether the private SCF payload works on the actual compact
route-owned PQS diagnostic fixture, beyond the tiny synthetic diagonal test.

Probe requirements:

- Build or reuse the compact one-center PQS source-box diagnostic path used in
  recent H1/J dry-run checks.
- Use explicit diagnostic inputs:
  - `electron_count = 4`
  - `fixture_role = :route_smoke`
- Use the private RHF input contract, initial-density payload, one-step payload,
  and SCF payload.
- Do not infer electron count from charges.
- Do not add a tracked test.

Report:

- Whether the probe reached RHF input contract availability.
- Whether initial density materialized.
- Whether SCF converged.
- If converged:
  - iteration count;
  - final total energy;
  - one-body energy;
  - two-body energy;
  - density trace;
  - final density change;
  - final energy change;
  - final one-step density match error.
- If not converged:
  - blocker/status;
  - iteration count;
  - last total energy if available;
  - last density change / energy change if available.
- Relation to existing H1/J diagnostic:
  - H1 energy;
  - H1/J self-Coulomb if available;
  - clarify that this is route-smoke diagnostic output, not physics acceptance.
- Runtime:
  - Use Julia-level `@elapsed` or `time()` inside the script.
  - If cold startup/build exceeds 60 seconds, that is acceptable for this
    probe, but report it explicitly.
  - Do not use `/usr/bin/time`.

Decision rules:

- If constructing the real compact route fixture is ambiguous, stop and report
  the exact missing call/field rather than inventing a new route.
- If SCF oscillates or fails, do not add damping/mixing in this pass.
- If the probe requires interactive approval/escalation, write `ATTENTION.md`
  and stop.

Validation/status:

- Run only the local probe and `git status --short --branch`.
- `git diff --check` is optional if no tracked files changed.

Deletion/shrinkage report:

- deleted:
- simplified:
- quarantined:
- not deleted because:
- exact remaining caller/blocker:

-- repo-manager@macmini

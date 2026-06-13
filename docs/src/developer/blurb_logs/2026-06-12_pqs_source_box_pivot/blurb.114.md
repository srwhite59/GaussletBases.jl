Pass 114 - local compact PQS RHF iteration trace, no code

Baseline:

- Current pushed HEAD should include `f62ce9e0 Record PQS compact RHF probe`.
- Pass 113 local probe reached private RHF SCF but blocked on
  `:scf_not_converged` after 25 iterations.
- Do not implement damping/mixing yet.

Task:

Run a local ignored iteration-trace probe for the same compact route-smoke PQS
RHF fixture.

Goal:

Collect enough iteration behavior to decide whether the failure looks like:

- slow monotone convergence;
- small oscillation / two-cycle;
- unstable divergence;
- or a bug in the SCF update/energy bookkeeping.

Allowed:

- Create or update ignored `tmp/work` probe artifacts.
- Reuse the pass-113 local probe script if convenient.
- Run one compact route-owned probe, even if cold elapsed time exceeds 60
  seconds; report elapsed time.

Do not:

- Change tracked source/tests/docs.
- Add damping/mixing.
- Change SCF controls in production code.
- Wire the route driver.
- Add report fields, public API, exports, artifacts, GTO, IDA/MWG, or fixture
  promotion.

Probe details:

- Same explicit diagnostic inputs:
  - `electron_count = 4`
  - `fixture_role = :route_smoke`
- Prefer `max_iterations = 50` if runtime is reasonable.
- Capture compact per-iteration fields from `scf.iteration_records`:
  - iteration
  - total_energy
  - density_change
  - energy_change
  - density_converged
  - energy_converged
  - converged
- Also report:
  - min/max/last density change;
  - min/max/last energy change excluding `nothing`;
  - whether total energy is monotone decreasing;
  - whether density change is generally decreasing;
  - any obvious two-cycle signal from the last several density/energy values.

Decision rules:

- If the existing pass-113 artifacts already contain enough iteration records,
  parse them without rerunning the cold route probe.
- If not, rerun a local ignored probe and write a trace summary file under
  `tmp/work/`.
- If the probe fails before producing iteration records, report the blocker
  exactly.
- Do not request interactive escalation/approval; if approval is truly needed,
  write `ATTENTION.md` and stop.

Validation/status:

- Run only the local probe/trace script and `git status --short --branch`.
- `git diff --check` is optional if no tracked files changed.

Report back:

- Whether this was a rerun or parsed existing artifacts.
- Elapsed time if rerun.
- SCF status/blocker.
- Iteration count.
- Trace summary and qualitative diagnosis.
- Recommended next pass:
  - controlled damping/mixing design;
  - longer-control probe;
  - or bug fix if a bookkeeping bug is evident.
- Deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini

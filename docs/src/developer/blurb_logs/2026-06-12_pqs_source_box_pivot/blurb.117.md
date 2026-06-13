Pass 117 - local compact RHF commutator/idempotency residual probe

Baseline:

- Current pushed HEAD should include `91377715 Record PQS RHF convention audit`.
- Pass 116 found no Fock/energy formula bug, but identified missing residual
  characterization.

Task:

Run a local ignored residual probe for the compact route-owned private RHF SCF
fixture.

Goal:

Measure ordinary-final-basis stationarity diagnostics alongside the existing
fixed-point density/energy changes, so we know whether the nonconverged SCF is
actually nonstationary or just failing the current fixed-point delta criterion.

Allowed:

- Create/update ignored `tmp/work` probe scripts and summaries.
- Reuse existing compact route fixture helpers/scripts.
- Rerun the compact route probe if needed; report elapsed time.

Do not:

- Change tracked source/tests/docs.
- Add damping/DIIS/acceleration to production code.
- Wire route driver/report fields.
- Add public API, exports, artifacts, GTO, IDA/MWG, or fixture promotion.

Diagnostics to compute:

For at least the undamped final state from a local compact run, and preferably
for one damped state such as `alpha = 0.25` or `alpha = 0.1` if convenient:

- fixed-point density change;
- energy change;
- total energy;
- density trace and trace error versus electron count;
- closed-shell idempotency error for `D = P / 2`, e.g.
  `norm(D * D - D, Inf)`;
- ordinary final-basis commutator residual, e.g.
  `norm(F * P - P * F, Inf)` or `norm(F * D - D * F, Inf)`;
- whether the occupied subspace from Fock diagonalization differs materially
  from the density subspace.

Use the same explicit diagnostic inputs:

- `electron_count = 4`
- `fixture_role = :route_smoke`

Decision rules:

- If commutator/idempotency residuals are small while fixed-point density delta
  is plateaued, recommend tightening labels/criteria rather than adding
  damping.
- If commutator residual is also large, recommend a controlled SCF-control
  design pass, likely acceleration/mixing beyond scalar damping.
- If residual computation reveals a bookkeeping bug, report the exact corrective
  task and do not patch in this pass.
- Do not request interactive escalation/approval; if approval is truly needed,
  write `ATTENTION.md` and stop.

Validation/status:

- Run only the local residual probe and `git status --short --branch`.

Report back:

- Artifact paths.
- Elapsed time if rerun.
- Residual table/summary.
- Qualitative diagnosis.
- Recommended next pass.
- Deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini

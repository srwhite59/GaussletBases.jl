Purpose:
  Probe the explicit PQS final-basis H1 seam after the direct retained
  overlap/kinetic/nuclear kernels, and identify what old oracle/helper pressure
  it can replace.

Context:
  Passes 036 and 037 changed the active retained one-body PQS construction so
  overlap, kinetic, and by-center nuclear blocks can be built directly in
  retained boundary source modes. Before adding a durable H1 test, run a narrow
  readiness probe and cleanup audit. The goal is to avoid adding another test
  unless it replaces or shrinks older oracle pressure.

Exact task:
  1. Build or update an ignored `tmp/work` probe that assembles the current
     explicit PQS final-basis H1 path using:

     - one `5 x 5 x 5` raw product source box;
     - the PQS boundary source-mode retained rule with retained count `98`;
     - direct retained overlap;
     - direct retained kinetic;
     - direct retained centered by-center nuclear;
     - `CartesianFinalBasisRealization` final-basis realization;
     - CPBM final by-center nuclear transfer and final one-electron Hamiltonian.

  2. The active construction path should not call
     `_pqs_current_route_safe_term_matrices(...)`.

  3. Compare the final Hamiltonian or lowest H1 eigenvalue against the available
     shell-support/oracle path if one is already reachable without broad new
     code. If not, report the exact missing oracle input instead of inventing a
     new oracle framework.

  4. Record a concise artifact under `tmp/work`, for example:

     `tmp/work/pqs_direct_retained_final_h1_probe_summary.txt`

     Include:
     - source dims/count and retained/final counts;
     - final overlap identity error;
     - whether retained one-body blocks used direct boundary construction;
     - whether raw source one-body blocks were materialized on the active path;
     - Hamiltonian symmetry/finite checks;
     - H1 eigenvalue if solved;
     - oracle comparison if available;
     - elapsed timings by coarse phase if cheap.

  5. Audit current tests/callers around `_pqs_current_route_safe_term_matrices`
     and the old shell-support projected helper. Report whether a compact H1
     gate could replace or shrink any of that coverage.

Do not:
  - add a permanent test in this pass unless it clearly replaces/shrinks older
    coverage and stays compact;
  - add IDA, density-density, RHF, driver wiring, exports/artifacts, or GTO
    changes;
  - create a broad route framework or new public API;
  - remove oracle helpers without a clear replacement.

Validation:
  - run the probe;
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`;
  - `git diff --check`.

Required response:
  - probe artifact path;
  - source/final dimensions and direct-path flags;
  - H1/eigenvalue/oracle comparison result if available;
  - current old helper/test pressure that could be replaced by a compact H1
    gate;
  - validation run and result;
  - deletion/shrinkage report:
      - whether anything was deleted or shrunk;
      - if no deletion, why this was probe-only;
      - exact proposed old coverage to shrink if the H1 gate is promoted next.

Continue the baton loop after writing `response.038.md`. Do not stop after one
pass unless blocked by a real design decision or a failure that cannot be
resolved without manager input. Do not request UI escalation; write
`.agent_handoffs/ATTENTION.md` if you are blocked.

-- repo-manager@macmini

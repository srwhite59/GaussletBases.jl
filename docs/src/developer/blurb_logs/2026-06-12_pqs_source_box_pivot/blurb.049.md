Purpose:
  Add the final IDA weight seam for the complete `223`-function PQS final
  basis, without density-density or RHF.

Context:
  Pass 048 stabilized the complete one-center H1 gate:

  ```text
  final dimension: 223
  H1:             -0.48047934800387226
  fixed oracle:   -0.48047920531279725
  ```

  The next physics prerequisite is final IDA weights. These must be the
  integrals of the actual final retained functions, not raw source weights and
  not boundary-shell diagnostic weights.

Exact task:
  Implement or probe the narrow final IDA weight path for
  `pqs_complete_core_shell_final_basis`.

  Intended convention:

  ```text
  final_ida_weights[j] =
      sum_support final_coefficients[support, j] * support_integral_weight[support]
  ```

  where support rows follow the complete basis support order:

  ```text
  core support rows first
  surrounding shell support rows second
  ```

Implementation guidance:
  - Add a compact helper only if useful, probably in
    `CartesianFinalBasisRealization`, such as:

    ```text
    pqs_complete_core_shell_final_ida_weights(final_basis, support_weights)
    ```

  - The helper should validate support-weight length, finite values, final
    coefficient shape, and return a small result/summary.
  - Use product PGDG integral weights for the same support states in the probe.
  - Compare against same-geometry fixed-block retained weights if reachable.
    Treat fixed-block weights as oracle only.
  - Report sign, min/max, sum, and any near-zero or negative weights.

Do not:
  - build density-density;
  - run RHF;
  - add GTO, driver wiring, exports, or artifacts;
  - use old fixed-block weights as active route data;
  - add broad tests. Use `tmp/work` first unless this replaces/shrinks old
    coverage.

Probe:
  Create ignored artifacts:

  ```text
  tmp/work/pqs_complete_core_shell_final_ida_weights_probe.jl
  tmp/work/pqs_complete_core_shell_final_ida_weights_probe_summary.txt
  ```

  Include:
  - final dimension;
  - support weight count;
  - final IDA weight count;
  - min/max/sum;
  - positivity/near-zero diagnostics;
  - fixed-block oracle comparison if reachable;
  - nonclaim flags for density-density/RHF.

Validation:
  - focused IDA-weight probe;
  - if source changed, run `test/nested/pqs_direct_retained_final_h1_runtests.jl`;
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`;
  - `git diff --check`.

Required response:
  - files edited and exact helper names;
  - final IDA weight numbers and oracle comparison/blocker;
  - validation run and result;
  - deletion/shrinkage report:
      - what old/fallback/oracle surface became less necessary;
      - what was deleted or simplified, if anything;
      - if nothing was deleted, why not;
      - whether any new test replaces/shrinks older coverage or is genuinely
        new live-contract coverage;
      - remaining stale or duplicate surfaces to retire next.

Continue the baton loop after writing `response.049.md`. Do not stop after one
pass unless blocked by a real design decision or a failure that cannot be
resolved without manager input. Do not request UI escalation; write
`.agent_handoffs/ATTENTION.md` if you are blocked.

-- repo-manager@macmini

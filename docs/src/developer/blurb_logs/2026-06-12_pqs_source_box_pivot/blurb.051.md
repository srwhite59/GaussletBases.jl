Purpose:
  Resolve the PQS density gauge convention before accepting a final
  density-density matrix.

Context:
  Pass 050 showed that signed final-gauge IDA weight division is unusable:

  ```text
  signed final-weight divided J = 18.132786333403647
  raw projected no-division J   = 0.56242162015939745
  fixed-block oracle J          = 0.6397857768997106
  exact H 1s J                  = 0.625
  ```

  The old nested fixed-block pair builders require positive retained weights
  and divide coefficient functions by those positive weights before raw-pair
  contraction. The current complete final basis has signed final function
  integrals after combined Lowdin cleanup, so it is probably the wrong gauge
  for IDA density weights.

Exact task:
  Audit whether the complete core/shell route has a localized/pre-final density
  gauge with positive weights that should own the density-density construction.

Required checks:
  1. Compute support weights in `core_then_shell` order.
  2. Compute weights for:

     ```text
     pre_final_coefficients
     final_coefficients
     ```

     from `pqs_complete_core_shell_final_basis`.

  3. Report min/max/sum/sign counts for both sets.
  4. If `pre_final_coefficients` have finite positive weights, build the
     weight-aware raw pair contraction in that pre-final gauge:

     ```text
     weighted_coefficients = pre_final_coefficients ./ pre_final_weights
     pair_matrix = weighted_coefficients' * raw_pair_numerator * weighted_coefficients
     ```

     Implement this explicitly for the complete support rows or reuse a trusted
     local kernel, but do not use fixed-block pair_sum as active authority.

  5. Compute H1 self-Coulomb J using:

     - pre-final positive-weight density matrix;
     - final signed-weight divided candidate from pass 050;
     - raw projected no-division candidate;
     - same-geometry fixed-block oracle.

  6. Decide whether the correct next convention is:

     - keep density-density in pre-final/localized positive-weight gauge and
       introduce a separate density-interaction object for RHF consumption;
     - transform a pre-final density interaction to final orbital coefficients
       in a specific reviewed way;
     - or block because no positive density gauge is available.

Do not:
  - run RHF;
  - add GTO, driver wiring, exports, or artifacts;
  - accept raw no-division only because it is closer;
  - accept signed final-weight division;
  - add permanent tests.

Probe:
  Create ignored artifacts:

  ```text
  tmp/work/pqs_complete_core_shell_density_gauge_probe.jl
  tmp/work/pqs_complete_core_shell_density_gauge_probe_summary.txt
  ```

  Include:
  - pre-final and final weight diagnostics;
  - J values for each convention;
  - fixed-block oracle J;
  - max matrix deltas where meaningful;
  - recommended next convention or exact blocker;
  - nonclaim flags for RHF/GTO/driver/export/artifact.

Validation:
  - focused density-gauge probe;
  - if source changed, run `test/nested/pqs_direct_retained_final_h1_runtests.jl`;
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`;
  - `git diff --check`.

Required response:
  - files edited and exact helper names, if any;
  - density gauge findings and J table;
  - recommended convention or exact blocker;
  - validation run and result;
  - deletion/shrinkage report:
      - what old/fallback/oracle surface became less necessary;
      - what was deleted or simplified, if anything;
      - if nothing was deleted, why not;
      - whether any new test replaces/shrinks older coverage or is genuinely
        new live-contract coverage;
      - remaining stale or duplicate surfaces to retire next.

Continue the baton loop after writing `response.051.md`. Do not stop after one
pass unless blocked by a real design decision or a failure that cannot be
resolved without manager input. Do not request UI escalation; write
`.agent_handoffs/ATTENTION.md` if you are blocked.

-- repo-manager@macmini

Purpose:
  Find the trusted same-geometry one-body oracle for the complete core/shell
  H1 route, or identify the exact missing comparison surface.

Context:
  Pass 046 showed the bad Z=1 H value is already present in support space:

  ```text
  support generalized H1 = -2.0639248059188007
  final H1               = -2.0638461028784740
  kinetic expectation    =  1.1012644410060735
  Z=1 nuclear            = -3.1651105438845457
  ```

  So the active blocker is the support/product one-body convention, especially
  the nuclear attraction construction. This is not final-basis transfer, charge
  double-counting, or generalized-solve behavior.

Exact task:
  Run a focused same-geometry one-body oracle audit. Do not add new physics
  features.

Required audit:
  1. Try to build a trusted old/fixed-block or existing route-global one-body
     oracle using the exact same basis/geometry as the complete PQS probe:

     ```text
     parent count = 7
     current_box = (1:7, 1:7, 1:7)
     inner_box = (2:6, 2:6, 2:6)
     raw_source_dims = (5, 5, 5)
     center charge = 1.0
     center location = (0,0,0)
     same Coulomb expansion
     ```

     Do not use the default White-Lindsey seed if it has a different mapping.
     If same-geometry oracle construction is impossible, report the exact
     missing adapter or constructor.

  2. Compare nuclear matrix construction variants on the same support rows:

     - current raw `gaussian_factor_matrices(base_layer)` product construction;
     - any available `pgdg_intermediate.gaussian_factor_terms` convention;
     - any reachable CCPM/PQS source-box nuclear helper convention, such as the
       `_source_box_nuclear_attraction_by_center` family or current CPBM PQS
       centered nuclear helper, without making them active route authority.

  3. For each reachable variant, report:

     ```text
     support generalized H1
     final H1 if cheap
     kinetic expectation
     nuclear expectation
     max matrix delta versus current construction
     sign/charge convention
     ```

  4. Decide whether the issue is:

     - a small local bug in the probe/operator construction;
     - use of the wrong Gaussian-factor convention;
     - a missing same-geometry source-box nuclear operator kernel;
     - or a broader PQS one-body convention gap.

Implementation rule:
  - If a small local fix is clear, implement it.
  - If not, do not patch around the bad H1 value. Record the exact blocker and
    keep the complete H1 route non-acceptance.

Do not:
  - add IDA weights, density-density, RHF, GTO, driver wiring, exports, or
    artifacts;
  - use old fixed-block matrices as active authority;
  - use `_pqs_current_route_safe_term_matrices(...)` as the active route;
  - add a generalized final solve acceptance path;
  - add permanent tests or broaden CPBM contract tests.

Probe:
  Update or create ignored artifacts:

  ```text
  tmp/work/pqs_complete_core_shell_one_body_oracle_probe.jl
  tmp/work/pqs_complete_core_shell_one_body_oracle_probe_summary.txt
  ```

Validation:
  - focused oracle/convention probe;
  - if source changed, rerun the complete H1 probe;
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`;
  - `git diff --check`.

Required response:
  - files edited and exact function/object names changed;
  - oracle(s) reached or exact oracle blocker;
  - convention comparison table;
  - whether a fix was made;
  - validation run and result;
  - deletion/shrinkage report:
      - what old/fallback/oracle surface became less necessary;
      - what was deleted or simplified, if anything;
      - if nothing was deleted, why not;
      - whether any new test replaces/shrinks older coverage or is genuinely
        new live-contract coverage;
      - remaining stale or duplicate surfaces to retire next.

Continue the baton loop after writing `response.047.md`. Do not stop after one
pass unless blocked by a real design decision or a failure that cannot be
resolved without manager input. Do not request UI escalation; write
`.agent_handoffs/ATTENTION.md` if you are blocked.

-- repo-manager@macmini

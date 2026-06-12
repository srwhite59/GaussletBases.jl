Purpose:
  Add the first complete PQS final-basis density-density diagnostic: H1
  one-orbital self-Coulomb J, before RHF.

Context:
  The complete PQS one-center route now has:

  ```text
  final dimension: 223
  H1:             -0.48047934800387226
  final IDA weights materialized
  final IDA weights match same-geometry fixed oracle after gauge alignment
  ```

  The density-density convention must follow the same rule that fixed the WL
  route:

  ```text
  project raw numerator first
  then divide by final_ida_weights[i] * final_ida_weights[j]
  ```

Exact task:
  Implement or probe a narrow final-basis density-density matrix for the
  complete core/shell final basis, then compute the H1 one-orbital
  self-Coulomb diagnostic J.

  Use:
  - complete final basis from `pqs_complete_core_shell_final_basis`;
  - final IDA weights from `pqs_complete_core_shell_final_ida_weights`;
  - raw pair numerator terms from the same PGDG bundle convention, not
    density-normalized pair terms as final authority.

Expected diagnostic:
  For exact hydrogen 1s, the self-Coulomb reference is:

  ```text
  J = 5/8 = 0.625
  ```

  The small `223` basis will not be exact, but the value should be physically
  sane and should be compared with any same-geometry fixed-block oracle if
  reachable.

Implementation guidance:
  - A compact helper may be added if useful, for example:

    ```text
    pqs_complete_core_shell_final_density_density_matrix(...)
    ```

  - It should clearly report:
    - raw numerator projected;
    - final IDA weight division applied;
    - final weights came from `pqs_complete_core_shell_final_ida_weights`;
    - density-density materialized;
    - RHF not materialized.
  - The H1 self-Coulomb contraction should use the same one-orbital convention
    as the decomposed WL He H1/J diagnostic. If that convention is not
    immediately reusable, report the exact contraction used.
  - Check symmetry, finite values, and positive J for the H1 orbital.

Do not:
  - run RHF;
  - add GTO, driver wiring, exports, or artifacts;
  - use fixed-block pair data as active authority;
  - use raw source weights or boundary diagnostic weights as final weights;
  - add broad tests. Use `tmp/work` unless replacing older coverage.

Probe:
  Create ignored artifacts:

  ```text
  tmp/work/pqs_complete_core_shell_final_density_j_probe.jl
  tmp/work/pqs_complete_core_shell_final_density_j_probe_summary.txt
  ```

  Include:
  - final dimension;
  - H1 energy;
  - final IDA weight diagnostics;
  - density-density symmetry/finite checks;
  - H1 self-Coulomb J;
  - exact `0.625` reference error;
  - fixed-block oracle comparison if reachable;
  - nonclaim flags for RHF/GTO/driver/export/artifact.

Validation:
  - focused density/J probe;
  - if source changed, run `test/nested/pqs_direct_retained_final_h1_runtests.jl`;
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`;
  - `git diff --check`.

Required response:
  - files edited and exact helper names;
  - J diagnostic and comparison/reference error;
  - validation run and result;
  - deletion/shrinkage report:
      - what old/fallback/oracle surface became less necessary;
      - what was deleted or simplified, if anything;
      - if nothing was deleted, why not;
      - whether any new test replaces/shrinks older coverage or is genuinely
        new live-contract coverage;
      - remaining stale or duplicate surfaces to retire next.

Continue the baton loop after writing `response.050.md`. Do not stop after one
pass unless blocked by a real design decision or a failure that cannot be
resolved without manager input. Do not request UI escalation; write
`.agent_handoffs/ATTENTION.md` if you are blocked.

-- repo-manager@macmini

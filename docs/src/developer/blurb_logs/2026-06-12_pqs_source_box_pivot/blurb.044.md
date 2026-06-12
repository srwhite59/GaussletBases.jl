Purpose:
  Build or expose the route-owned complete-shell final-basis realization for
  the one-center `5 x 5 x 5` core plus one surrounding shell fixture.

Context:
  Pass 043 showed that the attempted `125 + 98` basis was rank deficient
  because the `98` shell directions were boundary modes inside the same
  `5 x 5 x 5` source space as the 125 core modes. That gives rank `125`, not a
  physical `223`-function basis.

  The intended shell is an independent surrounding layer outside the core. The
  old White-Lindsey seed provides oracle geometry only:

  ```text
  working box:          (1:7, 1:7, 1:7)
  direct core:          (2:6, 2:6, 2:6), count 125
  surrounding support:  7^3 - 5^3 = 218 support points
  retained shell funcs: 98
  total retained dim:   223
  ```

Exact task:
  Implement the narrowest route-owned complete-shell final-basis realization,
  or stop with the exact smaller blocker.

  Preferred route-owned target:

  ```text
  direct core modes from the 5 x 5 x 5 inner core
  + projected/Lowdin-cleaned surrounding shell sector
  -> combined final overlap with rank 223 and identity after cleanup
  ```

Implementation guidance:
  - First focus on final-basis realization only. Do not add H1 operator
    placement unless the independent `223` final basis is already clean.
  - Check whether the projected q-shell machinery can provide the shell sector
    with:

    ```text
    current_box = (1:7, 1:7, 1:7)
    inner_box = (2:6, 2:6, 2:6)
    raw_source_dims = (5, 5, 5)
    ```

    That would match the desired 218 surrounding support points and 98 retained
    shell functions. Use it as the likely route-owned candidate if it works.
  - Keep the old White-Lindsey low-order seed/fixed-block packet only as an
    oracle for geometry, counts, and optional final overlap comparison. Do not
    take active operator matrices from it.
  - Put final-basis realization code in the final-basis realization area if a
    source edit is needed. Do not grow CPBM unless the helper is truly pair/block
    materialization.
  - Prefer a compact object or helper for the combined basis. Avoid a giant
    route report.

Do not:
  - use `_pqs_current_route_safe_term_matrices(...)`;
  - use old fixed-block matrices as active authority;
  - use a generalized-overlap final solve;
  - add overlap/kinetic/nuclear block placement unless the independent final
    basis is already materialized and the extra step is small;
  - add IDA, density-density, RHF, GTO, driver wiring, exports, or artifacts;
  - add broad tests or expand the CPBM contract file.

Probe:
  Update or create ignored probe artifacts:

  ```text
  tmp/work/pqs_complete_core_shell_final_basis_probe.jl
  tmp/work/pqs_complete_core_shell_final_basis_probe_summary.txt
  ```

  Report:
  - core support/count and shell support/count;
  - shell retained count and total final dimension;
  - combined pre-cleanup overlap rank and eigenvalue range;
  - final overlap identity error if materialized;
  - whether the shell support is disjoint from the core support;
  - whether old fixed-block matrices were avoided as active authority;
  - exact blocker if not materialized.

Validation:
  - focused probe;
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`;
  - `git diff --check`.

Required response:
  - files edited and exact function/object names added or changed;
  - whether the independent `223`-function final basis materializes;
  - probe artifact path and key numbers;
  - validation run and result;
  - deletion/shrinkage report:
      - what old/fallback/oracle surface became less necessary;
      - what was deleted or simplified, if anything;
      - if nothing was deleted, why not;
      - whether any new test replaces/shrinks older coverage or is genuinely
        new live-contract coverage;
      - remaining stale or duplicate surfaces to retire next.

Continue the baton loop after writing `response.044.md`. Do not stop after one
pass unless blocked by a real design decision or a failure that cannot be
resolved without manager input. Do not request UI escalation; write
`.agent_handoffs/ATTENTION.md` if you are blocked.

-- repo-manager@macmini

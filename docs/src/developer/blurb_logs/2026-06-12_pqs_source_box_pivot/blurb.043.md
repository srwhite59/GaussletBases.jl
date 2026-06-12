Purpose:
  Add the narrow combined core/shell retained-basis surface needed for the
  complete one-center PQS H1 route.

Context:
  Pass 042 found that a `223`-mode core-plus-shell object already exists only
  as the private White-Lindsey low-order seed/fixed-block packet:

  ```text
  core retained range:   1:125
  shell retained range:  126:223
  total retained dim:    223
  source:                :nested_fixed_block
  status:                :private_development_seed
  ```

  That packet is an oracle/reference only. It must not become the active PQS
  H1 route authority.

Exact target:
  Implement a route-owned combined direct-core plus PQS boundary-shell retained
  one-body assembly for the one-center H / Z=1 H1 probe.

  The needed active block families are:

  ```text
  core-core
  core-shell / shell-core
  shell-shell
  ```

  for:

  ```text
  overlap
  kinetic
  uncharged electron-nuclear by-center
  ```

Implementation guidance:
  - Add a compact combined retained layout/object if needed. It should identify
    core range `1:125`, shell range `126:223`, and final/retained dimension
    `223`.
  - Keep the existing `98`-function boundary-shell final-basis path as the
    shell-sector component, not as a physical endpoint.
  - Build core-core operator blocks directly from source/product-mode factors.
  - Build shell-shell blocks using the direct retained boundary product kernels
    already added.
  - Build core-shell blocks directly from product-mode factors and the boundary
    shell retained-mode list. If this is the hard missing kernel, implement it
    narrowly or stop with the exact function/blocker.
  - Transform/assemble into the combined final-basis H1 path only if the
    route-owned shell projection/Lowdin data can be represented cleanly.
  - The old fixed-block packet may be used only for oracle comparison and shape
    checks, not for active operator matrices.

Do not:
  - use `_pqs_current_route_safe_term_matrices(...)` as the active path;
  - use the private nested fixed-block packet as active matrix authority;
  - use a generalized-overlap final solve;
  - add IDA weights, density-density, RHF, GTO, driver wiring, exports, or
    artifacts;
  - add broad tests or grow the CPBM contract file as a notebook.

Preferred probe:
  Update or create ignored probe artifacts:

  ```text
  tmp/work/pqs_complete_core_shell_h1_probe.jl
  tmp/work/pqs_complete_core_shell_h1_probe_summary.txt
  ```

  If the path runs, report:
  - retained/final dimension, expected `223`;
  - core/shell ranges;
  - final overlap identity error;
  - Hamiltonian finite/symmetric checks;
  - H / Z=1 H1 value and error versus exact `-0.5`;
  - oracle delta versus fixed-block packet if available;
  - flags proving direct route-owned blocks were used and old fixed-block
    matrices were not active.

  If the path does not run, report the smaller blocker precisely, for example:
  - `:missing_core_shell_direct_retained_product_kernel`;
  - `:missing_combined_core_shell_final_basis_realization`;
  - `:missing_route_owned_shell_projection_lowdin_for_combined_basis`.

Validation:
  - focused probe or focused test;
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`;
  - `git diff --check`.

Required response:
  - files edited and exact function/object names added or changed;
  - whether the `223`-function complete H1 route runs, or the smaller blocker;
  - probe artifact path and key numbers;
  - validation run and result;
  - deletion/shrinkage report:
      - what old/fallback/oracle surface became less necessary;
      - what was deleted or simplified, if anything;
      - if nothing was deleted, why not;
      - whether any new test replaces/shrinks older coverage or is genuinely
        new live-contract coverage;
      - remaining stale or duplicate surfaces to retire next.

Continue the baton loop after writing `response.043.md`. Do not stop after one
pass unless blocked by a real design decision or a failure that cannot be
resolved without manager input. Do not request UI escalation; write
`.agent_handoffs/ATTENTION.md` if you are blocked.

-- repo-manager@macmini

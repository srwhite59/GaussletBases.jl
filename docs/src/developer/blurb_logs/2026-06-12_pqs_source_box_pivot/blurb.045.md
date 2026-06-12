Purpose:
  Add the one-center H1 operator placement on top of the new route-owned
  complete core/surrounding-shell final basis.

Context:
  Pass 044 materialized the independent `223`-function final basis:

  ```text
  core support count:        125
  surrounding shell support: 218
  shell retained count:       98
  final dimension:           223
  final overlap identity:    ~7e-14
  ```

  The next missing piece is not final-basis independence anymore. It is
  one-body operator transfer for this combined support/final basis.

Exact target:
  Build the focused H / Z=1 H1 probe using route-owned product-factor operator
  construction and the new complete final-basis realization.

  Required terms:

  ```text
  overlap
  kinetic
  uncharged electron-nuclear by-center
  ```

  Required output:

  ```text
  final overlap
  final kinetic
  final by-center nuclear
  final one-electron Hamiltonian
  ordinary symmetric H1 eigenvalue
  ```

Implementation guidance:
  - Use the complete final-basis object from pass 044, especially its
    `support_row_order`, `core_support_indices`, `shell_support_indices`,
    `pre_final_coefficients`, and `final_coefficients`.
  - Build the operator over the route-owned combined support rows or equivalent
    explicit core-core, core-shell, and shell-shell product blocks from the 1D
    overlap/kinetic/Gaussian factor data.
  - Transform to the final basis with the final coefficient matrix.
  - Preserve the by-center uncharged convention; apply Z=1 charge only in
    Hamiltonian assembly.
  - Use the private fixed-block packet only as an optional oracle comparison,
    not as active operator authority.
  - If a compact helper is useful, put final-basis operator transfer in the
    final-basis realization area or a clearly scoped CPBM helper only if it
    genuinely materializes pair/block operators.

Do not:
  - use `_pqs_current_route_safe_term_matrices(...)`;
  - use old fixed-block matrices as active authority;
  - use a generalized-overlap final solve;
  - add IDA weights, density-density, RHF, GTO, driver wiring, exports, or
    artifacts;
  - add a permanent test unless it replaces/shrinks older oracle pressure;
  - expand broad CPBM contract tests.

Probe:
  Update or create ignored probe artifacts:

  ```text
  tmp/work/pqs_complete_core_shell_h1_probe.jl
  tmp/work/pqs_complete_core_shell_h1_probe_summary.txt
  ```

  Report:
  - final dimension and core/shell ranges;
  - final overlap identity error;
  - Hamiltonian finite/symmetric checks;
  - H / Z=1 H1 value and error versus exact `-0.5`;
  - oracle delta versus fixed-block packet if available;
  - flags proving route-owned product operators were active and fixed-block
    matrices were not active;
  - coarse timings if cheap.

Validation:
  - focused H1 probe;
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`;
  - `git diff --check`.

Required response:
  - files edited and exact function/object names added or changed;
  - whether the `223`-function complete H1 route runs;
  - probe artifact path and key numbers;
  - validation run and result;
  - deletion/shrinkage report:
      - what old/fallback/oracle surface became less necessary;
      - what was deleted or simplified, if anything;
      - if nothing was deleted, why not;
      - whether any new test replaces/shrinks older coverage or is genuinely
        new live-contract coverage;
      - remaining stale or duplicate surfaces to retire next.

Continue the baton loop after writing `response.045.md`. Do not stop after one
pass unless blocked by a real design decision or a failure that cannot be
resolved without manager input. Do not request UI escalation; write
`.agent_handoffs/ATTENTION.md` if you are blocked.

-- repo-manager@macmini

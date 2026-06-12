Purpose:
  Implement the smallest complete one-center PQS/final-basis H1 path for the
  preferred `5 x 5 x 5` inner core plus one surrounding PQS boundary shell.

Context:
  Pass 041 selected the right physical target but found the current production
  surfaces cannot express it:

  ```text
  5 x 5 x 5 inner core:     125 modes
  one PQS boundary shell:    98 modes
  combined complete basis:  223 modes
  exact blocker: :missing_combined_direct_core_boundary_shell_final_basis_realization
  ```

  The previous `98`-function boundary-shell-only H1 path is mechanically useful
  but physically incomplete. It must not become a physical H acceptance target.
  The all-source `125`-mode fallback is diagnostic only, not the preferred
  route.

Exact task:
  Add the narrowest route-owned implementation needed to materialize the
  combined direct-core plus boundary-shell one-center H1 path, or stop with the
  first smaller exact blocker if the implementation cannot be kept narrow.

  The intended dataflow is:

  ```text
  complete one-center source/final basis
  -> combined core/shell retained overlap, kinetic, and by-center nuclear blocks
  -> final-basis transform / Lowdin cleanup
  -> final one-electron Hamiltonian
  -> ordinary symmetric H1 eigensolve
  ```

Implementation constraints:
  - Keep this as one-center H / Z=1 H1 work only.
  - Prefer a compact object or helper for the combined final-basis realization;
    do not add a broad PQS route framework.
  - The combined basis should represent the `125 + 98 = 223` core-plus-shell
    target, not a single all-source-mode fallback.
  - Use route-owned shell projection/Lowdin data if it is already available.
    If it is not available, stop with that exact blocker rather than building a
    loose oracle-only path.
  - Build or place the needed direct retained operator blocks for core-core,
    core-shell, and shell-shell sectors. If one of those block families is the
    real blocker, name it precisely.
  - Keep shell-support projection as oracle/comparison only, not the production
    operator path.

Do not:
  - add final IDA weights, density-density, RHF, GTO, driver wiring, exports, or
    artifacts;
  - use `_pqs_current_route_safe_term_matrices(...)` as the active path;
  - use a generalized-overlap final solve;
  - add a broad permanent test unless it replaces/shrinks older oracle pressure;
  - grow `test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
    as a development notebook.

Probe:
  If runnable, write an ignored probe/artifact such as:

  ```text
  tmp/work/pqs_complete_core_shell_h1_probe.jl
  tmp/work/pqs_complete_core_shell_h1_probe_summary.txt
  ```

  Report:
  - retained/final dimension, expected `223`;
  - physical extent;
  - final overlap identity error;
  - Hamiltonian symmetry/finite checks;
  - H / Z=1 H1 value and exact-reference error versus `-0.5`;
  - shell-support oracle delta if available;
  - direct-path flags and nonclaim flags;
  - coarse timings if cheap.

Validation:
  - run the focused probe or focused test;
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`;
  - `git diff --check`.

Required response:
  - files edited and exact object/function names added or changed;
  - whether the `223`-function complete H1 path runs, or the smaller exact
    blocker if not;
  - probe artifact path and key numbers if runnable;
  - validation run and result;
  - deletion/shrinkage report:
      - what CPBM/CFBR/source helper surface became simpler or more precise;
      - whether any old oracle-only path became less necessary;
      - whether any new test replaces/shrinks older coverage or is genuinely
        new live-contract coverage;
      - remaining stale or duplicate surfaces to retire next.

Continue the baton loop after writing `response.042.md`. Do not stop after one
pass unless blocked by a real design decision or a failure that cannot be
resolved without manager input. Do not request UI escalation; write
`.agent_handoffs/ATTENTION.md` if you are blocked.

-- repo-manager@macmini

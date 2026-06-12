Purpose:

Audit the remaining cold cost inside the new decomposed WL one-electron
matrix-set phase before doing any broader atom+GTO route-result cleanup.

Why now:

The one-electron matrix-set pass improved the Be S+P route:

```text
cold route elapsed    29.012344166s -> 25.465580292s
warm route elapsed     0.575187292s ->  0.495644334s
```

The previous separate overlap/kinetic/nuclear phases are gone from the top-level
timing table, but the new fused phase is now the largest cold route phase:

```text
decomposed_wl_one_electron_matrix_set   cold 4.432942625s   warm 0.054274709s
```

The current implementation in
`src/cartesian_pair_block_materialization/route_global_one_body_adapter.jl`
still builds the matrix set by calling internal report-shaped helpers:

```text
_route_global_decomposed_wl_factorized_one_body_matrix(..., :overlap)
_route_global_decomposed_wl_factorized_one_body_matrix(..., :kinetic)
_route_global_decomposed_wl_factorized_electron_nuclear_by_center_matrix(...)
```

That may still be too report-shaped for the hot atom+GTO path. Measure before
changing it.

Exact task:

1. Audit the subphases inside
   `_route_global_decomposed_wl_factorized_one_electron_matrix_set(...)`.
   Determine whether the cold cost is mainly:
   - overlap helper result construction;
   - kinetic helper result construction;
   - electron-nuclear helper result construction;
   - repeated factorized retained-basis/axis-table setup;
   - first-call compilation in old factorized kernels;
   - or something else.

2. If the audit shows helper report construction or repeated setup is the
   issue, add the smallest inner hot-path cleanup:
   - compute overlap, kinetic, and electron-nuclear matrices directly from one
     shared factorized retained-basis sidecar and shared axis tables where
     practical;
   - return the same `DecomposedWLOneElectronMatrixSet` shape;
   - keep the existing factorized helper reports and public wrappers available
     for direct term-level callers and reference/debug use.

3. Do not change the by-center electron-nuclear convention:
   - centers remain separated;
   - nuclear charges are recorded but not applied;
   - centers are not summed until Hamiltonian assembly.

4. Do not clean up the outer atom+GTO route result in this pass unless the audit
   shows it is the dominant source. The legacy overlap/kinetic/nuclear payload
   slots can remain temporarily awkward.

Trust boundary:

No public API/export/default changes. No PQS, ECP, Be2, H2, high-l Be, driver
defaults, acceptance fixtures, raw GTO final density-density, generalized final
solve, full-parent CPB fallback, direct Cartesian fallback, or ordinary
Cartesian IDA fallback.

Test policy:

Do not add tests by default. Use `tmp/work` timing probes and existing focused
tests. Add a test only if a genuinely new live module contract is introduced
and cannot be covered by an existing focused test or the Be physics probe.

Validation:

- rerun `tmp/work/be_atom_sp_decomposed_final_basis_warm_cold_timing_probe.jl`;
- run `julia --project=. -e 'using GaussletBases; println("load ok")'`;
- run `git diff --check`;
- if source changes touch factorized one-body machinery, run
  `test/nested/cartesian_wl_factorized_backend_equivalence_runtests.jl`;
- if source changes touch combined GTO assembly, run
  `test/nested/cartesian_route_global_combined_gto_layout_runtests.jl`.

Deletion/shrinkage report required:

- what old code, test, metadata, or compatibility path became unnecessary;
- what was deleted or simplified;
- if nothing was deleted, why no existing surface was made obsolete;
- whether any new test replaces/shrinks older coverage or is genuinely new
  live-contract coverage;
- any remaining stale or duplicate surfaces to retire next.

Report back:

- files changed;
- audit result for the one-electron matrix-set subphases;
- before/after cold and warm timing table for the matrix-set phase and total
  route;
- whether any helper report construction was bypassed in the hot path;
- physics comparison to old nested/QW Be S+P oracle;
- validation run;
- deletion/shrinkage report.

-- repo-manager@macmini

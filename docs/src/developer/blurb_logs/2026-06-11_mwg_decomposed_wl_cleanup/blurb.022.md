Purpose:

Audit and, if justified, reduce the next atom+GTO cold compile-pressure source:
the separate decomposed WL one-electron route calls for overlap, kinetic, and
electron-nuclear by center.

Why now:

The residual-moment pass succeeded locally:

```text
residual_moment_matrices cold   6.912517084s -> 1.074904167s
residual_moment_matrices warm   0.00290625s  -> 0.002948792s
```

The largest remaining cold block with the same route-shape smell is the
gausslet one-electron construction in
`src/cartesian_pair_block_materialization/route_global_atom_gto_final_basis_route.jl`.
The active atom+GTO route still calls these as separate report-shaped route
results:

```text
route_global_decomposed_wl_overlap_matrix
route_global_decomposed_wl_kinetic_matrix
route_global_decomposed_wl_electron_nuclear_by_center_matrices
```

Current Be S+P cold timings for those phases:

```text
electron_nuclear_by_center   4.311094166s
overlap                      1.940648875s
kinetic                      1.421223375s
combined                     7.672966416s
```

Exact task:

1. Audit before editing. Confirm whether these cold costs are mainly:
   - three separate decomposed WL one-body route calls;
   - report-shaped one-body result payloads;
   - repeated factorized retained-basis/axis table setup;
   - electron-nuclear-specific first-call compilation;
   - or a different source.

2. If the audit supports it, add the smallest route-owned compute seam for the
   active atom+GTO path. A good shape would be a narrow decomposed WL
   one-electron matrix set carrying:

   ```text
   overlap matrix
   kinetic matrix
   electron-nuclear by-center matrices
   status/blocker/provenance flags
   ```

   This should be a stable compute noun, not a giant route report object.

3. Reuse the existing factorized retained-basis backend and existing
   electron-nuclear by-center convention. Centers must remain separated and
   uncharged until Hamiltonian assembly. Do not change the physics convention.

4. Keep the existing public/compatibility helpers working:

   ```text
   route_global_decomposed_wl_overlap_matrix
   route_global_decomposed_wl_kinetic_matrix
   route_global_decomposed_wl_electron_nuclear_by_center_matrices
   ```

   They can remain wrappers/reference surfaces. The atom+GTO hot route should
   use the lean matrix-set seam only if it is clearly validated.

5. Do not generalize beyond the measured atom+GTO path in this pass. Do not
   touch GTO/GTO, MWG density-density, RHF, PQS, ECP, Be2, H2, driver defaults,
   or acceptance fixtures.

Physics and route contract to preserve:

- Be S+P RHF total matches old nested/QW oracle to roundoff.
- Final-basis solve remains ordinary/final-basis, not generalized-overlap.
- No raw GTO density-density accepted as a final operator.
- No full-parent CPB fallback.
- No direct Cartesian fallback.
- No ordinary Cartesian IDA fallback.

Test policy:

Do not add tests by default. Use the Be `tmp/work` probe and existing focused
tests. Add a test only if a genuinely new live module contract is introduced
and cannot be covered by an existing focused test or the physics probe.

Validation:

- rerun `tmp/work/be_atom_sp_decomposed_final_basis_warm_cold_timing_probe.jl`;
- run `julia --project=. -e 'using GaussletBases; println("load ok")'`;
- run `git diff --check`;
- if source changes touch one-body/factorized machinery, run
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
- audit result for overlap/kinetic/electron-nuclear one-electron phases;
- before/after cold and warm timing table for those phases and total route;
- whether any struct was introduced, and why it is a stable compute concept;
- physics comparison to old nested/QW Be S+P oracle;
- validation run;
- deletion/shrinkage report.

-- repo-manager@macmini

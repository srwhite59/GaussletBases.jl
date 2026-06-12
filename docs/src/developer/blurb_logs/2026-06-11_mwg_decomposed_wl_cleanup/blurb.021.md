Purpose:

Reduce the next measured atom+GTO cold compile-pressure source without turning
the cleanup into a broad framework pass.

Why now:

The inventory-summary cleanup removed the concrete side-7/side-15 tuple shape
and the Be S+P cold route timing improved substantially:

```text
previous cold route     159.006798084s
current cold route       30.843591209s
current warm route        0.492626583s
```

The largest remaining cold phase is now:

```text
residual_moment_matrices    cold 6.912517084s    warm 0.00290625s
```

The current seam is visible in
`src/cartesian_pair_block_materialization/route_global_atom_gto_final_basis_route.jl`.
Inside `_white_lindsey_decomposed_atom_gto_final_basis_route(...)`, the
`residual_moment_matrices` phase calls six separate decomposed WL moment
builders:

```text
route_global_decomposed_wl_position_x_matrix
route_global_decomposed_wl_position_y_matrix
route_global_decomposed_wl_position_z_matrix
route_global_decomposed_wl_x2_x_matrix
route_global_decomposed_wl_x2_y_matrix
route_global_decomposed_wl_x2_z_matrix
```

It then passes those six report-shaped results into
`route_global_combined_gto_residual_moment_matrices(...)` in
`src/cartesian_pair_block_materialization/route_global_combined_gto_matrix_assembly.jl`.

Exact task:

1. Audit the residual-moment phase shape before editing. Confirm whether the
   cold cost is mainly:
   - six separate decomposed one-body route calls;
   - report-shaped moment result payloads;
   - `route_global_combined_gto_residual_moment_matrices(...)` itself;
   - or a different first-call compile source.

2. If the audit confirms the current shape is the issue, add the smallest
   route-owned compute seam for the active atom+GTO path. Good directions:
   - a fused decomposed WL moment-matrix-set builder for position/x/y/z and
     x2/x/y/z that reuses the factorized retained-basis sidecar once;
   - a compact residual-moment compute object such as a position/x2 axis matrix
     triplet, if that makes the hot route carry matrices rather than six
     report-shaped route objects.

3. Do not replace `NamedTuple`s everywhere. Use a struct only for a stable
   compute noun that is passed across phases. Compact summaries may remain
   `NamedTuple`s.

4. Keep the existing public/compatibility APIs working. The existing six
   `route_global_decomposed_wl_*` moment helpers may remain as wrappers or
   reference surfaces. The atom+GTO hot path should use the leaner compute seam
   if it exists.

5. Preserve physics and final-basis rules:
   - Be S+P RHF total must still match the old nested/QW oracle to roundoff;
   - final overlap remains the ordinary final-basis contract;
   - no generalized final solve;
   - no raw GTO density-density accepted as a final operator;
   - no full-parent CPB, direct Cartesian fallback, or ordinary Cartesian IDA
     fallback.

Trust boundary:

No public exports unless already required by local module organization. No
PQS, ECP, Be2, H2, high-l Be, driver default changes, new acceptance fixtures,
new broad metadata layers, raw GTO final density-density, generalized final
solve, full-parent CPB fallback, direct Cartesian fallback, or ordinary
Cartesian IDA fallback.

Test policy:

Do not add tests by default. If source changes need validation, prefer an
existing focused module test plus the Be `tmp/work` probe. Add a new test only
if there is a new live module contract that cannot be covered by an existing
focused test or the physics probe.

Validation:

- rerun `tmp/work/be_atom_sp_decomposed_final_basis_warm_cold_timing_probe.jl`;
- rerun `tmp/work/atom_gto_specialization_shape_audit.jl` if result shapes
  change;
- run `julia --project=. -e 'using GaussletBases; println("load ok")'`;
- run `git diff --check`;
- if source changes are made, run the smallest directly affected focused test.

Deletion/shrinkage report required:

- what old code, test, metadata, or compatibility path became unnecessary;
- what was deleted or simplified;
- if nothing was deleted, why no existing surface was made obsolete;
- whether any new test replaces/shrinks older coverage or is genuinely new
  live-contract coverage;
- any remaining stale or duplicate surfaces to retire next.

Report back:

- files changed;
- audit result for the residual-moment phase;
- before/after cold and warm timing table, especially
  `residual_moment_matrices`;
- whether any struct was introduced, and why it is a stable compute concept
  rather than a giant report object;
- physics comparison to old nested/QW Be S+P oracle;
- validation run;
- deletion/shrinkage report.

-- repo-manager@macmini

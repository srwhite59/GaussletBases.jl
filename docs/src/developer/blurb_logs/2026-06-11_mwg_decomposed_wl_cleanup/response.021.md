Reduced the Be atom+GTO residual-moment cold compile pressure with a narrow decomposed WL moment matrix-set compute seam.

Files changed:

- `src/cartesian_pair_block_materialization/route_global_one_body_adapter.jl`
- `src/cartesian_pair_block_materialization/route_global_combined_gto_matrix_assembly.jl`
- `src/cartesian_pair_block_materialization/route_global_atom_gto_final_basis_route.jl`
- `docs/src/developer/blurb_logs/2026-06-11_mwg_decomposed_wl_cleanup/response.021.md`

Audit result:

- The residual-moment phase was dominated by six separate decomposed WL moment route calls:
  - `route_global_decomposed_wl_position_x_matrix`
  - `route_global_decomposed_wl_position_y_matrix`
  - `route_global_decomposed_wl_position_z_matrix`
  - `route_global_decomposed_wl_x2_x_matrix`
  - `route_global_decomposed_wl_x2_y_matrix`
  - `route_global_decomposed_wl_x2_z_matrix`
- Each helper could use the factorized retained-basis backend, but each call still rebuilt route-shaped one-body result payloads and repeated shared axis-overlap table setup.
- `route_global_combined_gto_residual_moment_matrices(...)` only needs the six gausslet/gausslet matrices, so carrying six full report-shaped one-body results through the atom+GTO hot path was unnecessary.

Implementation:

- Added internal `DecomposedWLMomentMatrixSet`.
- Added `route_global_decomposed_wl_moment_matrix_set(...)` as a route-owned compute seam for the active decomposed WL atom+GTO path.
- The factorized path computes all six retained gausslet moment matrices from one factorized retained-basis sidecar and shared axis-overlap tables:
  - `position_x`, `position_y`, `position_z`
  - `x2_x`, `x2_y`, `x2_z`
- The existing six `route_global_decomposed_wl_*` moment helpers remain available and still work as wrappers/reference surfaces.
- `route_global_combined_gto_residual_moment_matrices(...)` now accepts either:
  - the new `moment_matrix_set`, or
  - the previous six report-shaped result keywords.
- The atom+GTO hot route now uses `route_global_decomposed_wl_moment_matrix_set(...)` and passes the matrix set into combined residual-moment assembly.

Why the struct is justified:

- `DecomposedWLMomentMatrixSet` is a stable compute concept: the six retained gausslet moment matrices needed by residual MWG construction.
- It is not a broad route report object. It carries status/blocker, the inventory reference, the six matrices, provenance flags, and metadata only.
- The existing combined residual-moment result remains the report object consumed by later route stages.

Be S+P timing comparison:

| metric | before | after | change |
| --- | ---: | ---: | ---: |
| cold route elapsed | `30.843591209s` | `29.012344166s` | `-1.831247043s` |
| cold total elapsed | `33.516908959s` | `31.604365541s` | `-1.912543418s` |
| warm route elapsed | `0.492626583s` | `0.575187292s` | `+0.082560709s` |
| warm total elapsed | `2.490025333s` | `2.574732084s` | `+0.084706751s` |
| residual moment cold | `6.912517084s` | `1.074904167s` | `-5.837612917s` |
| residual moment warm | `0.00290625s` | `0.002948792s` | `+0.000042542s` |

The target phase improved substantially: residual moment cold time dropped from about `6.91s` to about `1.07s`. Overall cold route time improved by about `1.83s`; other cold phases in this fresh run were somewhat noisier/slower, so the phase-local improvement is larger than the top-level route delta.

Warm route time was slightly slower in this one measurement, but the residual-moment warm phase itself was unchanged at about `0.003s`. The warm-route increase is not from the new residual-moment seam; the largest warm phase remains `mixed_gto_blocks` at about `0.365s`, and RHF still dominates total same-process time.

Physics comparison to old nested/QW Be S+P oracle:

- `rhf_total_energy = -14.574514244574639`
- old nested/QW oracle total = `-14.574514244574694`
- delta = `5.5067062021407764e-14 Ha`
- final dimension = `636`
- retained gausslet dimension = `615`
- units/pairs = `131 / 8646`
- no full-parent CPB, direct Cartesian fallback, ordinary Cartesian IDA fallback, raw GTO final density-density, or generalized final-basis solve.

Shape audit:

- Reran `tmp/work/atom_gto_specialization_shape_audit.jl`.
- The inventory shape cleanup remains intact:
  - `unit_keys`: `Vector{Symbol}` for side-7 and side-15.
  - `unit_summaries`: `Vector{NamedTuple{...}}` for side-7 and side-15.
  - `pair_summaries`: compact count/status `NamedTuple` for side-7 and side-15.
- The audit script does not materialize the residual-moment route, so it does not report `DecomposedWLMomentMatrixSet` directly.

Validation run:

- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `julia --project=. test/nested/cartesian_route_global_combined_gto_layout_runtests.jl`
- `julia --project=. test/nested/cartesian_wl_factorized_backend_equivalence_runtests.jl`
- `julia --project=. tmp/work/be_atom_sp_decomposed_final_basis_warm_cold_timing_probe.jl`
- `julia --project=. tmp/work/atom_gto_specialization_shape_audit.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Deletion/shrinkage report:

- The atom+GTO hot route no longer stages six full report-shaped decomposed WL moment result objects just to extract their matrices.
- Removed 44 lines of repeated moment-helper calls from `_white_lindsey_decomposed_atom_gto_final_basis_route(...)`.
- The previous six public/local helper surfaces were not deleted because they remain compatibility/reference APIs and are still useful for direct term-level callers.
- No tests were added. Existing focused module tests and the Be physics/performance probe cover the live contract.
- Remaining stale or duplicate surfaces to retire next:
  - residual-moment assembly still returns a report-shaped `NamedTuple`; if it becomes a compile source, split compute matrices from report summary there too.
  - atom+GTO staged route result still mixes compute payloads with audit/report fields.
  - parent-axis bundle objects still ride through metadata instead of a route-owned parent-axis context.
  - one-body wrappers still build rich reports for direct term-level calls; keep them as reference/compatibility unless another hot path uses them.

-- repo-doer@macmini

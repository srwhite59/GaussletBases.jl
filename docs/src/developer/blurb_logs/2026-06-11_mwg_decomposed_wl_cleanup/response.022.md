Reduced the Be atom+GTO one-electron cold compile pressure with a narrow decomposed WL one-electron matrix-set compute seam.

Files changed:

- `src/cartesian_pair_block_materialization/route_global_one_body_adapter.jl`
- `src/cartesian_pair_block_materialization/route_global_combined_gto_matrix_assembly.jl`
- `src/cartesian_pair_block_materialization/route_global_atom_gto_final_basis_route.jl`
- `docs/src/developer/blurb_logs/2026-06-11_mwg_decomposed_wl_cleanup/response.022.md`

Audit result:

- The active atom+GTO route still made three separate decomposed WL one-electron route calls:
  - `route_global_decomposed_wl_overlap_matrix(...)`
  - `route_global_decomposed_wl_kinetic_matrix(...)`
  - `route_global_electron_nuclear_by_center_matrices(...)`
- Those calls already used the factorized backend on the `UnitPairIndexTable` route, but each call still built a report-shaped route result and repeated factorized one-body setup.
- The combined GTO one-electron assembly only needs:
  - gausslet retained overlap matrix;
  - gausslet retained kinetic matrix;
  - separated uncharged gausslet retained electron-nuclear by-center matrices.
- The measured cold smell was therefore the same as the residual-moment pass: hot route carrying report-shaped one-body results where a small compute payload is enough.

Implementation:

- Added internal `DecomposedWLElectronNuclearByCenterMatrix`.
- Added internal `DecomposedWLOneElectronMatrixSet`.
- Added `route_global_decomposed_wl_one_electron_matrix_set(...)`.
- The factorized implementation computes:
  - overlap through `_route_global_decomposed_wl_factorized_one_body_matrix(...)`;
  - kinetic through `_route_global_decomposed_wl_factorized_one_body_matrix(...)`;
  - each separated by-center nuclear matrix through `_route_global_decomposed_wl_factorized_electron_nuclear_by_center_matrix(...)`.
- The old wrapper path remains as a fallback/reference path for non-index-table inventories.
- `route_global_combined_gto_one_electron_matrices(...)` now accepts either:
  - the new `one_electron_matrix_set`, or
  - the previous `overlap_result`, `kinetic_result`, and `electron_nuclear_by_center_results` arguments.
- `_white_lindsey_decomposed_atom_gto_final_basis_route(...)` now uses the matrix set on the hot Be atom+GTO path.

Why the structs are justified:

- `DecomposedWLOneElectronMatrixSet` is a stable compute concept: the gausslet retained one-electron matrices needed by combined atom+GTO assembly.
- `DecomposedWLElectronNuclearByCenterMatrix` is the separated by-center nuclear payload with the existing uncharged convention.
- These are not giant report objects. They carry matrices, status/blocker, provenance/nonfallback flags, and compact metadata only.

Before/after timing:

| metric | before | after | change |
| --- | ---: | ---: | ---: |
| cold route elapsed | `29.012344166s` | `25.465580292s` | `-3.546763874s` |
| cold total elapsed | `31.604365541s` | `28.130448667s` | `-3.473916874s` |
| warm route elapsed | `0.575187292s` | `0.495644334s` | `-0.079542958s` |
| warm total elapsed | `2.574732084s` | `2.494674792s` | `-0.080057292s` |
| old cold overlap + kinetic + nuclear | `7.672966416s` | fused away | n/a |
| new cold one-electron matrix set | n/a | `4.432942625s` | about `-3.240023791s` versus old sum |
| old warm overlap + kinetic + nuclear | about `0.055876957s` | fused away | n/a |
| new warm one-electron matrix set | n/a | `0.054274709s` | effectively unchanged |

Current measured cold phase table for the touched path:

| phase | cold seconds | warm seconds |
| --- | ---: | ---: |
| `decomposed_wl_one_electron_matrix_set` | `4.432942625` | `0.054274709` |
| `combined_one_electron_assembly` | `2.074959666` | `0.009049458` |

Physics comparison to old nested/QW Be S+P oracle:

- `rhf_total_energy = -14.574514244574639`
- old nested/QW oracle total = `-14.574514244574694`
- delta = `5.5067062021407764e-14 Ha`
- final dimension = `636`
- retained gausslet dimension = `615`
- units/pairs = `131 / 8646`
- no full-parent CPB, direct Cartesian fallback, ordinary Cartesian IDA fallback, raw GTO final density-density, or generalized final-basis solve.

Validation run:

- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `julia --project=. test/nested/cartesian_wl_factorized_backend_equivalence_runtests.jl`
- `julia --project=. test/nested/cartesian_route_global_combined_gto_layout_runtests.jl`
- `julia --project=. tmp/work/be_atom_sp_decomposed_final_basis_warm_cold_timing_probe.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Deletion/shrinkage report:

- The atom+GTO hot route no longer stages three full decomposed WL one-electron report objects just to feed combined GTO one-electron assembly.
- Simplified `_white_lindsey_decomposed_atom_gto_final_basis_route(...)` by replacing the separate overlap, kinetic, and nuclear route calls with one matrix-set phase.
- The existing public/compatibility helpers were not deleted because direct term-level callers and tests still use them as route-global one-body APIs.
- No tests were added. Existing focused factorized/combined-GTO tests plus the Be physics/performance probe covered the live contract.
- Remaining stale or duplicate surfaces to retire next:
  - the atom+GTO route result still mixes compute payloads with audit/report fields;
  - the private route still stores the one-electron matrix set in the legacy overlap/kinetic/nuclear payload slots for now;
  - direct one-body helpers still build rich reports, which is fine for compatibility but should not be used in hot atom+GTO paths;
  - parent-axis bundle objects still ride through metadata instead of a route-owned parent-axis context.

-- repo-doer@macmini

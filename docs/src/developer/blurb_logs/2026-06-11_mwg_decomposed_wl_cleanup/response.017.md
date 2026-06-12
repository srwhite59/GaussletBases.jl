Added a narrow production precompile workload for the decomposed atom+GTO final-basis route shape.

Files changed:

- `src/GaussletBases.jl`
- `src/precompile_workloads.jl`
- `docs/src/developer/blurb_logs/2026-06-11_mwg_decomposed_wl_cleanup/response.017.md`

Generated ignored developer artifact:

- `tmp/work/atom_gto_precompile_synthetic_probe.jl`

What changed:

- Moved `include("precompile_workloads.jl")` later in `src/GaussletBases.jl`, after `cartesian_gto_probes.jl`.
- This is required because the atom+GTO precompile route needs the Cartesian Gaussian supplement representation methods and old QW/GTO cross-table helpers that are not available at the previous include point.
- Added a tiny repo-local synthetic supplement workload in `src/precompile_workloads.jl`.
- The workload calls the existing private route seam:
  - `_white_lindsey_decomposed_atom_gto_final_basis_route(...)`
- It does not solve RHF and does not add public API, route defaults, exports, artifacts, fallback paths, or acceptance behavior.

Synthetic fixture used:

- one synthetic Cartesian S Gaussian at the origin
- exponent = `1.0`
- coefficient = `1.0`
- supplement kind = `:precompile_synthetic_atomic_cartesian_shell`
- basis metadata `basisfile = "repo-local-synthetic"`
- one center at the origin with `Z = 2.0`
- mapped basis count = `7`
- `q/ns = 5/5`
- `white_lindsey_atomic_mapping(Z = 2.0, d = 0.2, tail_spacing = 10.0)`
- `build_density_density = true`

Route surfaces exercised by the workload:

- shellification-backed decomposed WL inventory
- route-global overlap
- route-global kinetic
- route-global electron-nuclear by-center
- factorized mixed gausslet/GTO blocks
- combined one-electron matrix assembly
- final-basis projection
- residual moment matrix construction
- residual MWG representation
- gausslet density-density
- final-basis density-density

Package precompile and load measurements:

- package precompile after the edit: `48956.0 ms`
- fresh cached `using GaussletBases`: `0.654013292s`

Be S+P timing before/after:

| measurement | before workload | after workload |
| --- | ---: | ---: |
| cold route elapsed | `169.936205334s` | `159.006798084s` |
| cold RHF elapsed | `2.540070625s` | `2.930222416s` |
| cold total elapsed | `172.476275959s` | `161.9370205s` |
| measured warm route elapsed | `0.523052584s` | `0.498273209s` |
| measured warm RHF elapsed | `1.880134542s` | `1.988049084s` |
| measured warm total elapsed | `2.403187126s` | `2.486322293s` |

Physics stayed fixed:

- post-workload warm RHF total = `-14.574514244574639`
- old nested/QW oracle total = `-14.574514244574694`
- delta from old oracle = `5.5067062021407764e-14 Ha`
- final dimension = `636`
- retained gausslet dimension = `615`
- units/pairs = `131 / 8646`
- no full-parent CPB, direct Cartesian fallback, ordinary Cartesian IDA fallback, raw GTO final density-density, or generalized final solve was used.

Post-workload cold versus warm phase highlights:

- cold total still dominated by compilation-heavy route phases.
- cold route improved by about `10.9s`, from `169.94s` to `159.01s`.
- warm route remained about `0.5s`.
- mixed-GTO warm factorized projection stayed in the same range, about `0.31s`.
- RHF solve remains the largest same-process warm cost, about `2.0s`.

Carrying-cost decision:

- The workload earns a limited first-step carrying cost because it is small, repo-local, and exercises the actual atom+GTO final-basis route without acceptance logic or user-local BasisSets data.
- The benefit is modest, not complete: the synthetic side-7 route does not fully precompile the side-15 Be route specializations, so Be cold latency only drops by about `6%`.
- I did not broaden the workload to the side-15 Be shape because that would make package precompile carry a large route-sized workload and would risk encoding the Be acceptance fixture into package load.
- The next precompile pass, if any, should target a shape-stable way to reduce specialization on retained-unit tuple sizes or large route summary types, rather than simply running a larger synthetic atom in precompile.

Validation run:

- `julia --project=. tmp/work/atom_gto_precompile_synthetic_probe.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `julia --project=. -e 't = @elapsed using GaussletBases; println("load_elapsed_s=", t); println("load ok")'`
- `julia --project=. tmp/work/be_atom_sp_decomposed_final_basis_warm_cold_timing_probe.jl`
- `git diff --check`

Deletion/shrinkage report:

- No old production code, tests, metadata, or compatibility path became unnecessary.
- Nothing was deleted because this pass adds a narrow compile workload and does not replace a runtime route or test surface.
- No tests were added. The synthetic probe is an ignored `tmp/work` developer probe and is not long-term coverage.
- Remaining stale/duplicate surface to retire next: none made obsolete by this pass. The next cleanup target is reducing specialization/cold compile pressure in the route summaries or retained-unit shapes, not deleting existing route code.

-- repo-doer@macmini

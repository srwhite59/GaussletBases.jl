Completed the Be S+P decomposed/final-basis warm/cold timing audit.

Files changed:

- Added `docs/src/developer/blurb_logs/2026-06-11_mwg_decomposed_wl_cleanup/response.016.md`

Generated ignored developer artifacts under `tmp/work`:

- `tmp/work/be_atom_sp_decomposed_final_basis_warm_cold_timing_probe.jl`
- `tmp/work/be_atom_sp_decomposed_final_basis_warm_cold_timing_summary.txt`
- `tmp/work/be_atom_sp_decomposed_final_basis_warm_cold_phase_timings.tsv`
- `tmp/work/be_atom_sp_decomposed_final_basis_warm_cold_mixed_gto_subphase_timings.tsv`

No production source code was changed.

Probe protocol:

- Fresh process cold run using the existing Be S+P decomposed final-basis probe.
- One same-process warmup run.
- One same-process measured warm run.
- Same fixture as the current Be S+P probe:
  - Be, `Z = 4`
  - `q/ns = 5/5`
  - Be `cc-pV5Z`, `lmax = 1`
  - BasisSets path: `/Users/srw/Library/CloudStorage/Dropbox/GaussletModules/BasisSets`
  - density-density and RHF enabled

Physics stayed fixed:

- cold RHF total = `-14.574514244574639`
- measured warm RHF total = `-14.574514244574639`
- warm minus cold = `0.0`
- old nested/QW oracle total = `-14.574514244574694`
- measured warm delta from oracle = `5.5067062021407764e-14 Ha`
- final dimension = `636`
- retained gausslet dimension = `615`
- units/pairs = `131 / 8646`
- no full-parent CPB, direct Cartesian fallback, ordinary Cartesian IDA fallback, raw GTO final density-density, or generalized final solve was used.

Top-level timing:

| label | route elapsed | RHF elapsed | total elapsed |
| --- | ---: | ---: | ---: |
| cold | `169.936205334s` | `2.540070625s` | `172.476275959s` |
| warmup | `0.549007292s` | `1.890630792s` | `2.439638084s` |
| measured warm | `0.523052584s` | `1.880134542s` | `2.403187126s` |

Cold versus measured warm phase timings:

| phase | cold | warm | cold/warm |
| --- | ---: | ---: | ---: |
| residual_moment_matrices | `36.155260708s` | `0.002820542s` | `12818.55x` |
| electron_nuclear_by_center | `28.065223s` | `0.049385625s` | `568.29x` |
| overlap | `20.711751125s` | `0.005225084s` | `3963.91x` |
| gausslet_density_density | `19.2412885s` | `0.004252083s` | `4525.14x` |
| kinetic | `12.651555542s` | `0.090058958s` | `140.48x` |
| mixed_gto_blocks | `9.199361667s` | `0.216573708s` | `42.48x` |
| final_basis_density_density | `8.028131375s` | `0.014046792s` | `571.53x` |
| combined_one_electron_assembly | `2.052166209s` | `0.009078292s` | `226.05x` |
| shellification_decomposed_inventory | `2.382645s` | `0.001274458s` | `1869.54x` |
| combined_gto_basis_layout | `1.567163708s` | `0.000020667s` | `75829.28x` |
| final_basis_projection | `0.118247125s` | `0.106098708s` | `1.11x` |
| residual_mwg_representation | `0.18844075s` | `0.007870542s` | `23.94x` |
| parent_axis_setup | `0.054773917s` | `0.016139833s` | `3.39x` |
| RHF solve | `2.540070625s` | `1.880134542s` | `1.35x` |

Mixed-GTO subphase timing:

| phase | cold | warm | cold/warm |
| --- | ---: | ---: | ---: |
| GTO/GTO self block construction | `0.765639792s` | `0.047742333s` | `16.04x` |
| factorized projection total | `0.338598417s` | `0.168615167s` | `2.01x` |

Attribution:

- The remaining long cold route timings are almost entirely compilation.
- The measured warm route construction is only about `0.52s`.
- The measured warm total is dominated by the probe-local RHF solve at about `1.88s`, not by route construction.
- Among warm route-construction phases, `mixed_gto_blocks` is the largest at about `0.217s`, with most of that in the factorized projection. This is not large enough to justify another algorithmic rewrite before precompile work.
- `final_basis_projection` is also a real warm cost at about `0.106s`, but it is modest relative to the RHF solve.

Precompile recommendation:

- A precompile workload extension is justified, but not added in this pass.
- Existing `src/precompile_workloads.jl` only compiles the small decomposed WL one-body route: overlap, kinetic, and electron-nuclear by-center.
- The cold Be data shows that a narrow production workload should target the existing atom+GTO final-basis route surfaces, especially:
  - shellification-backed decomposed WL overlap/kinetic/electron-nuclear for the larger atom+GTO route shape;
  - factorized mixed-GTO block construction;
  - residual moment matrix construction;
  - gausslet density-density and final-basis density-density.
- That workload should not depend on the user-local GaussletModules `BasisSets` path. The right next step is to use or add a tiny repo-local/synthetic supplement fixture that exercises the same code paths without carrying acceptance logic into precompile code.

Algorithmic optimization recommendation:

- No immediate algorithmic rewrite is justified by the warm data.
- If a later pass optimizes warm construction anyway, the first target should be the factorized mixed-GTO projection path, because it is the largest remaining warm route phase.
- The present evidence says cold latency, not warm route runtime, is the user-facing problem.

Validation run:

- `julia --project=. tmp/work/be_atom_sp_decomposed_final_basis_warm_cold_timing_probe.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Deletion/shrinkage report:

- No old production code, tests, metadata, or compatibility path became unnecessary in this timing-only pass.
- Nothing was deleted or simplified because this pass only attributed cold versus warm cost and intentionally avoided route-shape changes.
- No test was added; the generated probe is an ignored developer artifact under `tmp/work`, not long-term coverage.
- Remaining stale/duplicate surface to retire next: none identified by this pass. The next concrete cleanup should be precompile-workload design, not deletion.

-- repo-doer@macmini

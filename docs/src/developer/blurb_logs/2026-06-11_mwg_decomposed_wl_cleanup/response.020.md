Measured the effect of the decomposed WL inventory-summary shape cleanup. No production source, tests, or docs were changed in this pass.

Files changed:

- `docs/src/developer/blurb_logs/2026-06-11_mwg_decomposed_wl_cleanup/response.020.md`

Generated/updated ignored developer artifacts:

- `tmp/work/atom_gto_specialization_shape_audit_summary.txt`
- `tmp/work/be_atom_sp_decomposed_final_basis_warm_cold_timing_summary.txt`
- `tmp/work/be_atom_sp_decomposed_final_basis_warm_cold_phase_timings.tsv`
- `tmp/work/be_atom_sp_decomposed_final_basis_warm_cold_mixed_gto_subphase_timings.tsv`
- the existing Be S+P probe summary/phase artifacts refreshed by the warm/cold probe include path

Specialization-shape audit:

| field | before side-7 | before side-15 | after side-7 | after side-15 |
| --- | ---: | ---: | ---: | ---: |
| parent axis counts | `(7, 7, 7)` | `(15, 15, 15)` | `(7, 7, 7)` | `(15, 15, 15)` |
| retained dimension | `223` | `615` | `223` | `615` |
| unit count | `27` | `131` | `27` | `131` |
| pair count | `378` | `8646` | `378` | `8646` |
| unit pairs storage | `UnitPairIndexTable` | `UnitPairIndexTable` | `UnitPairIndexTable` | `UnitPairIndexTable` |
| retained units storage | `Vector{RetainedUnitRecord}` | `Vector{RetainedUnitRecord}` | `Vector{RetainedUnitRecord}` | `Vector{RetainedUnitRecord}` |
| `unit_keys` type | `NTuple{27, Symbol}` | `NTuple{131, Symbol}` | `Vector{Symbol}` | `Vector{Symbol}` |
| `unit_summaries` type | 27-element tuple | 131-element tuple | `Vector{NamedTuple{...}}` | `Vector{NamedTuple{...}}` |
| `pair_summaries` type | `NTuple{378, ...}` | `Symbol` | compact count/status `NamedTuple` | compact count/status `NamedTuple` |
| factorized sidecar result type | stable `NamedTuple` | stable `NamedTuple` | stable `NamedTuple` | stable `NamedTuple` |
| coefficient matrix size | `(343, 223)` | `(3375, 615)` | `(343, 223)` | `(3375, 615)` |

The concrete inventory-size tuple mismatch is gone for the three fields targeted by the cleanup. Remaining expected shape differences are physical route size differences: parent product size, retained dimension, unit count, pair count, and coefficient matrix size.

Remaining non-size-stable or broad shapes visible in the audit:

- `metadata_type` still carries full typed parent-axis bundle objects.
- The factorized sidecar is still a report-style `NamedTuple`, although its type is stable between side-7 and side-15 in this audit.
- The atom+GTO route result still mixes compute payloads with report fields outside this shape audit.

Be S+P warm/cold timing comparison:

| metric | previous post-precompile | after inventory cleanup | change |
| --- | ---: | ---: | ---: |
| cold route elapsed | `159.006798084s` | `30.843591209s` | `-128.163206875s` |
| cold total elapsed | `161.9370205s` | `33.516908959s` | `-128.420111541s` |
| warm route elapsed | `0.498273209s` | `0.492626583s` | `-0.005646626s` |
| warm total elapsed | `2.486322293s` | `2.490025333s` | `+0.00370304s` |

The inventory cleanup materially improved measured cold route latency in this fresh-process probe, while warm route timing stayed effectively unchanged. The warm total difference is RHF/noise-sized and not a route regression.

Physics comparability from the measured warm run:

- `rhf_total_energy = -14.574514244574639`
- old nested/QW oracle total = `-14.574514244574694`
- delta = `5.5067062021407764e-14 Ha`
- final dimension = `636`
- retained gausslet dimension = `615`
- units/pairs = `131 / 8646`
- no full-parent CPB, direct Cartesian fallback, ordinary Cartesian IDA fallback, raw GTO final density-density, or generalized final-basis solve.

Current cold phase timings, ranked:

| phase | cold seconds | warm seconds | cold/warm |
| --- | ---: | ---: | ---: |
| residual moment matrices | `6.912517084` | `0.00290625` | `2378.5` |
| electron-nuclear by center | `4.182564875` | `0.049564792` | `84.39` |
| gausslet density-density | `2.232777208` | `0.004201333` | `531.44` |
| mixed GTO blocks | `1.936833792` | `0.365753083` | `5.30` |
| overlap | `1.903553875` | `0.004397875` | `432.83` |
| combined one-electron assembly | `1.871207333` | `0.011082084` | `168.85` |
| final-basis density-density | `1.441699584` | `0.014010791` | `102.90` |
| kinetic | `1.366693209` | `0.000685709` | `1993.11` |
| shellification inventory | `0.7078105` | `0.001227` | `576.86` |

Ranked next compile-pressure candidates:

1. Residual moment bundles / residual-moment matrix construction. This is now the largest cold phase by a wide margin and is almost pure first-call compilation in this probe.
2. Giant atom+GTO route result and staged payload/report coupling. Several phases with tiny warm times still compile independently, suggesting full staged result shapes are still being carried into helpers that need narrower payloads.
3. One-body route method specialization for electron-nuclear, overlap, and kinetic. Warm runtime is fine, but cold ratios remain high for repeated operator surfaces.
4. Final-basis density-density assembly shape. It is no longer the largest cost, but still has a large cold/warm ratio.
5. Parent-axis metadata. It remains visible as a large typed metadata object in the shape audit, but the measured parent/inventory cold phases are not the dominant cost after the tuple cleanup.
6. Factorized sidecar. Its result type is stable across side-7 and side-15 in this audit, so it is not the next obvious target.
7. `pair_keys` small/large shape split. It was not implicated by this audit after the hot tuple cleanup and should stay lower priority unless a later trace shows it on the cold path.
8. Timing closures. Possible background pressure, but not a first target from these measurements.

Decision:

- Do not revert the inventory cleanup; it removed the bad type shape and the measured cold route latency improved substantially.
- Do not broaden precompile in this pass. The next useful cleanup should likely narrow residual moment bundle/result shapes or split compute payloads from audit reports in the atom+GTO route.
- Do not optimize warm kernels from this evidence; warm route time remains about `0.49s`, and the route is no longer the measured user-facing bottleneck in the same-process run.

Validation run:

- `julia --project=. tmp/work/atom_gto_specialization_shape_audit.jl`
- `julia --project=. tmp/work/be_atom_sp_decomposed_final_basis_warm_cold_timing_probe.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Deletion/shrinkage report:

- No old code, tests, metadata, or compatibility path became newly unnecessary in this measurement-only pass.
- Nothing was deleted or simplified because the blurb asked for evidence before another struct or precompile pass, and the evidence did not justify a tiny safe cleanup inside this pass.
- No tests were added. Existing ignored `tmp/work` probes were used.
- Remaining stale or duplicate surfaces to retire next:
  - residual moment matrix bundle/result fields that likely compile as report-shaped payloads;
  - atom+GTO staged route results that mix compute data with audit/report fields;
  - parent-axis bundle objects carried through metadata rather than a route-owned parent-axis context;
  - possible small/large `pair_keys` summary split if a later trace shows it still contributes.

-- repo-doer@macmini

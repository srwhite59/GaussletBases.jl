Audited the remaining cold cost inside the decomposed WL one-electron matrix-set phase. No production source was changed.

Files changed:

- `docs/src/developer/blurb_logs/2026-06-11_mwg_decomposed_wl_cleanup/response.023.md`

Generated ignored developer artifacts:

- `tmp/work/be_atom_sp_one_electron_matrix_set_subphase_probe.jl`
- `tmp/work/be_atom_sp_one_electron_matrix_set_subphase_summary.txt`
- `tmp/work/be_atom_sp_one_electron_matrix_set_subphase_timings.tsv`
- refreshed Be S+P warm/cold probe artifacts under `tmp/work`

Audit result:

- The committed one-electron matrix-set seam exists and currently calls:
  - `_route_global_decomposed_wl_factorized_one_body_matrix(..., :overlap)`
  - `_route_global_decomposed_wl_factorized_one_body_matrix(..., :kinetic)`
  - `_route_global_decomposed_wl_factorized_electron_nuclear_by_center_matrix(...)`
- I added a developer-only subphase probe that reconstructs the Be S+P decomposed inventory and calls those current private helpers directly.
- The probe used the same `q/ns = 5/5`, `Z = 4`, parent side count `15`, retained dimension `615`, `131` units, and `8646` unit pairs as the Be S+P route.
- The second probe run used the same one-field metadata shape as the live Be probe.

Subphase timings from the direct helper probe:

| subphase | cold seconds | warm seconds |
| --- | ---: | ---: |
| factorized retained basis lookup/extract | `0.038492417` | `0.0000015` |
| overlap helper | `0.047759` | `0.000300583` |
| kinetic helper | `0.000734333` | `0.000687667` |
| electron-nuclear helper | `0.359268209` | `0.046163792` |
| full matrix set after prior helper calls | `0.2198505` | `0.125867125` |

Interpretation:

- Direct helper result construction is not reproducing the full route’s `~4.36s` cold one-electron phase.
- Repeated factorized retained-basis setup is not the measured issue in the direct helper probe; the cached factorized sidecar is effectively free after the first lookup.
- Overlap and kinetic helper construction are not meaningful warm bottlenecks.
- Electron-nuclear remains the only nontrivial local subphase, but its direct cold/warm timings are far smaller than the full route cold phase.
- The evidence does not justify replacing the helper internals in this pass. The remaining full-route cold cost is more likely first-call specialization/compilation from the route context, TimeG-wrapped call graph, or staged atom+GTO route payload shape than simple inner helper report construction.

Be S+P route timing rerun:

| metric | previous response.022 | this audit rerun | change |
| --- | ---: | ---: | ---: |
| cold route elapsed | `25.465580292s` | `24.996941042s` | `-0.46863925s` |
| cold total elapsed | `28.130448667s` | `27.644617375s` | `-0.485831292s` |
| warm route elapsed | `0.495644334s` | `0.502428875s` | `+0.006784541s` |
| warm total elapsed | `2.494674792s` | `2.413790167s` | `-0.080884625s` |
| one-electron matrix-set cold | `4.432942625s` | `4.365419666s` | `-0.067522959s` |
| one-electron matrix-set warm | `0.054274709s` | `0.059145583s` | `+0.004870874s` |

Physics comparison to old nested/QW Be S+P oracle:

- `rhf_total_energy = -14.574514244574639`
- old nested/QW oracle total = `-14.574514244574694`
- delta = `5.5067062021407764e-14 Ha`
- final dimension = `636`
- retained gausslet dimension = `615`
- units/pairs = `131 / 8646`
- no full-parent CPB, direct Cartesian fallback, ordinary Cartesian IDA fallback, raw GTO final density-density, or generalized final-basis solve.

Decision:

- No production cleanup was made in this pass.
- The subphase audit did not support a narrow direct-matrix rewrite inside `_route_global_decomposed_wl_factorized_one_electron_matrix_set(...)`.
- The next useful compile-pressure target should be a more precise compile-attribution probe around the full atom+GTO route context, not another unmeasured inner helper rewrite.

Validation run:

- `julia --project=. tmp/work/be_atom_sp_one_electron_matrix_set_subphase_probe.jl`
- `julia --project=. tmp/work/be_atom_sp_decomposed_final_basis_warm_cold_timing_probe.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Deletion/shrinkage report:

- No old production code, tests, metadata, or compatibility path became unnecessary in this audit-only pass.
- Nothing was deleted or simplified because the measurements did not justify a source change.
- No tests were added. The only new artifact is an ignored `tmp/work` developer timing probe.
- Remaining stale or duplicate surfaces to retire next:
  - atom+GTO route result still mixes compute payloads with audit/report fields;
  - the one-electron matrix set still sits in legacy overlap/kinetic/nuclear payload slots in the private route result;
  - parent-axis bundle objects still ride through metadata instead of a route-owned parent-axis context;
  - if future compile traces implicate route context shape, split compute payloads from report summaries before adding broader precompile workloads.

-- repo-doer@macmini

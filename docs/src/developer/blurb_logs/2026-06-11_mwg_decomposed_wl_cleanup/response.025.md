Attributed the cold compile pressure inside the density-enabled Be atom+GTO
route. No production source was changed.

Files changed:

- `docs/src/developer/blurb_logs/2026-06-11_mwg_decomposed_wl_cleanup/response.025.md`

Generated ignored developer artifacts:

- `tmp/work/be_atom_sp_density_route_compile_attribution_probe.jl`
- `tmp/work/be_atom_sp_density_route_compile_attribution_summary.txt`
- `tmp/work/be_atom_sp_density_route_compile_attribution.tsv`

Probe setup:

- The probe runs fresh child Julia processes for cold subphase timings.
- It then runs warmup and measured warm variants in one process.
- Fixture is the current Be S+P atom+GTO decomposed route:
  - `Z = 4`
  - `q/ns = 5/5`
  - final dimension `636`
  - retained gausslet dimension `615`
  - raw/retained supplement count `21 / 21`
  - `131` units and `8646` unit pairs
- The probe first exposed that the one-electron route omits mixed GTO moment
  blocks by design. I updated the probe to build the moment-capable mixed block
  bundle as density-path setup data rather than treating the one-electron route
  as sufficient for residual-MWG moments.

Density subphase timing table:

| variant | cold elapsed | warm elapsed | property count | type-name length | shape |
| --- | ---: | ---: | ---: | ---: | --- |
| one-electron prerequisite route | `10.676180542s` | `0.685628916s` | `46` | `24945` | `unavailable` |
| mixed moment blocks only | `1.338955583s` | `0.368781584s` | `29` | `21712` | `unavailable` |
| residual moment matrices only | `7.178579542s` | `0.002445s` | `28` | `1152` | `unavailable` |
| density prerequisite sequence | `7.858426875s` | `0.296454875s` | `28` | `1152` | `unavailable` |
| gausslet density only | `2.308185875s` | `0.004224958s` | `47` | `2055` | `(615, 615)` |
| residual MWG representation only | `0.141823209s` | `0.007962708s` | `30` | `1054` | `(21, 3)` |
| residual MWG components only | `1.270469125s` | `0.015270917s` | `32` | `1191` | `(636, 636)` |
| residual MWG old-kernel call only | `0.264626459s` | `0.013451375s` | `5` | `805` | `(636, 636)` |
| component result only | `0.364438459s` | `0.000474291s` | `32` | `1191` | `(636, 636)` |
| density result only | `0.050661875s` | `0.000007334s` | `38` | `1761` | `(636, 636)` |
| final-basis density only | `1.517116s` | `0.014050333s` | `38` | `1761` | `(636, 636)` |
| gausslet plus residual | `2.449776458s` | `0.012074917s` | `2` | `3144` | `unavailable` |
| full density sequence | `4.061640709s` | `0.026057s` | `38` | `1761` | `(636, 636)` |

Result shape/type-size notes:

- The isolated final density-density result is not the giant type:
  - property count `38`
  - type-name length `1761`
  - final matrix shape `(636, 636)`
- The residual MWG representation is compact:
  - property count `30`
  - type-name length `1054`
  - residual centers/widths shape `(21, 3)`
- The residual MWG component result is also compact:
  - property count `32`
  - type-name length `1191`
- The large density-path object in this probe is the moment-capable mixed block
  bundle:
  - property count `29`
  - type-name length `21712`
- The one-electron prerequisite route remains large:
  - property count `46`
  - type-name length `24945`

Interpretation:

- The final density-density result constructor is not the cold bottleneck.
  `density_result_only` is about `0.051s` cold and effectively zero warm.
- The old residual-MWG analytic kernel is not the main cold bottleneck in this
  route context. The direct
  `_qwrg_final_residual_mwg_component_blocks(...)` call is about `0.265s` cold,
  while the full final-basis density call is about `1.52s` cold.
- `gausslet_density_density` is a real cold cost at about `2.31s`, but it is
  not enough to explain the `~13.5s` full-route density-enabled delta.
- The largest density-enabled cold attribution found here is the residual
  moment prerequisite path:
  - residual moment matrices only: about `7.18s` cold;
  - mixed moment block setup plus residual moment matrices: about `7.86s`
    cold.
- The mixed moment block bundle also has a large concrete type-name length
  (`21712`) and a noticeable warm cost (`0.369s`), so it is both a compile and
  runtime candidate.
- The full density sequence after prerequisites is only about `4.06s` cold and
  `0.026s` warm. This suggests the density-enabled full-route delta is caused
  more by residual moment prerequisites and moment-capable mixed block shape
  than by the final density matrix result wrapper.

Ranked next target:

1. Residual moment matrix construction and moment-capable mixed GTO block
   setup. Split `route_global_decomposed_wl_moment_matrix_set`,
   `route_global_combined_gto_residual_moment_matrices`, and the
   `include_moment_blocks = true` mixed-GTO path.
2. Moment-capable mixed block result shape. The type-name length is large and
   warm cost is still visible; this may need a lean moment payload or fused
   moment matrix-set path.
3. Gausslet density-density cold compile. It is about `2.31s` cold, compact in
   result shape, and probably a secondary precompile or specialization target.
4. Final residual-MWG component wrapper. It adds about `1.25s` cold over a
   `0.265s` direct old-kernel call, but the result object itself is compact.
5. Final density result construction. Low priority from this evidence.

Fallback/final-basis checks:

- Materialized route objects report:
  - no full-parent CPB;
  - no direct Cartesian product assembly;
  - no ordinary Cartesian IDA fallback;
  - no raw GTO density-density accepted as final operator;
  - no generalized final-basis solve.

Validation run:

- `julia --project=. tmp/work/be_atom_sp_density_route_compile_attribution_probe.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Deletion/shrinkage report:

- No production code, tests, metadata, or compatibility path became unnecessary
  in this measurement-only pass.
- Nothing was deleted or simplified because the evidence points to a next
  attribution target, not a tiny safe cleanup.
- No tests were added. The only new executable artifact is an ignored
  `tmp/work` developer timing probe.
- Remaining stale or duplicate surfaces to retire next:
  - moment-capable mixed GTO block payloads with large concrete result shape;
  - residual moment matrix construction that may carry more route context than
    downstream density-density needs;
  - giant full atom+GTO route compute/report result;
  - one-electron prerequisite route objects carried into density-only probes
    when compact matrices and moment inputs may be enough.

-- repo-doer@macmini

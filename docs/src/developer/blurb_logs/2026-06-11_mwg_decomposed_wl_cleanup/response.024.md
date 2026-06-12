Attributed the remaining Be atom+GTO cold compile pressure at the full-route level. No production source was changed.

Files changed:

- `docs/src/developer/blurb_logs/2026-06-11_mwg_decomposed_wl_cleanup/response.024.md`

Generated ignored developer artifacts:

- `tmp/work/be_atom_sp_full_route_compile_attribution_probe.jl`
- `tmp/work/be_atom_sp_full_route_compile_attribution_summary.txt`
- `tmp/work/be_atom_sp_full_route_compile_attribution.tsv`

Probe setup:

- The probe runs fresh child Julia processes for cold route timings.
- It then runs warmup and measured warm variants in one process.
- Fixture is the current Be S+P atom+GTO route:
  - `Z = 4`
  - `q/ns = 5/5`
  - parent side count `15`
  - retained gausslet dimension `615`
  - final dimension `636`
  - `131` units and `8646` unit pairs
- No production code, public API, fallback, PQS, ECP, Be2, H2, driver, or acceptance fixture path was changed.

Variant timing table:

| variant | cold route | warm route | cold total | warm total | status |
| --- | ---: | ---: | ---: | ---: | --- |
| current timing, full density | `24.095634292s` | `0.494972917s` | `26.691211876s` | `2.574068958s` | full density/RHF |
| timing disabled, full density | `24.585734292s` | `0.492720542s` | `27.169000209s` | `2.407433459s` | full density/RHF |
| one-electron only, timing | `10.57373825s` | `0.437666292s` | `10.57373825s` | `0.437666292s` | one-electron route |
| one-electron only, no timing | `10.5277665s` | `0.444984875s` | `10.5277665s` | `0.444984875s` | one-electron route |
| minimal metadata, full density | `28.330674s` | `0.504395416s` | `30.948890125s` | `2.512836875s` | full density/RHF |

Cold top phases with timing enabled:

| variant | top cold phases |
| --- | --- |
| current full density | `decomposed_wl_one_electron_matrix_set:4.123655`, `gausslet_density_density:2.127648`, `combined_one_electron_assembly:1.961618`, `mixed_gto_blocks:1.825074`, `final_basis_density_density:1.344373` |
| one-electron only | `decomposed_wl_one_electron_matrix_set:4.236753`, `mixed_gto_blocks:1.468092`, `combined_one_electron_assembly:0.801271`, `shellification_decomposed_inventory:0.742028`, `combined_gto_basis_layout:0.179342` |
| minimal metadata full density | `decomposed_wl_one_electron_matrix_set:6.390043`, `gausslet_density_density:2.204341`, `combined_one_electron_assembly:2.045183`, `mixed_gto_blocks:1.847031`, `shellification_decomposed_inventory:1.522698` |

Route-result shape summary:

| variant family | route property count | route type-name length |
| --- | ---: | ---: |
| full density | `46` | `43929` with current metadata, `43788` with minimal metadata |
| one-electron only | `46` | `24945` |

Physics/fallback checks:

- Full-density RHF total in all full-density variants: `-14.574514244574639`.
- Old nested/QW oracle total: `-14.574514244574694`.
- Delta: `5.5067062021407764e-14 Ha`.
- Fallback/final-basis flags stayed clean:
  - no full-parent CPB;
  - no direct Cartesian product assembly;
  - no ordinary Cartesian IDA fallback;
  - no raw GTO final density-density accepted;
  - no generalized final-basis solve.

Interpretation:

- Disabling timing sinks does not materially reduce cold route time:
  - current full density: `24.095634292s`;
  - no timing full density: `24.585734292s`.
  - Timing wrappers are not the next cleanup target.
- Minimal metadata does not reduce cold route time:
  - current metadata full density: `24.095634292s`;
  - minimal metadata full density: `28.330674s`.
  - Parent-axis metadata cleanup may still be conceptually useful, but this probe does not make it the next cold-latency target.
- `build_density_density = false` removes most cold route time:
  - full density current: `24.095634292s`;
  - one-electron only current: `10.57373825s`.
  - This points to the final density-density / residual MWG route shape as the next compile-pressure target.
- Route result type size is large:
  - full-density type-name length about `44k`;
  - one-electron-only type-name length about `25k`.
  - This supports a later compute/result split, but the biggest variant delta in this probe is still density-density enabled versus disabled.

Ranked next target:

1. Final density-density / residual MWG route shape and compile attribution.
   The full-density route adds about `13.5s` cold route time over the one-electron-only route.
2. Giant full atom+GTO route result type.
   The route result type remains very large and should eventually be split into lean compute payload plus audit summary.
3. Full-route one-electron matrix-set compile context.
   It remains the largest single top-level cold phase, but direct helper subphase timing did not justify an inner helper rewrite.
4. Parent-axis metadata and route context objects.
   Still a stale surface, but minimal metadata did not improve timing in this probe.
5. Timing wrapper specialization.
   Low priority from this evidence.

Validation run:

- `julia --project=. tmp/work/be_atom_sp_full_route_compile_attribution_probe.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Deletion/shrinkage report:

- No old production code, tests, metadata, or compatibility path became unnecessary in this measurement-only pass.
- Nothing was deleted or simplified because the probe identified the next target but did not support a tiny safe source cleanup.
- No tests were added. The only new artifact is an ignored `tmp/work` developer timing probe.
- Remaining stale or duplicate surfaces to retire next:
  - final density-density / residual MWG route result shapes;
  - giant atom+GTO compute/report result;
  - one-electron payloads stored in legacy overlap/kinetic/nuclear route slots;
  - parent-axis bundle objects carried through metadata.

-- repo-doer@macmini

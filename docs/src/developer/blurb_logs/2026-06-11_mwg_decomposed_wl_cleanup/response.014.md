What changed:

- Added mixed-GTO subphase timing to the private route-global mixed block path.
- Threaded a `mixed_gto_subphase_timings` sink through the private Be atom+GTO
  seam and the Be S+P probe.
- Replaced the old per-unit full CPB/GTO local bundle call in the active
  mixed-GTO route with narrower private helpers:
  - one route-level GTO/GTO self-block construction;
  - one per-unit mixed CPB/GTO local block construction.
- GTO/GTO self overlap, kinetic, position, x2, and nuclear by-center blocks are
  now built once per `route_global_mixed_gto_blocks_from_decomposed_units(...)`
  call and reused across retained units.
- Mixed CPB/GTO local blocks are still built per retained unit, because they
  depend on the unit source CPB.
- Updated artifacts:
  - `tmp/work/be_atom_sp_decomposed_final_basis_summary.txt`
  - `tmp/work/be_atom_sp_decomposed_final_basis_phase_timings.tsv`
  - `tmp/work/be_atom_sp_mixed_gto_subphase_timings.tsv`

Files changed:

- `src/cartesian_pair_block_materialization/route_global_atom_gto_final_basis_route.jl`
- `src/cartesian_pair_block_materialization/route_global_mixed_gto_blocks.jl`
- `tmp/work/be_atom_sp_decomposed_final_basis_probe.jl`
- `tmp/work/be_atom_sp_decomposed_final_basis_summary.txt`
- `tmp/work/be_atom_sp_decomposed_final_basis_phase_timings.tsv`
- `tmp/work/be_atom_sp_mixed_gto_subphase_timings.tsv`
- `docs/src/developer/blurb_logs/2026-06-11_mwg_decomposed_wl_cleanup/response.014.md`

Be S+P result:

- Final probe status:
  `:materialized_final_basis_be_sp_rhf_probe`.
- Route status:
  `:materialized_decomposed_atom_gto_final_basis_density_density_route`.
- RHF status:
  `:converged_final_basis_be_sp_restricted_hartree_fock`.
- RHF iterations: `27`.
- Final dimension: `636`.
- Final overlap identity error: `1.0125523569644675e-10`.
- Final density-density symmetry error: `0.0`.
- New RHF one-electron energy: `-19.066200470580668`.
- Old oracle one-electron energy: `-19.06620047058102`.
- New RHF electron-electron energy: `4.4916862260060055`.
- Old oracle electron-electron energy: `4.491686226006327`.
- New RHF total: `-14.574514244574662`.
- Old oracle RHF total: `-14.574514244574694`.
- RHF total delta from old oracle: `3.197442310920451e-14`.
- Fallback flags remain false:
  full-parent CPB, direct Cartesian product assembly, ordinary Cartesian IDA
  operators, raw GTO final density-density, and generalized final solve.

Before/after timing:

- Previous full Be S+P probe total: `357.444079292` seconds.
- New full Be S+P probe total: `342.53537725` seconds.
- Previous route elapsed: `354.822137417` seconds.
- New route elapsed: `339.913994709` seconds.
- Previous `mixed_gto_blocks`: `188.576837208` seconds.
- New `mixed_gto_blocks`: `177.196382542` seconds.
- Mixed-GTO improvement: about `11.38` seconds.

Mixed-GTO subphase timing:

```text
phase                                      elapsed_s          count
gto_gto_self_block_construction            0.823398625        1
per_unit_total                             169.936329462      131
unit_coefficient_construction              0.010550254        131
per_unit_provider_local_block_construction 168.066350504      131
support_coefficient_construction           0.000976289        131
retained_contraction                       0.00206471         131
row_placement_coverage                     0.000661169        131
```

Top-level phase timing after the change:

```text
phase                              elapsed_s
parent_axis_setup                  0.056173792
shellification_decomposed_inventory 2.449104834
combined_gto_basis_layout          1.597603459
overlap                            20.859601375
kinetic                            13.91740175
electron_nuclear_by_center         27.774019208
mixed_gto_blocks                   177.196382542
combined_one_electron_assembly     2.105359084
final_basis_projection             0.023578458
residual_moment_matrices           36.011697958
residual_mwg_representation        0.193488125
gausslet_density_density           19.338124042
final_basis_density_density        8.286646042
rhf_solve                          2.621382541
```

Interpretation:

- Repeated GTO/GTO self-block construction was removed from the hot loop, but it
  was not the main bottleneck.
- The remaining dominant subphase is per-unit mixed CPB/GTO local block
  construction: about `168.07` seconds over `131` units.
- Retained contraction, support coefficient construction, and row placement are
  negligible at this fixture size.
- The next replacement target should be a factorized or projected mixed-GTO
  path that avoids rebuilding full mixed local blocks for each retained unit.
  The route needs something closer to a retained-unit projection of reusable
  axis/cross tables rather than per-unit CPB-local mixed block materialization.

Validation:

- `julia --project=. tmp/work/be_atom_sp_decomposed_final_basis_probe.jl`
  passed and produced the RHF result above.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.
- No existing local test was run. The changed surface is the private Be
  mixed-GTO route-global path, and the Be probe directly exercises the modified
  full-density/RHF route. The available H/H2+ GTO acceptance test is useful but
  not cheap enough to run as an automatic extra gate for this timing pass.

Deletion/shrinkage:

- The active route no longer uses the old per-unit full
  `cpb_gto_supplement_local_operator_bundle(...)` path for route-global
  mixed-GTO assembly. That path remains available at the provider level for
  provider-bundle tests and direct provider consumers.
- The old route-local helpers that always bundled mixed and GTO/GTO self blocks
  together were replaced by narrower route-local helpers separating reusable
  GTO/GTO self data from per-unit mixed data.
- No tests were added.
- No compatibility shim, fallback route, public API, export, or broad metadata
  layer was added.
- No H/H2+ fixture-local wiring was deleted; it still owns fixture-specific
  acceptance assertions. It can be revisited after a shared route helper can
  replace the setup without weakening those scientific contracts.
- Remaining stale/duplicate surface to retire next:
  per-unit mixed CPB/GTO local block materialization as the production route for
  larger atom+supplement fixtures. It is correct but now measured as the main
  runtime bottleneck.

-- repo-doer@macmini

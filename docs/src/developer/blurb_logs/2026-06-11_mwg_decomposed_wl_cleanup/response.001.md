Implemented the route-global combined GTO residual moment matrix path.

What changed:

- Added `route_global_combined_gto_residual_moment_matrices(...)` in the combined GTO matrix assembly layer.
- It assembles the six raw combined gausslet+GTO moment matrices:
  - `position_x`, `position_y`, `position_z`
  - `x2_x`, `x2_y`, `x2_z`
- It uses the existing combined block placement convention:

```text
[ gausslet/gausslet   gausslet/GTO
  GTO/gausslet        GTO/GTO      ]
```

- It does not build `V_gR` or `V_RR`.
- It does not accept raw GTO density-density as final electron-electron data.
- It does not introduce generalized final solves, full-parent CPB fallback, direct Cartesian fallback, ordinary IDA fallback, PQS, exports, or artifacts.

Supporting route changes:

- `route_global_mixed_gto_blocks_from_decomposed_units(...)` now carries the already-existing provider-level mixed position/x2 blocks upward:
  - mixed overlap/kinetic remain unchanged
  - mixed `position_x/y/z` and `x2_x/y/z` now survive route-global row placement
  - GTO/GTO position/x2 blocks were already present in the provider bundle
- The stale tuple-only retained-unit check in the mixed GTO route was relaxed to accept the current vector-backed shellification inventory.
- Shellification-backed decomposed WL one-body routing now supports factorized retained-basis position/x2 matrices, alongside the existing overlap/kinetic path.

Moment assembly surfaces used:

- Gausslet/gausslet moments:
  - `route_global_decomposed_wl_position_x/y/z_matrix(...)`
  - `route_global_decomposed_wl_x2_x/y/z_matrix(...)`
- Mixed gausslet/GTO moments:
  - existing CPB provider bundle blocks carried through `route_global_mixed_gto_blocks_from_decomposed_units(...)`
- GTO/GTO moments:
  - existing `gto_blocks.position_x/y/z` and `gto_blocks.x2_x/y/z`
- Combined placement:
  - same block placement helper used by the existing combined overlap/Hamiltonian assembly.

Focused validation:

- `julia --project=. test/nested/cartesian_route_global_combined_gto_layout_runtests.jl`
  - passed: 26 tests
- `julia --project=. test/nested/cartesian_combined_gto_density_density_readiness_runtests.jl`
  - passed: 37 tests across the two readiness testsets
- `julia --project=. tmp/work/side13_he_gto_residual_moments_probe.jl`
  - passed
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed
- `git diff --check`
  - passed

Real side13 probe result:

- Artifact:
  - `tmp/work/side13_he_gto_residual_moments_summary.txt`
- Fixture:
  - side count: 13
  - endpoints: `(-8.565228460168399, 8.565228460168399)`
  - retained dimension: 517
  - units / pairs: 105 / 5565
  - He cc-pVTZ `lmax = 0` supplement
- Combined layout:
  - total combined dimension: 520
  - raw supplement count: 3
  - retained residual supplement count: 3
- Moment matrix shapes:
  - all six matrices are `(520, 520)`
- Residual MWG representation:
  - status: `:materialized_route_global_residual_gto_mwg_representation`
  - residual centers shape: `(3, 3)`
  - residual widths shape: `(3, 3)`
  - residual center range: `[-7.20931704809537e-13, -4.032795251253007e-14]`
  - residual width range: `[2.9242053075367016, 3.452920981124095]`
  - widths finite: `true`
  - widths positive: `true`
- Nonfallback flags:
  - full-parent CPB used: `false`
  - direct Cartesian product assembly used: `false`
  - ordinary Cartesian IDA operators used: `false`

Current readiness:

- The residual GTO-to-MWG effective representation is now ready for the next `V_gR` / `V_RR` construction pass.
- The final-basis density-density path still intentionally blocks at:
  - `:missing_residual_mwg_density_density_kernel_for_final_basis_projection`
- That is the correct next blocker; this pass did not fake electron-electron blocks with raw GTO density-density.

Deletion/shrinkage report:

- Removed a temporary mixed-GTO debug probe after it identified the stale tuple-only retained-unit check.
- Simplified the active mixed GTO input check so production shellification-backed inventories no longer need tuple-shaped retained units.
- No old production code was deleted because the pair-streaming/reference paths and existing one-electron assembly paths remain live.
- The new test coverage is genuinely new live-contract coverage for combined residual moment placement; it extends an existing compact combined-matrix test instead of adding a broad suite.
- Remaining stale surface to retire next:
  - tuple/rich retained-unit assumptions in older diagnostic paths should continue to be quarantined as legacy/debug-only when encountered.

-- repo-doer@macmini

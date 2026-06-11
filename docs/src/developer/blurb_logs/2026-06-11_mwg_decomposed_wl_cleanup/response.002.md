Implemented the first residual MWG density-density assembly path for final-basis He + GTO.

What changed:

- `route_global_combined_gto_final_basis_density_density_matrix(...)` now accepts the route-owned sources needed to build residual MWG electron-electron blocks:
  - `decomposed_inventory`
  - `parent_axis_counts`
  - `parent_axis_bundle_object`
  - `coulomb_expansion`
- When those sources and a materialized residual MWG representation are supplied, it now builds:
  - `V_gg` from the supplied decomposed WL density-density matrix
  - `V_gR` / `V_Rg` through the residual MWG component kernel
  - `V_RR` through the residual MWG component kernel
  - final `V_final = [V_gg V_gR; V_Rg V_RR]`
- It still blocks, with a sharper source blocker, when residual MWG data is available but route-owned inventory/axis sources are not supplied.

MWG kernel surfaces used:

- `ParentGaussletBases._qwrg_final_residual_mwg_component_blocks(...)`
- internally, that uses `_qwrg_mwg_interaction_components(...)`

Parent-product-to-final contraction source:

- The contraction matrix comes from `_white_lindsey_decomposed_factorized_retained_basis(...)`.
- That sidecar now exposes the dense coefficient matrix it already builds from shellification-backed retained-unit coefficient maps.
- Source recorded as:
  - `:shellification_retained_unit_coefficient_maps`
- Old fixed-block matrices are not route authority and are not used as the contraction source.

Real side13 validation:

- Artifact:
  - `tmp/work/side13_he_gto_residual_moments_summary.txt`
- Fixture:
  - side count: 13
  - endpoints: `(-8.565228460168399, 8.565228460168399)`
  - retained dimension: 517
  - units / pairs: 105 / 5565
  - He cc-pVTZ `lmax = 0` supplement
- Combined basis:
  - total combined dimension: 520
  - raw supplement count: 3
  - retained residual supplement count: 3
  - final dimension: 520
- Density-density status:
  - decomposed WL `V_gg`: materialized
  - final-basis density-density: `:materialized_route_global_combined_gto_final_basis_density_density_matrix`
  - blocker: `nothing`
- Shapes:
  - `V_final`: `(520, 520)`
  - `V_gR`: `(517, 3)`
  - `V_RR`: `(3, 3)`
- Checks:
  - supplied `V_gg` block matches the final matrix gausslet block: `true`
  - `V_gR` materialized: `true`
  - `V_RR` materialized: `true`
  - `V_final` symmetric: `true`
  - `V_final` finite: `true`
- Nonfallback flags:
  - raw GTO density-density final operator: `false`
  - full-parent CPB used: `false`
  - direct Cartesian product assembly used: `false`
  - ordinary Cartesian IDA operators used: `false`

Current readiness:

- The final-basis density-density matrix is ready for a side13 He + GTO RHF diagnostic pass.
- This pass did not run or add He + GTO RHF acceptance.
- The next implementation step can consume `V_final` in the existing final orthonormal basis and run a clearly labeled side13 He + GTO RHF probe.

Validation:

- `julia --project=. test/nested/cartesian_combined_gto_density_density_readiness_runtests.jl`
  - passed
- `julia --project=. tmp/work/side13_he_gto_residual_moments_probe.jl`
  - passed
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed
- `git diff --check`
  - passed

Deletion/shrinkage report:

- No old production code was deleted in this pass.
- Reason: the new path fills a previously real missing block (`V_gR` / `V_RR`) rather than replacing an existing production implementation.
- One stale blocked-contract assertion was simplified: the readiness test now expects the sharper missing-source blocker `:missing_decomposed_wl_density_density_inventory` instead of the old generic missing-kernel blocker once residual MWG data is present.
- The new source exposure is not a compatibility adapter; it exposes the coefficient matrix already built by the shellification-backed factorized sidecar.
- No broad metadata test was added.
- Remaining stale/duplicate surfaces to retire next:
  - blocked text and docs that still describe final-basis density-density as missing the residual MWG kernel should be updated after the next RHF probe confirms consumption.
  - older fixed-block MWG helpers remain oracle/kernel references only.

-- repo-doer@macmini

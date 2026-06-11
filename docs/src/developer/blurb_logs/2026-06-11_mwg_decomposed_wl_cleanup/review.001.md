Review result:

Accepted. The pass implements the requested combined GTO residual moment matrix
path and keeps the trust boundary: it assembles position/x2 raw combined
matrices for MWG moment extraction, but does not build `V_gR` / `V_RR`, does
not run He + GTO RHF, and does not accept raw GTO density-density as final
electron-electron data.

Corrections made:

None. The new test extension is compact and checks live block-placement
behavior rather than helper vocabulary. The response includes the required
deletion/shrinkage report.

Validation:

- `julia --project=. test/nested/cartesian_route_global_combined_gto_layout_runtests.jl`
- `julia --project=. test/nested/cartesian_combined_gto_density_density_readiness_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Commit/push:

Committed as `c20a5e70 Add combined GTO residual moment matrices`; push handled
by repo-manager after review.

Next target:

Implement the residual MWG density-density kernel for the final-basis
`V_gR` / `V_RR` blocks, using the now-materialized residual MWG centers and
widths.

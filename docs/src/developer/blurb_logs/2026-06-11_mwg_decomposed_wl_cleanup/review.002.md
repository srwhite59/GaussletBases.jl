Review result:

Accepted. The pass materializes the final-basis residual MWG density-density
matrix blocks for side13 He + GTO without accepting raw GTO density-density as
final electron-electron data.

The implementation obtains the parent-product-to-final gausslet contraction
from the shellification-backed decomposed retained-unit coefficient maps via
the existing factorized retained-basis sidecar. Old fixed-block matrices are not
used as route authority.

Corrections made:

None. The tracked tests remain focused. The real side13 validation is a
developer probe under `tmp/work`, not a routine acceptance gate.

Validation:

- `julia --project=. test/nested/cartesian_combined_gto_density_density_readiness_runtests.jl`
- `julia --project=. test/nested/cartesian_route_global_combined_gto_layout_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Reviewed artifact:

- `tmp/work/side13_he_gto_residual_moments_summary.txt`

Commit/push:

Committed as `b3964be1 Add residual MWG density blocks`; push handled by
repo-manager after review.

Next target:

Run a side13 He + GTO final-basis RHF diagnostic using the now-materialized
final-basis density-density matrix. After that, shrink the He acceptance/default
nested test bloat as a separate cleanup pass.

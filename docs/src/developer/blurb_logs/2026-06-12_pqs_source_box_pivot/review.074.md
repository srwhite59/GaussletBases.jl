Review 074: accepted.

The correction addresses the pass-073 blocker. The axis-layer path now converts
`gaussian_factor_matrices(...)` output into term-first `Array{Float64,3}`
factor arrays before support-space electron-nuclear contraction, and validates
term count and matrix shape. The focused H1 gate now exercises that path at the
origin and compares it against the explicit `pgdg_intermediate` Gaussian-factor
path.

Manager validation:

- `julia --project=. -e 't = @elapsed include("test/nested/pqs_direct_retained_final_h1_runtests.jl"); println("elapsed_s=", t)'`
  passed, 32/32, elapsed about 4.67s in the warm process.
- `git diff --check` passed before the review file was added.

Status:

- pass 073 and pass 074 may be committed together;
- the support electron-nuclear helper is accepted as a route-owned support
  operator helper;
- off-origin use is now tied to explicit `axis_layers`, not origin PGDG factors;
- no H1/RHF/IDA/density/driver/export behavior was added.

Deletion/shrinkage:

- accepted the deletion of the test-local `_pqs_h1_support_nuclear_matrix`
  helper from pass 073;
- no additional deletion was expected from pass 074 because it corrected a live
  helper path before acceptance;
- the added assertion is live-contract coverage for the corrected axis-layer
  convention, not broad helper-vocabulary coverage.

-- repo-manager@macmini

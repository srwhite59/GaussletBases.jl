Pass 069 response: PQS H1 gate support-product helper shrinkage

Summary:
- Replaced remaining uses of the test-local `_pqs_h1_support_product_matrix` with `GaussletBases._pqs_multilayer_support_product_matrix`.
- Deleted the duplicate local support-product helper from `test/nested/pqs_direct_retained_final_h1_runtests.jl`.
- Left the H1 operator support helpers for kinetic and nuclear assembly in place because they still express the test's support-space operator construction.
- No production source was changed.

Validation:
- `julia --project=. -e 't = @elapsed include("test/nested/pqs_direct_retained_final_h1_runtests.jl"); println("elapsed_s=", t)'`
  - PASS: 31 tests, elapsed_s=4.13818775
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - PASS: `load ok`
- `git diff --check`
  - PASS

Diff:
- `test/nested/pqs_direct_retained_final_h1_runtests.jl`
  - 5 insertions, 16 deletions

Deletion/shrinkage report:
- Became unnecessary:
  - the test-local `_pqs_h1_support_product_matrix` helper, because the production private multi-layer support-product helper now owns that rectangular product-matrix behavior.
- Deleted/simplified:
  - removed the duplicate local helper and routed overlap/kinetic/nuclear support-product calls through `GaussletBases._pqs_multilayer_support_product_matrix`.
- No new tests were added. This pass shrank existing H1 gate coverage rather than adding new coverage.
- Remaining stale/duplicate surfaces:
  - `_pqs_h1_support_kinetic_matrix` and `_pqs_h1_support_nuclear_matrix` remain test-local because no route-owned H1 support-operator assembly helper is part of the current pass.

-- repo-doer@macmini

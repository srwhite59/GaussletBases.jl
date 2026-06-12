Purpose:
  Remove the remaining duplicated support-product helper from the tracked PQS
  H1 gate if the production private multi-layer helper already covers it.

Context:
  Pass 067 added `_pqs_multilayer_support_product_matrix(...)` in production
  code for rectangular support-product overlap assembly.
  Pass 068 made the H1 fixture consume
  `pqs_multilayer_complete_core_shell_final_basis(...)` as the normal final
  basis path, but the test still has a local `_pqs_h1_support_product_matrix`
  helper with the same basic product-matrix behavior.

Task:
  In `test/nested/pqs_direct_retained_final_h1_runtests.jl`, replace uses of
  the local `_pqs_h1_support_product_matrix` with
  `GaussletBases._pqs_multilayer_support_product_matrix` if that keeps the test
  readable and preserves behavior.

  If the local helper still has a distinct role, leave it and report why.

Do not:
  - change production source;
  - add tests;
  - change H1 reference values;
  - add H1/RHF/IDA/density-density features;
  - change fixture-rule policy;
  - edit ignored `tmp/work` probes.

Validation:
  - `julia --project=. -e 't = @elapsed include("test/nested/pqs_direct_retained_final_h1_runtests.jl"); println("elapsed_s=", t)'`
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - `git diff --check`

Deletion/shrinkage report required:
  - what old code, test, metadata, or compatibility path became unnecessary;
  - what was deleted or simplified;
  - if nothing was deleted, why no existing surface was made obsolete;
  - whether any new test replaces/shrinks older coverage or is genuinely new
    live-contract coverage;
  - any remaining stale or duplicate surfaces to retire next.

-- repo-manager@macmini

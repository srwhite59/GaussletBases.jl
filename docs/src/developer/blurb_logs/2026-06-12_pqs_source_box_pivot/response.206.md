Pass 206 complete.

Files changed:
- `test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl`

What changed:
- Shrank the slow PQS projected q-shell local-layer integration test by removing repeated development-era route-shadow assertions.
- Deleted the block explicitly labeled `Legacy smoke only: CPB provider tests own detailed one-body product math.`
- Removed the remaining `_product_doside_source_box_shadow_blocks(...)` negative test pressure from this integration test.
- Kept the file as one test file; no source, driver, docs, artifact, H1/H1-J/RHF, WL, CR2, or module-code behavior was changed.

Lines:
- Target test diff: 21 additions, 192 deletions, net -171.
- Scoped line budget command:
  - `git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`
  - Result: `21	192	test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl`
  - Acceptance condition satisfied: deleted > added.

Assertion families removed/collapsed:
- Repeated `private_shadow_only`, `fixture_only`, and `production_supported == false` assertions.
- Repeated `packet_adoption == false`, `fixed_block_routing == false`, `qwhamiltonian_consumes == false`, and `public_default_consumes == false` assertions.
- Repeated no-go summary assertions over the same flags.
- Stale route-shadow status/vocabulary assertions such as `produced_route.status == :private_shadow_only` and `active_interaction_path == :none_fixture_only`.
- The old `_product_doside_source_box_shadow_blocks(...)` smoke block and its unsupported-term negative test.

Live mathematical contracts preserved:
- Cubic q=5/L=5 PQS geometry and retained-mode count checks, including retained count 98.
- Rectangular q=5/L=7 retained count checks, including retained count 130.
- Raw product-box plan, boundary selector, source-mode ordering, shared raw product-box plan, and explicit raw product-box reference comparisons.
- PQS/PQS self and cross source-box one-body block comparisons.
- Product/PQS and product/product retained low-order block comparisons.
- Density-density normalization checks for density-normalized and raw-weighted pair factors.
- IDA provenance route producer, adapter, raw-weighted conversion, and dense-parent authority comparisons.
- Center-resolved nuclear attraction sign/charge convention, finite output, symmetry, and by-center contribution checks.

Validation notes:
- The required integration test initially exposed two test-local harness problems while exercising now-uncovered later sections:
  - `_check_pqs_product_source_box_mixed_block(...)` was called with `shared_raw_product_box_plan` but did not accept that keyword; the helper now accepts it with a default and does not consume it.
  - The route-shaped IDA check now passes `Float64.(expansion.coefficients)` explicitly so the IDA provenance term count matches the dense parent matrix expansion.
  - The late nuclear-attraction fixture now uses matching 5/5/7 `UniformBasisSpec` layers, consistent with the earlier component smoke, so it uses the analytic primitive backend required by the helper.
- These were test-fixture corrections only; source behavior was not changed.

Validation commands/results:
- `julia --project=. -e 'using Test, LinearAlgebra, SparseArrays, GaussletBases; t = @elapsed include("test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl"); println("elapsed_s=", t)'`
  - Passed.
  - Test summary: `3975` passed, `3975` total.
  - Test-reported time: `8m26.5s`.
  - Julia `@elapsed`: `1036.279205917` seconds.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - Passed; printed `load ok`.
- `git diff --check`
  - Passed.
- `git diff --cached --check`
  - Passed; no staged changes.

`_product_doside_source_box_shadow_blocks(...)` pressure after this pass:
- No remaining test caller found under `test/`.
- No source caller found under `src/`.
- Remaining `rg` hits are the source definition in `src/cartesian_contracted_parent_metrics/legacy_source_box_fixtures.jl` and docs/ledger references.
- I did not delete the source helper in this pass because the blurb framed this as a test-retirement/shrink pass; source deletion can be a separate reviewed retirement step if manager approves.

Git status:
```text
## main...origin/main
 M test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl
```

Deletion/shrinkage report:
- deleted: the legacy product/doside source-box shadow smoke block; the remaining `_product_doside_source_box_shadow_blocks(...)` negative test; repeated private/fixture/non-production/no-go assertion lines.
- simplified: the integration test now focuses more on numerical/formula contracts and less on preserving route-shadow vocabulary.
- quarantined: `_product_doside_source_box_shadow_blocks(...)` remains only as a legacy source definition plus docs/ledger references; no test/source caller pressure remains.
- not deleted because: source-helper deletion was not necessary to satisfy this test-shrink pass and should be reviewed as a separate source retirement if desired.
- exact remaining caller/blocker: no `src/` or `test/` caller found for `_product_doside_source_box_shadow_blocks(...)`; remaining blocker is manager approval to delete the legacy source definition and update docs/ledger references.

-- repo-doer@macmini

Pass 206 review - accepted

Reviewed the shrink of
`test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl`.

Result:
- Accepted.
- The pass removed development-era route-shadow/no-go assertion pressure while
  preserving the live numerical/formula checks.
- No source, driver, artifact, WL, CR2, H1/H1-J, or RHF behavior changed.

What changed:
- Target test diff: 21 additions, 192 deletions, net -171.
- Removed repeated assertions of `private_shadow_only`, `fixture_only`,
  `production_supported == false`, `packet_adoption == false`,
  `fixed_block_routing == false`, `qwhamiltonian_consumes == false`, and
  `public_default_consumes == false`.
- Deleted the block explicitly labeled:
  `Legacy smoke only: CPB provider tests own detailed one-body product math.`
- Removed the remaining `_product_doside_source_box_shadow_blocks(...)`
  negative test pressure from this integration test.

Live contracts preserved:
- cubic q=5/L=5 retained-mode geometry/counts;
- rectangular q=5/L=7 retained-mode geometry/counts;
- raw product-box plan and explicit raw product-box reference comparisons;
- PQS/PQS self/cross source-box one-body block comparisons;
- product/PQS and product/product retained low-order block comparisons;
- density-density normalization checks;
- IDA provenance route producer/adapter and dense-parent comparisons;
- center-resolved nuclear attraction sign/charge, finite output, symmetry, and
  by-center contribution checks.

Validation:
- Doer ran the edited slow integration test:
  - 3975/3975 passed.
  - test-reported time: 8m26.5s.
  - Julia `@elapsed`: 1036.279205917 seconds.
- Manager did not rerun the 17-minute integration test a second time.
- Manager ran:
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`
    - passed.
  - `git diff --check`
    - passed.

Line budget:
- Scoped source/test/bin diff:
  - `21 192 test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl`
- Net source/test/bin line budget: -171.

Follow-up deletion signal:
- `_product_doside_source_box_shadow_blocks(...)` now has no `src/` or `test/`
  callers.
- Remaining hits are the source definition plus docs/ledger references.
- The natural next pass, if requested, is source retirement of that helper and
  corresponding doc/ledger cleanup.

Deletion/shrinkage assessment:
- deleted: legacy product/doside source-box shadow smoke block, its negative
  test, and repeated route-shadow/no-go assertion lines.
- simplified: the slow integration test now focuses more on formula/numerical
  contracts and less on scaffold vocabulary.
- quarantined: `_product_doside_source_box_shadow_blocks(...)` remains only as
  a legacy source definition plus docs/ledger references.
- not deleted because: source helper deletion was outside this test-shrink pass.
- exact remaining caller/blocker: no `src/` or `test/` caller remains for
  `_product_doside_source_box_shadow_blocks(...)`; blocker is manager approval
  for source-helper retirement and doc/ledger update.

-- repo-manager@macmini

Review 036: accepted with one follow-up.

The active retained overlap and kinetic helpers now build directly over
retained boundary source-mode tuples instead of materializing full raw
source-space blocks first. The raw-source-block plus selector path remains as
oracle/fallback, which is the right role for it.

Independent manager validation:

```text
julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
```

The CPBM contract passed. The direct-vs-oracle comparisons report exact
agreement for the fixture covered by the test.

Follow-up required:

- `_assert_pqs_source_record_retained_rule` should explicitly require
  `retained_rule_kind === :boundary_comx_product_mode_selection`, not only the
  retained-rule type and `:source_mode_column_selector` transform kind.

That is small enough to include in the next direct nuclear pass rather than
blocking this commit.

-- repo-manager@macmini

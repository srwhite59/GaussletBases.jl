Purpose:

Add the first compact PQS source-mode retained-rule object. This is the first
implementation step for the source-box-first PQS route.

Why now:

The restart audit identified the missing first concept:

```text
RawProductBoxPlan + retained source-mode boundary rule
```

The first rule should be source-mode boundary selection, not shell rows, not
shell projection, and not Lowdin. For the first fixture:

```text
source modes: 5 x 5 x 5 = 125
interior modes: 3 x 3 x 3 = 27
boundary retained modes: 98
```

Exact task:

Implement a compact PQS boundary source-mode retained rule under
`CartesianRawProductSources`, then attach it to PQS retained-unit transform
contracts.

Primary code surfaces:

```text
src/cartesian_raw_product_sources/records.jl
src/cartesian_raw_product_sources/summaries.jl
src/cartesian_raw_product_sources/CartesianRawProductSources.jl
src/cartesian_retained_unit_transform_contracts/unit_contracts.jl
src/cartesian_retained_unit_transform_contracts/summaries.jl
```

Recommended shape:

```text
PQSRetainedRule or similarly named compact internal record:
    source_key
    source_mode_dims
    source_mode_ordering
    retained_rule_kind = :boundary_comx_product_mode_selection
    retained_mode_indices
    retained_column_indices
    retained_count
    transform_kind = :source_mode_column_selector
    shell_realization_materialized = false
    lowdin_cleanup_used = false
```

Constructor:

```text
pqs_boundary_product_mode_retained_rule(raw_plan)
```

or, if cleaner:

```text
pqs_boundary_product_mode_retained_rule(source_mode_dims;
    source_mode_ordering = :x_major_y_major_z_fast,
    source_key = :pqs_raw_product_source,
)
```

Selection rule:

- use the source mode ordering already supported by
  `source_mode_indices(...)`;
- select modes where at least one local axis index is `1` or the maximum index
  on that axis;
- for `(5, 5, 5)`, retained count must be `98`;
- retain deterministic column order matching source-mode order.

Attach to retained-unit transform contracts:

- in `_raw_product_source_contract_metadata(unit)` for
  `:pqs_shell_retained_unit`, after `raw_product_source_plan` is available,
  add the retained rule and a compact retained-rule summary;
- keep shell realization and Lowdin as planned/not materialized metadata only.

Trust boundary:

Do not implement pair-block contraction yet. Do not implement shell projection,
Lowdin, support-row contraction, electron-nuclear, IDA, density-density,
Hamiltonian, RHF, driver adoption, exports, artifacts, full-parent fallback, or
direct Cartesian fallback.

Test policy:

One compact module-contract test is allowed if needed. It should check the live
contract, not helper vocabulary:

- `(5, 5, 5)` gives retained count `98`;
- retained mode indices are deterministic;
- retained columns match boundary source-mode selection;
- summary says no shell realization and no Lowdin;
- PQS retained-unit transform-contract metadata carries the rule/summary.

Do not add broad metadata tests or driver tests.

Validation:

- run the smallest affected raw-product-source / retained-transform contract
  tests;
- run `julia --project=. -e 'using GaussletBases; println("load ok")'`;
- run `git diff --check`.

Deletion/shrinkage report required:

- what old code, test, metadata, or compatibility path became unnecessary;
- what was deleted or simplified;
- if nothing was deleted, why no existing surface was made obsolete yet;
- whether any new test replaces/shrinks older coverage or is genuinely new
  live-contract coverage;
- any remaining stale or duplicate surfaces to retire next.

Report back:

- write `.agent_handoffs/response.002.md.tmp`, then atomically rename to
  `.agent_handoffs/response.002.md`;
- also write the curated copy to
  `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.002.md`;
- include files changed;
- include the retained-rule shape and `(5,5,5) -> 98` evidence;
- include validation run;
- include deletion/shrinkage report;
- sign `-- repo-doer@macmini`.

-- repo-manager@macmini

Pass 264 response - retire legacy source-box fixture builders

Caller audit result:
- Live `src/test/bin` callers before deletion:
  - none outside `src/cartesian_contracted_parent_metrics/legacy_source_box_fixtures.jl`.
- Historical docs/log references:
  - references remain in developer docs, baton logs, archived handoffs, and the pass blurb.
  - I did not chase those historical references.
- Caller proof after deletion:
  - `rg -n "_pqs_contact_cap_product_doside_unit|_pqs_outer_mismatch_product_doside_units|_pqs_atom_box_support_dense_units|_pqs_contact_cap_parent_coefficient_matrix|_pqs_support_local_parent_coefficient_matrix|_pqs_route_retained_unit_fact|_pqs_route_retained_unit_fact_audit|_pqs_route_direct_or_support_body_mismatches|_pqs_route_product_box_construction_rule|_pqs_route_slab_piece_rule" src test bin --glob '!*.tmp'`
  - returned no matches.

Functions deleted:
- `_pqs_contact_cap_product_doside_unit`
- `_pqs_outer_mismatch_product_doside_units`
- `_pqs_atom_box_support_dense_units`
- `_pqs_contact_cap_parent_coefficient_matrix`
- `_pqs_support_local_parent_coefficient_matrix`
- `_pqs_route_direct_or_support_body_mismatches`
- `_pqs_route_policy_region_for_build`
- `_pqs_route_realization_descriptor_for_build`
- `_pqs_route_retained_unit_classification`
- `_pqs_route_coefficient_scope`
- `_pqs_route_safe_term_capability`
- `_pqs_route_slab_piece_rule`
- `_pqs_route_product_box_construction_rule`
- `_pqs_route_retained_unit_fact`
- `_pqs_route_product_box_slab_rule_count`
- `_pqs_route_retained_unit_fact_audit`

Functions preserved:
- active product/doside source-box pair-plan, density-density, nuclear/local-Gaussian, and raw-plan conversion diagnostics earlier in `legacy_source_box_fixtures.jl`;
- current atom-growth/contact-cap/outer-mismatch route code outside this legacy fixture file.

Validation:
- `git diff --check`
  - passed.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed after precompile: `load ok`
- caller proof command above
  - passed with no matches.

Scoped line count for `src + test + bin`:
- `src/cartesian_contracted_parent_metrics/legacy_source_box_fixtures.jl`: +0 / -999
- total: +0 / -999

Git status:
```text
## main...origin/main
 M src/cartesian_contracted_parent_metrics/legacy_source_box_fixtures.jl
```

Deletion/shrinkage result:
- deleted: 999 lines from the now-unblocked legacy retained-unit fixture builder tail.
- simplified: the file no longer carries the old contact-cap, outer-mismatch, atom-box support-dense fixture builders or their route-fact audit/helper scaffolding.
- quarantined: independent H2 PQS route/source-plan/final-basis/H1/H1-J/RHF/support-partition/provider code, provider blocks, supplement values, CR2/export, HamV6, public API, and unrelated source-box pair-plan/density/nuclear helpers were untouched.
- not deleted because: active product/doside source-box pair-plan and raw-plan diagnostic helpers still have live source roles in this legacy fixture file.
- exact remaining caller/blocker: no remaining `src/test/bin` caller found for the deleted function/helper names; remaining docs/log mentions are historical.

-- repo-doer@macmini

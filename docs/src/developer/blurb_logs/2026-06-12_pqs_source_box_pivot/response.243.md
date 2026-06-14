Pass 243 complete.

Files changed:
- `src/cartesian_contracted_parent_metrics/legacy_source_box_fixtures.jl`
- `docs/src/developer/raw_product_source_retained_transform_policy.md`
- `docs/src/developer/projected_q_shell_policy.md`

Deleted obsolete wrapper families from `legacy_source_box_fixtures.jl`:
- Contact cap:
  - `_PQS_CONTACT_CAP_SAFE_TERM_OPERATOR_TERMS`
  - `_PQS_CONTACT_CAP_UNSUPPORTED_OPERATOR_TERMS`
  - `_pqs_contact_cap_axis_factor_terms`
  - `_pqs_contact_cap_direct_support_oracle_entries`
  - `_pqs_contact_cap_direct_support_oracle_block`
  - `_pqs_contact_cap_safe_term_operator_comparison`
- Outer mismatch:
  - `_PQS_OUTER_MISMATCH_SAFE_TERM_OPERATOR_TERMS`
  - `_PQS_OUTER_MISMATCH_UNSUPPORTED_OPERATOR_TERMS`
  - `_pqs_outer_mismatch_axis_factor_terms`
  - `_pqs_outer_mismatch_direct_support_oracle_entries`
  - `_pqs_outer_mismatch_direct_support_oracle_block`
  - `_pqs_outer_mismatch_local_column_range`
  - `_pqs_outer_mismatch_product_block`
  - `_pqs_outer_mismatch_safe_term_operator_comparison`
- Atom box:
  - `_PQS_ATOM_BOX_SAFE_TERM_OPERATOR_TERMS`
  - `_PQS_ATOM_BOX_UNSUPPORTED_OPERATOR_TERMS`
  - `_pqs_atom_box_axis_factor_terms`
  - `_pqs_atom_box_local_column_range`
  - `_pqs_atom_box_support_dense_unit_entries`
  - `_pqs_atom_box_direct_support_oracle_entries`
  - `_pqs_atom_box_support_local_block`
  - `_pqs_atom_box_safe_term_operator_comparison`

Live builders preserved:
- `_pqs_contact_cap_product_doside_unit`
- `_pqs_outer_mismatch_product_doside_units`
- `_pqs_atom_box_support_dense_units`

Caller audit:
- Pre-delete `rg` found no nonlocal source/test callers for the deleted names.
- Post-delete `rg` found no deleted symbols in `src`, `test`, or the updated policy docs.
- Preserved retained-unit builders still have active current-route metadata callers in `src/cartesian_contracted_parent_metrics/current_route_metadata_export.jl`:
  - `_pqs_outer_mismatch_product_doside_units(construction; audit = audit)`
  - `_pqs_atom_box_support_dense_units(construction; audit = audit)`
  - `_pqs_contact_cap_product_doside_unit(construction; audit = audit)`

Docs:
- Updated `raw_product_source_retained_transform_policy.md` and `projected_q_shell_policy.md` to say the old `*_safe_term_operator_comparison` wrappers are retired.
- The docs now point to the retained-unit fixture builders plus active current-route metadata export, CPB, and source-box contracts.
- No new design-note layer added.

Validation:
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed
  - one dependency precompiled; output included `GaussletBases` precompile in 56350.8 ms and `load ok`
- `git diff --check`
  - passed
- Current-route metadata smoke:
  - No non-slow current-route metadata/report smoke was available that actually exercises `current_route_metadata_export.jl` without a Be2 artifact directory.
  - The available real-artifact acceptance test is opt-in through `BE2_PQS_Q5_ARTIFACT_DIR`; without that artifact it only skips and does not exercise the path.
  - Used the active-caller audit above instead, per blurb fallback.

Line budget:
- `git diff --numstat`:
  - `19  48  docs/src/developer/projected_q_shell_policy.md`
  - `20  53  docs/src/developer/raw_product_source_retained_transform_policy.md`
  - `0   749 src/cartesian_contracted_parent_metrics/legacy_source_box_fixtures.jl`
- Total: 39 added, 850 deleted, net -811.
- Source deletion: 749 lines.
- Target of at least 500 net deleted met.

Git status:
```text
## main...origin/main
 M docs/src/developer/projected_q_shell_policy.md
 M docs/src/developer/raw_product_source_retained_transform_policy.md
 M src/cartesian_contracted_parent_metrics/legacy_source_box_fixtures.jl
```

Deletion/shrinkage report:
- deleted: the contact-cap, outer-mismatch, and atom-box safe-term/operator comparison wrapper families listed above.
- simplified: policy docs now describe those wrappers as retired and point to retained-unit builders plus active route metadata/CPB/source-box contracts.
- quarantined: none.
- not deleted because: retained-unit builders remain live and are active inputs to `current_route_metadata_export.jl`.
- exact remaining caller/blocker: `src/cartesian_contracted_parent_metrics/current_route_metadata_export.jl` still calls `_pqs_outer_mismatch_product_doside_units`, `_pqs_atom_box_support_dense_units`, and `_pqs_contact_cap_product_doside_unit` for current-route retained-unit inventory.

-- repo-doer@macmini

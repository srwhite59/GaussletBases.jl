# Pass 243 blurb - delete uncalled PQS route-shadow comparison oracles

Role: repo-doer.

Read before editing:

- `AGENTS.md`
- `BlurbStyle.md`
- `docs/src/developer/pqs_manager_running_log.md`
- `docs/src/developer/old_flat_cartesian_retirement_audit_2026-06-14.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.242.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/review.242.md`

Task type: cleanup/deletion after pass-241 line-budget exception.

Purpose:

Delete uncalled legacy comparison/oracle wrappers from the split
`CartesianContractedParentMetrics` legacy fixture file. These wrappers are old
route-shadow validation scaffolds. The current route metadata lane still needs
the retained-unit fixture builders, but it does not need the support-local
operator-comparison oracle families.

Deletion candidate:

```text
file:
  src/cartesian_contracted_parent_metrics/legacy_source_box_fixtures.jl
```

Subagent/read-only audit summary:

```text
risk class:
  green, if only the uncalled comparison/oracle wrapper families are removed.

expected line savings:
  about 700 source lines if all three families are safely deleted.

live builders to preserve:
  _pqs_contact_cap_product_doside_unit(...)
  _pqs_outer_mismatch_product_doside_units(...)
  _pqs_atom_box_support_dense_units(...)

current live source callers:
  src/cartesian_contracted_parent_metrics/current_route_metadata_export.jl
    calls/names those three builders for retained-unit/source metadata.
```

Required deletion target:

Delete the uncalled comparison/oracle families:

```text
contact cap comparison family:
  _PQS_CONTACT_CAP_SAFE_TERM_OPERATOR_TERMS
  _PQS_CONTACT_CAP_UNSUPPORTED_OPERATOR_TERMS
  _pqs_contact_cap_axis_factor_terms(...)
  _pqs_contact_cap_direct_support_oracle_entries(...)
  _pqs_contact_cap_direct_support_oracle_block(...)
  _pqs_contact_cap_safe_term_operator_comparison(...)

outer mismatch comparison family:
  _PQS_OUTER_MISMATCH_SAFE_TERM_OPERATOR_TERMS
  _PQS_OUTER_MISMATCH_UNSUPPORTED_OPERATOR_TERMS
  _pqs_outer_mismatch_axis_factor_terms(...)
  _pqs_outer_mismatch_direct_support_oracle_entries(...)
  _pqs_outer_mismatch_direct_support_oracle_block(...)
  _pqs_outer_mismatch_local_column_range(...)
  _pqs_outer_mismatch_product_block(...)
  _pqs_outer_mismatch_safe_term_operator_comparison(...)

atom box comparison family:
  _PQS_ATOM_BOX_SAFE_TERM_OPERATOR_TERMS
  _PQS_ATOM_BOX_UNSUPPORTED_OPERATOR_TERMS
  _pqs_atom_box_axis_factor_terms(...)
  _pqs_atom_box_local_column_range(...)
  _pqs_atom_box_support_dense_unit_entries(...)
  _pqs_atom_box_direct_support_oracle_entries(...)
  _pqs_atom_box_support_local_block(...)
  _pqs_atom_box_safe_term_operator_comparison(...)
```

Preserve:

- `_pqs_contact_cap_product_doside_unit(...)`;
- `_pqs_outer_mismatch_product_doside_units(...)`;
- `_pqs_atom_box_support_dense_units(...)`;
- current-route metadata export behavior;
- fake-PQS guardrails and independent H2 PQS work.

Before deleting:

- Re-run `rg` for the candidate function/constant names.
- If any candidate has a live nonlocal source caller, stop and write
  `ATTENTION.md` with the exact caller.
- Local calls inside the same soon-to-be-deleted comparison family do not block
  deletion.

Docs:

If developer docs still mention the deleted `*_safe_term_operator_comparison`
helpers, update those lines to say the comparison wrappers were retired and the
current authority is the retained-unit fixture builder plus active route
metadata/CPB/source-box contracts. Do not add a new long design note.

Forbidden:

- do not delete the retained-unit fixture builders listed above;
- do not delete raw product-box route producer or all-pairs inventory machinery;
- do not touch H2 independent source-plan/final-basis/H1 work;
- do not touch fake-PQS endpoint coverage;
- do not rework `current_route_metadata_export.jl` except for a trivial import/
  reference cleanup if deletion makes it necessary;
- do not run CR2, Hamiltonian assembly, H1, H1-J, RHF, full Be2 artifact
  acceptance, or slow old projected-Q-shell integration gates.

Line budget:

- This pass must be net-negative in `src + test + bin`.
- Target at least `500` net deleted source/test/bin lines if the audit is
  accurate.
- If you cannot safely delete at least `300` net lines, write `ATTENTION.md`
  explaining the blocker and exact remaining candidates.

Minimum validation:

Run:

```text
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
```

Also run the smallest current-route metadata/report smoke you can identify that
exercises `current_route_metadata_export.jl` without entering the slow old
integration gates. If the only available smoke is slow or protects the stale
comparison wrappers themselves, do not run it; explain the active-caller check
instead.

Do not run:

- `test/nested/integration_runtests.jl`;
- `test/nested/cartesian_pair_stage_low_order_policy_runtests.jl`;
- old projected-Q-shell integration gates;
- full Be2/CR2/Hamiltonian/HF/RHF routes.

Report:

- exact lines added/deleted in `src + test + bin`;
- function families deleted;
- live builders preserved;
- caller audit result;
- docs updated, if any;
- validation and timings;
- deleted/simplified/quarantined/not-deleted accounting;
- exact remaining blocker/caller, if any.

-- repo-manager@macmini

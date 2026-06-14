# Pass 236 response - independent H2 PQS retained-rule/source-plan seam audit

No source, test, or task-doc files were edited. This was a read-only audit.

## Direct Answers

1. `:atom_contact_core`

The old/fake retained count `251` is not justified by an independent PQS rule
in the inspected code. The generated support region is a union of two
atom-local `5^3` cores plus one `5*5*1` midpoint slab, total support `275`.
Existing terminal lowering classifies these direct regions with direct
identity-like contracts:

- `src/cartesian_terminal_lowering/region_contracts.jl`:
  `_direct_terminal_contract(...)`
- retained rule: `:direct_source_modes`
- retained-unit path:
  `src/cartesian_retained_units/lower_contract_units.jl`:
  `_direct_retained_unit(...)`

That authority would retain all direct source modes, so the honest independent
count from existing direct contracts is `275`, not `251`. A reduced
atom-contact count would require a new private retained rule with explicit
physics/algorithmic justification. The current missing object name is accurate:
`:independent_pqs_atom_contact_core_retained_rule`.

2. `:shared_shell_1`

The standard q=5 PQS boundary product-mode rule does explain retained count
`98`. The exact reusable authority is:

- `src/cartesian_raw_product_sources/records.jl`:
  `CartesianRawProductSources.raw_product_box_plan(...)`
- `src/cartesian_raw_product_sources/records.jl`:
  `CartesianRawProductSources.pqs_boundary_product_mode_retained_rule(...)`

The rule is metadata-only and selects source modes where any axis mode is on
the first or last boundary. For `(5, 5, 5)`, it retains `98` modes. Terminal
lowering already names the same rule for PQS complete shells:

- `src/cartesian_terminal_lowering/region_contracts.jl`:
  `_pqs_complete_shell_contract(...)`
- retained rule: `:pqs_boundary_comx_product_modes`
- source-mode shape: `(5, 5, 5)` for `PQSLowering(q = 5)`

3. `:shared_shell_2`

I found no current independent PQS authority for the old/fake retained count
`114`. The standard q=5 boundary product-mode rule gives `98`, not `114`.
Arithmetically, `114` matches a boundary product-mode count for a non-cubic
`(5, 5, 6)` source shape, but the inspected complete-shell PQS lowering path
does not generate a `(5, 5, 6)` source for `:shared_shell_2`.

Current constraints:

- `_pqs_complete_shell_contract(...)` records cubic `source_mode_shape =
  (q, q, q)`.
- `pqs_multilayer_shell_source_plan(...)` currently requires cubic raw source
  dimensions for shell layers.
- `distorted_product_box_comx` can carry non-cubic `source_mode_shape`, but it
  is attached to distorted product/gap regions, not to the shared molecular
  complete-shell regions produced here.

So the independent route should not claim `(251, 98, 114)`. With existing
authority, the plausible generated retained counts are `(275, 98, 98)`, with
final dimension `471`, or the route should remain blocked until a real
`shared_shell_2` retained rule is defined.

4. Next pass boundary

Do not build a partial `shared_shell_1`-only source plan first. That would add
a materialized fragment while leaving the controlling atom-contact decision
unresolved.

Smallest honest implementation pass:

- introduce a private independent H2 retained-rule/readiness summary that
  derives:
  - `:atom_contact_core => 275` from direct source modes;
  - `:shared_shell_1 => 98` from
    `pqs_boundary_product_mode_retained_rule((5,5,5))`;
  - `:shared_shell_2 => 98` from the same standard q=5 rule, or else keeps
    `shared_shell_2` blocked if manager wants to preserve a blocker instead of
    accepting final dimension `471`;
- update target readiness to stop carrying the old fake `(251,98,114)` as an
  independent-PQS expectation.

This should be metadata/readiness only. No coefficient matrices, final basis,
H1, H1-J, RHF, supplements, or fake/WL transforms.

5. Source-plan ownership

Use a small private H2 route source-plan/readiness object as the owner of the
cross-region contract. It should compose existing lower-level authorities:

- `CartesianShellification`: support-region geometry and terminal regions;
- `CartesianTerminalLowering`: direct versus PQS complete-shell contracts;
- `CartesianRawProductSources`: raw product source-mode facts and
  `pqs_boundary_product_mode_retained_rule(...)`;
- `CartesianRetainedUnits`: compact retained-unit summaries from selected
  lowering contracts.

Do not make `CartesianRawProductSources` own the whole H2 source plan; its own
docstring says it owns raw source facts and optional boundary retained rules,
not shell realization or route assembly. Do not use the current
`pqs_multilayer_shell_source_plan(...)` as the first H2 route authority; it is
valuable precedent, but it assumes a one-core multi-layer complete-shell shape
and cubic shell source dimensions.

Suggested missing private object/function:

```julia
_PQSIndependentH2CoreShellRetainedRulePlan
_pqs_source_box_route_driver_independent_h2_retained_rule_plan(...)
```

Fields should stay compact: status, blocker, support_order, support_counts,
retained_order, retained_counts, expected_final_dimension, per-unit rule
summary, authority/provenance, missing_objects.

6. Remaining deletion/shrink candidates

After pass 235, the cleanest near-term deletion candidates are smaller:

- route-skeleton helper-name/report pressure around `helper_by_pair_family`
  once the independent H2 retained-rule summary stops needing those aliases;
- remaining excess diagnostic assertions in
  `test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl`,
  preserving only compact boundary count, density normalization, and nuclear
  sign checks;
- `current_route_metadata_export.jl` migration scaffolding that is still only
  test/report vocabulary, after caller search proves no route-driver source
  caller remains.

Not deletion targets for the next implementation pass:

- fake-PQS H2 463 golden regression;
- active layered route spine modules;
- lower-level product/PQS density and nuclear convention helpers until their
  remaining source/test callers are replaced by compact module-contract tests.

## Recommended Next Implementation Pass

Add the private retained-rule/readiness plan only, with no materialized source
coefficients:

```text
support:  (275, 578, 362)
retained: (275, 98, 98) if manager accepts direct atom-contact + standard q=5 shells
final_dimension: 471
authority:
  atom_contact_core = :direct_source_modes
  shared_shell_1    = :pqs_boundary_comx_product_modes
  shared_shell_2    = :pqs_boundary_comx_product_modes
```

If manager does not want to accept `471` yet, implement the same plan with
`shared_shell_2` blocked by
`:missing_independent_pqs_shared_shell_2_retained_rule` and explicitly report
that `114` remains fake/WL-derived, not independent PQS.

## Validation

Read-only inspection only. No Julia commands and no tests were run.

Current status:

```text
## main...origin/main [ahead 1]
```

## No-Edit Confirmation

No source, test, or task-doc files were edited during the audit. Only this
baton response will be written.

-- repo-doer@macmini

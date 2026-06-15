# Pass 246 blurb - shrink legacy-default low-order policy vocabulary tests

Role: repo-doer.

Read before editing:

- `AGENTS.md`
- `BlurbStyle.md`
- `docs/src/developer/pqs_manager_running_log.md`
- `docs/src/developer/old_flat_cartesian_retirement_audit_2026-06-14.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.245.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/review.245.md`

Task type: cleanup/test shrink.

Purpose:

Shrink tests that preserve old default `:legacy_diatomic_source` low-order
policy vocabulary. This is a tests-only cleanup pass. The old legacy diatomic
source path may still exist in source, but these tests should stop preserving
large clouds of private default/deferred/not-selected vocabulary as if that
were the active public architecture.

Deletion candidate:

```text
primary files:
  test/nested/cartesian_shell_stage_low_order_policy_runtests.jl
  test/nested/cartesian_unit_stage_low_order_policy_runtests.jl
  test/nested/cartesian_transform_stage_low_order_policy_runtests.jl
  test/nested/cartesian_assembly_stage_low_order_policy_runtests.jl
  test/nested/cartesian_report_stage_low_order_policy_runtests.jl

avoid unless you have a very small safe edit:
  test/nested/cartesian_pair_stage_low_order_policy_runtests.jl
```

Audit summary:

```text
risk class:
  green for tests-only deletion of default-legacy vocabulary assertions when
  source behavior is unchanged.

expected line savings:
  about 150-300 test lines, depending on how much report-stage legacy default
  vocabulary can be safely collapsed.

obsolete pressure:
  :default_legacy_diatomic_source
  :legacy_diatomic_source_low_order_units
  :legacy_diatomic_source_low_order_transforms
  :legacy_diatomic_source_low_order_pairs
  :legacy_diatomic_source_pair_terms
  :legacy_diatomic_source_low_order_assembly
  :deferred_legacy_diatomic_source_*
  :not_selected_legacy_source_pairs

replacement/current authority:
  active atom-growth and terminal-shellification policy tests;
  independent H2 PQS artifact/readiness checks;
  old legacy source path remains reference/migration-only, not public route
  authority.
```

Prioritize:

1. Delete or collapse default-legacy test blocks whose main purpose is asserting
   exact `:legacy_diatomic_source*` default statuses across multiple stages.
2. Keep one compact smoke, if needed, that the old default path is still
   recognized as legacy/reference only.
3. Preserve active atom-growth and terminal-shellification policy coverage.
4. In the report-stage file, avoid expanding the stale broad gate; shrink only
   the legacy-default vocabulary block if it is clearly separate.

Forbidden:

- do not edit source code;
- do not change default driver behavior;
- do not delete active atom-growth or terminal-shellification tests;
- do not touch H2 independent source-plan/final-basis/H1 work;
- do not touch fake-PQS endpoint coverage;
- do not add new tests or new helper layers.

Line budget:

- This pass must be net-negative in `src + test + bin`.
- Target at least `150` net deleted lines.
- If you cannot safely delete at least `80` net lines, write `ATTENTION.md`
  explaining the blocker and exact remaining candidates.

Validation:

Run focused validation for each edited test file if it is expected to complete
under 60 seconds.

Recommended starting set if edited:

```text
julia --project=. test/nested/cartesian_shell_stage_low_order_policy_runtests.jl
julia --project=. test/nested/cartesian_unit_stage_low_order_policy_runtests.jl
julia --project=. test/nested/cartesian_transform_stage_low_order_policy_runtests.jl
julia --project=. test/nested/cartesian_assembly_stage_low_order_policy_runtests.jl
git diff --check
```

For `cartesian_report_stage_low_order_policy_runtests.jl`, do not force the
full stale broad gate. It is already known to fail on unrelated terminal/report
assertions. If you edit only its default-legacy block, use package load plus
caller/diff review and explain why the broad gate is not a valid focused
validation.

Do not run:

- `test/nested/integration_runtests.jl`;
- `test/nested/cartesian_pair_stage_low_order_policy_runtests.jl` as a gate;
- full Be2/CR2/Hamiltonian/HF/RHF routes.

Report:

- exact lines added/deleted;
- files/test blocks deleted or collapsed;
- active policy coverage preserved;
- validation and timings;
- deleted/simplified/quarantined/not-deleted accounting;
- exact remaining blocker/caller.

-- repo-manager@macmini

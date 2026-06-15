Pass 246 complete.

Files changed:
- `test/nested/cartesian_shell_stage_low_order_policy_runtests.jl`
- `test/nested/cartesian_unit_stage_low_order_policy_runtests.jl`
- `test/nested/cartesian_transform_stage_low_order_policy_runtests.jl`
- `test/nested/cartesian_assembly_stage_low_order_policy_runtests.jl`
- `test/nested/cartesian_report_stage_low_order_policy_runtests.jl`

Exact lines added/deleted:
- `git diff --numstat -- src test bin`
  - `0  31  test/nested/cartesian_assembly_stage_low_order_policy_runtests.jl`
  - `0  89  test/nested/cartesian_report_stage_low_order_policy_runtests.jl`
  - `0  16  test/nested/cartesian_shell_stage_low_order_policy_runtests.jl`
  - `0  20  test/nested/cartesian_transform_stage_low_order_policy_runtests.jl`
  - `0  20  test/nested/cartesian_unit_stage_low_order_policy_runtests.jl`
- Total `src + test + bin`: `0` added / `176` deleted, net `-176`.

Files/test blocks deleted or collapsed:
- Collapsed default-legacy blocks in shell, unit, transform, assembly, and report stage tests.
- Deleted exact assertion pressure for:
  - `:default_legacy_diatomic_source`
  - `:legacy_diatomic_source_low_order_units`
  - `:legacy_diatomic_source_low_order_transforms`
  - `:legacy_diatomic_source_low_order_pairs`
  - `:legacy_diatomic_source_pair_terms`
  - `:legacy_diatomic_source_low_order_assembly`
  - `:deferred_legacy_diatomic_source_*`
  - `:not_selected_legacy_source_pairs`
- In each edited default block, kept only compact legacy/reference smoke plus basic staged-data consistency where it was already present.
- Did not edit active atom-growth or terminal-shellification sections.

Active policy coverage preserved:
- Shell stage still checks selected atom-growth policy and terminal-shellification policy sections.
- Unit stage still checks active atom-growth and terminal unit inventory sections outside the collapsed default block.
- Transform stage still checks active atom-growth and terminal transform sections outside the collapsed default block.
- Assembly stage still checks active atom-growth and terminal assembly sections outside the collapsed default block.
- Report stage still keeps atom-growth RouteCore/CRC print coverage and the existing terminal/report coverage outside the collapsed default-legacy block.

Validation and timings:
- `julia --project=. test/nested/cartesian_shell_stage_low_order_policy_runtests.jl`
  - passed: 129 tests, 27.5s reported by Julia.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed, printed `load ok`.
- `git diff --check`
  - passed.
- Legacy-vocabulary `rg` check:
  - no remaining occurrences in the five edited files for `default_legacy_diatomic_source`, `legacy_diatomic_source_low_order*`, `legacy_diatomic_source_pair_terms`, `deferred_legacy_diatomic_source*`, or `not_selected_legacy_source_pairs`.
- Attempted focused gates for edited unit/transform/assembly files:
  - `julia --project=. test/nested/cartesian_unit_stage_low_order_policy_runtests.jl`
  - `julia --project=. test/nested/cartesian_transform_stage_low_order_policy_runtests.jl`
  - `julia --project=. test/nested/cartesian_assembly_stage_low_order_policy_runtests.jl`
  - These exceeded 60s and failed on unrelated stale terminal-shellification exact-field assertions outside the edited default-legacy blocks.
  - Unit failure example: missing `terminal_shellification_units_selected` field on `terminal_units`.
  - Transform failure example: missing `terminal_shellification_transforms_selected` field on `terminal_transforms`.
  - Assembly failure examples: missing `terminal_shellification_assembly_selected` and terminal lowering-contract fields on the current assembly object.
  - Per blurb fallback, I did not broaden this pass to repair those stale terminal gates. I terminated the remaining assembly process with `pkill -f cartesian_assembly_stage_low_order_policy_runtests.jl` after it continued streaming stale failures.
- Did not run `cartesian_report_stage_low_order_policy_runtests.jl`; it is already known stale from pass 244 and this pass only edited its default-legacy block.

Git status:
```text
## main...origin/main
 M test/nested/cartesian_assembly_stage_low_order_policy_runtests.jl
 M test/nested/cartesian_report_stage_low_order_policy_runtests.jl
 M test/nested/cartesian_shell_stage_low_order_policy_runtests.jl
 M test/nested/cartesian_transform_stage_low_order_policy_runtests.jl
 M test/nested/cartesian_unit_stage_low_order_policy_runtests.jl
```

Deletion/shrinkage accounting:
- deleted: exact default legacy/deferred/not-selected vocabulary assertions across five stage-policy tests.
- simplified: default path coverage is now compact legacy/reference smoke plus basic staged-data consistency.
- quarantined: none.
- not deleted because: active atom-growth and terminal-shellification policy sections were outside this pass scope and remain untouched.
- exact remaining blocker/caller: unit, transform, assembly, and report stage broad gates still contain unrelated stale terminal-shellification exact-field assertions; those are the remaining blockers to using those files as per-pass validation gates.

-- repo-doer@macmini

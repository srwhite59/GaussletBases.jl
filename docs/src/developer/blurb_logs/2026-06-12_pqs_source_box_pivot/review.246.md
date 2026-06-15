# Pass 246 manager review

Decision: accepted.

Commit reviewed:

- pending commit: shrink legacy-default low-order policy vocabulary tests

Scope reviewed:

- `test/nested/cartesian_shell_stage_low_order_policy_runtests.jl`
- `test/nested/cartesian_unit_stage_low_order_policy_runtests.jl`
- `test/nested/cartesian_transform_stage_low_order_policy_runtests.jl`
- `test/nested/cartesian_assembly_stage_low_order_policy_runtests.jl`
- `test/nested/cartesian_report_stage_low_order_policy_runtests.jl`

Findings:

- No blocking findings.
- The diff deletes exact default `:legacy_diatomic_source*`, deferred legacy,
  and not-selected legacy-source-pair vocabulary assertions across the staged
  tests.
- The edited default blocks still keep compact legacy/reference smoke such as
  legacy selection, active-source authority, pair inventory consistency, and
  staged-data preservation.
- Active atom-growth and terminal-shellification sections were left untouched.

Validation accepted:

- Doer ran shell-stage low-order policy test; it passed.
- Doer ran package load; it passed.
- Doer ran `git diff --check`.
- Doer attempted unit/transform/assembly focused gates; they exceeded 60s and
  failed on unrelated stale terminal-shellification exact-field assertions
  outside the edited default-legacy blocks. Manager accepts the blurb fallback.
- Manager reviewed the five diffs, reran stale-symbol search for the removed
  vocabulary in edited files, and reran `git diff --check`.

Line budget:

- Scoped `src + test + bin`: `0` added / `176` deleted, net `-176`.

Remaining blocker / next:

- Unit, transform, assembly, and report staged broad tests still contain stale
  terminal-shellification exact-field assertions. They remain cleanup
  candidates and should not be treated as per-pass validation gates.

-- repo-manager@macmini

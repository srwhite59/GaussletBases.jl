# Pass 248 manager review

Decision: accepted.

Commit reviewed:

- pending commit: materialize independent H2 PQS final basis

Scope reviewed:

- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `src/pqs_source_box_route_driver_helpers.jl`
- `src/pqs_source_box_route_driver_reporting.jl`
- staged low-order terminal-shellification assertion deletions

Findings:

- No blocking findings.
- The final-basis retained-count gate now follows the route-owned source-plan
  summary instead of the stale `(251, 98, 114)` tuple.
- The focused route now materializes an independent H2 PQS final basis with
  dimension `471`, retained counts `(275, 98, 98)`, full rank, and final overlap
  identity error about `1.3e-13`.
- The stale `source_backed_adapter = true` success metadata was replaced by a
  source-plan label lookup, preserving the fake/source-backed quarantine.
- Final-basis overlap/rank data are exposed as diagnostics. The pass did not
  introduce generalized-overlap transfer semantics or downstream working `S`.
- H1, H1-J, RHF, supplements, CR2, export, and public API remain off.
- The deletion offset removes more exact terminal-shellification mirror
  assertions from broad staged tests while preserving compact active terminal
  route/scaffold/deferred-materialization smoke.

Validation accepted:

- Doer ran package load; it passed in about 59 seconds.
- Doer ran the focused independent H2 PQS driver/artifact check with
  `run_final_basis=true`; it passed in about 78 seconds.
- Doer ran `git diff --check`; it passed.
- Manager reviewed the source/report/test diffs, reran `git diff --check`, and
  accepted the focused driver validation instead of duplicating the slow run.

Line budget:

- Scoped `src + test + bin`: `87` added / `100` deleted, net `-13`.

Remaining blocker / next:

- Independent H2 PQS has a route-owned source plan and final basis. The next
  physics seam is H1 one-body construction, still without H1-J, RHF,
  supplements, CR2, export, or public API.

-- repo-manager@macmini

# Pass 247 manager review

Decision: accepted.

Commit reviewed:

- pending commit: assemble independent H2 PQS source-plan payload

Scope reviewed:

- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `src/pqs_source_box_route_driver_helpers.jl`
- `src/pqs_source_box_route_driver_reporting.jl`
- stale terminal-shellification assertion deletions in staged low-order tests

Findings:

- No blocking findings.
- The independent H2 PQS path now assembles an available route-owned source-plan
  payload with support counts `(275, 578, 362)`, retained counts `(275, 98, 98)`,
  and final dimension `471`.
- The new source-plan summary and route artifact fields preserve the required
  guardrails: `fake_pqs_enabled = false`,
  `source_backed_fixed_source_oracle_used = false`, and
  `retained_transform_authority = :pqs_source_box_construction`.
- Final basis and physics remain blocked. This is important because the existing
  downstream final-basis helper still carries older physical-gausslet
  assumptions, including the old `(251, 98, 114)` retained-count expectation and
  stale source-backed metadata in the final-basis result path. That is the next
  blocker, not a pass-247 defect.
- The deletion offset removes exact terminal-shellification lowering/selected
  contract mirror assertions from broad staged tests while preserving compact
  terminal route/scaffold/materialization smoke.

Validation accepted:

- Doer ran package load; it passed in about 60 seconds.
- Doer ran the focused independent H2 PQS driver/artifact assertion; it passed
  in about 78 seconds and asserted source-plan availability, fake/source-backed
  guards, and continued final-basis/H1 blocking.
- Doer ran `git diff --check`; it passed.
- Manager reviewed the source/report/test diffs, reran `git diff --check`, and
  accepted the focused driver validation rather than duplicating the slow run.

Line budget:

- Scoped `src + test + bin`: `267` added / `287` deleted, net `-20`.

Remaining blocker / next:

- Next pass should review and correct the independent H2 PQS final-basis seam.
  Do not blindly enable final basis until the helper is updated away from the
  fake/WL-era retained-count assumptions and source-backed labels.

-- repo-manager@macmini

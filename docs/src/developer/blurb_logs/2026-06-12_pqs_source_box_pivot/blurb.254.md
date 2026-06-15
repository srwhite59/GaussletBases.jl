Pass 254 - shrink terminal report-stage low-order alias mirrors

Context:
- Current HEAD should include `aa8b96b4 Add independent H2 PQS private RHF input`.
- Passes 250-253 completed the independent H2 PQS H1-J/private-RHF diagnostic
  sequence and cleaned up the input taxonomy.
- The next accepted pass will trigger the required medium-term checkpoint for
  passes 250-254.
- This pass is deliberately cleanup-only. Do not continue physics, supplements,
  CR2/export, or public API work here.

Target:
- `test/nested/cartesian_report_stage_low_order_policy_runtests.jl`
- Focus on the terminal-shellification report-stage block, roughly the section
  that compares `terminal_summary` and `terminal_report` flat
  `low_order_terminal_shellification_*` aliases.

Task:
Delete stale report-stage alias/mirror assertions that preserve old flat
terminal-shellification vocabulary rather than a live contract.

Good deletion candidates:
- Exact object identity mirrors, for example:
  - `terminal_summary.terminal_shellification_scaffold === ...`
  - `terminal_summary.terminal_shellification_unit_inventory === ...`
  - `terminal_report.low_order_terminal_shellification_scaffold === ...`
  - `terminal_report.low_order_terminal_shellification_unit_inventory === ...`
- Flat alias equality mirrors where the compact summary/count/status assertion
  already remains nearby, for example:
  - `low_order_terminal_shellification_unit_keys`
  - `low_order_terminal_shellification_unit_roles`
  - `low_order_terminal_shellification_unit_kinds`
  - `low_order_terminal_shellification_unit_support_counts`
  - repeated alias-to-summary central count mirrors
- Deferred `low_order_route_core_*` terminal alias/status mirrors in this same
  terminal section, if they only duplicate the compact terminal summary and are
  not protecting an active route-core behavior.

Keep:
- The fact that the terminal report-stage summary exists and is selected.
- Compact terminal counts/statuses that still describe the active smoke:
  - terminal region count;
  - unit count if not redundant after deletion;
  - central gap/midpoint/distorted-product counts;
  - materialization required/status/blocker;
  - no Hamiltonian/operator/pair blocks materialized;
  - pair inventory deferred.
- Atom-growth route-core checks earlier in the file unless you can prove they
  are the same stale terminal alias pressure. The intended target is the
  terminal-shellification flat alias block.

Strict exclusions:
- Do not edit source.
- Do not change driver inputs, H2/PQS physics, RHF, H1-J, supplements, CR2,
  export, public API, or manifest files.
- Do not run broad stale integration gates.
- Do not use `test/nested/cartesian_pair_stage_low_order_policy_runtests.jl` as
  validation.

Validation:
- Run `git diff --check`.
- If the focused report-stage test is expected to complete within 60 seconds,
  run:
  `julia --project=. -e 'using Test; include("test/nested/cartesian_report_stage_low_order_policy_runtests.jl")'`.
- If that test is known or observed to be slow/stale, do not burn minutes on it;
  report why and use the smallest available syntax/load validation instead.
- No package load is required unless source is touched.

Line budget:
- This should be strictly net-negative in `src + test + bin`.
- Do not add replacement assertions for deleted flat aliases.

Report:
- Exact assertion families deleted.
- Any assertions deliberately kept and why.
- Validation command(s) and elapsed time if a Julia test was run.
- Scoped line count for `src + test + bin`.
- Deletion/shrinkage result:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini

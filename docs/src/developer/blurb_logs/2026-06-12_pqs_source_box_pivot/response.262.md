Pass 262 response - pay down support-partition exception by shrinking terminal assembly flat mirrors

Implemented the requested tests-only shrink in `test/nested/cartesian_assembly_stage_low_order_policy_runtests.jl`.

Lines deleted from the test file:
- `test/nested/cartesian_assembly_stage_low_order_policy_runtests.jl`: +6 / -51
- net: -45 lines

Stale assertion families removed:
- duplicated `terminal_shellification_*_available` flat mirrors across assembly and summary;
- exact object identity mirrors for terminal scaffold and unit inventory;
- deferred pair inventory/materialization status aliases duplicated between assembly and summary;
- duplicated matrix/materialization false flags on assembly and summary;
- duplicated pair inventory source/count/family helper mirrors;
- empty helper-by-family mirror checks.

Compact terminal-policy checks retained:
- terminal low-order policy route/source selection remains checked;
- terminal route remains summary/deferred only and not active source authority;
- region count and unit count are still checked against the staged pairs object;
- midpoint slab count, central gap region count, and distorted product box count remain checked;
- unit inventory is still present and its region-unit count is checked;
- owned supports and shellification regions remain asserted non-CPB;
- assembled stage still carries `shells`, `units`, `transforms`, and `pairs`.

Validation:
- `julia --project=. -e 'using GaussletBases; Meta.parseall(read("test/nested/cartesian_assembly_stage_low_order_policy_runtests.jl", String)); println("load and parse ok")'`
  - passed: `load and parse ok`
- `git diff --check`
  - passed.

I did not run the edited stale low-order integration test file as a gate. The pass requested lightweight syntax/package-load validation and explicitly excluded broad stale low-order integration gates.

Scoped line count for `src + test + bin`:
- `test/nested/cartesian_assembly_stage_low_order_policy_runtests.jl`: +6 / -51
- total: +6 / -51, net -45

Git status:
```text
## main...origin/main
 M test/nested/cartesian_assembly_stage_low_order_policy_runtests.jl
```

Deletion/shrinkage result:
- deleted: 51 lines of stale terminal-shellification flat mirror assertions.
- simplified: terminal policy coverage now emphasizes route selection, deferred status, compact inventory facts, non-CPB support facts, and major stage object preservation.
- quarantined: no source, driver input, support-partition/provider, running-log, or broad low-order policy coverage changes.
- not deleted because: compact live checks for terminal route selection, physical inventory facts, non-CPB supports, and major stage objects are still part of the pass contract.
- exact remaining caller/blocker: mature terminal-shellification implementation still exposes flat assembly/report mirrors; this pass only reduces test pressure on those mirrors and does not remove source fields.

-- repo-doer@macmini

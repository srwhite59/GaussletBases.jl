# Pass 254 response - shrink terminal report-stage low-order alias mirrors

Cleanup-only pass completed. I edited only:

```text
test/nested/cartesian_report_stage_low_order_policy_runtests.jl
```

No source, driver inputs, manifest, H2/PQS physics, RHF, H1-J, supplements,
CR2/export, or public API files were changed.

## Deleted assertion families

Deleted stale terminal-shellification report-stage alias/mirror assertions:

```text
terminal_summary.terminal_shellification_scaffold === staged scaffold
terminal_summary.terminal_shellification_unit_inventory === staged inventory
terminal_report.low_order_terminal_shellification_scaffold === summary scaffold
terminal_report.low_order_terminal_shellification_unit_inventory === summary inventory
terminal_report.low_order_terminal_shellification_unit_keys == summary keys
terminal_report.low_order_terminal_shellification_unit_roles == summary roles
terminal_report.low_order_terminal_shellification_unit_kinds == summary kinds
terminal_report.low_order_terminal_shellification_unit_support_counts == summary support counts
terminal_report low_order central gap/midpoint/distorted-product count mirrors
terminal_report.low_order_route_core_* deferred/status/count mirrors
```

## Deliberately kept

Kept compact active smoke for:

```text
terminal report summary exists and is selected
terminal region and unit counts
terminal scaffold/inventory availability
pair inventory deferred
assembly/materialization required/status/blocker
no Hamiltonian/operator/pair blocks materialized
central gap/midpoint/distorted-product counts
terminal print-line content
```

Atom-growth RouteCore checks earlier in the file were not touched.

## Line budget

Scoped `src + test + bin`:

```text
0 insertions, 36 deletions
net -36
```

Deletion/shrinkage report:

```text
deleted:
- stale terminal report-stage flat alias and object-identity mirrors
- stale terminal report-stage low-order RouteCore deferred/status mirrors

simplified:
- terminal report-stage coverage now leans on compact selected/available/count/
  deferral smoke instead of exact private alias clouds

quarantined:
- terminal-shellification report-stage aliases remain test-only compatibility
  surface, not physics or route authority

not deleted because:
- compact terminal summary/report smoke and print-line checks still protect the
  active visible report behavior
- atom-growth RouteCore checks are outside this terminal cleanup target

exact remaining caller/blocker:
- remaining terminal aliases are still referenced by compact smoke and print
  checks in this file
```

## Validation

Commands run:

```text
julia --project=. -e 'Meta.parseall(read("test/nested/cartesian_report_stage_low_order_policy_runtests.jl", String)); println("parse ok")'
git diff --check
git status --short --branch
```

Results:

```text
parse smoke: passed
git diff --check: passed
```

I did not run the full report-stage test. This loop has already observed and
recorded that broad report-stage low-order policy tests are slow/stale and not
appropriate as per-pass gates for these terminal alias cleanup passes.

Final git status:

```text
## main...origin/main
 M test/nested/cartesian_report_stage_low_order_policy_runtests.jl
```

-- repo-doer@macmini

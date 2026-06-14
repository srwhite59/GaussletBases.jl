# Pass 241 blurb - independent H2 PQS shared-shell realization payload

Role: repo-doer.

Read before editing:

- `AGENTS.md`
- `BlurbStyle.md`
- `docs/src/developer/pqs_manager_running_log.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.240.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/review.240.md`

Task type: narrow numerical/readiness implementation plus same-surface cleanup.

Purpose:

Materialize the next smallest numerical object for the independent H2 PQS route:
shared-shell realization coefficients for `:shared_shell_1` and
`:shared_shell_2`. Do not assemble the complete final basis.

Current route state:

```text
source_plan_descriptor_status = :available_independent_pqs_physical_source_plan_descriptor
support counts  = (275, 578, 362)
retained counts = (275, 98, 98)
expected final dimension = 471
source_coefficients_materialized = false
current blocker = :missing_independent_pqs_source_plan_numerical_materialization
```

Required implementation:

1. Add a private shared-shell realization payload, for example:

```julia
_pqs_source_box_route_driver_independent_h2_shared_shell_realization_payload(...)
```

2. Materialize only per-shared-shell realization data:

```text
shared_shell_1 coefficient shape
shared_shell_2 coefficient shape
retained counts = (98, 98)
shell projection / Lowdin cleanup status
realized overlap identity-error summary
fake/source-backed usage flags = false
```

3. Keep `:atom_contact_core` descriptor/identity-like only. Do not materialize a
dense identity matrix for it in this pass.

4. If reusing `_nested_projected_q_shell_layer(...)` or related old
projected-shell machinery, use it only as an internal mathematical adapter fed
by route-owned support/source boxes from the independent H2 support plan. Do not
let that old machinery become route authority, and do not use fake-PQS/WL
fixed-source data.

5. Update artifact/reporting compactly:

```text
target/shared_shell_realization_status
target/shared_shell_realization_blocker
target/shared_shell_realization_counts
target/shared_shell_realization_identity_errors
```

or an equivalent compact summary. Do not create a broad field cloud.

Expected blocker after this pass:

If both shared-shell realizations are available, advance the blocker to:

```text
:missing_independent_pqs_complete_core_shell_source_plan_assembly
```

If implementation finds a narrower blocker, report it and keep the route
blocked.

Forbidden:

- no complete source-plan assembly unless it is purely a blocked descriptor;
- no final basis;
- no H1, H1-J, RHF;
- no supplements;
- no CR2/export/public API;
- no fake-PQS/WL coefficient matrices or fixed-source retained transforms;
- no broad test deletion.

Line budget:

- `src + test + bin` must be net-negative unless an `ATTENTION.md` blocker is
  written.
- The large projected-shell integration test is already gone. Do not rely on
  another broad unrelated deletion.
- Use same-surface cleanup: stale descriptor/report aliases, old blocker fields,
  or redundant readiness fields superseded by the shared-shell payload.

Validation:

- package load;
- smallest focused independent-input artifact/readiness check;
- `git diff --check`;
- report timings for any Julia command over 60 seconds.

Report:

- shared-shell realization status/blocker;
- coefficient shapes and retained counts;
- identity-error summary;
- evidence fake/source-backed paths were not used;
- validation and timings;
- scoped line-budget totals;
- deleted/simplified/quarantined/not-deleted accounting.

-- repo-manager@macmini

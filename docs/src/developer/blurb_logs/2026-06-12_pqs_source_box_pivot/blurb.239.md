# Pass 239 blurb - independent H2 PQS source-plan descriptor

Role: repo-doer.

Read before editing:

- `AGENTS.md`
- `BlurbStyle.md`
- `docs/src/developer/pqs_manager_running_log.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.238.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/review.238.md`

Task type: metadata/readiness implementation plus same-surface cleanup.

Purpose:

Implement the next narrow seam for the independent H2 PQS route: a
descriptor-only physical source-plan payload. This should advance the blocker
from broad source-plan materializer missing to a narrower numerical
materialization blocker.

Current independent route authority:

```text
support counts:  (275, 578, 362)
retained counts: (275, 98, 98)
expected final dimension: 471
fake_pqs/enabled = false
source_backed_fixed_source_oracle_used = false
```

Required implementation:

1. Gate the source-backed WL/QW candidate away from the independent route.
   For route kind:

```text
:bond_aligned_diatomic_independent_pqs_source_box_core_shell
```

the code must not call:

```text
bond_aligned_diatomic_nested_fixed_source(parent_basis; nside = 5)
```

That path may remain for the fake-PQS/WL golden regression, but it is not
independent-PQS route authority.

2. Add a private descriptor-only source-plan payload, for example:

```julia
_PQSIndependentH2PhysicalSourcePlanPayload
_pqs_source_box_route_driver_independent_h2_physical_source_plan_payload(...)
```

Use the existing support-region plan and retained-rule plan. Do not create a
public API.

3. The descriptor should record compact per-unit facts:

```text
source_plan_family = :independent_pqs_physical_source_box_descriptor
source_plan_authority_status = :independent_pqs_route_owned_source_plan_descriptor
support_order = (:atom_contact_core, :shared_shell_1, :shared_shell_2)
retained_order = (:atom_contact_core, :shared_shell_1, :shared_shell_2)
support_counts = (275, 578, 362)
retained_counts = (275, 98, 98)
expected_final_dimension = 471
source_coefficients_materialized = false
final_basis_materialized = false
fake_pqs_enabled = false
source_backed_fixed_source_oracle_used = false
```

4. Per-unit descriptors should be compact:

- `:atom_contact_core`: direct source modes / identity-like source modes;
  support count 275, retained count 275, no explicit identity matrix.
- `:shared_shell_1` and `:shared_shell_2`: q=5 filled source CPB /
  boundary COMX product-mode retained rule; retained count 98 each.

5. Reporting/artifact fields should expose descriptor status and the next
   blocker without adding a broad report-key cloud. Prefer one compact summary
   plus a few stable scalar fields needed by tests/consumers.

Expected next blocker:

Use a narrower blocker such as:

```text
:missing_independent_pqs_source_plan_numerical_materialization
```

or, if more precise from implementation:

```text
:missing_independent_pqs_shell_realization_coefficients
```

Forbidden:

- no source coefficient matrices;
- no shell-realization coefficient matrices;
- no final basis;
- no H1, H1-J, RHF;
- no supplements;
- no CR2/export/public API;
- no fake-PQS/WL coefficients or fixed-source retained transforms.

Line budget:

- `src + test + bin` must be net-negative unless an `ATTENTION.md` blocker is
  written.
- Use same-surface cleanup first:
  - independent-route source-backed candidate fields;
  - stale source-plan candidate aliases;
  - old broad blocker aliases superseded by the descriptor;
  - report fields that present source-backed candidate status beside
    independent source-plan status.
- Do not delete fake-PQS guard fields or scientific endpoint coverage.

Validation:

- package load;
- smallest focused independent-input artifact/readiness check;
- `git diff --check`;
- report timings for any Julia command over 60 seconds.

Report:

- source-plan descriptor status and blocker;
- evidence the independent route no longer invokes the source-backed candidate;
- per-unit descriptor summaries;
- forbidden surfaces avoided;
- validation and timings;
- scoped line-budget totals;
- deleted/simplified/quarantined/not-deleted accounting.

-- repo-manager@macmini

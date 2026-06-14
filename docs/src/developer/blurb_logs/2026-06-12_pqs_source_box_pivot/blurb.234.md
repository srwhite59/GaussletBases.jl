# Pass 234 blurb - implement independent H2 PQS support-region plan

Role: repo-doer.

Read before editing:

- `AGENTS.md`
- `BlurbStyle.md`
- `docs/src/developer/pqs_manager_running_log.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.232.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.233.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/review.233.md`

Task:

Implement the independent H2 PQS support-region materializer/fingerprint only.
This pass should move the independent target from target support constants to
generated support-region authority if, and only if, the route geometry can
generate the expected support units.

Current route:

```julia
:bond_aligned_diatomic_independent_pqs_source_box_core_shell
```

Expected generated support grouping:

```text
:atom_contact_core = two atom-local 5^3 cores + 5*5*1 midpoint/contact slab
                   = 125 + 125 + 25 = 275

:shared_shell_1    = outer shared molecular shell
                   = 9*9*15 - 7*7*13 = 578

:shared_shell_2    = inner shared molecular shell
                   = 7*7*13 - 5*5*11 = 362
```

The shared-shell order must be outside-in, even if raw terminal geometry emits
inside-out.

Allowed implementation:

- Add a compact private support-plan helper/object near the independent target
  payload path, for example:

```julia
_PQSIndependentH2PhysicalSupportRegionPlan
_pqs_source_box_route_driver_independent_h2_support_region_plan(...)
```

- Use existing route geometry/shellification/lowering support APIs where
  possible:
  - `CartesianShellification.raw_terminal_geometry(...)`
  - `CartesianShellification.shellify(...)`
  - `CartesianRouteCore.OwnedSupport` / support summaries
  - current equivalent APIs if names differ
- Group primitive terminal regions into the three route support units above.
- Record support order, generated support counts, parent-axis counts, bond
  metadata, coverage/disjointness status, and authority/provenance.

Required success semantics:

- `target/support_plan_status = :available_independent_pqs_support_region_plan`
  or a similarly explicit available status.
- `target/support_plan_authority = :cartesian_shellification_route_geometry`
  or similarly explicit route-geometry/shellification authority.
- `target/support_counts_generated = true`.
- `target/support_counts_source` must no longer be
  `:target_constants_pending_support_region_materializer`.
- Support counts must be `(275, 578, 362)` in support order.
- Guard fields remain:
  - `fake_pqs/enabled = false`
  - `source_backed_fixed_source_oracle_used = false`
  - `physics/endpoint_ready = false`
- Retained counts remain empty/blocked/unknown.
- Expected final dimension remains absent/unknown.

If implementation is blocked:

- Do not add a second hard-coded source of truth.
- Keep the blocked support-plan status.
- Report the exact missing geometry/API field, not a generic blocker.

Forbidden:

- no atom-contact retained rule;
- no `shared_shell_2` retained rule;
- no source-plan authority beyond support-region generation;
- no fake-PQS/WL coefficient matrices or fixed-source retained transforms;
- no final basis, H1, H1-J, RHF, supplements, CR2, export, public API;
- no WL/QW physics comparison.

Line-budget/deletion requirement:

- Keep `src/ + test/ + bin/` net-negative.
- Delete the pass-232 hard-coded blocked `support_plan` field cloud once the
  generated plan lives in the target payload.
- Remove duplicate support-count/blocker aliases that become derivable from the
  support-plan summary.
- Do not delete accepted scientific endpoint tests or fake-PQS guard fields.
- If honest net-negative implementation is blocked, write `ATTENTION.md` and
  stop.

Validation:

- Run the smallest focused independent-input artifact/readiness check.
- Run package load.
- Run `git diff --check`.
- Report Julia-level timings, especially if the driver check exceeds 60s.

Report:

- support-plan status, authority, counts, and coverage/disjointness result;
- whether counts are generated or still target constants;
- exact source APIs used;
- exact forbidden surfaces avoided;
- validation and timings;
- scoped line-budget totals;
- deleted/simplified/quarantined/not-deleted carrying-cost accounting.

-- repo-manager@macmini

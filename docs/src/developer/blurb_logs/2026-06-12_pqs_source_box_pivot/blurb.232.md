# Pass 232 blurb - independent H2 PQS support plan authority

Role: repo-doer.

Read before editing:

- `AGENTS.md`
- `BlurbStyle.md`
- `docs/src/developer/pqs_manager_running_log.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.230.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.231.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/review.231.md`

Task:

Advance the new independent H2 PQS route by one layer only: turn the support
metadata from pass 231 into a route-owned support/region plan, or report the
exact blocker that prevents doing so.

This is still not a retained-transform or physics pass.

Current target:

```text
input:
  test/driver_inputs/h2_pqs_q5_independent_source_box_r4.jl

route kind:
  :bond_aligned_diatomic_independent_pqs_source_box_core_shell

intended support vocabulary:
  :atom_contact_core => 275
  :shared_shell_1    => 578
  :shared_shell_2    => 362
```

Allowed work:

- Add a compact private support-plan object or summary if the existing route
  skeleton/report style needs one.
- Use route geometry, parent axes, shellification/lowering-backed region facts,
  or existing PQS multilayer/terminal-lowering machinery to produce support
  regions.
- Record support order, support counts, parent-axis counts, center/bond
  metadata, disjointness/coverage status, and support-plan authority.
- If the support plan cannot be generated honestly, keep it blocked and record
  a precise blocker such as
  `:missing_independent_pqs_support_region_materializer`.

Required guard fields:

- `fake_pqs/enabled = false`
- `source_backed_fixed_source_oracle_used = false`
- support-plan authority must not be WL/QW fixed-source coefficient data
- `physics/endpoint_ready = false`
- no retained counts except empty/blocked/unknown
- no expected final dimension unless a real retained transform exists

Forbidden in this pass:

- no atom-contact retained rule;
- no `shared_shell_2` retained rule;
- no source-plan authority that imports fake-PQS/WL coefficient matrices;
- no final-basis materialization;
- no H1, H1-J, RHF, supplement, CR2, export, public API;
- no comparison to WL/QW physics values;
- no mutation of the fake-PQS route except shared-helper reuse that preserves
  its quarantine fields.

Deletion/line-budget target:

- Keep `src/ + test/ + bin/` net-negative again.
- Measure with `git diff --numstat -- src test bin`, and count any new driver
  input/test files.
- First deletion candidates should be redundant readiness/report default clouds,
  duplicate blocker aliases, or stale assertions in the same target/readiness
  surface. Do not delete accepted scientific endpoint tests or fake-PQS guard
  fields.
- If support-plan work cannot stay net-negative without unsafe deletion, write
  `ATTENTION.md` and stop.

Validation:

- Run the smallest focused driver/readiness check for the independent input.
- Run package load.
- Run `git diff --check`.
- Report Julia-level timings for any command over 60 seconds.

Report:

- support-plan status and blocker;
- whether support counts are generated or still target constants;
- authority/provenance fields;
- exact forbidden surfaces avoided;
- validation and timings;
- scoped line-budget totals;
- deleted/simplified/quarantined/not-deleted carrying-cost accounting.

-- repo-manager@macmini

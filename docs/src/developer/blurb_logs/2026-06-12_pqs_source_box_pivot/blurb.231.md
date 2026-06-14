# Pass 231 blurb - independent H2 PQS target/readiness surface

Role: repo-doer.

Read before editing:

- `AGENTS.md`
- `BlurbStyle.md`
- `docs/src/developer/pqs_manager_running_log.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.230.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/review.230.md`

Task:

Create the first separate independent-H2-PQS target/readiness surface. This is
not a physics materialization pass.

The target should be distinct from the fake-PQS H2 463 source-backed WL/QW
reproduction. It should make future work start from an explicitly
`fake_pqs=false` route whose retained-transform authority is intended to be PQS
source-box construction.

Add or expose a small driver input:

```text
test/driver_inputs/h2_pqs_q5_independent_source_box_r4.jl
```

Use a route kind such as:

```julia
:bond_aligned_diatomic_independent_pqs_source_box_core_shell
```

If a slightly different name fits the existing route-kind style better, use it,
but the name must include independent/PQS/source-box meaning and must not reuse
the fake-PQS route kind.

Required artifact/readiness semantics:

- `fake_pqs/enabled = false`
- `source_backed_fixed_source_oracle_used = false`
- `retained_transform_authority = :pqs_source_box_construction`
- `physics/endpoint_ready = false`
- primary blocker:
  `:missing_independent_pqs_atom_contact_core_retained_rule`
- secondary blocker:
  `:missing_independent_pqs_shared_shell_2_retained_rule`
- broader source-plan blocker:
  `:missing_independent_pqs_physical_source_plan_materializer`

Record the common H2 physical support vocabulary as target metadata only:

```text
support order/counts:
  :atom_contact_core => 275
  :shared_shell_1    => 578
  :shared_shell_2    => 362
```

Do not claim the fake-PQS/WL retained counts `(251, 98, 114)` as independently
generated PQS retained counts. It is okay to report:

- `shared_shell_1` has a plausible q=5 boundary retained-rule count of `98`;
- `atom_contact_core` retained rule is missing;
- `shared_shell_2` retained rule is missing for the fake-route count `114`.

But keep those as readiness facts, not final source-plan authority.

Testing:

- Add one thin explicit test only if needed to protect this new readiness
  surface. Keep it out of the default runner.
- The test should run or inspect the driver artifact/readiness contract and
  assert the guard fields above.
- Do not reconstruct the fake-PQS route and do not run H1/H1-J/RHF.

Forbidden in this pass:

- no final-basis materialization;
- no H1, H1-J, RHF, supplements, CR2, export/public API work;
- no provider-block work;
- no mutation of the fake-PQS route except harmless shared-helper reuse;
- no comparison to WL/QW physics values;
- no broad scalar report-key cloud if a compact target/readiness object or
  summary fits the existing style.

Line-budget rule:

- Keep `src/ + test/ + bin/` net-negative for this pass.
- Measure with `git diff --numstat -- src test bin`.
- Do not satisfy the rule by deleting accepted scientific endpoint tests,
  fake-PQS quarantine guard fields, or user-facing docs.
- If this target/readiness surface cannot be added while staying net-negative,
  stop and write `ATTENTION.md` explaining the exact blocker and the deletion
  candidates you found.

Validation:

- Run the smallest focused validation for the new readiness surface.
- Run package load.
- Run `git diff --check`.
- Report timings for any Julia command expected to exceed 60 seconds.

Report:

- files changed;
- new route/input name;
- artifact/readiness fields added;
- exact blocker values;
- validation run and timing;
- `git diff --numstat -- src test bin` totals;
- deleted/simplified/quarantined/not-deleted carrying-cost accounting.

-- repo-manager@macmini

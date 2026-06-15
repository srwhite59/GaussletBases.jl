# Pass 252 manager review

Decision: accepted.

Commit reviewed:

- pending commit: materialize independent H2 PQS private RHF diagnostic

Scope reviewed:

- `src/pqs_source_box_route_driver_helpers.jl`
- `src/pqs_source_box_route_driver_reporting.jl`
- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `src/pqs_multilayer_complete_core_shell_rhf.jl`
- `test/nested/cartesian_assembly_stage_low_order_policy_runtests.jl`

Findings:

- No blocking findings.
- The pass keeps the private RHF work on the independent H2 PQS route, with
  `fake_pqs_enabled=false`, `source_backed_fixed_source_oracle_used=false`, and
  `retained_transform_authority=:pqs_source_box_construction`.
- The new role classifier is compact and correctly covers the independent H2
  readiness, final-basis, H1, and H1-J diagnostic artifact roles.
- The RHF input-contract and shared RHF helper now accept the independent route
  kind only through that route/role gate.
- The test deletion removes stale terminal RouteCore mirror assertions while
  retaining compact terminal assembly/scaffold/count checks.

Validation accepted:

- Doer ran package load; it passed.
- Doer ran the focused independent H2 PQS private RHF driver/artifact probe; it
  passed in about 84 seconds.
- Doer ran `git diff --check`; it passed.
- Manager reviewed the diff, reran package load, and accepted doer validation
  without duplicating the slow focused RHF route run.

Line budget:

- Scoped `src + test + bin`: `29` added / `31` deleted, net `-2`.

Remaining blocker / next:

- The independent H2 PQS private RHF diagnostic is now materialized, but the
  route remains private/diagnostic-only. Supplements, CR2/export, public API,
  and public solver readiness remain blocked.

-- repo-manager@macmini

Purpose:
Record the PQS source-box fixture role policy as a tracked docs-only note.

Why now:
Pass 103 clarified that compact H1/H1-J fixtures should not become physics
acceptance gates by inertia. Before any RHF input-contract implementation, we
need a short tracked policy that distinguishes route smoke, convention
diagnostic, oracle/debug, and physics endpoint fixtures.

Governing framework:
Use `docs/src/developer/pqs_source_box_operator_framework.md`,
`docs/src/developer/successor_handoff_2026-06-12_pqs_source_box_pivot.md`, and
the pass 103 response.

Loop rule:
Do not request interactive approval/escalation during the baton loop. Use
approved commands and focused docs validation. If approval would be required,
write `ATTENTION.md` with the exact command, reason, and blocker, then stop.

Exact task:
Add a short tracked fixture policy note, preferably:

`docs/src/developer/pqs_source_box_fixture_policy.md`

Keep it concise and actionable. Include:

1. Fixture roles:
   - route smoke;
   - convention diagnostic;
   - oracle/debug;
   - physics endpoint.
2. Current compact fixture facts:
   - tracked Z=1 H1 seam;
   - direct structured H1/J convention probe;
   - one-center source-box driver H1/J dry-run.
3. Nonclaims:
   - compact H1/J route materialization is not RHF readiness;
   - `final_dimension == 223` is not physics acceptance;
   - self-Coulomb alone is not endpoint validation;
   - shell/support-row and fixed-block paths are oracle/debug, not route
     authority.
4. Parameter families that must move together:
   `Z`, electron count, charge/spin, spacing/distortion/radius, `q`, `n_s`,
   core/source side, shell depth, retained rule, Coulomb expansion, density
   gauge, and fixture family.
5. Before-RHF policy:
   - decide route-smoke-only versus physics endpoint;
   - if physics endpoint, choose He/Z/electron count/closed-shell fixture size
     and reference/error/timing thresholds;
   - if route smoke, assert only compact status/non-promotion facts.

Optionally add a one-line pointer from
`docs/src/developer/pqs_source_box_operator_framework.md` to the new note, but
do not do a broad doc rewrite.

Do not edit source or tests.
Do not implement RHF.
Do not add fixtures.
Do not run Julia.
Do not add GTO, IDA/MWG, exports, artifacts, fixture promotion, or production
route behavior.

Validation:

`git diff --check`
`LC_ALL=C rg -n "[^[:ascii:]]" docs/src/developer/pqs_source_box_fixture_policy.md docs/src/developer/pqs_source_box_operator_framework.md`
`git status --short --branch`

Report back:

- files changed;
- brief summary of the policy;
- whether a framework pointer was added;
- validation commands/results;
- git status;
- deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

After writing `.agent_handoffs/response.104.md`, continue polling for
`blurb.105.md`, `STOP.md`, or `ATTENTION.md`.

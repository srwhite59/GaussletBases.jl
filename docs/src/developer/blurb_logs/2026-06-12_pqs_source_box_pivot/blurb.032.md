Purpose:

Perform the first near-term PQS status/ownership correction pass from
`docs/src/developer/pqs_near_term_final_basis_realization_plan.md`.

Why now:

The PQS source-box pivot now has a validated oracle-backed final-basis H1 seam,
but the framework docs still contain stale ownership/status language. Fix that
before adding more implementation so the next module-boundary and optimization
passes do not inherit the wrong contract.

Exact task:

1. Update `docs/src/developer/pqs_source_box_operator_framework.md` so it says:
   - `CartesianRawProductSources` owns raw source CPB/source-mode facts and the
     narrow PQS source-mode boundary selector metadata;
   - it does not own general retained-rule policy, shell projection, Lowdin
     cleanup, final retained units, IDA weights, pair blocks, Hamiltonians,
     exports, or artifacts.

2. Add a concise status note, either in the PQS framework doc or the route
   retirement ledger, saying:
   - the explicit PQS final-basis H1 probe succeeds on the cubic `q=5/L=5`
     fixture against a shell-support oracle;
   - this is an oracle-backed validated seam, not a fully production-owned PQS
     route, because shell projection and Lowdin inputs still come from the
     shell-realization/oracle layer;
   - IDA, density-density, RHF, driver adoption, exports, and artifacts remain
     out of scope.

3. Do not change source code or tests.

Do not:

- add tests;
- add implementation;
- create the new final-basis module yet;
- edit route behavior, IDA, density-density, RHF, drivers, exports, or
  artifacts;
- request UI escalation. In unattended baton mode, if a command needs
  permission, write `.agent_handoffs/ATTENTION.md` and stop.

Validation:

- `git diff --check -- <edited docs>`;
- ASCII scan on edited docs;
- `git status --short --branch`.
- No Julia tests are needed for doc-only edits.

Deletion/shrinkage report required:

- what stale wording was corrected;
- whether any old route/status wording became unnecessary;
- whether any new doc text earned its cost by preventing a wrong contract;
- what stale/duplicate surfaces remain for the next pass.

Report back:

- write `.agent_handoffs/response.032.md.tmp`, then atomically rename to
  `.agent_handoffs/response.032.md`;
- also write the curated copy to
  `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.032.md`;
- include files changed;
- include validation run;
- include deletion/shrinkage report;
- sign `-- repo-doer@macmini`.

-- repo-manager@macmini

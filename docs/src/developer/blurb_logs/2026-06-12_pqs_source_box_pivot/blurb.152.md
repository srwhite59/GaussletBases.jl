Pass 152 - audit Be2 PQS Ham handoff contract

Purpose:

Do a no-edit audit for the current Be2/PQS readiness blocker:

```text
:missing_diatomic_ham_consumer_contract
```

The route now has private source plan, final basis, H1, and Ham-input/density
interaction data. The next question is what object should be built so a
downstream CR2 agent can eventually inspect or compare Be2 WL versus PQS.

Do not implement the handoff yet. Define the contract first.

Read first:

- `src/pqs_source_box_route_driver_helpers.jl`
- `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
- `docs/src/reference/export.md`
- `src/fullida_dense_export.jl`
- `docs/src/developer/pqs_source_box_operator_framework.md`
- `docs/src/developer/pqs_source_box_fixture_policy.md`
- `~/Dropbox/chatarchive/handoff/software_packets/gausslet-methods-fundamentals.md`

Optional local downstream inspection:

- If readable without approval, inspect only small docs or entry files under:
  - `~/Dropbox/codexhome/work/cr2`
  - `~/Dropbox/codexhome/work/hfdmrg`
- Do not run downstream jobs.
- Do not edit downstream repos.
- If sandbox/read permissions block this, report that and continue from
  GaussletBases-local export/reference docs.

Audit questions:

1. What minimal object should Be2/PQS expose next for downstream comparison?
   Options to evaluate:
   - private dense handoff object with `H1` and a dense/factorized
     electron-electron representation;
   - private factor handoff object carrying `H1`, density interaction, support
     weights, raw pair factors, and ordering metadata;
   - compatibility adapter toward an existing export/HamV6-like shape;
   - blocked payload with sharper missing export/handoff convention.

2. What exactly does the current Ham-input payload already provide?
   Include:
   - final basis/final dimension
   - H1/final Hamiltonian reference
   - pre-final density interaction
   - density gauge
   - raw-pair convention
   - support/pre-final ordering
   - nuclear metadata availability

3. What is still missing before CR2 should consume it?
   Pay attention to:
   - whether a dense `Vee` is required or a factor/pre-final interaction is
     acceptable;
   - how electron-electron data maps to the solver basis/order;
   - whether nuclear repulsion / center metadata must be included;
   - whether electron count/charge metadata belongs here;
   - whether the handoff should be private inspect-only or exportable.

4. What should the next implementation pass build?
   Prefer the smallest private object that:
   - makes the Be2/PQS Hamiltonian-constructor output inspectable;
   - preserves all ordering/convention labels;
   - does not claim public/export/HamV6/hfdmrg/CR2 readiness prematurely.

5. Propose exact status/blocker names.

Trust boundary:

- No source edits.
- No commits.
- No H1-J scalar diagnostic.
- No full public export.
- No RHF/SCF/DIIS work.
- No WL payload implementation.
- No CR2 or hfdmrg execution.
- No edits outside this repo.
- Do not promote shell/support-row contraction, raw product-box probes, or old
  WL adapter paths to route authority.
- Do not reinterpret retained diagnostic/self-integral weights as
  IDA/quadrature weights.
- Do not ask for interactive approval during unattended baton work. If approval
  would be required, write `.agent_handoffs/ATTENTION.md` and stop.

Validation:

- Read-only audit only.
- `git status --short --branch`
- No Julia commands are required.

Report back:

- Recommended handoff object shape.
- What current Be2/PQS payloads already satisfy.
- What remains missing or intentionally nonclaimed.
- Recommended next implementation pass.
- Proposed status/blocker labels.
- Git status.
- Deletion/shrinkage forecast:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini

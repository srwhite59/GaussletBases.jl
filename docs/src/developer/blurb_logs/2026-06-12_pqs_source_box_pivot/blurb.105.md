Purpose:
Add the first private complete core/shell RHF input-contract object.

Why now:
The fixture policy is documented, and pass 102 accepted the RHF boundary:
RHF must be separate from H1/J. The first code step should only validate and
summarize RHF inputs and occupation policy. It must not build Fock matrices,
run SCF, or wire the driver.

Governing framework:
Use `docs/src/developer/pqs_source_box_operator_framework.md` and
`docs/src/developer/pqs_source_box_fixture_policy.md`.

Keep these boundaries sharp:

- source-box-first PQS is the algorithmic framing;
- shell/support-row contraction is oracle/debug;
- retained diagnostic weights are not IDA/quadrature weights;
- compact fixtures are route-smoke/convention diagnostics unless explicitly
  promoted;
- RHF starts as private diagnostic/prototype only.

Loop rule:
Do not request interactive approval/escalation during the baton loop. Use
approved commands and focused validation. If approval would be required, write
`ATTENTION.md` with the exact command, reason, and blocker, then stop.

Exact task:
Add a private RHF input-contract helper/object only.

Suggested surface:

- new file if local include style fits:
  `src/pqs_multilayer_complete_core_shell_rhf.jl`
- include it near `src/pqs_multilayer_complete_core_shell_h1.jl`;
- no export and no public API claim.

Suggested private helper name:

`_pqs_multilayer_complete_core_shell_rhf_input_contract(...)`

It should validate and summarize:

- available source plan;
- available complete core/shell final basis;
- materialized complete core/shell H1 payload/final H1 matrix;
- available route-owned density inputs;
- Coulomb expansion presence/term count compatibility if that information is
  available from inputs;
- explicit `electron_count`;
- closed-shell occupation policy:
  - positive integer electron count;
  - even electron count;
  - `nocc = electron_count ÷ 2`;
  - occupancy `2`;
- fixture role, initially requiring an explicit role such as `:route_smoke` or
  `:physics_endpoint`.

Return a compact object or NamedTuple with:

- `object_kind`;
- `status`;
- `blocker`;
- `missing_inputs`;
- `electron_count`;
- `occupation`;
- `fixture_role`;
- `summary`;
- `metadata`.

For blocked cases, keep precise blockers such as:

- `:missing_source_plan`;
- `:missing_final_basis`;
- `:missing_h1_payload`;
- `:missing_density_inputs`;
- `:missing_electron_count`;
- `:unsupported_open_shell_rhf_input`;
- `:missing_fixture_role`.

Do not:

- build a Fock matrix;
- run SCF;
- compute RHF energy;
- wire route driver/report fields;
- add new driver options;
- infer electron count from nuclear charge;
- add GTO, IDA/MWG, exports, artifacts, fixture promotion, or production route
  behavior;
- add scalar report-field clouds.

Tests:
Add one small focused contract test if needed. Prefer synthetic compact inputs
over the heavy source-box dry-run. The test should catch non-obvious contract
bugs:

- available closed-shell route-smoke input contract;
- blocked missing electron count;
- blocked odd electron count or unsupported open-shell input;
- blocked missing fixture role.

Do not run broad tests.

Validation:

`julia --project=. <focused new test if added>`
`julia --project=. -e 'using GaussletBases; println("load ok")'`
`git diff --check`

Report back:

- files changed;
- helper/object name and fields;
- blocked/available statuses and blockers;
- tests added and what bug they catch;
- validation commands/results;
- git status;
- deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

After writing `.agent_handoffs/response.105.md`, continue polling for
`blurb.106.md`, `STOP.md`, or `ATTENTION.md`.

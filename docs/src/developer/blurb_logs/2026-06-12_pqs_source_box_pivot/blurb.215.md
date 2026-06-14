Pass 215 - add or block the physical H2 private RHF input contract

Role:
You are `repo-doer@macmini` implementing one bounded solver-boundary pass for
GaussletBases. Follow `AGENTS.md`, `JuliaStyle.md`, `BlurbStyle.md`, and the
unattended baton rules in `.agent_handoffs/RUN.md`.

Loop/approval rule:
- Unattended baton mode is active.
- Do not ask the user for permission through an escalation prompt.
- If a needed command requires approval, write `.agent_handoffs/ATTENTION.md`
  with the exact command, why it is needed, and the blocking condition, then
  stop.

Current state:
- Head before this pass should be:
  `d06e8dcc Add H2 physical gausslet H1-J`
- The H2 physical q5 gausslet-only driver path now has:
  - source plan available;
  - final basis available;
  - H1 available;
  - density/H1-J seam available;
  - final dimension `463`;
  - H1 lowest energy about `-0.7946609179724647`;
  - H1-J/self-Coulomb scalar about `0.45696639804337114`;
  - density gauge `:pre_final_localized_positive_weight`;
  - raw pair-factor convention `:raw_numerator`;
  - support weight count `1215`.
- The current endpoint blocker is:
  `:missing_physical_gausslet_rhf_or_solver_contract`.

Physics target:
H2 at R = 4.0, gausslet-only, q = n_s = 5, no supplement. This remains the
active PQS physics target. H2 221 remains source-box diagnostic only. He 419 is
the already validated atom endpoint.

Architectural guardrail:
WL and PQS should share the same physical shell/support decomposition wherever
possible. For H2 R=4 q=5, the route-level physical inventory is:

```text
support counts  = (275, 578, 362)
retained counts = (251, 98, 114)
retained order  = (:atom_contact_core, :shared_shell_1, :shared_shell_2)
final dimension = 463
```

WL and PQS differ in retained transform and fast operator construction, not in
the physical support/shell vocabulary. Do not reintroduce the 221 diagnostic
route as a solver target.

Purpose:
Add the private physical-H2 RHF input contract boundary, or block it precisely.
This pass should decide whether the already materialized physical H2 H1 and
density/H1-J payloads contain enough reviewed data to request the existing
private RHF machinery. It should not run SCF unless the current code already
has a clean, small adapter and the contract work is trivial. Default expectation
is input contract only, no SCF.

Important context:
The existing private RHF helpers in `src/pqs_multilayer_complete_core_shell_rhf.jl`
were written for the old one-center complete-core/shell payloads. They currently
gate on object kinds like:

```text
:pqs_multilayer_shell_source_plan
:pqs_complete_core_shell_final_basis
:pqs_multilayer_complete_core_shell_h1_payload
:pqs_complete_core_shell_pre_final_density_interaction
```

The physical H2 payloads use physical labels instead:

```text
:pqs_diatomic_physical_gausslet_core_shell_source_plan
:pqs_physical_gausslet_final_basis
:available_pqs_physical_gausslet_h1_payload
:pqs_physical_gausslet_pre_final_density_interaction
```

Do not fake old object kinds. Either add a small physical-H2 input contract, or
refactor the existing private RHF input-contract validators to accept both
reviewed payload families by explicit predicates/accessors.

Task:
Implement the smallest private contract that consumes the existing physical H2
objects:

```text
_PQSDiatomicPhysicalGaussletCoreShellSourcePlan
_PQSDiatomicPhysicalGaussletFinalBasisPayload
_PQSDiatomicPhysicalGaussletH1Payload
_PQSDiatomicPhysicalGaussletH1JPayload
```

The contract should validate and summarize:

```text
route kind and fixture role
electron count = 2
closed-shell occupation count = 1
final dimension = 463
H1 matrix available/finite/symmetric
density interaction available
density gauge = :pre_final_localized_positive_weight
raw pair-factor convention = :raw_numerator
final-to-pre-final transform available
pre-final pair matrix available/finite/symmetric
private diagnostic only
```

Expected success transition if the input contract is available:

```text
route/private_rhf_input_contract_status =
  :available_pqs_physical_gausslet_rhf_input_contract
private_rhf/input_contract_available = true
private_rhf/electron_count = 2
private_rhf/occupation = closed shell, nocc = 1
physics/endpoint_ready = false
physics/endpoint_blocker = :missing_physical_gausslet_private_rhf_execution
```

Expected blocker examples:

```text
:missing_physical_gausslet_h1_payload
:missing_physical_gausslet_density_interaction
:missing_physical_gausslet_final_to_pre_final_transform
:missing_physical_gausslet_pre_final_pair_matrix
:unsupported_physical_gausslet_electron_count
:unsupported_physical_gausslet_fixture_role
:physical_gausslet_rhf_input_contract_unreviewed
```

If you find the existing private RHF helpers can be safely reused after adding
explicit physical-object predicates, stop at the input-contract boundary. Do not
promote RHF to a product and do not add robust SCF/DIIS tuning.

Do not:
- implement new SCF strategy or DIIS tuning;
- run HFDMRG, DMRG, or CR2;
- add GTO/MWG supplement;
- compare to supplemented WL/QW H2 references;
- build public API/export/artifact-ready HamV6;
- pass along final-basis self-overlap as working downstream data;
- mutate or reuse the H2 221 diagnostic route;
- add broad/default-runner tests.

Artifact/test expectations:
- Extend the H2 physical driver artifact with compact private-RHF input-contract
  fields only. Do not add a report-field cloud.
- Keep `private_rhf/materialized = false` unless an existing clean execution path
  is genuinely invoked in this pass.
- Update
  `test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl`
  to assert the compact contract fields and the new remaining blocker.
- Preserve the active He and H2 endpoint tests.

Line-count rule:
The active source/test/bin line-count rule still applies. For edits under
`src`, `test`, `bin`, and the CR2 generator script, final tracked diff must be
net-negative:

```sh
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Docs and blurb logs do not count. If source lines are added, pay for them by
deleting stale source/test pressure in the same pass. Do not delete the active
He or H2 driver endpoint tests.

Line-budget deletion candidates to inspect first:
- `test/nested/cartesian_ham_builder_diatomic_config_smoke_runtests.jl`
  is about 603 lines and may now be superseded by focused driver endpoints plus
  the low-order materialization extraction.
- `test/nested/cartesian_ham_builder_one_center_config_smoke_runtests.jl`
  is about 195 lines; do not delete if it is still the only live one-center
  driver smoke.
- `test/nested/pqs_explicit_core_spacing_parent_axis_probe_runtests.jl`
  and `test/nested/pqs_raw_product_box_plan_probe_runtests.jl` are old probe
  scaffolds; only shrink/delete if their active contracts are already covered by
  H2 endpoint or smaller module tests.

If none of these can be safely deleted or shrunk, report exact live callers and
the reason. Do not satisfy line budget by deleting active scientific endpoint
coverage.

Validation:
Run:

```sh
julia --project=. -e 'using GaussletBases; println("load ok")'
```

Run the H2 463 driver-artifact test with Julia-level timing:

```sh
julia --project=. -e 'using Test; t = @elapsed include("test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl"); println("elapsed_s=", t)'
```

This test was about 78 seconds after pass 214. If it grows substantially,
report the likely phase. Do not start a broader suite unless a failure requires
a narrowly related check.

Always run:

```sh
git diff --check
git diff --cached --check
```

after staging if you stage.

Response file:
Write `.agent_handoffs/response.215.md` and also copy it to:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.215.md
```

Report:
- whether a physical-H2 RHF input contract is available or blocked;
- statuses/blocker;
- electron count and occupation policy;
- final dimension;
- H1 and density-interaction readiness facts;
- whether any existing RHF helper was adapted or left untouched;
- artifact fields changed;
- source/test/bin scoped line budget added/deleted/net;
- validation commands and elapsed times;
- deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

Continue polling after writing the response. Manager will review, commit/push,
and publish the next blurb unless a real blocker appears.

-- repo-manager@macmini

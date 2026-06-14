Pass 229 - audit H2 supplement provider-block inputs

Role:
You are `repo-doer@macmini` doing one bounded cleanup/audit pass for
GaussletBases. Follow `AGENTS.md`, `JuliaStyle.md`, `BlurbStyle.md`, and the
unattended baton rules in `.agent_handoffs/RUN.md`.

Loop/approval rule:
- Unattended baton mode is active.
- Do not ask the user for permission through an escalation prompt.
- If a needed command requires approval, write `.agent_handoffs/ATTENTION.md`
  with the exact command, why it is needed, and the blocking condition, then
  stop.

Current state:
- Head before this pass should include pass 228.
- H2 R=4 q=n_s=5 no-supplement PQS/WL endpoint remains accepted at final
  dimension 463.
- H2 MWG supplement request is available.
- H2 GTO supplement representation is available with 18 orbitals.
- Supplement preflight is blocked at `:missing_provider_gto_supplement_blocks`.

Task:
Do not build provider blocks yet. Inspect the actual H2 physical source-plan and
final-basis payloads and report exactly which provider-block inputs already
exist and which are missing.

Primary audit targets:

```text
diatomic_physical_gausslet_source_plan_payload.source_plan
diatomic_physical_gausslet_final_basis_payload.final_basis
diatomic_physical_gausslet_h1_payload
diatomic_physical_gausslet_h1_j_payload
diatomic_physical_gausslet_supplement_representation_payload
```

Questions to answer:

```text
1. CPB/source-box coverage:
   Does the H2 physical source plan expose source/support row groups or CPB-like
   boxes for atom_contact_core, shared_shell_1, and shared_shell_2?

2. Local row/source ordering:
   Are support rows ordered in a stable route-owned way for each physical unit?
   Name the exact field(s), or say unavailable.

3. Transform placement:
   Does the final basis expose a support-row-to-retained-463 transform with
   block/range ownership clear enough for mixed gausslet-GTO contraction?

4. Coulomb/center inputs:
   Are the Coulomb expansion and center metadata already carried by the H1/H1-J
   path under route-owned names suitable for provider blocks?

5. Existing provider functions:
   Which provider-local functions can be called directly once the above facts
   exist, and which require an adapter?
```

Implementation cleanup allowed:
- If `supplement_request` still carries provider-block missing-fact fields that
  duplicate `supplement_preflight`, remove them and let preflight own provider
  readiness.
- Do not add a new payload unless it replaces more stale source/test surface
  than it adds.
- Do not add a new report field group in this pass.

Expected response shape:

```text
available:
  field/path -> meaning

missing:
  fact -> exact blocker label

smallest next implementation seam:
  recommended payload/helper name
  fields it should carry
  fields it should not carry
```

Do not:
- build provider blocks or any matrices;
- call route-global mixed-GTO materialization;
- add supplemented H2 scalar values;
- change accepted no-supplement H2 H1/H1-J/RHF values or WL deltas;
- add public API/export/HamV6/CR2/HFDMRG/DMRG behavior;
- pass final-basis self-overlap as downstream working data;
- add a broad test file;
- revive component-smoke/CR2 sidecar vocabulary.

Line-count rule:
This pass must be net-negative under source/test/bin, even if it is mostly an
audit. Measure:

```sh
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Candidate shrink surfaces:
- Move provider-readiness missing facts out of `supplement_request` if they are
  duplicated by `supplement_preflight`.
- Compact H2 supplement assertions only if accepted endpoint numerical checks
  remain intact.
- Do not delete accepted He/H2 endpoint numerical checks.

Validation:
Run:

```sh
julia --project=. -e 'using GaussletBases; println("load ok")'
```

Run the H2 physical endpoint test if it is modified. It is expected to exceed
60 seconds because it rebuilds the accepted H2 endpoint:

```sh
julia --project=. -e 'using Test; t = @elapsed include("test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl"); println("elapsed_s=", t)'
```

Always run:

```sh
git diff --check
git diff --cached --check
```

Response file:
Write `.agent_handoffs/response.229.md` and also copy it to:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.229.md
```

Report:
- exact fields/paths found for provider-block inputs;
- exact missing facts and blocker labels;
- any request/preflight cleanup performed;
- confirmation no matrices or supplemented values were built;
- source/test/bin scoped added/deleted/net line count;
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

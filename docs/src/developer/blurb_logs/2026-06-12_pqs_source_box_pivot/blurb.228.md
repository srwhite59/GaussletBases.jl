Pass 228 - tighten supplement representation authority and audit provider blocks

Role:
You are `repo-doer@macmini` doing one cleanup/audit pass for GaussletBases.
Follow `AGENTS.md`, `JuliaStyle.md`, `BlurbStyle.md`, and the unattended baton
rules in `.agent_handoffs/RUN.md`.

Loop/approval rule:
- Unattended baton mode is active.
- Do not ask the user for permission through an escalation prompt.
- If a needed command requires approval, write `.agent_handoffs/ATTENTION.md`
  with the exact command, why it is needed, and the blocking condition, then
  stop.

Current state:
- Head before this pass should include pass 227.
- H2 R=4 q=n_s=5 no-supplement PQS/WL endpoint remains accepted at final
  dimension 463.
- The H2 MWG supplement request is route-owned.
- The H2 GTO supplement representation is available and matrix-free:

```text
object kind = :cartesian_gaussian_shell_supplement_representation
center count = 2
orbital count = 18
provider blocks materialized = false
```

Physics/architecture guardrail:
WL and PQS share the same physical support/shell decomposition. For this H2
target:

```text
support counts  = (275, 578, 362)
retained counts = (251, 98, 114)
retained order  = (:atom_contact_core, :shared_shell_1, :shared_shell_2)
gausslet final dimension = 463
```

The supplement policy sits above the WL/PQS retained-transform distinction.
Do not create a PQS-only supplement theory.

Task:
Do not build provider blocks yet. First tighten the route/artifact boundary and
audit the exact provider-block seam.

Implementation cleanup:
- Make `supplement_representation` the authority for representation status.
- Remove duplicate representation-status/object-kind fields from the
  supplement request summary/artifact if they are now only transitional
  duplication.
- Keep `supplement_request` as request metadata:

```text
status
blocker
fixture_label
supplement_policy
basis_name
lmax
atom_symbols
nuclear_charges
bond_axis
bond_length
required_provider_blocks
missing_fact_labels
matrices_materialized = false
```

- Keep `supplement_representation` as the representation metadata group:

```text
status
blocker
object_kind
basis_name
lmax
atom_symbols
center_count
orbital_count
matrices_materialized = false
provider_blocks_materialized = false
```

Audit output:
In the response, report the exact live provider-block seam for the next coding
pass:

```text
current available inputs:
  parent / axis bundles
  H2 physical source plan
  H2 physical final basis
  supplement representation
  Coulomb expansion / center metadata, if available

provider-local functions to reuse:
  cpb_gto_supplement_local_operator_bundle
  cpb_mixed_gto_overlap_block
  cpb_mixed_gto_kinetic_operator_block
  cpb_mixed_gto_position_operator_block
  cpb_mixed_gto_x2_operator_block
  cpb_mixed_gto_nuclear_by_center_block
  cpb_gto_overlap_operator_block
  cpb_gto_nuclear_by_center_block

missing route facts before blocks can materialize:
  CPB/source-box coverage for the H2 physical support units
  local row/source ordering for each support unit
  coefficient transform from support rows to retained 463 basis
  placement/accumulation rule for mixed gausslet-GTO blocks
  GTO/GTO self-block rule
  raw moment matrices needed by MWG residualization
```

If the audit proves one of these facts is already available in a named payload,
name the exact field and owner. If a fact is not available, keep the blocker
precise. Do not add a broad new field cloud to report this.

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
This pass must be net-negative under source/test/bin. Measure:

```sh
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Candidate shrink surfaces:
- Delete the duplicate request representation fields in source/reporting/tests.
- Further compact H2 request/representation/preflight assertions only if the
  accepted endpoint numerical checks remain intact.
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
Write `.agent_handoffs/response.228.md` and also copy it to:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.228.md
```

Report:
- exactly what duplicate request/report fields were removed or kept;
- current supplement request/representation/preflight statuses;
- provider-block seam audit with available facts and missing facts;
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

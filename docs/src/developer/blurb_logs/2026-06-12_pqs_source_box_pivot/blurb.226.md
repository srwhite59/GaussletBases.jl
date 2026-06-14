Pass 226 - add H2 PQS supplement request payload

Role:
You are `repo-doer@macmini` implementing one bounded request/metadata seam for
GaussletBases. Follow `AGENTS.md`, `JuliaStyle.md`, `BlurbStyle.md`, and the
unattended baton rules in `.agent_handoffs/RUN.md`.

Loop/approval rule:
- Unattended baton mode is active.
- Do not ask the user for permission through an escalation prompt.
- If a needed command requires approval, write `.agent_handoffs/ATTENTION.md`
  with the exact command, why it is needed, and the blocking condition, then
  stop.

Current state:
- Head before this pass should be the commit that publishes pass 225.
- H2 R=4 q=n_s=5 gausslet-only PQS/WL physical endpoint is accepted and
  comparison-ready at final dimension 463.
- Pass 224 added a metadata-only supplement preflight boundary.
- Pass 225 audited existing WL/MWG/GTO machinery and recommended a route-owned
  supplement request seam before provider-block materialization.

Physics target:
H2 R=4, q=n_s=5, common physical support plan:

```text
support counts  = (275, 578, 362)
retained counts = (251, 98, 114)
retained order  = (:atom_contact_core, :shared_shell_1, :shared_shell_2)
gausslet final dimension = 463
```

Task:
Add a compact, private, matrix-free H2 physical supplement request payload.

Suggested internal name:

```text
_PQSDiatomicPhysicalGaussletSupplementRequestPayload
```

Suggested fields:

```text
status
blocker
route_family
route_kind
fixture_label
supplement_policy
atom_symbols
nuclear_charges
atom_locations
bond_axis
bond_length
basis_name
lmax
uncontracted
residual_keep_policy
residual_drop_tolerance
representation_status
representation_object_kind
required_provider_blocks
available_fact_labels
missing_fact_labels
summary
metadata
```

Default request facts for this H2 target should be explicit and conservative:

```text
supplement_policy = :mwg_residual_gto
basis_name = "H/cc-pVTZ"
lmax = 1
fixture_label = :h2_r4_physical_gausslet_q5
residual policy = route-private MWG residual GTO preflight only
```

Decision rule:
- For `supplement_policy = :none`, request status should be `:not_requested`.
- For `supplement_policy = :mwg_residual_gto`, create the route-owned request
  payload if the H2 physical target inventory is available.
- Do not build provider blocks or any matrices.
- If creating a real `CartesianGaussianShellSupplementRepresentation3D` is
  straightforward using existing route-neutral code and does not build matrices,
  it is allowed, but not required.
- If the representation is not built, set a precise blocker such as
  `:missing_gto_supplement_representation` and make the preflight blocker point
  there before `:missing_provider_gto_supplement_blocks`.
- If the representation is built, keep the next blocker as
  `:missing_provider_gto_supplement_blocks`.

Artifact/reporting:
Write compact fields only. A small `supplement_request` group is acceptable:

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
representation_status
required_provider_blocks
missing_fact_labels
matrices_materialized = false
```

Do not copy full center tables or large objects into the artifact.

Do not:
- build GTO/GTO, mixed gausslet/GTO, or MWG residual matrices;
- add supplemented H2 scalar values;
- change accepted no-supplement H2 H1/H1-J/RHF values or WL deltas;
- add public API/export/HamV6/CR2/HFDMRG/DMRG behavior;
- pass final-basis self-overlap as downstream working data;
- add a new broad test file;
- revive component-smoke/CR2 sidecar vocabulary.

Line-count rule:
This pass must be net-negative under source/test/bin. Measure:

```sh
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Candidate shrink surfaces:
- In `test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl`,
  replace repeated `supplement_preflight/*` field assertions with one compact
  request/preflight fingerprint plus endpoint-critical fields.
- In `test/nested/pqs_source_box_route_driver_report_runtests.jl`, continue
  deleting stale route-report metadata assertions that only preserve
  transitional field names.
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
Write `.agent_handoffs/response.226.md` and also copy it to:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.226.md
```

Report:
- request payload status/blocker behavior for `:none` and
  `:mwg_residual_gto`;
- whether a representation object was created or explicitly left blocked;
- artifact fields written;
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

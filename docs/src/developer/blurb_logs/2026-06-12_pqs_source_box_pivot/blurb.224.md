Pass 224 - add H2 physical supplement preflight boundary

Role:
You are `repo-doer@macmini` implementing one bounded metadata/readiness seam
for GaussletBases. Follow `AGENTS.md`, `JuliaStyle.md`, `BlurbStyle.md`, and
the unattended baton rules in `.agent_handoffs/RUN.md`.

Loop/approval rule:
- Unattended baton mode is active.
- Do not ask the user for permission through an escalation prompt.
- If a needed command requires approval, write `.agent_handoffs/ATTENTION.md`
  with the exact command, why it is needed, and the blocking condition, then
  stop.

Current state:
- Head before this pass should be:
  `23990c8c Retire component smoke sidecar surface`
- H2 R=4 q=n_s=5 gausslet-only PQS/WL physical endpoint is accepted and
  comparison-ready.
- The old component-smoke/CR2 sidecar route-shadow surface has been deleted.

Physics target:
H2 R=4, q=n_s=5, physical support plan:

```text
support counts  = (275, 578, 362)
retained counts = (251, 98, 114)
retained order  = (:atom_contact_core, :shared_shell_1, :shared_shell_2)
gausslet-only final dimension = 463
```

Architectural guardrail:
WL and PQS share the physical support/intermediate gausslet plan. They differ in
retained transform and fast operator construction. GTO/MWG supplement policy
should sit above that distinction:

```text
common physical support/intermediate gausslet plan
+ retained_transform_kind = :pqs
+ supplement_policy = :mwg_residual_gto
```

Purpose:
Add the first private H2 physical supplement preflight/payload boundary. This
pass is metadata/readiness only. It should make the next required facts explicit
without building GTO matrices or changing the accepted no-supplement endpoint.

Suggested implementation surface:
- `src/pqs_source_box_diatomic_complete_core_shell.jl` for the private physical
  supplement preflight payload if that is where adjacent H2 physical target,
  source-plan, final-basis, H1, H1-J, and RHF payloads live.
- `src/pqs_source_box_route_driver_reporting.jl` only if needed to write compact
  artifact fields.
- `bin/cartesian_ham_builder.jl` only if a missing visible-driver input default
  is needed.
- Existing H2 physical endpoint test only if you can add a compact override case
  without making it broader.

Preflight behavior:
For the accepted H2 physical route:

```text
supplement_policy = :none
  -> supplement_preflight/status = :not_requested
  -> accepted no-supplement endpoint values remain unchanged

supplement_policy = :mwg_residual_gto
  -> supplement_preflight/status = :blocked_pqs_physical_gausslet_mwg_residual_gto_preflight
  -> supplement_preflight/blocker = first precise missing fact
```

Use compact status/blocker vocabulary. Expected missing facts are likely:

```text
:missing_gto_supplement_request
:missing_provider_gto_supplement_blocks
:missing_mixed_gausslet_gto_blocks
:missing_gto_gto_blocks
:missing_combined_raw_moment_matrices
:missing_residual_mwg_representation
:missing_combined_density_density_readiness
```

Use the exact names that fit the existing code, but keep them specific.

The preflight should record compact facts, not matrices:

```text
route kind
system/fixture label
support counts
retained counts
retained order
retained_transform_kind = :pqs
gausslet_final_dimension = 463
supplement_policy
required fact labels
available fact labels
status/blocker
```

Do not:
- build GTO/GTO, mixed gausslet/GTO, or MWG residual matrices;
- write supplemented H2 comparison values;
- change the accepted no-supplement H2 endpoint values or deltas;
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
- `test/nested/pqs_source_box_route_driver_report_runtests.jl` still carries
  many route-report metadata assertions. Shrink only stale report-field
  assertions that are superseded by the accepted H2 physical endpoint and the
  current compact artifact checks.
- `test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl`
  remains very large. If you need line budget, remove stale route-shadow or
  no-go assertion blocks only when a live mathematical check remains.

Do not delete accepted He/H2 endpoint tests to satisfy line budget.

Validation:
Run:

```sh
julia --project=. -e 'using GaussletBases; println("load ok")'
```

Run the smallest focused test that exercises the new preflight. If you add the
check to the H2 physical endpoint test and it still takes over 60 seconds,
report that explicitly with elapsed time:

```sh
julia --project=. -e 'using Test; t = @elapsed include("test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl"); println("elapsed_s=", t)'
```

Always run:

```sh
git diff --check
git diff --cached --check
```

Response file:
Write `.agent_handoffs/response.224.md` and also copy it to:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.224.md
```

Report:
- preflight object/status/blocker behavior for `:none` and `:mwg_residual_gto`;
- exact fields written, if any;
- confirmation no supplemented matrices/values were built;
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

Pass 210 - add checked H2 physical source-plan candidate adapter

Role:
You are `repo-doer@macmini` implementing one bounded H2 physical-route
readiness seam for GaussletBases. Follow `AGENTS.md`, `JuliaStyle.md`,
`BlurbStyle.md`, and the unattended baton rules in `.agent_handoffs/RUN.md`.

Loop/approval rule:
- Unattended baton mode is active.
- Do not ask the user for permission through an escalation prompt.
- If a needed command requires approval, write `.agent_handoffs/ATTENTION.md`
  with the exact command, why it is needed, and the blocking condition, then
  stop.

Current state:
- Head before this pass should be:
  `1d8a0d55 Audit H2 physical source plan inputs`
- The H2 463 physical source-plan seam is blocked on:
  `:missing_atom_contact_core_support_rows`
- Pass 209 found that candidate rows/coefs exist in old/source-backed and
  atom-growth materializer objects, but neither is yet route-owned PQS source
  authority.

Physics target:
H2 at R = 4.0, gausslet-only, q = n_s = 5, no supplement, physical target:

```text
support counts:  (275, 578, 362)
retained counts: (251, 98, 114)
expected final dimension: 463
support/retained order:
  (:atom_contact_core, :shared_shell_1, :shared_shell_2)
```

Purpose:
Add a checked candidate adapter/readiness summary that can locate and verify the
physical H2 rows/coefs, without yet declaring the physical PQS source plan
materialized unless the route-authority question is genuinely resolved.

Task:
Add one compact private candidate payload/helper near the physical source-plan
payload, for example:

```text
_PQSDiatomicPhysicalGaussletSourcePlanCandidatePayload
_pqs_source_box_route_driver_diatomic_physical_gausslet_source_plan_candidate_payload
```

The candidate should attempt to locate one of the known construction sources:

1. Preferred as an oracle/adapter candidate:
   `bond_aligned_diatomic_nested_fixed_source(...; nside = 5)` or the
   equivalent source-backed object for the H2 R=4 basis.

2. Alternative candidate:
   the route-configured atom-growth complete-rectangular materializer output,
   if it is already available in the driver path and can be checked without
   broad rewiring.

Do not silently choose a different route. The candidate summary must say which
source was used:

```text
:source_backed_fixed_source_oracle
:atom_growth_low_order_materializer
:not_available
```

Required checks:
- support order exactly `(:atom_contact_core, :shared_shell_1, :shared_shell_2)`;
- retained order exactly `(:atom_contact_core, :shared_shell_1, :shared_shell_2)`;
- support counts exactly `(275, 578, 362)`;
- retained counts exactly `(251, 98, 114)`;
- final dimension exactly `463`;
- no residual/GTO supplement;
- no H2 221 diagnostic source-plan object kind;
- source-backed/low-order authority flag is explicit.

Expected status behavior:
- If the candidate cannot be built, keep the existing source-plan payload
  blocked, but use a sharper blocker such as:
  `:missing_physical_gausslet_source_plan_materializer`.
- If the candidate can be built and all counts/order checks pass, report:

```text
candidate_status = :available_physical_gausslet_source_plan_candidate
candidate_source = ...
candidate_counts_match = true
```

but keep the actual source plan blocked unless the candidate is route-owned
authority. A good blocker in that case is:

```text
:source_plan_candidate_not_route_authority
```

Do not call this a materialized source plan merely because an oracle candidate
exists.

Artifact/report:
- Extend the existing H2 physical target artifact with compact candidate fields
  only, for example:

```text
target/source_plan_candidate_status
target/source_plan_candidate_source
target/source_plan_candidate_counts_match
target/source_plan_authority_status
```

- The artifact should still show `physics/endpoint_ready = false`.
- If candidate exists but is not authority, `physics/endpoint_blocker` should be
  `:source_plan_candidate_not_route_authority`.

Do not:
- implement final basis, H1, H1-J, density interaction, RHF, HFDMRG, CR2, or
  exports;
- compare to supplemented WL/QW H2 references;
- reuse the H2 221 source-box diagnostic plan;
- promote low-order/WL materialization to route authority by naming alone;
- add a broad new test file;
- add many scalar report aliases.

Line-count rule:
The active source/test/bin line-count rule still applies. For edits under
`src`, `test`, `bin`, and the CR2 generator script, final tracked diff must be
net-negative:

```sh
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Docs and blurb logs do not count. If the candidate adapter adds source lines,
pay for them by deleting stale source/test pressure in the same pass.

Suggested shrink candidates:
- The old route-configured atom-growth materializer/report path still has many
  flat diagnostic aliases in `src/pqs_source_box_low_order_materialization.jl`;
  if any are source-unreferenced or only preserve report vocabulary, shrink
  there carefully.
- Old non-default route-shadow tests may contain repeated receipt/diagnostic
  assertions. Delete only stale internal-vocabulary checks, not geometry,
  dimension, finite-output, symmetry, or physics checks.

Validation:
Run:

```sh
julia --project=. -e 'using GaussletBases; println("load ok")'
```

Run the H2 463 artifact-readiness test with Julia-level timing:

```sh
julia --project=. -e 'using Test; t = @elapsed include("test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl"); println("elapsed_s=", t)'
```

If constructing the candidate makes that test substantially slower than the
recent 60-65 second range, report the timing and the phase likely responsible.

Always run:

```sh
git diff --check
git diff --cached --check
```

after staging if you stage.

Response file:
Write `.agent_handoffs/response.210.md` and also copy it to:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.210.md
```

Report:
- candidate status/source/counts;
- whether source plan remains blocked or became available;
- exact blocker labels;
- artifact fields added;
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

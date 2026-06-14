Pass 211 - materialize private H2 physical source-plan object from checked candidate

Role:
You are `repo-doer@macmini` implementing one bounded H2 physical-route source
plan seam for GaussletBases. Follow `AGENTS.md`, `JuliaStyle.md`,
`BlurbStyle.md`, and the unattended baton rules in `.agent_handoffs/RUN.md`.

Loop/approval rule:
- Unattended baton mode is active.
- Do not ask the user for permission through an escalation prompt.
- If a needed command requires approval, write `.agent_handoffs/ATTENTION.md`
  with the exact command, why it is needed, and the blocking condition, then
  stop.

Current state:
- Head before this pass should be:
  `804fa551 Add H2 physical source plan candidate`
- The H2 463 candidate is available:
  - candidate source `:source_backed_fixed_source_oracle`
  - candidate counts match `true`
  - support counts `(275, 578, 362)`
  - retained counts `(251, 98, 114)`
  - final dimension `463`
- The actual source plan remains blocked on:
  `:source_plan_candidate_not_route_authority`

Physics target:
H2 at R = 4.0, gausslet-only, q = n_s = 5, no supplement. The source plan must
represent:

```text
atom_contact_core support/coefs + shared_shell_1 support/coefs + shared_shell_2 support/coefs
```

with core-then-shared order and expected dimension `463`.

Purpose:
Turn the verified source-backed candidate into an explicit private route-owned
physical source-plan object, if and only if the adapter can carry the actual
support rows/states/coefficient matrices and preserve the source-backed
provenance honestly.

Task:
Add a private source-plan object for the physical H2 route. Use explicit naming,
for example:

```text
_PQSDiatomicPhysicalGaussletCoreShellSourcePlan
:pqs_diatomic_physical_gausslet_core_shell_source_plan
```

This object should carry the actual source-plan data by reference/value:

- parent/bundles or metrics needed by later final-basis construction;
- atom-contact-core support indices and states;
- shared-shell support indices and states;
- shell coefficient matrix for the shared layers;
- core coefficient matrix or enough core sequence data to reconstruct the
  retained core columns;
- support order `(:atom_contact_core, :shared_shell_1, :shared_shell_2)`;
- retained order `(:atom_contact_core, :shared_shell_1, :shared_shell_2)`;
- retained ranges / final dimension;
- convention labels:
  - `source_plan_family = :physical_gausslet_core_shell_source_plan`
  - `source_backed_adapter = true`
  - `source_backed_candidate_source = :source_backed_fixed_source_oracle`
  - `route_owned_authority = true`
  - `supplement_policy = :none`
  - `h2_221_diagnostic_source_plan_reused = false`
  - final basis/H1/H1-J/RHF/export/public flags false.

The route-owned authority claim is allowed only for this private adapter if:

- the candidate is available and count-checked;
- support and retained order match exactly;
- support/retained counts match exactly;
- final dimension is exactly `463`;
- no supplement is present;
- the constructed object does not claim the old 221 diagnostic object kind.

Expected status transition:

```text
source_plan_status:
  :blocked_pqs_diatomic_physical_gausslet_core_shell_source_plan
    -> :available_pqs_diatomic_physical_gausslet_core_shell_source_plan

source_plan_authority_status:
  :candidate_not_route_authority
    -> :private_source_backed_adapter_authority
```

After this pass, the H2 physical artifact should still have:

```text
physics/endpoint_ready = false
```

but the blocker should advance to the next seam, for example:

```text
:missing_physical_gausslet_final_basis_builder
```

Do not:
- implement final basis;
- implement H1, H1-J, density interaction, RHF, HFDMRG, CR2, exports, or public
  API;
- compare to supplemented WL/QW H2 references;
- reuse or mutate the H2 221 source-box diagnostic plan;
- add broad tests or default-runner includes.

Line-count rule:
The active source/test/bin line-count rule still applies. For edits under
`src`, `test`, `bin`, and the CR2 generator script, final tracked diff must be
net-negative:

```sh
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Docs and blurb logs do not count. If source lines are added, pay for them by
deleting stale source/test pressure in the same pass. Prefer deleting
transitional metadata/report-vocabulary checks over touching endpoint or
physics checks.

Validation:
Run:

```sh
julia --project=. -e 'using GaussletBases; println("load ok")'
```

Run the H2 463 artifact-readiness test with Julia-level timing:

```sh
julia --project=. -e 'using Test; t = @elapsed include("test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl"); println("elapsed_s=", t)'
```

If the test runtime remains around 75 seconds, report that. If it grows
substantially, report the likely phase.

Always run:

```sh
git diff --check
git diff --cached --check
```

after staging if you stage.

Response file:
Write `.agent_handoffs/response.211.md` and also copy it to:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.211.md
```

Report:
- source-plan object/status/blocker;
- whether authority status advanced to
  `:private_source_backed_adapter_authority`;
- support/retained counts and final dimension;
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

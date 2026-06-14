Pass 219 - add H2 WL gausslet-only reference-candidate readiness

Role:
You are `repo-doer@macmini` implementing one bounded no-matrix comparison
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
  `bff7aa43 Publish PQS blurb 218`
- Pass 218 found the honest next blocker:
  `:missing_route_configured_wl_h2_gausslet_only_463_reference_candidate`.
- The current H2 physical PQS endpoint is:
  - H2 R=4, q=n_s=5, supplement none;
  - final dimension 463;
  - H1/H1-J/private RHF diagnostic materialized;
  - still not comparison-ready.

Physics target:
H2 at R = 4.0, gausslet-only, q = n_s = 5, no supplement. H2 221 remains
source-box diagnostic only. He 419 is already validated.

Architectural guardrail:
WL and PQS share the physical shell/support decomposition. They differ in the
retained transform and fast operator construction path. For the H2 R=4 q5
gausslet-only comparison candidate, require:

```text
parent axis counts = (9, 9, 15)
support counts     = (275, 578, 362)
retained counts    = (251, 98, 114)
retained order     = (:atom_contact_core, :shared_shell_1, :shared_shell_2)
final dimension    = 463
supplement policy  = :none
```

Purpose:
Add a small route-owned no-matrix readiness object for the matching WL/old-QW
gausslet-only H2 reference candidate. This is not a scalar comparison and not a
WL Hamiltonian materialization. It only records whether the current driver/report
has enough reviewed facts to name the matching WL gausslet-only reference route.

Task:
Add a compact candidate/readiness seam, preferably near the existing physical H2
target/reporting code:

```text
_pqs_source_box_route_driver_h2_wl_gausslet_only_reference_candidate(...)
```

or a clearer local name if the surrounding code suggests one.

The candidate should be available only when all reviewed conditions hold:
- current route is `:bond_aligned_diatomic_physical_gausslet_core_shell_pqs`;
- geometry is H2 R=4 with centers equivalent to `(0,0,-2)` and `(0,0,2)`;
- parent axis counts are `(9, 9, 15)`;
- target support counts are `(275, 578, 362)`;
- target retained counts are `(251, 98, 114)`;
- retained order is `(:atom_contact_core, :shared_shell_1, :shared_shell_2)`;
- supplement policy is `:none`;
- final dimension is `463`;
- the reference label is explicit, e.g.
  `WL/QW H2 R=4 gausslet-only 463`;
- retained transform kind is explicitly WL/old-QW, not PQS.

If any condition fails, return a blocked candidate with a precise blocker. Use
one blocker if practical, for example:

```text
:wl_h2_gausslet_only_reference_candidate_mismatch
```

and include compact mismatch facts in the summary.

Artifact/report expectation:
Expose only compact comparison-readiness fields, likely under the existing
`comparison` group or the physical-target summary:

```text
comparison/wl_reference_candidate_status
comparison/wl_reference_candidate_blocker
comparison/wl_reference_final_dimension
comparison/wl_reference_retained_transform_kind
comparison/wl_reference_supplement_policy
comparison/wl_reference_label
```

Do not add large matrices, self-overlap matrices, source arrays, or full
metadata objects to the artifact. If the existing save path makes a different
compact field grouping cleaner, use it, but avoid a broad field cloud.

Expected endpoint transition:

If the no-matrix candidate is available, keep:

```text
comparison/ready = false
physics/endpoint_ready = false
```

and move the endpoint blocker to a sharper next blocker such as:

```text
:missing_wl_h2_gausslet_only_reference_values
```

Do not mark comparison ready until actual WL gausslet-only H1/H1-J/RHF or other
reviewed scalar values exist.

Do not:
- build WL H1/H1-J/RHF matrices in this pass;
- add GTO/MWG supplement;
- compare to supplemented WL/QW H2 values;
- alter the H2 221 diagnostic route;
- add public API, export, HamV6, CR2, HFDMRG, or DMRG behavior;
- pass final-basis self-overlap as downstream working data;
- hide driver behavior behind an opaque `run_driver(config)`;
- add a broad/default-runner test.

Test expectation:
Update only the existing H2 physical endpoint test to assert the compact
candidate fields and the sharper blocker. Do not add a new test file.

Line-count rule:
The active source/test/bin line-count rule applies. Final tracked diff under
`src`, `test`, `bin`, and the CR2 generator script must be net-negative:

```sh
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Candidate stale-test deletion surfaces to audit first:
- `test/nested/cartesian_ham_builder_one_center_config_smoke_runtests.jl`
  remains listed as a one-center WL driver smoke. Delete/shrink only if He
  driver endpoints fully supersede it.
- `test/nested/cartesian_pair_block_route_state_global_safe_terms_runtests.jl`
  may be old route-global staging pressure. Delete/shrink only if live source
  callers and active endpoint tests no longer need it.
- Do not delete active He/H2 driver endpoint tests or compact CPB/local provider
  module-contract tests.

If no safe deletion exists, stop with exact live callers instead of adding a
net-positive source/test diff.

Validation:
Run:

```sh
julia --project=. -e 'using GaussletBases; println("load ok")'
```

Run the focused H2 endpoint with Julia-level timing:

```sh
julia --project=. -e 'using Test; t = @elapsed include("test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl"); println("elapsed_s=", t)'
```

The focused H2 endpoint was about 83 seconds after pass 217. If it grows
substantially, report why.

Always run:

```sh
git diff --check
git diff --cached --check
```

after staging if you stage.

Response file:
Write `.agent_handoffs/response.219.md` and also copy it to:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.219.md
```

Report:
- candidate status/blocker;
- exact candidate fields exposed;
- whether the old supplemented WL/QW scalar references remain blocked;
- endpoint blocker before and after;
- whether any actual WL matrix/scalar value was materialized, which should be
  no in this pass;
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

Pass 208 - first H2 physical gausslet source-plan seam

Role:
You are `repo-doer@macmini` implementing one bounded H2 physical-route seam for
GaussletBases. Follow `AGENTS.md`, `JuliaStyle.md`, `BlurbStyle.md`, and the
unattended baton rules in `.agent_handoffs/RUN.md`.

Loop/approval rule:
- Unattended baton mode is active again.
- Do not ask the user for permission through an escalation prompt.
- If a needed command requires approval, write `.agent_handoffs/ATTENTION.md`
  with the exact command, why it is needed, and the blocking condition, then
  stop.

Current state:
- Head before this pass should be:
  `abb2c443 Clarify driver endpoint taxonomy`
- The driver endpoint taxonomy is now explicit:
  - He 419: physical atom endpoint.
  - H2 221: source-box diagnostic only, not a physics endpoint.
  - H2 463: physical gausslet-only H2 target inventory, currently blocked on
    `:missing_physical_gausslet_source_plan`.

Physics target:
H2 at R = 4.0, gausslet-only, q = n_s = 5, no residual/GTO supplement, physical
route shape:

```text
(:atom_contact_core, :shared_shell_1, :shared_shell_2)
```

The reviewed target inventory is:

```text
support units:   (:atom_contact_core, :shared_shell_1, :shared_shell_2)
support counts:  (275, 578, 362)
retained counts: (251, 98, 114)
expected final dimension: 463
retained atom-core interiors: true
supplement_policy: :none
```

Purpose:
Start replacing the H2 463 target-inventory blocker with a real, route-owned
source-plan seam. This pass should either materialize an honest physical
gausslet source-plan object or report the exact missing construction fact. Do
not fake success by adapting the old 221 source-box diagnostic plan.

Relevant current surfaces:
- Physical target skeleton:
  `src/pqs_source_box_route_driver_skeletons.jl`
  - `_pqs_source_box_route_driver_physical_gausslet_core_shell_skeleton`
- Physical target payload:
  `src/pqs_source_box_diatomic_complete_core_shell.jl`
  - `_PQSDiatomicPhysicalGaussletCoreShellTargetPayload`
  - `_pqs_source_box_route_driver_diatomic_physical_gausslet_target_payload`
- Existing diagnostic source-plan path, for reference only:
  `src/pqs_source_box_diatomic_complete_core_shell.jl`
  - `_PQSDiatomicCompleteCoreShellSourcePlan`
  - `_pqs_source_box_route_driver_diatomic_complete_core_shell_source_plan`
  - `_pqs_source_box_route_driver_diatomic_complete_core_shell_source_plan_payload`
- Driver assembly currently calls both:
  `src/pqs_source_box_route_driver_helpers.jl`
  - `cartesian_assembly(...)`
- Physical H2 input/test:
  `test/driver_inputs/h2_pqs_q5_physical_gausslet_r4.jl`
  `test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl`
- Endpoint manifest:
  `docs/src/developer/cartesian_driver_endpoint_manifest.md`

Task:
Add the narrowest private source-plan seam for the H2 463 physical gausslet
target.

Preferred implementation shape:
1. Add a compact private physical source-plan payload/object near the existing
   diatomic complete-core/shell code. Good names would be explicit, for example:

```text
_PQSDiatomicPhysicalGaussletCoreShellSourcePlanPayload
_pqs_source_box_route_driver_diatomic_physical_gausslet_source_plan_payload
```

2. It may reuse shared structs/helpers only where the semantics really match.
   It must not claim the 221 diagnostic object kind unless it is truly the same
   source-plan contract. If a new object kind is needed, use one such as:

```text
:pqs_diatomic_physical_gausslet_core_shell_source_plan
```

3. The available path, if it can be built honestly, should carry:
   - target payload status;
   - parent axis counts;
   - support order `(:atom_contact_core, :shared_shell_1, :shared_shell_2)`;
   - retained order `(:atom_contact_core, :shared_shell_1, :shared_shell_2)`;
   - support counts `(275, 578, 362)`;
   - retained counts `(251, 98, 114)`;
   - expected/precleanup final dimension `463`;
   - retained atom-core interiors `true`;
   - source-plan role `:atom_contact_core_plus_pqs_shared_shells`;
   - supplement policy `:none`;
   - compact convention labels showing this is physical H2 gausslet-only,
     private route-owned, not supplemented, not RHF, not export/public.

4. If the actual support arrays, support states, coefficient matrices, or
   shared-shell contractions cannot yet be constructed from existing route-owned
   objects, do not synthesize placeholders. Return a blocked payload with the
   sharpest blocker, for example:

```text
:missing_atom_contact_core_support_rows
:missing_shared_shell_1_coefficients
:missing_shared_shell_2_coefficients
:missing_physical_gausslet_source_plan_construction
```

5. Wire the payload into assembly/report/artifact only as a compact status and
   summary. The H2 physical target artifact should move from the broad blocker
   `:missing_physical_gausslet_source_plan` to either:
   - available physical source-plan status, with final basis still not
     requested/materialized; or
   - a sharper physical source-plan blocker from this pass.

6. Update the H2 463 target test to assert the new status/blocker. Keep it an
   artifact-readiness test. It should not require final basis, H1, H1-J, RHF, or
   comparison readiness.

Do not:
- mutate the H2 221 source-box diagnostic route;
- use the H2 221 diagnostic plan as the H2 463 physical source plan;
- compare to supplemented WL/QW H2 references;
- implement final basis, H1, H1-J, density interaction, RHF, HFDMRG, CR2, or
  exports/artifacts beyond the existing driver artifact;
- add a broad new test file;
- add public API;
- add a flat field cloud of many report aliases.

Line-count rule:
The active source/test/bin line-count rule still applies. For edits under
`src`, `test`, `bin`, and the CR2 generator script, final tracked diff must be
net-negative:

```sh
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Docs and blurb logs do not count. If this implementation needs added source
lines, pay for them in the same pass by deleting or shrinking stale source/test
surfaces. It is acceptable to delete old scaffolding outside the touched H2
test when it is clearly obsolete and only preserves transitional vocabulary.
Do not delete live physics endpoint checks or scientific references to game the
budget.

Suggested shrink candidates if you need budget:
- Remove the non-contract `println` from
  `test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl`.
- Look for stale route-shadow metadata assertions in old non-default tests such
  as `test/diatomic/runtests.jl`; delete only clearly redundant
  report-text/internal-vocabulary checks, not geometry or physics contract
  checks.
- If a superseded helper is now source-unreferenced and only test-referenced,
  delete the helper and its stale test pressure together.

Validation:
Run:

```sh
julia --project=. -e 'using GaussletBases; println("load ok")'
```

Run the H2 463 target inventory/readiness test. It has recently run around 60
seconds because of driver compilation, so it is acceptable to run it with a
Julia-level elapsed timer and report the elapsed time:

```sh
julia --project=. -e 'using Test; t = @elapsed include("test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl"); println("elapsed_s=", t)'
```

Run a smaller focused check if you add a module-level source-plan helper and
can validate it without the full driver. Do not add a new broad test.

Always run:

```sh
git diff --check
git diff --cached --check
```

after staging if you stage.

Response file:
Write `.agent_handoffs/response.208.md` and also copy it to:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.208.md
```

Report:
- whether the physical H2 source-plan seam is available or blocked;
- exact status/blocker labels;
- support/retained counts and dimension observed;
- files edited/deleted;
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

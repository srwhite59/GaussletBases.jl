Pass 222 - retire H2 221 diagnostic scaffold and audit supplement seam

Role:
You are `repo-doer@macmini` doing one bounded cleanup plus design-audit pass
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
  `c486288d Wire H2 WL comparison values`
- The H2 R=4 q=n_s=5 gausslet-only physical endpoint is now comparison-ready:
  PQS and reviewed WL no-supplement values agree at about `1e-12` or better.
- The old H2 221 route remains a source-box diagnostic, not a physics endpoint.

Architectural guardrail:
Treat WL and PQS as sharing the same physical shell/support decomposition where
possible. For H2 R=4 q=5, the physical gausslet-only support plan is:

```text
support counts  = (275, 578, 362)
retained counts = (251, 98, 114)
retained order  = (:atom_contact_core, :shared_shell_1, :shared_shell_2)
final dimension = 463
```

WL and PQS should differ in retained transform and fast operator construction
path, not in the route-level physical shell/support vocabulary.

Purpose:
1. Remove stale slow test pressure from the superseded H2 221 source-box
   diagnostic route if caller audit confirms it is only development scaffold.
2. Audit the next GTO/MWG supplement seam for the H2 physical route without
   implementing supplement support yet.

Task A - retire H2 221 diagnostic scaffold:
1. Audit references to:

```text
test/driver_inputs/h2_pqs_q5_source_box_diagnostic_r4.jl
test/nested/cartesian_ham_builder_h2_pqs_q5_source_box_diagnostic_r4_runtests.jl
```

2. If no live source, default-runner, or current physics workflow depends on
   them, delete both files.
3. If deletion is blocked by a live caller, do not add an adapter. Report the
   exact caller and shrink the test/input as much as possible instead.
4. Do not delete the 221 route implementation itself in this pass. This is a
   test/input scaffold retirement pass, not route surgery.

Task B - supplement seam audit:
Write a concise audit in the response identifying the current repo surfaces for
the existing WL/MWG or GTO supplement machinery and what the PQS H2 physical
route would need next.

The audit should answer:

```text
1. Where does the current WL/MWG residual Gaussian supplement path live?
2. What objects/matrices does it require before the final transform?
3. Which of those facts does the H2 physical PQS route already have?
4. What is the smallest next implementation seam for PQS supplement support?
5. What should remain explicitly forbidden until the seam exists?
```

Use the unified view:

```text
common physical support/intermediate gausslet plan
+ retained transform kind = :wl or :pqs
+ supplement policy = :none or :mwg_residual_gto
```

Do not implement:
- GTO/MWG supplement support;
- H2 supplemented comparison values;
- new matrices or artifacts;
- public API/export/HamV6/CR2/HFDMRG/DMRG behavior;
- changes to the accepted H2 physical endpoint values;
- changes to the visible driver shape beyond deleting stale input/test files.

Line-count rule:
This pass must be net-negative under source/test/bin. Measure:

```sh
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

If deleting the H2 221 diagnostic input/test is allowed, that should be the
line-budget source. Do not satisfy line budget by deleting accepted He/H2
physics endpoint tests.

Validation:
Run the smallest checks that validate the cleanup:

```sh
rg -n "h2_pqs_q5_source_box_diagnostic_r4" test bin src docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
git diff --cached --check
```

If the `rg` command only finds old blurb/review history, that is acceptable.
Do not rerun the slow deleted H2 221 diagnostic test. Do not rerun the H2
physical endpoint unless you edited it.

Response file:
Write `.agent_handoffs/response.222.md` and also copy it to:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.222.md
```

Report:
- whether the H2 221 diagnostic input/test were deleted or why not;
- reference audit result;
- supplement seam audit answers;
- source/test/bin scoped added/deleted/net line count;
- validation commands and results;
- deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

Continue polling after writing the response. Manager will review, commit/push,
and publish the next blurb unless a real blocker appears.

-- repo-manager@macmini

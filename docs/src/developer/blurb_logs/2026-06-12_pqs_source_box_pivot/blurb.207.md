Pass 207 - clarify driver endpoint taxonomy and H2 diagnostic naming

Role:
You are `repo-doer@macmini` implementing one bounded driver/test cleanup pass
for GaussletBases. Follow `AGENTS.md`, `JuliaStyle.md`, and `BlurbStyle.md`.

Loop/approval rule:
- Unattended baton mode is active.
- Do not ask the user for permission through an escalation prompt.
- If a needed command requires approval, write `.agent_handoffs/ATTENTION.md`
  with the exact command, why it is needed, and the blocking condition, then
  stop.

Current state:
- Head before this pass should be:
  `01c3b9b2 Retire product doside shadow blocks`
- The important route taxonomy is:
  - He 419: real physical atom endpoint.
  - H2 221: source-box diagnostic only, not a physics endpoint.
  - H2 463: physical gausslet-only H2 target inventory, currently blocked on
    `:missing_physical_gausslet_source_plan`.

Problem:
The H2 221 diagnostic input/test still use the ambiguous name:

```text
test/driver_inputs/h2_pqs_q5_gausslet_only_r4.jl
test/nested/cartesian_ham_builder_h2_pqs_q5_gausslet_only_r4_readiness_runtests.jl
```

Those names sound like a gausslet-only physics endpoint, but the artifact now
correctly says:

```text
artifact_role = :source_box_diagnostic
physics_endpoint_ready = false
physics_endpoint_blocker = :retained_atom_core_interiors_missing
retained_atom_core_interiors = false
source_plan_role = :boundary_source_box_diagnostic
```

Goal:
Make the driver endpoint taxonomy impossible to confuse before any physical H2
source-plan work resumes.

Allowed work:
1. Rename the H2 221 diagnostic input/test to use `source_box_diagnostic` in
   the filename, for example:

```text
test/driver_inputs/h2_pqs_q5_source_box_diagnostic_r4.jl
test/nested/cartesian_ham_builder_h2_pqs_q5_source_box_diagnostic_r4_runtests.jl
```

2. Update references in the renamed test and any non-log docs/tests that point
   to the old filenames.

3. Add a compact driver endpoint manifest doc, for example:

```text
docs/src/developer/cartesian_driver_endpoint_manifest.md
```

The manifest should be a short table, not a framework. Include at least:
- input file;
- test file;
- route role;
- expected/current final dimension;
- endpoint-ready status;
- current blocker;
- whether it is in a default runner.

Rows should include:
- He PQS q5 WL-map physical atom endpoint, dimension 419;
- He PQS q5 WL-map RHF diagnostic/endpoint test, dimension 419;
- H2 PQS q5 source-box diagnostic, dimension 221, blocker
  `:retained_atom_core_interiors_missing`;
- H2 PQS q5 physical gausslet target, expected dimension 463, blocker
  `:missing_physical_gausslet_source_plan`.

4. Make the renamed H2 221 diagnostic test assert role/blocker fields before
   scalar details if it does not already. Keep it clear that H1/RHF on the 221
   route is not a physics endpoint.

Do not:
- implement a physical H2 source plan;
- add final basis/H1/H1-J/RHF behavior;
- compare H2 gausslet-only values to supplemented WL/QW references;
- change driver semantics;
- add new source modules;
- add broad tests or default-runner includes.

Line-count rule:
The active source/test/bin line-count rule still applies. For edits under
`src`, `test`, `bin`, and the CR2 generator script, final tracked diff must be
net-negative:

```sh
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Docs do not count for this budget. If you rename test/input files, make a small
test/input simplification so the scoped diff is line-negative. Do not delete
live endpoint assertions to game the rule; remove duplication or comments only
where the role assertions make them redundant.

Validation:
Run:

```sh
julia --project=. -e 'using GaussletBases; println("load ok")'
```

Run the renamed H2 221 diagnostic test. If it is expected to exceed 60 seconds,
state why it is necessary and report elapsed time using Julia-level timing.

Run the H2 463 target inventory test unless it is expected to exceed 60 seconds;
if skipped, report the exact command and reason.

Check that no non-log references to the old ambiguous filename remain:

```sh
rg -n "h2_pqs_q5_gausslet_only_r4" test bin docs/src/developer --glob '!docs/src/developer/blurb_logs/**'
```

Always run:

```sh
git diff --check
git diff --cached --check
```

after staging if you stage.

Response file:
Write `.agent_handoffs/response.207.md` and also copy it to:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.207.md
```

Report:
- files renamed/created/edited;
- endpoint manifest contents in summary form;
- source/test/bin scoped line budget added/deleted/net;
- validation commands and elapsed times;
- old filename reference search result;
- deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

Stop after writing the response. Manager will review, commit, push, and pause.

-- repo-manager@macmini

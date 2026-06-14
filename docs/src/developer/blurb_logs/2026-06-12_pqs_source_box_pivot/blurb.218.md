Pass 218 - audit the H2 gausslet-only WL reference comparison path

Role:
You are `repo-doer@macmini` doing one bounded comparison-path audit for
GaussletBases. Follow `AGENTS.md`, `JuliaStyle.md`, `BlurbStyle.md`, and the
unattended baton rules in `.agent_handoffs/RUN.md`.

Loop/approval rule:
- Unattended baton mode is active.
- Do not ask the user for permission through an escalation prompt.
- If a needed command requires approval, write `.agent_handoffs/ATTENTION.md`
  with the exact command, why it is needed, and the blocking condition, then
  stop.

Current state:
- Head before this pass should be:
  `f6aa0b4a Add H2 physical RHF diagnostic execution`
- The physical H2 q5 gausslet-only PQS driver endpoint now has:
  - source plan available;
  - final basis available;
  - H1 available;
  - density/H1-J available;
  - private RHF input contract available from explicit
    `private_rhf_electron_count = 2`;
  - private RHF execution materialized/converged;
  - endpoint still blocked on
    `:missing_h2_gausslet_only_reference_comparison`.

Physics target:
H2 at R = 4.0, gausslet-only, q = n_s = 5, no supplement. The current physical
PQS endpoint has final dimension 463. H2 221 remains source-box diagnostic only.
He 419 is the already validated atom endpoint.

Architectural guardrail:
Treat WL and PQS as sharing the same physical shell/support decomposition
wherever possible. For H2 R=4 q=5, the common physical gausslet-only target is:

```text
support counts  = (275, 578, 362)
retained counts = (251, 98, 114)
retained order  = (:atom_contact_core, :shared_shell_1, :shared_shell_2)
final dimension = 463
supplement      = none
```

WL and PQS should differ in retained transform and fast operator construction
path, not in physical shell/support vocabulary. Do not require the final bases
to be identical as bases; they are different contractions. The comparison path
must be honest about the shared parent/support plan and the different retained
transforms.

Purpose:
Do not implement a new comparison yet. Audit exactly what WL/old-QW route can
serve as the matching gausslet-only H2 reference for the current PQS endpoint,
or conclude that no matching route is currently available. The next code pass
should be based on this audit, not on a guessed reference.

Task:
Read-only audit first. Locate and report:

1. The old WL/QW H2 R=4 gausslet-only reference source, including:
   - where the 463 gausslet-only fixed block came from;
   - whether it was truly supplement-free;
   - whether it used the common support inventory above.
2. Any current visible-driver or low-order materializer path that can reproduce
   the same gausslet-only physical H2 target without supplement.
3. Whether current route-configured WL materialization can be constrained to:
   - R=4 H2;
   - q=n_s=5;
   - no residual GTO/MWG supplement;
   - final dimension 463;
   - the common support/retained inventory above.
4. Whether the existing old WL/QW HF scalar references are supplemented and
   therefore forbidden as direct comparison values for this pass.
5. The exact next implementation seam if a matching WL route exists.
6. The exact blocker if it does not.

Important distinction:
Shared `q=5` and `n_s=5` are not enough. The audit must check the physical
support decomposition, supplement policy, and retained dimension. If the WL
route produces a different dimension or uses residual supplements, it is not the
matching gausslet-only comparison route.

Do not:
- edit `src`, `test`, or `bin` unless you find a clearly stale deletion that is
  safe and still leaves the pass net-negative;
- add a new driver input;
- add a test;
- add supplement support;
- compare to supplemented WL/QW H2 HF totals;
- run HFDMRG, DMRG, CR2, or external solvers;
- mutate the H2 221 diagnostic route;
- hide driver behavior behind an opaque `run_driver(config)` wrapper;
- pass final-basis self-overlap as downstream working data.

Allowed:
- use `rg`, `sed`, `git log`, `git show`, and small local read-only Julia probes
  if needed;
- inspect old reports/docs/handoffs and current driver/materializer code;
- create no tracked source/test changes in this pass.

Line-count rule:
This pass should normally be no-source/no-test. If you do edit tracked files
under `src`, `test`, `bin`, or the CR2 generator script, the scoped diff must be
net-negative:

```sh
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Docs and blurb logs do not count. Do not add source/test just to preserve audit
metadata.

Validation:
This is an audit pass. No broad tests are required. If you run any Julia probe
expected to exceed 60 seconds, explain why first in the response and use
Julia-level timing.

Response file:
Write `.agent_handoffs/response.218.md` and also copy it to:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.218.md
```

Report:
- files/docs inspected;
- old WL/QW gausslet-only reference source and values found;
- whether the old references are supplemented or gausslet-only;
- whether a current visible-driver/materializer path can produce the matching
  WL gausslet-only H2 reference;
- required comparability conditions:
  - same geometry;
  - same physical shell/support plan;
  - no supplement;
  - final dimension 463;
  - explicit retained-transform kind;
  - endpoint comparison label;
- recommended next pass:
  - implementation seam if available, or
  - blocker if not;
- source/test/bin scoped line budget if any source/test/bin files changed;
- deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

Continue polling after writing the response. Manager will review, commit/push,
and publish the next blurb unless a real blocker appears.

-- repo-manager@macmini

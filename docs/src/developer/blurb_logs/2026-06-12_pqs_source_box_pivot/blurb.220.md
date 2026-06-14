Pass 220 - local probe for H2 WL gausslet-only reference values

Role:
You are `repo-doer@macmini` doing one bounded local probe/audit for
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
  `5e864d69 Add H2 WL reference candidate readiness`
- The H2 physical PQS endpoint has a no-matrix WL reference candidate available.
- The remaining endpoint blocker is:
  `:missing_wl_h2_gausslet_only_reference_values`.

Physics target:
H2 at R = 4.0, gausslet-only, q = n_s = 5, no supplement. H2 221 remains
source-box diagnostic only.

Purpose:
Find out whether the repo can currently compute or recover the matching
WL/old-QW gausslet-only 463 reference values without supplement. Do not wire the
values into the driver yet. This pass should produce evidence for the next
implementation decision.

Task:
Use a local ignored probe under `tmp/work/`, or a read-only audit if a probe is
not possible, to answer:

1. Can existing WL/old-QW code build the no-supplement H2 R=4 q5 fixed block
   with final dimension 463?
2. Can it produce gausslet-only WL H1, H1-J/self-Coulomb, and private RHF
   diagnostic values for that exact 463-dimensional no-supplement basis?
3. If yes, report the values with provenance and timing.
4. If no, report the exact missing callable/object/transform needed.

Required comparability conditions:

```text
geometry: H2 R=4, centers equivalent to (0,0,-2) and (0,0,2)
parent axis counts: (9, 9, 15)
support counts: (275, 578, 362)
retained counts: (251, 98, 114)
retained order: (:atom_contact_core, :shared_shell_1, :shared_shell_2)
final dimension: 463
supplement policy: :none
retained transform kind: WL/old-QW gausslet retained transform
```

Do not:
- edit `src`, `test`, or `bin` in this pass unless you are deleting stale
  surfaces and staying net-negative;
- wire values into the driver/artifact;
- add tests;
- add GTO/MWG supplement;
- compare to supplemented WL/QW H2 values;
- run HFDMRG, DMRG, CR2, or external solvers;
- mutate the H2 221 diagnostic route;
- pass final-basis self-overlap as downstream working data;
- hide the driver behind an opaque wrapper.

Allowed:
- create or update an ignored `tmp/work/h2_wl_gausslet_only_463_reference_probe.jl`
  or similarly named probe;
- run narrow Julia commands using `julia --project=.` with Julia-level timing;
- inspect old docs/reports/source.

Runtime rule:
If a probe is expected to exceed 60 seconds, state that in the response and
preserve phase timing from the probe. Do not run a broad suite.

Line-count rule:
This pass should normally make no tracked `src`, `test`, or `bin` edits. If it
does, the scoped diff must be net-negative:

```sh
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Ignored `tmp/work` probe files do not count.

Response file:
Write `.agent_handoffs/response.220.md` and also copy it to:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.220.md
```

Report:
- probe path if created;
- exact commands run;
- whether final dimension 463 was reproduced;
- whether support/retained inventory matched the candidate;
- WL H1/H1-J/RHF values if obtained;
- timing by phase if available;
- blocker if values were not obtained;
- recommended next pass:
  - add values to the driver artifact if values are reviewed, or
  - implement the missing WL route seam if not;
- tracked source/test/bin line budget if any tracked files changed;
- deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

Continue polling after writing the response. Manager will review, commit/push,
and publish the next blurb unless a real blocker appears.

-- repo-manager@macmini

Purpose:
Design the RHF contract boundary after driver-owned private H1/J materialization.

Why now:
The source-box driver path now reaches a private H1/J diagnostic with route-owned
density inputs. Before any RHF implementation, we need a separate route-owned
density/HF contract so the H1/J diagnostic helper does not turn into an
acceptance driver by inertia.

Governing framework:
Use `docs/src/developer/pqs_source_box_operator_framework.md`,
`docs/src/developer/successor_handoff_2026-06-12_pqs_source_box_pivot.md`, and
the methods background packet:
`~/Dropbox/chatarchive/handoff/software_packets/gausslet-methods-fundamentals.md`.

Keep these boundaries sharp:

- source-box-first PQS is the algorithmic framing;
- shell/support-row contraction is oracle/debug;
- retained diagnostic weights are not IDA/quadrature weights;
- H1/J is diagnostic/private unless explicitly promoted;
- RHF/SCF/Fock must be a separate design boundary, not an extension of H1/J by
  accident.

Loop rule:
Do not request interactive approval/escalation during the baton loop. Use
approved commands and read-only inspection. If approval would be required, write
`ATTENTION.md` with the exact command, reason, and blocker, then stop.

Exact task:
Do a no-edit RHF contract design audit.

Inspect:

- `src/pqs_multilayer_complete_core_shell_h1.jl`
- `src/pqs_multilayer_support_density.jl`
- `src/pqs_source_box_route_driver_helpers.jl`
- `docs/src/developer/pqs_source_box_operator_framework.md`
- `docs/src/developer/successor_handoff_2026-06-12_pqs_source_box_pivot.md`
- recent blurb log responses/reviews around passes 097-101.

Answer:

1. What density gauge should RHF use, and how does it relate to the current
   pre-final positive-weight H1/J diagnostic gauge?
2. What orbital/final-to-prefinal map is needed for Fock construction?
3. Where should Fock construction live: H1/J helper, new PQS RHF seam, support
   density module, or another module?
4. What are the required RHF inputs and outputs?
5. What electron count / occupation convention should be explicit before code?
6. What energy components must be reported?
7. What convergence/SCF object should exist before any route adoption?
8. What conventions are forbidden?
9. What should remain diagnostic/private?
10. What is the smallest safe first implementation pass after design, if any?

Do not edit files.
Do not implement RHF.
Do not add Fock construction.
Do not add tests.
Do not run Julia unless needed for a narrow read-only query.
Do not add GTO, IDA/MWG, exports, artifacts, fixture promotion, or production
route behavior.

Validation:
Read-only inspection only. Run:

`git status --short --branch`

Report back:

- proposed RHF contract;
- proposed module/seam ownership;
- explicit density/orbital/occupation conventions;
- forbidden shortcuts;
- smallest next implementation pass or reason to stop before implementation;
- git status;
- deletion/shrinkage forecast:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

After writing `.agent_handoffs/response.102.md`, continue polling for
`blurb.103.md`, `STOP.md`, or `ATTENTION.md`.

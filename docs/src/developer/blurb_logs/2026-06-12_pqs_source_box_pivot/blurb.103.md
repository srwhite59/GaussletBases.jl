Purpose:
Design the fixture/science rule boundary before any RHF implementation.

Why now:
The source-box driver now materializes a private H1/J diagnostic, and pass 102
defined an RHF contract boundary. Before coding RHF, we need to decide what the
current compact fixture means scientifically and what remains only route-smoke
or oracle/debug evidence.

Governing framework:
Use `docs/src/developer/pqs_source_box_operator_framework.md`,
`docs/src/developer/successor_handoff_2026-06-12_pqs_source_box_pivot.md`, and
`~/Dropbox/chatarchive/handoff/software_packets/gausslet-methods-fundamentals.md`.

Keep these boundaries sharp:

- source-box-first PQS is the algorithmic framing;
- shell/support-row contraction is oracle/debug;
- retained diagnostic weights are not IDA/quadrature weights;
- compact H1/J fixture results are diagnostics unless explicitly promoted;
- RHF/SCF/Fock remains unimplemented in this pass.

Loop rule:
Do not request interactive approval/escalation during the baton loop. Use
approved commands and read-only inspection. If approval would be required, write
`ATTENTION.md` with the exact command, reason, and blocker, then stop.

Exact task:
Do a no-edit fixture-rule design audit.

Inspect:

- `test/nested/pqs_direct_retained_final_h1_runtests.jl`
- `src/pqs_source_box_route_driver_helpers.jl`
- recent pass logs/reviews around 095-102
- any existing docs that discuss fixture acceptance, side13/q-ladder, H1/J, or
  PQS source-box route validation.

Answer:

1. What is the current compact fixture?
   Include `Z`, `q`, `n_s`, parent axis counts, core/source dimensions, final
   dimension, H1 energy, self-Coulomb, and density gauge.
2. Which facts are route-smoke facts?
3. Which facts are scientific endpoint/physics facts?
4. Which facts are oracle/debug only?
5. What parameter families must move together for physics comparison?
   Consider `Z`, `d`/spacing, distortion `s`, `q`, parent radius/axis count,
   shell depth, core side/source box size, and side13 versus compact fixtures.
6. What should not become an acceptance gate by inertia?
7. What fixture policy should exist before RHF implementation?
8. What is the smallest next implementation/design pass after this audit?

Do not edit files.
Do not implement RHF.
Do not add tests.
Do not run broad tests.
Do not add GTO, IDA/MWG, exports, artifacts, fixture promotion, or production
route behavior.

Validation:
Read-only inspection only. Run:

`git status --short --branch`

Report back:

- fixture-rule proposal;
- route-smoke versus physics endpoint distinction;
- what should remain oracle/debug;
- what must be decided before RHF;
- smallest next pass recommendation;
- git status;
- deletion/shrinkage forecast:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

After writing `.agent_handoffs/response.103.md`, continue polling for
`blurb.104.md`, `STOP.md`, or `ATTENTION.md`.

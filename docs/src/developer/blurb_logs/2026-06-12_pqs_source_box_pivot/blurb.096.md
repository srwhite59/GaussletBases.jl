Purpose:
Audit the exact H1/J missing-input seam before implementing density/J diagnostic
inputs.

Why now:
The complete core/shell diagnostic route now has a compact internal payload
boundary. The next handoff task is to move H1/J from blocked diagnostic state
toward materialized diagnostic state, but the density inputs are scientifically
sensitive and must use the right convention.

Governing framework:
Use `docs/src/developer/pqs_source_box_operator_framework.md`.

Keep these boundaries sharp:

- source-box-first PQS is the algorithmic framing;
- shell/support-row contraction is oracle/debug;
- retained diagnostic weights are not IDA/quadrature weights;
- H1/J remains diagnostic/private until explicitly promoted.

Loop rule:
Do not request interactive approval/escalation during the baton loop. Use
approved commands and narrow inspection. If approval would be required, stop and
report `ATTENTION` instead.

Exact task:
Do a no-edit audit of the H1/J missing-input seam.

Inspect:

- `src/pqs_source_box_route_driver_helpers.jl`
- `src/pqs_multilayer_complete_core_shell_h1.jl`
- `test/nested/pqs_direct_retained_final_h1_runtests.jl`
- `docs/src/developer/pqs_source_box_operator_framework.md`
- `docs/src/developer/successor_handoff_2026-06-12_pqs_source_box_pivot.md`

Answer:

1. Where are `:axis_weights` and `:raw_pair_factor_terms` currently required?
2. What exact existing objects already contain the needed axis/support weights?
3. What exact existing objects already contain the raw support pair numerator
   terms, if any?
4. Which convention is required for the H1/J diagnostic density path?
5. What is the smallest implementation seam to pass those inputs into
   `_pqs_source_box_route_driver_complete_core_shell_h1_j_diagnostic_payload`?
6. What must remain blocked if any input is not already available as structured
   route-owned data?

Do not edit files.
Do not commit.
Do not run broad tests.
Do not implement density/J inputs yet.
Do not add RHF/SCF/Fock.
Do not add GTO, IDA/MWG, exports, artifacts, fixture promotion, or production
route behavior.
Do not use retained diagnostic weights as quadrature/IDA weights.
Do not promote shell/support-row contraction beyond oracle/debug.

Validation:
Read-only inspection only. If a focused Julia probe is truly needed, explain why
and keep it local under `tmp/work`.

Report back:

- concise inventory of required inputs and current sources;
- convention statement for the diagnostic density path;
- smallest recommended implementation pass;
- exact blocker if structured data is missing;
- git status;
- deletion/shrinkage forecast:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

After writing `.agent_handoffs/response.096.md`, continue polling for
`blurb.097.md`, `STOP.md`, or `ATTENTION.md`.

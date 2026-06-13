Purpose:
Tighten the internal H1/J diagnostic payload contract after density-input
materialization.

Why now:
Pass 097 made the private H1/J diagnostic materialize from route-owned density
inputs. The density-input object is currently used transiently inside the
complete core/shell diagnostic route payload builder, but the compact route
payload should carry it as internal route-owned state before we move to probe
cleanup or broader review.

Governing framework:
Use `docs/src/developer/pqs_source_box_operator_framework.md`.

Keep these boundaries sharp:

- source-box-first PQS is the algorithmic framing;
- shell/support-row contraction is oracle/debug;
- retained diagnostic weights are not IDA/quadrature weights;
- H1/J remains diagnostic/private until explicitly promoted.

Loop rule:
Do not request interactive approval/escalation during the baton loop. Use
approved commands and focused validation. If approval would be required, write
`ATTENTION.md` with the exact command, reason, and blocker, then stop.

Exact task:
Make a tiny internal cleanup in `src/pqs_source_box_route_driver_helpers.jl`.

1. Add `density_inputs` as an internal field on
   `_PQSCompleteCoreShellDiagnosticRoutePayload`.
2. Return the density-input object from
   `_pqs_source_box_route_driver_complete_core_shell_diagnostic_route_payload(...)`
   along with the existing source payload, final basis, H1 payload, and H1/J
   payload.
3. Keep report-facing fields and aliases unchanged.
4. Do not add scalar report-field clouds for density inputs.
5. If a compact metadata flag is needed, make it clearly diagnostic/private,
   for example `private_h1_j_diagnostic_materialized`, but do not rename or
   reinterpret existing report aliases in this pass.

Also inspect whether any docstring or nearby comment now falsely implies that
`driver_route_materialized == true` means public route/global matrix
materialization. If such wording exists in the changed area, correct it
narrowly. Do not do a broad documentation pass.

Trust boundary:
Do not change numerical behavior.
Do not change H1/J values.
Do not add public API.
Do not add new driver options.
Do not add RHF/SCF/Fock.
Do not add GTO, IDA/MWG, exports, artifacts, or fixture promotion.
Do not change lattice size.
Do not add caching/checkpointing.
Do not add tests unless an existing focused assertion must be adjusted for the
new internal field.

Validation:
Run:

`julia --project=. -e 'using GaussletBases; println("load ok")'`
`git diff --check`

If the edit touches only the struct/return shape and no behavior, do not rerun
the long route-driver dry-run. If you see behavior changes, stop and report.

Report back:

- files changed;
- exact field added and where it is populated;
- confirmation report-facing fields are unchanged;
- confirmation numerical behavior was not changed;
- validation commands/results;
- git status;
- deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

After writing `.agent_handoffs/response.098.md`, continue polling for
`blurb.099.md`, `STOP.md`, or `ATTENTION.md`.

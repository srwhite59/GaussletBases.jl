Purpose:
Implement the narrow private density-input adapter for the complete core/shell
H1/J diagnostic route.

Why now:
Pass 096 confirmed that the route-owned source plan already carries the parent
axis bundles needed for the diagnostic density inputs. The remaining H1/J
blockers are `:axis_weights` and `:raw_pair_factor_terms`. The next step is to
feed those inputs from structured source-plan provenance, not from scalar report
aliases or shell/support-row authority.

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
In `src/pqs_source_box_route_driver_helpers.jl`, add a narrow private helper
beside the complete core/shell diagnostic payload helpers, for example:

`_pqs_source_box_route_driver_complete_core_shell_density_inputs(...)`

The helper should:

1. Accept the route-owned source plan and Coulomb expansion, directly or through
   the existing source-plan payload.
2. Require an available `:pqs_multilayer_shell_source_plan`.
3. Require structured source-plan `bundles`.
4. Use the existing validated `_pqs_source_box_ida_factor_provenance(...)`
   path with `expected_term_count = length(coulomb_expansion.coefficients)`.
5. Return only a compact private payload/NamedTuple with:
   - `status`
   - `blocker`
   - `axis_weights`
   - `raw_pair_factor_terms`
   - `missing_inputs`
   - `summary`
   - `metadata` if needed locally

Then wire `_pqs_source_box_route_driver_complete_core_shell_diagnostic_route_payload(...)`
to pass:

- `axis_weights = density_inputs.axis_weights`
- `raw_pair_factor_terms = density_inputs.raw_pair_factor_terms`

into `_pqs_source_box_route_driver_complete_core_shell_h1_j_diagnostic_payload(...)`
when the helper is available.

If density inputs are unavailable, keep H1/J blocked. Do not fall back to
density-normalized pair terms, report scalar aliases, retained diagnostic
weights, fixed-block data, or shell/support-row oracle authority.

Behavior expectation:
For the current shellification-backed source-box driver path, H1/J may now
materialize as a private diagnostic if all structured inputs are present. Route
global materialization, RHF/SCF/Fock, exports, artifacts, and public route
promotion must remain false/unimplemented.

Trust boundary:
Do not add public API.
Do not add new driver options.
Do not add RHF/SCF/Fock.
Do not add GTO, IDA/MWG, exports, artifacts, or fixture promotion.
Do not change lattice size.
Do not add caching/checkpointing.
Do not add scalar report-field clouds.
Preserve existing report aliases unless a focused assertion must update because
the private H1/J diagnostic is now materialized.

Validation:
Run:

`julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl`

Run the smallest focused route-driver dry-run or existing route-driver test that
exercises `cartesian_assembly` and the complete core/shell H1/J diagnostic. Do
not run broad suites.

Run:

`julia --project=. -e 'using GaussletBases; println("load ok")'`
`git diff --check`

Report back:

- files changed;
- helper name and returned fields;
- H1/J diagnostic status before/after for the focused route-driver dry-run;
- final dimension, H1 energy, and self-Coulomb/J diagnostic value if
  materialized;
- confirmation that route/global/RHF/export/artifact flags remain false;
- validation commands/results;
- git status;
- deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

After writing `.agent_handoffs/response.097.md`, continue polling for
`blurb.098.md`, `STOP.md`, or `ATTENTION.md`.

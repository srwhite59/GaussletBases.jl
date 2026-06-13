Pass 098 response

Implemented the tiny internal payload-shape cleanup.

Files changed:

- `src/pqs_source_box_route_driver_helpers.jl`
- `.agent_handoffs/response.098.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.098.md`

Exact field added:

- Added `density_inputs` to `_PQSCompleteCoreShellDiagnosticRoutePayload`.

Where populated:

- `_pqs_source_box_route_driver_complete_core_shell_diagnostic_route_payload(...)`
  now passes the local `density_inputs` object into the compact route payload
  constructor immediately before `h1_j_payload`.

Report-facing fields:

- Unchanged.
- No density-input scalar report-field cloud was added.
- Existing report aliases and H1/J report summary fields were not renamed or
  reinterpreted.

Nearby wording/comment inspection:

- Inspected the changed area for wording implying that
  `driver_route_materialized == true` means public/global route matrix
  materialization.
- No nearby docstring or comment with that false implication was present, so no
  wording change was made.

Numerical behavior:

- Not changed. The edit only stores an already-built private diagnostic object
  on the internal compact route payload.
- No H1/J values, density inputs, source-plan construction, or report aliases
  were changed.

Validation commands/results:

- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed: `load ok`
  - precompilation ran because the source file changed
- `git diff --check`
  - passed
- Long route-driver dry-run was not rerun because the edit only touches
  internal struct/return shape and load/diff checks passed.

Git status before writing this response:

```text
## main...origin/main
 M src/pqs_source_box_route_driver_helpers.jl
```

Deletion/shrinkage report:

- deleted: none
- simplified: compact complete core/shell route payload now carries
  `density_inputs` as route-owned internal state instead of leaving it as a
  transient local only.
- quarantined: none
- not deleted because: current report aliases remain compatibility surfaces and
  the pass was explicitly no-behavior-change.
- exact remaining caller/blocker: no deletion blocker introduced; the remaining
  route-driver materialization blocker from pass 097 is unchanged:
  `cartesian_assembly` does not yet provide an available terminal
  shellification/lowering source-plan path for the focused
  `route_family=:pqs_source_box` route-driver probe.

-- repo-doer@macmini

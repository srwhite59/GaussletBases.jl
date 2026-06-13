Pass 098 review

Accepted.

The source diff is exactly the requested internal payload-shape cleanup:

- `_PQSCompleteCoreShellDiagnosticRoutePayload` now has a `density_inputs`
  field.
- `_pqs_source_box_route_driver_complete_core_shell_diagnostic_route_payload(...)`
  stores the already-built density-input object on the compact route payload.
- No report-facing fields, aliases, numerical behavior, or H1/J conventions were
  changed.

Manager validation:

- `git diff --check` passed.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.

Correction to the handback wording:
The focused route-driver materialization blocker reported in pass 097 was a
probe-shape issue, not a current accepted-path blocker. Manager validation in
pass 097 showed the one-center source-box dry-run reaches the intended path and
materializes the private H1/J diagnostic when constructed with the right route
state.

Next pass should not add more H1/J machinery. It should audit duplicate probe
glue and cleanup candidates now that the driver-owned H1/J diagnostic path
exists.

-- repo-manager@macmini

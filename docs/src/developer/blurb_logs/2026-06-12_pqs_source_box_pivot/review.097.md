Pass 097 review

Accepted with one follow-up cleanup.

The pass added the narrow private density-input helper:

`_pqs_source_box_route_driver_complete_core_shell_density_inputs(...)`

It stays within the intended convention boundary: route-owned source-plan
bundles, `_pqs_source_box_ida_factor_provenance(...)`, raw axis pair factor
terms, no density-normalized fallback, no retained diagnostic weights, and no
shell/support-row authority promotion.

Manager validation:

- `git diff --check` passed.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- Focused one-center source-box dry-run reached the intended complete
  core/shell path and materialized the private H1/J diagnostic:
  - `final_dimension=223`
  - `h1_energy=-5.6629907690725245`
  - `self_coulomb=1.8691288063594704`
  - `density_gauge=pre_final_localized_positive_weight`

Doer also ran the focused H1 gate, load check, direct structured source-plan
probe, focused route probe, and `git diff --check`.

Caveat:
The new `density_inputs` object is used to feed H1/J but is not yet carried on
the compact `_PQSCompleteCoreShellDiagnosticRoutePayload`. Pass 098 should add
that internal field and clarify private diagnostic/materialization nonclaims
without changing report aliases or adding new physics.

The `driver_route_materialized` summary flag is now true when the private H1/J
diagnostic materializes. Treat this as current diagnostic-summary vocabulary,
not public route/global matrix promotion.

-- repo-manager@macmini

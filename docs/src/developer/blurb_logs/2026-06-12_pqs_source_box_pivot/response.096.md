Pass 096 response

No-edit audit completed for the complete core/shell H1/J missing-input seam.

Files changed:

- none

Read-only inspection surfaces:

- `src/pqs_source_box_route_driver_helpers.jl`
- `src/pqs_multilayer_complete_core_shell_h1.jl`
- `src/pqs_multilayer_support_density.jl`
- `src/pqs_multilayer_shell_source_plan.jl`
- `src/CartesianContractedParentMetrics.jl`
- `test/nested/pqs_direct_retained_final_h1_runtests.jl`
- `docs/src/developer/pqs_source_box_operator_framework.md`
- `docs/src/developer/successor_handoff_2026-06-12_pqs_source_box_pivot.md`

Required-input inventory:

1. `:axis_weights` and `:raw_pair_factor_terms` are required at the driver
   diagnostic boundary:
   `_pqs_source_box_route_driver_complete_core_shell_h1_j_missing_inputs(...)`
   marks them missing, and
   `_pqs_source_box_route_driver_complete_core_shell_h1_j_diagnostic_payload(...)`
   forwards them to `pqs_multilayer_complete_core_shell_h1_j_payload(...)`.
   The current compact diagnostic route payload does not yet supply either
   input, so H1/J remains blocked on exactly those density/J inputs.

2. The H1/J payload consumes those inputs through support-density helpers:
   `pqs_multilayer_support_weights(plan; axis_weights)` builds support-row
   product weights in `core_support_states` then `shell_support_states` order,
   and `pqs_multilayer_support_pair_raw_numerator_matrix(plan;
   raw_pair_factor_terms, coulomb_expansion)` builds the raw support pair
   numerator matrix. The raw numerator helper explicitly does not divide by
   weights and does not use density-normalized `pair_factor_terms` as
   authority.

3. Existing route-owned objects already contain the axis/support weight source:
   the shellification-backed `pqs_multilayer_shell_source_plan` carries
   `bundles` and `metrics`. The parent axis bundle PGDG intermediates carry
   `weights`; the existing IDA provenance extractor
   `_pqs_source_box_ida_factor_provenance(bundles)` exposes them as
   `axis_weights = (x = ..., y = ..., z = ...)`.

4. Existing objects also already contain raw support pair numerator terms:
   the same parent axis bundle PGDG intermediates carry
   `pair_factor_terms_raw`; `_pqs_source_box_ida_factor_provenance(bundles)`
   validates and exposes them as `raw_axis_pair_factor_terms = (x = ..., y =
   ..., z = ...)`. The density-normalized `axis_pair_factor_terms` are present
   too, but they are not the correct authority for this H1/J diagnostic path.

Convention statement:

- The diagnostic density path must use the localized pre-final positive-weight
  gauge. The support-row raw numerator matrix and positive support-row weights
  are passed to
  `CartesianFinalBasisRealization.pqs_complete_core_shell_pre_final_density_interaction(...)`.
  That helper forms pre-final weights from `pre_final_coefficients' *
  support_weights`, divides pre-final coefficients by those positive pre-final
  weights, and builds the pre-final pair matrix. It must not use signed final
  weight division, raw no-division density, density-normalized pair terms as
  authority, retained diagnostic weights, RHF/SCF/Fock, GTO, exports,
  artifacts, or production route behavior.

Smallest recommended implementation pass:

- Add a narrow private density-input helper beside the existing complete
  core/shell diagnostic payload helpers, for example
  `_pqs_source_box_route_driver_complete_core_shell_density_inputs(source_plan)`.
- It should accept the available route-owned `source_plan`, require an
  available `source_plan.bundles`, call or mirror the existing validated
  `_pqs_source_box_ida_factor_provenance(source_plan.bundles; expected_term_count
  = length(coulomb_expansion.coefficients))`, and return only:
  `status`, `blocker`, `axis_weights`, `raw_pair_factor_terms`,
  `missing_inputs`, and a compact `summary`.
- Then wire `_pqs_source_box_route_driver_complete_core_shell_diagnostic_route_payload(...)`
  to pass those two fields into
  `_pqs_source_box_route_driver_complete_core_shell_h1_j_diagnostic_payload(...)`.
- Preserve the current report aliases and keep H1/J private diagnostic-only.

What must remain blocked:

- If `source_plan` is absent, blocked, or does not carry route-owned `bundles`,
  keep blocking on `:pqs_multilayer_shell_source_plan` or a specific
  `:source_plan_bundles`/axis-bundle blocker.
- If any axis PGDG intermediate lacks `weights`, keep blocking on
  `:axis_weights` or the specific missing axis weight source.
- If any axis PGDG intermediate lacks `pair_factor_terms_raw`, keep blocking
  on `:raw_pair_factor_terms` or the specific missing axis raw terms.
- If term counts do not match the Coulomb expansion, or raw/density-normalized
  reconstruction/shape validation fails in the provenance extractor, keep H1/J
  blocked rather than falling back to density-normalized pair terms or shell/
  support-row oracle authority.

Validation:

- Read-only inspection only.
- No Julia probe was needed.
- `git status --short --branch`:

```text
## main...origin/main
```

Deletion/shrinkage forecast:

- deleted: none in the next pass unless a local missing-input adapter becomes
  unreachable.
- simplified: the diagnostic route payload can stop manually leaving
  `axis_weights` and `raw_pair_factor_terms` absent when `source_plan.bundles`
  already owns validated PGDG provenance.
- quarantined: shell/support-row contraction remains diagnostic/oracle;
  density-normalized pair terms remain non-authority for this path; H1/J
  remains private diagnostic.
- not deleted because: current report aliases and focused H1 tests still need
  compatibility until the H1/J diagnostic is materialized and reviewed.
- exact remaining caller/blocker: `_pqs_source_box_route_driver_complete_core_shell_diagnostic_route_payload(...)`
  currently calls `_pqs_source_box_route_driver_complete_core_shell_h1_j_diagnostic_payload(...)`
  without `axis_weights` or `raw_pair_factor_terms`; those two fields are the
  remaining structured-input blocker.

-- repo-doer@macmini

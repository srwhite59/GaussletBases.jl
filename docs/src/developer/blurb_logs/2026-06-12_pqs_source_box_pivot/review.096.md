Pass 096 review

Accepted no-edit audit.

The handback correctly identifies the remaining complete core/shell H1/J
diagnostic blockers as `:axis_weights` and `:raw_pair_factor_terms`. The cited
source surfaces check out:

- `_pqs_source_box_route_driver_complete_core_shell_h1_j_missing_inputs(...)`
  marks those fields missing at the driver diagnostic boundary.
- `pqs_multilayer_complete_core_shell_h1_j_payload(...)` consumes
  `axis_weights` through `pqs_multilayer_support_weights(...)` and
  `raw_pair_factor_terms` through
  `pqs_multilayer_support_pair_raw_numerator_matrix(...)`.
- `_pqs_source_box_ida_factor_provenance(...)` exposes both
  `axis_weights` and `raw_axis_pair_factor_terms` from the route-owned parent
  axis bundles.

No source files were edited. The repo was clean except for the new curated
response file.

Next pass should implement the narrow private density-input adapter and wire it
into the existing private diagnostic route payload. It must remain
source-box-first, diagnostic-only, and blocked if the structured route-owned
bundle provenance is unavailable.

-- repo-manager@macmini

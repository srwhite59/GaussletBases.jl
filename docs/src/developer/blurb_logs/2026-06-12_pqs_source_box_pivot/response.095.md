Pass 095 response

Implemented the driver-owned complete core/shell final-basis and H1 seam for
the one-center `:pqs_source_box` route.

Files changed:

- `src/pqs_source_box_route_driver_helpers.jl`

What changed:

- Added a local H1 construction payload in `cartesian_assembly(...)` that reuses
  the pass-094 complete core/shell source-plan payload.
- Built `pqs_multilayer_complete_core_shell_final_basis(source_plan; ...)` from
  the route-owned shellification/lowering source plan.
- Built `pqs_multilayer_complete_core_shell_h1_payload(...)` from route-owned
  center records, route-owned Coulomb expansion, and axis layers obtained from
  the parent PGDG axis bundle.
- Fed the resulting final basis and H1 payload into the existing H1/J diagnostic
  helper.
- Kept dense final-basis/H1 objects internal to assembly; the public report only
  carries the existing compact summary fields.

Focused dry-run result:

```text
h1j_status=blocked_missing_complete_core_shell_h1_j_route_inputs
h1j_blocker=missing_complete_core_shell_h1_j_route_inputs
h1j_summary=(status = :blocked_missing_complete_core_shell_h1_j_route_inputs,
             blocker = :missing_complete_core_shell_h1_j_route_inputs,
             final_dimension = 223,
             h1_energy = -5.6629907690725245,
             self_coulomb = nothing,
             density_gauge = nothing,
             missing_inputs = (:axis_weights, :raw_pair_factor_terms),
             support_density_input_source = :not_materialized,
             h1_orbital_source = :not_materialized,
             signed_final_weight_division_used = false,
             raw_no_division_used = false,
             density_normalized_pair_terms_used_as_authority = false,
             driver_route_materialized = false,
             ida_data_materialized = false,
             density_density_materialized = false,
             rhf_materialized = false,
             gto_materialized = false,
             exports_materialized = false,
             artifacts_materialized = false)
h1_energy=-5.6629907690725245
final_dimension=223
shellization_source=terminal_cartesian_shellification_geometry
terminal_selected=true
driver_route_materialized=false
```

Status:

- Driver-owned final basis is now available for the one-center PQS source-box
  dry-run.
- Driver-owned H1 payload is now available and summarized as final dimension
  `223` with H1 energy `-5.6629907690725245`.
- H1/J remains correctly blocked before J/RHF on only the density/J inputs:
  `(:axis_weights, :raw_pair_factor_terms)`.
- No explicit-box bridge, fixed-block authority, RHF/SCF, GTO, export, or
  artifact path was added.

Validation:

- Focused one-center PQS source-box dry-run: passed with H1/J blocked only on
  `axis_weights` and `raw_pair_factor_terms`.
- `julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl`:
  passed, 46/46.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`: passed.
- `git diff --check`: passed.

Deletion/shrinkage report:

- No old code path became obsolete in this pass; this was a connection seam from
  the already-built route-owned source plan to the existing final-basis/H1
  helpers.
- No tests were added. The existing focused H1 helper test remains the live
  contract check.
- No compatibility adapter was added for the explicit-box bridge or fixed-block
  route.
- Remaining stale/duplicate surface to retire next: none introduced here. The
  next real missing route-owned inputs are the retained density weights and raw
  pair numerator terms for the J diagnostic.

-- repo-doer@macmini

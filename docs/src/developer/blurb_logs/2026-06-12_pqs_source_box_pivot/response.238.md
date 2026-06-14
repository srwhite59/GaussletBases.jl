# Pass 238 response - independent H2 PQS source-plan materializer audit

No files were edited. This was a read-only audit.

## Recommended Source-Plan Object Shape

Introduce a private route-owned descriptor payload in
`src/pqs_source_box_diatomic_complete_core_shell.jl`, not a public API:

```text
_PQSIndependentH2PhysicalSourcePlanPayload
_pqs_source_box_route_driver_independent_h2_physical_source_plan_payload(...)
```

At this stage it should be metadata/readiness only. It should not carry
coefficient matrices.

Recommended fields:

```text
status
blocker
route_family
route_kind
source_plan_family = :independent_pqs_physical_source_box_descriptor
source_plan_authority_status = :independent_pqs_route_owned_source_plan_descriptor
support_plan
retained_rule_plan
support_order = (:atom_contact_core, :shared_shell_1, :shared_shell_2)
retained_order = (:atom_contact_core, :shared_shell_1, :shared_shell_2)
support_counts = (275, 578, 362)
retained_counts = (275, 98, 98)
expected_final_dimension = 471
unit_source_descriptors
unit_retained_rule_descriptors
source_backed_fixed_source_oracle_used = false
fake_pqs_enabled = false
source_coefficients_materialized = false
final_basis_materialized = false
missing_objects
summary
metadata
```

The important split is:

- source-plan descriptor: available after support regions and retained rules are
  described;
- source coefficient / shell realization materialization: still blocked and not
  in the same object as if it were available.

## Exact Implementation Seam

The next implementation should happen in the current driver assembly seam:

```text
cartesian_assembly(...)
  -> _pqs_source_box_route_driver_diatomic_physical_gausslet_target_payload(...)
  -> independent source-plan descriptor payload
  -> existing source-plan payload/reporting
```

But first, gate the fake/source-backed candidate away from independent PQS.
Currently `cartesian_assembly` always builds:

```text
_pqs_source_box_route_driver_diatomic_physical_gausslet_source_plan_candidate_payload(...)
```

and that helper calls:

```text
bond_aligned_diatomic_nested_fixed_source(parent_basis; nside = 5)
```

This is the source-backed WL/QW candidate path. It should not run for
`:bond_aligned_diatomic_independent_pqs_source_box_core_shell`. The independent
route should instead call a new independent descriptor payload and leave the
fake/source-backed candidate for the fake-PQS regression path only.

## Atom-Contact Core

Direct source modes are enough for the next source-plan pass, provided they are
represented as route-owned descriptors and not as coefficient matrices.

The atom-contact unit is a route aggregation of three direct support pieces:

```text
left atom-local core
right atom-local core
midpoint slab
```

Existing authority:

- `CartesianShellification.raw_terminal_geometry(...)` produces the direct
  regions and support counts.
- `CartesianTerminalLowering._direct_terminal_contract(...)` marks direct
  regions with:

```text
retained_rule = :direct_source_modes
realization_rule = :direct_or_trivial_embedding
metadata.identity_like = true
metadata.source_cpb_equals_owned_support = true
```

- `CartesianRetainedUnits._direct_retained_unit(...)` is the metadata precedent
  for an identity-like direct retained unit.

The route-owned H2 source-plan descriptor should own the aggregation into one
`:atom_contact_core` unit. The lower modules own the direct-region facts; they
should not own the H2 cross-region aggregation.

Minimal descriptor for `:atom_contact_core`:

```text
unit_key = :atom_contact_core
source_family = :direct_terminal_source_modes
source_region_roles = (:atom_local_core, :atom_local_core, :midpoint_slab)
source_cpbs = direct contracts' source CPBs
support_count = 275
source_mode_count = 275
retained_rule = :direct_source_modes
retained_count = 275
transform_kind = :identity_source_modes
coefficient_matrix_materialized = false
```

No explicit identity matrix is needed in the next pass.

## Shared Shells

For `:shared_shell_1` and `:shared_shell_2`, use the layered route modules, not
the old fixed-source route.

Existing reusable authority:

- `CartesianShellification.raw_terminal_geometry(...)` / `shellify(...)`
  produce the shared molecular complete-shell regions and support counts.
- `CartesianTerminalLowering.lower_terminal_regions(..., PQSLowering(q = 5))`
  selects `_pqs_complete_shell_contract(...)` for complete-shell regions.
- `_pqs_complete_shell_contract(...)` supplies:

```text
lowering_kind = :pqs_filled_source_cpb
retained_rule = :pqs_boundary_comx_product_modes
realization_rule = :shell_projection_lowdin
metadata.source_mode_shape = (5, 5, 5)
source_cpbs = (filled outer source CPB,)
```

- `CartesianRawProductSources.raw_product_box_plan(...)` should describe each
  shared shell's filled source CPB and source-mode ordering.
- `CartesianRawProductSources.pqs_boundary_product_mode_retained_rule(...)`
  should supply the q=5 boundary COMX product-mode retained-rule descriptor;
  for `(5, 5, 5)` it gives retained count `98`.
- `CartesianRetainedUnits._pqs_shell_retained_unit(...)` is useful precedent
  for the metadata-only PQS shell retained unit, but the H2 source-plan
  descriptor should carry only compact per-unit summaries, not whole final-basis
  materialization data.

`pqs_multilayer_shell_source_plan(...)` is useful precedent for a
shellification-backed source-plan entry point, but it currently materializes
shell projection / Lowdin realization data (`shell_final_coefficients`). That
is too far for the next independent H2 source-plan descriptor pass.

## Owner Boundaries

Recommended ownership:

- `CartesianShellification`: support geometry and owned support regions only.
- `CartesianTerminalLowering`: direct/PQS lowering contracts and source CPBs.
- `CartesianRawProductSources`: raw product-box source facts and boundary
  product-mode retained-rule descriptors.
- `CartesianRetainedUnits`: compact retained-unit metadata from selected
  lowering contracts.
- `src/pqs_source_box_diatomic_complete_core_shell.jl`: H2 route-owned
  composition of atom-contact + shared-shell descriptors into one physical
  source-plan payload.

Do not move H2 route composition into `CartesianRawProductSources`. That module
intentionally does not own shell realization, route assembly, artifacts, or
cross-region physical source plans.

## Missing Primitives / Blockers

The next pass does not need a new numerical primitive if it stays descriptor
only.

The remaining missing objects after descriptor creation should be named more
narrowly than the current broad blocker, for example:

```text
:missing_independent_pqs_source_plan_numerical_materialization
:missing_independent_pqs_shell_realization_coefficients
```

depending on how manager wants to phrase the next blocker.

The true missing numerical work is later:

- source-axis transform coefficients if needed beyond metadata facts;
- shell projection and Lowdin cleanup for the shared shells;
- final basis assembly;
- H1/H1-J/RHF materialization.

## Forbidden Paths Confirmed Avoided

Do not use these for the independent source-plan descriptor:

- fake-PQS H2 463 source-backed WL/QW route;
- `bond_aligned_diatomic_nested_fixed_source(...)`;
- `_pqs_source_box_route_driver_diatomic_physical_gausslet_source_plan_candidate_payload(...)` for independent PQS;
- `_pqs_source_box_route_driver_physical_gausslet_source_plan_from_candidate(...)`;
- old fixed-source coefficient matrices;
- `:wl_qw_source_backed_retained_transform`;
- `:fake_pqs_private_source_backed_adapter_authority`;
- dense parent or shell-row oracle paths as route authority.

Those paths may remain for the fake-PQS golden regression, but they should be
gated out of the independent route.

## Deletion / Shrink Candidates

Best same-surface candidates for the next implementation:

- Make the source-backed candidate payload conditional so independent PQS never
  calls `bond_aligned_diatomic_nested_fixed_source(...)`; this should reduce
  runtime and remove conceptual pressure from the independent path.
- Split or rename the current generic
  `_pqs_source_box_route_driver_diatomic_physical_gausslet_source_plan_candidate_payload`
  so its fake/source-backed role is explicit.
- Retire independent-route report aliases that exist only to explain the old
  missing atom-core/shared-shell retained-rule blockers.
- After the descriptor exists, shrink candidate/report fields that present
  source-backed candidate status beside independent source-plan status.

Separate cleanup pool, not necessarily the same pass:

- route-shaped density-density producer/consumer sections in
  `src/cartesian_contracted_parent_metrics/current_route_metadata_export.jl`;
- corresponding route-shadow pressure in
  `test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl`,
  preserving only compact density/nuclear convention checks.

Do not delete the fake-PQS H2 463 golden regression yet.

## Next-Pass Recommendation

Implement a source-plan descriptor payload next, but do it as a narrow
metadata-only pass:

1. Gate the fake/source-backed candidate away from independent PQS.
2. Add the independent H2 source-plan descriptor payload in
   `pqs_source_box_diatomic_complete_core_shell.jl`.
3. Populate per-unit descriptors from the existing support plan, retained-rule
   plan, terminal lowering contracts, raw product source plans, and retained
   rule descriptors.
4. Update reporting/artifact fields to show descriptor status and the next
   numerical blocker.
5. Do not build source coefficients, shell realization coefficients, final
   basis, H1, H1-J, RHF, supplements, CR2, exports, or public API.

Stop earlier only if the implementation cannot prevent the independent route
from invoking the source-backed candidate without a broader driver refactor.
That would be an `ATTENTION.md` blocker.

## Validation

Read-only inspection only. No Julia commands and no tests were run.

Current status:

```text
## main...origin/main
```

-- repo-doer@macmini

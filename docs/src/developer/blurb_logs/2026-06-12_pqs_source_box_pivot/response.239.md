# Pass 239 response - independent H2 PQS source-plan descriptor

Implemented the descriptor-only independent H2 PQS physical source-plan
payload and gated the source-backed WL/QW candidate away from the independent
route.

## Files changed

- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `src/pqs_source_box_route_driver_helpers.jl`
- `src/pqs_source_box_route_driver_reporting.jl`
- `test/nested/integration_runtests.jl`

Deleted:

- `test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl`

## Source-Plan Descriptor

Added private descriptor helper:

```text
_pqs_source_box_route_driver_independent_h2_physical_source_plan_descriptor(...)
```

Descriptor status/blocker:

```text
source_plan_descriptor_status = :available_independent_pqs_physical_source_plan_descriptor
source_plan_status            = :blocked_pqs_diatomic_physical_gausslet_core_shell_source_plan
source_plan_blocker           = :missing_independent_pqs_source_plan_numerical_materialization
source_plan_authority_status  = :independent_pqs_route_owned_source_plan_descriptor
source_plan_family            = :independent_pqs_physical_source_box_descriptor
```

The descriptor records:

```text
support_order            = (:atom_contact_core, :shared_shell_1, :shared_shell_2)
retained_order           = (:atom_contact_core, :shared_shell_1, :shared_shell_2)
support_counts           = (275, 578, 362)
retained_counts          = (275, 98, 98)
expected_final_dimension = 471
source_coefficients_materialized = false
final_basis_materialized = false
fake_pqs_enabled = false
source_backed_fixed_source_oracle_used = false
```

Per-unit descriptor summaries:

```text
:atom_contact_core
  source_family = :direct_terminal_source_modes
  source_region_roles = (:atom_local_core, :atom_local_core, :midpoint_slab)
  support_count = 275
  retained_count = 275
  transform_kind = :identity_source_modes
  coefficient_matrix_materialized = false

:shared_shell_1
  source_family = :pqs_filled_source_cpb
  source_mode_dims = (5, 5, 5)
  support_count = 578
  retained_count = 98
  retained_rule = :boundary_comx_product_mode_selection
  transform_kind = :source_mode_column_selector
  coefficient_matrix_materialized = false

:shared_shell_2
  source_family = :pqs_filled_source_cpb
  source_mode_dims = (5, 5, 5)
  support_count = 362
  retained_count = 98
  retained_rule = :boundary_comx_product_mode_selection
  transform_kind = :source_mode_column_selector
  coefficient_matrix_materialized = false
```

## Source-Backed Candidate Gate

`cartesian_assembly(...)` now skips
`_pqs_source_box_route_driver_diatomic_physical_gausslet_source_plan_candidate_payload(...)`
when:

```text
recipe.route_kind === :bond_aligned_diatomic_independent_pqs_source_box_core_shell
```

That prevents the independent route from calling:

```text
bond_aligned_diatomic_nested_fixed_source(parent_basis; nside = 5)
```

The fake-PQS/WL golden regression path is still allowed to use the source-backed
candidate.

## Reporting/Artifact

Added compact target artifact fields:

```text
target/source_plan_descriptor_status
target/source_plan_descriptor_blocker
target/source_plan_family
target/source_coefficients_materialized
```

The independent endpoint blocker now advances to:

```text
:missing_independent_pqs_source_plan_numerical_materialization
```

## Forbidden Surfaces Avoided

No source coefficient matrices, shell-realization coefficient matrices, final
basis, H1, H1-J, RHF, supplements, CR2, exports, public API, fake-PQS/WL
coefficients, or fixed-source retained transforms were added.

## Validation

Package load:

```text
julia --project=. -e 't = @elapsed begin using GaussletBases; println("load ok") end; println("elapsed_s=", t)'
```

Result:

```text
load ok
elapsed_s=0.670290375
```

Focused independent-input artifact/readiness check:

```text
julia --project=. -e 'using JLD2; ... include("bin/cartesian_ham_builder.jl") ...'
```

The check captured driver stdout and asserted:

```text
!occursin("diatomic.fixed_source", driver_output)
fake_pqs/enabled == false
route/source_backed_fixed_source_oracle_used == false
physics/endpoint_ready == false
physics/endpoint_blocker == :missing_independent_pqs_source_plan_numerical_materialization
target/support_counts == (275, 578, 362)
target/retained_counts == (275, 98, 98)
target/expected_final_dimension == 471
target/source_plan_status == :blocked_pqs_diatomic_physical_gausslet_core_shell_source_plan
target/source_plan_blocker == :missing_independent_pqs_source_plan_numerical_materialization
target/source_plan_authority_status == :independent_pqs_route_owned_source_plan_descriptor
target/source_plan_descriptor_status == :available_independent_pqs_physical_source_plan_descriptor
target/source_plan_family == :independent_pqs_physical_source_box_descriptor
target/source_coefficients_materialized == false
```

Result:

```text
source_plan_descriptor_status=available_independent_pqs_physical_source_plan_descriptor
source_plan_blocker=missing_independent_pqs_source_plan_numerical_materialization
source_plan_authority_status=independent_pqs_route_owned_source_plan_descriptor
physics_endpoint_blocker=missing_independent_pqs_source_plan_numerical_materialization
source_backed_candidate_output_seen=false
elapsed_s=128.857161458
```

The focused check exceeded 60 seconds because it included source edit
precompilation plus first-call driver compilation.

Whitespace:

```text
git diff --check
```

passed.

## Scoped Line Budget

For `src + test + bin`:

```text
185 added / 5671 deleted, net -5486
```

`git diff --numstat -- src test bin`:

```text
145	0	src/pqs_source_box_diatomic_complete_core_shell.jl
32	0	src/pqs_source_box_route_driver_helpers.jl
8	0	src/pqs_source_box_route_driver_reporting.jl
0	1	test/nested/integration_runtests.jl
0	5670	test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl
```

## Git Status

```text
## main...origin/main
 M src/pqs_source_box_diatomic_complete_core_shell.jl
 M src/pqs_source_box_route_driver_helpers.jl
 M src/pqs_source_box_route_driver_reporting.jl
 M test/nested/integration_runtests.jl
 D test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl
```

## Deletion/Shrinkage Accounting

deleted:
- `test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl`
- its include from `test/nested/integration_runtests.jl`

simplified:
- independent route no longer invokes the source-backed WL/QW fixed-source
  candidate path;
- independent endpoint blocker is now the narrower numerical-materialization
  blocker;
- reporting exposes a compact descriptor status/family rather than promoting
  fake/source-backed candidate status beside independent-PQS status.

quarantined:
- fake-PQS H2 463 golden regression remains separate and unchanged;
- descriptor remains metadata/readiness only.

not deleted because:
- fake-PQS guard fields and endpoint coverage are still needed until an
  independent final-basis/H1/H1-J/RHF endpoint exists;
- source coefficient and shell-realization coefficient materialization remain
  future work;
- public/export/supplement paths remain out of scope.

exact remaining caller/blocker:
- `:missing_independent_pqs_source_plan_numerical_materialization`

-- repo-doer@macmini

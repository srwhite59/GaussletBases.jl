Pass 228 response

Implemented the supplement representation authority cleanup without building provider blocks or supplement matrices.

Files changed:
- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `src/pqs_source_box_route_driver_helpers.jl`
- `src/pqs_source_box_route_driver_reporting.jl`
- `test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl`

What changed:
- Removed duplicate representation status/object-kind fields from `_PQSDiatomicPhysicalGaussletSupplementRequestPayload`.
- Removed the duplicate request-side representation summary/report aliases:
  - `supplement_request_representation_status`
  - `supplement_request_representation_object_kind`
- Removed `supplement_request/representation_status` from the saved artifact and the H2 target test request fingerprint.
- Kept `supplement_representation` as the sole authority for representation status/object kind.
- Kept request metadata focused on request/preflight facts:
  - status
  - blocker
  - fixture_label
  - supplement_policy
  - basis_name
  - lmax
  - atom_symbols
  - nuclear_charges
  - bond_axis
  - bond_length
  - required_provider_blocks
  - missing_fact_labels
  - matrices_materialized=false

Current status ownership:
- `supplement_policy = :none`
  - request status: `:not_requested`
  - representation status: `:not_requested`
  - preflight status: `:not_requested`
- `supplement_policy = :mwg_residual_gto`
  - request status: `:available_pqs_physical_gausslet_supplement_request`
  - representation status: `:available_pqs_physical_gausslet_gto_supplement_representation`
  - representation object kind: `:cartesian_gaussian_shell_supplement_representation`
  - representation center count: 2
  - representation orbital count: 18
  - preflight status: `:blocked_pqs_physical_gausslet_mwg_residual_gto_preflight`
  - preflight blocker: `:missing_provider_gto_supplement_blocks`

Provider-block seam audit:
- Available route inputs now:
  - parent object from the driver assembly call
  - axis bundles from the H2 physical source plan
  - H2 physical source plan via `diatomic_physical_gausslet_source_plan_payload.source_plan`
  - H2 physical final basis via `diatomic_physical_gausslet_final_basis_payload.final_basis`
  - matrix-free supplement representation via `diatomic_physical_gausslet_supplement_representation_payload.representation`
  - center metadata via the parent contract/parent center table
  - Coulomb expansion is available in existing physical H1/J helper paths, but is not yet packaged as a supplement provider-block request input
- Provider-local functions already present in `src/CartesianCPBBlockProviders.jl`:
  - `cpb_gto_supplement_local_operator_bundle`
  - `cpb_mixed_gto_overlap_block`
  - `cpb_mixed_gto_kinetic_operator_block`
  - `cpb_mixed_gto_position_operator_block`
  - `cpb_mixed_gto_x2_operator_block`
  - `cpb_mixed_gto_nuclear_by_center_block`
  - `cpb_gto_overlap_operator_block`
  - `cpb_gto_nuclear_by_center_block`
- Still missing before provider blocks can be built:
  - CPB/source-box coverage for H2 physical support units
  - local row/source ordering for each support unit
  - coefficient transform from support rows to retained 463 basis
  - placement/accumulation rule for mixed gausslet-GTO blocks
  - GTO/GTO self-block rule
  - raw moment matrices for MWG residualization
- Existing route-global mixed-GTO machinery is still WL/decomposed-inventory shaped and blocked on facts such as decomposed WL inventories, row coverage, provider-level GTO bundle, and raw moment matrices. It was not promoted to H2 PQS route authority.

No matrices/supplemented values built:
- No provider blocks were materialized.
- No GTO supplement matrices were materialized.
- `matrices_materialized=false` remains on the request and representation surfaces.
- `provider_blocks_materialized=false` remains on the representation surface.

Net line count:

```text
$ git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
0	19	src/pqs_source_box_diatomic_complete_core_shell.jl
0	5	src/pqs_source_box_route_driver_helpers.jl
0	4	src/pqs_source_box_route_driver_reporting.jl
0	4	test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl
```

Net: -32 lines over the requested source/test/bin scope.

Validation:

```text
$ julia --project=. -e 'using GaussletBases; println("load ok")'
load ok

$ julia --project=. -e 'using Test; t = @elapsed include("test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl"); println("elapsed_s=", t)'
Test Summary:                                                        | Pass  Total     Time
cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4 target artifact |  115    115  2m23.2s
elapsed_s=143.57978475

$ git diff --check
passed

$ git diff --cached --check
passed
```

Git status:

```text
$ git status --short --branch
## main...origin/main
 M src/pqs_source_box_diatomic_complete_core_shell.jl
 M src/pqs_source_box_route_driver_helpers.jl
 M src/pqs_source_box_route_driver_reporting.jl
 M test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl
```

Deletion/shrinkage report:
- deleted: duplicate request-side representation status/object-kind fields, aliases, artifact field, and test expectation.
- simplified: `supplement_request` now reports request/preflight facts only; `supplement_representation` is the representation authority.
- quarantined: provider-block materialization remains blocked/private; existing WL/decomposed route-global mixed-GTO machinery was not promoted.
- not deleted because: request still needs to carry required provider-block names and missing provider facts; representation still needs to carry the matrix-free GTO supplement object for the next seam; the H2 endpoint test remains the active artifact gate.
- exact remaining caller/blocker: `:missing_provider_gto_supplement_blocks`; missing route facts are CPB/source-box coverage, local row ordering, support-row-to-retained transform placement, mixed block accumulation, GTO/GTO self-block rule, and raw moment matrices.

-- repo-doer@macmini

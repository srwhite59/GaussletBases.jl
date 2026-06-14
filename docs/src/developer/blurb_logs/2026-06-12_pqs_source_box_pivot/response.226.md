Pass 226 response - H2 PQS supplement request payload

Files changed:
- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `src/pqs_source_box_route_driver_helpers.jl`
- `src/pqs_source_box_route_driver_reporting.jl`
- `test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl`
- `test/nested/pqs_source_box_route_driver_report_runtests.jl`

Implemented:
- Added private payload `_PQSDiatomicPhysicalGaussletSupplementRequestPayload`.
- Wired `cartesian_assembly` to build the request payload before the existing supplement preflight payload.
- Added compact artifact group `supplement_request`.
- Kept all behavior matrix-free: no GTO/GTO, mixed gausslet/GTO, MWG residual, density-density, supplemented scalar, public API, export, CR2, or HFDMRG behavior was added.

Request/preflight behavior:
- `supplement_policy = :none`
  - request status: `:not_requested`
  - request blocker: `nothing`
  - preflight status: `:not_requested`
  - preflight blocker: `nothing`
- `supplement_policy = :mwg_residual_gto`
  - request status: `:blocked_pqs_physical_gausslet_mwg_residual_gto_request`
  - request blocker: `:missing_gto_supplement_representation`
  - preflight status: `:blocked_pqs_physical_gausslet_mwg_residual_gto_preflight`
  - preflight blocker: `:missing_gto_supplement_representation`
  - next downstream blocker remains `:missing_provider_gto_supplement_blocks` after representation exists.

Representation:
- No `CartesianGaussianShellSupplementRepresentation3D` object is created in this pass.
- `representation_status` is explicitly:
  - `:not_requested` for `:none`
  - `:not_materialized_pqs_physical_gausslet_gto_supplement_representation` for `:mwg_residual_gto`

Default H2 request facts:
- `fixture_label = :h2_r4_physical_gausslet_q5`
- `basis_name = "H/cc-pVTZ"`
- `lmax = 1`
- `uncontracted = false`
- atom symbols: `("H", "H")`
- nuclear charges: `(1.0, 1.0)`
- bond axis/length: `:z`, `4.0`
- residual policy: `:route_private_mwg_residual_gto_preflight_only`

Artifact fields written in `supplement_request`:
- `status`
- `blocker`
- `fixture_label`
- `supplement_policy`
- `basis_name`
- `lmax`
- `atom_symbols`
- `nuclear_charges`
- `bond_axis`
- `bond_length`
- `representation_status`
- `required_provider_blocks`
- `missing_fact_labels`
- `matrices_materialized = false`

Test/shrink work:
- Updated the H2 physical endpoint test to check `supplement_request/*`.
- Updated MWG preflight expectations to block first on `:missing_gto_supplement_representation`.
- Compactified repeated H2 `supplement_preflight/*` and target/route metadata assertions while keeping accepted H2 numerical endpoint checks.
- Deleted stale route-driver materialization/report-artifact assertion blocks from `test/nested/pqs_source_box_route_driver_report_runtests.jl`; these preserved transitional private metadata and were not endpoint physics checks.

Source/test/bin scoped line count:

```text
211	1	src/pqs_source_box_diatomic_complete_core_shell.jl
47	0	src/pqs_source_box_route_driver_helpers.jl
40	0	src/pqs_source_box_route_driver_reporting.jl
38	29	test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl
0	751	test/nested/pqs_source_box_route_driver_report_runtests.jl
```

Totals: 336 added, 781 deleted, net -445.

Validation:
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed
  - precompiled `GaussletBases` in about 57 seconds
- `julia --project=. -e 'using Test; t = @elapsed include("test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl"); println("elapsed_s=", t)'`
  - passed: 121/121
  - elapsed: `143.205909792` seconds
- `git diff --check`
  - passed
- `git diff --cached --check`
  - passed
- `julia --project=. -e 'Meta.parse("begin\n" * read("test/nested/pqs_source_box_route_driver_report_runtests.jl", String) * "\nend"); println("route report test parse ok")'`
  - passed
  - parse-only check for the shrunk slow integration test

Not run:
- `test/nested/pqs_source_box_route_driver_report_runtests.jl` was not executed; it is marked as an integration/slow test and this pass only deleted stale transitional assertion pressure there. I ran a parse-only check instead.

Git status:

```text
## main...origin/main
 M src/pqs_source_box_diatomic_complete_core_shell.jl
 M src/pqs_source_box_route_driver_helpers.jl
 M src/pqs_source_box_route_driver_reporting.jl
 M test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl
 M test/nested/pqs_source_box_route_driver_report_runtests.jl
```

Deletion/shrinkage report:
- deleted:
  - 751 lines of stale route-driver materialization/report-artifact assertion pressure from `test/nested/pqs_source_box_route_driver_report_runtests.jl`.
- simplified:
  - repeated H2 supplement preflight assertions replaced with compact request/preflight checks.
  - redundant target/route metadata assertions reduced in the H2 endpoint test.
- quarantined:
  - GTO representation/provider-block materialization remains explicitly blocked and private.
- not deleted because:
  - accepted H2 numerical endpoint checks and WL/PQS deltas remain the live physics gate.
  - provider-local GTO block tests remain needed for the future provider seam.
- exact remaining caller/blocker:
  - `:missing_gto_supplement_representation` blocks the request/preflight first.
  - after representation exists, the next blocker is `:missing_provider_gto_supplement_blocks`.

-- repo-doer@macmini

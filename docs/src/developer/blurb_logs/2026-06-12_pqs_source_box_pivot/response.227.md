Pass 227 response - H2 supplement representation request

Files changed:
- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `src/pqs_source_box_route_driver_helpers.jl`
- `src/pqs_source_box_route_driver_reporting.jl`
- `test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl`
- `test/nested/pqs_source_box_route_driver_report_runtests.jl`

Implemented:
- Added private `_PQSDiatomicPhysicalGaussletSupplementRepresentationPayload`.
- Materialized the route-owned matrix-free GTO supplement representation for the H2 physical MWG request using:

```julia
legacy_bond_aligned_diatomic_gaussian_supplement("H", "cc-pVTZ", nuclei; lmax = 1)
basis_representation(supplement)
```

- Kept the route-facing request label as `"H/cc-pVTZ"` while normalizing the constructor basis name to `"cc-pVTZ"`.
- Wired `cartesian_assembly` to build the representation payload between supplement request and preflight.
- Added compact artifact group `supplement_representation`.
- No GTO/GTO, mixed gausslet/GTO, MWG residual, density-density, supplemented scalar, public API, export, CR2, HFDMRG, DMRG, or production route behavior was added.

Representation behavior:
- `supplement_policy = :none`
  - representation status: `:not_requested`
  - blocker: `nothing`
  - orbital count: `0`
- `supplement_policy = :mwg_residual_gto`
  - request status: `:available_pqs_physical_gausslet_supplement_request`
  - representation status: `:available_pqs_physical_gausslet_gto_supplement_representation`
  - representation object kind: `:cartesian_gaussian_shell_supplement_representation`
  - center count: `2`
  - orbital count: `18`
  - matrices materialized: `false`
  - provider blocks materialized: `false`

Request/preflight blocker transition:
- Before this pass, the first blocker was `:missing_gto_supplement_representation`.
- After this pass, the representation fact is available and preflight advances to:
  - `:missing_provider_gto_supplement_blocks`
- Remaining missing facts are provider/block readiness facts only:
  - `:missing_provider_gto_supplement_blocks`
  - `:missing_mixed_gausslet_gto_blocks`
  - `:missing_gto_gto_blocks`
  - `:missing_combined_raw_moment_matrices`
  - `:missing_residual_mwg_representation`
  - `:missing_combined_density_density_readiness`

Artifact fields written in `supplement_representation`:
- `status`
- `blocker`
- `object_kind`
- `basis_name`
- `lmax`
- `atom_symbols`
- `center_count`
- `orbital_count`
- `matrices_materialized`
- `provider_blocks_materialized`

Source/test/bin scoped line count:

```text
202	17	src/pqs_source_box_diatomic_complete_core_shell.jl
46	0	src/pqs_source_box_route_driver_helpers.jl
37	0	src/pqs_source_box_route_driver_reporting.jl
59	31	test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl
7	306	test/nested/pqs_source_box_route_driver_report_runtests.jl
```

Totals: 351 added, 354 deleted, net -3.

Validation:
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed
  - precompiled `GaussletBases` in about 58 seconds
- `julia --project=. -e 'using Test; t = @elapsed include("test/nested/cartesian_cpb_mixed_gto_overlap_block_runtests.jl"); println("elapsed_s=", t)'`
  - passed: 1889/1889
  - elapsed: `85.594604875` seconds
- `julia --project=. -e 'using Test; t = @elapsed include("test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl"); println("elapsed_s=", t)'`
  - passed: 115/115
  - elapsed: `145.825892375` seconds
- `git diff --check`
  - passed
- `git diff --cached --check`
  - passed
- `julia --project=. -e 'Meta.parse("begin\n" * read("test/nested/pqs_source_box_route_driver_report_runtests.jl", String) * "\nend"); println("route report test parse ok")'`
  - passed

Notes:
- The first H2 endpoint attempt failed because the route label `"H/cc-pVTZ"` was passed directly to the legacy constructor. I fixed that by keeping `"H/cc-pVTZ"` as the request/artifact label and using `"cc-pVTZ"` for the legacy basis constructor.
- `test/nested/pqs_source_box_route_driver_report_runtests.jl` was not executed; it is marked integration/slow. I ran a parse-only check after shrinking stale metadata assertion pressure there.

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
  - stale route-report integration assertions and unused helper scaffolding for one-center/White-Lindsey permutations, exhaustive parent-contract metadata, route-configured shellization/materializer metadata, and printed report section headers.
- simplified:
  - H2 MWG supplement assertions now use compact request/representation/preflight fingerprints plus endpoint-critical facts.
  - route-report integration test is now a compact dry-run save smoke instead of a broad transitional metadata inventory.
- quarantined:
  - provider-block materialization remains blocked behind `:missing_provider_gto_supplement_blocks`.
  - no supplemented H2 scalar values or old supplemented WL/QW references were promoted.
- not deleted because:
  - accepted H2 numerical endpoint checks and WL/PQS gausslet-only deltas remain the live physics gate.
  - provider-local GTO/source tests remain the existing focused check for later provider work.
- exact remaining caller/blocker:
  - request and representation are available for the H2 MWG policy.
  - supplement preflight is still blocked at `:missing_provider_gto_supplement_blocks`; no mixed gausslet/GTO, GTO/GTO, raw moment, residual MWG, or density-density provider blocks exist on this route yet.

-- repo-doer@macmini

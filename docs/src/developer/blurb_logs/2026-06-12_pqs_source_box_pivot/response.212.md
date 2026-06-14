Pass 212 response - H2 physical gausslet final-basis seam

Final-basis seam materialized.

- final-basis object kind: `:pqs_physical_gausslet_final_basis`
- final-basis status: `:available_pqs_physical_gausslet_final_basis`
- final-basis blocker: `nothing`
- endpoint remains not ready: `physics/endpoint_ready = false`
- endpoint blocker advanced to `:missing_physical_gausslet_h1_builder`

Implementation notes:

- Fixed the private physical source-plan object to carry the merged
  atom-contact core support from `source.sequence.core_indices/core_states`,
  not only the first child sequence.
- Added private payload `_PQSDiatomicPhysicalGaussletFinalBasisPayload`.
- Added a private physical final-basis builder that:
  - verifies support/retained order
    `(:atom_contact_core, :shared_shell_1, :shared_shell_2)`;
  - verifies counts `(275, 578, 362)` -> `(251, 98, 114)`;
  - restricts parent-row coefficient maps to route-owned support rows;
  - assembles the support-order transform;
  - computes overlap diagnostics from axis overlap factors;
  - applies the final Lowdin cleanup and records only identity-error
    diagnostics for the artifact/report path.
- Did not reuse the complete-core-shell identity-core helper.
- Did not implement H1, H1-J, density interaction, RHF, supplement, CR2,
  exports, or public API.

Counts and diagnostic:

- support counts: `(275, 578, 362)`
- retained counts: `(251, 98, 114)`
- final dimension: `463`
- `basis/final_overlap_identity_error = 1.6181500583911657e-13`

Artifact fields changed:

- `config/run_final_basis = true`
- `target/source_plan_blocker = nothing`
- `route/final_basis_status = :available_pqs_physical_gausslet_final_basis`
- `basis/final_dimension = 463`
- `basis/final_overlap_identity_error` is written
- `physics/endpoint_blocker = :missing_physical_gausslet_h1_builder`
- endpoint manifest row updated to the H1-builder blocker

Files changed:

- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `src/pqs_source_box_route_driver_helpers.jl`
- `src/pqs_source_box_route_driver_reporting.jl`
- `test/driver_inputs/h2_pqs_q5_physical_gausslet_r4.jl`
- `test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl`
- `test/nested/integration_runtests.jl`
- `test/nested/bond_aligned_diatomic_high_order_recipe_opt_in_source_construction_integration_runtests.jl`
- `docs/src/developer/cartesian_driver_endpoint_manifest.md`

Source/test/bin scoped line budget:

```text
404  5     src/pqs_source_box_diatomic_complete_core_shell.jl
27   0     src/pqs_source_box_route_driver_helpers.jl
29   5     src/pqs_source_box_route_driver_reporting.jl
2    2     test/driver_inputs/h2_pqs_q5_physical_gausslet_r4.jl
0    1644  test/nested/bond_aligned_diatomic_high_order_recipe_opt_in_source_construction_integration_runtests.jl
7    8     test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl
0    1     test/nested/integration_runtests.jl
```

Totals for scoped files: 469 added, 1665 deleted, net -1196.

Validation:

```sh
julia --project=. -e 'using GaussletBases; println("load ok")'
```

Result: passed; precompile reported `56683.1 ms`, then printed `load ok`.

```sh
julia --project=. -e 'using Test; t = @elapsed include("test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl"); println("elapsed_s=", t)'
```

Result: passed, 40/40, `elapsed_s=75.662340541`.

Runtime remains around the expected 75 seconds. The fixed-source phase printed
`diatomic.fixed_source.total: 11.42 seconds`; the added final-basis seam did
not introduce an obvious new dominant runtime phase.

Focused artifact scalar probe:

```sh
julia --project=. -e 'using JLD2; driver = normpath(joinpath(pwd(), "bin", "cartesian_ham_builder.jl")); input = normpath(joinpath(pwd(), "test", "driver_inputs", "h2_pqs_q5_physical_gausslet_r4.jl")); outfile = "/private/tmp/h2_pqs_q5_physical_gausslet_r4_pass212.jld2"; tsvfile = "/private/tmp/h2_pqs_q5_physical_gausslet_r4_pass212.tsv"; saved_args = copy(ARGS); empty!(ARGS); append!(ARGS, [input, "outfile=$(repr(outfile))", "tsvfile=$(repr(tsvfile))"]); t = @elapsed try include(driver) finally empty!(ARGS); append!(ARGS, saved_args) end; jldopen(outfile, "r") do file; println("artifact_probe_elapsed_s=", t); println("final_overlap_identity_error=", file["basis/final_overlap_identity_error"]); println("final_dimension=", file["basis/final_dimension"]); println("route_final_basis_status=", file["route/final_basis_status"]); println("physics_endpoint_blocker=", file["physics/endpoint_blocker"]); end'
```

Result:

```text
artifact_probe_elapsed_s=73.85325
final_overlap_identity_error=1.6181500583911657e-13
final_dimension=463
route_final_basis_status=available_pqs_physical_gausslet_final_basis
physics_endpoint_blocker=missing_physical_gausslet_h1_builder
```

```sh
git diff --check
git diff --cached --check
```

Result: both passed.

Git status:

```text
## main...origin/main
 M docs/src/developer/cartesian_driver_endpoint_manifest.md
 M src/pqs_source_box_diatomic_complete_core_shell.jl
 M src/pqs_source_box_route_driver_helpers.jl
 M src/pqs_source_box_route_driver_reporting.jl
 M test/driver_inputs/h2_pqs_q5_physical_gausslet_r4.jl
 D test/nested/bond_aligned_diatomic_high_order_recipe_opt_in_source_construction_integration_runtests.jl
 M test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl
 M test/nested/integration_runtests.jl
```

Deletion/shrinkage report:

- deleted:
  `test/nested/bond_aligned_diatomic_high_order_recipe_opt_in_source_construction_integration_runtests.jl`
- simplified:
  removed its include from `test/nested/integration_runtests.jl`; the active
  H2 endpoint guard now validates the physical source-plan/final-basis seam
  directly through the driver artifact
- quarantined:
  physical final-basis remains private/internal and carries no H1, H1-J,
  density, RHF, export, supplement, or public API claims
- not deleted because:
  the H2 driver-owned artifact test remains the active physics seam guard; the
  source-backed candidate remains needed as private adapter input until H1/H1-J
  route builders exist
- exact remaining caller/blocker:
  `:missing_physical_gausslet_h1_builder`

-- repo-doer@macmini

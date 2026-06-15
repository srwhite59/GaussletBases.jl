# Pass 247 response - independent H2 PQS source-plan payload

Implemented the complete route-owned source-plan payload assembly for the
independent H2 PQS target. No final-basis, H1, H1-J, RHF, supplement, CR2,
export, or public API work was enabled.

## Files changed

- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `src/pqs_source_box_route_driver_helpers.jl`
- `src/pqs_source_box_route_driver_reporting.jl`
- `test/nested/cartesian_assembly_stage_low_order_policy_runtests.jl`
- `test/nested/cartesian_report_stage_low_order_policy_runtests.jl`
- `test/nested/cartesian_transform_stage_low_order_policy_runtests.jl`
- `test/nested/cartesian_unit_stage_low_order_policy_runtests.jl`

## Source-plan payload

The independent H2 PQS route now assembles an available
`_PQSDiatomicPhysicalGaussletCoreShellSourcePlan` from:

- route-owned atom-contact core support rows derived from the atom-local core
  and midpoint slab boxes;
- pass-241 shared-shell realization payloads for `shared_shell_1` and
  `shared_shell_2`;
- sparse parent-row coefficient maps for the direct atom-contact core selector
  and shared-shell realization coefficients.

Focused artifact result:

```text
source_plan_status=available_pqs_diatomic_physical_gausslet_core_shell_source_plan
source_plan_blocker=nothing
support_counts=(275, 578, 362)
retained_counts=(275, 98, 98)
expected_final_dimension=471
source_plan_authority_status=independent_pqs_route_owned_source_plan
```

Recorded route/source-plan facts include:

```text
support_order = (:atom_contact_core, :shared_shell_1, :shared_shell_2)
retained_order = (:atom_contact_core, :shared_shell_1, :shared_shell_2)
support_counts = (275, 578, 362)
retained_counts = (275, 98, 98)
final_dimension = 471
source_backed_fixed_source_oracle_used = false
fake_pqs_enabled = false
retained_transform_authority = :pqs_source_box_construction
source_plan_authority_status = :independent_pqs_route_owned_source_plan
```

## Fake/source-backed guard evidence

The focused artifact check asserted:

```text
fake_pqs_enabled=false
source_backed_fixed_source_oracle_used=false
route/retained_transform_authority=pqs_source_box_construction
```

The independent driver path still skips the fake/source-backed WL/QW candidate
payload.

## Remaining blocked surfaces

The endpoint remains non-physics-ready:

```text
physics_endpoint_blocker=missing_pqs_physical_gausslet_final_basis
final_basis_status=not_materialized_pqs_physical_gausslet_final_basis
route/h1_materialized=false
route/h1_j_materialized=false
```

No H1, H1-J, density interaction, RHF, supplement provider, MWG/GTO, CR2,
export, artifact-public API readiness, or fake-PQS promotion was added.

## Line budget

Scoped `src + test + bin`:

```text
233   2   src/pqs_source_box_diatomic_complete_core_shell.jl
28    6   src/pqs_source_box_route_driver_helpers.jl
6     2   src/pqs_source_box_route_driver_reporting.jl
0     67  test/nested/cartesian_assembly_stage_low_order_policy_runtests.jl
0     111 test/nested/cartesian_report_stage_low_order_policy_runtests.jl
0     60  test/nested/cartesian_transform_stage_low_order_policy_runtests.jl
0     39  test/nested/cartesian_unit_stage_low_order_policy_runtests.jl
```

Total: `267 added / 287 deleted`, net `-20`.

Deletion offset used:

- removed stale exact terminal-shellification lowering/selected-contract mirror
  assertions in the broad staged unit/transform/assembly/report tests;
- kept compact terminal route/scaffold/materialization assertions.

## Validation

Package load:

```text
julia --project=. -e 't = @elapsed begin using GaussletBases; println("load ok") end; println("elapsed_s=", t)'
```

Result:

```text
load ok
elapsed_s=59.697142584
```

Focused independent H2 PQS driver/artifact assertion:

```text
julia --project=. -e 'using JLD2; outfile = "/private/tmp/h2_independent_pqs_247.jld2"; stdout_path = "/private/tmp/h2_independent_pqs_247_stdout.txt"; old_args = copy(ARGS); empty!(ARGS); append!(ARGS, ["test/driver_inputs/h2_pqs_q5_independent_source_box_r4.jl", "outfile=\"$(outfile)\"", "save_artifact=true", "save_tsv=false"]); elapsed = @elapsed open(stdout_path, "w") do out; redirect_stdout(out) do; include("bin/cartesian_ham_builder.jl"); end; end; empty!(ARGS); append!(ARGS, old_args); jldopen(outfile, "r") do file; @assert file["target/source_plan_status"] === :available_pqs_diatomic_physical_gausslet_core_shell_source_plan; @assert file["target/source_plan_blocker"] === nothing; @assert Tuple(file["target/support_counts"]) == (275, 578, 362); @assert Tuple(file["target/retained_counts"]) == (275, 98, 98); @assert file["target/expected_final_dimension"] == 471; @assert file["route/source_plan_status"] === :available_pqs_diatomic_physical_gausslet_core_shell_source_plan; @assert file["route/source_backed_fixed_source_oracle_used"] == false; @assert file["route/fake_pqs_enabled"] == false; @assert file["route/retained_transform_authority"] === :pqs_source_box_construction; @assert file["target/source_plan_authority_status"] === :independent_pqs_route_owned_source_plan; @assert file["target/shared_shell_realization_status"] === :available_independent_pqs_shared_shell_realization_payload; @assert Tuple(file["target/shared_shell_realization_counts"]) == (98, 98); @assert file["target/source_coefficients_materialized"] == true; @assert file["fake_pqs/enabled"] == false; @assert file["physics/endpoint_ready"] == false; @assert file["physics/endpoint_blocker"] in (:physical_gausslet_final_basis_request_not_enabled, :missing_pqs_physical_gausslet_final_basis); @assert file["route/final_basis_status"] === :not_materialized_pqs_physical_gausslet_final_basis; @assert file["route/h1_materialized"] == false; @assert file["route/h1_j_materialized"] == false; println("source_plan_status=", file["target/source_plan_status"]); println("source_plan_blocker=", file["target/source_plan_blocker"]); println("support_counts=", Tuple(file["target/support_counts"])); println("retained_counts=", Tuple(file["target/retained_counts"])); println("expected_final_dimension=", file["target/expected_final_dimension"]); println("source_plan_authority_status=", file["target/source_plan_authority_status"]); println("fake_pqs_enabled=", file["fake_pqs/enabled"]); println("source_backed_fixed_source_oracle_used=", file["route/source_backed_fixed_source_oracle_used"]); println("physics_endpoint_blocker=", file["physics/endpoint_blocker"]); println("final_basis_status=", file["route/final_basis_status"]); end; println("stdout_path=", stdout_path); println("elapsed_s=", elapsed)'
```

Result:

```text
source_plan_status=available_pqs_diatomic_physical_gausslet_core_shell_source_plan
source_plan_blocker=nothing
support_counts=(275, 578, 362)
retained_counts=(275, 98, 98)
expected_final_dimension=471
source_plan_authority_status=independent_pqs_route_owned_source_plan
fake_pqs_enabled=false
source_backed_fixed_source_oracle_used=false
physics_endpoint_blocker=missing_pqs_physical_gausslet_final_basis
final_basis_status=not_materialized_pqs_physical_gausslet_final_basis
stdout_path=/private/tmp/h2_independent_pqs_247_stdout.txt
elapsed_s=78.41146575
```

This exceeded 60 seconds because it is the only focused driver artifact check
that exercises the independent H2 source-plan assembly seam.

Whitespace:

```text
git diff --check
```

passed.

## Git status

```text
## main...origin/main
 M src/pqs_source_box_diatomic_complete_core_shell.jl
 M src/pqs_source_box_route_driver_helpers.jl
 M src/pqs_source_box_route_driver_reporting.jl
 M test/nested/cartesian_assembly_stage_low_order_policy_runtests.jl
 M test/nested/cartesian_report_stage_low_order_policy_runtests.jl
 M test/nested/cartesian_transform_stage_low_order_policy_runtests.jl
 M test/nested/cartesian_unit_stage_low_order_policy_runtests.jl
```

## Deletion/shrinkage report

deleted:
- stale exact terminal-shellification lowering and selected-contract mirror
  assertion blocks in broad staged tests.

simplified:
- independent H2 source-plan status now comes from one route-owned source-plan
  payload instead of descriptor-plus-blocker adjacency;
- endpoint blocker merging no longer lets unsupported optional H1/RHF fixture
  diagnostics overwrite an existing final-basis blocker.

quarantined:
- fake-PQS/source-backed WL/QW route remains separate and unused by this
  independent path.

not deleted because:
- final-basis builder still has the older physical-gausslet retained-count
  assumptions and was intentionally not changed in this pass;
- broad staged tests were not run because the blurb explicitly excluded them
  as stale validation gates.

exact remaining caller/blocker:
- `:missing_pqs_physical_gausslet_final_basis`

-- repo-doer@macmini

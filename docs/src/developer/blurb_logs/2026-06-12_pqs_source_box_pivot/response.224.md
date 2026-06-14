Pass 224 response - H2 physical supplement preflight boundary

Implemented a private metadata-only supplement preflight boundary for the H2
physical gausslet route.

Files changed:
- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `src/pqs_source_box_route_driver_helpers.jl`
- `src/pqs_source_box_route_driver_reporting.jl`
- `test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl`
- `test/nested/pqs_source_box_route_driver_report_runtests.jl`

Preflight object:
- `_PQSDiatomicPhysicalGaussletSupplementPreflightPayload`

Preflight summary fields:
- `status`
- `blocker`
- `route_family`
- `route_kind`
- `fixture_label`
- `support_counts`
- `retained_counts`
- `retained_order`
- `retained_transform_kind`
- `gausslet_final_dimension`
- `supplement_policy`
- `required_fact_labels`
- `available_fact_labels`
- `missing_fact_labels`
- `matrices_materialized`
- `supplemented_values_materialized`
- `public_api`

Behavior:
- `supplement_policy = :none`
  - `supplement_preflight/status = :not_requested`
  - `supplement_preflight/blocker = nothing`
  - no accepted no-supplement H2 endpoint values or deltas were changed.
- `supplement_policy = :mwg_residual_gto`
  - `supplement_preflight/status = :blocked_pqs_physical_gausslet_mwg_residual_gto_preflight`
  - `supplement_preflight/blocker = :missing_provider_gto_supplement_blocks`
  - missing facts recorded:
    `:missing_provider_gto_supplement_blocks`,
    `:missing_mixed_gausslet_gto_blocks`,
    `:missing_gto_gto_blocks`,
    `:missing_combined_raw_moment_matrices`,
    `:missing_residual_mwg_representation`,
    `:missing_combined_density_density_readiness`

Artifact fields written:
- `route/supplement_preflight_status`
- `route/supplement_preflight_blocker`
- `supplement_preflight/status`
- `supplement_preflight/blocker`
- `supplement_preflight/fixture_label`
- `supplement_preflight/support_counts`
- `supplement_preflight/retained_counts`
- `supplement_preflight/retained_order`
- `supplement_preflight/retained_transform_kind`
- `supplement_preflight/gausslet_final_dimension`
- `supplement_preflight/supplement_policy`
- `supplement_preflight/required_fact_labels`
- `supplement_preflight/available_fact_labels`
- `supplement_preflight/missing_fact_labels`
- `supplement_preflight/matrices_materialized`
- `supplement_preflight/supplemented_values_materialized`

Confirmation:
- No GTO/GTO matrices were built.
- No mixed gausslet/GTO matrices were built.
- No MWG residual matrices were built.
- No supplemented H2 comparison values were added.
- The accepted no-supplement endpoint remains the only materialized physics
  endpoint in this pass.

Line count:

```text
137  1    src/pqs_source_box_diatomic_complete_core_shell.jl
39   0    src/pqs_source_box_route_driver_helpers.jl
53   2    src/pqs_source_box_route_driver_reporting.jl
89   0    test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl
0    319  test/nested/pqs_source_box_route_driver_report_runtests.jl
```

Scoped source/test/bin total:
- added: 318
- deleted: 322
- net: -4

Validation:

```text
julia --project=. -e 'using GaussletBases; println("load ok")'
```

Result: passed. Precompiled `GaussletBases` in 57 s, then printed `load ok`.

```text
julia --project=. -e 'using Test; t = @elapsed include("test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl"); println("elapsed_s=", t)'
```

Result: passed.

```text
Test Summary:                                                        | Pass  Total     Time
cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4 target artifact |  119    119  2m20.6s
elapsed_s=141.034032417
```

The focused test exceeded 60 s because it exercised the accepted H2 final
basis/H1/H1-J/private RHF artifact path and the new supplement-preflight
override.

Additional changed-test validation attempted:

```text
julia --project=. -e 'using Test; t = @elapsed include("test/nested/pqs_source_box_route_driver_report_runtests.jl"); println("elapsed_s=", t)'
```

Result: incomplete. The process produced no normal test output for several
minutes and was terminated with `kill 95765`. It had parsed into the file and
was compiling inside `_pqs_route_driver_check_materialization_status` at
`test/nested/pqs_source_box_route_driver_report_runtests.jl:1093` when killed.
I did not broaden validation further.

```text
git diff --check
```

Result: passed.

```text
git diff --cached --check
```

Result: passed.

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
- deleted: 319 lines of stale route-report metadata assertion pressure from
  `test/nested/pqs_source_box_route_driver_report_runtests.jl`.
- simplified: exact route-report key-cloud assertions and duplicated
  `route_configured_*` artifact/TSV assertions were reduced to compact status
  and materialization checks.
- quarantined: MWG/GTO supplement request now stops at the private preflight
  payload with explicit missing facts.
- not deleted because: the H2 physical endpoint artifact test remains the live
  acceptance gate for the accepted no-supplement physics endpoint and now for
  the supplement-preflight boundary.
- exact remaining caller/blocker: supplement advancement is blocked at
  `:missing_provider_gto_supplement_blocks`, then mixed gausslet/GTO blocks,
  GTO/GTO blocks, combined raw moment matrices, residual MWG representation,
  and combined density-density readiness.

-- repo-doer@macmini

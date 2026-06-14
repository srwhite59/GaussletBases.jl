Pass 210 response - checked H2 physical source-plan candidate adapter

Implemented a compact private checked candidate/readiness seam for the H2 463
physical gausslet source-plan path.

Candidate result:
- candidate status: `:available_physical_gausslet_source_plan_candidate`
- candidate source: `:source_backed_fixed_source_oracle`
- candidate counts match: `true`
- source-plan authority status: `:candidate_not_route_authority`
- actual source plan remains blocked
- actual source-plan blocker: `:source_plan_candidate_not_route_authority`

Candidate checks:
- support order: `(:atom_contact_core, :shared_shell_1, :shared_shell_2)`
- retained order: `(:atom_contact_core, :shared_shell_1, :shared_shell_2)`
- support counts: `(275, 578, 362)`
- retained counts: `(251, 98, 114)`
- final dimension: `463`
- no supplement: true
- H2 221 diagnostic source-plan object reused: false

Implementation:
- Added private payload:
  `_PQSDiatomicPhysicalGaussletSourcePlanCandidatePayload`
- Added helper:
  `_pqs_source_box_route_driver_diatomic_physical_gausslet_source_plan_candidate_payload`
- The helper builds the source-backed fixed-source oracle candidate from the
  carried H2 parent QW basis with `bond_aligned_diatomic_nested_fixed_source`;
  it checks counts/order and stores only compact status/count fields in the
  route/report payload.
- The existing physical source-plan payload now consumes the candidate payload
  and keeps `source_plan = nothing` with blocker
  `:source_plan_candidate_not_route_authority`.

Artifact fields added:
- `target/source_plan_candidate_status`
- `target/source_plan_candidate_source`
- `target/source_plan_candidate_counts_match`
- `target/source_plan_authority_status`

Artifact behavior:
- `physics/endpoint_ready` remains `false`.
- `physics/endpoint_blocker` is now
  `:source_plan_candidate_not_route_authority`.
- No final basis, H1, H1-J, density interaction, RHF, HFDMRG, CR2, export, or
  public API behavior was added.

Files changed:
- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `src/pqs_source_box_route_driver_helpers.jl`
- `src/pqs_source_box_route_driver_reporting.jl`
- `test/driver_inputs/h2_pqs_q5_physical_gausslet_r4.jl`
- `test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl`
- `test/nested/bond_aligned_diatomic_high_order_recipe_opt_in_source_construction_integration_runtests.jl`

Source/test/bin scoped line budget:

```text
143  0    src/pqs_source_box_diatomic_complete_core_shell.jl
17   0    src/pqs_source_box_route_driver_helpers.jl
12   0    src/pqs_source_box_route_driver_reporting.jl
1    1    test/driver_inputs/h2_pqs_q5_physical_gausslet_r4.jl
0    213  test/nested/bond_aligned_diatomic_high_order_recipe_opt_in_source_construction_integration_runtests.jl
9    2    test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl
```

Totals: 182 added, 216 deleted, net -34.

Validation:
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed
  - precompiled `GaussletBases` in about 56 s, then printed `load ok`
- `julia --project=. -e 'using Test; t = @elapsed include("test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl"); println("elapsed_s=", t)'`
  - passed: 40/40
  - elapsed: `74.511691083` s
  - timing impact: source-backed candidate construction added about 11.3 s,
    dominated by `diatomic.fixed_source.source_assembly` at about 10.1 s;
    assembly total was about 27.2 s
- `git diff --check`
  - passed
- `git diff --cached --check`
  - passed

Git status:

```text
## main...origin/main
 M src/pqs_source_box_diatomic_complete_core_shell.jl
 M src/pqs_source_box_route_driver_helpers.jl
 M src/pqs_source_box_route_driver_reporting.jl
 M test/driver_inputs/h2_pqs_q5_physical_gausslet_r4.jl
 M test/nested/bond_aligned_diatomic_high_order_recipe_opt_in_source_construction_integration_runtests.jl
 M test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl
```

Deletion/shrinkage report:
- deleted:
  - none
- simplified:
  - removed 213 lines of old non-default high-order opt-in source-construction
    test pressure that checked detailed projected-q-shell descriptor,
    source-mode dimension, and metadata vocabulary
  - kept the test's active construction, dimension, support coverage,
    finite-output, symmetry, overlap, and route dimension checks
- quarantined:
  - the candidate remains explicitly oracle/adapter-only through
    `:source_backed_fixed_source_oracle` and
    `:candidate_not_route_authority`
- not deleted because:
  - the H2 463 artifact-readiness test remains the active guard for target
    inventory and candidate blocker behavior
  - source-backed fixed-source construction is still needed as the checked
    oracle candidate until a route-owned physical source-plan producer exists
- exact remaining caller/blocker:
  - source plan remains blocked on
    `:source_plan_candidate_not_route_authority`; the next real producer must
    turn the verified source-backed rows/coefs into a route-owned physical PQS
    source-plan object before final-basis or H1 work starts

-- repo-doer@macmini

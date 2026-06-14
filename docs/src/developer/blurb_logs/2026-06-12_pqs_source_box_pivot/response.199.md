Pass 199 response - label current H2 221 route as diagnostic-only.

Files changed:

- `bin/cartesian_ham_builder.jl`
- `src/pqs_source_box_route_driver_helpers.jl`
- `src/pqs_source_box_route_driver_reporting.jl`
- `test/driver_inputs/h2_pqs_q5_gausslet_only_r4.jl`
- `test/nested/cartesian_ham_builder_h2_pqs_q5_gausslet_only_r4_readiness_runtests.jl`
- deleted `test/nested/cartesian_pair_block_one_body_placement_plan_runtests.jl`

Exact artifact labels added:

- `route/artifact_role = :source_box_diagnostic`
- `physics/endpoint_ready = false`
- `physics/endpoint_blocker = :retained_atom_core_interiors_missing`
- `basis/retained_atom_core_interiors = false`
- `basis/source_plan_role = :boundary_source_box_diagnostic`

Driver/input surface:

- Added optional driver route-input fields:
  - `artifact_role`
  - `physics_endpoint_ready`
  - `physics_endpoint_blocker`
  - `retained_atom_core_interiors`
  - `source_plan_role`
- The H2 input sets those fields to the diagnostic-only labels above.
- The fields are carried through `cartesian_recipe`, durable recipe metadata, and the diatomic readiness artifact writer.

Diagnostic-only confirmation:

- Current 221-dimensional route construction is unchanged.
- H2 readiness test still asserts:
  - `basis/final_dimension == 221`
  - `comparison/ready == false`
  - `route/h1_j_materialized == false`
  - `private_rhf/materialized == false`
- No H1-J, density interaction, private RHF, supplemented comparison, source-plan redesign, export, or public solver work was added.

Validation and timing:

- `julia --project=. -e 't = @elapsed include("test/nested/cartesian_ham_builder_h2_pqs_q5_gausslet_only_r4_readiness_runtests.jl"); println("h2_diagnostic_label_elapsed_s=", t)'`
  - passed after wiring labels through durable recipe metadata
  - 52 tests passed
  - `h2_diagnostic_label_elapsed_s=153.783415166`
  - precompilation took about 57 seconds
  - route phases remained compilation dominated:
    - transform about 0.27s, about 99.9% compilation
    - pair terms about 1.04s, about 99.97% compilation
    - assembly about 33.5s, about 99.2% compilation
- `julia --project=. test/nested/cartesian_ham_builder_h2_pqs_q5_gausslet_only_r4_readiness_runtests.jl`
  - passed
  - 52 tests passed
  - test summary time `1m34.9s`
  - no package precompilation on this run; assembly was still about 33.1s and about 99.4% compilation
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed: `load ok`
- `git diff --check`
  - passed with no output
- `rg -n "cartesian_pair_block_one_body_placement_plan_runtests|CartesianPairBlockMaterialization local one-body placement plan" test src bin`
  - no output; exit 1 from no matches
- `git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`
  - `7	0	bin/cartesian_ham_builder.jl`
  - `16	0	src/pqs_source_box_route_driver_helpers.jl`
  - `6	1	src/pqs_source_box_route_driver_reporting.jl`
  - `5	0	test/driver_inputs/h2_pqs_q5_gausslet_only_r4.jl`
  - `7	0	test/nested/cartesian_ham_builder_h2_pqs_q5_gausslet_only_r4_readiness_runtests.jl`
  - `0	301	test/nested/cartesian_pair_block_one_body_placement_plan_runtests.jl`
- `git status --short --branch`
  - `## main...origin/main`
  - modified:
    - `bin/cartesian_ham_builder.jl`
    - `src/pqs_source_box_route_driver_helpers.jl`
    - `src/pqs_source_box_route_driver_reporting.jl`
    - `test/driver_inputs/h2_pqs_q5_gausslet_only_r4.jl`
    - `test/nested/cartesian_ham_builder_h2_pqs_q5_gausslet_only_r4_readiness_runtests.jl`
  - deleted:
    - `test/nested/cartesian_pair_block_one_body_placement_plan_runtests.jl`

Scoped line-budget arithmetic:

- Insertions: 41
- Deletions: 302
- Net: -261 lines

Deletion/shrinkage report:

- deleted:
  - `test/nested/cartesian_pair_block_one_body_placement_plan_runtests.jl`
- simplified:
  - removed a non-runner private one-body placement-plan scaffold while preserving the smaller mixed one-body consumer/materialization coverage.
  - H2 artifact now plainly marks the 221-dimensional route as diagnostic-only.
- quarantined:
  - none
- not deleted because:
  - the H2 readiness test remains the active driver-owned artifact gate;
  - explicit He/H2 driver endpoints, WL/H2 references, and scientific endpoints were preserved;
  - route construction was intentionally unchanged in this labeling pass.
- exact remaining caller/blocker:
  - H2 physics endpoint remains blocked by `:retained_atom_core_interiors_missing`; current 221-dimensional route is a source-box diagnostic only.

-- repo-doer@macmini

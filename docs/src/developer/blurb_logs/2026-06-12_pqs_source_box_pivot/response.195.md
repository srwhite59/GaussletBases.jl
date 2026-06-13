Pass 195 complete.

Files changed:

- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `src/pqs_source_box_route_driver_reporting.jl`
- `test/driver_inputs/h2_pqs_q5_gausslet_only_r4.jl`
- `test/nested/cartesian_ham_builder_h2_pqs_q5_gausslet_only_r4_readiness_runtests.jl`
- deleted `test/nested/cartesian_pair_block_one_body_dispatch_runtests.jl`

Parent/axis-bundle path connected:

- Reused the existing parent-axis probe path:
  `_pqs_explicit_core_spacing_parent_axis_probe(...; carry_objects = true)`.
- The H2 driver input now explicitly sets:
  - `parent_axis_counts = nothing`
  - `probe_parent_axis_construction = true`
- This lets the existing `_cartesian_parent_object_carry` path carry the
  transient parent basis, wrapped QW basis, and axis-bundle objects into the
  parent contract.
- No second parent constructor path was added.

Artifact/report additions:

- Added a small `parent/*` group to the H2 readiness artifact writer.
- The explicit H2 artifact test now verifies:
  - `parent/parent_axis_counts == (x = 9, y = 9, z = 15)`
  - `parent/parent_axis_counts_source === :constructed_parent_axis_probe`
  - `parent/parent_materialization_blocker === nothing`
  - `parent/parent_basis_object_available == true`
  - `parent/parent_qw_basis_object_available == true`
  - `parent/parent_axis_bundle_object_available == true`
  - `parent/parent_basis_object_type_label == "CartesianParentGaussletBasis3D"`
  - `parent/parent_qw_basis_object_type_label == "BondAlignedDiatomicQWBasis3D"`
  - `parent/parent_axis_bundle_object_type_label == "_CartesianNestedAxisBundles3D"`

Axis-count result:

- The constructed parent-axis probe reports `(x = 9, y = 9, z = 15)`.
- These counts match the current H2 target extents:
  - transverse extent `4.0`
  - parallel extent `6.0`
  - `core_spacing = 0.5`
- The old manual fixture count is no longer used as route authority for this
  input.

Route/source-plan readiness:

- Enabling the real parent/axis-bundle carry made the already-existing private
  source-plan path available.
- Current artifact route status:
  - `route/readiness_status === :blocked_diatomic_complete_core_shell_ham_readiness`
  - `route/readiness_blocker === :missing_diatomic_complete_core_shell_final_basis_consumer`
  - `route/source_plan_status === :available_pqs_diatomic_complete_core_shell_source_plan`
  - `route/final_basis_status === :not_materialized_diatomic_complete_core_shell_final_basis`
  - `route/h1_status === :not_materialized_diatomic_complete_core_shell_h1`
- I added a narrow final-basis request guard so this readiness-only H2 input
  does not cascade into final-basis/H1 materialization while `run_h1`,
  `run_h1_j`, and `run_private_rhf` are all false.
- No H1-J, density interaction, private RHF, supplemented comparison, export,
  or public solver behavior was added.

Scoped line-budget arithmetic:

```text
9	0	src/pqs_source_box_diatomic_complete_core_shell.jl
1	0	src/pqs_source_box_route_driver_reporting.jl
2	2	test/driver_inputs/h2_pqs_q5_gausslet_only_r4.jl
16	2	test/nested/cartesian_ham_builder_h2_pqs_q5_gausslet_only_r4_readiness_runtests.jl
0	261	test/nested/cartesian_pair_block_one_body_dispatch_runtests.jl
```

Total for the requested scoped diff: 28 added, 265 deleted, net `-237`.

Validation results:

- `julia --project=. test/nested/cartesian_ham_builder_h2_pqs_q5_gausslet_only_r4_readiness_runtests.jl`
  - passed: 38/38
  - elapsed: 1m31.0s on final run
  - final run had no package precompile; assembly was about 31s and
    compilation-dominated for the carried parent/object path
  - artifact path printed by the final run:
    `/Users/srw/dmrgtmp/jl_iNulGb/h2_pqs_q5_gausslet_only_r4.jld2`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed: `load ok`
- `git diff --check`
  - passed
- `git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`
  - passed and shown above
- `git status --short --branch`
  - `## main...origin/main`
  - ` M src/pqs_source_box_diatomic_complete_core_shell.jl`
  - ` M src/pqs_source_box_route_driver_reporting.jl`
  - ` M test/driver_inputs/h2_pqs_q5_gausslet_only_r4.jl`
  - ` M test/nested/cartesian_ham_builder_h2_pqs_q5_gausslet_only_r4_readiness_runtests.jl`
  - ` D test/nested/cartesian_pair_block_one_body_dispatch_runtests.jl`

Deletion/shrinkage report:

- deleted:
  - `test/nested/cartesian_pair_block_one_body_dispatch_runtests.jl`
- simplified:
  - removed an unrunner-included metadata-only mixed one-body dispatch scaffold
  - kept the active mixed one-body consumer smoke and accessor/plan-batch
    coverage intact
  - removed the stale manual H2 parent-axis-count authority from the input by
    setting `parent_axis_counts = nothing`
- quarantined:
  - nothing
- not deleted because:
  - the H2 readiness test is the current driver-owned endpoint-readiness gate
  - He driver endpoints and scientific/reference tests remain live
- exact remaining caller/blocker:
  - H2 PQS route readiness now blocks at
    `:missing_diatomic_complete_core_shell_final_basis_consumer`
  - the source plan is available, but final-basis/H1 remain not materialized
    because this input has `run_h1=false`, `run_h1_j=false`, and
    `run_private_rhf=false`
  - supplemented WL/QW H2 comparison remains blocked until supplement policy is
    represented explicitly in the driver/PQS path

-- repo-doer@macmini

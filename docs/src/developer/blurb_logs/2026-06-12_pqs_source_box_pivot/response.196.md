Pass 196 complete.

Files changed:

- `bin/cartesian_ham_builder.jl`
- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `src/pqs_source_box_route_driver_helpers.jl`
- `src/pqs_source_box_route_driver_reporting.jl`
- `test/driver_inputs/h2_pqs_q5_gausslet_only_r4.jl`
- `test/nested/cartesian_ham_builder_h2_pqs_q5_gausslet_only_r4_readiness_runtests.jl`
- deleted `test/nested/cartesian_pair_block_one_body_plan_dispatch_runtests.jl`

Final-basis request surface:

- Added explicit visible driver input `run_final_basis`.
- Default driver value is `nothing`; `cartesian_recipe` resolves `nothing` to
  the existing downstream requests:
  `run_h1 || run_h1_j || run_private_rhf`.
- The H2 input sets:
  - `run_final_basis = true`
  - `run_h1 = false`
  - `run_h1_j = false`
  - `run_private_rhf = false`
- The H2 artifact now records `config/run_final_basis = true`.

Final-basis materialization:

- Reused the existing private diatomic complete-core/shell final-basis payload
  path.
- Added a final-basis request gate so final basis can be requested without
  using `run_h1` as a backdoor.
- Added a separate H1 request gate so final-basis-only H2 does not solve H1.
- No second final-basis implementation path was added.

Final-basis result:

- `route/source_plan_status === :available_pqs_diatomic_complete_core_shell_source_plan`
- `route/final_basis_status === :available_pqs_complete_core_shell_final_basis`
- `basis/final_dimension = 221`
- `basis/final_overlap_identity_error = 2.510171831978689e-13`
- `route/h1_status === :not_materialized_diatomic_complete_core_shell_h1`
- `route/h1_j_materialized == false`
- `private_rhf/materialized == false`
- `comparison/ready == false`

New route blocker:

- `route/readiness_status === :blocked_diatomic_complete_core_shell_ham_readiness`
- `route/readiness_blocker === :missing_diatomic_complete_core_shell_h1_consumer`

Validation results and timing:

- `julia --project=. test/nested/cartesian_ham_builder_h2_pqs_q5_gausslet_only_r4_readiness_runtests.jl`
  - passed: 41/41
  - elapsed: 2m30.9s
  - runtime was compilation dominated:
    - package precompile: about 57s
    - assembly/final-basis path: about 32s, 99.94% compilation time
  - artifact path printed by the test:
    `/Users/srw/dmrgtmp/jl_5AtRuk/h2_pqs_q5_gausslet_only_r4.jld2`
- Focused readback command:
  - `julia --project=. -e 'using JLD2; ... include("bin/cartesian_ham_builder.jl"); ...'`
  - confirmed:
    - `final_dimension=221`
    - `final_overlap_identity_error=2.510171831978689e-13`
    - `readiness_blocker=missing_diatomic_complete_core_shell_h1_consumer`
  - readback artifact:
    `/Users/srw/dmrgtmp/jl_SfoNH5/h2_pqs_final_basis_readback.jld2`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed: `load ok`
- `git diff --check`
  - passed
- `git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`
  - passed:

```text
2	1	bin/cartesian_ham_builder.jl
12	1	src/pqs_source_box_diatomic_complete_core_shell.jl
14	4	src/pqs_source_box_route_driver_helpers.jl
1	1	src/pqs_source_box_route_driver_reporting.jl
1	0	test/driver_inputs/h2_pqs_q5_gausslet_only_r4.jl
5	2	test/nested/cartesian_ham_builder_h2_pqs_q5_gausslet_only_r4_readiness_runtests.jl
0	343	test/nested/cartesian_pair_block_one_body_plan_dispatch_runtests.jl
```

Scoped total: 35 added, 352 deleted, net `-317`.

- `git status --short --branch`
  - `## main...origin/main`
  - ` M bin/cartesian_ham_builder.jl`
  - ` M src/pqs_source_box_diatomic_complete_core_shell.jl`
  - ` M src/pqs_source_box_route_driver_helpers.jl`
  - ` M src/pqs_source_box_route_driver_reporting.jl`
  - ` M test/driver_inputs/h2_pqs_q5_gausslet_only_r4.jl`
  - ` M test/nested/cartesian_ham_builder_h2_pqs_q5_gausslet_only_r4_readiness_runtests.jl`
  - ` D test/nested/cartesian_pair_block_one_body_plan_dispatch_runtests.jl`

Deletion/shrinkage report:

- deleted:
  - `test/nested/cartesian_pair_block_one_body_plan_dispatch_runtests.jl`
- simplified:
  - removed an unrunner-included plan-level mixed one-body dispatch metadata
    scaffold
  - kept active mixed one-body consumer smoke and endpoint/readiness tests
  - kept final-basis request explicit instead of overloading `run_h1`
- quarantined:
  - nothing
- not deleted because:
  - H2 readiness test is the live driver-owned H2 endpoint-readiness gate
  - He driver endpoint and scientific/reference tests remain live
- exact remaining caller/blocker:
  - H2 PQS route readiness now blocks at
    `:missing_diatomic_complete_core_shell_h1_consumer`
  - final basis is available, but H1/H1-J/RHF remain not materialized because
    the H2 input explicitly leaves those requests false
  - supplemented WL/QW H2 comparison remains blocked until supplement policy is
    represented explicitly in the driver/PQS path

-- repo-doer@macmini

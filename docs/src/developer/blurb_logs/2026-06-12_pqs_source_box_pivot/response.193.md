Pass 193 complete.

Implemented a driver-visible H2 gausslet-only readiness artifact without
making a supplemented WL/QW physics comparison and without promoting H1/H1-J or
private RHF for H2.

Files changed:

- `bin/cartesian_ham_builder.jl`
  - added visible/default override keys for `bond_axis`, `bond_length`,
    `xmax_parallel`, `xmax_transverse`, `supplement_policy`,
    `comparison_ready`, `comparison_blocker`, `run_h1`, and `run_h1_j`
  - kept the visible staged construction intact
- `src/pqs_source_box_route_driver_helpers.jl`
  - carried the new driver metadata into system/recipe metadata
  - added compact report-facing diatomic complete core/shell readiness fields
    sourced from the existing private diatomic payloads
- `src/pqs_source_box_route_driver_reporting.jl`
  - added a blocker-aware H2/diatomic readiness artifact writer gated on
    two-center PQS with `supplement_policy = :none`
  - left the existing He H1/H1-J artifact writer unchanged
- `test/driver_inputs/h2_pqs_q5_gausslet_only_r4.jl`
  - new explicit H2 gausslet-only driver input
- `test/nested/cartesian_ham_builder_h2_pqs_q5_gausslet_only_r4_readiness_runtests.jl`
  - new explicit/manual driver readiness artifact test, not added to default
    or integration runners
- deleted:
  - `test/nested/cartesian_pair_block_one_body_block_set_preflight_runtests.jl`
  - `test/nested/cartesian_pair_block_one_body_batch_summary_runtests.jl`

Artifact path produced by focused test:

```text
/Users/srw/dmrgtmp/jl_Pl7mmH/h2_pqs_q5_gausslet_only_r4.jld2
```

H2 artifact fields verified:

- `system/atom_symbols == ("H", "H")`
- `system/nuclear_charges == (1, 1)`
- `system/atom_locations == ((0.0, 0.0, -2.0), (0.0, 0.0, 2.0))`
- `system/bond_axis == :z`
- `system/bond_length == 4.0`
- `config/route_family == :pqs_source_box`
- `config/route_kind == :bond_aligned_diatomic_fixed_q_complete_core_shell`
- `config/q == 5`
- `config/n_s == 5`
- `config/supplement_policy == :none`
- `config/comparison_ready == false`
- `comparison/ready == false`
- `comparison/blocker == :supplemented_reference_not_comparable_to_gausslet_only`
- no `comparison/wl_rhf_total`
- no `comparison/delta_rhf`
- no private RHF materialization
- route/global/export/public flags checked false where exposed

Route/materialization readiness recorded:

- `route/readiness_status == :blocked_diatomic_complete_core_shell_ham_readiness`
- `route/readiness_blocker == :missing_diatomic_complete_core_shell_source_plan_producer`
- `route/source_plan_status == :not_materialized_diatomic_complete_core_shell_source_plan`
- `route/final_basis_status == :not_materialized_diatomic_complete_core_shell_final_basis`
- `route/h1_status == :not_materialized_diatomic_complete_core_shell_h1`

This is intentional readiness reporting. The artifact honestly records that
the H2 gausslet-only driver target is configured but not yet a materialized H2
H1/H1-J/RHF endpoint. It does not compare to the documented supplemented WL/QW
H2 HF/ED references.

Scoped line-budget arithmetic:

- requested tracked numstat:
  - `14	4	bin/cartesian_ham_builder.jl`
  - `144	0	src/pqs_source_box_route_driver_helpers.jl`
  - `79	0	src/pqs_source_box_route_driver_reporting.jl`
  - `0	206	test/nested/cartesian_pair_block_one_body_batch_summary_runtests.jl`
  - `0	276	test/nested/cartesian_pair_block_one_body_block_set_preflight_runtests.jl`
- requested `git diff --numstat` subtotal: `237` added, `486` deleted, net `-249`
- new untracked files not shown by `git diff --numstat`:
  - `32` lines in `test/driver_inputs/h2_pqs_q5_gausslet_only_r4.jl`
  - `72` lines in `test/nested/cartesian_ham_builder_h2_pqs_q5_gausslet_only_r4_readiness_runtests.jl`
- real scoped budget including new files: `341` added, `486` deleted, net `-145`

Validation results:

- `julia --project=. test/nested/cartesian_ham_builder_h2_pqs_q5_gausslet_only_r4_readiness_runtests.jl`
  - passed: `29` assertions, `43.9s` on final run
  - earlier first run took `1m42.5s`, dominated by about `57s` package
    precompilation
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed: `load ok`
- `git diff --check`
  - passed
- `git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`
  - passed and shown above
- `git status --short --branch`
  - `## main...origin/main`
  - ` M bin/cartesian_ham_builder.jl`
  - ` M src/pqs_source_box_route_driver_helpers.jl`
  - ` M src/pqs_source_box_route_driver_reporting.jl`
  - ` D test/nested/cartesian_pair_block_one_body_batch_summary_runtests.jl`
  - ` D test/nested/cartesian_pair_block_one_body_block_set_preflight_runtests.jl`
  - `?? test/driver_inputs/h2_pqs_q5_gausslet_only_r4.jl`
  - `?? test/nested/cartesian_ham_builder_h2_pqs_q5_gausslet_only_r4_readiness_runtests.jl`

Deletion/shrinkage report:

- deleted:
  - `test/nested/cartesian_pair_block_one_body_block_set_preflight_runtests.jl`
  - `test/nested/cartesian_pair_block_one_body_batch_summary_runtests.jl`
- simplified:
  - mixed one-body metadata scaffold carrying cost reduced by 482 deleted lines
  - active mixed one-body routine coverage remains in the consumer smoke and
    accessor/plan-batch style tests
- quarantined: nothing
- not deleted because:
  - He driver endpoint tests, WL H/H2 reference tests, and runner-included
    tests remain live
- exact remaining caller/blocker:
  - H2 PQS endpoint materialization remains blocked on
    `:missing_diatomic_complete_core_shell_source_plan_producer`; the current
    readiness artifact also records the lower source/final/H1 statuses as not
    materialized.
  - Direct comparison to old WL/QW H2 HF/ED totals remains blocked until the
    H/cc-pVTZ S/P residual supplement policy is represented in the driver/PQS
    path.

-- repo-doer@macmini

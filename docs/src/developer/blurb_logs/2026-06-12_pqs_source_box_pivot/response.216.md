Pass 216 response - H2 physical private RHF execution boundary

Result:
- Private RHF execution was not attempted.
- The physical-H2 RHF input contract remains available, but execution is blocked precisely on the missing physical execution adapter.
- `route/private_rhf_execution_status = :blocked_pqs_physical_gausslet_private_rhf_execution`
- `private_rhf/executed = false`
- `private_rhf/materialized = false`
- `private_rhf/converged = false`
- `physics/endpoint_blocker = :missing_physical_gausslet_rhf_execution_adapter`

Implemented:
- Added private `_PQSDiatomicPhysicalGaussletRHFExecutionPayload`.
- The execution payload consumes:
  - `_PQSDiatomicPhysicalGaussletRHFInputContractPayload`
  - `_PQSDiatomicPhysicalGaussletH1Payload`
  - `_PQSDiatomicPhysicalGaussletH1JPayload`
- It records that the input contract is available, but blocks before SCF because no reviewed physical-H2 execution adapter exists.
- `cartesian_assembly` now carries this execution payload.
- `cartesian_report`/artifact output now exposes compact private-RHF execution facts.
- Existing complete-core/shell private RHF helpers were left untouched; they remain label-gated to the old complete-core/shell payload family. No old object kinds were faked.

Execution/readiness facts:
- execution attempted: `false`
- execution status: `:blocked_pqs_physical_gausslet_private_rhf_execution`
- execution blocker: `:missing_physical_gausslet_rhf_execution_adapter`
- input contract available: `true`
- electron count source: `:explicit_private_rhf_electron_count`
- electron count: `2`
- occupation policy: `:closed_shell_rhf`
- occupation `nocc`: `1`
- materialized/converged: `false` / `false`
- total energy: `nothing`
- iteration count: `0`
- final density / one-step consistency status:
  `:not_evaluated_missing_physical_gausslet_rhf_execution_adapter`

Artifact fields changed:
- `route/private_rhf_execution_status`
- `private_rhf/executed`
- `private_rhf/execution_status`
- `private_rhf/execution_blocker`
- `private_rhf/converged`
- `private_rhf/total_energy`
- `private_rhf/one_body_energy`
- `private_rhf/two_body_energy`
- `private_rhf/iteration_count`
- `private_rhf/density_trace`
- `private_rhf/idempotency_residual`
- `private_rhf/commutator_residual`
- `private_rhf/energy_delta`
- `private_rhf/final_density_one_step_consistency_status`
- `physics/endpoint_blocker`

Files changed:
- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `src/pqs_source_box_route_driver_helpers.jl`
- `src/pqs_source_box_route_driver_reporting.jl`
- `test/driver_inputs/h2_pqs_q5_physical_gausslet_r4.jl`
- `test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl`
- `test/nested/integration_runtests.jl`
- `docs/src/developer/cartesian_driver_endpoint_manifest.md`
- deleted `test/nested/pqs_explicit_core_spacing_parent_axis_probe_runtests.jl`
- deleted `test/nested/pqs_raw_product_box_plan_probe_runtests.jl`

Source/test/bin scoped line budget:
- added: `183`
- deleted: `310`
- net: `-127`

Validation:
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed; precompile/load output ended with `load ok`
  - elapsed reported by Julia precompile output: `56876.3 ms`
- `julia --project=. -e 'using Test; t = @elapsed include("test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl"); println("elapsed_s=", t)'`
  - passed, `78/78`
  - `elapsed_s=80.925129167`
  - Assembly reported `31.412742` seconds with about `95%` compilation share.
- `git diff --check`
  - passed
- `git diff --cached --check`
  - passed
- `git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`
  - confirms net `-127` in scoped files.

Deletion/shrinkage report:
- deleted:
  - `test/nested/pqs_explicit_core_spacing_parent_axis_probe_runtests.jl`
  - `test/nested/pqs_raw_product_box_plan_probe_runtests.jl`
- simplified:
  - `test/nested/integration_runtests.jl` no longer includes those old probe scaffolds.
  - The H2 physical endpoint now covers the parent-axis probe path used by the active target and asserts the private-RHF execution boundary.
- quarantined:
  - RHF remains blocked at private physical execution adapter.
  - No SCF execution, DIIS tuning, HFDMRG/DMRG/CR2, GTO/MWG supplement, public API, export, HamV6 artifact, or reference comparison was added.
- not deleted because:
  - Active He and H2 endpoint tests were preserved.
  - `test/nested/cartesian_ham_builder_one_center_config_smoke_runtests.jl` was left intact because it remains the one-center WL driver smoke, distinct from the He PQS endpoint.
- exact remaining caller/blocker:
  - `physics/endpoint_blocker = :missing_physical_gausslet_rhf_execution_adapter`

Git status at response time:
- dirty worktree with the pass-216 source/test/doc edits above.
- no commits made.

-- repo-doer@macmini

Pass 215 response - H2 physical private RHF input contract

Result:
- Physical-H2 private RHF input contract is available.
- `route/private_rhf_input_contract_status = :available_pqs_physical_gausslet_rhf_input_contract`
- `private_rhf/input_contract_available = true`
- `private_rhf/materialized = false`
- `physics/endpoint_blocker = :missing_physical_gausslet_private_rhf_execution`

Manager acceptance correction:
- The first response inferred electron count from nuclear charges. That was
  corrected before acceptance.
- The H2 driver input now supplies explicit `private_rhf_electron_count = 2`.
- The physical-H2 RHF input contract now reads only that explicit route/private
  RHF electron-count input and blocks with
  `:missing_physical_gausslet_electron_count` if it is absent.

Implemented:
- Added private `_PQSDiatomicPhysicalGaussletRHFInputContractPayload`.
- The contract consumes the physical route-owned payloads:
  - `_PQSDiatomicPhysicalGaussletCoreShellSourcePlan`
  - `_PQSDiatomicPhysicalGaussletFinalBasisPayload`
  - `_PQSDiatomicPhysicalGaussletH1Payload`
  - `_PQSDiatomicPhysicalGaussletH1JPayload`
- The contract validates:
  - route kind and fixture role;
  - explicit route/private-RHF electron count `2`;
  - closed-shell RHF occupation with `nocc = 1`;
  - final dimension `463`;
  - H1 matrix availability, finiteness, and symmetry;
  - density interaction availability;
  - density gauge `:pre_final_localized_positive_weight`;
  - raw pair-factor convention `:raw_numerator`;
  - final-to-pre-final transform availability;
  - pre-final pair matrix availability, finiteness, and symmetry.
- Existing complete-core/shell private RHF helpers were left untouched. They remain label-gated to the old complete-core/shell object family; this pass adds only the physical-H2 input contract boundary and does not run SCF.

Artifact fields changed:
- `route/private_rhf_input_contract_status`
- `private_rhf/input_contract_status`
- `private_rhf/input_contract_blocker`
- `private_rhf/input_contract_available`
- `private_rhf/electron_count`
- `private_rhf/occupation_policy`
- `private_rhf/occupation_nocc`
- `private_rhf/h1_matrix_available`
- `private_rhf/h1_matrix_finite`
- `private_rhf/h1_matrix_symmetry_error`
- `private_rhf/density_interaction_available`
- `private_rhf/final_to_pre_final_transform_available`
- `private_rhf/pre_final_pair_matrix_available`
- `private_rhf/pre_final_pair_matrix_finite`
- `private_rhf/pre_final_pair_matrix_symmetry_error`
- `physics/endpoint_blocker`

Observed readiness facts:
- input contract status: `:available_pqs_physical_gausslet_rhf_input_contract`
- input contract blocker: `nothing`
- electron count: `2`
  - electron count source: explicit `private_rhf_electron_count`
- occupation policy: `:closed_shell_rhf`
- occupation `nocc`: `1`
- final dimension: `463`
- H1 matrix available/finite/symmetric: true, true, symmetry error `< 1e-8`
- density interaction available: true
- final-to-pre-final transform available: true
- pre-final pair matrix available/finite/symmetric: true, true, symmetry error `< 1e-8`
- private RHF execution/materialization: false

Files changed:
- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `src/pqs_source_box_route_driver_helpers.jl`
- `src/pqs_source_box_route_driver_reporting.jl`
- `test/driver_inputs/h2_pqs_q5_physical_gausslet_r4.jl`
- `test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl`
- `test/nested/integration_runtests.jl`
- `docs/src/developer/cartesian_driver_endpoint_manifest.md`
- deleted `test/nested/cartesian_ham_builder_diatomic_config_smoke_runtests.jl`

Source/test/bin scoped line budget:
- added: `382`
- deleted: `608`
- net: `-226`

Validation:
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed; precompile/load output ended with `load ok`
  - elapsed reported by Julia precompile output: `56808.1 ms`
- `julia --project=. -e 'using Test; t = @elapsed include("test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl"); println("elapsed_s=", t)'`
  - passed, `70/70`
  - `elapsed_s=80.608992292`
  - Runtime grew slightly from pass 214. Assembly reported `31.859427` seconds with about `95%` compilation share.
  - Manager rerun after the explicit electron-count correction also passed
    `70/70`; that run included package precompile and printed
    `elapsed_s=138.575093583`, with assembly `31.933407` seconds.
- `git diff --check`
  - passed
- `git diff --cached --check`
  - passed
- `git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`
  - confirms net `-226` in scoped files.

Deletion/shrinkage report:
- deleted:
  - `test/nested/cartesian_ham_builder_diatomic_config_smoke_runtests.jl`
- simplified:
  - `test/nested/integration_runtests.jl` no longer includes that stale slow Be2 route-configured WL driver smoke.
  - The active H2 physical endpoint test now carries the compact private-RHF input-contract assertions.
- quarantined:
  - RHF work remains an input contract only.
  - No SCF execution, DIIS tuning, HFDMRG/DMRG/CR2, GTO/MWG supplement, public API, export, or HamV6 artifact work was added.
- not deleted because:
  - The active H2 physical endpoint test remains the physics seam guard.
  - `test/nested/cartesian_ham_builder_one_center_config_smoke_runtests.jl` was left intact as the one-center driver smoke.
- exact remaining caller/blocker:
  - `physics/endpoint_blocker = :missing_physical_gausslet_private_rhf_execution`

Git status at response time:
- dirty worktree with the pass-215 source/test/doc edits above.
- no commits made.

-- repo-doer@macmini

Pass 214 response - H2 physical density/H1-J input seam

Result:
- Density/H1-J seam materialized for the H2 physical q5 gausslet-only final basis.
- `route/h1_j_status = :materialized_pqs_physical_gausslet_h1_j_payload`
- `route/h1_j_materialized = true`
- `physics/endpoint_blocker = :missing_physical_gausslet_rhf_or_solver_contract`

Implemented:
- Added private `_PQSDiatomicPhysicalGaussletH1JPayload`.
- The payload consumes:
  - `_PQSDiatomicPhysicalGaussletCoreShellSourcePlan`
  - `_PQSDiatomicPhysicalGaussletFinalBasisPayload`
  - `_PQSDiatomicPhysicalGaussletH1Payload`
- Added physical-labeled private helpers for:
  - density provenance from physical axis bundles;
  - support weights in `(:atom_contact_core, :shared_shell_1, :shared_shell_2)` order;
  - support raw pair numerator matrix;
  - physical pre-final density interaction in the localized positive-weight gauge;
  - compact H1-J/self-Coulomb diagnostic.
- `cartesian_assembly` now builds and carries the physical H1-J payload after the physical H1 payload.
- `cartesian_report`/artifact output now exposes compact density/H1-J facts.
- Turned on `run_h1_j = true` for `test/driver_inputs/h2_pqs_q5_physical_gausslet_r4.jl`.
- Updated the H2 physical driver endpoint test and endpoint manifest for the new remaining blocker.

Observed artifact facts:
- final dimension: `463`
- H1 lowest energy: `-0.7946609179724647`
- H1-J/self-Coulomb scalar: `0.45696639804337114`
- density interaction status:
  `:materialized_pqs_physical_gausslet_pre_final_density_interaction`
- density gauge: `:pre_final_localized_positive_weight`
- raw pair-factor convention: `:raw_numerator`
- support weight count: `1215`
- support weights all positive: `true`
- support raw pair shape: `(1215, 1215)`
- pre-final pair matrix shape: `(463, 463)`
- pre-final pair matrix finite: `true`
- pre-final pair matrix symmetry error: `8.881784197001252e-16`

Artifact fields changed:
- `route/h1_j_status`
- `route/h1_j_materialized`
- `density_interaction/status`
- `density_interaction/density_gauge`
- `density_interaction/raw_pair_factor_convention`
- `density_interaction/pre_final_pair_matrix_shape`
- `density_interaction/pre_final_pair_matrix_finite`
- `density_interaction/pre_final_pair_matrix_symmetry_error`
- `density_interaction/support_weight_count`
- `density_interaction/support_weights_all_positive`
- `density_interaction/support_raw_pair_shape`
- `density_interaction/support_raw_pair_finite`
- `density_interaction/h1_j_self_coulomb`
- `physics/endpoint_blocker`

Files changed:
- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `src/pqs_source_box_route_driver_helpers.jl`
- `src/pqs_source_box_route_driver_reporting.jl`
- `test/driver_inputs/h2_pqs_q5_physical_gausslet_r4.jl`
- `test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl`
- `test/nested/integration_runtests.jl`
- `docs/src/developer/cartesian_driver_endpoint_manifest.md`
- deleted `test/nested/cartesian_route_diatomic_materializer_probe_runtests.jl`

Source/test/bin scoped line budget:
- added: `671`
- deleted: `673`
- net: `-2`

Validation:
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed; precompile/load output ended with `load ok`
  - elapsed reported by Julia precompile output: `56915.3 ms`
- `julia --project=. -e 'using Test; t = @elapsed include("test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl"); println("elapsed_s=", t)'`
  - passed, `55/55`
  - `elapsed_s=78.957880083`
  - Runtime grew only slightly from the pass-213 ~78 second focused run; assembly still reported high compilation share.
- focused artifact scalar probe:
  - `artifact_probe_elapsed_s=77.535774625`
  - scalar/status lines listed in observed artifact facts above.
- `git diff --check`
  - passed
- `git diff --cached --check`
  - passed
- `git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`
  - confirms net `-2` in scoped files.

Deletion/shrinkage report:
- deleted:
  - `test/nested/cartesian_route_diatomic_materializer_probe_runtests.jl`
- simplified:
  - `test/nested/integration_runtests.jl` no longer includes that stale route-configured materializer probe.
  - The H2 physical endpoint test now checks compact artifact facts instead of requiring a separate development probe.
- quarantined:
  - Density/H1-J remains private/internal on the H2 physical route.
  - No RHF/HF/DMRG, GTO/MWG supplement, CR2/export, final dense Ham contract, or public API work was added.
- not deleted because:
  - The active H2 physical endpoint test remains the physics seam guard.
  - The physical source-plan, final-basis, H1, and H1-J payloads are needed for the next solver-contract decision.
- exact remaining caller/blocker:
  - `physics/endpoint_blocker = :missing_physical_gausslet_rhf_or_solver_contract`

Git status at response time:
- dirty worktree with the pass-214 source/test/doc edits above.
- no commits made.

-- repo-doer@macmini

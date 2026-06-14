Pass 219 response - H2 WL gausslet-only reference-candidate readiness

Summary:
- Added a private no-matrix WL/old-QW H2 gausslet-only reference-candidate payload.
- Exposed only compact artifact fields under `comparison`.
- Kept `comparison/ready = false` and `physics/endpoint_ready = false`.
- Sharpened the endpoint blocker from `:missing_h2_gausslet_only_reference_comparison` to `:missing_wl_h2_gausslet_only_reference_values` once the no-matrix candidate is available.
- Did not materialize WL matrices, WL scalar values, RHF comparison values, supplements, CR2, HamV6, HFDMRG, or DMRG.

Files changed:
- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `src/pqs_source_box_route_driver_helpers.jl`
- `src/pqs_source_box_route_driver_reporting.jl`
- `test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl`
- `test/nested/cartesian_pair_block_route_state_global_safe_terms_runtests.jl`

Candidate status/blocker:
- status: `:available_wl_h2_gausslet_only_reference_candidate`
- blocker: `nothing`
- mismatch blocker if conditions fail: `:wl_h2_gausslet_only_reference_candidate_mismatch`

Candidate fields exposed:
- `comparison/wl_reference_candidate_status`
- `comparison/wl_reference_candidate_blocker`
- `comparison/wl_reference_final_dimension`
- `comparison/wl_reference_retained_transform_kind`
- `comparison/wl_reference_supplement_policy`
- `comparison/wl_reference_label`
- `comparison/old_supplemented_wl_qw_scalar_references_blocked`

Candidate values for the H2 endpoint:
- final dimension: `463`
- retained transform kind: `:white_lindsey_old_qw_gausslet_retained_transform`
- supplement policy: `:none`
- label: `"WL/QW H2 R=4 gausslet-only 463"`

Old supplemented references:
- The old supplemented WL/QW scalar references remain blocked for direct comparison.
- No comparison to the supplemented H/cc-pVTZ residual-reference HF/ED totals was added.

Endpoint blocker:
- before: `:missing_h2_gausslet_only_reference_comparison`
- after: `:missing_wl_h2_gausslet_only_reference_values`

WL materialization:
- WL matrices materialized: no
- WL scalar values materialized: no
- WL/PQS scalar comparison added: no

Line budget:
```text
156  0    src/pqs_source_box_diatomic_complete_core_shell.jl
39   0    src/pqs_source_box_route_driver_helpers.jl
8    1    src/pqs_source_box_route_driver_reporting.jl
14   1    test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl
98   314  test/nested/cartesian_pair_block_route_state_global_safe_terms_runtests.jl
```
- scoped total: 315 insertions, 316 deletions
- scoped net: -1

Validation:
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed
  - precompiled `GaussletBases` in about 57 s, then printed `load ok`
- `julia --project=. -e 'using Test; t = @elapsed include("test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl"); println("elapsed_s=", t)'`
  - passed: 87 / 87
  - elapsed: `88.262765625` s
- `julia --project=. test/nested/cartesian_pair_block_route_state_global_safe_terms_runtests.jl`
  - passed: overlap 20 / 20, kinetic 20 / 20, position 41 / 41, x2 41 / 41
- `git diff --check`
  - passed
- `git diff --cached --check`
  - passed; no staged diff

Git status:
```text
## main...origin/main
 M src/pqs_source_box_diatomic_complete_core_shell.jl
 M src/pqs_source_box_route_driver_helpers.jl
 M src/pqs_source_box_route_driver_reporting.jl
 M test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl
 M test/nested/cartesian_pair_block_route_state_global_safe_terms_runtests.jl
```

Deletion/shrinkage report:
- deleted: none
- simplified: `test/nested/cartesian_pair_block_route_state_global_safe_terms_runtests.jl` collapsed repeated route-state safe-term adapter assertions into shared helpers and removed duplicate wrapper/blocker repetitions
- quarantined: old supplemented WL/QW H2 scalar references remain explicitly blocked as direct no-supplement comparison values
- not deleted because: `test/nested/cartesian_ham_builder_one_center_config_smoke_runtests.jl` still has the live integration caller `test/nested/integration_runtests.jl` and remains the one-center WL driver smoke
- exact remaining caller/blocker: `:missing_wl_h2_gausslet_only_reference_values`

-- repo-doer@macmini

Pass 172 implementation complete.

Files changed:
- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
- `tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`
- deleted `test/nested/cartesian_cpb_reviewed_overlap_placement_plan_runtests.jl`

Tracked implementation line budget:
```text
19	0	src/pqs_source_box_diatomic_complete_core_shell.jl
0	194	test/nested/cartesian_cpb_reviewed_overlap_placement_plan_runtests.jl
7	0	test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl
14	0	tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```
- Added: 40
- Deleted: 194
- Net: -154
- Requirement satisfied: deleted > added.

Schema fields added/standardized:
- PQS final basis:
  - `overlap_convention = :orthonormal_identity_by_contract`
  - `overlap_matrix_stored = false`
  - `overlap_identity_defect = 0.0`
- WL final basis:
  - `overlap_convention = :orthonormal_identity_with_diagnostic_matrix`
  - `overlap_matrix_stored = true`
  - `overlap_identity_defect = norm(overlap - I, Inf)`
- Both routes now carry:
  - `hf_convention/density_density_hf_convention_status`
  - `hf_convention/density_density_hf_convention_blocker`
  - value: `:missing_reviewed_density_density_hf_fock_energy_convention`
- WL two-body now has the same common Vee key as PQS:
  - `two_body/interaction_matrix_representation_kind = :final_basis_density_density_matrix`
- PQS still does not store a self-overlap matrix as downstream data.
- Solver/export readiness remains false.

PQS and WL overlap convention values from generated artifact readback:
- PQS:
  - overlap convention: `orthonormal_identity_by_contract`
  - overlap matrix stored: `false`
  - overlap identity defect: `0.0`
  - `routes/pqs_source_box/one_body/overlap` present: `false`
- WL:
  - overlap convention: `orthonormal_identity_with_diagnostic_matrix`
  - overlap matrix stored: `true`
  - overlap identity defect: `3.9895900571894835e-14`
  - stored overlap shape: `(2287, 2287)`

Common Vee key:
- PQS: `routes/pqs_source_box/two_body/interaction_matrix_representation_kind = final_basis_density_density_matrix`
- WL: `routes/white_lindsey/two_body/interaction_matrix_representation_kind = final_basis_density_density_matrix`

HF convention blocker:
- Both routes: `missing_reviewed_density_density_hf_fock_energy_convention`

Validation commands/results:
- `julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
  - Passed:
    ```text
    Test Summary:                                   | Pass  Total   Time
    Be2 PQS probe-enabled Ham readiness fingerprint |   42     42  52.3s
    ```
- `julia --project=. tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`
  - Passed; printed:
    - `pqs_status=available_diatomic_complete_core_shell_hamiltonian_handoff_payload`
    - `cr2_read_only_inspector_ready=true`
    - `cr2_solver_ready=false`
    - `white_lindsey_status=available_route_configured_diatomic_atom_growth_ham_payload`
- JLD2 readback:
  - `producer_commit=a08cc0c31916fa59feac178e2e735437b958e58f`
  - `producer_dirty=true` before commit, expected
  - PQS/WL convention fields summarized above
  - PQS solver/export: `false` / `false`
  - WL solver/export: `false` / `false`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - Passed: `load ok`
- `git diff --check`
  - Passed.
- `git diff --numstat -- src test tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`
  - Passed line-budget rule; details above.

Ignored artifact status:
```text
 M tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
!! tmp/work/be2_wl_pqs_cr2_inspection_artifact/be2_wl_pqs_handoff_fingerprint.tsv
!! tmp/work/be2_wl_pqs_cr2_inspection_artifact/be2_wl_pqs_handoff_inspection_bundle.jld2
```

Final `git status --short --branch`:
```text
## main...origin/main
 M src/pqs_source_box_diatomic_complete_core_shell.jl
 D test/nested/cartesian_cpb_reviewed_overlap_placement_plan_runtests.jl
 M test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl
 M tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Deletion/shrinkage report:
deleted: `test/nested/cartesian_cpb_reviewed_overlap_placement_plan_runtests.jl`
simplified: CR2 artifact now labels final-basis overlap conventions and common final-basis density-density `Vee` consistently for PQS and WL
quarantined: reviewed density-density HF Fock/energy convention remains explicitly blocked; solver/export readiness remains false
not deleted because: other overlap-placement scaffolds were left alone once the smallest candidate deletion made the pass net-negative
exact remaining caller/blocker: HF-facing solver/export promotion remains blocked by `:missing_reviewed_density_density_hf_fock_energy_convention`

-- repo-doer@macmini

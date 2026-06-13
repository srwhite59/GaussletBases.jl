Pass 149 response - add private diatomic H1 payload

Commit:
- Not committed; blurb did not request a commit.

Files changed:
- `src/pqs_source_box_route_driver_helpers.jl`
- `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`

Exact payload/helper names:
- `_PQSDiatomicCompleteCoreShellH1Payload`
- `_pqs_source_box_route_driver_diatomic_complete_core_shell_h1_payload`
- Private support helpers:
  - `_pqs_source_box_route_driver_diatomic_support_kinetic_matrix`
  - `_pqs_source_box_route_driver_diatomic_support_electron_nuclear_by_center_matrices`
- New `cartesian_assembly` field:
  - `diatomic_complete_core_shell_h1_payload`

Implementation summary:
- Added a private Be2/PQS complete core/shell H1 payload after the private diatomic final-basis payload.
- Did not call `pqs_multilayer_complete_core_shell_h1_payload(...)`.
- Left old `pqs_multilayer_support_*` and old one-center H1 object-kind guards unchanged.
- Built support kinetic from the private diatomic source-plan metrics and `core_support_states` then `shell_support_states`.
- Built support electron-nuclear by-center matrices from parent center records, axis layers, Coulomb Gaussian expansion data, and existing lower Gaussian/product primitives.
- Reused lower `CartesianFinalBasisRealization` final one-body transfer, final one-electron Hamiltonian assembly, and H1 solve.

Default Be2/PQS H1 status:
- `diatomic_complete_core_shell_h1_payload.status = :blocked_diatomic_complete_core_shell_h1_payload`
- `blocker = :missing_pqs_diatomic_complete_core_shell_source_plan`
- `source_plan_status = :not_available`
- `final_basis_status = :not_available`
- support/final/H1 statuses remain not materialized.
- H1-J, density-density, Ham, RHF, public/export/artifact flags remain false.

Probe-enabled Be2/PQS H1 status:
- `diatomic_complete_core_shell_h1_payload.status = :available_diatomic_complete_core_shell_h1_payload`
- `blocker = nothing`
- `source_plan_status = :available_pqs_diatomic_complete_core_shell_source_plan`
- `final_basis_status = :available_pqs_complete_core_shell_final_basis`
- `support_kinetic_status = :materialized_diatomic_complete_core_shell_support_kinetic_matrix`
- `support_electron_nuclear_status = :materialized_diatomic_complete_core_shell_support_electron_nuclear_by_center_matrix_set`
- `final_kinetic_status = :materialized_pqs_complete_core_shell_final_one_body_matrix`
- `final_electron_nuclear_status = :materialized_diatomic_final_electron_nuclear_by_center`
- `final_hamiltonian_status = :materialized_pqs_complete_core_shell_final_one_electron_hamiltonian`
- `h1_status = :materialized_pqs_complete_core_shell_final_h1_solve`

Observed H1 facts:
- `final_dimension = 221`
- `lowest_energy = -0.27746109235228694`
- `center_count = 2`
- final-basis support row order remains `:core_then_shell`
- source-plan support row order remains `:core_product_then_shell_left_right_pqs`

Readiness/Ham fingerprint:
- Probe-enabled readiness now reports:
  - `:diatomic_complete_core_shell_h1_payload` in `available_objects`
  - `blocker = :missing_diatomic_complete_core_shell_h1_j_consumer`
  - `:diatomic_complete_core_shell_h1_j_consumer` in `missing_objects`
- Probe-enabled readiness no longer reports `:diatomic_complete_core_shell_h1_consumer` as missing once private H1 is available.
- The older generic complete-core/shell Ham payload remains blocked; this pass did not claim a full Ham payload.

Boundary confirmation:
- No H1-J, density-density, electron-electron Ham payload, full Ham payload, RHF/SCF/DIIS, WL payload, public API, export, artifact, hfdmrg, or CR2 behavior was added.
- No shell/support-row contraction, raw product-box probe, or old WL adapter path was promoted to route authority.
- No retained diagnostic/self-integral weights were reinterpreted as IDA/quadrature weights.

Validation:
- `git diff --check`
  - passed.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed.
  - precompilation: `56682.1 ms  ✓ GaussletBases`
  - printed `load ok`.
- `julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
  - passed.
  - `Be2 PQS Ham payload readiness fingerprint | Pass 182 Total 182 Time 13.3s`
  - `Be2 PQS probe-enabled Ham readiness fingerprint | Pass 259 Total 259 Time 40.1s`
  - focused test body time was about 53.4s, below 60s after precompilation.
- Focused energy probe:
  - `h1_status=available_diatomic_complete_core_shell_h1_payload`
  - `final_dimension=221`
  - `lowest_energy=-0.27746109235228694`
  - `readiness_blocker=missing_diatomic_complete_core_shell_h1_j_consumer`

Git status:
```text
## main...origin/main
 M src/pqs_source_box_route_driver_helpers.jl
 M test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl
```

Deletion/shrinkage report:
- deleted: none; this pass added the private H1 seam.
- simplified: probe-enabled Be2/PQS readiness now advances past the missing private H1 consumer and reports the sharper H1-J consumer blocker.
- quarantined: old `:pqs_multilayer_shell_source_plan` claim; old one-center H1 helper guard; H1-J/density-density/electron-electron Ham payload; full Ham payload; RHF/SCF/DIIS; WL payloads; public API/exports/artifacts/hfdmrg/CR2.
- not deleted because: old one-center H1 helpers remain the validated compact one-center fixture path, and no H1-J/electron-electron consumer has been approved yet.
- exact remaining caller/blocker: `cartesian_assembly(...).diatomic_complete_core_shell_ham_readiness_payload` now blocks on `:missing_diatomic_complete_core_shell_h1_j_consumer`.

-- repo-doer@macmini

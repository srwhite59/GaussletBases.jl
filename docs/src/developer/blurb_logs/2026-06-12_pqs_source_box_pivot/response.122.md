Pass 122 response - bounded private RHF Fock-DIIS controls

Files changed:
- `src/pqs_multilayer_complete_core_shell_rhf.jl`
- `test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl`

Helper/control names:
- Added private control helper:
  `_pqs_multilayer_complete_core_shell_rhf_scf_control_payload(...)`
- Added private loop helpers:
  - `_pqs_multilayer_complete_core_shell_rhf_density_from_fock(...)`
  - `_pqs_multilayer_complete_core_shell_rhf_diis_error_vector(...)`
  - `_pqs_multilayer_complete_core_shell_rhf_diis_coefficients(...)`
  - `_pqs_multilayer_complete_core_shell_rhf_mixed_fock(...)`
- Updated `_pqs_multilayer_complete_core_shell_rhf_scf_payload(...)` to accept
  either `scf_control_payload` or existing keyword controls.

Control payload fields:
- `object_kind`
- `status`
- `blocker`
- `mixing_kind`
- `max_iterations`
- `density_atol`
- `energy_atol`
- `residual_atol`
- `trace_atol`
- `idempotency_atol`
- `max_history`
- `diis_start_iteration`
- `diis_regularization`
- `diis_coefficient_max_abs`
- residual/metric labels
- `summary`
- `metadata`

Status/blocker labels:
- Control payload available:
  `:available_pqs_multilayer_complete_core_shell_rhf_scf_control_payload`
- Control payload blocked:
  `:blocked_pqs_multilayer_complete_core_shell_rhf_scf_control_payload`
- Unsupported mixing:
  `:unsupported_scf_mixing_kind`
- Invalid controls:
  `:invalid_scf_controls`
- SCF max-iteration failure remains:
  `:blocked_pqs_multilayer_complete_core_shell_rhf_scf_payload` with
  `:scf_not_converged`

Implemented control behavior:
- Supported `mixing_kind`:
  - `:fixed_point`
  - `:fock_diis`
- Default behavior remains fixed-point for existing callers.
- `:fock_diis` defaults:
  - `max_history = 6`
  - `diis_start_iteration = 2`
  - `diis_regularization = 1.0e-12`
  - `diis_coefficient_max_abs = 25.0`
- `:scalar_density_damping` is not implemented and blocks as unsupported.
- DIIS uses symmetrized Fock matrices and commutator error vectors
  `vec(F * D - D * F)` with `D = P / occupancy`.
- The selected Fock is diagonalized to rebuild an idempotent closed-shell
  spin-summed density.
- Convergence now requires density, energy, commutator residual, trace, and
  idempotency checks. The first-iteration missing-energy rule remains explicit.

DIIS fallback/pathology fields:
- Iteration records now include:
  - `mixing_kind`
  - `diis_used`
  - `diis_fallback`
  - `diis_fallback_reason`
  - `commutator_residual`
  - `trace_error`
  - `idempotency_error`
  - `update_density_change`
  - convergence booleans for density, energy, residual, trace, and idempotency
- SCF summaries include compact counters:
  - `diis_used_count`
  - `diis_fallback_count`
  - `diis_solve_failure_count`
  - `diis_coefficient_pathology_count`

Tests:
- Updated the focused synthetic SCF test to cover:
  - fixed-point control payload defaults;
  - Fock-DIIS control payload defaults;
  - unsupported mixing kind blocker;
  - preserved fixed-point zero-interaction convergence;
  - Fock-DIIS on the tiny self-consistent fixture;
  - trace/idempotency/residual labels and nonclaim flags.
- No tracked compact route-smoke test was added.

Validation:
- `julia --project=. test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl`
  - Passed: `88 / 88`
  - Note: run precompiled `GaussletBases`.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - Passed: `load ok`
- `git diff --check`
  - Passed.

Git status:
```text
## main...origin/main
 M src/pqs_multilayer_complete_core_shell_rhf.jl
 M test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl
```

Deletion/shrinkage report:
- deleted: none.
- simplified: loose SCF control semantics are now carried by one compact
  private control payload.
- quarantined: no new tracked quarantine; ignored local DIIS/residual probes
  remain under `tmp/work`.
- not deleted because: the real compact route behavior still needs local
  validation through the new tracked private Fock-DIIS path before any cleanup
  of historical probes or route-facing decisions.
- exact remaining caller/blocker:
  `_pqs_multilayer_complete_core_shell_rhf_scf_payload(...)` now has bounded
  private Fock-DIIS controls, but no compact route rerun has validated the
  tracked implementation yet; route-driver integration remains intentionally
  blocked.

-- repo-doer@macmini

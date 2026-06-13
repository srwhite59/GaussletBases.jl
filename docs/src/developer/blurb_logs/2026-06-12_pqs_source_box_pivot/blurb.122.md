Pass 122 - implement bounded private RHF Fock-DIIS controls

Baseline:

- Current pushed HEAD should include `277feae0 Record PQS RHF DIIS experiment`.
- Local pass 121 showed bounded Fock DIIS can reduce the compact route-smoke
  commutator residual, but no tracked implementation exists yet.

Task:

Implement a bounded private SCF-control payload and Fock-DIIS path inside the
private RHF diagnostic helper. Keep it private; no route/report/public behavior.

Implementation surface:

- `src/pqs_multilayer_complete_core_shell_rhf.jl`
- `test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl`

Required pieces:

1. Private control payload/helper

- Add a compact helper such as
  `_pqs_multilayer_complete_core_shell_rhf_scf_control_payload(...)`.
- Supported `mixing_kind` values:
  - `:fixed_point`
  - `:fock_diis`
- Optional `:scalar_density_damping` may be recognized as blocked/unsupported
  unless it is trivial to keep honest; do not implement density damping unless
  needed for compatibility.
- Include compact fields:
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
- Defaults:
  - fixed-point behavior should preserve existing callers.
  - for `:fock_diis`: `max_history = 6`,
    `diis_start_iteration = 2`, `diis_regularization = 1.0e-12`,
    `diis_coefficient_max_abs = 25.0`.
- Invalid controls should block with precise labels such as
  `:invalid_scf_controls` or `:unsupported_scf_mixing_kind`.

2. SCF loop integration

- Update `_pqs_multilayer_complete_core_shell_rhf_scf_payload(...)` to consume
  either an explicit control payload or existing keyword controls.
- Preserve existing default fixed-point behavior for current tests/callers.
- For `:fock_diis`:
  - store bounded history of symmetrized Fock matrices and commutator error
    vectors;
  - solve the constrained DIIS system with regularization;
  - if solve fails or coefficients exceed `diis_coefficient_max_abs`, fall back
    for that iteration to the current Fock and increment compact fallback
    counters;
  - diagonalize the selected/mixed Fock and rebuild idempotent closed-shell
    density from occupied orbitals.
- Keep iteration records compact. Add fields only if useful, such as:
  - `mixing_kind`
  - `diis_used`
  - `diis_fallback`
  - `commutator_residual`
  - `trace_error`
  - `idempotency_error`
- Convergence should require density, energy, residual, trace, and idempotency
  checks. Keep the first-iteration energy rule explicit.

3. Tests

- Update the existing focused SCF test.
- Preserve current fixed-point zero-interaction convergence behavior.
- Add synthetic control tests for:
  - control payload accepts `:fixed_point`;
  - control payload accepts `:fock_diis` with default history/regularization;
  - unsupported mixing kind blocks;
  - `:fock_diis` on the tiny self-consistent fixture still converges and keeps
    trace/idempotency/residual labels valid.
- Do not add a tracked compact route-smoke test.

Exclusions:

- Do not run the compact route probe in this pass.
- Do not add route-driver wiring.
- Do not add report aliases/options/public API.
- Do not touch GTO, IDA/MWG, exports, artifacts, or production route behavior.
- Do not promote fixtures to physics acceptance.
- Do not request interactive escalation/approval; if approval is truly needed,
  write `ATTENTION.md` and stop.

Validation:

- `julia --project=. test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Report back:

- Helper/control names.
- Status/blocker labels.
- DIIS fallback/pathology fields.
- Validation commands/results.
- Git status.
- Deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini

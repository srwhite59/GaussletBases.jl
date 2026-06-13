Pass 120 response - private RHF SCF-control/acceleration contract design

No code was changed.

Recommended first implementation pass:
- First do a local-only Fock-DIIS experiment under `tmp/work`, using the
  residual-aware private SCF payload and the compact route-smoke fixture.
- Do not add production/private DIIS code until that experiment shows the
  commutator residual moves down from the current `~1.3e-5` scale.
- Reason: scalar damping already failed to converge within 100 iterations, and
  density mixing reduced the reported step size without reducing the underlying
  candidate density/subspace or commutator residual. Density mixing also makes
  the mixed density slightly non-idempotent. Fock DIIS keeps the accelerated
  object on the Fock side; diagonalizing the mixed Fock then rebuilds an
  idempotent closed-shell density with the correct trace.

Control object:
- Add a private compact controls payload first, not public API:
  `_pqs_multilayer_complete_core_shell_rhf_scf_control_payload(...)`.
- Keep it in `src/pqs_multilayer_complete_core_shell_rhf.jl` near the existing
  SCF controls validation.
- Proposed fields:
  - `object_kind = :pqs_multilayer_complete_core_shell_rhf_scf_control_payload`
  - `status`
  - `blocker`
  - `mixing_kind`
  - `max_iterations`
  - `density_atol`
  - `energy_atol`
  - `residual_atol`
  - `trace_atol`
  - `idempotency_atol`
  - `mixing_alpha`
  - `max_history`
  - `diis_start_iteration`
  - `diis_regularization`
  - `density_change_rule`
  - `residual_metric`
  - `idempotency_rule`
  - `orbital_metric`
  - `summary`
  - `metadata`
- Required caller-visible configuration after defaults:
  - `mixing_kind`
  - `max_iterations`
  - `density_atol`
  - `energy_atol`
  - `residual_atol`
- Optional/defaulted fields:
  - `trace_atol = density_atol`
  - `idempotency_atol = density_atol`
  - `mixing_alpha = 1.0`
  - `max_history = 6` for DIIS, `0` or `nothing` for fixed-point/scalar damping
  - `diis_start_iteration = 2`
  - `diis_regularization = 1.0e-12`
  - labels default to the existing residual-aware names:
    - `density_change_rule = :fixed_point_spin_summed_density_inf_norm`
    - `residual_metric = :ordinary_final_basis_commutator_inf_norm`
    - `idempotency_rule = :closed_shell_spatial_density_idempotency`
    - `orbital_metric = :ordinary_orthonormal_final_basis`
- Supported `mixing_kind` values should start small:
  - `:fixed_point`
  - `:scalar_density_damping`
  - `:fock_diis`
- Unsupported values should block with `:unsupported_scf_mixing_kind`.

First acceleration strategy:
- Choose local-only Fock DIIS first.
- Do not make residual-aware scalar damping the first implementation; the pass
  115/117 evidence says scalar damping changes the displayed density step but
  does not materially reduce the residual.
- Do not choose density DIIS first; it mixes density matrices directly and can
  lose idempotency unless followed by projection/diagonalization rules that
  become a second algorithm.
- Fock DIIS is the smallest next experiment because the existing SCF loop
  already builds a symmetric ordinary final-basis Fock and then diagonalizes it
  to rebuild `P = occupancy * C_occ*C_occ'`.

Residual/history definition:
- Use the ordinary final-basis commutator as the DIIS error vector:
  `vec(F * D - D * F)` where `D = P / occupancy`.
- Keep the public-facing residual label aligned with the existing payload:
  `:ordinary_final_basis_commutator_inf_norm`.
- Store histories locally/privately:
  - symmetrized Fock matrices for extrapolation;
  - flattened commutator vectors for DIIS errors;
  - compact per-iteration summaries only in returned payloads.
- History size:
  - default `max_history = 6`;
  - minimum useful DIIS history `2`;
  - for the compact route-smoke fixture, test `4` and `6` locally before
    increasing.
- Preserve trace/idempotency by mixing Fock matrices, diagonalizing the mixed
  Fock, and rebuilding the density from the occupied orbitals. Do not accept a
  mixed density as final without idempotency and trace checks.

Convergence labels:
- `rhf_converged = true` should require all of:
  - fixed-point spin-summed density delta `<= density_atol`;
  - energy delta `<= energy_atol`, with the current first-iteration missing
    energy rule kept only for cases where density/residual are already below
    tolerance;
  - ordinary final-basis commutator residual `<= residual_atol`;
  - density trace error `<= trace_atol`;
  - closed-shell spatial-density idempotency error `<= idempotency_atol`.
- A compact `convergence_diagnostics` or extended residual summary can record
  booleans:
  - `density_converged`
  - `energy_converged`
  - `residual_converged`
  - `trace_converged`
  - `idempotency_converged`
- Max-iteration failure should keep:
  - status `:blocked_pqs_multilayer_complete_core_shell_rhf_scf_payload`
  - blocker `:scf_not_converged`
- Add a compact blocker tuple only if useful, e.g.
  `convergence_blockers = (:density_delta_above_atol, :residual_above_atol)`.

Implementation boundary:
- First tracked implementation surface, after local Fock-DIIS experiment:
  `src/pqs_multilayer_complete_core_shell_rhf.jl`.
- Touch first:
  - `_pqs_multilayer_complete_core_shell_rhf_scf_controls_blocker(...)`
    or replace it with the compact control payload helper;
  - `_pqs_multilayer_complete_core_shell_rhf_scf_payload(...)` loop internals.
- Keep route-driver integration blocked.
- Keep default tracked tests synthetic only:
  - update `test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl`
    for control payload validation and zero-interaction convergence labels;
  - add a tiny synthetic DIIS/control case only if it remains small and does
    not encode route-smoke fixture details.
- Real compact route validation should remain ignored/local until the private
  SCF control is stable.

Explicit exclusions:
- No route-driver wiring.
- No public API, report aliases, exports, artifacts, GTO, IDA/MWG, or fixture
  promotion.
- No RHF production claim.
- No final-basis generalized-overlap SCF path; ordinary final-basis metric
  remains the assumption.
- No shell/support-row contraction promotion.

Git status:
```text
## main...origin/main
```

Deletion/shrinkage forecast:
- deleted: none.
- simplified: future controls should replace loose `max_iterations`,
  `density_atol`, and `energy_atol` handling with one compact private controls
  payload.
- quarantined: real compact route acceleration experiments should remain
  ignored `tmp/work` artifacts until stable.
- not deleted because: current private SCF helper is still the active
  diagnostic seam and is the right place to test residual-controlled
  convergence before route adoption.
- exact remaining caller/blocker:
  `_pqs_multilayer_complete_core_shell_rhf_scf_payload(...)` still blocks the
  compact real fixture with `:scf_not_converged` and commutator residual around
  `1.34e-5`; the next blocker is proving a private Fock-DIIS/control contract
  can reduce that residual without changing route behavior.

-- repo-doer@macmini

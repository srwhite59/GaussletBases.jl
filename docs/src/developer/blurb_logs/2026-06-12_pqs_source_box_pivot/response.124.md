Pass 124 response - tracked Fock-DIIS compact route-smoke validation

Scope:
- Local ignored `tmp/work` probe only.
- No tracked source/test/doc edits.
- No route wiring/report fields/public API/exports/artifacts/GTO/IDA/MWG.
- No fixture promotion.

Artifact paths:
- Probe script:
  `tmp/work/pqs_compact_route_owned_rhf_tracked_fock_diis_probe.jl`
- Summary:
  `tmp/work/pqs_compact_route_owned_rhf_tracked_fock_diis_probe_summary.txt`

Command/results:
- `julia --project=. tmp/work/pqs_compact_route_owned_rhf_tracked_fock_diis_probe.jl`
  - elapsed: `107.736029333` seconds
- `git status --short --branch`
  - `## main...origin/main`

Probe summary:
```text
control_status                         available_pqs_multilayer_complete_core_shell_rhf_scf_control_payload
control_blocker                        nothing
mixing_kind                            fock_diis
max_history                            6
diis_start_iteration                   2
diis_regularization                    9.9999999999999998e-13
diis_coefficient_max_abs               25
max_iterations                         100
density_atol                           1e-08
energy_atol                            1e-10
residual_atol                          1e-08
scf_status                             blocked_pqs_multilayer_complete_core_shell_rhf_scf_payload
scf_blocker                            scf_not_converged
rhf_converged                          false
iteration_count                        100
final_or_last_total_energy             -10.032119189981447
last_iteration_total_energy            -10.032119189981199
density_change                         1.1434019359102265e-07
update_density_change                  3.2634555308108659e-08
energy_change                          2.4158453015843406e-13
residual_status                        materialized_pqs_multilayer_complete_core_shell_rhf_scf_residual_diagnostics
residual_blocker                       nothing
commutator_residual                    2.9183295714024782e-08
spatial_commutator_residual            1.4591647857012391e-08
trace_error                            4.4408920985006262e-16
idempotency_error                      5.5511151231257827e-17
diis_used_count                        99
diis_fallback_count                    0
diis_solve_failure_count               0
diis_coefficient_pathology_count       0
last_diis_used                         true
last_diis_fallback                     false
last_diis_fallback_reason              nothing
last_density_converged                 false
last_energy_converged                  true
last_residual_converged                false
last_trace_converged                   true
last_idempotency_converged             true
final_convergence_blocked_by_recomputed_diagnostics false
driver_route_materialized              false
route_report_materialized              false
exports_materialized                   false
artifacts_materialized                 false
public_api                             false
private_diagnostic_only                true
```

Diagnosis:
- The tracked private Fock-DIIS path reproduces the pass-121 local experiment
  scale.
- It does not converge under the strict `1e-8` density/residual gates within
  100 iterations.
- It does substantially reduce the residual compared with fixed-point:
  - previous fixed-point commutator residual scale: `~1.3e-5`;
  - tracked Fock-DIIS commutator residual: `2.9183295714024782e-08`.
- Trace and idempotency are at machine precision.
- DIIS itself behaved cleanly in this run:
  - `diis_used_count = 99`;
  - no solve failures;
  - no coefficient pathologies;
  - no fallbacks.
- Final convergence was not blocked by the new final recomputed diagnostic gate;
  the last iteration itself still missed density and residual tolerances.

Recommended next pass:
- Do not loosen tolerances yet.
- Run a local-only asymptote/parameter probe for the tracked private Fock-DIIS
  path:
  - same `max_history = 6`, `regularization = 1e-12`;
  - increase `max_iterations` to 200;
  - optionally compare `max_history = 8` only locally;
  - preserve strict `1e-8` gates for measurement.
- If the residual crosses tolerance with more iterations, the next tracked
  change should be only a private max-iteration/default-control adjustment.
- If it plateaus near `3e-8`, the next design question is whether this compact
  diagnostic path needs a realistic residual floor/tolerance policy before any
  route-facing work.

Git status:
```text
## main...origin/main
```

Deletion/shrinkage report:
- deleted: none.
- simplified: this validates that the tracked private Fock-DIIS path replaces
  the separate local DIIS implementation for this measurement.
- quarantined: tracked-Fock-DIIS probe script and summary remain ignored
  `tmp/work` artifacts.
- not deleted because: prior ignored residual/DIIS probes remain useful
  comparison artifacts until the convergence/tolerance policy is settled.
- exact remaining caller/blocker:
  `_pqs_multilayer_complete_core_shell_rhf_scf_payload(...)` with tracked
  `:fock_diis` remains blocked by `:scf_not_converged` at 100 iterations:
  density delta `1.14e-7`, commutator residual `2.92e-8`, strict gate `1e-8`.

-- repo-doer@macmini

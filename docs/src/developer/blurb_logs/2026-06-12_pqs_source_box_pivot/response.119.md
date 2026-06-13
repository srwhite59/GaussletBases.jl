Pass 119 response - local compact PQS probe with residual-aware SCF payload

Scope:
- Local ignored `tmp/work` probe only.
- No tracked source/test/doc edits.
- No damping, DIIS, acceleration, route wiring, report fields, public API,
  exports, artifacts, GTO, IDA/MWG, or fixture promotion.

Artifact paths:
- Probe script:
  `tmp/work/pqs_compact_route_owned_rhf_residual_payload_probe.jl`
- Summary:
  `tmp/work/pqs_compact_route_owned_rhf_residual_payload_probe_summary.txt`

Command/results:
- `julia --project=. tmp/work/pqs_compact_route_owned_rhf_residual_payload_probe.jl`
  - elapsed: `99.980484375000003` seconds
- `git status --short --branch`
  - `## main...origin/main`

Probe results:
```text
electron_count                         4
fixture_role                           route_smoke
scf_status                             blocked_pqs_multilayer_complete_core_shell_rhf_scf_payload
scf_blocker                            scf_not_converged
iteration_count                        50
last_iteration_total_energy            -10.032175867411672
final_total_energy                     -10.032176042852385
density_change                         4.5438494082283842e-05
energy_change                          1.750599221139737e-07
final_one_step_recomputed              true
final_one_step_payload_role            recomputed_final_density_diagnostic
residual_status                        materialized_pqs_multilayer_complete_core_shell_rhf_scf_residual_diagnostics
residual_blocker                       nothing
density_trace_error                    3.5527136788005009e-15
closed_shell_idempotency_error         9.0205620750793969e-17
commutator_residual                    1.3415879049821777e-05
spatial_commutator_residual            6.7079395249108886e-06
density_change_rule                    fixed_point_spin_summed_density_inf_norm
residual_metric                        ordinary_final_basis_commutator_inf_norm
idempotency_rule                       closed_shell_spatial_density_idempotency
orbital_metric                         ordinary_orthonormal_final_basis
rhf_materialized                       false
rhf_converged                          false
driver_route_materialized              false
route_report_materialized              false
exports_materialized                   false
artifacts_materialized                 false
public_api                             false
private_diagnostic_only                true
```

Diagnosis:
- The updated private SCF payload exposes the expected residual diagnostics
  directly on the nonconverged real compact route-smoke fixture.
- The values match the pass-117 local residual scale:
  - commutator residual is about `1.34e-5`;
  - spatial commutator residual is about `6.71e-6`;
  - trace and closed-shell idempotency are near machine precision.
- The blocked SCF state is therefore a real residual/stationarity issue, not a
  missing residual-reporting issue.

Recommended next pass:
- Design-only private SCF-control contract for acceleration/mixing.
- Keep the contract private and diagnostic-only; specify what residual is used
  for control, how density/Fock history is bounded, and what convergence labels
  are allowed before implementing DIIS or another accelerator.
- Do not route-wire RHF or promote report/public behavior yet.

Git status:
```text
## main...origin/main
```

Deletion/shrinkage report:
- deleted: none.
- simplified: no tracked code changed; this pass confirmed the new compact
  residual payload removes dependence on the separate residual script for this
  measurement.
- quarantined: local probe script and summary remain ignored `tmp/work`
  artifacts.
- not deleted because: ignored historical probes are still useful comparison
  artifacts until the private SCF-control design is settled.
- exact remaining caller/blocker:
  `_pqs_multilayer_complete_core_shell_rhf_scf_payload(...)` now reports
  materialized residual diagnostics, but remains blocked by
  `:scf_not_converged` with `commutator_residual ≈ 1.34e-5`; the next blocker is
  private acceleration/mixing contract design.

-- repo-doer@macmini

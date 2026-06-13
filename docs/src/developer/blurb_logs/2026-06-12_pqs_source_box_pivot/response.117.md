Pass 117 response - compact RHF commutator/idempotency residual probe

Scope:
- Added/updated ignored local `tmp/work` probe only.
- No tracked source, tracked tests, tracked docs, route wiring, public API,
  exports, artifacts, damping/DIIS, or production behavior changed.

Artifacts:
- Probe script:
  `tmp/work/pqs_compact_route_owned_rhf_residual_probe.jl`
- Summary:
  `tmp/work/pqs_compact_route_owned_rhf_residual_probe_summary.txt`
- Table:
  `tmp/work/pqs_compact_route_owned_rhf_residual_probe.tsv`

Commands/results:
- `julia --project=. tmp/work/pqs_compact_route_owned_rhf_residual_probe.jl`
  - First attempt failed after fixture build due to a typo in the ignored local
    script. I fixed only that local probe typo.
  - Successful rerun elapsed: `99.661981999999995` seconds.
- `git status --short --branch`
  - `## main...origin/main`

Residual table:

```text
label               alpha  status                                      iter  total_energy        last_dE              final_dE_vs_last     fixed_point_dP      candidate_dP        trace               trace_err            idempotency_err      comm_p              comm_d              occupied_projector_delta  subspace_differs
undamped            1.0    blocked_local_pqs_complete_core_shell_rhf_scf 50 -10.032176042852385 1.750599221139737e-7 1.7544071262420857e-7 4.5470083933824013e-5 4.5470083933824013e-5 4.0000000000000036 3.5527136788005009e-15 9.0205620750793969e-17 1.3415879049821777e-5 6.7079395249108886e-6 2.2735041966912006e-5 true
damped_alpha_0p25   0.25   blocked_local_pqs_complete_core_shell_rhf_scf 50 -10.032169325493694 2.1194997401607907e-7 1.6907443978197989e-7 1.1187935222756185e-5 4.4751740890969227e-5 4.0000000000000018 1.7763568394002505e-15 1.7013735073567737e-8 1.3003822598410542e-5 6.5019112992052708e-6 2.2375870445484614e-5 true
```

Interpretation:
- Trace is fine for both states: errors are around `1e-15`.
- Undamped closed-shell idempotency is fine: about `9e-17`.
- The damped mixed density is slightly non-idempotent, about `1.7e-8`, as
  expected for scalar density mixing.
- The ordinary final-basis commutator residual is not small relative to the
  requested `1e-8` density tolerance:
  - undamped `norm(F*P - P*F, Inf) ≈ 1.34e-5`;
  - damped alpha 0.25 `norm(F*P - P*F, Inf) ≈ 1.30e-5`.
- The occupied projector/subspace delta is also material:
  - undamped about `2.27e-5`;
  - damped alpha 0.25 about `2.24e-5`.
- Damping alpha 0.25 reduces the damped fixed-point step to about `1.12e-5`,
  but the underlying candidate density change remains about `4.48e-5`.

Qualitative diagnosis:
- This is not just a bookkeeping failure of the current fixed-point density
  criterion. The final densities are not stationary with respect to the
  recomputed ordinary-final-basis Fock operator.
- The evidence still does not point to a factor-of-two Fock/energy convention
  bug: trace and idempotency are controlled, and pass 116 verified the
  Fock/energy derivative convention.
- The nonconvergence is therefore best classified as a real private diagnostic
  SCF-control problem. Scalar damping can reduce the reported step size, but it
  does not remove the commutator/subspace residual at this fixture scale.

Recommended next pass:
- Do a small controlled SCF-control design/implementation pass for the private
  diagnostic path only:
  - carry explicit residual labels in the private SCF summary, especially
    `density_change_rule` and `residual_metric`;
  - add commutator/idempotency residual reporting to the private diagnostic
    payload or local probe path before route wiring;
  - then test one acceleration/mixing strategy beyond scalar damping, such as a
    compact DIIS-style density/Fock residual experiment.
- Do not wire the route driver or promote this to production behavior until the
  private residual-controlled SCF behavior is understood.

Bookkeeping/formula bugs:
- Production formula bug found: none.
- Residual computation bookkeeping bug found: none in production. The only bug
  encountered was a typo in the ignored local probe script, fixed locally.

Git status:
```text
## main...origin/main
```

Deletion/shrinkage report:
- deleted: none.
- simplified: none in this local probe pass.
- quarantined: residual probe, table, and summary remain ignored `tmp/work`
  artifacts.
- not deleted because: the private RHF one-step/SCF helpers remain the active
  diagnostic seam and now need residual-aware SCF control before shrinkage.
- exact remaining caller/blocker:
  `_pqs_multilayer_complete_core_shell_rhf_scf_payload(...)` still blocks with
  `:scf_not_converged`; measured ordinary-final-basis commutator residual is
  about `1.3e-5`, so the next blocker is private SCF control/residual handling,
  not route wiring.

-- repo-doer@macmini

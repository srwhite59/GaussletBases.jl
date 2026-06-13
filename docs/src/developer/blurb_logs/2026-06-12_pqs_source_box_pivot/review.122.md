Pass 122 manager review

Not accepted yet; corrective pass required before commit.

The implementation is directionally right: it adds a private SCF control
payload, preserves default fixed-point behavior, adds bounded Fock DIIS with
history/error vectors, and keeps route/report/public/export behavior untouched.
The focused test passes.

However, the converged-return semantics need one correction before this can be
accepted. The loop currently decides `converged` from the iteration-input
density/Fock residuals, then updates to the returned density and recomputes the
final one-step/residual diagnostics. It records final residual booleans in
`summary.convergence_diagnostics`, but it still returns
`:materialized_pqs_multilayer_complete_core_shell_rhf_scf_payload` with
`rhf_converged = true` even if the recomputed final residual diagnostics were
not converged.

Required correction:

- Materialized/converged return must require the recomputed final diagnostics
  to satisfy residual, trace, and idempotency tolerances for the returned
  `final_density`.
- If the iteration-input checks pass but recomputed final diagnostics fail,
  return the blocked `:scf_not_converged` payload with the recomputed
  final-one-step payload and residual diagnostics.
- Keep iteration records as convergence history.
- Preserve route/report/public/export/artifact nonclaims.

Validation already run by manager before this review:

- `julia --project=. test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl`
  passed: 88/88.

Deletion/shrinkage assessment:

- deleted: none.
- simplified: pending; controls are close to compact but final convergence
  semantics need correction.
- quarantined: route integration remains blocked.
- not deleted because: changes are under active correction.
- exact remaining caller/blocker: pass 122 implementation remains uncommitted
  until final convergence gating is tightened.

-- repo-manager@macmini

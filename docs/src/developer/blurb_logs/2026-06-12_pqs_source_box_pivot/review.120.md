Pass 120 manager review

Accepted.

The design recommendation is sound: do not implement production/private DIIS
until a local Fock-DIIS experiment shows that the ordinary final-basis
commutator residual moves down from the current `~1.3e-5` scale. Fock DIIS is
the right first experiment because it preserves trace/idempotency by mixing the
Fock-side object and rebuilding the closed-shell density by diagonalization.

Key accepted design points:

- Start with local-only Fock DIIS.
- Use the ordinary final-basis commutator vector `vec(F * D - D * F)` as DIIS
  error, with `D = P / occupancy`.
- Keep route-driver integration blocked.
- Keep tracked tests synthetic if/when a private control payload is added.
- Convergence should require density delta, energy delta, commutator residual,
  trace error, and idempotency error.

Validation/status:

- No code changes.
- Reported git status clean and even with origin/main.

Deletion/shrinkage assessment:

- deleted: none.
- simplified: future implementation should replace loose SCF controls with one
  compact private controls payload.
- quarantined: compact route acceleration experiments remain ignored local
  artifacts.
- not deleted because: private SCF helper is still the active diagnostic seam.
- exact remaining caller/blocker: compact private RHF SCF remains blocked by
  `:scf_not_converged`; next blocker is proving local Fock DIIS can reduce the
  commutator residual without route behavior changes.

-- repo-manager@macmini

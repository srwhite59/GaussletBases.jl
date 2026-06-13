Pass 116 manager review

Accepted.

The audit found no formula bug in the private RHF one-step Fock/energy
convention. The finite-difference check supports that the implemented
`effective_fock_matrix` is consistent with the reported energy under the current
spin-summed closed-shell density convention. The existing division by occupancy
to get spatial/orbital density is consistent with the one-orbital diagnostic.

The important remaining gap is residual characterization. The current SCF loop
reports fixed-point density and energy changes, but does not measure an
ordinary final-basis stationarity residual such as `[F, P]`, nor idempotency and
trace diagnostics for the final density. That should be probed locally before
adding damping/DIIS or route wiring.

Validation/status:

- Local ignored algebraic check ran and matched finite-difference/Fock
  directional derivatives to about `1e-10` absolute error.
- Reported git status clean and even with origin/main, apart from ignored local
  artifacts.

Deletion/shrinkage assessment:

- deleted: none.
- simplified: none; read-only/local audit.
- quarantined: algebraic check remains local ignored evidence.
- not deleted because: private RHF helpers remain the active diagnostic seam.
- exact remaining caller/blocker: compact real private SCF remains blocked by
  `:scf_not_converged`; next need is ordinary-final-basis commutator/idempotency
  residual evidence.

-- repo-manager@macmini

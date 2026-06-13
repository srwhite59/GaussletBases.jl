Pass 119 manager review

Accepted as a local ignored residual-aware probe.

The updated private SCF payload now reports materialized residual diagnostics on
the real compact route-smoke fixture directly:

- commutator residual: `1.3415879049821777e-05`
- spatial commutator residual: `6.7079395249108886e-06`
- density trace error: `3.55e-15`
- closed-shell idempotency error: `9.02e-17`

This matches the pass-117 local residual scale and confirms the nonconverged
state is a real stationarity/control issue, not missing reporting.

Decision:

- Do not route-wire RHF.
- Do not implement DIIS/mixing immediately.
- Next pass should be design-only for a private SCF-control/acceleration
  contract: residual used for control, bounded history, allowed controls,
  convergence labels, and failure labels.

Validation/status:

- Local ignored residual-aware probe elapsed: `99.980484375000003` seconds.
- Reported git status clean and even with origin/main, apart from ignored local
  artifacts.

Deletion/shrinkage assessment:

- deleted: none.
- simplified: the separate local residual script is no longer needed for this
  measurement once the private payload is used, though ignored artifacts remain
  useful evidence for now.
- quarantined: local probe script and summary remain ignored artifacts.
- not deleted because: historical local probes are useful until the SCF-control
  design settles.
- exact remaining caller/blocker: private compact RHF SCF remains blocked by
  `:scf_not_converged` with commutator residual around `1.34e-5`; next blocker
  is private acceleration/mixing contract design.

-- repo-manager@macmini

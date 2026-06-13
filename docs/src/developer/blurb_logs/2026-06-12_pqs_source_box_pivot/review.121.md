Pass 121 manager review

Accepted as a local ignored Fock-DIIS experiment.

The experiment provides enough evidence for a bounded private implementation
pass:

- best run: `max_history = 6`, `diis_regularization = 1.0e-12`
- commutator residual reduced to `2.9183295714024782e-8`
- fixed-point density delta reduced to `1.1441678879342554e-7`
- energy change reduced to `2.4158453015843406e-13`
- trace and idempotency remained at machine precision

The history-4/regularization-1e-12 case showed coefficient pathologies, so any
implementation must include bounded history, explicit regularization, and a
coefficient-pathology fallback.

Decision:

- Implement private diagnostic Fock DIIS controls next.
- Keep route-driver integration blocked.
- Keep tracked tests synthetic and focused.
- Re-run the compact route-smoke fixture only after the private implementation
  exists.

Validation/status:

- Local ignored DIIS experiment elapsed: `114.94605325000001` seconds after a
  local-script typo fix.
- Reported git status clean and even with origin/main, apart from ignored local
  artifacts.

Deletion/shrinkage assessment:

- deleted: none.
- simplified: the implementation target is now bounded Fock DIIS, not generic
  mixing.
- quarantined: DIIS experiment artifacts remain ignored local diagnostics.
- not deleted because: local experiments are useful comparison evidence until
  private SCF controls stabilize.
- exact remaining caller/blocker: private SCF has no built-in control payload or
  DIIS history yet.

-- repo-manager@macmini

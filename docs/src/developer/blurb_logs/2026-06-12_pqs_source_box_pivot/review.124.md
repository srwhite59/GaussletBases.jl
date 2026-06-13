Pass 124 manager review

Accepted as a local ignored validation probe.

The tracked private Fock-DIIS path reproduces the pass-121 local experiment
scale on the compact route-smoke fixture:

- status: `:blocked_pqs_multilayer_complete_core_shell_rhf_scf_payload`
- blocker: `:scf_not_converged`
- iterations: 100
- density delta: `1.1434019359102265e-07`
- commutator residual: `2.9183295714024782e-08`
- energy change: `2.4158453015843406e-13`
- trace/idempotency errors: machine precision
- DIIS used count: 99
- DIIS fallback/solve/pathology counts: 0

The corrected final-diagnostic gate is not the blocker; the last iteration
still misses density and residual tolerances under the strict `1e-8` gates.

Decision:

- Do not route-wire RHF.
- Do not loosen tolerances yet.
- Next pass should be a local ignored asymptote/parameter probe, primarily
  `max_iterations = 200` for the same history/regularization, optionally
  `max_history = 8`, to see whether the strict gate is reachable or whether a
  residual-floor/tolerance policy is needed.

Validation/status:

- Local ignored tracked-Fock-DIIS probe elapsed: `107.736029333` seconds.
- Reported git status clean and even with origin/main, apart from ignored local
  artifacts.

Deletion/shrinkage assessment:

- deleted: none.
- simplified: tracked private Fock-DIIS now replaces the separate local DIIS
  implementation for this measurement.
- quarantined: local probe script and summary remain ignored artifacts.
- not deleted because: prior ignored probes are useful comparison artifacts
  until convergence/tolerance policy is settled.
- exact remaining caller/blocker: tracked private Fock-DIIS remains blocked by
  `:scf_not_converged` at strict `1e-8` gates.

-- repo-manager@macmini

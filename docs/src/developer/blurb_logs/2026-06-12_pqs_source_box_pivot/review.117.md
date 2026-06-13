Pass 117 manager review

Accepted as a local ignored residual probe.

The residual probe establishes that the compact route-owned private RHF SCF
state is not merely failing a fixed-point delta label. The ordinary final-basis
commutator residual is about `1.3e-5`, and the occupied subspace delta is about
`2.2e-5`, while trace and idempotency are controlled. This supports the pass 116
conclusion: no factor-of-two Fock/energy formula bug is apparent, but the SCF
state is genuinely nonstationary at this fixture scale.

Decision:

- Do not wire RHF into the route driver.
- Do not add production damping/DIIS yet.
- Next pass should add residual reporting to the private SCF payload itself:
  ordinary final-basis commutator residual, trace error, idempotency error, and
  explicit residual metric labels. That gives future local probes and any later
  private controls a stable diagnostic contract.

Validation/status:

- Local ignored residual probe reran the compact route and completed in
  `99.661981999999995` seconds after one local-script typo was fixed.
- Reported git status clean and even with origin/main, apart from ignored local
  artifacts.

Deletion/shrinkage assessment:

- deleted: none.
- simplified: none; this was evidence gathering.
- quarantined: residual probe artifacts remain ignored local diagnostics.
- not deleted because: the private RHF SCF helper remains the active diagnostic
  seam.
- exact remaining caller/blocker: private compact RHF SCF remains blocked by
  `:scf_not_converged`; next blocker is missing residual-aware private payload
  diagnostics before any SCF-control implementation.

-- repo-manager@macmini

Pass 123 manager review

Accepted. This accepts the corrected combined pass-122/123 source changes.

The SCF helper now gates materialized convergence on the recomputed final
diagnostics for the returned final density. If iteration-input convergence
passes but the recomputed final residual/trace/idempotency checks fail, the
helper returns the blocked `:scf_not_converged` payload rather than claiming
`rhf_converged = true`.

The private Fock-DIIS controls remain diagnostic-only:

- no route wiring;
- no report aliases/options;
- no public API;
- no exports/artifacts;
- no fixture promotion.

Validation repeated by manager:

- `julia --project=. test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl`
  passed: 91/91.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  passed.
- `git diff --check` passed.

Deletion/shrinkage assessment:

- deleted: none.
- simplified: loose SCF control semantics are now carried by one compact
  private control payload, and convergence authority is the recomputed final
  residual diagnostics.
- quarantined: route integration remains blocked; compact route validation
  remains local/ignored.
- not deleted because: historical local probes remain useful until tracked
  private Fock-DIIS is validated on the compact route-smoke fixture.
- exact remaining caller/blocker: private Fock-DIIS exists in tracked code but
  has not yet been validated on the compact route-owned fixture.

-- repo-manager@macmini

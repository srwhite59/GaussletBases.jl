Pass 128 manager review

Accepted as a local ignored validation probe.

The default-control path is now confirmed:

- omitted Fock-DIIS `max_history` resolves to 8;
- compact route-smoke fixture converged under strict gates;
- iteration: 34;
- total energy: `-10.032119189804888`;
- density delta: `8.3320104149464669e-9`;
- commutator residual: `2.9455127052713248e-9`;
- spatial commutator residual: `1.4727563526356624e-9`;
- trace/idempotency errors: machine precision;
- coefficient-pathology fallbacks: 3, at iterations 8, 9, and 10;
- solve failures: 0;
- final recomputed diagnostics did not block convergence.

Decision:

- The SCF-control issue is settled enough for the next step.
- Next pass should be no-edit route-RHF adoption audit: identify the exact
  private driver slot, request object, required fields, and nonclaim/blocker
  vocabulary before any implementation.
- Do not route-wire RHF yet.
- Do not treat this private PQS seam as serious HF. Serious HF remains an
  `hfdmrg` concern, and CR2 downstream validation should come later through the
  CR2 agent when this line is actually ready.

Validation/status:

- Local ignored probe elapsed: `107.467666458` seconds.
- Reported git status was clean and even with origin/main.

Deletion/shrinkage assessment:

- deleted: none.
- simplified: callers can omit Fock-DIIS history and get the validated private
  route-smoke control default.
- quarantined: default-history probe script/table/summary remain ignored
  `tmp/work` artifacts.
- not deleted because: local probe artifacts document the private-control
  evidence chain.
- exact remaining caller/blocker: the private RHF helper converges when called
  by local probe code; the route driver does not yet have an audited private
  request/slot for RHF diagnostics.

-- repo-manager@macmini

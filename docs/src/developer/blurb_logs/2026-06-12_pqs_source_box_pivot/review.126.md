Pass 126 manager review

Accepted as a local ignored validation probe.

The confirmatory history-8 probe reproduced pass 125:

- converged: true
- iteration: 34
- total energy: `-10.032119189804888`
- density delta: `8.3320104149464669e-9`
- commutator residual: `2.9455127052713248e-9`
- spatial commutator residual: `1.4727563526356624e-9`
- trace/idempotency errors: machine precision
- DIIS used count: 30
- coefficient-pathology fallbacks: 3, at iterations 8, 9, and 10
- solve failures: 0
- final recomputed diagnostics did not block convergence

Decision:

- The next pass may make the small tracked private-control default change:
  Fock-DIIS `max_history` from 6 to 8.
- Keep strict tolerances, regularization, and coefficient guard unchanged.
- Do not route-wire RHF, promote fixtures, add public/report behavior, or treat
  this as serious HF. Serious HF remains an `hfdmrg`/downstream CR2 concern;
  this PQS seam remains private route-smoke/control diagnostics.

Validation/status:

- Local ignored probe elapsed: `105.134415291` seconds.
- Reported git status was clean and even with origin/main.

Deletion/shrinkage assessment:

- deleted: none.
- simplified: the remaining private-control decision is now a single default
  value, not tolerance policy or route design.
- quarantined: history-8 probe script/table/summary remain ignored `tmp/work`
  artifacts.
- not deleted because: prior ignored probes are useful comparison artifacts
  until the private default change is committed and reviewed.
- exact remaining caller/blocker: explicit history 8 converges on the compact
  route-smoke fixture; default Fock-DIIS still resolves to history 6 until the
  next tracked pass changes it.

-- repo-manager@macmini

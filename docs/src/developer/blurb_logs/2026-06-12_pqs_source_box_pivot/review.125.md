Pass 125 manager review

Accepted as a local ignored validation probe.

The asymptote probe answered the immediate control question:

- `max_history = 6` still does not converge by 200 iterations under the strict
  `1e-8` density/residual gates.
- `max_history = 8` converges by iteration 34 under the same strict gates.
- The converged history-8 case has commutator residual `2.9455127052713248e-9`,
  density delta `8.3320104149464669e-9`, energy change
  `2.8421709430404007e-14`, and trace/idempotency errors at machine precision.
- The history-8 case had 3 coefficient-pathology fallbacks and no solve
  failures; this is bounded enough to justify one confirmatory run but not yet
  enough to change defaults.

Decision:

- Do not loosen tolerances.
- Do not route-wire RHF.
- Do not promote the compact route-smoke fixture to a physics endpoint.
- Next pass should confirm `max_history = 8` with `max_iterations = 100`,
  record the fallback iterations, and preserve the same strict gates.

The user also reminded us that serious HF belongs to `hfdmrg`. That reinforces
the boundary here: this PQS RHF seam remains a private route-smoke/control
diagnostic, not the repo's serious HF package or physics-reference path.

Validation/status:

- Local ignored probe elapsed: `103.647220416` seconds.
- Reported git status was clean and even with origin/main.

Deletion/shrinkage assessment:

- deleted: none.
- simplified: the next choice narrows to confirming history 8 rather than
  loosening tolerances or changing route behavior.
- quarantined: asymptote probe script/table/summary remain ignored `tmp/work`
  artifacts.
- not deleted because: prior ignored probes remain useful comparison artifacts
  until the private control default is settled.
- exact remaining caller/blocker: default `max_history = 6` still blocks at the
  strict `1e-8` gates; history 8 has one successful local result and needs one
  confirmatory run before any default change.

-- repo-manager@macmini

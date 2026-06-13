Pass 140 manager review

Accepted.

Plain-language state:

- Be2/PQS now has a private support-window/order payload.
- It derives compact source-box windows and support counts from route-owned
  source boxes and parent dimensions.
- It records the key order mismatch explicitly:
  retained order is `(:pqs_left, :pqs_right, :product)`, while the candidate
  complete-core/shell support order is `(:product, :pqs_left, :pqs_right)`.
- It does not materialize raw product-box plans, support states, coefficients,
  source plans, final basis, H1, H1/J, Ham data, RHF, exports, or artifacts.

Decision:

- The next missing facts are raw product-box plan objects and PQS axis-local
  coefficients.
- Before wrapping those as route-owned data, add a focused fingerprint of the
  existing private raw-box route producer on the same Be2 route.
- This should be test-only: prove what the private shadow producer returns and
  whether it contains the objects needed by the future source-realization
  contract.

Deletion/shrinkage assessment:

- deleted: none.
- simplified: source-window/order facts are now one compact route-owned payload.
- quarantined: raw product-box probe output, raw product-box plans,
  coefficients, source-plan/final-basis/H1/H1-J/Ham materialization, RHF/SCF,
  WL payloads, public APIs, exports, artifacts, hfdmrg, and CR2 execution remain
  out of scope.
- not deleted because: support-window/source-plan/Ham readiness payloads are now
  the active private blocker chain.
- exact remaining caller/blocker: no route-owned diatomic source realization
  currently carries raw product-box plan objects, PQS axis-local coefficients,
  support-order permutation provenance, and an honest source-plan consumer
  shape.

-- repo-manager@macmini

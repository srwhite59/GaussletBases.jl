Pass 138 manager review

Accepted.

Plain-language state:

- Be2/PQS now has a private diatomic complete-core/shell source-plan payload.
- In the default path it blocks on the missing parent axis-bundle object.
- In the probe-enabled path it has the parent axis bundle and blocks more
  narrowly on the missing diatomic source-realization contract.
- The payload does not claim to return a `:pqs_multilayer_shell_source_plan`,
  and it does not build final basis, H1, H1/J, Ham data, RHF, exports, or
  artifacts.

Decision:

- The next blocker is no longer payload shape. It is the source-realization
  content: support states, support indices, retained/source coefficients,
  source-box ordering, and factor provenance for the diatomic
  `(:pqs_left, :product, :pqs_right)` route.
- Before coding a materializer, audit whether existing retained-unit records,
  raw product-box probe machinery, or source-box density-density helpers already
  carry the needed realization facts.
- Do not make raw product-box probe output route authority unless the audit
  finds a clean structured path and manager approves it.

Deletion/shrinkage assessment:

- deleted: none.
- simplified: the source-plan blocker now lives in a route-owned private
  payload instead of only in Ham-readiness missing-object lists.
- quarantined: fake `:pqs_multilayer_shell_source_plan` objects, one-center
  shellification semantics, WL adapters, final-basis/H1/H1-J/Ham
  materialization, RHF/SCF, public API, exports, artifacts, hfdmrg, and CR2
  execution remain out of scope.
- not deleted because: the default and probe-enabled Be2 fingerprints still
  document distinct route states that matter for the next producer decision.
- exact remaining caller/blocker: no diatomic source-realization contract
  currently supplies the support states/indices and retained/source
  coefficients needed by the existing final-basis/H1 source-plan consumer
  shape.

-- repo-manager@macmini

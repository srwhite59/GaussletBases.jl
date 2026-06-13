Pass 139 manager review

Accepted.

Plain-language state:

- Be2/PQS retained units currently carry route metadata, retained ranges, counts,
  provenance, and weight semantics, but not retained/source coefficient
  matrices.
- Be2/PQS source boxes do carry enough x/y/z window information to define
  support windows against a parent axis bundle.
- Existing raw product-box probe plumbing can build richer raw plan objects
  internally, but the current probe result is metadata-only and explicitly
  private/probe-only.
- Forcing Be2 into the existing `core_then_shell` source-plan consumer shape
  would require an explicit order/permutation contract, because the natural
  retained order `(:pqs_left, :pqs_right, :product)` differs from the plausible
  support order `(:product, :pqs_left, :pqs_right)`.

Decision:

- Add a compact private support-window/order payload next.
- It should derive source-box windows and support counts from route-owned source
  boxes and parent dimensions.
- It should record the retained order, candidate support order, and that a
  retained-to-support permutation is required.
- It should not use raw product-box probe output as route authority and should
  not materialize coefficients, final basis, H1, H1/J, Ham data, RHF, WL, or
  exports.

Deletion/shrinkage assessment:

- deleted: none.
- simplified: next pass should move source-window/order facts out of ad hoc
  route-skeleton interpretation into one compact route-owned private payload.
- quarantined: raw product-box probe results, shell/support-row contraction,
  final-basis/H1/H1-J/Ham materialization, density-density materialization,
  RHF/SCF/Fock, WL payloads, public APIs, exports, artifacts, hfdmrg, and CR2
  execution remain out of scope.
- not deleted because: current Be2 fingerprints and source-plan payload remain
  active compatibility checks while the route-owned support-window boundary is
  being defined.
- exact remaining caller/blocker: no diatomic source materializer currently
  converts source boxes plus parent axis bundles into coefficient blocks,
  support-order provenance, and an honest source-plan consumer shape.

-- repo-manager@macmini

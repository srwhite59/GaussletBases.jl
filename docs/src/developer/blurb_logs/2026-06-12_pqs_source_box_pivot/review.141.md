Pass 141 manager review

Accepted.

Plain-language state:

- The existing private raw-box route producer is directly callable for the
  probe-enabled Be2/PQS fixture.
- It returns raw product-box plans for left/right PQS units, raw PQS plan views
  with axis-local coefficient arrays, the product/doside unit, descriptor, and
  pair inventory.
- It is still explicitly `:private_shadow_only` and not route authority.
- This proves the missing raw-box/coefficient objects are constructible, but not
  yet wrapped as a route-owned private payload for the driver.

Decision:

- Add a private raw-box route payload next.
- It should wrap the underlying raw-box route producer, not the metadata-only
  raw product-box probe result.
- It should be private/candidate-only, compactly summarized, and wired into the
  existing private assembly payload chain.
- It should not claim to materialize a source plan, final basis, H1, H1/J, Ham
  data, RHF, exports, public API, hfdmrg, or CR2 handoff.

Deletion/shrinkage assessment:

- deleted: none.
- simplified: producer availability is now proven by one focused fingerprint
  instead of inferred from raw-box probe metadata.
- quarantined: the producer remains private-shadow evidence; it is not yet
  route authority or a source-plan materializer.
- not deleted because: the fingerprint is useful until the private raw-box
  payload exists and can replace direct test calls.
- exact remaining caller/blocker: no route-owned private payload currently
  carries raw product-box plans, PQS axis-local coefficient objects, product
  unit summary, pair inventory, and support-order provenance for downstream
  diatomic source realization.

-- repo-manager@macmini

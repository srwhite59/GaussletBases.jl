Pass 136 manager review

Accepted and committed by the manager as:

```text
c102916c Fingerprint probe-enabled Be2 readiness
```

Plain-language state:

- The default Be2/PQS readiness path still documents the unprobed route state:
  no parent axis-bundle object and no diatomic complete-core/shell source-plan
  producer.
- With `probe_parent_axis_construction = :auto`, the same focused Be2/PQS route
  can expose the parent axis-bundle object through existing structured carry.
- In both paths, the private Ham payload remains correctly blocked.
- The narrower remaining blocker is now the missing route-owned diatomic
  complete-core/shell source-plan producer.

Decision:

- Do not make probe-enabled behavior the default yet.
- Do not jump straight into producer implementation.
- Next pass should be a no-edit producer-contract audit: decide what a diatomic
  PQS complete-core/shell source-plan producer should consume and return, and
  whether it should be a new source-plan object or a carefully bounded adapter
  around existing source-box route skeleton objects.

Deletion/shrinkage assessment:

- deleted: none.
- simplified: the parent-axis-bundle blocker is separated from the source-plan
  producer blocker.
- quarantined: probe-enabled parent-axis construction remains a focused
  fingerprint path, not a default route behavior; final basis, H1, H1/J, RHF,
  WL payloads, exports, artifacts, hfdmrg, and CR2 execution remain out of
  scope.
- not deleted because: both default and probe-enabled fingerprints are currently
  useful to document the two route states.
- exact remaining caller/blocker: Be2/PQS needs a route-owned diatomic
  complete-core/shell source-plan producer for the
  `(:pqs_left, :product, :pqs_right)` source-box skeleton before final-basis/H1
  and Ham payload materialization can proceed.

-- repo-manager@macmini

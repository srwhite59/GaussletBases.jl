Pass 142 manager review

Accepted.

Plain-language state:

- Be2/PQS now has a private raw-box route payload wired through
  `cartesian_assembly(...)`.
- In the default route it blocks on the missing parent axis-bundle object.
- In the probe-enabled route it carries the private raw-box producer, raw
  product-box plan summaries, raw PQS plan summaries, product/doside unit
  summary, and pair-inventory summary.
- The direct test-only producer calls from pass 141 were removed in favor of
  assertions through the route-owned private payload.
- It still does not materialize a source plan, final basis, H1, H1/J, Ham data,
  RHF, exports, artifacts, hfdmrg, or CR2 handoff.

Decision:

- The next question is the source-realization mapping, not availability of raw
  objects.
- Audit how the raw-box payload could satisfy the existing
  `:pqs_multilayer_shell_source_plan` consumer shape:
  product unit as core/body sector, left/right raw PQS plans as shell/source
  sectors, explicit ordering/permutation, support states/indices, coefficient
  matrix shape, bundles, metrics, and provenance.
- Do not implement the materializer until that mapping is checked.

Deletion/shrinkage assessment:

- deleted: direct/ad hoc raw-box producer assertions from the focused test.
- simplified: tests now prefer the private route-owned `diatomic_raw_box_route_payload`.
- quarantined: raw-box payload is private candidate data, not public route
  authority or a source-plan materializer.
- not deleted because: support-window/source-plan/Ham readiness payloads remain
  active blockers and review surfaces.
- exact remaining caller/blocker: no source-realization materializer currently
  converts the raw-box payload plus support-order/permutation facts into an
  honest `:pqs_multilayer_shell_source_plan` consumer shape.

-- repo-manager@macmini

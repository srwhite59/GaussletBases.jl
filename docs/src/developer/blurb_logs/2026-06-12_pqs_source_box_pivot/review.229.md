Pass 229 post-hoc manager review

The user redirected pass 229 away from the queued supplement-provider audit and
toward quarantining the H2 463 fake-PQS endpoint.

Accepted commits:
- `208df0ad Relabel H2 source-backed route as fake PQS`
- `26ecc737 Quarantine fake PQS endpoint readiness`

Review judgment:
- The first commit made the fake status visible in filenames and under the
  `fake_pqs/*` artifact group, but still allowed stable `route/*` and
  `physics/*` fields to look like an independent endpoint.
- The corrective commit fixed the important leak by making
  `physics/endpoint_ready = false` and placing fake/independent markers
  directly in `route/*` and `physics/*`.

Remaining caveat:
- Some internal type/status names still contain `physical_gausslet` because
  the fake reproduction reuses that implementation path. That is tolerable
  for now because the artifact-facing contract is quarantined, but it remains
  cleanup pressure once an independent PQS H2 route exists.

Line-count note:
- The corrective relabel was intentionally allowed to be line-positive. This
  was a safety correction to prevent downstream misclassification, not a normal
  feature pass. The line-reduction rule resumes after this audit/correction.

State:
- `main` is clean/even at `26ecc737`.
- No loop continuation was started from this audit entry.

-- repo-manager@macmini

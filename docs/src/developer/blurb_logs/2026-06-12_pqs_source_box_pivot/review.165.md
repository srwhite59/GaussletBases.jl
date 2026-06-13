Pass 165 review - accepted

Verdict: accepted.

The audit answered the right question before implementation: CR2 wants a
durable read-only inspection artifact, not HamV6, not solver readiness, and not
a public API. The proposed first schema is appropriate because it stores plain
arrays and compact scalar metadata under route groups rather than private route
payload objects.

Important conclusions:

- JLD2 is already a project dependency, so the first artifact can be JLD2
  without dependency churn.
- JSON is not established in the main source path, so TSV is the conservative
  first companion fingerprint.
- PQS can fill the useful read-only inspection pieces now:
  final H1, pre-final density interaction, final-to-pre-final coefficients,
  pre-final/support weights, raw pair numerator/provenance, nuclear/electron
  metadata, ordering labels, and readiness flags.
- White-Lindsey can fill a related final-basis density-density inspection shape
  through existing low-order/final-basis export paths, but it does not have PQS
  pre-final/source-box support-row internals. Those fields should be
  `not_applicable` or `unavailable`, not synthesized.
- The first implementation should write the fixed schema with PQS populated and
  WL placeholders if needed, rather than changing the schema later.

The next implementation pass should be private, artifact-only, and line-budget
negative across `src` + `test`. It should add the smallest JLD2+TSV inspection
writer for the current Be2/PQS handoff and shrink older duplicate readiness
assertions in the focused Be2 fingerprint test enough to pay for the new code.

-- repo-manager@macmini

Pass 110 manager review

Accepted as a read-only route-integration audit, but not accepted as the next
implementation direction.

The audit usefully identifies the eventual route seam: the diagnostic route
payload already carries source payload, final basis, H1 payload, density inputs,
and H1/J diagnostic data; the missing route-facing ingredients are explicit
`electron_count`, explicit `fixture_role`, and compact SCF controls. The audit
also correctly rejects inferring electron count from center charges and rejects
scalar report-field clouds.

However, after pass 109 the immediate risk is the private SCF payload contract
itself, not route adoption. The next pass should not wire a route request object
yet. First inspect the private SCF loop against the pass-102 RHF contract,
especially final-density/final-one-step consistency and convergence semantics.
Then run a local `tmp/work` compact real PQS SCF probe if the SCF contract review
does not find a bug.

Validation:

- No code changes.
- No Julia commands.
- Reported git status: clean and even with origin/main.

Deletion/shrinkage assessment:

- deleted: none.
- simplified: none; read-only audit.
- quarantined: route RHF adoption remains blocked.
- not deleted because: no implementation occurred.
- exact remaining caller/blocker: route-driver RHF caller remains blocked on
  explicit diagnostic request design, but before that we need a private SCF
  contract review and real compact local probe.

-- repo-manager@macmini

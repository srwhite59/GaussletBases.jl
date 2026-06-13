Pass 104 response

Docs-only fixture policy pass completed.

Files changed:

- `docs/src/developer/pqs_source_box_fixture_policy.md`
- `docs/src/developer/pqs_source_box_operator_framework.md`

Policy summary:

- Added a short tracked fixture policy note that separates:
  - route smoke;
  - convention diagnostic;
  - oracle/debug;
  - physics endpoint.
- Recorded current compact fixture facts for:
  - tracked Z=1 H1 seam;
  - direct structured H1/J convention probe;
  - one-center source-box driver H1/J dry-run.
- Recorded nonclaims:
  - compact H1/J materialization is not RHF readiness;
  - final dimension `223` is not physics acceptance;
  - self-Coulomb alone is not endpoint validation;
  - shell/support-row, fixed-block, and explicit-box paths remain oracle/debug.
- Listed coupled fixture parameters that must move as reviewed fixture
  families before physics comparison.
- Added the before-RHF policy: choose route-smoke-only versus physics endpoint
  first; if physics endpoint, choose target, electron count, closed-shell rule,
  fixture size, reference/error thresholds, and timing threshold before code.

Framework pointer:

- Added a two-line pointer from
  `docs/src/developer/pqs_source_box_operator_framework.md` to
  `docs/src/developer/pqs_source_box_fixture_policy.md`.
- No broad framework rewrite was done.

Validation:

- `git diff --check`
  - passed
- `LC_ALL=C rg -n "[^[:ascii:]]" docs/src/developer/pqs_source_box_fixture_policy.md docs/src/developer/pqs_source_box_operator_framework.md`
  - passed: no non-ASCII matches
- `git status --short --branch`

```text
## main...origin/main
 M docs/src/developer/pqs_source_box_operator_framework.md
?? docs/src/developer/pqs_source_box_fixture_policy.md
```

Deletion/shrinkage report:

- deleted: none.
- simplified: fixture-role policy is now centralized in one tracked note,
  reducing pressure for future RHF blurbs/tests to restate route-smoke versus
  physics-endpoint boundaries.
- quarantined: compact H1/J fixture facts remain route-smoke/convention
  diagnostics unless a later fixture card promotes one; side13/q-ladder probe
  history remains ignored/historical evidence.
- not deleted because: this was a docs-only policy pass; tracked H1 seam tests
  and H1/J report aliases remain active until a replacement smoke/policy-backed
  test is explicitly assigned.
- exact remaining caller/blocker: before RHF implementation, manager still must
  choose route-smoke-only versus physics endpoint and, for endpoint work, choose
  target system, electron count, fixture size, reference/error threshold, and
  timing threshold.

-- repo-doer@macmini

Pass 109 manager review

Accepted.

The implementation adds the first private closed-shell RHF SCF payload and
keeps it diagnostic/prototype-only. It consumes the existing input contract,
initial-density payload, and one-step payload; it does not wire the route
driver, add report aliases/options, export API, or promote any fixture. The
focused synthetic test is small and deterministic.

Validation repeated by manager:

- `julia --project=. test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl`
  passed: 30/30.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  passed.
- `git diff --check` passed.

Follow-up before route wiring:

- The private RHF input contract intentionally requires explicit
  `electron_count` and `fixture_role`. The route driver currently should not
  infer these from nuclei or silently promote compact route-smoke fixtures.
  Before adding any driver integration, audit the exact seam for carrying those
  explicit diagnostic inputs without adding a scalar report-field cloud or a
  public route option.

Deletion/shrinkage assessment:

- deleted: none.
- simplified: RHF iteration logic now has a private payload instead of pressure
  to grow the one-step or H1/J helpers.
- quarantined: SCF remains private diagnostic/prototype behavior, with route,
  report, public API, export, and artifact flags false.
- not deleted because: no previous SCF payload existed to replace.
- exact remaining caller/blocker: no route-driver caller. Route integration is
  blocked on an explicit diagnostic-input seam for electron count and fixture
  role, plus a manager decision about whether to run the cold route smoke.

-- repo-manager@macmini

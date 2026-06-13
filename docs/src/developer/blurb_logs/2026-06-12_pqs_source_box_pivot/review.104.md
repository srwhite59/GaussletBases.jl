Pass 104 review

Accepted docs-only fixture policy.

The new `docs/src/developer/pqs_source_box_fixture_policy.md` records the
fixture roles and nonclaims we need before RHF work:

- compact fixtures are route-smoke/convention diagnostics unless explicitly
  promoted;
- `final_dimension == 223`, H1/J materialization, and self-Coulomb alone are
  not physics acceptance;
- shell/support-row, fixed-block, and explicit-box paths remain oracle/debug;
- RHF requires an explicit route-smoke versus physics-endpoint decision.

The framework pointer is narrow and appropriate. Validation passed:

- `git diff --check`
- ASCII scan of the fixture policy and framework doc.

Next pass may add the first private RHF input-contract object, but only as a
blocked/available input contract. No Fock, no SCF, no driver wiring.

-- repo-manager@macmini

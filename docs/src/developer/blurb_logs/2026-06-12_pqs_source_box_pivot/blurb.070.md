Purpose:
  Audit the remaining H1 support-operator assembly seam before promoting more
  probe/test logic into production. This is read-only/docs-only unless you find
  a trivial stale wording fix.

Context:
  The tracked PQS H1 gate still has local helpers for:

  - support-space kinetic assembly;
  - support-space electron-nuclear assembly from PGDG Gaussian factor terms and
    Coulomb coefficients.

  The rectangular support-product primitive is now production-owned by
  `_pqs_multilayer_support_product_matrix(...)`, and the complete final-basis
  assembly is production-owned by
  `pqs_multilayer_complete_core_shell_final_basis(...)`.

Task:
  Audit whether the next production seam should be route-owned support
  one-body assembly for the multi-layer PQS plan.

  Answer:

  1. Is support kinetic assembly safe to promote now as a helper over
     `plan.metrics` and `plan.core/shell_support_states`?

  2. Is support electron-nuclear assembly safe to promote now, or does it need
     a clearer by-center convention first?

  3. For electron-nuclear, what should be the convention?
     Explicitly address:
     - uncharged by-center matrix versus charged Hamiltonian-stage matrix;
     - sign convention;
     - center location and off-origin center readiness;
     - source of Gaussian factor terms;
     - whether old fixed-block or WL code is oracle only.

  4. What helper, if any, should be implemented next?
     Name its proposed inputs/outputs and nonclaims.

  5. What existing test/probe-local code would shrink if that helper existed?

Deliverable:
  Add a concise note to
  `docs/src/developer/pqs_near_term_final_basis_realization_plan.md` under the
  driver/assembly section, or add a small adjacent developer note if that file
  would become too crowded.

Do not:
  - implement the helper in this pass;
  - change tests;
  - add H1/RHF/IDA/density-density features;
  - change fixture-rule policy;
  - edit ignored `tmp/work` probes;
  - add exports or driver wiring.

Validation:
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - `git diff --check`

Deletion/shrinkage report required:
  - what old code, test, metadata, or compatibility path became unnecessary;
  - what was deleted or simplified;
  - if nothing was deleted, why no existing surface was made obsolete;
  - whether any new test replaces/shrinks older coverage or is genuinely new
    live-contract coverage;
  - any remaining stale or duplicate surfaces to retire next.

-- repo-manager@macmini

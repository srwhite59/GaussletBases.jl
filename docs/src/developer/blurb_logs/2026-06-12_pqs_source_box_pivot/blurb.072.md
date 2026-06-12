Purpose:
  Nail down the route-owned PQS support electron-nuclear convention before any
  implementation. This is a convention/design audit, not a coding pass.

Context:
  The H1 gate still has a local `_pqs_h1_support_nuclear_matrix` helper. It
  currently handles only the centered/origin case using
  `pgdg_intermediate.gaussian_factor_terms` and Coulomb expansion coefficients.

  Pass 070 said electron-nuclear promotion should wait until these are explicit:

  - separated by-center output;
  - uncharged matrix versus Hamiltonian-stage charge application;
  - sign convention;
  - center location/off-origin readiness;
  - source of centered Gaussian factor terms;
  - old fixed-block/WL surfaces as oracle only.

Task:
  Audit the existing PQS retained by-center nuclear helpers and WL by-center
  convention, then write the exact support electron-nuclear helper contract that
  should be implemented next.

  Answer:

  1. What should the helper be named?
     Candidate: `pqs_multilayer_support_electron_nuclear_by_center_matrices`.

  2. What should its inputs be?
     Include plan, Coulomb expansion, center records, and axis/source factor
     data as needed.

  3. What should it return?
     Prefer one support-space matrix per center, ordered over core support rows
     followed by shell support rows.

  4. What is the sign convention?
     Decide whether each uncharged center matrix is `-1/r_A` or `+1/r_A`, and
     state how charge application in final H assembly should work.

  5. How should centered/origin and off-origin centers get Gaussian factor
     terms?
     Identify the existing centered factor helper to reuse if available.

  6. What should remain oracle/reference only?

  7. What exact test-local code would shrink if the helper lands?

Deliverable:
  Add a concise convention note to
  `docs/src/developer/pqs_near_term_final_basis_realization_plan.md`, or an
  adjacent small developer note if clearer.

Do not:
  - implement the nuclear helper;
  - change tests;
  - change H1/RHF/IDA/density-density behavior;
  - change fixture-rule policy;
  - add exports or driver wiring;
  - edit ignored `tmp/work` probes.

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

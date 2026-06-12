Purpose:
  Promote only the safe support-kinetic assembly helper for multi-layer PQS
  plans, then shrink the tracked H1 gate's local kinetic helper.

Context:
  Pass 070 concluded:

  - support kinetic is safe: it only needs plan support states and
    `plan.metrics.{x,y,z}.overlap/kinetic`;
  - support electron-nuclear is not ready for promotion until separated
    by-center, sign, charge, centered/off-origin Gaussian factor, and oracle
    conventions are explicit.

Task:
  Add a narrow helper near the multi-layer plan code, with a name like:

  ```julia
  pqs_multilayer_support_kinetic_matrix(plan)
  ```

  It should:

  - require a `:pqs_multilayer_shell_source_plan`;
  - require `plan.status == :available_pqs_multilayer_shell_source_plan`;
  - build the support-space kinetic matrix over `vcat(plan.core_support_states,
    plan.shell_support_states)`;
  - use the standard three-term Cartesian product form:
    `Kx*Sy*Sz + Sx*Ky*Sz + Sx*Sy*Kz`;
  - return either the matrix directly or a compact result with matrix plus
    status/provenance, following the local style;
  - make no final-basis transfer, H1, nuclear, IDA, RHF, driver, export, or
    artifact claims.

  Then update `test/nested/pqs_direct_retained_final_h1_runtests.jl` so it uses
  the route-owned kinetic helper instead of local `_pqs_h1_support_kinetic_matrix`.
  Delete the local kinetic helper if it becomes unnecessary.

Do not:
  - implement support electron-nuclear in this pass;
  - change nuclear sign/charge convention;
  - add H1/RHF/IDA/density-density features;
  - change fixture-rule policy;
  - add exports or driver wiring;
  - edit ignored `tmp/work` probes.

Validation:
  - `julia --project=. -e 't = @elapsed include("test/nested/pqs_direct_retained_final_h1_runtests.jl"); println("elapsed_s=", t)'`
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - `git diff --check`

Docs:
  Update `docs/src/developer/pqs_near_term_final_basis_realization_plan.md`
  only to mark support kinetic helper implemented. Leave the nuclear wait note.

Deletion/shrinkage report required:
  - what old code, test, metadata, or compatibility path became unnecessary;
  - what was deleted or simplified;
  - if nothing was deleted, why no existing surface was made obsolete;
  - whether any new test replaces/shrinks older coverage or is genuinely new
    live-contract coverage;
  - any remaining stale or duplicate surfaces to retire next.

-- repo-manager@macmini

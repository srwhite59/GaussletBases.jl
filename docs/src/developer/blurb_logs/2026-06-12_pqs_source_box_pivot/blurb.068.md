Purpose:
  Shrink the tracked PQS H1 gate now that
  `pqs_multilayer_complete_core_shell_final_basis(plan; ...)` exists. The test
  should consume the route-owned assembly helper as the normal path instead of
  keeping manual final-basis assembly as a parallel authority.

Context:
  Pass 067 added the helper and a comparison in
  `test/nested/pqs_direct_retained_final_h1_runtests.jl`. The fixture still
  constructs `fixture.final_basis` manually by building `core_overlap`,
  `core_shell_overlap`, and `shell_overlap`, then calling
  `CartesianFinalBasisRealization.pqs_complete_core_shell_final_basis(...)`.

Task:
  Refactor the H1 test fixture so the final basis comes from:

  ```julia
  plan = pqs_multilayer_shell_source_plan(...)
  final_basis = pqs_multilayer_complete_core_shell_final_basis(plan; ...)
  ```

  Delete or simplify the duplicated manual final-basis assembly in the test
  fixture. Keep only the support matrix helpers still needed for final overlap,
  kinetic, nuclear, and H1 checks.

  The final test should still verify the same physics/workflow facts:

  - complete core/shell count facts;
  - final overlap identity;
  - final H finite/symmetric;
  - ordinary H1 solve;
  - H1 reference value and fixed-block oracle agreement;
  - no generalized-overlap solve, no IDA/RHF/driver/export/artifact claim.

Do not:
  - change production source unless the refactor exposes a trivial bug;
  - add new tests or broaden assertions;
  - add H1/RHF/IDA/density-density features;
  - change fixture-rule policy;
  - edit ignored `tmp/work` probes in this pass.

Validation:
  - `julia --project=. -e 't = @elapsed include("test/nested/pqs_direct_retained_final_h1_runtests.jl"); println("elapsed_s=", t)'`
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

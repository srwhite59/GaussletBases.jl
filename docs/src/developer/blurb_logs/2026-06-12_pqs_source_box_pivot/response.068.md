Pass 068 response: PQS H1 gate assembly-helper shrinkage

Summary:
- Refactored `test/nested/pqs_direct_retained_final_h1_runtests.jl` so `_pqs_h1_complete_fixture()` now builds:
  - `plan = GaussletBases.pqs_multilayer_shell_source_plan(...)`
  - `final_basis = GaussletBases.pqs_multilayer_complete_core_shell_final_basis(plan; ...)`
- Removed the duplicated manual final-basis construction from the H1 fixture.
- Kept the test-local support product/kinetic/nuclear helpers because the final overlap, kinetic, nuclear, and H1 checks still consume support-space matrices.
- Preserved the active H1 gate facts:
  - complete core/shell counts;
  - final overlap identity;
  - finite symmetric final H;
  - ordinary H1 solve;
  - H1 reference value and fixed-block oracle agreement;
  - no generalized-overlap, IDA, RHF, driver, export, or artifact claim.

Validation:
- `julia --project=. -e 't = @elapsed include("test/nested/pqs_direct_retained_final_h1_runtests.jl"); println("elapsed_s=", t)'`
  - PASS: 31 tests, elapsed_s=4.188590334
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - PASS: `load ok`
- `git diff --check`
  - PASS

Diff:
- `test/nested/pqs_direct_retained_final_h1_runtests.jl`
  - 17 insertions, 99 deletions

Deletion/shrinkage report:
- Became unnecessary:
  - the test-local `PQSH1CCPM` alias;
  - `_pqs_h1_axis_metrics`;
  - `_pqs_h1_core_indices_and_states`;
  - manual staged layer/descriptor/shell-plan construction inside the fixture;
  - manual `core_overlap`, `core_shell_overlap`, and `shell_overlap` assembly for final-basis construction;
  - the separate helper-vs-manual final-basis comparison block in the test body.
- Deleted/simplified:
  - the H1 fixture now uses the route-owned multi-layer source plan and final-basis helper as its normal path;
  - the test body checks that the fixture final basis came from the multi-layer source plan/helper instead of comparing two authorities.
- No new tests were added. This pass shrank existing coverage rather than adding new scaffolding.
- Remaining duplicate/stale surface:
  - support-space product helpers remain in this file only because the H1 gate still assembles final one-body operator checks from support matrices.

-- repo-doer@macmini

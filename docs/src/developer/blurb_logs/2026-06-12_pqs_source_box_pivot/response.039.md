Pass 039 complete.

Summary:
- Added one compact durable H1 gate:
  - `test/nested/pqs_direct_retained_final_h1_runtests.jl`
- Added it to the default nested runner:
  - `test/nested/runtests.jl`
- Shrunk old helper-vocabulary assertions in:
  - `test/nested/bond_aligned_diatomic_high_order_recipe_opt_in_source_construction_integration_runtests.jl`

H1 gate result:
- Direct command:
  - `julia --project=. -e 't = @elapsed include("test/nested/pqs_direct_retained_final_h1_runtests.jl"); println("elapsed_s=", t)'`
- Result:
  - `PQS direct-retained final H1 gate | 29 passed`
  - testset time: `4.0s`
  - elapsed: `4.961041333s`
- The gate checks:
  - source dims/count `5 x 5 x 5 / 125`;
  - boundary/final retained count `98 / 98`;
  - final overlap identity error below tolerance;
  - retained overlap, kinetic, and centered by-center nuclear direct-boundary flags;
  - no raw source one-body block materialization on the active retained path;
  - final Hamiltonian finite/symmetric;
  - final Hamiltonian and H1 eigenvalue match shell-support oracle within `1.0e-12`;
  - ordinary symmetric H1 diagnostic, not generalized-overlap solve;
  - no IDA, density-density, RHF, driver, export, or artifact claim.

Old integration shrink:
- Removed the broad per-term safe-term matrix/oracle assertion loop over `_pqs_current_route_safe_term_matrices(...)`.
- Removed detailed safe-term diagnostic metadata assertions such as whole-route matrix consumer, retained dimension, pair count, support-local oracle fields, and source-box availability flags.
- Removed the broad per-term authority-comparison loop over `_pqs_current_route_safe_term_authority_comparison(...)`.
- Removed detailed authority-comparison metadata assertions such as authority source/field/shape and compared-term metadata.
- Retained a small private-diagnostic smoke:
  - helper object kind/status;
  - expected term tuple;
  - global/max authority error bounds;
  - finite-output flags;
  - empty unavailable terms.

Validation:
- New H1 gate passed directly.
- The edited slow integration file was not run through the full integration harness because it requires the top-level fixture harness and is broad/slow. A direct include without the harness fails before reaching the edited block on missing `bond_aligned_homonuclear_qw_basis`, which confirms it is harness-bound rather than standalone.
- Narrow validation for the edited integration file:
  - `julia --project=. -e 'Meta.parseall(read("test/nested/bond_aligned_diatomic_high_order_recipe_opt_in_source_construction_integration_runtests.jl", String)); println("integration parse ok")'` passed.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.
- Extra whitespace scan for the new untracked test file found no trailing whitespace.

Deletion/shrinkage report:
- Tracked integration diff removes `33` lines of old helper-vocabulary assertions.
- The new gate is `375` lines and adds `29` live-contract assertions.
- Runtime assertion pressure in the old helper-vocabulary section drops by about `95` checks; the new gate adds `29`, so this coverage area shifts from old implementation vocabulary to a smaller scientific/workflow contract.
- The H1 gate replaces old coverage pressure rather than duplicating it: it directly checks retained overlap/kinetic/nuclear construction, final-basis transfer, Hamiltonian assembly, and H1 oracle agreement.
- No private CCPM helper was deleted in this pass.
- Private CCPM surfaces left for future deletion/quarantine:
  - `_pqs_current_route_safe_term_matrices(...)`
  - `_pqs_current_route_safe_term_authority_comparison(...)`
  - related support-local oracle/debug vocabulary in `CartesianContractedParentMetrics.jl`

-- repo-doer@macmini

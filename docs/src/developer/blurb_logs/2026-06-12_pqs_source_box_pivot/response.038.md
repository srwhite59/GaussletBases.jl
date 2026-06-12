Pass 038 complete.

Summary:
- Updated and ran the ignored developer probe `tmp/work/pqs_final_basis_h1_probe.jl`.
- Wrote the requested concise artifact:
  - `tmp/work/pqs_direct_retained_final_h1_probe_summary.txt`
- The probe assembles the explicit PQS final-basis H1 seam from:
  - one `5 x 5 x 5` raw product source box;
  - PQS boundary source-mode retained rule with retained count `98`;
  - direct retained overlap;
  - direct retained kinetic;
  - direct retained centered by-center nuclear;
  - `CartesianFinalBasisRealization` final-basis realization;
  - CPBM final by-center nuclear transfer;
  - CPBM final one-electron Hamiltonian assembly.
- The active probe path does not call `_pqs_current_route_safe_term_matrices(...)`.

Probe result:
- `source_mode_dims = (5, 5, 5)`
- `source_mode_count = 125`
- `boundary_count = 98`
- `shell_support_count = 98`
- `final_retained_count = 98`
- `final_overlap_identity_error = 2.2360255762214364e-14`
- retained direct-boundary flags:
  - overlap: `true`
  - kinetic: `true`
  - centered by-center nuclear: `true`
- raw source one-body blocks materialized on active retained path:
  - overlap: `false`
  - kinetic: `false`
  - centered by-center nuclear: `false`
- source-space input used by active retained path:
  - overlap: `false`
  - kinetic: `false`
  - centered by-center nuclear: `false`
- Hamiltonian checks:
  - shape: `(98, 98)`
  - finite: `true`
  - symmetry error: `2.7755575615628914e-17`
- H1:
  - direct final-basis lowest eigenvalue: `-0.08171962129085239`
  - shell-support oracle lowest eigenvalue: `-0.081719621290852015`
  - eigenvalue delta: `3.7470027081099033e-16`
- Final Hamiltonian versus shell-support oracle max error:
  - `1.3322676295501878e-15`
- Solve convention:
  - ordinary symmetric eigensolve: `true`
  - generalized overlap solve: `false`

Coarse probe timings from the recorded run:
- total elapsed: `2.986572s`
- parent axis setup: `0.053107s`
- shell layer: `1.166869s`
- shell realization plan: `0.033901s`
- final basis: `0.076325s`
- direct retained overlap block: `0.062749s`
- direct retained kinetic block: `0.076376s`
- direct retained centered nuclear block: `0.405413s`
- final overlap transfer: `0.017613s`
- final kinetic transfer: `0.000053s`
- final nuclear transfer: `0.022356s`
- final Hamiltonian assembly: `0.096599s`
- shell-support oracle construction: `0.075016s`
- H1 eigensolve: `0.000559s`

Old helper/test pressure audit:
- `_pqs_current_route_safe_term_matrices(...)` is still defined in `src/CartesianContractedParentMetrics.jl` as private diagnostic/oracle machinery and feeds `_pqs_current_route_safe_term_authority_comparison(...)`.
- The only live `test/` caller found is in `test/nested/bond_aligned_diatomic_high_order_recipe_opt_in_source_construction_integration_runtests.jl`, where it still checks old current-route safe-term matrix vocabulary and support-local oracle comparisons.
- `pqs_source_shell_projected_one_body_matrix(...)` remains a compact CFBR oracle/helper. It is still useful for shell-support oracle comparisons and should not be deleted by this pass.
- If the H1 probe is promoted to a compact gate, the best shrink target is the slow integration section around the `_pqs_current_route_safe_term_matrices(...)` call in `bond_aligned_diatomic_high_order_recipe_opt_in_source_construction_integration_runtests.jl`. The H1 gate can cover overlap/kinetic/nuclear retained-to-final scientific behavior more directly, while leaving the old helper as private oracle/debug until a caller-driven deletion pass removes or quarantines its authority-comparison use.

Validation:
- `julia --project=. tmp/work/pqs_final_basis_h1_probe.jl` passed.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.

Deletion/shrinkage report:
- Nothing was deleted or shrunk in tracked source/tests because this pass was explicitly probe/readiness first and the H1 gate has not yet been promoted.
- The probe updated an ignored `tmp/work` script and wrote an ignored artifact only.
- No permanent test was added.
- No new public API, driver route, IDA, density-density, RHF, GTO, export, or artifact path was added.
- Exact proposed old coverage to shrink next, if the H1 gate is promoted: the `_pqs_current_route_safe_term_matrices(...)` assertions in `test/nested/bond_aligned_diatomic_high_order_recipe_opt_in_source_construction_integration_runtests.jl`, especially the old whole-route safe-term matrix and shell-row oracle vocabulary now superseded by direct retained-boundary one-body and final H1 seam checks.
- Remaining stale/duplicate surface: `_pqs_current_route_safe_term_matrices(...)` and `_pqs_current_route_safe_term_authority_comparison(...)` remain private CCPM oracle/debug vocabulary until that caller-driven shrink pass happens.

-- repo-doer@macmini

Pass 048 response: PQS complete core-shell H1 durable gate stabilization

Summary:
- Replaced the old permanent 98-function boundary-shell-only H1 gate with a compact complete core/surrounding-shell H1 gate in `test/nested/pqs_direct_retained_final_h1_runtests.jl`.
- Updated `docs/src/developer/pqs_near_term_final_basis_realization_plan.md` so the durable H1 gate is the 223-function complete route and the old 98-function path is mechanical/oracle-only.

Files edited:
- `test/nested/pqs_direct_retained_final_h1_runtests.jl`
- `docs/src/developer/pqs_near_term_final_basis_realization_plan.md`

Permanent test result:
- Testset: `PQS complete core-shell final H1 gate`
- Route:
  - `current_box = (1:7, 1:7, 1:7)`
  - `inner_box = (2:6, 2:6, 2:6)`
  - `raw_source_dims = (5, 5, 5)`
  - `core_support_count = 125`
  - `shell_support_count = 218`
  - `shell_final_retained_count = 98`
  - `final_retained_count = 223`
  - nuclear factor source: `pgdg_intermediate.gaussian_factor_terms`
  - solve kind: `ordinary_symmetric`
- H1 value pinned by the permanent test:
  - `-0.48047934800387226`
- Same-geometry fixed-block oracle:
  - `-0.48047920531279725`
- Oracle delta:
  - `1.4269107501130307e-7`

Old boundary-shell-only pressure removed/shrunk:
- Removed the synthetic retained-unit / terminal-lowering / unit-pair / pair-materialization setup that supported the old 98-function boundary-only test.
- Removed assertions that presented the 98-function retained boundary shell as the durable H1 gate.
- The replacement test uses the complete route-owned final-basis helper:
  - `pqs_complete_core_shell_final_basis`
  - `pqs_complete_core_shell_final_one_body_matrix`
  - `pqs_complete_core_shell_final_one_electron_hamiltonian`
  - `pqs_complete_core_shell_final_h1_solve`
- The test file shrank from 375 lines to 274 lines.

Validation:
- `julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl`
  - passed: 28 checks in 2.8s
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed
- `git diff --check`
  - passed
- Docs ASCII scan:
  - `LC_ALL=C rg -n "[^\\x00-\\x7F]" docs/src/developer/pqs_near_term_final_basis_realization_plan.md`
  - no matches

Deletion/shrinkage report:
- Old/fallback/oracle surface made less necessary:
  - the 98-function boundary-shell-only H1 gate is no longer the permanent H1 coverage target.
  - the CPBM/raw-source retained boundary scaffolding in this test is no longer needed for the live H1 contract.
- Deleted/simplified:
  - deleted old synthetic retained-unit planning helpers from the H1 test;
  - deleted old pair-materialization record setup from the H1 test;
  - deleted boundary-only source CPB/raw-plan retained-rule assertions from the H1 gate;
  - relabeled the old 98-function path in docs as mechanical/source-box coverage or oracle aid, not a physical one-center acceptance basis.
- New test coverage status:
  - no new test file was added;
  - the existing focused H1 test was replaced/shrunk and now protects the live 223-function complete core/shell route.
- Remaining stale or duplicate surfaces to retire next:
  - tmp/work retained-source/PQS H1 probes or handoff text that still describe the boundary-shell-only route as physical acceptance;
  - older docs that mention the 98-function H1 path without the mechanical/nonphysical qualification;
  - any future complete H1 probe that still uses raw `gaussian_factor_matrices(base_layer)` instead of `pgdg_intermediate.gaussian_factor_terms`.

-- repo-doer@macmini

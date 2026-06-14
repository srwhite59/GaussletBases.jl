Pass 222 response - H2 221 diagnostic scaffold retirement and supplement seam audit

Task A result:
- Deleted `test/driver_inputs/h2_pqs_q5_source_box_diagnostic_r4.jl`.
- Deleted `test/nested/cartesian_ham_builder_h2_pqs_q5_source_box_diagnostic_r4_runtests.jl`.
- Did not delete or alter the 221 route implementation.

Reference audit:
- Command: `rg -n "h2_pqs_q5_source_box_diagnostic_r4" test bin src docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot`
- Before deletion, live hits were the diagnostic test pointing at its own input plus historical logs.
- After deletion, hits are only historical blurb/review/response log text under `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot`.
- No live source, `bin`, default-runner, or current physics workflow reference remains.

Supplement seam audit:

1. Current WL/MWG residual Gaussian supplement path:
- Ordinary WL/QW residual supplement machinery lives around `ordinary_cartesian_qiu_white_operators(..., supplement; interaction_treatment = :mwg)` in `src/ordinary_qw_operator_assembly.jl`, with residual centers/widths carried on `QiuWhiteResidualGaussianOperators`.
- Legacy supplement constructors live in `src/legacy_basis_adapter.jl`, including `legacy_atomic_gaussian_supplement(...)`.
- Final-residual MWG component kernels are `_qwrg_mwg_interaction_components(...)` and `_qwrg_final_residual_mwg_component_blocks(...)`.
- Route-global combined GTO/MWG preflight surfaces live in `src/cartesian_pair_block_materialization/route_global_combined_gto_*.jl`: basis layout, mixed GTO blocks, one-electron/moment matrix assembly, final-basis projection, residual MWG representation, and combined density-density readiness.

2. Required objects/matrices before the final transform:
- A common gausslet retained sector/inventory with retained dimension and ranges.
- Provider-level GTO supplement orbitals and mixed gausslet/GTO plus GTO/GTO blocks for overlap, kinetic, by-center nuclear, position, and x2 moments.
- Combined raw overlap and one-electron Hamiltonian matrices over gausslet plus raw supplement columns.
- Combined raw moment matrices `position_x/y/z` and `x2_x/y/z` to derive residual MWG centers and widths after final-basis projection.
- Final-basis residual projection transform and retained supplement count.
- For density-density, a gausslet density-density matrix, residual MWG representation, decomposed inventory/factorized retained contraction, parent axis counts/bundles, and Coulomb expansion.

3. H2 physical PQS route facts already present:
- The accepted H2 physical route has the common physical support/intermediate plan facts: support counts `(275, 578, 362)`, retained counts `(251, 98, 114)`, retained order `(:atom_contact_core, :shared_shell_1, :shared_shell_2)`, and final dimension `463`.
- It has parent axis counts/bundles, source plan, final PQS basis, H1, H1-J/density interaction, and private RHF diagnostic execution for the gausslet-only route.
- It does not yet have provider-level GTO supplement blocks, combined raw GTO matrices, raw moment matrices, a residual GTO MWG representation, or a route-owned bridge from the PQS physical final-basis payload into the existing CPBM combined-GTO supplement surfaces.

4. Smallest next implementation seam:
- Add a private H2 physical supplement preflight/payload boundary that keeps the visible driver shape and reports:
  `common physical support/intermediate gausslet plan + retained_transform_kind = :pqs + supplement_policy = :mwg_residual_gto`.
- The first implementation should be metadata/readiness only: validate the common physical gausslet plan, normalize the supplement request, identify required provider-level GTO blocks and raw moment matrices, and return explicit blockers for missing facts.
- Only after that boundary exists should matrix materialization wire into the existing CPBM combined-GTO layout/projection/MWG representation surfaces.

5. Explicitly forbidden until the seam exists:
- No supplemented H2 comparison values.
- No raw GTO/GTO or fixed/raw-GTO MWG blocks treated as final density-density data.
- No public API/export/HamV6/CR2/HFDMRG/DMRG behavior.
- No QW/Hamiltonian adoption or old fixed-block route authority.
- No IDA/MWG semantic change and no retained/source/final-residual weight division.
- No generalized-overlap final-basis transfer path.

Scoped line-count result:
- Command: `git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`
- Result: `0 added`, `68 deleted`, net `-68`.

Validation:
- `rg -n "h2_pqs_q5_source_box_diagnostic_r4" test bin src docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot`
  - passed; only historical log hits remain.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed: `load ok`.
- `git diff --check`
  - passed.
- `git diff --cached --check`
  - passed.

Git status before writing this response:
```text
## main...origin/main
 D test/driver_inputs/h2_pqs_q5_source_box_diagnostic_r4.jl
 D test/nested/cartesian_ham_builder_h2_pqs_q5_source_box_diagnostic_r4_runtests.jl
```

Deletion/shrinkage report:
- deleted: `test/driver_inputs/h2_pqs_q5_source_box_diagnostic_r4.jl`; `test/nested/cartesian_ham_builder_h2_pqs_q5_source_box_diagnostic_r4_runtests.jl`.
- simplified: removed the slow H2 221 diagnostic test pressure now that the H2 463 physical gausslet-only endpoint owns the reviewed WL/PQS comparison gate.
- quarantined: GTO/MWG supplement support remains audit-only; supplemented WL/QW scalar references remain outside the accepted H2 gausslet-only endpoint.
- not deleted because: the underlying 221 route implementation was explicitly out of scope for this pass.
- exact remaining caller/blocker: no remaining live caller for the deleted scaffold; next supplement work is blocked on a private route-owned H2 physical supplement preflight/payload boundary for `supplement_policy = :mwg_residual_gto`.

-- repo-doer@macmini

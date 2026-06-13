Pass 167 response - revised CR2 inspection artifact schema

Status: complete
Role: repo-doer@macmini
Commit: uncommitted; repo-manager owns commits.

Files changed:
- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
- deleted `test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl`

Corrections made:
- White-Lindsey route group is now schema-shaped with the same top-level families as PQS where practical:
  - `routes/white_lindsey/route`
  - `routes/white_lindsey/readiness`
  - `routes/white_lindsey/system`
  - `routes/white_lindsey/final_basis`
  - `routes/white_lindsey/one_body`
  - `routes/white_lindsey/two_body`
  - `routes/white_lindsey/validation`
- WL remains explicitly unavailable/not-applicable. No PQS pre-final/source-box data is synthesized for WL.
- Restored compact PQS route-spine semantic guards in the focused Be2 test:
  - source plan object kind is `:pqs_diatomic_complete_core_shell_source_plan`;
  - support order is `(:product, :pqs_left, :pqs_right)`;
  - route retained order is `(:pqs_left, :pqs_right, :product)`;
  - final basis dimension is `221`;
  - final support row order is `:core_then_shell`.
- Kept H1 dimension/finite/symmetry/lowest-energy checks and the artifact smoke.
- Made the artifact schema block more readable with grouped NamedTuple/Dict entries instead of over-compressed scalar key listing.

Deletion/running status:
- Deleted `test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl` under the explicit pass-167 permission.
- `rg -n "pqs_multilayer_complete_core_shell_rhf_scf_runtests" test/runtests.jl test/nested/runtests.jl test` returned no matches, so no default or nested runner included it.

Validation results:
- `julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl` passed: 25 passed, 51.3 s.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed: `load ok`.
- `git diff --check` passed.
- `git diff --numstat -- src test`:
  - `177  17  src/pqs_source_box_diatomic_complete_core_shell.jl`
  - `0    96  test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl`
  - `73   163 test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
  - total added 250, deleted 276; net-negative by 26 lines.

Git status:
`## main...origin/main`
` M src/pqs_source_box_diatomic_complete_core_shell.jl`
` D test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl`
` M test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`

Readability notes:
- Kept compact fixture argument listings where they were already local and clear.
- Expanded the schema/control block enough to show route/system/final-basis/one-body/two-body/validation group structure directly.

Deletion/shrinkage report:
- deleted: standalone private synthetic RHF SCF validator `test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl`, which no runner included.
- simplified: focused Be2 test now guards route-spine semantics and artifact contract without restoring the old assertion wall.
- quarantined: WL route remains placeholder-only in the CR2 inspection artifact; no solver/export/HamV6/HFDMRG readiness claim was added.
- not deleted because: H1 scientific checks and Be2 artifact smoke are the live contract for this pass.
- exact remaining caller/blocker: CR2 solver/export handoff remains blocked by `:missing_cr2_solver_handoff_format`; downstream solver/HamV6/HFDMRG readiness remains blocked by `:missing_hfdmrg_density_density_contract` and downstream contract objects.

-- repo-doer@macmini

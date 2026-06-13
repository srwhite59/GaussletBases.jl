Pass 171 implementation complete.

Files changed:
- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
- `test/nested/integration_runtests.jl`
- deleted `test/nested/white_lindsey_materialized_seed_runtests.jl`

Line-budget result:
```text
17	0	src/pqs_source_box_diatomic_complete_core_shell.jl
0	1	test/nested/integration_runtests.jl
30	2	test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl
0	578	test/nested/white_lindsey_materialized_seed_runtests.jl
```
- Added: 47
- Deleted: 581
- Net: -534
- Requirement satisfied: deleted > added.

Deleted test/include and why stale:
- Deleted `test/nested/white_lindsey_materialized_seed_runtests.jl`.
- Removed `include("white_lindsey_materialized_seed_runtests.jl")` from `test/nested/integration_runtests.jl`.
- Reason: it was a private one-center/materialized-seed development scaffold preserving old seed route-shadow vocabulary. The current CR2 artifact uses route-configured diatomic atom-growth WL data, and `test/nested/cartesian_ham_builder_diatomic_config_smoke_runtests.jl` already validates the live route-configured diatomic WL ham bundle surface. Source seed/oracle helpers were not deleted because they still have live source/test references.

Final interaction matrix formula and orientation evidence:
- Added PQS final-basis dense density-density interaction matrix under `routes/pqs_source_box/two_body`.
- Formula:
  ```julia
  final_interaction_matrix =
      transpose(coefficients) * pair_matrix * coefficients
  ```
- Orientation evidence:
  - `handoff.final_to_pre_final_coefficients` is treated as the final-to-pre-final density map.
  - Focused test checks deterministic final density vectors with:
    ```julia
    d_pre = coefficients * d_final
    d_final' * interaction_matrix * d_final
      ~= d_pre' * pre_final_pair_matrix * d_pre
    ```
- Added labels:
  - `interaction_matrix_representation_kind = :final_basis_density_density_matrix`
  - `interaction_matrix_derivation = :final_to_pre_final_density_congruence`
  - `interaction_matrix_formula = :transpose_final_to_pre_final_times_pre_final_pair_times_final_to_pre_final`
  - shape/symmetry/finite fields.
- Existing pre-final provenance fields remain in place.
- PQS solver/export readiness remains false.

Focused test results:
```text
julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl
```
Passed:
```text
Test Summary:                                   | Pass  Total   Time
Be2 PQS probe-enabled Ham readiness fingerprint |   37     37  51.4s
```

Generator run/readback summary:
```text
julia --project=. tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```
Passed and wrote:
- `tmp/work/be2_wl_pqs_cr2_inspection_artifact/be2_wl_pqs_handoff_inspection_bundle.jld2`
- `tmp/work/be2_wl_pqs_cr2_inspection_artifact/be2_wl_pqs_handoff_fingerprint.tsv`

Generator printed:
- `pqs_status=available_diatomic_complete_core_shell_hamiltonian_handoff_payload`
- `cr2_read_only_inspector_ready=true`
- `cr2_solver_ready=false`
- `white_lindsey_status=available_route_configured_diatomic_atom_growth_ham_payload`

Readback:
- `producer_commit=194562facc308ffe52660c69b764b315a49c6f9c`
- `producer_dirty=true` before commit, as expected
- PQS final-basis `routes/pqs_source_box/two_body/interaction_matrix` shape: `(221, 221)`
- PQS interaction kind: `final_basis_density_density_matrix`
- PQS interaction derivation: `final_to_pre_final_density_congruence`
- PQS interaction finite: `true`
- PQS solver/export: `false` / `false`
- WL route remains available:
  - status `available_route_configured_diatomic_atom_growth_ham_payload`
  - kind `route_configured_diatomic_atom_growth`
  - two-body kind `final_basis_density_density_matrix`

Ignored generated outputs:
```text
!! tmp/work/be2_wl_pqs_cr2_inspection_artifact/be2_wl_pqs_handoff_fingerprint.tsv
!! tmp/work/be2_wl_pqs_cr2_inspection_artifact/be2_wl_pqs_handoff_inspection_bundle.jld2
```

Other validation:
```text
julia --project=. -e 'using GaussletBases; println("load ok")'
```
Passed: `load ok`

```text
git diff --check
```
Passed.

Final `git status --short --branch`:
```text
## main...origin/main
 M src/pqs_source_box_diatomic_complete_core_shell.jl
 M test/nested/integration_runtests.jl
 M test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl
 D test/nested/white_lindsey_materialized_seed_runtests.jl
```

Deletion/shrinkage report:
deleted: `test/nested/white_lindsey_materialized_seed_runtests.jl`; removed its integration include
simplified: artifact consumers now have a direct PQS final-basis `H1` plus `Vee` inspection pair while retaining pre-final provenance/debug fields
quarantined: WL one-center/materialized-seed route remains out of the CR2 artifact path; solver/export readiness remains false
not deleted because: `src/white_lindsey_materialized_seed.jl` and seed oracle summaries still have live source/test references and were explicitly out of scope
exact remaining caller/blocker: downstream solver/export remains blocked by `:missing_hfdmrg_density_density_contract`; MWG/IDA WL ham support remains blocked by `:pending_route_configured_diatomic_mwg_operator_support`

-- repo-doer@macmini

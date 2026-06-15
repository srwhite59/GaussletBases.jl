Pass 256 response - independent H2 PQS supplement staging audit

Summary:
- No source, test, bin, or docs edits.
- Existing supplement request/representation/preflight helpers can mechanically attach to the independent H2 PQS target because `_pqs_source_box_route_driver_diatomic_physical_gausslet_target_payload` accepts `:bond_aligned_diatomic_independent_pqs_source_box_core_shell` and normalizes the generated support/retained plans to `:available_physical_gausslet_core_shell_target_inventory`.
- The seam is not clean enough for provider-block work yet: helper/type/status names and fixture labels are still `physical_gausslet`/old-target shaped, and the fake-PQS preflight test remains schema/history only.

Direct audit answers:

1. Existing payload attachment:
   - `_pqs_source_box_route_driver_diatomic_physical_gausslet_supplement_request_payload(parent, target_payload)` is target-driven, not fake-source-backed-driven. For an independent target with `supplement_policy = :mwg_residual_gto`, it would produce an available request, `H/cc-pVTZ`, `lmax = 1`, and missing provider-block facts.
   - `_pqs_source_box_route_driver_diatomic_physical_gausslet_supplement_representation_payload(request_payload)` is route-neutral after the request exists. It constructs the H/cc-pVTZ diatomic supplement representation from parent nuclei and should still produce 18 orbitals for H2/lmax=1.
   - `_pqs_source_box_route_driver_diatomic_physical_gausslet_supplement_preflight_payload(target_payload, request_payload, representation_payload)` is target-driven and would use independent target counts/dimension if the independent input requested `:mwg_residual_gto`.
   - However, these payloads are still effectively old physical-target shaped in vocabulary: function/type/status names, blockers such as `:missing_physical_gausslet_target_inventory`, fixture label `:h2_r4_physical_gausslet_q5`, and fact labels such as `:physical_gausslet_core_shell_target_inventory`. They are not fake-PQS shaped in data dependency, but they are not yet clean independent-PQS supplement names either.

2. First independent supplement input:
   - Recommended later input name: `test/driver_inputs/h2_pqs_q5_independent_source_box_r4_supplement_preflight.jl`.
   - Base it on `h2_pqs_q5_independent_source_box_r4_h1_j.jl`, not the private-RHF input. The supplement preflight needs the fake-free independent final basis and H1-J density-interaction seam as route context, but should not run private RHF or supplemented physics.
   - Set `artifact_role = :independent_h2_pqs_supplement_preflight_diagnostic`, `supplement_policy = :mwg_residual_gto`, `run_final_basis = true`, `run_h1 = true`, `run_h1_j = true`, `run_private_rhf = false`, and keep all supplemented-value flags off.

3. Artifact facts the first independent supplement preflight should expose:
   - `fake_pqs/enabled = false`
   - `route/source_backed_fixed_source_oracle_used = false`
   - `route/retained_transform_authority = :pqs_source_box_construction`
   - support counts `(275, 578, 362)`
   - retained counts `(275, 98, 98)`
   - gausslet final dimension `471`
   - supplement request/representation basis `"H/cc-pVTZ"`
   - supplement `lmax = 1`
   - supplement representation orbital count `18`
   - preflight status blocked on provider work, with current honest top-level blocker `:missing_provider_gto_supplement_blocks`
   - preflight missing facts should include at least `:missing_mixed_gausslet_gto_blocks`, `:missing_gto_gto_blocks`, `:missing_combined_raw_moment_matrices`, `:missing_residual_mwg_representation`, and `:missing_combined_density_density_readiness`.

4. Existing provider/kernel objects:
   - Reusable as provider/kernel references: `CartesianCPBBlockProviders.cpb_mixed_gto_overlap_block`, `cpb_mixed_gto_position_operator_block`, `cpb_mixed_gto_x2_operator_block`, `cpb_mixed_gto_kinetic_operator_block`, and whole-supplement variants covered by `cartesian_cpb_mixed_gto_overlap_block_runtests.jl`.
   - Reusable as GTO/GTO one-body providers: `cpb_gto_overlap_operator_block`, `cpb_gto_position_operator_block`, `cpb_gto_x2_operator_block`, and `cpb_gto_kinetic_operator_block`, covered by `cartesian_cpb_gto_supplement_one_body_block_runtests.jl`.
   - Reusable as by-center nuclear provider references: `cpb_mixed_gto_nuclear_by_center_block`, `cpb_gto_nuclear_by_center_block`, and `cpb_gto_supplement_local_operator_bundle`, covered by `cartesian_cpb_gto_nuclear_by_center_block_runtests.jl`.
   - Partly reusable as readiness/reference vocabulary: residual MWG representation and combined density-density readiness in `cartesian_combined_gto_density_density_readiness_runtests.jl`.
   - Still missing for the independent route: a route-owned provider-block payload that maps the independent final gausslet basis/CPB rows into mixed and GTO/GTO provider blocks; route-global combined one-electron matrices; route-global raw moment matrices for residual MWG; residual MWG representation from independent final-basis moments; combined density-density readiness that does not depend on fake-PQS or WL/QW scalar references.

5. Must not copy from fake-PQS/WL preflight:
   - 463 final dimension.
   - Retained counts `(251, 98, 114)`.
   - `:source_backed_fixed_source_oracle` source labels.
   - `:fake_pqs_private_source_backed_adapter_authority`.
   - fake-PQS comparison role or WL/QW reproduction labels.
   - old supplemented WL/QW scalar reference comparisons.
   - source-backed retained-transform authority or any implication that fake-PQS preflight proves independent-PQS provider readiness.

6. Smallest safe next pass:
   - Do input/preflight only.
   - Add an explicit independent supplement-preflight artifact role to the compact independent-H2 role classifier if needed.
   - Add the small input file named above.
   - If a focused driver artifact test is added, it should assert only preflight facts and blocker facts, not provider matrices or supplemented energies.
   - Do not start provider-block object implementation in the next pass. The current names/status vocabulary should be quarantined or renamed toward independent supplement preflight before materialization code is attached.

Recommended next-pass boundaries:
- Allowed:
  - `test/driver_inputs/h2_pqs_q5_independent_source_box_r4_supplement_preflight.jl`
  - the compact independent-H2 artifact-role classifier, only if required for the new artifact role
  - a narrow artifact preflight test, only if it replaces or avoids copying fake-PQS assertions
- Forbidden:
  - provider block materialization
  - mixed gausslet/GTO matrix assembly
  - residual MWG representation
  - supplemented H1/H1-J/RHF values
  - supplemented WL/QW scalar references
  - CR2/export/public API/HamV6 readiness
  - fake-PQS evidence as independent-PQS evidence

Validation:
- `git diff --check` passed.
- `git status --short --branch` before response: `## main...origin/main`.
- No Julia command run; this was a no-edit inspection pass and the blurb prohibited slow H2/broad gates.

Scoped line count:
- `src + test + bin`: 0.

Deletion/shrinkage result:
- deleted: none.
- simplified: none.
- quarantined: fake-PQS preflight remains schema/history only; not used as independent-PQS evidence.
- not deleted because: this was an audit-only pass, and the supplement helpers are still live staging surfaces for the next input/preflight seam.
- exact remaining caller/blocker: `cartesian_assembly` still builds the supplement request/representation/preflight payloads for the target payload; first blocker is missing independent supplement preflight input/artifact role, then `:missing_provider_gto_supplement_blocks` for actual matrix/provider work.

-- repo-doer@macmini

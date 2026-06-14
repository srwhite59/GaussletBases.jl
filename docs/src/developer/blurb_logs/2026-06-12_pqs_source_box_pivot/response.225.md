Pass 225 response - audit H2 PQS MWG/GTO supplement provider seam

Scope:
- No source/test/bin/doc edits were made besides this response and the curated response copy.
- No Julia validation was run; this was a read-only implementation-readiness audit.

Read-only commands run:
- `sed -n '1,220p' AGENTS.md`
- `sed -n '1,220p' .agent_handoffs/RUN.md`
- `sed -n '1,220p' .agent_handoffs/DOER_STARTUP.md`
- `sed -n '1,220p' .agent_handoffs/state.md`
- `sed -n '1,260p' .agent_handoffs/blurb.225.md`
- `sed -n '1,260p' docs/src/developer/pqs_source_box_operator_framework.md`
- `git log --oneline -5`
- `git status --short --branch`
- `rg` searches over `residual`, `MWG`, `GTO`, `supplement`, `cc-pVTZ`, `mixed`, `gaussian`, `basis_representation`, and route/provider function names in `src`, `test`, `docs/src/developer`, and `tmp/work` while ignoring `*.tmp`. One broad search included a missing `src/backends` path and exited 2 after printing useful matches; no files were changed.
- Targeted `sed -n` reads of the files listed below.

Files/functions audited:
- `src/CartesianCPBBlockProviders.jl`
  - `cpb_mixed_gto_overlap_block`
  - `cpb_gto_overlap_operator_block`
  - `cpb_gto_nuclear_by_center_block`
  - `cpb_gto_supplement_local_operator_bundle`
- `src/cartesian_gto_probes.jl`
  - `basis_representation` conversion for legacy supplements
  - `gto_overlap_matrix`
  - `_cartesian_final_gto_cross_overlap_handoff`
  - `_pqs_source_box_gto_cross_overlap_shadow`
- `src/cartesian_basis_representation.jl`
  - `CartesianGaussianShellSupplementRepresentation3D`
- `src/ordinary_qw_residuals.jl`
  - `_qwrg_residual_space_analysis`
  - `_qwrg_residual_space`
  - `_qwrg_residual_space_by_owner`
  - `_qwrg_residual_center_data`
  - `_qwrg_residual_moment_data`
  - `_qwrg_final_residual_mwg_component_blocks`
- `src/ordinary_cartesian_ida.jl`
  - `ordinary_cartesian_ida_operators(...; interaction_treatment = :residual_gaussian_mwg)`
  - `_hybrid_residual_gaussian_pair_factors_mwg`
- `src/cartesian_pair_block_materialization/route_global_combined_gto_basis_layout.jl`
  - `route_global_combined_gto_basis_layout`
- `src/cartesian_pair_block_materialization/route_global_mixed_gto_blocks.jl`
  - `route_global_mixed_gto_blocks_from_decomposed_units`
- `src/cartesian_pair_block_materialization/route_global_combined_gto_matrix_assembly.jl`
  - `route_global_combined_gto_one_electron_matrices`
  - `route_global_combined_gto_residual_moment_matrices`
- `src/cartesian_pair_block_materialization/route_global_combined_gto_final_basis.jl`
  - `route_global_combined_gto_final_basis_projection`
- `src/cartesian_pair_block_materialization/route_global_combined_gto_density_density.jl`
  - `route_global_residual_gto_mwg_representation`
  - `route_global_combined_gto_final_basis_density_density_matrix`
- `src/cartesian_pair_block_materialization/route_global_atom_gto_final_basis_route.jl`
  - `_white_lindsey_decomposed_atom_gto_final_basis_route`
- `src/pqs_source_box_diatomic_complete_core_shell.jl`
  - `_pqs_source_box_route_driver_diatomic_physical_gausslet_target_payload`
  - `_pqs_source_box_route_driver_diatomic_physical_gausslet_supplement_preflight_payload`
  - `_pqs_source_box_route_driver_diatomic_physical_gausslet_source_plan_payload`
  - `_pqs_source_box_route_driver_diatomic_physical_gausslet_final_basis_payload`
  - physical H1/H1-J/RHF payload helpers
- `src/pqs_source_box_route_driver_helpers.jl`
  - `cartesian_assembly` physical target/supplement preflight wiring
- `src/pqs_source_box_route_driver_reporting.jl`
  - H2 physical artifact `supplement_preflight/*` reporting
- Tests audited:
  - `test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl`
  - `test/nested/cartesian_cpb_mixed_gto_overlap_block_runtests.jl`
  - `test/nested/cartesian_cpb_gto_supplement_one_body_block_runtests.jl`
  - `test/nested/cartesian_cpb_gto_supplement_local_operator_bundle_runtests.jl`
  - `test/nested/cartesian_mixed_gto_whole_supplement_source_runtests.jl`
  - `test/nested/cartesian_route_global_combined_gto_layout_runtests.jl`
  - `test/nested/cartesian_combined_gto_density_density_readiness_runtests.jl`
  - `test/ordinary/mwg_residual_component_helper_runtests.jl`
  - relevant old ordinary/diatomic QW supplement sections in `test/ordinary/runtests.jl` and `test/diatomic/runtests.jl`

Existing WL/MWG/GTO supplement provider surfaces:
- Provider-local blocks live in `CartesianCPBBlockProviders`:
  - CPB/GTO mixed overlap, kinetic, coordinate moments, and by-center nuclear blocks.
  - GTO/GTO overlap, kinetic, coordinate moments, and by-center nuclear blocks.
  - `cpb_gto_supplement_local_operator_bundle` packages the mixed and GTO/GTO provider blocks with uncharged by-center nuclear blocks.
- Supplement request/source conversion lives in the legacy supplement constructors plus `basis_representation`, producing `CartesianGaussianShellSupplementRepresentation3D` with explicit orbitals, centers, angular powers, primitive exponents/coefficients, basis metadata, and nuclei.
- Route-global combined-GTO machinery lives in `CartesianPairBlockMaterialization`:
  - `route_global_combined_gto_basis_layout` defines gausslet/GTO block ranges and provider source names.
  - `route_global_mixed_gto_blocks_from_decomposed_units` consumes decomposed WL retained units and provider blocks to produce retained-row mixed blocks plus GTO/GTO blocks.
  - `route_global_combined_gto_one_electron_matrices` assembles combined overlap/H1 from gausslet route-global matrices plus provider mixed and GTO/GTO blocks.
  - `route_global_combined_gto_final_basis_projection` residualizes GTOs against an orthonormal gausslet final sector.
  - `route_global_residual_gto_mwg_representation` extracts residual MWG centers/widths from combined overlap/position/x2 moments.
  - `route_global_combined_gto_final_basis_density_density_matrix` blocks until residual MWG representation and density-density components exist; it explicitly does not accept raw GTO density-density as the final operator.
- Old QW/MWG machinery in `ordinary_qw_residuals.jl` and `ordinary_cartesian_ida.jl` owns residual orthogonalization, residual owner/moment extraction, and final residual MWG component construction. It is useful as a kernel/convention source, but not as PQS route authority.

Reusable versus WL-specific:
- Route-independent/reusable:
  - `CartesianGaussianShellSupplementRepresentation3D` and legacy supplement-to-representation conversion.
  - Provider-local CPB/GTO and GTO/GTO blocks, provided the PQS route supplies the correct CPBs/support rows and center records.
  - Combined layout/matrix/projection helpers conceptually: they only need gausslet retained dimension, gausslet-sector matrices, provider bundle blocks, and ranges. Their names and some status labels say decomposed WL, but the algebra is not inherently WL-specific.
  - Residual MWG moment extraction and component kernels after a final residual supplement basis exists.
- WL-specific or not route authority for PQS:
  - `route_global_mixed_gto_blocks_from_decomposed_units` as currently written requires `:available_white_lindsey_decomposed_unit_pair_inventory`, `retained_units`, `source_cpbs`, WL unit coefficient maps, and for the factorized fast path currently requires one origin-centered center. That is not directly usable for H2 PQS without an adapter.
  - `_white_lindsey_decomposed_atom_gto_final_basis_route` is a one-center WL stitching seam, not the H2 PQS route.
  - Old ordinary/nested QW supplemented acceptance paths are comparison/oracle history. They should not become PQS construction authority.
  - Raw GTO MWG blocks are explicitly residual-construction/debug inputs only, not final two-body route payloads.

Exact PQS H2 physical data already available:
- Geometry/centers:
  - H2 R=4 center match at `(0.0, 0.0, -2.0)` and `(0.0, 0.0, 2.0)`.
  - two nuclear charges, both `1.0`.
  - parent/system center table through the route parent and by-center H1 payload.
- Parent/axis facts:
  - parent axis counts `(9, 9, 15)`.
  - parent axis bundle object available through the physical source plan (`source_plan.axis_bundles`).
  - parent basis available through `source_plan.parent_basis`.
- Physical target inventory:
  - support order/counts: `(:atom_contact_core, :shared_shell_1, :shared_shell_2)` with `(275, 578, 362)`.
  - retained order/counts: `(:atom_contact_core, :shared_shell_1, :shared_shell_2)` with `(251, 98, 114)`.
  - expected/actual gausslet final dimension `463`.
- Source/final route payloads:
  - source support indices and states for atom-contact core plus shared shell groups.
  - core coefficient matrix, shared-shell coefficient matrices, retained ranges, final dimension.
  - final basis, pre-final coefficients, final coefficients, cleanup/transform data, and final overlap identity diagnostics.
- Operator data already built for the no-supplement endpoint:
  - support kinetic and by-center nuclear blocks.
  - final kinetic and by-center nuclear blocks.
  - final H1/H1 solve.
  - H1-J diagnostic/private density path, pre-final density interaction, density gauge `:pre_final_localized_positive_weight`, and raw pair factor convention `:raw_numerator`.
  - private RHF input/execution payloads for the gausslet-only endpoint.

Missing facts before `:provider_gto_supplement_blocks` can become available:
- A route-owned supplement request payload for this exact H2 case:
  - `supplement_policy = :mwg_residual_gto`.
  - atom symbols/charges/locations bound to the route geometry.
  - basis family/name, likely `H/cc-pVTZ` unless manager selects otherwise.
  - `lmax`, contraction/uncontracted policy, residual keep/drop tolerances, and whether owner-local residualization is required.
  - supplement provenance and fixture label.
- Conversion from that request to a `CartesianGaussianShellSupplementRepresentation3D` with two H-centered orbital sets at z = `-2.0` and `2.0`.
- A PQS provider-block adapter that maps the PQS final/pre-final gausslet sector to provider inputs:
  - CPB/support source for each support group or a factorized projection path that does not require WL retained units.
  - mixed gausslet/GTO block orientation `gausslet_rows_by_gto_columns`.
  - full row coverage for all 463 PQS final gausslet rows or a documented pre-final-to-final projection convention.
  - GTO/GTO overlap, kinetic, coordinate moments, and uncharged by-center nuclear blocks.
  - center-count agreement between route by-center nuclear records and provider by-center blocks.
- A PQS-facing combined layout wrapper:
  - current layout status/kinds say decomposed WL; a PQS wrapper should identify the gausslet sector as `:pqs_physical_gausslet_final_basis` while reusing the same range/orientation algebra.
- Raw combined moment matrices if proceeding past one-electron assembly:
  - position/x2 gausslet-sector matrices and mixed/GTO moment blocks in combined raw basis order.
- Residual MWG representation and final density-density readiness:
  - residual centers/widths from combined moments.
  - residual MWG density-density component construction using PQS contraction data rather than decomposed WL inventory.

Smallest recommended implementation pass:
- First pass should be a compact supplement request/preflight payload, not provider-block materialization.
- Suggested internal name:
  - `_PQSDiatomicPhysicalGaussletSupplementRequestPayload`
- Suggested fields:
  - `status`
  - `blocker`
  - `route_family`
  - `route_kind`
  - `fixture_label`
  - `supplement_policy`
  - `atom_symbols`
  - `nuclear_charges`
  - `atom_locations`
  - `bond_axis`
  - `bond_length`
  - `basis_name`
  - `lmax`
  - `uncontracted`
  - `residual_keep_policy`
  - `residual_drop_tolerance`
  - `representation`
  - `representation_status`
  - `required_provider_blocks`
  - `available_fact_labels`
  - `missing_fact_labels`
  - `summary`
  - `metadata`
- Behavior should remain blocked and matrix-free:
  - if request/repr is not implemented, sharpen the blocker to the request-level missing fact;
  - if request/repr is implemented, keep the next blocker as `:missing_provider_gto_supplement_blocks`.
- This is the smallest useful seam because it binds policy, geometry, and supplement identity before choosing a provider adapter. Without it, a provider-block adapter would have to infer chemistry and basis choice from scattered recipe/system fields.

Candidate deletion/shrink surfaces for the first coding pass:
- deleted:
  - none in this no-edit pass.
- simplified:
  - After the request payload exists, shrink `supplement_preflight` report construction so it reads compact request/preflight summaries instead of carrying repeated flat support/retained/request facts.
  - In `test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl`, replace repeated supplement preflight field assertions with a compact request/preflight fingerprint plus the endpoint-critical fields.
  - In route-global combined GTO layout/readiness tests, replace exhaustive metadata flag assertions with compact status/blocker/range/orientation checks where the compact payload now owns the invariant.
- quarantined:
  - Old ordinary/nested QW supplemented scalar paths should remain comparison/oracle history until a driver-owned supplemented H2 endpoint exists. They should not be extended for PQS.
  - `tmp/work` Be2 supplement reports are useful historical audits only; do not promote them as live H2 route tests.
- not deleted because:
  - Provider-local CPB/GTO and GTO/GTO tests still protect the analytic block conventions needed by the provider seam.
  - Residual MWG component helper test still protects the key convention that raw GTO rows are residual-construction inputs only.
  - H2 physical endpoint/preflight test is the live driver artifact contract.
- exact remaining caller/blocker:
  - Current preflight blocker is `:missing_provider_gto_supplement_blocks`.
  - The more precise first missing implementation fact is route-owned H2 supplement request/representation, then a PQS adapter that provides mixed gausslet/GTO and GTO/GTO blocks in the route-global combined-GTO orientation.

Git status:

```text
## main...origin/main
```

-- repo-doer@macmini

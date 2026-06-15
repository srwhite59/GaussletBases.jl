Pass 259 response - independent H2 PQS supplement provider-block seam audit

Summary:
- No source, test, bin, docs, or driver-input edits.
- The first provider-block step should be a private route-owned local-provider payload, not route-global supplemented operators.
- Key seam hazard: existing CPB provider functions operate on rectangular CPBs. The independent H2 shared-shell support is outer-box minus inner-box, and `source_cpb` is the filled outer source box. A provider-block implementation must therefore carry explicit support CPB tiling/row ownership or an equivalent row map before calling provider blocks. It must not blindly call provider blocks on the filled shared-shell `source_cpb`.

1. First route-owned payload:
- Recommended private payload name:
  - struct: `_PQSIndependentH2PQSSupplementProviderBlocksPayload`
  - helper: `_pqs_source_box_route_driver_independent_h2_pqs_supplement_provider_blocks_payload`
- It should live beside the current supplement request/representation/preflight helpers in `src/pqs_source_box_diatomic_complete_core_shell.jl`.
- Keep it private and route-owned. Do not export it and do not make it a public supplement API.
- It may wrap/reuse existing `physical_gausslet` request/representation payloads, but the new provider payload name should say independent-H2-PQS to avoid extending the old physical/fake vocabulary.

2. Required inputs:
- `parent`: the `CartesianParentGaussletBasis3D` and axis/factor data used by CPB providers.
- `target_payload`: route kind, support counts `(275, 578, 362)`, retained counts `(275, 98, 98)`, final dimension `471`, independent authority labels.
- `source_plan_payload`: source plan status, support order, support states/indices, sparse coefficient maps, parent basis/axis bundles, and retained ranges.
- `final_basis_payload`: final-basis status and final/pre-final coefficient maps. Needed for later projection/fingerprints and to prove the provider payload is attached to the independent final basis, even if first-pass matrices stay local/provider-level.
- `supplement_request_payload`: policy `:mwg_residual_gto`, basis `"H/cc-pVTZ"`, `lmax = 1`, geometry and required provider facts.
- `supplement_representation_payload`: the actual supplement and representation, with 18 orbitals.
- Parent Coulomb/nuclear expansion data for by-center GTO nuclear blocks. Existing code commonly uses `coulomb_gaussian_expansion(doacc = false)`.
- Center records with center key/index, nuclear charge, and location for both H nuclei.
- A route-owned CPB support partition or equivalent row-ownership map:
  - atom-contact core should be represented as its actual atom-local-core plus midpoint-slab pieces, not just a flat support vector;
  - shared shells need rectangular support tiles or row filters for `outer_box minus inner_box`;
  - every tile needs a mapping back to parent rows and route support-unit rows.

3. Existing provider functions to call/wrap:
- Prefer wrapping the existing bundle first:
  - `CartesianCPBBlockProviders.cpb_gto_supplement_local_operator_bundle(parent, cpb, supplement; expansion, center_records)`
- That bundle already calls the local providers for:
  - mixed gausslet/GTO overlap, position x/y/z, x2 x/y/z, kinetic;
  - GTO/GTO supplement self overlap, position x/y/z, x2 x/y/z, kinetic;
  - mixed by-center nuclear blocks;
  - GTO/GTO by-center nuclear blocks.
- If the implementation needs finer control, the direct functions are:
  - `cpb_mixed_gto_overlap_block`
  - `cpb_mixed_gto_position_operator_block`
  - `cpb_mixed_gto_x2_operator_block`
  - `cpb_mixed_gto_kinetic_operator_block`
  - `cpb_gto_overlap_operator_block`
  - `cpb_gto_position_operator_block`
  - `cpb_gto_x2_operator_block`
  - `cpb_gto_kinetic_operator_block`
  - `cpb_mixed_gto_nuclear_by_center_block`
  - `cpb_gto_nuclear_by_center_block`

4. First-pass payload products:
- Produce a compact provider-block payload with:
  - status/blocker;
  - target/source/final/supplement status fingerprints;
  - support partition coverage summary: unit keys, tile count, parent-row coverage count, duplicate count, missing count, outside count;
  - provider bundle summaries per tile, not report-level scalar clouds;
  - term availability counts for mixed one-body, GTO self one-body, mixed nuclear by-center, and GTO self nuclear by-center;
  - local dense block handles retained inside provider block/bundle objects only;
  - compact block fingerprints such as term, status, blocker, shape, finite flag, support count, orbital count, center count.
- Matrices should remain local/provider-level in the first implementation pass.
- Do not assemble route-global mixed/GTO matrices.
- Do not project through final coefficients.
- Do not build combined raw moment matrices or residual MWG representation.
- Artifact/reporting should expose only status/count/fingerprint summaries, not dense matrices or broad scalar aliases.

5. Artifact fields/status after provider blocks exist:
- `route/supplement_preflight_blocker` should advance from `:missing_provider_gto_supplement_blocks` to the next honest blocker, likely `:missing_combined_raw_moment_matrices`.
- `supplement_preflight/available_fact_labels` should add:
  - `:provider_gto_supplement_blocks`
  - `:mixed_gausslet_gto_blocks`
  - `:gto_gto_blocks`
- `supplement_preflight/missing_fact_labels` should remove provider/mixed/GTO block blockers and retain:
  - `:missing_combined_raw_moment_matrices`
  - `:missing_residual_mwg_representation`
  - `:missing_combined_density_density_readiness`
- Add a compact artifact group only if needed, e.g. `supplement_provider_blocks`, with:
  - `status`
  - `blocker`
  - `object_kind`
  - `support_partition_status`
  - `support_tile_count`
  - `support_coverage_complete`
  - `mixed_one_body_term_count`
  - `gto_one_body_term_count`
  - `mixed_by_center_nuclear_count`
  - `gto_by_center_nuclear_count`
  - `provider_level_local_blocks_materialized`
  - `route_global_matrix_materialized = false`
  - `supplemented_values_materialized = false`
  - `public_api = false`
- Keep `fake_pqs/enabled = false`, `route/source_backed_fixed_source_oracle_used = false`, and `route/retained_transform_authority = :pqs_source_box_construction`.

6. Must remain blocked after provider blocks:
- combined raw moment matrices;
- residual MWG representation;
- combined density-density readiness;
- supplemented H1/H1-J/RHF values;
- CR2/export/HamV6/public API readiness;
- WL/QW supplemented scalar comparisons;
- fake-PQS/source-backed evidence.

7. Smallest validation fixture:
- Do not copy the large CPB provider tests.
- Use a compact route-owned smoke that exercises:
  - independent H2 support-region/retained-rule facts;
  - H/cc-pVTZ lmax-1 supplement representation;
  - one atom-contact tile and one shared-shell support tile or partition fixture;
  - the provider bundle summary, with local matrix handles present but no route-global matrices.
- The smoke should assert:
  - fake-free independent route authority;
  - support tile row coverage/fingerprint;
  - mixed/GTO one-body and by-center nuclear term counts/statuses;
  - `route_global_matrix_materialized == false`;
  - next blocker remains `:missing_combined_raw_moment_matrices`.
- Avoid the full slow H2 preflight route unless the report/artifact plumbing itself is changed and cannot be validated by a smaller smoke.

Recommended exact next implementation pass:
- Implement only the private independent-H2 provider-block payload and a compact support-partition summary.
- Use `cpb_gto_supplement_local_operator_bundle` for local provider-level block handles once CPB tile ownership is explicit.
- Wire it after supplement representation/preflight in `cartesian_assembly`, gated to the independent-H2 supplement-preflight diagnostic role and `supplement_policy = :mwg_residual_gto`.
- Update preflight/report summaries only enough to mark provider blocks available and advance the blocker to combined raw moment matrices.
- Add a small route-owned smoke for the provider payload. Do not run or add supplemented-value endpoint tests.

Validation:
- `git diff --check` passed.
- No Julia command was run; this was a no-edit audit pass.

Scoped line count:
- `src + test + bin`: 0.

Git status:
```text
## main...origin/main
```

Deletion/shrinkage result:
- deleted: none.
- simplified: none.
- quarantined: provider blocks remain unimplemented; route-global matrices, residual MWG, combined density readiness, supplemented values, CR2/export/HamV6/public API remain blocked.
- not deleted because: audit-only pass; no stale mirror deletion was needed to offset code.
- exact remaining caller/blocker: `cartesian_assembly` currently builds supplement request/representation/preflight metadata only; the next implementation blocker is explicit support CPB tiling/row ownership for independent H2 provider-local blocks, followed by `:missing_combined_raw_moment_matrices`.

-- repo-doer@macmini

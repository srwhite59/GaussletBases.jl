Pass 260 response - independent H2 PQS supplement support-tiling audit

Summary:
- No source, test, bin, docs, or driver-input edits.
- A support-partition summary is ready to implement before provider blocks.
- The next implementation should expose route-owned support tiling/row ownership first, then provider blocks in a later pass.

1. Where support rows/unit keys/source boxes/parent row indices are available:
- `_pqs_source_box_route_driver_independent_h2_support_region_plan(parent, route_skeleton, recipe)` constructs the independent support facts.
- Returned plan fields already include:
  - `support_order = (:atom_contact_core, :shared_shell_1, :shared_shell_2)`
  - `support_counts = (atom_contact_core = 275, shared_shell_1 = 578, shared_shell_2 = 362)`
  - route coverage facts: `coverage_complete`, `duplicate_count`, `missing_count`, `outside_count`
  - `atom_contact_core_descriptor.support_indices`
  - `atom_contact_core_descriptor.support_states`
  - `shared_shell_descriptors.shared_shell_1.support_indices`
  - `shared_shell_descriptors.shared_shell_2.support_indices`
  - shared shell `current_box`, `inner_box`, and filled `source_cpb`
- `_PQSDiatomicPhysicalGaussletCoreShellSourcePlan` carries the source-plan row data after materialization:
  - `atom_contact_core_support_indices`
  - `atom_contact_core_support_states`
  - `shared_shell_support_indices`
  - `shared_shell_support_states`
  - `support_order`
  - `retained_order`
  - `retained_ranges`
  - sparse parent-row coefficient maps.
- Final-basis payloads add:
  - `support_row_ranges`
  - `retained_ranges`
  - `pre_final_coefficients`
  - `final_coefficients`
  - final dimension and overlap diagnostics.

2. Atom-contact core decomposition:
- Conceptually yes: the current construction builds atom-contact support from `atom_regions..., midpoint_regions...`.
- The code computes:
  - `atom_regions = regions with role == :atom_local_core`
  - `midpoint_regions = regions with role == :midpoint_slab`
  - `atom_contact_regions = (atom_regions..., midpoint_regions...)`
  - `atom_contact_core_support_indices = reduce(vcat, atom_contact_support(region) for region in atom_contact_regions)`
- Expected conceptual counts are two atom-local cores plus one midpoint/contact slab, totaling `275`. For the current H2 q5 geometry this is consistent with `125 + 125 + 25`.
- Current blocker: the returned `atom_contact_core_descriptor` preserves only aggregate `source_region_roles`, aggregate support indices/states/count, and direct-selector coefficient representation. It does not preserve per-piece boxes, per-piece CPBs, per-piece support counts, or per-piece row ranges.
- Therefore atom-contact provider tiling is implementable, but first the support-partition payload should retain per-region atom/contact tile descriptors instead of reconstructing them later from a flattened support vector.

3. Shared-shell decomposition:
- The shared-shell descriptors already expose the critical facts:
  - `role`
  - `terminal_region_key`
  - `current_box` from `raw.outer_box`
  - `inner_box` from `raw.inner_exclusion_box`
  - filled `source_cpb`
  - actual `support_indices`
  - actual `support_states`
  - `support_count`
- This is enough to avoid treating filled `source_cpb` as support. The support is `current_box \ inner_box`, not `source_cpb`.
- Existing lower-level utilities already encode the shell-difference idea:
  - `_nested_diatomic_shell_support_indices(outer_box, inner_box, dims)` gives actual parent rows for `outer_box \ inner_box`.
  - `CartesianRouteCore.complete_shell_support(outer, inner)` represents shell-owned support as `OwnedSupport`, not a CPB.
  - `CartesianCPB.complete_shell_boundary_strata(outer, inner)` provides a CPB stratum tiling for one-layer complete shells.
- Caveat: `complete_shell_boundary_strata` has a one-layer complete-shell invariant. The next pass must verify whether each independent shared shell's `current_box`/`inner_box` satisfies that invariant. If not, implement a narrow rectangular-difference tiler for the existing boxes and validate against `support_indices`.

4. First support-partition row-map/fingerprint:
- Recommended payload fields:
  - `status`
  - `blocker`
  - `route_family`
  - `route_kind`
  - `support_order`
  - `retained_order`
  - `support_counts`
  - `retained_counts`
  - `retained_ranges`
  - `unit_partitions`
  - `total_tile_count`
  - `total_support_count`
  - `coverage_complete`
  - `duplicate_parent_row_count`
  - `missing_parent_row_count`
  - `outside_parent_row_count`
  - `source_backed_fixed_source_oracle_used = false`
  - `fake_pqs_enabled = false`
  - `retained_transform_authority = :pqs_source_box_construction`
- Per-unit partition fingerprint:
  - `unit_key`
  - `support_count`
  - `retained_range`
  - `tile_count`
  - `tile_support_counts`
  - `parent_row_count`
  - `parent_row_min`
  - `parent_row_max`
  - `duplicate_parent_row_count`
  - `missing_parent_row_count`
  - `outside_parent_row_count`
  - `coverage_complete`
- Per-tile fingerprint:
  - `tile_key`
  - `unit_key`
  - `source_region_role`
  - `cpb_role`
  - `intervals`
  - `support_count`
  - `parent_row_ids` or a compact parent-row fingerprint plus optional row ids for the route-owned payload
  - `unit_local_row_range` or explicit `unit_local_rows`
  - `retained_range`
  - `provider_tile_ready`

5. Is a support-partition summary implementable before provider blocks?
- Yes.
- Recommended private object/helper:
  - struct: `_PQSIndependentH2PQSSupplementSupportPartitionPayload`
  - helper: `_pqs_source_box_route_driver_independent_h2_pqs_supplement_support_partition_payload`
- Inputs should be `target_payload`, `source_plan_payload`, and the support-region plan or enough parent/route data to recover it.
- The next pass should add no provider block calls. It should only materialize tile descriptors and coverage/fingerprint summaries.
- Exact missing data to add or preserve:
  - atom-contact per-region boxes/counts/row ids are currently flattened away in the returned support plan;
  - shared shell support tiles are not yet represented as CPBs or row-owned tiles, only as `current_box`, `inner_box`, and support rows.

6. Smallest validation smoke:
- No full H2 H1/H1-J route should be needed for support partition.
- Build only the independent H2 readiness/source-plan context:
  - use the existing independent input geometry;
  - construct system/recipe/parent/shells enough to obtain the route skeleton and target/source-plan payloads;
  - call the new private support-partition helper.
- Assert:
  - support order `(:atom_contact_core, :shared_shell_1, :shared_shell_2)`
  - support counts `(275, 578, 362)`
  - retained counts `(275, 98, 98)`
  - atom-contact tile counts/piece roles include two `:atom_local_core` and one `:midpoint_slab`
  - shared-shell tiles cover exactly the stored support rows and do not include inner-box rows
  - total duplicate/missing/outside counts are zero
  - fake/source-backed flags are false
  - no provider blocks or route-global matrices are materialized.

Recommended exact next implementation pass:
- Implement `_PQSIndependentH2PQSSupplementSupportPartitionPayload` and `_pqs_source_box_route_driver_independent_h2_pqs_supplement_support_partition_payload`.
- Extend or wrap the independent support-region plan just enough to retain atom-contact per-piece descriptors and shared-shell tile descriptors.
- Use existing `CartesianCPB`/`CartesianRouteCore` support concepts where they fit; if shared shells are not one-layer complete shells, add only a narrow rectangular-difference tiler and validate it against stored `support_indices`.
- Wire only a compact summary into assembly/reporting if needed for the next provider-block payload.
- Do not call CPB provider functions yet.

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
- quarantined: provider calls, provider blocks, mixed matrices, residual MWG, combined density readiness, supplemented values, CR2/export/HamV6/public API remain blocked.
- not deleted because: audit-only pass; no stale mirror deletion was needed.
- exact remaining caller/blocker: support rows are available in the independent support/source-plan payloads, but atom-contact per-piece descriptors and shared-shell support tiles are not yet carried as a compact route-owned partition; that partition is the next blocker before provider-block payload materialization.

-- repo-doer@macmini

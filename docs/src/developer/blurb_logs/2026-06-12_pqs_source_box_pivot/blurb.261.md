Pass 261 - independent H2 PQS supplement support-partition payload

Context:
- Current HEAD should include
  `2a4df13f Record independent H2 PQS support tiling audit`.
- Pass 260 found that provider blocks must wait for a route-owned
  support-partition payload. The key hazard is that shared-shell support is
  `outer_box \ inner_box`; the filled shared-shell `source_cpb` is not itself
  the support and must not be passed blindly to CPB provider helpers.
- The independent H2 PQS route is fake-free and currently has:
  - support counts `(275, 578, 362)`;
  - retained counts `(275, 98, 98)`;
  - final dimension `471`;
  - supplement preflight blocked by `:missing_provider_gto_supplement_blocks`.

Task:
Implement only the private route-owned support-partition payload needed before
independent H2 PQS supplement provider blocks.

Suggested object/helper names:
- `_PQSIndependentH2PQSSupplementSupportPartitionPayload`
- `_pqs_source_box_route_driver_independent_h2_pqs_supplement_support_partition_payload`

Implementation scope:
1. Extend or wrap the independent support-region/source-plan facts just enough
   to expose row-owned support partitions:
   - atom-contact core split into the two `:atom_local_core` pieces and one
     `:midpoint_slab` piece;
   - shared-shell pieces represented as outer-minus-inner support tiles or row
     maps, never as filled outer source CPBs.
2. The payload/summary should expose compact facts:
   - `status`, `blocker`;
   - `route_family`, `route_kind`;
   - `support_order`, `retained_order`;
   - `support_counts`, `retained_counts`, `retained_ranges`;
   - `unit_partitions`;
   - `total_tile_count`, `total_support_count`;
   - duplicate/missing/outside parent-row counts;
   - `coverage_complete`;
   - `source_backed_fixed_source_oracle_used = false`;
   - `fake_pqs_enabled = false`;
   - `retained_transform_authority = :pqs_source_box_construction`.
3. Per-unit partition fingerprints should include:
   - `unit_key`;
   - `support_count`;
   - `retained_range`;
   - `tile_count`;
   - `tile_support_counts`;
   - parent-row coverage/min/max and duplicate/missing/outside counts;
   - `coverage_complete`.
4. Per-tile fingerprints should include:
   - `tile_key`;
   - `unit_key`;
   - `source_region_role`;
   - `cpb_role` or `support_tile_role`;
   - `intervals` when rectangular;
   - `support_count`;
   - compact parent-row fingerprint, plus row ids if the route-owned payload
     needs them for the later provider-block pass;
   - `unit_local_row_range` or explicit local rows;
   - `retained_range`;
   - `provider_tile_ready`.
5. If existing `CartesianCPB.complete_shell_boundary_strata(outer, inner)` fits
   a shared shell, use it. If it does not fit the current shared-shell boxes,
   add only a narrow rectangular-difference tiler for these boxes and validate
   it against the stored support row ids.
6. Wire only a compact summary into `cartesian_assembly`/reporting if needed for
   the next provider-block payload. Do not create a broad report-key cloud.

Strict exclusions:
- Do not call CPB provider functions.
- Do not implement provider blocks.
- Do not build mixed gausslet/GTO matrices, GTO/GTO matrices, combined raw
  moment matrices, residual MWG representation, combined density readiness,
  supplemented H1/H1-J/RHF values, CR2/export, HamV6, or public API.
- Do not modify fake-PQS/WL source-backed routes except if a tiny shared helper
  must preserve existing behavior.
- Do not compare to fake-PQS/WL scalar results.
- Do not run broad stale low-order integration tests as a per-pass gate.

Validation:
- `git diff --check`.
- Run the smallest focused smoke that constructs the independent H2
  readiness/source-plan context and calls the new private support-partition
  helper.
- Assertions should cover:
  - support order `(:atom_contact_core, :shared_shell_1, :shared_shell_2)`;
  - support counts `(275, 578, 362)`;
  - retained counts `(275, 98, 98)`;
  - atom-contact tile roles include two `:atom_local_core` tiles and one
    `:midpoint_slab` tile;
  - shared-shell support tiles cover exactly stored support rows and exclude
    inner-box rows;
  - duplicate/missing/outside counts are zero;
  - fake/source-backed flags are false;
  - provider blocks and route-global matrices are not materialized.
- If the focused smoke requires a full driver route, explain why and report
  elapsed seconds. Do not run RHF or supplemented-value routes.

Line budget:
- Keep scoped `src + test + bin` net-negative if feasible by deleting a mature
  stale assertion/metadata candidate.
- If the narrow support-partition payload is slightly net-positive and no mature
  deletion is safe in the same pass, report the exact positive diff and why it
  should be treated as a temporary implementation exception. Do not broaden the
  pass just to force arithmetic.

Report:
- Files/functions changed.
- Whether a support-partition payload is available or precisely blocked.
- Tile counts and coverage counts by support unit.
- Validation command(s) and timing if a route smoke was run.
- Scoped line count for `src + test + bin`.
- Deletion/shrinkage result:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini

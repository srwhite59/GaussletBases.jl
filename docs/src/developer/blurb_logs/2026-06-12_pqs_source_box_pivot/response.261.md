Pass 261 response - independent H2 PQS supplement support-partition payload

Implemented the private support-partition boundary for the independent H2 PQS supplement preflight path.

Files/functions changed:
- `src/pqs_source_box_diatomic_complete_core_shell.jl`
  - added `_PQSIndependentH2PQSSupplementSupportPartitionPayload`;
  - added support-tile helpers for atom-contact rectangular pieces and shared-shell `outer_box \ inner_box` tiles;
  - extended the independent H2 support-region plan descriptors with route-owned support tiles;
  - added `_pqs_source_box_route_driver_independent_h2_pqs_supplement_support_partition_payload`.
- `src/pqs_source_box_route_driver_helpers.jl`
  - `cartesian_assembly` now builds and carries `independent_h2_pqs_supplement_support_partition_payload`.

Payload status:
- available: `:available_independent_h2_pqs_supplement_support_partition`
- support order: `(:atom_contact_core, :shared_shell_1, :shared_shell_2)`
- support counts: `(275, 578, 362)`
- retained counts: `(275, 98, 98)`
- total support count: `1215`
- total tile count: `55`
- duplicate/missing/outside parent-row counts: `0 / 0 / 0`
- coverage complete: `true`
- `source_backed_fixed_source_oracle_used = false`
- `fake_pqs_enabled = false`
- `retained_transform_authority = :pqs_source_box_construction`
- provider/global matrices remain unmaterialized.

Tile and coverage counts:
- `:atom_contact_core`: 3 tiles, support count 275, tile support counts `(125, 125, 25)`, zero duplicate/missing/outside rows.
- `:shared_shell_1`: 26 tiles, support count 578, tile support counts `(91, 91, 91, 91, 49, 49, 13, 13, 13, 13, 7, 7, 7, 7, 7, 7, 7, 7, 1, 1, 1, 1, 1, 1, 1, 1)`, zero duplicate/missing/outside rows.
- `:shared_shell_2`: 26 tiles, support count 362, tile support counts `(55, 55, 55, 55, 25, 25, 11, 11, 11, 11, 5, 5, 5, 5, 5, 5, 5, 5, 1, 1, 1, 1, 1, 1, 1, 1)`, zero duplicate/missing/outside rows.

Validation:
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed.
- `julia --project=. -e 't = @elapsed include("tmp/work/h2_independent_pqs_support_partition_smoke.jl"); println("elapsed_s=", t)'`
  - passed.
  - elapsed: `164.176072458` seconds.
  - This focused smoke used the full independent H2 driver context because the support partitions are assembled from route-owned target/source-plan facts. It did not run H1, H1-J, RHF, provider blocks, or supplemented-value routes.
- `git diff --check`
  - passed.

Scoped line count for `src + test + bin`:
- `src/pqs_source_box_diatomic_complete_core_shell.jl`: +429 / -0
- `src/pqs_source_box_route_driver_helpers.jl`: +6 / -0
- total: +435 / -0
- This is a temporary implementation exception: the pass introduced a new private authority payload and row/tile fingerprints before provider blocks. I did not find a mature stale deletion target inside this narrow boundary without broadening into unrelated cleanup.

Git status:
```text
## main...origin/main
 M src/pqs_source_box_diatomic_complete_core_shell.jl
 M src/pqs_source_box_route_driver_helpers.jl
```

Deletion/shrinkage result:
- deleted: none.
- simplified: support partition authority is now concentrated in one private payload instead of leaving the next provider-block pass to infer support rows from filled shared-shell source CPBs or adjacent route facts.
- quarantined: provider blocks, route-global matrices, mixed/GTO matrices, combined raw moments, residual MWG, supplemented values, RHF, CR2/export, HamV6, public API, and fake/WL comparison paths remain untouched.
- not deleted because: no mature stale assertion/helper/report field was safely deletable inside the pass-261 boundary; forcing deletion would have widened the task.
- exact remaining caller/blocker: the next provider-block payload should consume `independent_h2_pqs_supplement_support_partition_payload`; supplement preflight remains blocked on `:missing_provider_gto_supplement_blocks`.

-- repo-doer@macmini

Pass 194 review

Accepted.

This pass fixed the main concern from pass 193: the diatomic readiness report
surface is compact again. The on-disk H2 readiness artifact remains stable, but
the report no longer carries the full scalar alias cloud.

Accepted shrink:

- Removed duplicated scalar report aliases such as
  `diatomic_complete_core_shell_readiness_status`,
  `diatomic_complete_core_shell_source_plan_status`,
  `diatomic_complete_core_shell_h1_status`,
  `diatomic_complete_core_shell_rhf_materialized`, and related fields.
- Kept `diatomic_complete_core_shell_readiness_summary` as the compact report
  authority.
- Moved optional diagnostic facts needed by the artifact writer into that
  summary.
- Updated the H2 readiness artifact writer to derive route, basis, physics,
  density-interaction, and private-RHF artifact fields from the compact summary.
- Preserved the artifact keys checked by the H2 readiness test.

Accepted blocker diagnosis:

- The live H2 blocker is produced by
  `_pqs_source_box_route_driver_diatomic_complete_core_shell_ham_readiness_payload`.
- The immediate readiness blocker remains
  `:missing_diatomic_complete_core_shell_source_plan_producer`.
- This is not primarily a route-kind or fixed-q rejection. The H2 route is
  classified as bond-aligned diatomic, and q/n_s retained metadata is present.
- The next materialization seam depends on real diatomic parent/axis-bundle
  materialization and a narrow source-plan producer that consumes the existing
  raw-box/source-realization payloads.

Validation reviewed:

- Focused H2 readiness test passed, 29/29, in 43.3s on my rerun.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.
- Removed-alias search has no live `src/test/bin` hits.
- Scoped line budget is `28 added / 53 deleted`, net `-25`.

Next manager decision:

Pass 195 should be a narrow materialization step, not H1/H1-J/RHF. The most
natural next seam is the reviewed diatomic parent/axis-bundle materializer for
the H2 driver path, because the source-plan producer cannot honestly become
available until it has real parent axes/axis bundles to consume.

-- repo-manager@macmini

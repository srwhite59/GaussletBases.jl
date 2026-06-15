# PQS Thin Route Demolition History

This archive records the compressed history of the
`demolition/pqs-thin-route` branch. The current doctrine now lives in:

- [Cartesian route migration](../cartesian/route_migration.md)

## Compact History

- The demolition branch was created to remove Cartesian/PQS report, status,
  readiness, probe, print-line, and helper-schema scaffolding without preserving
  temporary green tests.
- The branch was allowed to break while old scaffolding was cut, but the driver
  entry point and input-driven workflow were later protected explicitly.
- The first cut deleted exact print-line tests, flat terminal route glue tests,
  staged low-order policy/status tests, and private field-cloud driver tests.
- Route-global and metadata demolition removed route-core sidecar adapters,
  legacy source-box fixture layers, and old White-Lindsey adapter-summary
  paths.
- Precompile workloads that exercised deleted WL/route-global adapter paths
  were removed instead of being repaired.
- A temporary overreach deleted too much driver/source-route code; that was
  rolled back, establishing the guardrail that `bin/cartesian_ham_builder.jl`
  and input-driven construction are protected.
- Old Cartesian nested helper tests, manifest mirror tests, endpoint key-list
  tests, and status/blocker/materialized field tests were deleted.
- The old low-order/WL route-configured materialization stack was collapsed
  before WL capability was restored compactly through driver-level smoke paths.
- Generic driver stages were then cleaned in execution order: system, recipe,
  parent, shells, units, transforms, pair terms, assembly, report,
  materialization, print, and save.
- At each stage, deleted upstream fields were not restored for downstream
  report code; downstream stages were simplified when execution reached them.
- The route stopped carrying broad fields such as `status`, `readiness`,
  `available`, `blocker`, `materialized`, `probe`, `pending`, `diagnostics`,
  and `object_kind` as driver-stage bookkeeping.
- The compact driver report now summarizes current route objects rather than
  old parent-axis/readiness/probe/raw-box scaffolding.
- Route-specific PQS payloads were normalized to keep real endpoint data for
  source plans, final basis, H1, H1-J, RHF, supplement preflight, and artifact
  paths.
- The full 2x2x2 driver ladder became the validation authority for the private
  route: WL/PQS by atomic/diatomic by pure/supplement policy.
- Temporary line ladders were added for `wl_atomic`, `wl_diatomic`,
  `pqs_atomic`, and `pqs_diatomic`.
- The old executable Cartesian `nested`, `ordinary`, and `diatomic` test groups
  were removed from the default test dispatcher.
- Dormant root-test Cartesian/QW/WL fixture factories were deleted after active
  radial/core fixture dependencies were preserved.
- After review, the thin route branch was fast-forwarded into `main` without a
  merge commit or squash, and the demolition branch name was deleted.
- The remaining old Cartesian code is now treated as feature-donor inventory,
  not dead code by default.
- Future migration should use compact driver ladders and feature-specific
  materialization checks, not restored helper/schema tests.

Pass 221 response - H2 no-supplement WL comparison artifact wiring

Implemented the bounded H2 physical gausslet-only comparison-artifact wiring.

Files changed:
- `bin/cartesian_ham_builder.jl`
- `src/pqs_source_box_route_driver_helpers.jl`
- `src/pqs_source_box_route_driver_reporting.jl`
- `test/driver_inputs/h2_pqs_q5_physical_gausslet_r4.jl`
- `test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl`
- `test/nested/cartesian_ham_builder_h2_pqs_q5_source_box_diagnostic_r4_runtests.jl`

Readiness before/after:
- Before: H2 physical input had `comparison_ready = false`, `comparison_blocker = :gausslet_only_reference_not_selected`, `physics_endpoint_ready = false`, `physics_endpoint_blocker = :missing_h2_gausslet_only_reference_comparison`; the artifact writer converted the reviewed WL candidate into `:missing_wl_h2_gausslet_only_reference_values`.
- After: H2 physical input has `comparison_ready = true`, `comparison_blocker = nothing`, `physics_endpoint_ready = true`, `physics_endpoint_blocker = nothing`, and `comparison_reference_label = "WL/QW H2 R=4 gausslet-only 463"`.

Comparison fields added to the physical H2 artifact:
- `comparison/wl_h1_lowest`
- `comparison/delta_h1`
- `comparison/wl_h1_self_coulomb`
- `comparison/delta_h1_j`
- `comparison/wl_rhf_electronic_energy`
- `comparison/delta_rhf_electronic_energy`
- `comparison/wl_rhf_nuclear_repulsion`
- `comparison/pqs_rhf_total_with_nuclear_repulsion`
- `comparison/wl_rhf_total_with_nuclear_repulsion`
- `comparison/delta_rhf_total_with_nuclear_repulsion`

The new scalar comparison fields are written only when the reviewed WL constants are present, so the H2 221 source-box diagnostic artifact does not acquire these physical endpoint comparison fields.

WL constants used:
- `wl_h1_lowest = -0.7946609179724673`
- `wl_h1_self_coulomb = 0.45696639804337047`
- `wl_rhf_one_electron_energy = -1.5611571934181985`
- `wl_rhf_electron_electron_energy = 0.40220533775308426`
- `wl_rhf_electronic_energy = -1.1589518556651142`
- `wl_rhf_nuclear_repulsion = 0.25`
- `wl_rhf_total_with_nuclear_repulsion = -0.9089518556651142`

PQS-vs-WL deltas:
- `delta_h1 = 2.55351295663786e-15`
- `delta_h1_j = 6.661338147750939e-16`
- `delta_rhf_electronic_energy = -3.2713831643604863e-12`
- `pqs_rhf_total_with_nuclear_repulsion = -0.9089518556683855`
- `delta_rhf_total_with_nuclear_repulsion = -3.2713831643604863e-12`

Supplemented WL/QW scalar references remain quarantined:
- `comparison/old_supplemented_wl_qw_scalar_references_blocked = true`
- no ambiguous molecular `comparison/wl_rhf_total` or `comparison/delta_rhf` was added for this H2 endpoint.

Scoped line-count result:
- Command: `git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`
- Added/deleted/net under the scoped paths: `102 added`, `105 deleted`, net `-3`.

Validation:
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed; final rerun printed `load ok` after precompile.
- `julia --project=. -e 'using Test; t = @elapsed include("test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl"); println("elapsed_s=", t)'`
  - passed: `91` assertions, `elapsed_s=89.513741708`.
- `julia --project=. test/nested/cartesian_ham_builder_h2_pqs_q5_source_box_diagnostic_r4_runtests.jl`
  - passed: `12` assertions, test summary time `1m39.8s`.
- `git diff --check`
  - passed.
- `git diff --cached --check`
  - passed.

Git status:
```text
## main...origin/main
 M bin/cartesian_ham_builder.jl
 M src/pqs_source_box_route_driver_helpers.jl
 M src/pqs_source_box_route_driver_reporting.jl
 M test/driver_inputs/h2_pqs_q5_physical_gausslet_r4.jl
 M test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl
 M test/nested/cartesian_ham_builder_h2_pqs_q5_source_box_diagnostic_r4_runtests.jl
```

Deletion/shrinkage report:
- deleted: no whole files.
- simplified: shrank `test/nested/cartesian_ham_builder_h2_pqs_q5_source_box_diagnostic_r4_runtests.jl` from development-era parent/materialization/source-plan assertions toward the live H2 221 diagnostic contract.
- quarantined: supplemented WL/QW scalar references remain blocked by `comparison/old_supplemented_wl_qw_scalar_references_blocked = true`; H2 221 diagnostic route remains source-box diagnostic only and has no new physical WL scalar comparison fields.
- not deleted because: the physical H2 endpoint test still owns the reviewed WL/PQS scalar gate and the H2 221 diagnostic test still protects the compact diagnostic contract after shrinkage.
- exact remaining caller/blocker: next cleanup can continue reducing stale diagnostic assertions once the H2 physical endpoint test is the accepted route-owned comparison gate; no source caller requires the removed diagnostic assertions.

-- repo-doer@macmini

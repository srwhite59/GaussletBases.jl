Pass 194 complete.

Files changed:

- `src/pqs_source_box_route_driver_helpers.jl`
- `src/pqs_source_box_route_driver_reporting.jl`

Report alias shrinkage:

- Removed the pass-193 scalar report aliases for the diatomic complete
  core/shell readiness surface:
  - `diatomic_complete_core_shell_readiness_status`
  - `diatomic_complete_core_shell_readiness_blocker`
  - `diatomic_complete_core_shell_source_plan_status`
  - `diatomic_complete_core_shell_final_basis_status`
  - `diatomic_complete_core_shell_h1_payload_status`
  - `diatomic_complete_core_shell_h1_status`
  - `diatomic_complete_core_shell_ham_input_status`
  - `diatomic_complete_core_shell_hamiltonian_handoff_status`
  - `diatomic_complete_core_shell_hamiltonian_consumer_contract_status`
  - `diatomic_complete_core_shell_final_dimension`
  - `diatomic_complete_core_shell_final_overlap_identity_error`
  - `diatomic_complete_core_shell_h1_lowest_energy`
  - `diatomic_complete_core_shell_density_gauge`
  - `diatomic_complete_core_shell_raw_pair_factor_convention`
  - `diatomic_complete_core_shell_nuclear_repulsion`
  - `diatomic_complete_core_shell_electron_count`
  - `diatomic_complete_core_shell_spin_sector`
  - `diatomic_complete_core_shell_final_basis_materialized`
  - `diatomic_complete_core_shell_h1_materialized`
  - `diatomic_complete_core_shell_h1_j_materialized`
  - `diatomic_complete_core_shell_ham_input_materialized`
  - `diatomic_complete_core_shell_hamiltonian_handoff_materialized`
  - `diatomic_complete_core_shell_hamiltonian_consumer_contract_materialized`
  - `diatomic_complete_core_shell_rhf_materialized`
  - `diatomic_complete_core_shell_public_api`
  - `diatomic_complete_core_shell_exports_materialized`
  - `diatomic_complete_core_shell_artifacts_materialized`
- Retained the compact report field
  `diatomic_complete_core_shell_readiness_summary`.
- Merged optional diagnostic facts needed by the artifact writer into that
  summary: final dimension, final-overlap identity error, H1 lowest energy,
  density gauge, raw pair-factor convention, nuclear repulsion, electron
  count, and spin sector.
- `rg` found no remaining live `src/test/bin` callers for the removed scalar
  aliases.

Artifact compatibility:

- Preserved the pass-193 on-disk JLD2 keys.
- The diatomic readiness writer now derives `route/*`, optional `basis/*`,
  optional `physics/*`, optional `density_interaction/*`, and
  `private_rhf/materialized` from
  `report.diatomic_complete_core_shell_readiness_summary`.
- The explicit H2 readiness artifact test still checks the same public artifact
  keys and passed unchanged.

Source-plan producer blocker diagnosis:

- The top-level blocker is produced by
  `_pqs_source_box_route_driver_diatomic_complete_core_shell_ham_readiness_payload`
  in `src/pqs_source_box_diatomic_complete_core_shell.jl`.
- In the current H2 driver run, that payload sees
  `source_plan_status != :available_pqs_diatomic_complete_core_shell_source_plan`,
  so its summary records:
  - `status = :blocked_diatomic_complete_core_shell_ham_readiness`
  - `blocker = :missing_diatomic_complete_core_shell_source_plan_producer`
- Existing private payloads before the blocker:
  - route skeleton and retained-unit/pair inventory metadata
  - diatomic support-window payload
  - diatomic raw-box route payload
  - diatomic source-realization payload
  - diatomic source-plan payload object
  - downstream final-basis, H1, ham-input, handoff, consumer-contract, and
    ham-readiness payload objects as blocked/not-materialized summaries
- Missing inputs for a real source-plan producer are still the parent/axis
  materialization path and source-plan materialization path:
  - current parent contract reports `parent_axis_bundle_object_available = false`
  - current parent materialization blocker is
    `:pending_reviewed_diatomic_parent_materializer`
  - current route-configured materializer readiness is
    `:blocked_missing_materializer_inputs`
  - the source-plan producer still needs real parent axis bundles/transforms,
    raw product/source realization wiring, and a producer that returns
    `:pqs_diatomic_complete_core_shell_source_plan`
- This is not primarily a route-kind rejection. The route is classified as
  `:bond_aligned_diatomic` under `:pqs_source_box`. It is also not a fixed-q
  policy rejection; q/n_s retained-unit metadata is already present. The
  immediate next seam is parent/axis-bundle materialization plus the narrow
  diatomic source-plan producer that consumes it.

Recommended pass-195 seam:

- Implement the smallest private route-owned diatomic source-plan producer
  path that turns the current blocked source-plan payload into an available
  `:pqs_diatomic_complete_core_shell_source_plan`.
- Keep it behind the existing private diatomic complete core/shell payloads.
- First consume or construct the real parent axis bundle needed by the current
  metadata path; then wire the raw-box/source-realization payload into the
  source-plan producer.
- Do not add H2 final-basis, H1, H1-J, RHF, supplemented comparison, exports,
  or public API in that pass.

Scoped line-budget arithmetic:

```text
13	46	src/pqs_source_box_route_driver_helpers.jl
15	7	src/pqs_source_box_route_driver_reporting.jl
```

Total for the requested scoped diff: 28 added, 53 deleted, net `-25`.

Validation results:

- `julia --project=. test/nested/cartesian_ham_builder_h2_pqs_q5_gausslet_only_r4_readiness_runtests.jl`
  - passed: 29/29
  - elapsed: 1m42.9s
  - over-60s runtime was driven by package precompilation:
    `GaussletBases` precompiled in about 57s
  - artifact path printed by the test:
    `/Users/srw/dmrgtmp/jl_8yMGod/h2_pqs_q5_gausslet_only_r4.jld2`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed: `load ok`
- `git diff --check`
  - passed
- `git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`
  - passed and shown above
- `git status --short --branch`
  - `## main...origin/main`
  - ` M src/pqs_source_box_route_driver_helpers.jl`
  - ` M src/pqs_source_box_route_driver_reporting.jl`

Deletion/shrinkage report:

- deleted:
  - no files
- simplified:
  - removed the duplicated diatomic readiness scalar report alias cloud
  - kept the artifact-facing scalar keys stable while moving report authority
    back to one compact summary
- quarantined:
  - nothing
- not deleted because:
  - the H2 driver input and explicit readiness artifact test are current
    endpoint-readiness coverage
  - the He driver endpoints and reference tests are live
- exact remaining caller/blocker:
  - H2 PQS materialization remains blocked at
    `:missing_diatomic_complete_core_shell_source_plan_producer`
  - lower current blockers are the reviewed diatomic parent/axis materializer
    and the source-plan producer that consumes real parent axis/source
    realization payloads
  - supplemented WL/QW H2 comparison remains blocked until supplement policy is
    represented explicitly in the driver/PQS path

-- repo-doer@macmini

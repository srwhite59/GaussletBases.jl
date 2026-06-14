Pass 207 complete.

Files renamed:
- `test/driver_inputs/h2_pqs_q5_gausslet_only_r4.jl`
  -> `test/driver_inputs/h2_pqs_q5_source_box_diagnostic_r4.jl`
- `test/nested/cartesian_ham_builder_h2_pqs_q5_gausslet_only_r4_readiness_runtests.jl`
  -> `test/nested/cartesian_ham_builder_h2_pqs_q5_source_box_diagnostic_r4_runtests.jl`

Files created:
- `docs/src/developer/cartesian_driver_endpoint_manifest.md`

Files edited:
- `test/nested/cartesian_ham_builder_h2_pqs_q5_source_box_diagnostic_r4_runtests.jl`

What changed:
- Renamed the H2 221 diagnostic input/test from ambiguous `gausslet_only_r4` naming to `source_box_diagnostic_r4`.
- Updated the renamed test's input path, testset name, temporary output names, and diagnostic artifact print label.
- Moved the H2 221 role/blocker assertions to the start of the artifact checks:
  - `route/artifact_role = :source_box_diagnostic`
  - `physics/endpoint_ready = false`
  - `physics/endpoint_blocker = :retained_atom_core_interiors_missing`
  - `basis/retained_atom_core_interiors = false`
  - `basis/source_plan_role = :boundary_source_box_diagnostic`
  - `comparison/ready = false`
  - `comparison/blocker = :supplemented_reference_not_comparable_to_gausslet_only`
- Removed one non-contract `println` from the renamed diagnostic test to keep the src/test/bin scoped line budget net-negative.
- Added a compact endpoint manifest table.

Endpoint manifest summary:
- He PQS q5 WL-map physical atom endpoint: `he_pqs_q5_wlmap.jl`, dimension 419, endpoint-ready, no blocker, not default-runner.
- He PQS q5 WL-map RHF diagnostic/endpoint check: `he_pqs_q5_wlmap.jl`, dimension 419, endpoint-ready, no blocker, not default-runner.
- H2 PQS q5 source-box diagnostic: `h2_pqs_q5_source_box_diagnostic_r4.jl`, dimension 221, not endpoint-ready, blocker `:retained_atom_core_interiors_missing`, not default-runner.
- H2 PQS q5 physical gausslet target inventory: `h2_pqs_q5_physical_gausslet_r4.jl`, expected dimension 463, not endpoint-ready, blocker `:missing_physical_gausslet_source_plan`, not default-runner.

Source/test/bin scoped line budget:
- Actual rename-aware accounting:
  - deleted: 140 lines
  - added: 139 lines
  - net: -1 line
- Old files:
  - `test/driver_inputs/h2_pqs_q5_gausslet_only_r4.jl`: 38 lines
  - `test/nested/cartesian_ham_builder_h2_pqs_q5_gausslet_only_r4_readiness_runtests.jl`: 102 lines
- New files:
  - `test/driver_inputs/h2_pqs_q5_source_box_diagnostic_r4.jl`: 38 lines
  - `test/nested/cartesian_ham_builder_h2_pqs_q5_source_box_diagnostic_r4_runtests.jl`: 101 lines
- Raw command output before staging:
  - `git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`
  - `0	38	test/driver_inputs/h2_pqs_q5_gausslet_only_r4.jl`
  - `0	102	test/nested/cartesian_ham_builder_h2_pqs_q5_gausslet_only_r4_readiness_runtests.jl`
  - The new renamed files are currently untracked, so this raw unstaged command lists the deleted old paths only.

Validation:
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - Passed; printed `load ok`.
- `julia --project=. -e 'using Test; t = @elapsed include("test/nested/cartesian_ham_builder_h2_pqs_q5_source_box_diagnostic_r4_runtests.jl"); println("elapsed_s=", t)'`
  - Passed.
  - Test summary: 52 passed / 52 total.
  - Test-reported time: 1m26.1s.
  - Julia `@elapsed`: 86.515485708 seconds.
- `julia --project=. -e 'using Test; t = @elapsed include("test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl"); println("elapsed_s=", t)'`
  - Passed.
  - Test summary: 37 passed / 37 total.
  - Test-reported time: 1m01.2s.
  - Julia `@elapsed`: 61.555269917 seconds.
  - I expected this to stay under 60 seconds because it skips final basis/H1, but it still ran just over 60 seconds due remaining driver compilation.
- `rg -n "h2_pqs_q5_gausslet_only_r4" test bin docs/src/developer --glob '!docs/src/developer/blurb_logs/**'`
  - No non-log hits.
- `git diff --check`
  - Passed.
- `git diff --cached --check`
  - Passed; no staged changes.

Git status:
```text
## main...origin/main
 D test/driver_inputs/h2_pqs_q5_gausslet_only_r4.jl
 D test/nested/cartesian_ham_builder_h2_pqs_q5_gausslet_only_r4_readiness_runtests.jl
?? docs/src/developer/cartesian_driver_endpoint_manifest.md
?? test/driver_inputs/h2_pqs_q5_source_box_diagnostic_r4.jl
?? test/nested/cartesian_ham_builder_h2_pqs_q5_source_box_diagnostic_r4_runtests.jl
```

Deletion/shrinkage report:
- deleted: old ambiguous H2 221 diagnostic filenames; one non-contract diagnostic artifact `println`.
- simplified: H2 221 diagnostic test naming and first assertions now match its source-box diagnostic role before scalar/H1 details.
- quarantined: H2 221 remains explicitly labeled `:source_box_diagnostic`, not a physics endpoint.
- not deleted because: H2 221 scalar/H1 diagnostic checks still validate the current driver artifact; H2 463 remains target inventory only until source-plan work is assigned.
- exact remaining caller/blocker: H2 physical gausslet endpoint remains blocked by `:missing_physical_gausslet_source_plan`; H2 221 diagnostic remains blocked as a physics endpoint by `:retained_atom_core_interiors_missing`.

-- repo-doer@macmini

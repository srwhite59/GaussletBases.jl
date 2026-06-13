Pass 190 response

Implemented the optional private RHF route-smoke diagnostic for the He q=5/n_s=5 PQS driver artifact.

Deletion verification:

- `test/nested/pqs_route_driver_one_center_materializer_probe_runtests.jl` had no live runner include or source caller in the pass search; deleted.
- The first deletion was not enough for the true line budget after adding the new endpoint, so I also checked `test/nested/cartesian_white_lindsey_adapter_timing_probe.jl`. It had no source/runner/doc caller and was an explicitly manual development timing probe; deleted.
- Other obvious materializer probe tests were not deleted because `test/nested/integration_runtests.jl` still includes them.

Driver inputs added:

- `run_private_rhf = false`
- `private_rhf_electron_count = nothing`
- `private_rhf_fixture_role = :route_smoke`
- `private_rhf_mixing_kind = :fock_diis`
- `private_rhf_max_iterations = 25`
- `private_rhf_density_atol = 1.0e-8`
- `private_rhf_energy_atol = 1.0e-10`
- `private_rhf_residual_atol = 1.0e-8`
- `private_rhf_trace_atol = private_rhf_density_atol`
- `private_rhf_idempotency_atol = private_rhf_density_atol`
- `private_rhf_max_history = nothing`
- `private_rhf_diis_start_iteration = 2`
- `private_rhf_diis_regularization = 1.0e-12`
- `private_rhf_diis_coefficient_max_abs = 25.0`
- `wl_rhf_total = nothing`

Where RHF is called:

- `bin/cartesian_ham_builder.jl` still exposes the staged driver spine.
- `cartesian_assembly` now calls `_pqs_source_box_route_driver_complete_core_shell_private_rhf_payload(...)` only when `recipe.private_rhf_inputs.run_private_rhf == true`.
- That helper consumes `assembly.complete_core_shell_diagnostic_route_payload` data: source plan, final basis, H1 payload, density inputs, Coulomb expansion, H1/J payload, and pre-final density interaction.
- Electron count uses explicit `private_rhf_electron_count` when supplied; otherwise it conservatively derives the one-center closed-shell He value from the one-center nuclear charge.
- The helper calls `_pqs_multilayer_complete_core_shell_rhf_input_contract(...)` and `_pqs_multilayer_complete_core_shell_rhf_scf_payload(...)`. No public API/export/product RHF path was added.

Artifact keys written when `run_private_rhf=true`:

- `private_rhf/status`
- `private_rhf/blocker`
- `private_rhf/total_energy`
- `private_rhf/iteration_count`
- `private_rhf/converged`
- `private_rhf/residual`
- `private_rhf/mixing_kind`
- `comparison/wl_rhf_total`
- `comparison/delta_rhf`

The existing default H1/H1-J artifact keys are unchanged when `run_private_rhf=false`; private RHF keys are written only for requested RHF runs.

Explicit RHF endpoint:

- Added `test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_rhf_runtests.jl`.
- It runs `bin/cartesian_ham_builder.jl` with `test/driver_inputs/he_pqs_q5_wlmap.jl` and ARGS overrides:
  - `run_private_rhf=true`
  - `wl_rhf_total=-2.85080350301779`
  - temp `outfile`
  - temp `tsvfile`
- It reads only the JLD2 artifact and checks:
  - `private_rhf/converged == true`
  - `private_rhf/total_energy ≈ -2.850817886618113`
  - `comparison/wl_rhf_total ≈ -2.85080350301779`
  - `comparison/delta_rhf ≈ -1.4383600322798173e-5`

Endpoint runtime and printed phases:

- `julia --project=. test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_rhf_runtests.jl`
- Passed: 5/5, `3m16.2s`
- Printed phases:
  - `Transforming: 12.340135 seconds`, 99.99% compilation
  - `Pair terms: 25.109986 seconds`, 100.00% compilation
  - `Assembly: 32.339243 seconds`, 34.87M allocations, 5.222 GiB, 92.49% compilation

Default-runner confirmation:

- `rg -n "he_pqs_q5_wlmap_rhf|pqs_route_driver_one_center_materializer_probe" test/nested/runtests.jl test/nested/integration_runtests.jl test/nested`
- Only hit was the new explicit endpoint file; no default nested or integration runner include.

Line budget:

- `git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`

```text
23	1	bin/cartesian_ham_builder.jl
134	1	src/pqs_source_box_route_driver_helpers.jl
30	0	src/pqs_source_box_route_driver_reporting.jl
0	166	test/nested/cartesian_white_lindsey_adapter_timing_probe.jl
0	167	test/nested/pqs_route_driver_one_center_materializer_probe_runtests.jl
```

`git diff --numstat` does not include the untracked new endpoint test. Counting `wc -l test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_rhf_runtests.jl` = 40 lines, true source/test/bin arithmetic is:

- added: 23 + 134 + 30 + 40 = 227
- deleted: 1 + 1 + 166 + 167 = 335
- net: -108 lines

Validation:

- `julia --project=. test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_rhf_runtests.jl` passed 5/5 in 3m16.2s.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.
- `rg -n "[ \t]+$" test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_rhf_runtests.jl` had no trailing-whitespace hits.
- `rg -n "he_pqs_q5_wlmap_rhf|pqs_route_driver_one_center_materializer_probe" test/nested/runtests.jl test/nested/integration_runtests.jl test/nested` confirmed no default-runner include.
- `git status --short --branch`:

```text
## main...origin/main
 M bin/cartesian_ham_builder.jl
 M src/pqs_source_box_route_driver_helpers.jl
 M src/pqs_source_box_route_driver_reporting.jl
 D test/nested/cartesian_white_lindsey_adapter_timing_probe.jl
 D test/nested/pqs_route_driver_one_center_materializer_probe_runtests.jl
?? test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_rhf_runtests.jl
```

Deletion/shrinkage report:

- deleted:
  - `test/nested/pqs_route_driver_one_center_materializer_probe_runtests.jl`
  - `test/nested/cartesian_white_lindsey_adapter_timing_probe.jl`
- simplified:
  - private RHF is route-owned through the existing complete core/shell diagnostic route payload instead of a separate hand-built test fixture.
  - default artifact behavior stays non-RHF; private RHF scalar artifact keys are opt-in.
- quarantined:
  - private RHF remains diagnostic/private and default-off.
- not deleted because:
  - included materializer/report integration tests remain live in `test/nested/integration_runtests.jl`.
  - no RHF matrices, densities, orbitals, Fock matrices, or histories were written to artifacts.
- exact remaining caller/blocker:
  - remaining RHF endpoint is explicit-only in `test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_rhf_runtests.jl`; it is not included by default runners.

-- repo-doer@macmini

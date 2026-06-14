Pass 217 response - physical H2 private RHF execution adapter

Implemented the physical H2 private RHF execution adapter.

Files changed:
- `src/pqs_multilayer_complete_core_shell_rhf.jl`
- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `test/driver_inputs/h2_pqs_q5_physical_gausslet_r4.jl`
- `test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl`
- deleted `test/nested/cartesian_pair_block_route_global_matrix_set_smoke_runtests.jl`
- deleted `test/nested/cartesian_pair_block_route_global_one_body_adapter_runtests.jl`

What changed:
- Generalized the existing private RHF helper predicates/accessors to accept the reviewed physical H2 object family explicitly:
  - physical RHF input contract status `:available_pqs_physical_gausslet_rhf_input_contract`;
  - physical H1 payload status `:available_pqs_physical_gausslet_h1_payload`;
  - physical density interaction kind/status `:pqs_physical_gausslet_pre_final_density_interaction` / `:materialized_pqs_physical_gausslet_pre_final_density_interaction`.
- Reused the existing private RHF lower math for:
  - H1-Aufbau initial density;
  - one-step Fock/energy;
  - Fock-DIIS SCF iteration;
  - final recomputed one-step/residual diagnostics.
- Wired `_pqs_source_box_route_driver_diatomic_physical_gausslet_rhf_execution_payload` to execute the reused SCF helper for the physical H2 input contract.
- Kept artifact output compact: status/blocker, executed/materialized/converged, electron count, occupation, energies, iteration count, density trace, residuals, energy delta, and final-density/one-step consistency status. No matrices were added.
- Updated the H2 endpoint test to assert execution and to accept the required materialized or nonconverged execution state.

Private RHF execution result:
- attempted: yes
- converged: yes
- status: `:materialized_pqs_physical_gausslet_private_rhf_execution`
- blocker: `nothing`
- endpoint blocker: `:missing_h2_gausslet_only_reference_comparison`
- electron count source: explicit `private_rhf_electron_count`
- electron count: 2
- occupation policy: closed-shell RHF, `nocc = 1`
- total energy: `-1.1589518556683855`
- one-body energy: `-1.5611567624163218`
- two-body energy: `0.40220490674793624`
- iteration count: 8
- density trace: `1.9999999999999987`
- idempotency residual: `4.163336342344337e-17`
- commutator residual: `3.2406053396272227e-9`
- energy delta: `2.2115642650533118e-13`
- final-density/one-step consistency status: `:reviewed_recomputed`

Validation:
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed, `load ok`
- `julia --project=. -e 'using Test; t = @elapsed include("test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl"); println("elapsed_s=", t)'`
  - passed, 79/79 tests
  - elapsed: `142.753179292` seconds
- Capture/readback command for compact RHF fields:
  - `julia --project=. bin/cartesian_ham_builder.jl test/driver_inputs/h2_pqs_q5_physical_gausslet_r4.jl 'outfile="/private/tmp/h2_pqs_q5_physical_gausslet_r4_pass217.jld2"' 'tsvfile="/private/tmp/h2_pqs_q5_physical_gausslet_r4_pass217.tsv"' save_artifact=true save_tsv=true`
  - `julia --project=. -e 'using JLD2; path="/private/tmp/h2_pqs_q5_physical_gausslet_r4_pass217.jld2"; ...'`
- `git diff --check`
  - passed
- `git diff --cached --check`
  - passed

Line budget:
- scoped command:
  - `git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`
- added: 236
- deleted: 1020
- net: -784

Deletion/shrinkage report:
- deleted:
  - `test/nested/cartesian_pair_block_route_global_matrix_set_smoke_runtests.jl`
  - `test/nested/cartesian_pair_block_route_global_one_body_adapter_runtests.jl`
- simplified:
  - physical H2 private RHF now reuses the existing private RHF SCF helper instead of remaining as a blocked adapter placeholder.
- quarantined:
  - none
- not deleted because:
  - `test/nested/cartesian_ham_builder_one_center_config_smoke_runtests.jl` remains active in `test/nested/integration_runtests.jl` as the one-center driver smoke.
- exact remaining caller/blocker:
  - physics endpoint remains blocked only by missing gausslet-only H2 reference comparison: `:missing_h2_gausslet_only_reference_comparison`.
  - Manager acceptance correction: `run_private_rhf` was changed to `true`
    before acceptance so the visible driver input explicitly requests the
    private diagnostic RHF execution. The path remains private/report-only and
    still does not mark solver/export or public HamV6 behavior ready.

Git status:
```text
## main...origin/main
 M src/pqs_multilayer_complete_core_shell_rhf.jl
 M src/pqs_source_box_diatomic_complete_core_shell.jl
 M test/driver_inputs/h2_pqs_q5_physical_gausslet_r4.jl
 M test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl
 D test/nested/cartesian_pair_block_route_global_matrix_set_smoke_runtests.jl
 D test/nested/cartesian_pair_block_route_global_one_body_adapter_runtests.jl
```

-- repo-doer@macmini

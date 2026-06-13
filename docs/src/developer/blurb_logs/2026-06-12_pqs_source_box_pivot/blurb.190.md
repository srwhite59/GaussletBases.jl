Pass 190 - add optional private RHF to the He driver artifact, paid by deleting stale probe pressure.

Purpose:

The He q=5/n_s=5 PQS gate is now driver-owned and explicit:

```text
test/driver_inputs/he_pqs_q5_wlmap.jl
-> bin/cartesian_ham_builder.jl
-> saved artifact
-> test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_runtests.jl
```

Now add the optional private RHF route-smoke diagnostic behind a default-off
driver flag. This is a Hamiltonian validator, not the product. HFDMRG remains
the serious HF/DMRG path outside this repo.

Physics target:

Same He 419 WL-mapped PQS endpoint, with optional private RHF:

```text
PQS RHF total ≈ -2.850817886618113
WL  RHF total ≈ -2.85080350301779
delta ≈ -1.4383600322798173e-5
```

Task type:

Implementation plus deletion. This pass must be source/test/bin/generator
net-negative.

Deletion target:

Before adding code, verify whether this file has any live runner include or
source caller:

```text
test/nested/pqs_route_driver_one_center_materializer_probe_runtests.jl
```

Current manager inspection found no include in `test/nested/runtests.jl` or
`test/nested/integration_runtests.jl`, and no source caller. It appears to be
stale one-center materializer probe pressure from the older route-driver
transition.

If that remains true, delete the file in this pass. That should pay most or all
of the line budget for the private RHF driver seam.

If it has a live caller that manager missed, do not delete it silently. Report
the exact caller. If you cannot keep the pass source/test/bin/generator
net-negative after preserving live tests, write `.agent_handoffs/ATTENTION.md`
and stop.

Implementation task:

1. Add default-off driver inputs in `bin/cartesian_ham_builder.jl`.

   Suggested names:

   ```julia
   run_private_rhf = false
   private_rhf_electron_count = nothing
   private_rhf_fixture_role = :route_smoke
   private_rhf_mixing_kind = :fock_diis
   private_rhf_max_iterations = 25
   private_rhf_density_atol = 1.0e-8
   private_rhf_energy_atol = 1.0e-10
   private_rhf_residual_atol = 1.0e-8
   private_rhf_trace_atol = private_rhf_density_atol
   private_rhf_idempotency_atol = private_rhf_density_atol
   private_rhf_max_history = nothing
   private_rhf_diis_start_iteration = 2
   private_rhf_diis_regularization = 1.0e-12
   private_rhf_diis_coefficient_max_abs = 25.0
   wl_rhf_total = nothing
   ```

   Keep defaults false/nonintrusive. Existing driver behavior must not change
   unless `run_private_rhf=true`.

2. Carry those controls to the route-owned assembly/report path.

   Use the existing private RHF helpers in:

   ```text
   src/pqs_multilayer_complete_core_shell_rhf.jl
   ```

   The driver path should build:

   ```text
   _pqs_multilayer_complete_core_shell_rhf_input_contract(...)
   _pqs_multilayer_complete_core_shell_rhf_scf_payload(...)
   ```

   using data already present in:

   ```text
   assembly.complete_core_shell_diagnostic_route_payload
   ```

   Required route-owned inputs are already there:

   ```text
   source_payload.source_plan
   source_payload.coulomb_expansion
   final_basis
   h1_payload
   density_inputs
   h1_j_payload
   complete_core_shell_ham_payload.density_interaction
   ```

   Electron count:

   - Prefer explicit `private_rhf_electron_count` when supplied.
   - For this He route, it is acceptable to derive `2` from the one-center
     nuclear charge if the input is `nothing`, but keep that derivation
     conservative and clearly metadata-labeled.
   - Do not infer general molecular occupation policy beyond closed-shell
     private diagnostics in this pass.

3. Add compact report/artifact fields only.

   Suggested artifact keys:

   ```text
   private_rhf/status
   private_rhf/blocker
   private_rhf/total_energy
   private_rhf/iteration_count
   private_rhf/converged
   private_rhf/residual
   private_rhf/mixing_kind
   comparison/wl_rhf_total
   comparison/delta_rhf
   ```

   Residual should use `summary.residual_diagnostics.commutator_residual` when
   available. If the RHF payload is blocked/unavailable, write compact status
   and blocker without fake energy values.

   Do not write large RHF matrices, densities, orbitals, Fock matrices, or
   iteration histories into the driver artifact in this pass.

4. Add an explicit endpoint test, not a default nested include.

   Add a small explicit test, likely:

   ```text
   test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_rhf_runtests.jl
   ```

   It should:

   - run the same input file with ARGS overrides:

     ```text
     run_private_rhf=true
     wl_rhf_total=-2.85080350301779
     outfile=<temp>
     tsvfile=<temp>
     ```

   - read the JLD2 artifact only;
   - assert compact RHF facts:

     ```text
     private_rhf/converged == true
     private_rhf/total_energy ≈ -2.850817886618113
     comparison/wl_rhf_total ≈ -2.85080350301779
     comparison/delta_rhf ≈ -1.4383600322798173e-5
     ```

   Do not include this test in `test/nested/runtests.jl`. It is an explicit
   physics endpoint, like the H1/H1-J driver test.

Trust boundary:

- Private RHF only. Do not add public solver API, exports, artifacts intended
  for downstream HF, or any HFDMRG integration.
- Do not make RHF part of default driver behavior.
- Do not add RHF to the default nested runner.
- Do not add H2, Be2, Cr2, DMRG, HFDMRG, GTO, ECP, Qiu-White correction, or
  production solver behavior.
- Do not resurrect the deleted hand-built He H1/H1-J test.
- Do not request interactive approval or sandbox escalation. If approval is
  genuinely required, write `.agent_handoffs/ATTENTION.md` and stop.

Line-budget rule:

This pass must be net-negative for:

```text
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Acceptance condition:

```text
sum(deleted) > sum(added)
```

Expected budget source:

```text
delete test/nested/pqs_route_driver_one_center_materializer_probe_runtests.jl
```

If that deletion is not safe or not enough, do not pad with arbitrary
formatting deletion. Either find a genuinely stale source/test surface and
report it, or write `ATTENTION.md`.

Validation:

Run focused validation:

```text
julia --project=. test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_rhf_runtests.jl
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
git status --short --branch
```

Also confirm the RHF endpoint is not in the default nested runner:

```text
rg -n "he_pqs_q5_wlmap_rhf|pqs_route_driver_one_center_materializer_probe" test/nested/runtests.jl test/nested/integration_runtests.jl test/nested
```

The explicit RHF endpoint may take several minutes. Report its runtime and
dominant printed phases if available. Do not run broad nested or integration
suites in this pass.

Report:

- deletion verification for the stale materializer probe test;
- exact driver input names added;
- where RHF is called from the route-owned payload path;
- artifact keys written;
- explicit RHF endpoint runtime and result values;
- confirmation the RHF endpoint is not default-runner included;
- line-budget arithmetic;
- validation results;
- deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

Write the result to `.agent_handoffs/response.190.md` and copy it to:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.190.md
```

-- repo-manager@macmini

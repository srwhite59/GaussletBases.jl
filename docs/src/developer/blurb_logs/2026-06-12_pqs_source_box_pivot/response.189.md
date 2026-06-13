Pass 189 complete.

Implementation change:

- Removed `include("cartesian_ham_builder_he_pqs_q5_wlmap_runtests.jl")` from `test/nested/runtests.jl`.
- Kept `test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_runtests.jl` as an explicit endpoint test.
- Changed that endpoint testset label to include `cartesian_ham_builder_he_pqs_q5_wlmap`, so the requested `rg` validation can confirm the endpoint remains without matching the default runner.

Default runner status:

- `test/nested/runtests.jl` no longer includes the multi-minute He driver endpoint.
- The explicit endpoint remains available directly:

```text
julia --project=. test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_runtests.jl
```

Private RHF seam audit:

1. Current private RHF helper to call:

- The future driver flag should build an input contract with `_pqs_multilayer_complete_core_shell_rhf_input_contract`.
- For an actual route-smoke RHF total, it should then call `_pqs_multilayer_complete_core_shell_rhf_scf_payload`.
- `_pqs_multilayer_complete_core_shell_rhf_initial_density_payload` and `_pqs_multilayer_complete_core_shell_rhf_one_step_payload` are supporting pieces; the driver should not expose those as separate user-facing route stages for the first flag.

2. Existing data already available:

- `assembly.complete_core_shell_diagnostic_route_payload` already carries the needed private route inputs:
  - `source_payload.source_plan`
  - `source_payload.coulomb_expansion`
  - `final_basis`
  - `h1_payload`
  - `density_inputs`
  - `h1_j_payload`
  - `complete_core_shell_ham_payload.density_interaction`
- Electron count can be derived from the one-center system/nuclear charges for He as `2`.
- The existing artifact already writes H1, H1-J, density gauge, final dimension, retained shell counts, and WL H1/J deltas.

3. Suggested input flags:

```julia
run_private_rhf = false
private_rhf_max_iterations = 25
private_rhf_density_atol = 1.0e-8
private_rhf_energy_atol = 1.0e-10
private_rhf_residual_atol = 1.0e-8
private_rhf_mixing_kind = :fock_diis
private_rhf_fixture_role = :route_smoke
wl_rhf_total = -2.85080350301779
```

Use `private_rhf_max_iterations` to match the helper's current keyword. A short alias `private_rhf_max_iter` could be accepted in the driver input, but should resolve to the helper spelling before calling the payload.

4. Future RHF artifact keys:

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

The residual should come from `summary.residual_diagnostics.commutator_residual` when available.

5. Later RHF test shape:

- Reuse the same driver input and override `run_private_rhf=true`.
- Read the artifact only.
- Check:
  - `private_rhf/status`
  - `private_rhf/converged`
  - `private_rhf/total_energy ≈ -2.850817886618113`
  - `comparison/wl_rhf_total ≈ -2.85080350301779`
  - `comparison/delta_rhf ≈ -1.4383600322798173e-5`
- Do not reconstruct source plans, final bases, density interactions, Fock matrices, or SCF internals in the test.

Recommended next implementation pass:

- Add the `run_private_rhf=false` defaults and RHF control inputs to `bin/cartesian_ham_builder.jl`.
- In `cartesian_assembly`, or a compact route-owned payload derived immediately from `complete_core_shell_diagnostic_route_payload`, conditionally call the RHF input contract and SCF payload only when `run_private_rhf=true`.
- Extend `cartesian_save` with the compact `private_rhf/*` and `comparison/delta_rhf` keys.
- Add a focused explicit endpoint test that runs the existing He input with `run_private_rhf=true`.
- Keep it out of `test/nested/runtests.jl`.

Validation:

```text
rg -n "cartesian_ham_builder_he_pqs_q5_wlmap" test/nested/runtests.jl test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_runtests.jl
```

returned only the explicit endpoint test:

```text
test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_runtests.jl:7:@testset "cartesian_ham_builder_he_pqs_q5_wlmap artifact" begin
```

```text
julia --project=. -e 'using GaussletBases; println("load ok")'
```

passed with `load ok`.

```text
git diff --check
```

passed with no whitespace errors.

```text
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

reported:

```text
1	1	test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_runtests.jl
0	1	test/nested/runtests.jl
```

Line budget:

- added: `1`
- deleted: `2`
- net: `-1`

```text
git status --short --branch
```

returned:

```text
## main...origin/main
 M test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_runtests.jl
 M test/nested/runtests.jl
```

I did not rerun the multi-minute driver endpoint because this pass did not change the endpoint's behavior or artifact assertions.

Deletion/shrinkage report:

- deleted: one default-runner include for the multi-minute He driver endpoint.
- simplified: default nested runner no longer pays the explicit physics endpoint cost.
- quarantined: the He driver endpoint is now explicit/manual, not default.
- not deleted because: the endpoint itself is the driver-owned He PQS physics acceptance gate and should remain available.
- exact remaining caller/blocker: no default nested caller remains; future RHF work is blocked only on implementing the optional `run_private_rhf` driver/materialization path.

-- repo-doer@macmini

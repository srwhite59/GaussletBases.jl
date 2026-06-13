Pass 196 - materialize the H2 final basis only.

Purpose:

The H2 gausslet-only driver path now has real diatomic parent axes/axis bundles
and an available source plan. The current readiness blocker is:

```text
:missing_diatomic_complete_core_shell_final_basis_consumer
```

This pass should materialize the H2 final basis through the driver artifact
path, but still not build H1, H1-J, density interaction, or private RHF.

Task type:

Narrow implementation plus focused test update.

Current state:

- H2 input:
  `test/driver_inputs/h2_pqs_q5_gausslet_only_r4.jl`
- H2 readiness test:
  `test/nested/cartesian_ham_builder_h2_pqs_q5_gausslet_only_r4_readiness_runtests.jl`
- Parent axes are now constructed, with counts `(x = 9, y = 9, z = 15)`.
- Source plan is available:
  `:available_pqs_diatomic_complete_core_shell_source_plan`.
- Final basis is not materialized yet because the H2 input requests no H1/H1-J/RHF.

Exact task:

1. Add the smallest explicit final-basis request surface needed by the visible
   driver, if one does not already exist.

   Preferred shape:

   ```julia
   run_final_basis = true
   ```

   or an existing equivalent. Do not use `run_h1 = true` as a backdoor just to
   force final-basis construction.

2. Materialize the H2 final basis when the source plan is available and the
   final-basis request is true.

   Reuse the existing diatomic complete-core/shell final-basis payload path.
   Do not add a second final-basis implementation path.

3. Extend the H2 artifact/test to assert the final-basis state.

   Expected assertions should include:

   ```text
   route/source_plan_status == :available_pqs_diatomic_complete_core_shell_source_plan
   route/final_basis_status indicates materialized/available final basis
   basis/final_dimension is present and positive
   basis/final_overlap_identity_error is present if available as a scalar diagnostic
   route/h1_status remains not materialized
   route/h1_j_materialized == false
   private_rhf/materialized == false
   comparison/ready == false
   ```

   Do not hard-code a final dimension unless the implementation has an
   independent reason to expect it. If the dimension is surprising, report it
   rather than forcing the test to match an old assumption.

4. Report the new blocker after final basis materializes.

   A likely next blocker is an H1 consumer, but do not assume the label. Preserve
   the actual blocker.

Trust boundary:

- No H2 H1 materialization.
- No H1-J/density interaction.
- No private RHF for H2.
- No supplemented WL/QW comparison.
- No supplement support.
- No Be2/Cr2 artifact work.
- No HFDMRG, DMRG, ECP, exports, public solver behavior, or Qiu-White
  correction work.
- Do not add this test to default runners.
- Preserve the visible staged driver style.
- Do not request interactive approval or sandbox escalation. If approval is
  genuinely required, write `.agent_handoffs/ATTENTION.md` and stop.

Line-budget rule:

This pass must be net-negative for:

```text
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Pay for new code/test assertions by shrinking stale readiness/scaffold code or
stale tests. Do not delete scientific endpoints, reference tests, explicit He
driver endpoints, or the H2 readiness test.

Validation:

```text
julia --project=. test/nested/cartesian_ham_builder_h2_pqs_q5_gausslet_only_r4_readiness_runtests.jl
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
git status --short --branch
```

If the H2 readiness test exceeds 60 seconds, report elapsed time and whether it
is compilation, parent construction, or final-basis construction dominated.

Report back:

- final-basis request surface used;
- final-basis status, dimension, and overlap identity diagnostic if available;
- route/source-plan/H1 readiness status after the pass;
- scoped line-budget arithmetic;
- validation results and timing;
- deletion/shrinkage report with exact remaining caller/blocker.

-- repo-manager@macmini

Pass 197 - materialize H2 H1 only.

Purpose:

The H2 gausslet-only driver path now has real parent axes, an available source
plan, and a materialized final basis:

```text
final_dimension = 221
final_overlap_identity_error ~= 2.51e-13
```

The current blocker is:

```text
:missing_diatomic_complete_core_shell_h1_consumer
```

This pass should materialize H2 H1 through the existing diatomic
complete-core/shell H1 payload path, and write/read only compact H1 diagnostics
from the driver artifact.

Task type:

Narrow implementation plus focused test update.

Exact task:

1. Update the H2 driver input to request H1:

   ```julia
   run_final_basis = true
   run_h1 = true
   run_h1_j = false
   run_private_rhf = false
   ```

2. Reuse the existing diatomic complete-core/shell H1 payload path.

   Do not add a second H1 implementation path.

3. Extend the H2 artifact/test to assert compact H1 diagnostics:

   Required:

   ```text
   route/source_plan_status == :available_pqs_diatomic_complete_core_shell_source_plan
   route/final_basis_status == :available_pqs_complete_core_shell_final_basis
   route/h1_status indicates available/materialized H1
   basis/final_dimension == 221
   basis/final_overlap_identity_error is finite and small
   physics/h1_lowest is present and finite
   route/h1_j_materialized == false
   private_rhf/materialized == false
   comparison/ready == false
   ```

   Also record/check a scalar H1 symmetry/finiteness diagnostic if the existing
   H1 payload exposes one. Do not save a broad matrix unless already part of the
   artifact contract.

4. Report the H1 lowest energy value, but do not compare it to the old
   supplemented WL/QW H2 HF/ED references. This is still gausslet-only PQS.

5. Preserve the actual next blocker after H1 materializes. A likely next
   blocker is H1-J/density interaction, but do not assume the label.

Trust boundary:

- No H1-J/density interaction.
- No private RHF for H2.
- No supplemented WL/QW comparison.
- No supplement support.
- No Be2/Cr2 artifact work.
- No HFDMRG, DMRG, ECP, exports, public solver behavior, or Qiu-White
  correction work.
- Do not add this test to default runners.
- Do not store final-basis self-overlap as downstream data; scalar identity
  diagnostics only.
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
is compilation, parent construction, final-basis construction, or H1
construction dominated.

Report back:

- H1 request surface used;
- H1 status and lowest energy;
- H1 scalar symmetry/finiteness diagnostics, if available;
- route/source-plan/final-basis/H1-J readiness status after the pass;
- scoped line-budget arithmetic;
- validation results and timing;
- deletion/shrinkage report with exact remaining caller/blocker.

-- repo-manager@macmini

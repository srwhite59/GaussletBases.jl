Pass 199 - label current H2 221 route as diagnostic-only.

Purpose:

Pass 198 proved that the current H2 221-dimensional route is not the intended
physical atom-core-plus-shell H2 target. It retains only atom-side boundary
modes plus a midpoint/product slab:

```text
98 + 98 + 25 = 221
```

Before any future H2 work, the artifact should say this plainly. This pass
should label the current route as a source-box diagnostic / route smoke, not a
physics endpoint.

Task type:

Small implementation plus focused test update. No new physics.

Exact task:

1. Add compact H2 artifact labels for the current target role.

   Suggested meaning:

   ```text
   route/artifact_role = :source_box_diagnostic
   physics/endpoint_ready = false
   physics/endpoint_blocker = :retained_atom_core_interiors_missing
   basis/retained_atom_core_interiors = false
   basis/source_plan_role = :boundary_source_box_diagnostic
   ```

   Use names that fit existing artifact conventions, but keep the payload
   compact. Do not add a large inventory table.

2. Add or set the corresponding H2 driver input keys if useful:

   ```julia
   artifact_role = :source_box_diagnostic
   physics_endpoint_ready = false
   physics_endpoint_blocker = :retained_atom_core_interiors_missing
   ```

3. Update the explicit H2 readiness test to assert the labels.

   The test should still assert:

   ```text
   comparison/ready == false
   route/h1_j_materialized == false
   private_rhf/materialized == false
   ```

4. Do not change the current 221-dimensional route construction in this pass.

5. Delete or shrink stale scaffold code/tests to keep the scoped line budget
   negative.

Trust boundary:

- No H2 source-plan/final-basis redesign.
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

Pay for new labels/assertions by shrinking stale readiness/scaffold code or
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
is compilation dominated.

Report back:

- exact artifact labels added;
- confirmation current 221 route remains diagnostic-only;
- scoped line-budget arithmetic;
- validation results and timing;
- deletion/shrinkage report with exact remaining caller/blocker.

-- repo-manager@macmini

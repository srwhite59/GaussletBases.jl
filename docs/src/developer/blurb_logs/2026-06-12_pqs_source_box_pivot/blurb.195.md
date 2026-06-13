Pass 195 - materialize the H2 diatomic parent/axis bundle only.

Purpose:

The H2 gausslet-only driver artifact is honest but still blocked before a real
source plan:

```text
:missing_diatomic_complete_core_shell_source_plan_producer
```

Pass 194 diagnosed the lower dependency: the current H2 run still has no real
diatomic parent basis/axis-bundle objects. Existing code already has a reviewed
probe path that can construct and carry bond-aligned diatomic parent objects:

```text
_pqs_explicit_core_spacing_parent_axis_probe(...; carry_objects = true)
bond_aligned_homonuclear_qw_basis
_qwrg_bond_aligned_axis_bundles
```

This pass should connect only that parent/axis-bundle materialization seam for
the H2 driver readiness path. Do not implement the source-plan producer yet.

Task type:

Narrow implementation plus focused test update.

Read/inspect:

```text
src/pqs_source_box_route_driver_helpers.jl
src/pqs_source_box_diatomic_complete_core_shell.jl
src/pqs_source_box_route_driver_reporting.jl
test/driver_inputs/h2_pqs_q5_gausslet_only_r4.jl
test/nested/cartesian_ham_builder_h2_pqs_q5_gausslet_only_r4_readiness_runtests.jl
test/nested/pqs_explicit_core_spacing_parent_axis_probe_runtests.jl
test/nested/cartesian_parent_contract_runtests.jl
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/review.194.md
```

Known relevant existing behavior:

- `_cartesian_parent_object_carry` already consumes a parent axis probe with
  carried `basis_object` and `axis_bundle_object`.
- `pqs_explicit_core_spacing_parent_axis_probe_runtests.jl` already proves
  `carry_objects = true` can construct:

  ```text
  basis_object_type_label == "BondAlignedDiatomicQWBasis3D"
  axis_bundle_object_type_label == "_CartesianNestedAxisBundles3D"
  ```

- `cartesian_parent_contract_runtests.jl` already shows a probed Be2 parent
  can report `parent_axis_bundle_object_available = true`.

Exact task:

1. Enable the H2 driver readiness path to request/carry the real diatomic
   parent basis and axis-bundle objects.

   Prefer reusing the existing explicit-core-spacing parent-axis probe path.
   Do not add a second parent constructor path if the existing probe can do the
   job.

   If this can be achieved by input/config changes plus a small driver
   connection, prefer that over a new helper.

2. Update the H2 readiness artifact/test to assert the parent/axis-bundle step
   is now real.

   Expected artifact/report signals should include, where currently available:

   ```text
   parent_basis_object_available = true
   parent_axis_bundle_object_available = true
   parent_basis_object_type_label = "CartesianParentGaussletBasis3D" or wrapped parent label
   parent_qw_basis_object_type_label = "BondAlignedDiatomicQWBasis3D"
   parent_axis_bundle_object_type_label = "_CartesianNestedAxisBundles3D"
   parent_materialization_blocker === nothing
   ```

   Use the actual existing report/artifact key names. If the artifact does not
   currently expose these, add the smallest stable diagnostic keys under a
   `parent/*` or existing parent-report group.

3. The H2 source-plan blocker should move forward or become more precise.

   Acceptable outcomes:

   ```text
   A. source plan still blocked, but parent axis/bundle availability is true and
      the blocker is now a narrower source-realization/source-plan producer
      blocker.

   B. source plan unexpectedly becomes available through existing wiring.
      If so, keep final basis/H1/H1-J/RHF assertions out of this pass and report
      the new status for the next pass.
   ```

   If parent/axis-bundle construction fails or produces axis counts inconsistent
   with the H2 target, do not paper over it. Preserve the blocker and report the
   mismatch.

Trust boundary:

- Do not implement a diatomic source-plan producer in this pass.
- Do not materialize H2 final basis, H1, H1-J, density interaction, or private
  RHF by requirement.
- Do not compare to WL/QW supplemented HF/ED references.
- Do not add supplement support.
- Do not add Be2/Cr2 artifact work.
- Do not add HFDMRG, DMRG, ECP, exports, public solver behavior, or Qiu-White
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

Pay for any new code/test assertions by shrinking stale readiness/scaffold code
or stale tests. Do not delete scientific endpoints, reference tests, explicit
He driver endpoints, or the new H2 readiness test.

Validation:

```text
julia --project=. test/nested/cartesian_ham_builder_h2_pqs_q5_gausslet_only_r4_readiness_runtests.jl
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
git status --short --branch
```

If the H2 readiness test exceeds 60 seconds, report elapsed time and whether it
is compilation or parent construction dominated.

Report back:

- what parent/axis-bundle path was connected;
- parent availability/type-label fields verified;
- whether axis counts match the H2 target;
- new route/source-plan readiness status and blocker;
- scoped line-budget arithmetic;
- validation results;
- deletion/shrinkage report with exact remaining caller/blocker.

-- repo-manager@macmini

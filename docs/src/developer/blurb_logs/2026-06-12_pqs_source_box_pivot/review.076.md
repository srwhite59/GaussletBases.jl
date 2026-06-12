Review 076: accepted.

The pass adds the intended shellification/lowering-backed boundary:

- `PQSMultilayerShellLayerRegion`;
- `PQSMultilayerShellRegionPlan`;
- `pqs_multilayer_shell_region_plan(shellification_plan, lowering_plan; ...)`;
- `pqs_multilayer_shell_source_plan(bundles, region_plan; ...)`.

The new region plan carries geometry/lowering facts only: direct core region,
ordered complete-shell layer records, selected PQS lowering contracts/source
CPBs, coverage and duplicate/disjointness fingerprints, and shellification /
terminal-lowering provenance. It does not carry PQS descriptors, Lowdin data,
support operators, final-basis transfer, H1, IDA, density-density, RHF, driver
data, exports, or artifacts.

The new source-plan overload still delegates to the explicit-box bridge for
descriptor/Lowdin source realization. That is acceptable for this pass because
the route authority for the boxes is now the region plan. The old explicit-box
entry point remains working as a bridge.

Manager validation:

- `julia --project=. -e 't = @elapsed include("test/nested/pqs_direct_retained_final_h1_runtests.jl"); println("elapsed_s=", t)'`
  passed, 44/44, elapsed about 6.01s.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.

Follow-up:

- The H1 gate now compares region-plan and explicit-box source plans, but the
  active H1 operator path still uses the explicit-box `fixture.plan`.
- Next pass should make the region-plan source plan the active H1 fixture path
  and keep the explicit-box path only as a compact bridge comparison.
- Keep support-space dense operator helpers bounded as H1 seam/oracle machinery;
  do not generalize them into scalable PQS operator construction here.

Deletion/shrinkage:

- Geometry authority moved out of ad hoc test/probe box arithmetic into the new
  region-plan object.
- No deletion yet because the explicit-box entry point is still the source
  realization bridge and was explicitly kept.
- The next shrinkage target is the test fixture's active use of explicit
  `core_box` / `outer_box` planning.

-- repo-manager@macmini

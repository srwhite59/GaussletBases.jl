Pass 076 complete.

Object/helper added:

- Added `PQSMultilayerShellLayerRegion`.
- Added `PQSMultilayerShellRegionPlan`.
- Added `pqs_multilayer_shell_region_plan(shellification_plan, lowering_plan; metadata = (;))`.
- Added `pqs_multilayer_shell_source_plan(bundles, region_plan; ...)`.

Ownership boundary:

- `pqs_multilayer_shell_region_plan(...)` owns only shellification/lowering facts:
  - direct core terminal region and core box;
  - ordered complete-shell layer records;
  - shell layer outer/inner boxes;
  - selected PQS lowering contracts;
  - selected source CPBs;
  - support coverage and duplicate/disjointness fingerprints;
  - provenance from `CartesianShellification` and `CartesianTerminalLowering`.
- It does not build PQS descriptors, Lowdin matrices, support operator blocks, final-basis transfer, H1, IDA, density-density, RHF, driver data, exports, or artifacts.
- The new `pqs_multilayer_shell_source_plan(bundles, region_plan; ...)` overload consumes that boundary and delegates through the existing explicit-box bridge for now. The geometry authority is now the region plan, while the old explicit-box entry point remains available as a bridge.

Comparison result versus explicit-box bridge:

- Extended `test/nested/pqs_direct_retained_final_h1_runtests.jl` with a compact comparison using the existing H1 fixture geometry:
  - parent count 7;
  - shellification center at `(4, 4, 4)`;
  - `OneCenterShellification(core_side = 5, q = 5)`;
  - `PQSLowering(q = 5)`;
  - direct core box `(2:6, 2:6, 2:6)`;
  - outer/current box `(1:7, 1:7, 1:7)`.
- The region-plan entry point matches the explicit-box bridge on:
  - source-plan availability;
  - layer count;
  - core support count;
  - shell support count;
  - shell final retained count;
  - complete final retained count;
  - final overlap identity check.

Validation run:

- `julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl`
  - passed: 44 tests.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed.
- `git diff --check`
  - passed.

Deletion/shrinkage report:

- Responsibility moved out of ad hoc PQS box arithmetic:
  - the new region plan now records the direct core and ordered complete-shell regions from shellification/lowering;
  - support coverage and duplicate/disjointness fingerprints are now available before PQS descriptor/Lowdin realization;
  - selected PQS source CPBs are carried from terminal lowering rather than inferred only from explicit source-plan boxes.
- Explicit-box/probe code remaining as bridge only:
  - `pqs_multilayer_shell_source_plan(bundles, core_box, outer_box; ...)` still performs the actual descriptor/Lowdin source realization and remains the compatibility bridge for current probes/tests;
  - the new region-plan overload delegates to that bridge for this pass.
- Tests:
  - no new broad test file was added;
  - the existing H1 gate was minimally extended with compact region-plan-vs-explicit bridge checks.
- Nothing was deleted in this pass because the explicit-box entry point is still the active realization bridge and the blurb explicitly asked to keep it working.
- Remaining stale/duplicate surface to retire next:
  - once the region-plan entry point owns more of the flow, the explicit-box bridge should be quarantined as legacy/probe-only or reduced to a compatibility wrapper;
  - box-depth and layer-box arithmetic inside the explicit-box implementation remains duplicated geometry authority until that follow-up.

-- repo-doer@macmini

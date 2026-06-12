Pass 078 complete.

Private realization helper:

- Added `_pqs_multilayer_realize_shell_source_plan(...)`.
- It consumes:
  - `bundles`;
  - `core_box`;
  - `outer_box`;
  - ordered layer specs with `current_box`, `inner_box`, and provenance;
  - `bond_axis`, `term_coefficients`, `source_kind`, and `metadata`.
- It owns the common source realization work:
  - support index/state construction;
  - projected-q shell layer construction;
  - descriptor extraction;
  - shell realization plan and Lowdin cleanup;
  - collapsed shell sector assembly;
  - compact source-plan summary.

Entry point routing:

- Explicit-box bridge:
  - `pqs_multilayer_shell_source_plan(bundles, core_box, outer_box; ...)`
  - now builds layer specs with `_pqs_multilayer_explicit_box_layer_specs(...)`;
  - then calls `_pqs_multilayer_realize_shell_source_plan(...)`.
- Region-plan path:
  - `pqs_multilayer_shell_source_plan(bundles, region_plan; ...)`
  - now builds layer specs from `region_plan.shell_layers` with `_pqs_multilayer_region_plan_layer_specs(...)`;
  - then calls `_pqs_multilayer_realize_shell_source_plan(...)` directly.

The region-plan entry point no longer calls the explicit-box public entry point.

Validation run:

- `julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl`
  - passed: 44 tests.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed.
- `git diff --check`
  - passed.

Deletion/shrinkage report:

- Explicit-box responsibility isolated as bridge-only:
  - `_pqs_multilayer_box_depth(...)` and `_pqs_multilayer_core_box_at_depth(...)` now feed only `_pqs_multilayer_explicit_box_layer_specs(...)`;
  - the region-plan path gets layer boxes from shellification/lowering-backed `region_plan.shell_layers`.
- Duplicate arithmetic remaining:
  - the explicit-box bridge still constructs layer specs from `core_box`/`outer_box` for compatibility;
  - the common realization helper still performs support coverage and duplicate checks after realization, because both entry points need the same source-plan validation.
- Test code:
  - no test code changed in this pass;
  - the existing H1 gate from pass 077 continued to validate the active region-plan path.
- Remaining stale/duplicate surface to retire next:
  - the explicit-box bridge remains available and still carries compatibility box arithmetic;
  - a later pass can quarantine it as legacy/probe-only once no active fixture needs to construct plans from explicit boxes.

-- repo-doer@macmini

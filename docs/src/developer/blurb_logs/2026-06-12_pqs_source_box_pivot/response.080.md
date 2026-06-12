Pass 080 complete.

Labels/docs changed:

- Updated the explicit-box `pqs_multilayer_shell_source_plan(bundles, core_box, outer_box; ...)` docstring to call it a compatibility/probe bridge.
- The docstring now directs new route work to the shellification/lowering-backed `pqs_multilayer_shell_source_plan(bundles, region_plan; ...)` entry point.
- Added compact source-plan summary labels:
  - `explicit_box_bridge`;
  - `shellification_backed_geometry`.
- The explicit-box entry point now passes:
  - `explicit_box_bridge = true`;
  - `shellification_backed_geometry = false`.
- Updated `docs/src/developer/pqs_near_term_final_basis_realization_plan.md` so the explicit-box entry point is described as compatibility/probe machinery, not route shellification authority.
- The docs now list `pqs_multilayer_shell_region_plan(...)` plus the region-plan source entry point as the driver-facing source-plan seam.

`source_kind` decision:

- Left `source_kind = :repeated_one_cell_projected_q_shell_layers` stable for the explicit-box bridge.
- Reason: no live source/test assertion depends on the old symbol, but historical tracked responses mention it. Metadata/summary labels were enough to downgrade authority without unnecessary churn.
- The region-plan path remains marked as `:shellification_backed_repeated_one_cell_projected_q_shell_layers`.

Validation run:

- `git diff --check`
  - passed.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed.

Deletion/shrinkage report:

- Conceptual authority downgraded:
  - explicit boxes are now labeled compatibility/probe bridge machinery;
  - shellification/lowering-backed region plans are documented as the route authority for new work.
- Explicit-box bridge remaining:
  - the public explicit-box entry point still exists for compatibility/probes;
  - it still carries explicit layer-spec construction internally, but is labeled as bridge-only.
- Test coverage:
  - no tests were changed or added;
  - no broad metadata assertion was added for the new labels.
- Remaining stale surface:
  - any future active caller using explicit boxes for new route behavior should move to `pqs_multilayer_shell_region_plan(...)` plus the region-plan source entry point.

-- repo-doer@macmini

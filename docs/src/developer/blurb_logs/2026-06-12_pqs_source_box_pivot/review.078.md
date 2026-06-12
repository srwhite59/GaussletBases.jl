Review 078: accepted.

The refactor extracts `_pqs_multilayer_realize_shell_source_plan(...)` as the
shared private source-realization helper. The explicit-box bridge now builds
layer specs with `_pqs_multilayer_explicit_box_layer_specs(...)` and calls that
helper. The region-plan entry point builds layer specs from
`region_plan.shell_layers` with `_pqs_multilayer_region_plan_layer_specs(...)`
and calls the same helper directly.

This satisfies the main requirement: the region-plan source path no longer
calls the public explicit-box entry point. Explicit `core_box` / `outer_box`
arithmetic remains isolated to the explicit-box bridge's layer-spec producer.

Manager validation:

- `julia --project=. -e 't = @elapsed include("test/nested/pqs_direct_retained_final_h1_runtests.jl"); println("elapsed_s=", t)'`
  passed, 44/44, elapsed about 6.26s.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.

Deletion/shrinkage:

- explicit-box box-depth/layer-box arithmetic is now bridge-local rather than
  on the region-backed path;
- no tests changed, because pass 077 already made the H1 gate exercise the
  region-backed route authority;
- the explicit-box bridge remains for compatibility and probe/oracle use.

Follow-up:

- A later pass can document or enforce bridge-only status for the explicit-box
  entry point if active callers no longer need it.
- Do not grow dense support-space one-body helpers into a general PQS operator
  algorithm without a separate scale/design pass.

-- repo-manager@macmini

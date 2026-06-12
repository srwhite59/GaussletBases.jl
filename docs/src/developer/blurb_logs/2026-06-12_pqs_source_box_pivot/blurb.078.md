Purpose:
  Continue shrinking explicit-box authority in multi-layer PQS source planning.
  Pass 077 made the H1 gate active path region-plan-backed, but the
  region-plan overload still delegates to the public explicit-box bridge:

  ```julia
  pqs_multilayer_shell_source_plan(bundles, region_plan; ...)
      -> pqs_multilayer_shell_source_plan(bundles, core_box, outer_box; ...)
  ```

  That keeps box-depth/layer-box arithmetic on the active internal path.

Task:
  Refactor `src/pqs_multilayer_shell_source_plan.jl` so explicit-box and
  region-plan entry points share a lower internal source-realization helper.

Implementation guidance:
  - Extract the common realization work into a private helper that consumes:
      - `bundles`;
      - `core_box`;
      - `outer_box`;
      - an ordered set of layer specifications containing `current_box`,
        `inner_box`, and optional provenance;
      - `bond_axis`, `term_coefficients`, and `metadata`.
  - The explicit-box entry point may remain responsible for building layer
    specifications from `_pqs_multilayer_box_depth(...)` and
    `_pqs_multilayer_core_box_at_depth(...)`.
  - The region-plan entry point should build/pass layer specifications from
    `region_plan.shell_layers`, not call the explicit-box public entry point.
  - Preserve current output shape and numerical behavior as much as possible.
  - Keep the explicit-box entry point working as bridge/compatibility.

Do not:
  - add new physics, H1/RHF/IDA/density-density, fixture-rule policy, driver
    wiring, exports, or artifacts;
  - change shellification or terminal lowering APIs unless absolutely needed;
  - add broad tests or metadata vocabulary assertions;
  - generalize support-space dense one-body helpers;
  - request UI escalation. In unattended baton mode, write
    `.agent_handoffs/ATTENTION.md` and stop if permission is genuinely needed.

Validation:
  - focused H1 gate;
  - load check;
  - `git diff --check`.

Report:
  - name of the private realization helper and how each entry point reaches it;
  - confirmation that region-plan entry no longer calls the explicit-box public
    entry point;
  - validation run;
  - deletion/shrinkage report:
      - what explicit-box responsibility was isolated as bridge-only;
      - what duplicate arithmetic remains;
      - whether any test code changed or stayed untouched.

-- repo-manager@macmini

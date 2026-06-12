Accepted.

This pass correctly stopped at the missing input boundary instead of hard-coding
the focused H1 fixture into the driver. The result is useful: the one-center
driver path can construct the parent axis bundle, but the active
`:pqs_source_box` route does not yet select terminal shellification and does
not carry the complete core/shell region/source/final-basis inputs needed to
feed the H1/J diagnostic slot.

The next implementation boundary is now sharper:

- shell stage must select/own the `CartesianShellification.shellify(...)`
  output for the complete core/shell PQS route;
- lowering/transform stages must own the `PQSLowering` plan,
  `pqs_multilayer_shell_region_plan(...)`, and
  `pqs_multilayer_shell_source_plan(...)`;
- parent or transform stage must surface route-owned Coulomb expansion, axis
  weights, and raw pair numerator terms for the H1/J diagnostic;
- assembly can then build final basis, H1 payload, and H1/J payload.

No production code or tests changed. That is appropriate for this audit-only
boundary pass.

Test-bloat guidance for follow-up:

- Do not grow `test/nested/pqs_direct_retained_final_h1_runtests.jl` into a PQS
  notebook. It is already carrying region-plan, explicit-box bridge, final
  basis, H1, axis-layer convention, fixed-block oracle, and nonclaim checks.
- Do not use the broad CPBM contract file, broad route-driver report file, or
  assembly/report low-order policy files as routine per-pass validation.
- A new PQS test should replace or shrink existing pressure, not add another
  layer of helper-vocabulary assertions.

Deletion/shrinkage review:

- Nothing became removable in this pass because the driver-owned input producer
  does not exist yet.
- No test was added, which is correct.
- The stale surface to retire later remains fixture/probe H1/J construction
  glue once the driver owns the source-plan/final-basis/H1/J input producer.

Next blurb should target the first missing owner: shell-stage selection of the
complete core/shell PQS shellification plan, with no final-basis/H1/J expansion
unless the shell/lowering seam is already present.

-- repo-manager@macmini

Pass 205 review - accepted

Reviewed the low-order/WL materialization extraction.

Result:
- Accepted.
- `src/pqs_source_box_route_driver_helpers.jl` no longer owns the
  route-configured White-Lindsey / low-order materialization implementation.
- The implementation now lives in private include:
  `src/pqs_source_box_low_order_materialization.jl`.
- `src/GaussletBases.jl` includes the new file immediately after
  `pqs_source_box_route_driver_helpers.jl`.

What changed:
- Moved the low-order/WL preflight, one-center/diatomic materializer probes,
  atom-growth materializer probe, basis/Ham adapters, low-order shellization
  policy helper, and `_pqs_source_box_route_driver_materialization(...)`.
- No intended behavior changes.
- No public API changes.
- No report field/return-shape changes.
- No H1/H1-J/RHF, artifact schema, GTO, IDA/MWG, export, parent, shell, unit,
  or transform behavior changes.

Line counts:
```text
src/pqs_source_box_route_driver_helpers.jl      8973
src/pqs_source_box_low_order_materialization.jl 2655
src/GaussletBases.jl                             715
```

Line budget:
- `src/GaussletBases.jl`: +1
- `src/pqs_source_box_route_driver_helpers.jl`: -2657
- new file counted separately: +2655
- scoped source/test/bin net: -1

Validation run by manager:
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed.
- `julia --project=. bin/cartesian_ham_builder.jl save_artifact=false save_tsv=false`
  - passed; printed `driver complete`.
  - This covered the default save-disabled driver materialization/report path.
- `git diff --check`
  - passed.

`:legacy_diatomic_source` audit:
- Still live.
- It remains the default low-order shellization policy when
  `low_order_shellization_policy = nothing`.
- It is still consumed/reported by focused low-order stage/report/materializer
  tests and by the default driver smoke.
- Not deletion-ready in this pass.

Deletion/shrinkage assessment:
- deleted: no behavior path deleted.
- simplified: the driver helper file no longer contains the low-order/WL
  materialization subsystem.
- quarantined: low-order/WL materialization is now an explicitly private
  include.
- not deleted because: `:legacy_diatomic_source` is still default policy and
  live in focused tests/default driver behavior.
- exact remaining caller/blocker: `cartesian_materialization(...)` still calls
  `_pqs_source_box_route_driver_materialization(...)`; the implementation is
  now private, but the low-order stage/report tests and default driver path
  still preserve `:legacy_diatomic_source` metadata.

-- repo-manager@macmini

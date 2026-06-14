Pass 205 complete.

Files changed:
- `src/GaussletBases.jl`
- `src/pqs_source_box_route_driver_helpers.jl`
- `src/pqs_source_box_low_order_materialization.jl`

What changed:
- Extracted the route-configured White-Lindsey / low-order materialization implementation from `src/pqs_source_box_route_driver_helpers.jl` into private include file `src/pqs_source_box_low_order_materialization.jl`.
- Added `include("pqs_source_box_low_order_materialization.jl")` immediately after `include("pqs_source_box_route_driver_helpers.jl")`.
- Left the public/staged driver surface unchanged. `cartesian_materialization(...)` still reaches `_pqs_source_box_route_driver_materialization(...)`; the implementation now lives in the private low-order materialization file.

Moved helper families:
- `_pqs_source_box_route_driver_white_lindsey_preflight_fixed_block`
- `_pqs_source_box_route_driver_white_lindsey_ham_preflight`
- `_pqs_source_box_route_driver_one_center_materializer_probe`
- `_pqs_source_box_route_driver_route_configured_one_center_report`
- `_pqs_source_box_route_driver_diatomic_materializer_probe`
- `_pqs_source_box_route_driver_diatomic_atom_growth_materializer_probe`
- `_pqs_source_box_route_driver_route_configured_diatomic_basis_adapter`
- `_pqs_source_box_route_driver_route_configured_diatomic_basis_adapter_summary`
- `_pqs_source_box_route_driver_diatomic_atom_growth_basis_adapter`
- `_pqs_source_box_route_driver_route_configured_diatomic_ham_adapter`
- `_pqs_source_box_route_driver_route_configured_diatomic_ham_adapter_summary`
- `_pqs_source_box_route_driver_low_order_shellization_policy`
- `_pqs_source_box_route_driver_materialization`

Behavior/report shape:
- Behavior intended unchanged.
- Report-facing fields, aliases, statuses, and materialization return shapes were not renamed or removed.
- No H1/H1-J/RHF, artifact behavior, GTO, IDA/MWG, exports, parent/shell/unit/transform splitting, or lattice-size changes were added.

`:legacy_diatomic_source` audit:
- Still live and not deleted.
- It remains the default resolved low-order shellization policy for `low_order_shellization_policy = nothing` in the moved `_pqs_source_box_route_driver_low_order_shellization_policy`.
- It is still used by the source helper for low-order shell/unit/transform/pair/assembly/report metadata and by focused tests such as:
  - `test/nested/cartesian_shell_stage_low_order_policy_runtests.jl`
  - `test/nested/cartesian_unit_stage_low_order_policy_runtests.jl`
  - `test/nested/cartesian_transform_stage_low_order_policy_runtests.jl`
  - `test/nested/cartesian_pair_stage_low_order_policy_runtests.jl`
  - `test/nested/cartesian_assembly_stage_low_order_policy_runtests.jl`
  - `test/nested/cartesian_report_stage_low_order_policy_runtests.jl`
  - `test/nested/cartesian_route_diatomic_materializer_probe_runtests.jl`
- The focused driver smoke also reported `materialization.low_order_shellization_policy_resolved = :legacy_diatomic_source` on the default save-disabled path.

Line/carrying-cost accounting:
- `src/pqs_source_box_route_driver_helpers.jl`: 11630 -> 8973 lines.
- `src/pqs_source_box_low_order_materialization.jl`: 2655 lines.
- `src/GaussletBases.jl`: 714 -> 715 lines.
- Scoped real line count: 12344 -> 12343, net -1.
- Tracked diff before staging: `src/GaussletBases.jl` +1, `src/pqs_source_box_route_driver_helpers.jl` -2657; new file is untracked and 2655 lines.

Validation:
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - Passed. It precompiled `GaussletBases` in about 57 seconds, then printed `load ok`.
- `julia --project=. bin/cartesian_ham_builder.jl save_artifact=false save_tsv=false`
  - Passed; printed `driver complete`.
  - This covered the default materialization/report path with artifact and TSV saving disabled.
- `git diff --check`
  - Passed.
- `git diff --cached --check`
  - Passed; no staged changes.

Git status:
```text
## main...origin/main
 M src/GaussletBases.jl
 M src/pqs_source_box_route_driver_helpers.jl
?? src/pqs_source_box_low_order_materialization.jl
```

Deletion/shrinkage report:
- deleted: no behavior path deleted; removed the low-order/WL materialization implementation block from `src/pqs_source_box_route_driver_helpers.jl`.
- simplified: route-driver helpers no longer own the low-order/WL materializer implementation; that implementation is isolated in private `src/pqs_source_box_low_order_materialization.jl`.
- quarantined: low-order/WL materialization remains private route support; `:legacy_diatomic_source` remains explicitly live and unpromoted.
- not deleted because: `:legacy_diatomic_source` is still the default low-order policy and is asserted by focused stage/report/materializer tests.
- exact remaining caller/blocker: `cartesian_materialization(...)` still depends on `_pqs_source_box_route_driver_materialization(...)`; the low-order stage/report tests and default driver path still consume/report `:legacy_diatomic_source` metadata.

-- repo-doer@macmini

Pass 179 response - fixed-q multi-shell PQS source modes for one-center He

Files changed:
- `src/cartesian_terminal_lowering/region_contracts.jl`
- `src/pqs_multilayer_shell_source_plan.jl`
- `test/nested/pqs_direct_retained_final_h1_runtests.jl`

What changed:
- `_pqs_complete_shell_contract(region, policy::PQSLowering)` now records `source_mode_shape = (q, q, q)` from the lowering policy and keeps the physical support CPB shape separately as `source_box_shape`.
- Region-backed `pqs_multilayer_shell_source_plan(bundles, region_plan; ...)` now passes each layer's `source_mode_shape` through to realization.
- `_pqs_multilayer_realize_shell_source_plan(...)` now uses a layer-provided source-mode shape when present, requiring it to be cubic for this pass.
- Explicit-box bridge behavior is preserved: explicit-box specs still do not carry `source_mode_shape`, so they continue to derive source shape/q from `inner_box` as before.
- The focused PQS test was narrowed to the fixed-q source/final-basis inventory contract and no longer builds duplicate explicit-box comparison data or H1/nonclaim scaffolding.

Observed fixed-q source/final-basis inventory:

```text
parent/current box        = (1:11, 1:11, 1:11)
direct core box           = (4:8, 4:8, 4:8)
core_support_count        = 125
shell_layer_count         = 3
shell_support_count       = 1206
shell raw_source_dims     = (5, 5, 5) for all three shell records
shell retained_count      = 98 for all three shell records
source_mode_shape_source  = :terminal_lowering_contract
shell_final_retained_count = 294
final_retained_count      = 419
final_overlap_identity_error < 1e-10
```

Explicit-box bridge behavior:
- Preserved by construction. `_pqs_multilayer_explicit_box_layer_specs(...)` still does not set `source_mode_shape`, and `_pqs_multilayer_realize_shell_source_plan(...)` falls back to `length.(inner_box)` when `source_mode_shape` is absent.
- The shellification/lowering-backed path now records `fixed_source_mode_shape_used = true` and `source_mode_shape_sources = (:terminal_lowering_contract,)`.

Source/test line budget:

```text
git diff --numstat -- src test
2       1       src/cartesian_terminal_lowering/region_contracts.jl
21      5       src/pqs_multilayer_shell_source_plan.jl
25      136     test/nested/pqs_direct_retained_final_h1_runtests.jl
```

Total added: 48
Total deleted: 142
Net: -94 lines

Validation:

```text
julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl
```

First run result:

```text
Precompiling packages...
  57430.8 ms  ✓ GaussletBases
  1 dependency successfully precompiled in 57 seconds. 50 already precompiled.
diatomic.pqs_layer.lowdin_cleanup: 0.1215 seconds
Test Summary:                         | Pass  Total  Time
PQS complete core-shell final H1 gate |   22     22  4.1s
```

After renaming the testset label, rerun result:

```text
diatomic.pqs_layer.lowdin_cleanup: 0.1073 seconds
Test Summary:                                  | Pass  Total  Time
PQS fixed-q complete core-shell inventory gate |   22     22  4.1s
```

The >60s wall time on the first validation was package precompilation from the source edit; the test body itself reports 4.1s.

```text
julia --project=. -e 'using GaussletBases; println("load ok")'
```

Result:

```text
load ok
```

```text
git diff --check
```

Result: passed with no output.

```text
git diff --numstat -- src test
```

Result: net-negative by 94 lines, shown above.

Git status:

```text
## main...origin/main
 M src/cartesian_terminal_lowering/region_contracts.jl
 M src/pqs_multilayer_shell_source_plan.jl
 M test/nested/pqs_direct_retained_final_h1_runtests.jl
```

Deletion/shrinkage report:
- deleted: duplicate explicit-box bridge construction/comparison from `test/nested/pqs_direct_retained_final_h1_runtests.jl`; local fixed-block H1 oracle helper; H1 operator construction/assertions and repeated nonclaim assertions that were outside this pass's inventory target.
- simplified: focused test now protects the fixed-q He source/final-basis inventory directly: side 11 parent, three fixed `(5,5,5)` PQS shell records, final dimension 419.
- quarantined: H1/H1-J/RHF/physics energy interpretation remains out of scope for this pass; explicit-box bridge remains compatibility/probe behavior, not route authority.
- not deleted because: `pqs_multilayer_shell_source_plan(bundles, core_box, outer_box; ...)` still has live explicit-box compatibility/probe use and was intentionally preserved.
- exact remaining caller/blocker: next physics step is H1/H1-J on the 419-dimensional fixed-q He basis, then only after review RHF; AHGBS residuals and Fig.8 supplemented comparison remain separate future work.

-- repo-doer@macmini

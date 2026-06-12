Pass 092 response: PQS H1/J route-owned input audit

Summary:
- I did not feed the H1/J diagnostic payload slot yet.
- The route driver already can construct a parent axis bundle for the one-center probe, but the active `:pqs_source_box` dry-run still does not select terminal shellification or produce the complete core/shell region/source/final-basis inputs.
- I made no production, test, or docs edits in this pass.

Focused dry-run result:

```
status=blocked_missing_complete_core_shell_h1_j_route_inputs
blocker=missing_complete_core_shell_h1_j_route_inputs
missing=(:pqs_multilayer_shell_region_plan, :pqs_multilayer_shell_source_plan, :pqs_multilayer_complete_core_shell_final_basis, :pqs_multilayer_complete_core_shell_h1_payload, :axis_weights, :raw_pair_factor_terms, :coulomb_expansion)
low_order_shellization_source=not_applicable
low_order_terminal_selected=false
parent_axis_bundle_available=true
```

Exact missing input boundary by driver stage:
- Parent stage: `parent.parent_axis_bundle_object` is available in the one-center probe path.
- Parent density/Coulomb inputs: axis weights and raw pair numerator terms are derivable from the PGDG parent bundle, but they are not currently surfaced as a named route-owned H1/J density input. The Coulomb expansion is used in parent-axis construction, but it is not carried as a named route input for the complete core/shell H1/J payload.
- Shell stage: `:pqs_source_box` reports `low_order_shellization_source = :not_applicable` and `low_order_terminal_shellification_selected = false`. The driver does not yet own/select the `CartesianShellification.shellify(...)` output needed by the complete core/shell path.
- Lowering/region/source stage: because there is no selected shellification plan, the driver cannot produce the `PQSLowering` terminal lowering, `pqs_multilayer_shell_region_plan(...)`, or `pqs_multilayer_shell_source_plan(...)` from route-owned objects.
- Assembly stage: without the source plan, the driver cannot produce `pqs_multilayer_complete_core_shell_final_basis(...)` or `pqs_multilayer_complete_core_shell_h1_payload(...)`.

Conclusion:
- The next implementation boundary is not the pass-091 payload slot itself. It is a driver-owned PQS source-box shellification/lowering/complete-core-shell input producer.
- I did not synthesize the `test/nested/pqs_direct_retained_final_h1_runtests.jl` fixture inside the driver, did not use the explicit-box bridge, and did not add report fields.

Validation:
- Ran focused driver dry-run probe above.
- Ran `julia --project=. -e 'using GaussletBases; println("load ok")'`: passed.
- Ran `git diff --check`: passed.
- Did not run `test/nested/pqs_direct_retained_final_h1_runtests.jl` because no H1/final-basis construction code was touched.

Deletion/shrinkage report:
- No old code, tests, metadata, or compatibility paths became unnecessary in this audit-only pass.
- Nothing was deleted or simplified because the route-owned input seam is still missing; deleting diagnostics now would remove the only precise blocker signal.
- No new tests were added.
- Remaining stale/duplicate surface to retire later: once a driver-owned complete core/shell input producer exists, any fixture-only H1/J glue that duplicates `pqs_direct_retained_final_h1_runtests.jl` construction should either become reference-only or be removed.

-- repo-doer@macmini

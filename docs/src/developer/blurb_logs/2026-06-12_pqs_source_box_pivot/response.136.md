Pass 136 response - probe-enabled Be2 readiness fingerprint

Files changed:
- `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`

Commit SHA:
- Not committed; pass 136 did not explicitly instruct the doer to commit.

Implementation summary:
- Kept the existing default Be2/PQS fixture behavior unchanged.
- Added a keyword to the focused fixture helper:
  - `_pqs_be2_ham_payload_fingerprint_assembly(; probe_parent_axis_construction = false)`
- Added a second focused testset using:
  - `probe_parent_axis_construction = :auto`
- No source files were changed.
- No scalar report aliases, production behavior, final-basis/H1/H1-J, density interaction, RHF/SCF, WL payload, exports, or artifacts were added.

Default Be2 readiness:
```text
parent_axis_bundle_object_available = false
status = :blocked_diatomic_complete_core_shell_ham_readiness
blocker = :missing_diatomic_complete_core_shell_source_plan_producer
missing_objects includes :parent_axis_bundle_object
missing_objects includes :diatomic_complete_core_shell_source_plan_producer
```

Probe-enabled Be2 readiness:
```text
probe_parent_axis_construction = :auto
parent_axis_bundle_object_available = true
status = :blocked_diatomic_complete_core_shell_ham_readiness
blocker = :missing_diatomic_complete_core_shell_source_plan_producer
available_objects includes :parent_axis_bundle_object
missing_objects does not include :parent_axis_bundle_object
missing_objects includes :diatomic_complete_core_shell_source_plan_producer
```

Conclusion:
- The parent axis-bundle object blocker is already solvable by the existing probe-enabled route option.
- The remaining Be2/PQS readiness blocker is the missing route-owned diatomic complete-core/shell source-plan producer.
- Probe-enabled behavior was not made default.

Private Ham payload blocker:
- Confirmed unchanged for both default and probe-enabled fixtures:
```text
ham.status = :blocked_complete_core_shell_ham_payload
ham.blocker = :missing_complete_core_shell_ham_inputs
ham.missing_inputs = (
    :pqs_multilayer_complete_core_shell_final_basis,
    :pqs_multilayer_complete_core_shell_h1_payload,
    :pqs_complete_core_shell_final_one_electron_hamiltonian,
    :complete_core_shell_density_inputs,
    :complete_core_shell_h1_j_diagnostic_payload,
    :pqs_complete_core_shell_pre_final_density_interaction,
)
```

Validation commands/results:
- `julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
  - Passed.
  - Default fixture: `46 / 46`, `11.6s`.
  - Probe-enabled fixture: `23 / 23`, `16.6s`.
  - No precompilation-dominated over-60s run in the final validation.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - Passed: `load ok`
- `git diff --check`
  - Passed.

Git status:
```text
## main...origin/main
 M test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl
```

Deletion/shrinkage report:
- deleted: none.
- simplified: the focused Be2 readiness test now separates parent-axis-bundle availability from the source-plan producer blocker.
- quarantined: probe-enabled route behavior remains a test-only fingerprint; final-basis/H1/H1-J, density interaction, RHF/SCF, WL payloads, exports, artifacts, and public API remain out of scope.
- not deleted because: the default fixture still documents the unprobed route path and existing Ham blocker compatibility.
- exact remaining caller/blocker: Be2/PQS still lacks the route-owned diatomic complete-core/shell source-plan producer. With `probe_parent_axis_construction = :auto`, the parent axis-bundle object is available, so the source-plan producer is now the narrower remaining blocker before final-basis/H1/density-input/Ham materialization can proceed.

-- repo-doer@macmini

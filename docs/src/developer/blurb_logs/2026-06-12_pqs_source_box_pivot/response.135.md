Pass 135 response - diatomic PQS Ham readiness payload

Files changed:
- `src/pqs_source_box_route_driver_helpers.jl`
- `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`

Commit SHA:
- Not committed; pass 135 did not explicitly instruct the doer to commit.

Exact payload/helper names:
- Payload type: `_PQSDiatomicCompleteCoreShellHamReadinessPayload`
- Main helper: `_pqs_source_box_route_driver_diatomic_complete_core_shell_ham_readiness_payload`
- Compact summary helpers:
  - `_pqs_source_box_route_driver_diatomic_center_summary`
  - `_pqs_source_box_route_driver_diatomic_source_box_summary`
  - `_pqs_source_box_route_driver_diatomic_retained_unit_summary`
  - `_pqs_source_box_route_driver_diatomic_pair_inventory_summary`
- Assembly field: `diatomic_complete_core_shell_ham_readiness_payload`

Implementation summary:
- Added one private compact readiness payload next to the existing complete-core/shell Ham payload seam.
- Wired it into `cartesian_assembly(...)` as a single object field, not as scalar report aliases.
- The helper returns compact not-applicable payloads for non-PQS or non-diatomic routes.
- The Be2/PQS path now records that route skeleton/source-box facts exist, but the route-owned diatomic complete-core/shell source-plan producer does not.
- No final basis, H1, H1/J, density interaction, RHF/SCF, WL payload, exports, artifacts, or public API were added.

Observed Be2 readiness status/blocker:
```text
status = :blocked_diatomic_complete_core_shell_ham_readiness
blocker = :missing_diatomic_complete_core_shell_source_plan_producer
```

Compact summary facts asserted by the Be2 focused test:
- `route_family = :pqs_source_box`
- `system_classification = :bond_aligned_diatomic`
- `bond_axis = :x`
- `center_count = 2`
- `nuclear_charges = (4, 4)`
- `parent_axis_bundle_object_available = false`
- `source_box_count = 3`
- `source_box_keys = (:pqs_left, :product, :pqs_right)`
- `retained_unit_count = 3`
- `retained_unit_kinds = (:pqs, :pqs, :product_doside)`
- `pair_count = 6`
- `pair_family_counts = (pqs_pqs = 3, pqs_product = 2, product_pqs = 0, product_product = 1)`
- `available_objects` includes `:route_skeleton`, `:source_boxes`, `:retained_units`, and `:pair_inventory`
- `missing_objects` includes `:diatomic_complete_core_shell_source_plan_producer` and `:parent_axis_bundle_object`
- non-promotion flags remain false for public API, final basis, H1, H1/J, Ham payload materialization, RHF, exports, and artifacts.

Existing Ham payload blocker:
- Behavior preserved.
- Focused test still asserts:
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
  - Passed: `46 / 46`
  - Final rerun test-set time: `12.2s`
  - Total command exceeded 60 seconds because `GaussletBases` precompiled for about `57s`; precompilation dominated wall time.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - Passed: `load ok`
- `git diff --check`
  - Passed.

Git status:
```text
## main...origin/main
 M src/pqs_source_box_route_driver_helpers.jl
 M test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl
```

Deletion/shrinkage report:
- deleted: none.
- simplified: Be2/PQS Ham readiness is now one route-owned object instead of a test-only inference from scattered downstream blockers.
- quarantined: WL payload implementation, shell/support-row contraction authority, final-basis/H1/H1-J materialization, RHF/SCF, exports, artifacts, and public API remain out of scope.
- not deleted because: existing report and Ham blocker compatibility still need to remain stable while the route-owned diatomic source-plan producer is missing.
- exact remaining caller/blocker: Be2/PQS still lacks a route-owned diatomic complete-core/shell source-plan producer for the `(:pqs_left, :product, :pqs_right)` source-box route. Until that exists, the complete-core/shell final-basis, H1, density-input, H1/J, and Ham payload path must remain blocked.

-- repo-doer@macmini

Pass 148 response - audit diatomic H1 consumer seam

Scope:
- No edits to source, tests, or docs.
- No Julia commands were run; the blurb only required read-only inspection.
- Read-only commands used `sed`, `nl`, `rg`, and `git status --short --branch`.

Existing H1/support-one-body consumer contract:
- `pqs_multilayer_support_kinetic_matrix(plan)` is old one-center/multilayer-plan shaped:
  - requires `object_kind === :pqs_multilayer_shell_source_plan`;
  - requires `status === :available_pqs_multilayer_shell_source_plan`;
  - builds states as `vcat(plan.core_support_states, plan.shell_support_states)`;
  - consumes `plan.metrics.x/y/z.overlap` and `.kinetic`;
  - returns a dense support matrix in core-then-shell row order.
  - References: `src/pqs_multilayer_support_one_body.jl:196`.
- `pqs_multilayer_support_electron_nuclear_by_center_matrices(plan; ...)` has the same old object-kind/status guard, then:
  - requires nonempty center records;
  - consumes `coulomb_expansion.coefficients`;
  - uses the same `vcat(core_support_states, shell_support_states)` order;
  - builds negative unit-charge by-center attraction matrices;
  - records nuclear charge and center metadata but does not apply charges or sum centers.
  - References: `src/pqs_multilayer_support_one_body.jl:240`.
- `pqs_complete_core_shell_final_one_body_matrix(final_basis, support_operator; term, ...)` is already final-basis generic:
  - requires a `NamedTuple` final basis with `object_kind === :pqs_complete_core_shell_final_basis`;
  - requires `status === :available_pqs_complete_core_shell_final_basis`;
  - requires `final_basis_materialized`;
  - requires `support_row_order === :core_then_shell`;
  - requires a finite symmetric support matrix of size `(core_support_count + shell_support_count)^2`;
  - applies `transpose(final_basis.final_coefficients) * support_operator * final_basis.final_coefficients`.
  - References: `src/cartesian_final_basis_realization/pqs_complete_core_shell_final_basis.jl:268` and `src/cartesian_final_basis_realization/pqs_complete_core_shell_final_basis.jl:790`.
- `pqs_multilayer_complete_core_shell_h1_payload(plan; final_basis, ...)` is a coordinating old-plan helper:
  - repeats the old `:pqs_multilayer_shell_source_plan` and `:available_pqs_multilayer_shell_source_plan` guard;
  - requires available complete core/shell final basis;
  - calls the two support-one-body helpers;
  - transfers kinetic and by-center nuclear matrices to the final basis;
  - calls final one-electron Hamiltonian assembly and H1 solve;
  - reports H1 materialized while keeping IDA/density/RHF/GTO/export/artifact false.
  - References: `src/pqs_multilayer_complete_core_shell_h1.jl:112`.
- Final Hamiltonian assembly applies nuclear charges and sums centers only after transfer:
  - requires final by-center metadata with recorded charge, uncharged nuclear input, and unsummed centers;
  - adds `charge .* nuclear` to the kinetic matrix.
  - References: `src/cartesian_final_basis_realization/pqs_complete_core_shell_final_basis.jl:330` and `src/cartesian_final_basis_realization/pqs_complete_core_shell_final_basis.jl:344`.

Requirements already satisfied by the private diatomic objects:
- `_PQSDiatomicCompleteCoreShellSourcePlan` already carries the support-one-body raw ingredients:
  - `bundles`;
  - `metrics`;
  - `core_support_indices`;
  - `core_support_states`;
  - `shell_support_indices`;
  - `shell_support_states`;
  - `shell_final_coefficients`;
  - `support_order`;
  - route/convention metadata.
  - References: `src/pqs_source_box_route_driver_helpers.jl:11741`.
- The diatomic source plan is constructed from route-owned raw-box/source-realization data and uses the same computational support order needed by the final-basis transfer:
  - core support states first;
  - shell support states second, with left PQS rows then right PQS rows;
  - duplicate core/shell support is checked;
  - metrics are derived from the parent axis bundle;
  - convention label says `support_row_order = :core_product_then_shell_left_right_pqs`.
  - References: `src/pqs_source_box_route_driver_helpers.jl:12634` and `src/pqs_source_box_route_driver_helpers.jl:12737`.
- `_PQSDiatomicCompleteCoreShellFinalBasisPayload` already carries:
  - `source_plan`;
  - `source_plan_status`;
  - `final_basis`;
  - `final_basis_status`;
  - available/missing objects;
  - summary/metadata.
  - References: `src/pqs_source_box_route_driver_helpers.jl:11763`.
- The diatomic final basis is built directly from the private diatomic source plan and the lower final-basis realization helper:
  - overlap blocks are assembled from the same source-plan support states and metrics;
  - final basis records `support_row_order == :core_then_shell`;
  - available payload currently advances the blocker to `:diatomic_complete_core_shell_h1_consumer`.
  - References: `src/pqs_source_box_route_driver_helpers.jl:13066`, `src/pqs_source_box_route_driver_helpers.jl:13106`, and `src/pqs_source_box_route_driver_helpers.jl:13158`.
- Current focused Be2 fixture confirms the available diatomic final-basis facts:
  - source plan available;
  - final basis available;
  - `core_support_count == 25`;
  - `shell_support_count == 250`;
  - `shell_final_retained_count == 196`;
  - `final_retained_count == 221`;
  - readiness blocker is `:missing_diatomic_complete_core_shell_h1_consumer`.
  - References: `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl:600` and `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl:635`.

Missing or ambiguous requirements:
- The direct blockers are the old object-kind/status guards:
  - old support helpers require `:pqs_multilayer_shell_source_plan`;
  - old H1 helper requires `:pqs_multilayer_shell_source_plan`;
  - the new route intentionally returns `:pqs_diatomic_complete_core_shell_source_plan` and records `returns_pqs_multilayer_shell_source_plan = false`.
  - References: `src/pqs_multilayer_support_one_body.jl:197`, `src/pqs_multilayer_complete_core_shell_h1.jl:121`, and `src/pqs_source_box_route_driver_helpers.jl:12778`.
- The support row-order labels are compatible but not identical:
  - support operators are built over `core_support_states` then `shell_support_states`;
  - diatomic source-plan convention is more specific: `core_product_then_shell_left_right_pqs`;
  - final-basis transfer only accepts the compact final label `:core_then_shell`.
  - This should be asserted in the diatomic H1 seam rather than silently relying on label similarity.
- `bundles`/`metrics` are structurally compatible with the old support helpers for kinetic and centered Gaussian factors, but the old helper names would still report `pqs_multilayer_*` sources. A diatomic wrapper should either preserve distinct metadata or call the lower product/factor primitives directly.
- Actual nuclear center records are available from `parent.center_table`, but the diatomic source-plan payload currently carries only `center_summary`, not full records. The existing helper `_pqs_source_box_route_driver_complete_core_shell_center_records(parent)` already produces the needed records.
  - Reference: `src/pqs_source_box_route_driver_helpers.jl:10659`.
- Actual Coulomb expansion data are needed for nuclear support matrices. The diatomic source-plan payload currently carries `coulomb_expansion_summary`, not the expansion object. The next H1 pass can call `coulomb_gaussian_expansion(doacc = false)` at the H1 seam or add a compact private expansion field if manager approves.
- Axis layers are available from `parent.parent_axis_bundle_object` through `_pqs_source_box_route_driver_complete_core_shell_axis_layers(parent_axis_bundle_object)`.
  - Reference: `src/pqs_source_box_route_driver_helpers.jl:10677`.
- The electron-nuclear sign convention is already clear: support by-center matrices are negative unit-charge attraction, nuclear charges are recorded but not applied, and Hamiltonian assembly applies charges and sums centers.
  - References: `src/pqs_multilayer_support_one_body.jl:233` and `src/cartesian_final_basis_realization/pqs_complete_core_shell_final_basis.jl:353`.
- H1 can safely reuse the lower final-basis transfer, final Hamiltonian assembly, and H1 solve. It should not reuse the top-level old `pqs_multilayer_complete_core_shell_h1_payload` unchanged because that would require broadening the old object-kind contract and would blur the new private diatomic source-plan boundary.

Recommended smallest implementation pass:
- Add a private route-owned payload:
  - `_PQSDiatomicCompleteCoreShellH1Payload`
  - helper: `_pqs_source_box_route_driver_diatomic_complete_core_shell_h1_payload(...)`
- Wire only the Be2/PQS `cartesian_assembly` diagnostic/readiness path to build this payload after the diatomic final-basis payload.
- Keep the old `pqs_multilayer_support_kinetic_matrix`, `pqs_multilayer_support_electron_nuclear_by_center_matrices`, and `pqs_multilayer_complete_core_shell_h1_payload` guards unchanged in the first pass.
- Inside the new private helper, assert the diatomic source-plan object kind/status and final-basis status, then use private/local support-one-body wrappers over the existing lower primitives:
  - support kinetic from source-plan `metrics` and `vcat(core_support_states, shell_support_states)`;
  - support electron-nuclear by-center from parent center records, axis layers, Coulomb expansion coefficients, and the existing centered Gaussian factor helpers;
  - final one-body transfer through `CartesianFinalBasisRealization.pqs_complete_core_shell_final_one_body_matrix`;
  - final Hamiltonian and H1 solve through existing `CartesianFinalBasisRealization` helpers.
- Proposed report fields:
  - `status`;
  - `blocker`;
  - `route_family`;
  - `source_plan`;
  - `source_plan_status`;
  - `final_basis`;
  - `final_basis_status`;
  - `support_kinetic_status`;
  - `support_electron_nuclear_status`;
  - `final_kinetic_status`;
  - `final_electron_nuclear_status`;
  - `final_hamiltonian_status`;
  - `h1_status`;
  - `h1_payload`;
  - `final_dimension`;
  - `lowest_energy` only when H1 is naturally solved;
  - `available_objects`;
  - `missing_objects`;
  - `summary`;
  - `metadata`.
- Nonclaim flags in the summary/metadata:
  - `h1_materialized`;
  - `h1_j_materialized = false`;
  - `ham_payload_materialized = false` unless the next blurb explicitly promotes the private route payload into the Ham readiness report;
  - `density_density_materialized = false`;
  - `rhf_materialized = false`;
  - `public_api = false`;
  - `exports_materialized = false`;
  - `artifacts_materialized = false`;
  - `shell_support_row_contraction_authority = false`;
  - `retained_diagnostic_weights_are_ida_weights = false`.

Why not broaden the old guards first:
- The raw numerical ingredients are compatible, but the old helper names and object-kind contract encode `:pqs_multilayer_shell_source_plan`.
- The current pivot is deliberately retiring that adapter claim for the Be2/PQS path.
- A private diatomic H1 payload keeps source-box-first route authority distinct while still reusing lower matrix/transfer primitives.

Git status:
```text
## main...origin/main
```

Deletion/shrinkage forecast:
- deleted: none in this no-edit audit.
- simplified: next implementation can shrink the current Ham readiness blocker from missing generic old multilayer H1 payloads to a private diatomic H1 payload status.
- quarantined: old `:pqs_multilayer_shell_source_plan` claim; shell/support-row contraction; H1-J/density/RHF; WL payloads; public API/exports/artifacts/hfdmrg/CR2.
- not deleted because: old one-center H1 helpers remain the validated compact one-center fixture path, and the diatomic H1 consumer has not yet been implemented.
- exact remaining caller/blocker: `cartesian_assembly(...).diatomic_complete_core_shell_ham_readiness_payload` currently blocks on `:missing_diatomic_complete_core_shell_h1_consumer` after source plan and final basis are available.

-- repo-doer@macmini

Pass 130 response - no-edit Hamiltonian constructor payload audit

Scope:
- No source/test/doc edits.
- No Julia commands run; read-only inspection was sufficient.
- Inspected the requested route-driver, H1/H1J, RHF-as-validator, reporting/export, framework, fixture-policy, hfdmrg, and CR2 context surfaces.

Source inventory:
- `src/pqs_source_box_route_driver_helpers.jl`
  - `_pqs_source_box_route_driver_complete_core_shell_diagnostic_route_payload(...)` already builds the route-owned complete core/shell diagnostic bundle.
  - `_PQSCompleteCoreShellDiagnosticRoutePayload` carries `source_payload`, `final_basis`, `h1_payload`, `density_inputs`, and `h1_j_payload`; this is the right object to consume, not report aliases.
  - `_pqs_source_box_route_driver_complete_core_shell_density_inputs(...)` carries axis weights, raw pair-factor terms, Coulomb term count/provenance, and explicit diagnostic/non-authority flags.
  - `_pqs_source_box_route_driver_complete_core_shell_h1_j_diagnostic_payload(...)` materializes the H1/J diagnostic route payload from real objects and stores the nested `density_interaction`.
  - `_pqs_source_box_route_driver_complete_core_shell_h1_j_report_fields(...)` exposes only summary/status/scalar aliases and should not become the Ham authority boundary.
- `src/pqs_multilayer_complete_core_shell_h1.jl`
  - `pqs_multilayer_complete_core_shell_final_basis(...)` builds the final basis and preserves route/source-plan metadata while explicitly not materializing IDA, density-density product export, RHF, driver wiring, exports, or artifacts.
  - `pqs_multilayer_complete_core_shell_h1_payload(...)` builds final kinetic, separated by-center nuclear terms, `final_hamiltonian`, and H1 solve. The one-body Ham matrix is already present as `h1_payload.final_hamiltonian.hamiltonian_matrix`.
  - `pqs_multilayer_complete_core_shell_h1_j_payload(...)` builds support weights, support raw pair numerator, and `density_interaction`; its self-Coulomb is a diagnostic.
- `src/cartesian_final_basis_realization/pqs_complete_core_shell_final_basis.jl`
  - `pqs_complete_core_shell_final_one_electron_hamiltonian(...)` returns the one-body Hamiltonian object with kinetic matrix, charged nuclear matrix, center summaries, symmetry/finiteness diagnostics, and nonclaim flags.
  - `pqs_complete_core_shell_pre_final_density_interaction(...)` returns the available electron-electron representation today: a two-index pre-final density-interaction matrix plus final-to-pre-final coefficients, pre-final positive-weight diagnostics, support row ordering, raw-pair numerator shape/symmetry, and convention flags.
- `src/pqs_multilayer_support_density.jl`
  - `pqs_multilayer_support_weights(...)` and `pqs_multilayer_support_pair_raw_numerator_matrix(...)` define the support-row density input convention: core-then-shell support order, raw numerator, no density-normalized pair-factor authority.
- `src/pqs_multilayer_complete_core_shell_rhf.jl`
  - RHF/SCF consumes the H1 payload and the pre-final density interaction as a private validator. It should not define the product surface.
- `src/pqs_source_box_route_driver_reporting.jl`
  - Current materialization/report fields around `ham_preflight_status`, `ham_operator_payload_status`, `ham_interaction_status`, `ham_bundle_export_status`, and artifact fields are broad private driver/report vocabulary, not a clean PQS Ham payload.
  - Durable report serialization elides heavy objects such as `route_materializer_payload`; this reinforces keeping the first Ham payload private and object-carrying, not artifact-facing.
- `src/fullida_dense_export.jl`
  - Existing dense export shape is `fullida_dense_v1`: dense `H1`, dense two-index IDA `Vee`, optional `Vps`, ordering metadata, and JLD2 writing.
  - It explicitly warns that this is density-density/two-index IDA, not a full four-index Coulomb Hamiltonian export.

Proposed private Ham payload:
- Helper/object name:
  - `_pqs_source_box_route_driver_complete_core_shell_ham_payload(...)`
  - payload slot: `complete_core_shell_ham_payload`
- Preferred placement:
  - Build it immediately after the existing complete core/shell diagnostic route payload has materialized or blocked `final_basis`, `h1_payload`, `density_inputs`, and `h1_j_payload`.
  - Add it as a compact slot on `_PQSCompleteCoreShellDiagnosticRoutePayload`, or build a sibling private assembly payload from that object. I prefer the slot because it prevents a new adjacent field cloud.
- Suggested fields:
  - `object_kind`
  - `status`
  - `blocker`
  - `route_family`
  - `source_payload`
  - `source_plan_summary`
  - `final_basis`
  - `final_basis_summary`
  - `h1_payload`
  - `one_body_hamiltonian`
  - `density_inputs`
  - `density_input_summary`
  - `h1_j_payload`
  - `density_interaction`
  - `electron_electron_representation`
  - `coulomb_expansion`
  - `coulomb_expansion_summary`
  - `center_summaries`
  - `dimension_summary`
  - `ordering_summary`
  - `convention_labels`
  - `missing_inputs`
  - `summary`
  - `metadata`
- The payload should carry heavy objects internally and expose compact summaries only. It should not flatten matrices or provenance into report fields.

Required convention labels:
- `route_family = :pqs_source_box`
- `fixture_role = :route_smoke` or equivalent private fixture label when known
- `final_basis_kind = :pqs_complete_core_shell_final_basis`
- `support_row_order = :core_then_shell` from final-basis/density-interaction objects
- `one_body_hamiltonian_kind = :pqs_complete_core_shell_final_one_electron_hamiltonian`
- `nuclear_charge_application_stage = :hamiltonian_assembly`
- `center_summation_stage = :hamiltonian_assembly`
- `electron_electron_representation = :pre_final_density_interaction`
- `interaction_model = :density_density_pre_final_gauge`
- `density_gauge = :pre_final_localized_positive_weight`
- `raw_pair_factor_convention = :raw_numerator`
- `weight_application_stage = :pre_final_density_interaction_boundary`
- `final_orbital_consumption_rule = :combined_lowdin_cleanup_times_final_coefficients`
- `signed_final_weight_division_used = false`
- `raw_no_division_used = false`
- `density_normalized_pair_terms_used_as_authority = false`
- `fixed_block_pair_data_authority_used = false`

Required nonclaims:
- No production HF solver.
- No RHF product surface.
- No public API.
- No export/artifact writing in the first payload pass.
- No IDA/MWG semantic promotion.
- No full four-index Coulomb tensor.
- No GTO behavior.
- No CR2 science acceptance.
- No serious-HF claim.
- No shell/support-row contraction promoted to route algorithm.
- No retained diagnostic weights promoted to final IDA/quadrature weights.

Existing objects to carry rather than flatten:
- Carry `_PQSCompleteCoreShellDiagnosticRoutePayload` or its constituent route-owned objects:
  - `source_payload`
  - `source_payload.source_plan`
  - `source_payload.coulomb_expansion`
  - `final_basis`
  - `h1_payload`
  - `h1_payload.final_hamiltonian`
  - `density_inputs`
  - `h1_j_payload`
  - `h1_j_payload.density_interaction`
- Do not copy report aliases such as `complete_core_shell_h1_j_h1_energy`, `complete_core_shell_h1_j_self_coulomb`, or `complete_core_shell_h1_j_density_gauge` into a new scalar cloud.
- Do not derive the Ham payload from `cartesian_report(...)`; it is summary/report-facing and loses the object authority downstream consumers need.

Missing data for downstream handoff:
- Already route-owned:
  - final basis object and summaries;
  - dense final-basis one-body Hamiltonian matrix;
  - center/nuclear summaries for the one-body term;
  - Coulomb expansion provenance;
  - source plan / retained transform provenance;
  - support-row density weights and raw pair-factor provenance;
  - pre-final density-interaction matrix and final-to-pre-final mapping.
- Available only through diagnostic helpers:
  - the H1/J `density_interaction` object and self-Coulomb diagnostic live under the H1/J diagnostic path. That is acceptable for the first private Ham payload if the payload labels the representation as diagnostic/private and pre-final gauge.
- Available only as report aliases:
  - summary scalars such as H1 energy, self-Coulomb, density gauge, and driver-route materialization status. These are not sufficient for a downstream Ham handoff.
- Truly missing:
  - a reviewed downstream-facing schema for consuming `pre_final_pair_matrix` plus `final_to_pre_final_coefficients`;
  - a full four-index final-basis Coulomb tensor, if a downstream consumer requires one;
  - final IDA/MWG weights or a promoted IDA density-density payload;
  - artifact/export writer for this PQS Ham payload;
  - CR2/HFDMRG acceptance validation against the new payload.

Smallest next implementation pass:
- Files:
  - `src/pqs_source_box_route_driver_helpers.jl`
  - one focused nested route-driver test, preferably extending the existing one-center H1/J route-smoke coverage if present; otherwise add a narrow test that calls the one-center PQS source-box driver path and inspects the private assembly payload.
- Implementation:
  1. Add `_pqs_source_box_route_driver_complete_core_shell_ham_payload(...)`.
  2. Consume the existing complete core/shell diagnostic route payload, not report fields.
  3. Return blocked payloads for non-PQS routes, missing final basis, missing H1 payload/final Hamiltonian, missing density inputs, missing H1/J payload, or missing/materialization-blocked density interaction.
  4. On success, carry `final_basis`, `h1_payload`, `h1_payload.final_hamiltonian`, `density_inputs`, `h1_j_payload`, `h1_j_payload.density_interaction`, Coulomb/source provenance, compact summaries, convention labels, and nonclaim flags.
  5. Keep `cartesian_report(...)`, artifact writing, exports, RHF/SCF, and public API unchanged.
- Focused validation target:
  - a route-smoke assertion that the private Ham payload is materialized on the one-center source-box route and that:
    - final dimension is 223 for the current compact fixture;
    - one-body Hamiltonian status is materialized;
    - density interaction status is materialized;
    - density gauge is `:pre_final_localized_positive_weight`;
    - no export/artifact/public/RHF product flags are promoted.
- This seam is clear. It should be a private Ham payload pass, not an export/report/public API pass.

Git status:
```text
## main...origin/main
```

Deletion/shrinkage forecast:
- deleted: none in this no-edit audit.
- simplified: once the Ham payload exists, downstream handoff discussions can refer to one compact route-owned object instead of H1/J report aliases plus RHF validator fields.
- quarantined: RHF remains a private validator; current report/artifact/export vocabulary remains outside the first PQS Ham payload.
- not deleted because: existing report/materialization/export fields are compatibility surfaces for `cartesian_ham_builder` and White-Lindsey/low-order paths, and this pass forbids implementation.
- exact remaining caller/blocker: no private PQS Ham payload caller exists yet; the next pass must add the payload and one focused route-smoke inspection before CR2/HFDMRG can consume it.

-- repo-doer@macmini

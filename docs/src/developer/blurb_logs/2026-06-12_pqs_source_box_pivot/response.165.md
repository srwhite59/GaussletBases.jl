Pass 165 response - CR2 inspection artifact schema audit

Status: complete
Role: repo-doer@macmini
Scope: no-edit implementation-readiness audit only; no source/test/doc/runtime edits; no tests run.

Commands run:
- `git status --short --branch`
- `sed -n '1,220p' AGENTS.md`
- `sed -n '1,220p' BlurbStyle.md`
- `sed -n '1,260p' docs/src/developer/pqs_source_box_operator_framework.md`
- `sed -n '1,240p' docs/src/developer/successor_handoff_2026-06-12_pqs_source_box_pivot.md`
- `sed -n '1,260p' /Users/srw/Dropbox/codexhome/work/cr2/reports/be2_wl_pqs_downstream_handoff_request_2026-06-13.md`
- `sed -n '3240,3335p' src/pqs_source_box_diatomic_complete_core_shell.jl`
- `sed -n '3435,3495p' src/pqs_source_box_diatomic_complete_core_shell.jl`
- `sed -n '160,235p' test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
- `sed -n '1,220p' src/cartesian_bundle_export.jl`
- `sed -n '650,830p' src/cartesian_bundle_export.jl`
- `sed -n '1,180p' src/fullida_dense_export.jl`
- `sed -n '1,220p' src/sliced_ham_export.jl`
- `rg -n "WhiteLindsey|white_lindsey|_WhiteLindseyLowOrderHamBundleAdapter|white_lindsey_low_order|ordinary_cartesian_qiu_white_operators|write_cartesian_basis_bundle_jld2|cartesian_basis_bundle_payload" src test/nested`
- `rg -n "Be2|be2|white_lindsey|qiu_white|ordinary_cartesian_qiu_white" test src | head -120`
- `sed -n '700,830p' src/cartesian_bundle_export.jl`
- `sed -n '740,820p' test/nested/pqs_source_box_route_driver_report_runtests.jl`
- `sed -n '840,925p' test/nested/pqs_source_box_route_driver_report_runtests.jl`
- `rg -n "struct CartesianBasisRepresentation3D|struct OrdinaryCartesianOperators3D|ordinary_cartesian_qiu_white_operators\(" src`
- `sed -n '65,145p' src/cartesian_basis_representation.jl`
- `sed -n '90,135p' src/ordinary_qw_types_and_bases.jl`
- `sed -n '3080,3160p' src/pqs_source_box_route_driver_helpers.jl`

Proposed compact JLD2 schema:
- `schema/name = "be2_wl_pqs_handoff_inspection_bundle"`
- `schema/version = 1`
- `bundle/purpose = "cr2_read_only_hamiltonian_inspection"`
- `bundle/routes = ["white_lindsey", "pqs_source_box"]`
- `producer/package`, `producer/repo_commit`, `producer/dirty`, `producer/generated_at`
- per route under `routes/<route_label>/...`:
  - `route/label`, `route/family`, `route/kind`, `route/status`, `route/blocker`, `route/nonclaims/*`
  - `system/atom_symbols`, `system/nuclear_charges`, `system/cartesian_coordinates`, `system/nuclear_repulsion`, `system/electron_count`, `system/spin_sector`, `system/bond_axis`, `system/bond_length`, `system/units`
  - `fixture/q`, `fixture/n_s`, `fixture/spacing/*`, `fixture/parent_axis_counts`, `fixture/parent_dimension`, `fixture/backend`, `fixture/provenance/*`
  - `final_basis/final_dimension`, `final_basis/orbital_labels`, `final_basis/basis_centers`, `final_basis/integral_weights`, `final_basis/order_label`, `final_basis/support_indices_status`
  - `pre_final/dimension`, `pre_final/support_weight_count`, `pre_final/support_order`, `pre_final/pre_final_weights`, `pre_final/support_weights`, `pre_final/retained_diagnostic_weights_are_ida_weights = false`
  - `one_body/overlap`, `one_body/one_body_hamiltonian`, `one_body/kinetic`, `one_body/nuclear_by_center/<i>`, plus availability keys for each optional matrix family
  - `two_body/representation_kind`, `two_body/density_gauge`, `two_body/raw_pair_factor_convention`, `two_body/pre_final_pair_matrix`, `two_body/final_to_pre_final_coefficients`, `two_body/pre_final_weights`, `two_body/support_weights`, `two_body/availability/*`
  - `ordering/final_basis_order`, `ordering/pre_final_order`, `ordering/support_order`, `ordering/matrix_storage = "julia_column_major"`
  - `validation/finite/*`, `validation/symmetry_defect/*`, `validation/dimension_checks/*`, `validation/readiness/*`

Schema rule: write only plain arrays, strings, numbers, booleans, and small tuple/NamedTuple-like metadata converted to JLD2-friendly values. Do not write private `_PQS...` structs, route reports, or opaque Julia objects.

Small fingerprint format:
- Prefer TSV first, not JSON, unless manager explicitly approves adding/reusing a JSON dependency. The current export surfaces are JLD2-centric and already have simple scalar/string JLD2 metadata patterns; I did not find an established project JSON writer in the inspected export examples.
- Suggested TSV columns:
  - `route_label`
  - `route_family`
  - `route_kind`
  - `status`
  - `blocker`
  - `repo_commit`
  - `dirty`
  - `final_dimension`
  - `pre_final_dimension`
  - `support_weight_count`
  - `one_body_shape`
  - `two_body_shape`
  - `h1_lowest`
  - `h1_symmetry_defect`
  - `one_body_finite`
  - `two_body_finite`
  - `density_gauge`
  - `raw_pair_factor_convention`
  - `cr2_read_only_inspector_ready`
  - `cr2_solver_ready`
  - `cr2_export_ready`
  - `cr2_handoff_blocker`
  - `hfdmrg_density_density_ready`
  - `hamv6_export_ready`
  - `nuclear_charges`
  - `nuclear_repulsion`
  - `electron_count`
  - `spin_sector`
  - `parent_axis_counts`
  - `missing_objects`

PQS source-box availability map:
- Route labels/status/readiness: available from the diatomic complete core/shell route payload and consumer contract readiness. Current route-level blocker remains `:missing_hfdmrg_density_density_contract`; CR2 read-only view reports `cr2_read_only_inspector_ready = true` and `cr2_handoff_blocker = :missing_cr2_solver_handoff_format`.
- System metadata: available in the Hamiltonian handoff summary for the focused Be2 path: nuclear charges `(4.0, 4.0)`, coordinates `((-2.0,0.0,0.0),(2.0,0.0,0.0))`, nuclear repulsion `4.0`, electron count `8`, spin `:closed_shell_singlet`. Atom symbols, bond axis, and fixture axis metadata are available upstream through parent/route construction rather than all being copied into the handoff.
- Fixture metadata: q/n_s/spacing/parent-axis information is available at the assembly/fixture level. It should be lifted directly by the artifact writer, not copied into every downstream readiness object.
- Final-basis metadata: final dimension and final labels/centers/weights are available through the route final basis and handoff surfaces. The focused Be2 Hamiltonian fingerprint currently has final dimension `221` for the PQS diagnostic path.
- Pre-final/support metadata: available for the PQS CR2 read-only view through the density-interaction handoff: `density_gauge = :pre_final_localized_positive_weight`, `raw_pair_factor_convention = :raw_numerator`, support weight count `275`, and pre-final pair matrix shape `(221, 221)`.
- One-body arrays: final H1 is available as `one_body_hamiltonian` in the handoff. H1 low spectrum, finite/symmetry checks, and kinetic/electron-nuclear components are available from the route/H1 diagnostic payloads or can be computed by the writer.
- Two-body inspection data: PQS can provide a read-only pre-final density-interaction representation. This is inspection/debug data, not a solver-ready HamV6 or final dense V export.
- Ordering contracts: final/pre-final/support order labels should be explicit in the artifact. Current objects expose enough order and dimension information, but the writer should not rely on implicit matrix shape alone.

White-Lindsey availability map:
- Existing route-configured one-center ham export already writes a JLD2 `ham` group with `model_kind = "white_lindsey_low_order"`, `route_family = "white_lindsey_low_order"`, final-basis `overlap`, `one_body_hamiltonian`, `interaction_matrix`, orbital labels, basis centers, final integral weights, expansion exponents/coefficients, nuclear charge, and source/status strings.
- `CartesianBasisRepresentation3D` carries final metadata, coefficient map, parent labels/centers, support indices/states, and parent data. `OrdinaryCartesianOperators3D` carries overlap, H1, interaction matrix, raw-to-final, nuclear charges, kinetic one-body, and by-center nuclear terms when built through the ordinary QW path.
- Existing seed-based WL materialization still has a blocked ham preflight without the route-configured low-order builder: blocker `:missing_pure_low_order_fixed_block_density_density_interaction_builder`. The route-configured one-center path can write a ham bundle when expansion and Z are supplied.
- WL has final-basis density-density inspection data, not PQS pre-final/source-box support-row internals. `pre_final/*`, `support_weights`, `final_to_pre_final_coefficients`, and PQS raw pair convention should be `not_applicable` or `unavailable` for WL, not synthesized.
- WL one-body component availability is path-dependent: the specialized low-order ham adapter writes total overlap/H/interaction; ordinary Cartesian operators can carry kinetic and nuclear-by-center matrices, but the inspected WL low-order adapter does not currently emit those components.

Smallest implementation seam I would introduce first:
- Add one private inspection-bundle builder/writer boundary, not a public API: `_pqs_source_box_route_driver_be2_cr2_inspection_bundle_payload(...)` or shorter `_be2_cr2_ham_inspection_payload(...)` if kept near the route driver.
- It should return plain grouped values for both routes, then reuse the existing JLD2 value-writing pattern from `cartesian_bundle_export.jl` to write groups. Keep private structs out of the artifact. If implementation starts with one route, write the full schema with unavailable/not-applicable statuses for the other route rather than changing schema later.
- The writer should consume existing route-owned payloads directly: PQS Hamiltonian handoff + consumer readiness + route final-basis/H1 diagnostic payloads; WL route-configured ham adapter or ordinary Cartesian operator payload where available. It should not add new readiness scalar clouds to the route report.

What it would replace or shrink:
- It would make the current focused test’s repeated CR2 readiness scalar assertions less valuable. The test can eventually assert a compact inspection fingerprint tuple or TSV row instead of separately checking `ham_input_payload.summary`, handoff summary, consumer readiness flags, and CR2 view duplicates.
- It would let the route report stop accumulating copied CR2/consumer diagnostic fields. The CR2 read-only view can remain the route-owned status gate; the artifact payload becomes the detailed inspection carrier.
- It would avoid promoting HamV6, sliced Vblocks, final dense V, or solver-specific CR2 objects just to answer a read-only comparison request.

Exact first implementation pass recommended:
1. Add a private in-memory schema builder for the PQS route only, plus the fixed top-level schema metadata and unavailable WL group placeholders.
2. Write a compact JLD2 artifact and TSV fingerprint from that plain payload using existing JLD2 writer patterns.
3. Validate only with the focused Be2 ham payload fingerprint route test or an equally narrow artifact smoke that checks group names, shapes, finite/symmetry scalars, and readiness/nonclaim values.
4. Do not add H1/J density inputs, HamV6, dense final V, solver export, cache/checkpointing, or public exports.

Deletion/shrinkage forecast:
- deleted: none in this no-edit audit.
- simplified: future pass can simplify duplicate CR2 readiness checks in `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl` once the TSV/JLD2 fingerprint is the active inspection contract.
- quarantined: WL final-basis density-density matrix should be labeled as a final-basis inspection representation; PQS support-row/pre-final contraction should remain debug/oracle/inspection, not route-production behavior.
- not deleted because: current pass is no-edit; existing scalar assertions still protect the only implemented CR2 read-only inspection readiness seam.
- exact remaining caller/blocker: the current CR2 artifact writer does not exist; WL route-configured ham export is separate from the requested Be2 WL-vs-PQS inspection bundle; PQS downstream readiness still blocks solver/export claims on `:missing_hfdmrg_density_density_contract` and CR2 solver handoff on `:missing_cr2_solver_handoff_format`.

Git status:
`## main...origin/main`

-- repo-doer@macmini

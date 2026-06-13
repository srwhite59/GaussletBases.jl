Pass 132 response - no-edit Be2 WL/PQS driver comparison readiness audit

Scope:
- No file edits.
- No Julia commands run; read-only inspection was sufficient.
- No broad report test, hfdmrg, CR2 execution, artifact/export writing, RHF/SCF work, or fixture promotion.

Plain-language driver status:
- The driver can now produce a private, object-carrying PQS Ham payload for the compact one-center complete core/shell route through `cartesian_assembly(...)`.
- The driver can produce Be2/diatomic route reports and route-shape inventories for PQS, but the source-box materialization branch still reports `:pending_source_box_retained_route`.
- The driver has White-Lindsey/low-order Be2 materializer and Ham-adapter machinery, including object-carrying adapter state under materialization/probe paths, but it does not yet expose a matching compact WL Ham payload boundary parallel to the new PQS one.

One-center PQS status:
- Current usable path:
  - staged driver functions through `cartesian_assembly(...)`;
  - `assembly.complete_core_shell_diagnostic_route_payload.complete_core_shell_ham_payload`.
- The payload carries final basis, H1 payload/final Hamiltonian, density inputs, H1/J diagnostic payload, nested density interaction, Coulomb/source provenance, summaries, and convention labels.
- Validation exists in `test/nested/pqs_source_box_route_driver_complete_core_shell_ham_payload_runtests.jl`.
- This is still private route-smoke, not public API, not export/artifact, not IDA/MWG promotion, and not a physics endpoint.

Be2/diatomic PQS status:
- The Be2 PQS dry-run/report path exists:
  - `_pqs_source_box_route_driver_dry_run(...)`;
  - `cartesian_system -> cartesian_recipe -> cartesian_parent -> cartesian_shells -> cartesian_units -> cartesian_transforms -> cartesian_pair_terms -> cartesian_assembly`;
  - `test/nested/pqs_source_box_route_driver_report_runtests.jl` checks Be2 PQS route report metadata and standard unit inventory.
- The Be2 PQS report has source-box route skeleton/inventory data:
  - retained dimension `221` in the route report test;
  - unit keys `(:pqs_left, :pqs_right, :product)`;
  - pair families `(:pqs_pqs, :pqs_product, :product_product)`;
  - output representation `:retained_two_index_density_density`.
- The Be2 PQS materialization branch is not ready:
  - `shellization_source = :pending_source_box_route_shellization`;
  - `materialized_shellization_stage = :pending_source_box_retained_route`;
  - `final_integral_weights_status = :pending_final_ida_weights`;
  - `one_body_operator_status = :pending_source_box_retained_blocks`;
  - `ham_preflight_status = :not_applicable_to_pqs_source_box_route`;
  - `ham_operator_payload_status = :pending_source_box_retained_operator_payload`;
  - `ham_interaction_status = :pending_source_box_retained_density_density_blocks`;
  - `ham_bundle_export_status = :pending_source_box_retained_route`.
- The new complete-core/shell PQS Ham payload is one-center-only in practice today. It relies on the complete core/shell source-plan/final-basis/H1/H1J path; there is no route-owned Be2/diatomic PQS final-basis/H1/Ham payload seam yet.

Be2 WL status:
- WL/low-order Be2 has materially more route materialization machinery than Be2 PQS:
  - `_pqs_source_box_route_driver_diatomic_materializer_probe(...)`;
  - `_pqs_source_box_route_driver_diatomic_atom_growth_materializer_probe(...)`;
  - `_pqs_source_box_route_driver_route_configured_diatomic_basis_adapter(...)`;
  - `_pqs_source_box_route_driver_route_configured_diatomic_ham_adapter(...)`;
  - `_pqs_source_box_route_driver_white_lindsey_ham_preflight(...)`.
- The atom-growth materializer probe can carry `basis_adapter`, `ham_adapter`, `ham_adapter_summary`, retained dimension, support count, coverage audit, final-integral-weight status, and Ham adapter status.
- The WL diatomic Ham adapter carries an `operators` object when available, plus retained dimension, matrix sizes, nuclear metadata, one-body/operator status, interaction status, interaction treatment, and final-integral-weight status.
- However, this is not yet a clean WL counterpart to the PQS Ham payload. Some paths expose only preflight/materialization statuses; artifact paths can write bundles, but the first CR2 comparison should not depend on artifact/export plumbing.

Comparable payload target:
- Smallest useful first comparison: H1-only fingerprint pair.
  - Compare final-basis dimensions/order labels, center metadata, one-body Hamiltonian status, one-body matrix shape/finiteness/symmetry, nuclear charge convention, and route/source provenance.
  - This is the least ambiguous first Be2 comparison because PQS Be2 density-interaction/Ham payload is not route-owned yet.
- Next comparison after H1: H1 plus density-interaction surrogate.
  - PQS one-center currently labels `:pre_final_density_interaction`; Be2 PQS does not yet have the analogous diatomic object.
  - WL has density-density interaction matrices in the Ham adapter/export-oriented path, but the convention needs a compact payload label rather than relying on `ham_interaction_status`.
- Fuller Ham payload/export should wait until both sides have private object-carrying payloads with explicit conventions.

Conventions that must match or be explicitly labeled:
- final basis ordering and dimensions;
- center metadata, nuclear charges, and nuclear one-body storage;
- one-body terms and H1 matrix ordering;
- Coulomb expansion identity and coefficient source;
- density/electron-electron representation;
- density gauge;
- pair-factor convention;
- final integral/IDA/MWG status;
- export/artifact status;
- route authority labels: source-box-first PQS versus WL/low-order adapter;
- nonclaim flags: public API false, artifacts false unless explicitly requested, RHF product false, serious-HF false.

Exact blockers for CR2 Be2 WL/PQS comparison:
- Missing PQS Be2 object:
  - no route-owned Be2/diatomic PQS final-basis/H1/Ham payload equivalent to `complete_core_shell_ham_payload`;
  - Be2 PQS materialization status remains `:pending_source_box_retained_route`;
  - Be2 PQS one-body operator status remains `:pending_source_box_retained_blocks`.
- Missing WL counterpart object:
  - WL has `route_configured_diatomic_ham_adapter` and summaries, but no compact WL Ham payload boundary with the same fields/convention vocabulary as PQS.
- Missing comparison handoff:
  - no private WL/PQS Be2 fingerprint pair that CR2 can consume without parsing broad report fields or artifact metadata.
- Missing reviewed interaction comparison level:
  - H1-only is feasible to specify first; H1 plus density interaction needs explicit mapping between PQS pre-final density gauge and WL density-density matrix conventions.

Recommended next passes:
1. PQS Be2/diatomic Ham readiness fingerprint pass.
   - No broad implementation first.
   - Build the Be2 staged PQS route to `cartesian_assembly(...)` in a focused probe/test and inspect the private payload/blocker surface.
   - Record exact missing inputs for Be2 PQS final-basis/H1/Ham, rather than assuming the one-center complete-core/shell helper generalizes.
2. WL object-carrying Ham payload boundary pass.
   - Add a private compact WL payload wrapper around `route_configured_diatomic_ham_adapter`/`operators`.
   - Carry object references and compact summaries; do not write artifacts or add public API.
   - Use labels matching the PQS payload where possible: final ordering, one-body Hamiltonian, density-density interaction, Coulomb expansion, center metadata, final weights, export/public nonclaims.
3. Private Be2 WL/PQS handoff fingerprint pass.
   - Once the PQS missing-object status and WL payload boundary are clear, add a private comparison/handoff payload that emits H1-only fingerprints for both sides.
   - Keep CR2-facing shape private and summary/fingerprint-only at first.
   - Defer H1+density comparison until the PQS Be2 density-interaction object exists and its gauge can be compared honestly with WL.

Git status:
```text
## main...origin/main
```

Deletion/shrinkage forecast:
- deleted: none in this no-edit audit.
- simplified: the recommended path avoids parsing broad report/materialization field clouds by creating compact private payload/fingerprint boundaries.
- quarantined: RHF, hfdmrg, CR2 execution, public API, artifacts/exports, IDA/MWG promotion, and physics endpoint claims remain out of scope.
- not deleted because: existing report/materialization/export fields are compatibility surfaces for `cartesian_ham_builder` and current WL low-order tests.
- exact remaining caller/blocker: CR2 has no comparable Be2 WL/PQS private handoff yet; the immediate blocker is the missing Be2/diatomic PQS final-basis/H1/Ham payload seam, followed by the missing compact WL Ham payload wrapper.

-- repo-doer@macmini

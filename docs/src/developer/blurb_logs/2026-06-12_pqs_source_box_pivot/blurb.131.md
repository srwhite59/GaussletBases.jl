Pass 131 - add private complete-core/shell Ham payload seam

Purpose:

Implement the first private Hamiltonian-constructor payload boundary for the
one-center PQS source-box complete core/shell route.

Why now:

Pass 130 showed the route already owns the needed objects: source payload,
final basis, H1 payload/final Hamiltonian, density inputs, H1/J diagnostic
payload, and density interaction. RHF has served as a validator; do not expand
RHF. The next product boundary is a clean Ham payload that downstream consumers
can later inspect or hand off.

Governing boundary:

- Source-box-first PQS remains the algorithmic framing.
- Shell/support-row data is diagnostic/oracle support, not route authority.
- Retained diagnostic weights are not IDA/quadrature weights.
- Serious HF belongs to `codexhome/work/hfdmrg`.
- Downstream Cr2 validation belongs to the CR2 agent after this line is ready.

Files/surfaces:

- `src/pqs_source_box_route_driver_helpers.jl`
  - `_PQSCompleteCoreShellDiagnosticRoutePayload`
  - `_pqs_source_box_route_driver_complete_core_shell_diagnostic_route_payload(...)`
  - `_pqs_source_box_route_driver_complete_core_shell_h1_j_diagnostic_payload(...)`
  - `cartesian_assembly(...)`
- Add one focused test file if cleaner, for example:
  - `test/nested/pqs_source_box_route_driver_complete_core_shell_ham_payload_runtests.jl`
  Do not add it to the default nested runner if it is a cold route-smoke test.

Exact implementation:

1. Add a private helper:

   `_pqs_source_box_route_driver_complete_core_shell_ham_payload(...)`

   It should consume the complete core/shell route-owned objects, not
   `cartesian_report(...)` aliases.

2. Add a compact slot to `_PQSCompleteCoreShellDiagnosticRoutePayload`, named
   `complete_core_shell_ham_payload` or `ham_payload`.

   Prefer `complete_core_shell_ham_payload` for clarity unless the local style
   strongly argues otherwise.

3. On success, the payload should carry object references, not copy matrices
   into a report field cloud:

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
   - compact `dimension_summary`
   - compact `ordering_summary`
   - compact `convention_labels`
   - `missing_inputs`
   - `summary`
   - `metadata`

   If this exact list is too broad in code, keep the object-carrying fields and
   a compact summary/convention label set; do not create scalar report aliases.

4. Required success labels/nonclaims:

   - `electron_electron_representation === :pre_final_density_interaction`
   - `density_gauge === :pre_final_localized_positive_weight`
   - `raw_pair_factor_convention === :raw_numerator`
   - `support_row_order === :core_then_shell`
   - `signed_final_weight_division_used === false`
   - `raw_no_division_used === false`
   - `density_normalized_pair_terms_used_as_authority === false`
   - `public_api === false`
   - `exports_materialized === false`
   - `artifacts_materialized === false`
   - `rhf_product_surface === false` or equivalent nonclaim
   - `serious_hf_claim === false` or equivalent nonclaim

5. Required blockers should cover at least:

   - non-PQS route;
   - missing final basis;
   - missing H1 payload or final Hamiltonian;
   - missing density inputs;
   - missing H1/J payload;
   - missing density interaction.

Test shape:

- Build the current one-center PQS route through the normal stage functions,
  ending at `cartesian_assembly(...)`, rather than inspecting report aliases.
- Use the same compact route-smoke family as the recent probes:
  one Be-like center at the origin, `nuclear_charges = (4,)`,
  parent counts `(x = 7, y = 7, z = 7)`, `q = 5`, `n_s = 5`,
  `route_family = :pqs_source_box`, route shape
  `(:pqs_left, :product, :pqs_right)`, and the standard source-box retained
  rules.
- Inspect:
  `assembly.complete_core_shell_diagnostic_route_payload.complete_core_shell_ham_payload`
  or the final chosen slot name.
- Assert:
  - status is materialized/available for the Ham payload;
  - final dimension is 223;
  - one-body Hamiltonian object is available and finite/symmetric through the
    existing H1 payload object;
  - density interaction is available;
  - density gauge is `:pre_final_localized_positive_weight`;
  - electron-electron representation is `:pre_final_density_interaction`;
  - export/artifact/public/RHF-product/serious-HF flags remain false.

Trust boundary:

- Do not add report fields or public API.
- Do not add export/artifact writing.
- Do not route-wire RHF or add SCF/DIIS work.
- Do not execute `hfdmrg` or CR2.
- Do not promote any fixture to physics endpoint.
- Do not add full four-index Coulomb tensors.
- Do not change IDA/MWG semantics.

Validation:

- Run the focused new/updated test directly. This route-smoke may take over
  60 seconds cold; that is acceptable for this pass because it validates the
  private Ham payload seam. Report elapsed time if available.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Decision rules:

- If the helper can be added with only the route helper file and one focused
  test, implement, validate, and commit.
- If the implementation requires report alias fields, public API, export
  plumbing, RHF route adoption, or scalar field clouds, stop and report the
  exact pressure point.
- If density interaction is unavailable from the H1/J payload, stop and report
  the exact missing object rather than synthesizing aliases.
- Do not request interactive escalation/approval; if approval is truly needed,
  write `ATTENTION.md` and stop.

Report back:

- Commit SHA if committed.
- Files changed.
- Validation commands/results, with elapsed time for any >60s route-smoke.
- Final payload/slot name and key status labels.
- Git status.
- Deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini

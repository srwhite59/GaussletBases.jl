Pass 257 response - independent H2 PQS supplement preflight input

Files changed:
- `test/driver_inputs/h2_pqs_q5_independent_source_box_r4_supplement_preflight.jl`
- `src/pqs_source_box_route_driver_helpers.jl`
- `docs/src/developer/cartesian_driver_endpoint_manifest.md`
- `test/nested/cartesian_report_stage_low_order_policy_runtests.jl`

New input:
- Path: `test/driver_inputs/h2_pqs_q5_independent_source_box_r4_supplement_preflight.jl`
- Includes: `h2_pqs_q5_independent_source_box_r4_h1_j.jl`
- Sets:
  - `artifact_role = :independent_h2_pqs_supplement_preflight_diagnostic`
  - `supplement_policy = :mwg_residual_gto`
  - `comparison_ready = false`
  - `comparison_blocker = :independent_pqs_supplement_preflight_only`
  - `physics_endpoint_ready = false`
  - `physics_endpoint_blocker = :missing_provider_gto_supplement_blocks`
  - `(run_final_basis, run_h1, run_h1_j, run_private_rhf) == (true, true, true, false)`

Classifier/manifest updates:
- Added `:independent_h2_pqs_supplement_preflight_diagnostic` to `_pqs_source_box_route_driver_independent_h2_pqs_artifact_role`.
- This keeps the independent route artifact behavior for the new role:
  - `source_backed_fixed_source_oracle_used = false`
  - `retained_transform_authority = :pqs_source_box_construction`
  - fake-PQS guard behavior unchanged.
- Added a manifest row for the independent H2 PQS q5 MWG/GTO supplement preflight, explicitly stating no provider blocks or supplemented values.

Blocked surfaces:
- Provider blocks were not implemented.
- Mixed gausslet/GTO matrices were not built.
- GTO/GTO matrices were not built.
- Residual MWG representation was not built.
- Combined density-density readiness was not built.
- Supplemented H1/H1-J/RHF values were not added.
- No supplemented WL/QW scalar references were added.
- No public/export/CR2/HamV6 readiness was added.
- Fake-PQS evidence was not used as independent-PQS evidence.

Cleanup offset:
- Deleted 24 stale report-stage `low_order_route_core_*` alias mirror assertion lines from `test/nested/cartesian_report_stage_low_order_policy_runtests.jl`.
- The active compact summary assertions immediately above remain, including route-core pair count/order/readiness and typed pair-operator blocker checks.

Validation:
- `julia --project=. -e 'include("test/driver_inputs/h2_pqs_q5_independent_source_box_r4_supplement_preflight.jl"); @assert route_family === :pqs_source_box; @assert route_kind === :bond_aligned_diatomic_independent_pqs_source_box_core_shell; @assert artifact_role === :independent_h2_pqs_supplement_preflight_diagnostic; @assert supplement_policy === :mwg_residual_gto; @assert comparison_ready == false; @assert comparison_blocker === :independent_pqs_supplement_preflight_only; @assert physics_endpoint_blocker === :missing_provider_gto_supplement_blocks; @assert (run_final_basis, run_h1, run_h1_j, run_private_rhf) == (true, true, true, false); println("supplement_preflight_input_smoke ok")'`
  - passed: `supplement_preflight_input_smoke ok`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed: `load ok`
- `julia --project=. -e 'Meta.parseall(read("test/nested/cartesian_report_stage_low_order_policy_runtests.jl", String)); println("report_stage_test_parse ok")'`
  - passed: `report_stage_test_parse ok`
- `julia --project=. -e 'using GaussletBases; @assert GaussletBases._pqs_source_box_route_driver_independent_h2_pqs_artifact_role(:independent_h2_pqs_supplement_preflight_diagnostic); @assert !GaussletBases._pqs_source_box_route_driver_independent_h2_pqs_artifact_role(:fake_pqs_source_backed_wl_reproduction); println("independent_h2_role_classifier_smoke ok")'`
  - passed: `independent_h2_role_classifier_smoke ok`
- `git diff --check`
  - passed

Scoped line count:
- `src + test + bin`: `+13/-24`, net `-11`.
- Breakdown: `src` `+1`, new driver input `+12`, stale report-stage test assertions `-24`.
- Manifest docs: `+1` outside scoped line budget.

Git status:
```text
## main...origin/main
 M docs/src/developer/cartesian_driver_endpoint_manifest.md
 M src/pqs_source_box_route_driver_helpers.jl
 M test/nested/cartesian_report_stage_low_order_policy_runtests.jl
?? test/driver_inputs/h2_pqs_q5_independent_source_box_r4_supplement_preflight.jl
```

Deletion/shrinkage result:
- deleted: 24 stale report-stage RouteCore alias mirror assertion lines.
- simplified: independent supplement preflight now has its own tiny include/override input and manifest row instead of overloading H1-J or private-RHF inputs.
- quarantined: provider-block and supplemented-value work remains blocked behind `:missing_provider_gto_supplement_blocks`; fake-PQS remains separate.
- not deleted because: the compact summary assertions still protect the active report-stage low-order contract, and supplement helpers remain live staging surfaces.
- exact remaining caller/blocker: `cartesian_assembly` still builds the supplement request/representation/preflight payloads; the new preflight role reaches independent-PQS artifact classification, but actual supplement materialization remains blocked by `:missing_provider_gto_supplement_blocks`.

-- repo-doer@macmini

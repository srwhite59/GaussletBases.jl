Pass 252 - independent H2 PQS private RHF diagnostic

Context:
- Current HEAD should include `2cb0a11e Clarify independent H2 PQS input taxonomy`.
- The independent H2 PQS route is now fake-free through H1-J:
  - support/source plan from PQS/source-box authority, not WL/QW fixed-source;
  - retained counts `(275, 98, 98)`;
  - final dimension `471`;
  - final basis materialized with small identity-overlap error;
  - H1 materialized with finite symmetric `(471, 471)` matrix;
  - H1-J/density diagnostic materialized with
    `density_gauge = :pre_final_localized_positive_weight`,
    `raw_pair_factor_convention = :raw_numerator`, and
    `h1_j_self_coulomb = 0.4569117646737236`.
- Pass 251 added stage-specific driver inputs:
  - `test/driver_inputs/h2_pqs_q5_independent_source_box_r4.jl`
  - `test/driver_inputs/h2_pqs_q5_independent_source_box_r4_final_basis.jl`
  - `test/driver_inputs/h2_pqs_q5_independent_source_box_r4_h1.jl`
  - `test/driver_inputs/h2_pqs_q5_independent_source_box_r4_h1_j.jl`
- The fake-PQS H2 463 route remains quarantined as a source-backed WL/QW
  reproduction only. Do not use it as evidence for independent PQS.

Known relevant code surfaces:
- `src/pqs_source_box_diatomic_complete_core_shell.jl`
  - `_pqs_source_box_route_driver_diatomic_physical_gausslet_rhf_input_contract`
  - `_pqs_source_box_route_driver_diatomic_physical_gausslet_rhf_execution_payload`
- `src/pqs_source_box_route_driver_helpers.jl`
  - private RHF report merge fields around the complete-core/shell report path
- `src/pqs_source_box_route_driver_reporting.jl`
  - writes the `private_rhf/*` artifact group and route RHF status fields
- Existing fake-PQS RHF test can be used only as a schema/example reference:
  - `test/nested/cartesian_ham_builder_h2_fake_pqs_wl_source_backed_r4_runtests.jl`

Task:
Run and, if narrowly necessary, repair the existing private RHF diagnostic seam
for the independent H2 PQS H1-J basis.

Preferred route:
1. Start from
   `test/driver_inputs/h2_pqs_q5_independent_source_box_r4_h1_j.jl`
   with driver overrides:
   - `run_private_rhf=true`
   - `private_rhf_electron_count=2`
   - a temporary `outfile`
   - `save_tsv=false`
2. If that works cleanly, add the smallest artifact/test coverage needed to lock
   the diagnostic contract.
3. If a tiny input variant is cleaner than command-line overrides, add:
   `test/driver_inputs/h2_pqs_q5_independent_source_box_r4_private_rhf.jl`
   as an include/override file. Offset any added lines in the same pass.

Required artifact facts if RHF materializes:
- `route/artifact_role` stays independent-H2-PQS diagnostic, not public endpoint.
- `route/fake_pqs_enabled == false`.
- `route/source_backed_fixed_source_oracle_used == false`.
- `route/retained_transform_authority == :pqs_source_box_construction`.
- `private_rhf/input_contract_status ==
  :available_pqs_physical_gausslet_rhf_input_contract`.
- `private_rhf/execution_status ==
  :materialized_pqs_physical_gausslet_private_rhf_execution`.
- `private_rhf/requested == true`.
- `private_rhf/executed == true`.
- `private_rhf/materialized == true`.
- `private_rhf/converged == true`.
- `private_rhf/electron_count == 2`.
- `private_rhf/occupation_policy == :closed_shell_rhf`.
- `private_rhf/total_energy`, `one_body_energy`, `two_body_energy`,
  `density_trace`, `idempotency_residual`, and `commutator_residual` are finite
  when present.
- Keep and report the existing final-density one-step consistency status if it
  is available.

If RHF does not materialize:
- Do not force it by changing the physics convention.
- Preserve the strongest available input-contract facts.
- Report the exact blocker from the artifact and code path.
- Leave the route diagnostic-only.

Strict exclusions:
- Do not implement a new solver.
- Do not add supplement provider blocks, CR2, export/HamV6, public API, or
  public solver readiness.
- Do not compare the independent PQS RHF value to supplemented WL/QW references.
- Do not use the fake-PQS source-backed route for independent-PQS evidence.
- Do not broaden report-key clouds. If a compact existing field is enough, use
  it.
- Do not run broad stale integration gates. In particular, do not use
  `test/nested/cartesian_pair_stage_low_order_policy_runtests.jl` as validation.

Line budget:
- Keep `src + test + bin` net-negative if source/test/bin files are touched.
- If this is mostly an artifact/test pass and needs a tiny positive source
  change, offset it with a mature deletion candidate from the old-flat audit.
- If the existing RHF seam works with no source edits, still try to take a small
  safe deletion offset, but do not delay the pass on a risky deletion.

Validation:
- Before any long route run, state why it is needed. This focused private RHF
  driver run is allowed to take more than 60 seconds because it is the current
  H2 physics diagnostic seam.
- Run one focused independent H2 PQS private RHF driver/artifact check.
- Run `julia --project=. -e 'using GaussletBases; println("load ok")'` if source
  code is touched.
- Run `git diff --check`.
- Do not run the full suite.

Report:
- RHF input-contract status and blocker.
- RHF execution status and blocker.
- Energies and convergence diagnostics if materialized:
  total, one-body, two-body, iteration count, density trace, idempotency
  residual, commutator residual, energy delta, final-density one-step status.
- Confirmation that fake-PQS/source-backed authority remains false for the
  independent route.
- Validation commands and elapsed time for the focused RHF run.
- Deletion/shrinkage result:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini

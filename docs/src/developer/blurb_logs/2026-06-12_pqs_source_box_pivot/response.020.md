Ran the real projected-q-shell shell-realization oracle probe for the new CPBM final-basis object.

Artifacts:

- `tmp/work/pqs_shell_final_basis_real_oracle_probe.jl`
- `tmp/work/pqs_shell_final_basis_real_oracle_probe_output.txt`

What the probe did:

- Built the existing cubic projected-q-shell fixture:
  - `current_box = (1:5, 1:5, 1:5)`
  - `inner_box = (2:4, 2:4, 2:4)`
  - `q = 5`, `L = 5`
- Obtained the old oracle shell plan from `CCPM._pqs_shell_realization_plan(...)`.
- Obtained shell support overlap from `CCPM._pqs_product_box_support_overlap_matrix(...)`.
- Built a CRPS raw source plan and PQS boundary retained rule with matching source dims/order.
- Called `CPBM.pqs_source_shell_realization_final_basis(...)` with:
  - `shell_support_indices`
  - `shell_overlap`
  - `shell_projection = shell_plan.shell_projection_matrix`
  - `lowdin_cleanup = shell_plan.lowdin_cleanup`
- Did not call `_pqs_current_route_safe_term_matrices(...)`.
- Did not materialize one-body operators, H1, IDA, density-density, RHF, drivers, exports, artifacts, full-parent fallback, direct Cartesian fallback, or old fixed-block matrix authority.

Key diagnostics:

- `status = :available_pqs_shell_realization_final_basis`
- `blocker = nothing`
- `descriptor_kind = :projected_q_shell`
- `source_mode_dims = (5, 5, 5)`
- `source_mode_count = 125`
- `retained_count = 98`
- `shell_support_count = 98`
- `final_retained_count = 98`
- `shell_projection_shape = (98, 98)`
- `lowdin_cleanup_shape = (98, 98)`
- `final_shell_coefficients_shape = (98, 98)`
- `projected_boundary_overlap_rank = 98`
- `final_overlap_rank = 98`
- `final_overlap_identity_error = 2.2360255762214361e-14`
- `shell_plan_isometry_error = 2.2360255762214361e-14`
- `coefficient_error = 0.0`
- `final_overlap_error = 0.0`
- `identity_error_delta = 0.0`
- `one_body_operator_blocker = :missing_pqs_shell_projected_one_body_operator_materialization`
- `one_body_operator_materialized = false`
- `h1_solve_materialized = false`
- `ida_data_materialized = false`
- `rhf_materialized = false`
- `driver_route_materialized = false`
- `artifacts_materialized = false`
- `old_shell_plan_oracle_used = true`
- `current_route_safe_term_matrices_called = false`
- `old_fixed_block_matrix_authority_used = false`

Result:

The new CPBM final-basis object matches the old shell-realization oracle exactly for:

- `final_shell_coefficients == shell_plan.shell_isometry_matrix`
- `final_overlap == shell_plan.realized_overlap`
- final identity error equal to `shell_plan.isometry_error`

The final-basis object is validated on a real projected-q-shell shell-realization fixture. The next exact blocker remains:

```text
:missing_pqs_shell_projected_one_body_operator_materialization
```

Validation:

- `julia --project=. tmp/work/pqs_shell_final_basis_real_oracle_probe.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Deletion/shrinkage report:

- No tracked source, docs, or tests changed in this pass; this was intentionally probe-only.
- Nothing was deleted yet because the final-basis object is validated, but no final shell-projected one-body operator path exists yet.
- The new ignored `tmp/work` probe earned its temporary carrying cost by validating the CPBM final-basis object against real old projected-q-shell shell-realization data before operator work.
- The identity-only pass-019 contract test remains useful as the small module contract; this probe is the real-fixture oracle check and should not become a permanent test unless a future regression needs it.
- Remaining stale surfaces to retire next: bridge/readiness tests that only preserve `:shell_realization_not_materialized`, old retained-source H1 probes, and current-route support-local safe-term oracle vocabulary should shrink after final shell-projected one-body materialization and a real final PQS H1 probe exist.

-- repo-doer@macmini

Pass 106 response - one-step RHF density contraction seam

Contraction decision:
- Unambiguous enough to implement.
- Existing `pqs_complete_core_shell_pre_final_density_interaction(...)` provides the pre-final positive-weight two-index interaction and `final_to_pre_final_coefficients`.
- Existing `_pqs_complete_core_shell_restricted_one_orbital_interaction_energy(...)` defines the direct-minus-exchange convention as:
  - symmetrized pre-final interaction;
  - orbital density `rho = c_pre * c_pre'`;
  - occupations from `diag(rho)`;
  - `2 * dot(n, V*n) - dot(vec(rho), vec(V .* rho))`.
- The new seam applies the same convention to an externally supplied spin-summed closed-shell final density by dividing by the explicit RHF occupancy before pre-final contraction.

Files changed:
- `src/pqs_multilayer_complete_core_shell_rhf.jl`
- `test/nested/pqs_multilayer_complete_core_shell_rhf_one_step_runtests.jl`
- `test/nested/runtests.jl`

Helper/object names:
- Added private helper:
  `_pqs_multilayer_complete_core_shell_rhf_one_step_payload(...)`
- Added private support helpers:
  `_pqs_multilayer_complete_core_shell_rhf_one_step_missing_inputs(...)`
  `_pqs_multilayer_complete_core_shell_rhf_one_step_missing_blocker(...)`
  `_pqs_multilayer_complete_core_shell_rhf_one_step_blocked_payload(...)`
  `_pqs_multilayer_complete_core_shell_rhf_density_interaction_payload(...)`
  `_pqs_multilayer_complete_core_shell_rhf_h1_payload(...)`

Payload shape:
- `object_kind`
- `status`
- `blocker`
- `missing_inputs`
- `final_density`
- `effective_fock_matrix`
- `fock_matrix`
- `one_body_energy`
- `two_body_energy`
- `total_energy`
- `summary`
- `metadata`

Status labels:
- available/materialized:
  `:materialized_pqs_multilayer_complete_core_shell_rhf_one_step_payload`
- blocked:
  `:blocked_pqs_multilayer_complete_core_shell_rhf_one_step_payload`

Blocker labels:
- missing/input blockers:
  `:missing_rhf_input_contract`
  `:missing_h1_payload`
  `:missing_density_interaction`
  `:missing_final_density`
  `:missing_rhf_one_step_inputs`
- validation blockers:
  `:nonsquare_final_density`
  `:final_density_dimension_mismatch`
  `:nonfinite_final_density`
  `:nonsymmetric_final_density`
  `:electron_trace_mismatch`
  `:h1_dimension_mismatch`
  `:nonfinite_h1_matrix`
  `:final_to_pre_final_dimension_mismatch`
  `:density_interaction_dimension_mismatch`
  `:nonfinite_density_interaction`

Notes:
- The helper is private and not exported.
- It does not add SCF, RHF convergence, inferred electron count, route-driver fields, exports, artifacts, IDA/MWG, GTO, or fixture promotion.
- The pass-105 input contract validates H1 availability but does not carry the large H1 matrix. Therefore the one-step helper accepts `h1_payload` explicitly, and can also extract it from a compact payload that carries `h1_payload`.
- Density interaction may be supplied directly or extracted from an H1/J-style payload carrying `density_interaction`.
- Summary and metadata carry the required nonclaims:
  `scf_materialized = false`
  `rhf_converged = false`
  `rhf_energy_materialized = false`
  `driver_route_materialized = false`
  `exports_materialized = false`
  `artifacts_materialized = false`

Validation commands/results:
- `julia --project=. test/nested/pqs_multilayer_complete_core_shell_rhf_one_step_runtests.jl`
  - passed: 31/31
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed: `load ok`
- `git diff --check`
  - passed

Git status:

```text
## main...origin/main
 M src/pqs_multilayer_complete_core_shell_rhf.jl
 M test/nested/runtests.jl
?? test/nested/pqs_multilayer_complete_core_shell_rhf_one_step_runtests.jl
```

Deletion/shrinkage report:
- deleted: none
- simplified: none yet; this pass creates the private one-step seam needed before future RHF/J density cleanup.
- quarantined: RHF remains private diagnostic-only through helper summary/metadata nonclaims.
- not deleted because: this pass was additive by assignment and did not retire existing H1/J or input-contract callers.
- exact remaining caller/blocker: no route-driver caller yet; next cleanup remains blocked on manager approval to wire or consume this private payload from a specific diagnostic path.

-- repo-doer@macmini

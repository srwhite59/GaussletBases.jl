Pass 112 response - SCF final one-step/final-density consistency

Files changed:
- `src/pqs_multilayer_complete_core_shell_rhf.jl`
- `test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl`

Change:
- Updated `_pqs_multilayer_complete_core_shell_rhf_scf_payload(...)` so the converged return recomputes `_pqs_multilayer_complete_core_shell_rhf_one_step_payload(...)` using the returned `final_density`.
- Returned `final_one_step_payload` is now the recomputed payload.
- `summary.final_total_energy` now comes from the recomputed final one-step payload.
- `iteration_records` remain unchanged as convergence history.
- The blocked/nonconverged summary now labels that its one-step payload is not recomputed or is the last attempted input-density payload.

SCF summary labels added/changed:
- `orbital_metric = :ordinary_orthonormal_final_basis`
- `final_one_step_recomputed = true` on converged returns
- `final_one_step_density_matches_final_density`
- `final_one_step_density_match_error`
- `converged_iteration_input_total_energy`
- blocked/nonconverged path labels:
  - `final_one_step_recomputed = false`
  - `final_one_step_payload_role = :not_recomputed_or_last_attempted_input_density`
  - `final_density_role = :not_materialized_or_last_post_diagonalization_density`

Test updates:
- Asserted `scf.final_one_step_payload.final_density ≈ scf.final_density`.
- Asserted `scf.summary.final_one_step_recomputed === true`.
- Asserted `scf.summary.orbital_metric === :ordinary_orthonormal_final_basis`.
- Asserted `scf.summary.final_total_energy ≈ scf.final_one_step_payload.total_energy`.
- Preserved route/report/export/artifact/public nonclaim assertions.

Validation commands/results:
- `julia --project=. test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl`
  - passed: 35/35
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed: `load ok`
- `git diff --check`
  - passed

Git status:

```text
## main...origin/main
 M src/pqs_multilayer_complete_core_shell_rhf.jl
 M test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl
```

Deletion/shrinkage report:
- deleted: none
- simplified: final SCF energy/payload labels now point to the returned final density instead of the converged-iteration input density.
- quarantined: SCF remains private diagnostic/prototype behavior; no route wiring, report aliases, public API, exports, or artifacts were added.
- not deleted because: this was a corrective consistency pass on an active private helper.
- exact remaining caller/blocker: no route-driver caller yet; real compact PQS probe is now blocked on manager approval for a local `tmp/work` probe or explicit request-object wiring.

-- repo-doer@macmini

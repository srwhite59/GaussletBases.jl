Pass 108 response - private RHF initial-density payload

Files changed:
- `src/pqs_multilayer_complete_core_shell_rhf.jl`
- `test/nested/pqs_multilayer_complete_core_shell_rhf_initial_density_runtests.jl`
- `test/nested/runtests.jl`

Helper/object names:
- Added private helper:
  `_pqs_multilayer_complete_core_shell_rhf_initial_density_payload(...)`
- Added private blocked-payload helper:
  `_pqs_multilayer_complete_core_shell_rhf_initial_density_blocked_payload(...)`
- Object kind:
  `:pqs_multilayer_complete_core_shell_rhf_initial_density_payload`

Behavior:
- Builds H1-Aufbau spin-summed final-basis density from an available RHF input contract and materialized H1 payload.
- Symmetrizes the H1 matrix before `eigen(Symmetric(...))`.
- Occupies the lowest `nocc` orbitals with occupancy 2.
- Returns occupied orbital coefficients, final density, eigenvalues, occupied eigenvalues, electron trace, compact summary, and metadata.
- Reuses the existing private H1-payload extraction helper so a compact H1/J/diagnostic payload carrying `h1_payload` can be used.
- Does not call `_pqs_multilayer_complete_core_shell_rhf_one_step_payload(...)`.
- Does not add SCF, RHF convergence, RHF energy, route-driver wiring, report aliases, exports, artifacts, GTO, IDA/MWG, or fixture promotion.

Status/blocker labels:
- available/materialized:
  `:materialized_pqs_multilayer_complete_core_shell_rhf_initial_density_payload`
- blocked:
  `:blocked_pqs_multilayer_complete_core_shell_rhf_initial_density_payload`
- blockers:
  `:missing_rhf_input_contract`
  `:missing_h1_payload`
  `:missing_rhf_initial_density_inputs`
  `:h1_dimension_mismatch`
  `:nonfinite_h1_matrix`
  `:insufficient_final_dimension_for_occupation`

Required nonclaims:
- `initial_density_source = :h1_aufbau`
- `scf_materialized = false`
- `rhf_converged = false`
- `rhf_energy_materialized = false`
- `driver_route_materialized = false`
- `exports_materialized = false`
- `artifacts_materialized = false`

Validation commands/results:
- `julia --project=. test/nested/pqs_multilayer_complete_core_shell_rhf_initial_density_runtests.jl`
  - passed: 27/27
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed: `load ok`
- `git diff --check`
  - passed

Git status:

```text
## main...origin/main
 M src/pqs_multilayer_complete_core_shell_rhf.jl
 M test/nested/runtests.jl
?? test/nested/pqs_multilayer_complete_core_shell_rhf_initial_density_runtests.jl
```

Deletion/shrinkage report:
- deleted: none
- simplified: none yet; this pass adds the initial-density boundary required before any future SCF loop.
- quarantined: initial-density payload remains private diagnostic-only with explicit nonclaims.
- not deleted because: no existing initial-density path was being replaced in this pass.
- exact remaining caller/blocker: no route-driver caller and no SCF loop; future use remains blocked on manager-approved private RHF iteration or diagnostic wiring.

-- repo-doer@macmini

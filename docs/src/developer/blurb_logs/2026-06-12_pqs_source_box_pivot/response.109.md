Pass 109 response - private RHF SCF payload

Files changed:
- `src/pqs_multilayer_complete_core_shell_rhf.jl`
- `test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl`
- `test/nested/runtests.jl`

Helper/object names:
- Added private helper:
  `_pqs_multilayer_complete_core_shell_rhf_scf_payload(...)`
- Added private support helpers:
  `_pqs_multilayer_complete_core_shell_rhf_scf_controls_blocker(...)`
  `_pqs_multilayer_complete_core_shell_rhf_scf_blocked_payload(...)`
- Object kind:
  `:pqs_multilayer_complete_core_shell_rhf_scf_payload`

SCF control defaults:
- `max_iterations = 25`
- `density_atol = 1.0e-8`
- `energy_atol = 1.0e-10`

Behavior:
- Starts from a spin-summed final density.
- Builds initial density internally with `_pqs_multilayer_complete_core_shell_rhf_initial_density_payload(...)` when one is not supplied.
- On each iteration, calls `_pqs_multilayer_complete_core_shell_rhf_one_step_payload(...)`, symmetrizes the effective Fock matrix, diagonalizes it in the ordinary final-basis metric, occupies the lowest `nocc` orbitals with occupancy 2, and builds the next spin-summed density.
- Iteration records are compact: iteration number, total energy, density change, energy change, density convergence, energy convergence, and combined convergence.
- Iteration 1 has no previous energy, so `energy_change = nothing` and energy convergence is treated as satisfied for that first comparison. Density convergence is still required.

Status/blocker labels:
- converged/materialized:
  `:materialized_pqs_multilayer_complete_core_shell_rhf_scf_payload`
- blocked:
  `:blocked_pqs_multilayer_complete_core_shell_rhf_scf_payload`
- blockers:
  `:missing_rhf_input_contract`
  `:missing_h1_payload`
  `:missing_density_interaction`
  `:missing_initial_density`
  `:missing_rhf_scf_inputs`
  `:invalid_scf_controls`
  `:scf_not_converged`
  plus propagated blockers from the private initial-density or one-step payloads.

Tiny-fixture convergence result:
- Fixture: 2x2 diagonal H1, zero density interaction, two electrons.
- Converged in iteration 1.
- Final density trace: 2.0.
- Final total energy: -2.0.
- Final density change: 0.0.
- First iteration energy change: `nothing` by contract.
- Nonclaims remain false for driver route, route report, exports, artifacts, and public API.

Validation commands/results:
- `julia --project=. test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl`
  - passed: 30/30
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed: `load ok`
- `git diff --check`
  - passed

Git status:

```text
## main...origin/main
 M src/pqs_multilayer_complete_core_shell_rhf.jl
 M test/nested/runtests.jl
?? test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl
```

Deletion/shrinkage report:
- deleted: none
- simplified: none yet; this pass intentionally adds the first private SCF diagnostic loop.
- quarantined: SCF remains private diagnostic/prototype behavior with `public_api = false`, no route report, no driver route, no exports, and no artifacts.
- not deleted because: no previous SCF payload existed to replace.
- exact remaining caller/blocker: no route-driver caller; production/report adoption remains blocked on manager-approved route wiring and physics validation.

-- repo-doer@macmini

Pass 107 response - RHF one-step convention cross-check

Scope:
- Test-only.
- No production code changes.
- No SCF, route wiring, report aliases, exports, artifacts, GTO, IDA/MWG, or fixture promotion.

Change:
- Updated `test/nested/pqs_multilayer_complete_core_shell_rhf_one_step_runtests.jl`.
- Extended the synthetic density interaction with the metadata fields consumed by the existing one-orbital diagnostic:
  `pre_final_weight_count` and `final_to_pre_final_reconstruction_error`.
- Added a one-orbital closed-shell density cross-check:
  - occupied final orbital `c = [1.0, 0.0]`;
  - spin-summed density `2.0 .* (c * c')`;
  - one-step helper `_pqs_multilayer_complete_core_shell_rhf_one_step_payload(...)`;
  - existing diagnostic `CartesianFinalBasisRealization.pqs_complete_core_shell_pre_final_orbital_self_coulomb(...)`.

Observed comparison:
- `one_step.two_body_energy = 2.0`
- `self_coulomb.self_coulomb = 2.0`
- They agree under `≈`.
- Density convention remains `:spin_summed_closed_shell_final_density`.
- Contraction rule remains `:pre_final_restricted_direct_minus_exchange_from_orbital_density`.
- `scf_materialized === false`.
- `rhf_converged === false`.

Validation commands/results:
- `julia --project=. test/nested/pqs_multilayer_complete_core_shell_rhf_one_step_runtests.jl`
  - passed: 38/38
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed: `load ok`
- `git diff --check`
  - passed

Git status:

```text
## main...origin/main
 M test/nested/pqs_multilayer_complete_core_shell_rhf_one_step_runtests.jl
```

Deletion/shrinkage report:
- deleted: none
- simplified: none; this was a focused convention cross-check.
- quarantined: RHF one-step remains private diagnostic-only through existing nonclaim assertions.
- not deleted because: existing one-step test is now the narrow convention guard for this private seam.
- exact remaining caller/blocker: no route-driver caller; further use remains blocked on a manager-approved private diagnostic wiring pass.

-- repo-doer@macmini

Pass 113 response - local compact PQS RHF SCF probe

Scope:
- No tracked source/test/doc edits.
- Added ignored local probe artifacts under `tmp/work/`:
  - `tmp/work/pqs_compact_route_owned_rhf_scf_probe.jl`
  - `tmp/work/pqs_compact_route_owned_rhf_scf_probe_summary.txt`

Command:
- `julia --project=. tmp/work/pqs_compact_route_owned_rhf_scf_probe.jl`

Runtime:
- Julia-level elapsed time from inside the script: `104.527284583` seconds.
- This exceeded 60 seconds, as expected and reported before launch, because the cold process compiled/built the compact route-owned final basis, H1/J density interaction, and private RHF SCF helpers.

Probe result:
- RHF input contract availability:
  `available_pqs_multilayer_complete_core_shell_rhf_input_contract`
- Initial density:
  `materialized_pqs_multilayer_complete_core_shell_rhf_initial_density_payload`
- Initial density trace:
  `3.9999999999999987`
- Initial one-step payload:
  `materialized_pqs_multilayer_complete_core_shell_rhf_one_step_payload`
- SCF status:
  `blocked_pqs_multilayer_complete_core_shell_rhf_scf_payload`
- SCF blocker:
  `scf_not_converged`
- Iterations:
  `25`
- Last total energy:
  `-10.032171604675229`
- Last one-body energy:
  `-14.21915637345403`
- Last two-body energy:
  `4.1869847687788013`
- Last density change:
  `4.4593642533563571e-05`
- Last energy change:
  `1.6560067273019285e-07`
- Final one-step density match error:
  not available because the SCF payload did not converge and did not perform the converged final recomputation.

Relation to H1/J diagnostic:
- H1 energy:
  `-5.6629907690725245`
- H1/J status:
  `materialized_pqs_multilayer_complete_core_shell_h1_j_payload`
- H1/J self-Coulomb:
  `1.8691288063594704`
- H1/J density gauge:
  `pre_final_localized_positive_weight`
- H1/J driver route materialized:
  `true`

Nonclaims:
- This is route-smoke diagnostic output, not physics acceptance.
- Explicit probe inputs were `electron_count = 4` and `fixture_role = :route_smoke`.
- Electron count was not inferred from charge.
- No damping/mixing was added after nonconvergence.
- No route-driver wiring, report fields, public API, exports, artifacts, GTO, IDA/MWG, or fixture promotion were added.

Validation/status:
- `git status --short --branch`

```text
## main...origin/main
```

Ignored artifact status:

```text
!! tmp/work/pqs_compact_route_owned_rhf_scf_probe.jl
!! tmp/work/pqs_compact_route_owned_rhf_scf_probe_summary.txt
```

Deletion/shrinkage report:
- deleted: none
- simplified: none; this was a local probe only.
- quarantined: RHF SCF remains private diagnostic/prototype and route-smoke-only.
- not deleted because: ignored probe artifacts are useful local evidence for the next corrective decision.
- exact remaining caller/blocker: private SCF does not converge on the compact route-owned PQS fixture within 25 iterations without damping/mixing; next decision is whether to audit iteration behavior, add a controlled damping/mixing seam, or adjust SCF controls in another local probe.

-- repo-doer@macmini

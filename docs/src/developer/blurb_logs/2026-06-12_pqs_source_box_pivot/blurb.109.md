Pass 109 - private closed-shell RHF SCF payload, tiny diagnostic loop

Baseline:

- Current pushed HEAD should include `7e9abf82 Add PQS RHF initial density payload`.
- Private RHF helpers live in `src/pqs_multilayer_complete_core_shell_rhf.jl`.
- Existing private helpers:
  - `_pqs_multilayer_complete_core_shell_rhf_input_contract(...)`
  - `_pqs_multilayer_complete_core_shell_rhf_initial_density_payload(...)`
  - `_pqs_multilayer_complete_core_shell_rhf_one_step_payload(...)`

Task:

Add the first private closed-shell RHF SCF payload over the existing private
diagnostic seams.

This is still private diagnostic/prototype behavior. It is not route adoption,
not a physics endpoint, and not a public API.

Preferred implementation surface:

- Stay in `src/pqs_multilayer_complete_core_shell_rhf.jl`.
- Add a private helper with a name like:
  `_pqs_multilayer_complete_core_shell_rhf_scf_payload(...)`.
- Inputs:
  - available RHF input contract;
  - materialized H1 payload, or extractable H1 payload via existing private
    helper;
  - materialized density interaction, or extractable density interaction via
    existing private helper;
  - optional initial-density payload;
  - small SCF controls such as `max_iterations`, `density_atol`,
    `energy_atol`, and optional metadata.
- If `initial_density_payload` is not supplied, it is acceptable to build one
  internally by calling the private initial-density helper.

SCF loop contract:

- Start from a spin-summed final density.
- On each iteration:
  - call `_pqs_multilayer_complete_core_shell_rhf_one_step_payload(...)`;
  - symmetrize the returned Fock/effective Fock matrix;
  - diagonalize it in the ordinary final-basis metric;
  - occupy the lowest `nocc` orbitals with occupancy 2;
  - build the next spin-summed density;
  - record compact iteration data: iteration number, total energy,
    density change, energy change, and convergence flag.
- Converged when density and energy tolerances are satisfied. Be explicit about
  how iteration 1 handles missing previous energy.
- Return a compact payload containing final density, final orbital coefficients
  or occupied coefficients, final one-step payload, compact iteration summary,
  and metadata.
- Summary must not duplicate large matrices beyond the top-level payload fields.

Expected status/blocker labels:

- converged/materialized:
  `:materialized_pqs_multilayer_complete_core_shell_rhf_scf_payload`
- blocked:
  `:blocked_pqs_multilayer_complete_core_shell_rhf_scf_payload`
- useful blockers:
  - `:missing_rhf_input_contract`
  - `:missing_h1_payload`
  - `:missing_density_interaction`
  - `:missing_initial_density`
  - `:invalid_scf_controls`
  - `:scf_not_converged`
  - or a blocker propagated from the private initial-density/one-step payloads.

Required nonclaims:

- `rhf_materialized = true` only for this private diagnostic payload if
  converged.
- `rhf_converged = true` only when the loop actually satisfies tolerances.
- `driver_route_materialized = false`
- `route_report_materialized = false`
- `exports_materialized = false`
- `artifacts_materialized = false`
- `public_api = false`
- fixture role carried through unchanged.

Tests:

- Add one focused synthetic test file or extend the RHF initial-density/one-step
  tests only if that stays clearer.
- Use a tiny diagonal or near-diagonal fixture where the H1-Aufbau density is
  already self-consistent, so the test is fast and deterministic.
- Cover:
  - converged private SCF payload;
  - final density trace equals electron count;
  - nonclaims remain false for driver/report/export/artifacts;
  - missing density interaction or missing contract blocker;
  - nonconvergence with `max_iterations = 1` only if easy and deterministic.

Exclusions:

- Do not wire the route driver.
- Do not add report aliases or driver options.
- Do not add public exports/API.
- Do not infer electron count from nuclei.
- Do not promote compact fixtures to physics acceptance.
- Do not touch GTO, IDA/MWG, exports, artifacts, or production route behavior.
- Do not run the heavy source-box dry-run.
- Do not request interactive escalation/approval; if approval is truly needed,
  write `ATTENTION.md` and stop.

Validation:

- `julia --project=. test/nested/<focused SCF test>.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Report back:

- Helper/object names.
- Status/blocker labels.
- SCF control defaults.
- Tiny-fixture convergence result.
- Validation commands/results.
- Git status.
- Deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini

Pass 112 - fix SCF final one-step/final-density consistency

Baseline:

- Current pushed HEAD should include `c3668502 Record PQS pass 111 audit`.
- Pass 111 found a labeling/consistency issue in
  `_pqs_multilayer_complete_core_shell_rhf_scf_payload(...)`.

Task:

Correct the private SCF payload so the returned final diagnostic energy and
`final_one_step_payload` correspond to the returned `final_density`.

Problem to fix:

- Inside the SCF loop, `one_step` is computed from the density at the start of
  the iteration.
- Then the loop diagonalizes the Fock matrix and updates `density` to
  `next_density`.
- The converged return currently returns `final_density = next_density` but
  `final_one_step_payload = one_step`, which was computed from the previous
  density. They are close under tolerance, but not exactly the same object
  convention.

Required implementation:

- In `src/pqs_multilayer_complete_core_shell_rhf.jl`, update
  `_pqs_multilayer_complete_core_shell_rhf_scf_payload(...)`.
- On converged return, after `density = next_density`, recompute
  `_pqs_multilayer_complete_core_shell_rhf_one_step_payload(...)` with
  `final_density = density`.
- Return that recomputed payload as `final_one_step_payload`.
- Use the recomputed payload for final energy fields such as
  `summary.final_total_energy`.
- Keep `iteration_records` unchanged as convergence history.
- Add compact summary fields such as:
  - `orbital_metric = :ordinary_orthonormal_final_basis`
  - `final_one_step_recomputed = true`
  - `final_one_step_density_matches_final_density = true`
  - optionally `converged_iteration_input_total_energy` if useful.
- For `:scf_not_converged`, either leave existing behavior but label it
  explicitly as last attempted/input one-step, or make an analogous compact
  label. Do not overbuild this path.

Tests:

- Update `test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl`.
- Assert:
  - `scf.final_one_step_payload.final_density ≈ scf.final_density`;
  - `scf.summary.final_one_step_recomputed === true`;
  - `scf.summary.orbital_metric === :ordinary_orthonormal_final_basis`;
  - final total energy equals `scf.final_one_step_payload.total_energy`;
  - route/report/export/artifact/public nonclaims remain false.

Exclusions:

- Do not change the RHF contraction formula.
- Do not add damping/mixing.
- Do not run a real compact PQS route probe yet.
- Do not wire the route driver.
- Do not add report aliases/options or public exports.
- Do not touch GTO, IDA/MWG, exports, artifacts, or production route behavior.
- Do not request interactive escalation/approval; if approval is truly needed,
  write `ATTENTION.md` and stop.

Validation:

- `julia --project=. test/nested/pqs_multilayer_complete_core_shell_rhf_scf_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Report back:

- Exact SCF summary labels added/changed.
- Validation commands/results.
- Git status.
- Deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini

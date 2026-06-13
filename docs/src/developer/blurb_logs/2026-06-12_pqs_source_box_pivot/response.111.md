Pass 111 response - private RHF SCF payload audit

Scope:
- No code changes.
- No Julia commands.
- Read-only code/contract inspection only.

Audit result:
- Not ready for a local real compact PQS probe yet.
- The private SCF loop is directionally consistent with the pass-102 RHF contract, but one payload-labeling issue should be corrected before using it on a real compact PQS fixture.

1. Energy and density consistency:
- `final_one_step_payload` is produced from the density at the start of the iteration.
- `final_density` is then updated from the eigensolve of that one-step Fock matrix.
- On a converged return, these densities are close by `density_atol`, so the current result is acceptable for the synthetic self-consistent fixture.
- For a real probe, the labels are not tight enough: `final_total_energy` is the energy of `final_one_step_payload.final_density`, not necessarily exactly the returned `final_density`.
- Smallest corrective pass: after convergence, either recompute one final one-step payload on the returned `final_density`, or explicitly label the existing payload as `converged_iteration_input_one_step_payload` and record the density mismatch. Recomputing is cleaner before a real probe.

2. Convergence semantics:
- The first-iteration rule is acceptable: `energy_change = nothing`, and energy convergence is treated as satisfied only because density convergence is still required.
- `:scf_not_converged` has the same labeling issue as the converged path: it returns the last attempted one-step payload plus the post-diagonalization density. That is useful, but should be labeled as last input one-step versus last output density, or recomputed if used as a final diagnostic.

3. Density/occupation convention:
- The loop preserves spin-summed final density mathematically: `P = 2 * C_occ * C_occ'`, so `tr(P) == electron_count` for the closed-shell input contract.
- Ordinary final-basis metric is assumed through `eigen(Symmetric(fock_matrix))`. It should be made explicit in the SCF summary before a real probe, e.g. `orbital_metric = :ordinary_orthonormal_final_basis`.
- Open-shell, odd electron count, fractional occupation, and inferred electron count remain blocked upstream by `_pqs_multilayer_complete_core_shell_rhf_input_contract(...)`.

4. Payload carrying cost:
- Large matrices are kept in top-level payload fields or nested top-level payloads, not duplicated in summaries/metadata.
- Iteration records are compact enough: iteration number, total energy, density change, energy change, density/energy convergence flags, and combined convergence.

5. Test pressure:
- The four current RHF test files are still justified as separate fast contract tests for now:
  - input contract;
  - initial density;
  - one-step Fock/energy convention;
  - SCF loop behavior.
- Do not consolidate yet. Consolidation becomes attractive after the route request object exists and the seams stabilize.

Ready/not ready:
- Not ready for a real compact local PQS probe until the final one-step payload/density labeling is corrected.

Exact corrective task:
- In `src/pqs_multilayer_complete_core_shell_rhf.jl`, adjust `_pqs_multilayer_complete_core_shell_rhf_scf_payload(...)` so the returned final energy and `final_one_step_payload` clearly correspond to the returned `final_density`.
- Preferred fix: after convergence, recompute `_pqs_multilayer_complete_core_shell_rhf_one_step_payload(...)` using the returned `final_density`, return that as `final_one_step_payload`, and keep iteration records unchanged as the convergence history.
- Add/adjust the focused SCF test to assert that `scf.final_one_step_payload.final_density ≈ scf.final_density`.
- Add summary fields for the ordinary final-basis metric assumption and the final one-step recomputation/consistency check.

Git status:

```text
## main...origin/main
```

Deletion/shrinkage forecast:
- deleted: none in the corrective pass.
- simplified: labels/summary semantics become clearer before route/request integration.
- quarantined: SCF remains private diagnostic/prototype, not route/report/public behavior.
- not deleted because: the four focused RHF tests still protect distinct active private contracts.
- exact remaining caller/blocker: real compact PQS probe is blocked by final SCF one-step/density consistency labeling.

-- repo-doer@macmini

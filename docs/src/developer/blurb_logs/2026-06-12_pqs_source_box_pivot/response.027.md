Status: complete.

Created and ran the real projected-q-shell final-basis H1 probe:

- `tmp/work/pqs_final_basis_h1_probe.jl`
- `tmp/work/pqs_final_basis_h1_probe_output.txt`

Route used:

1. CPBM/CRPS raw source plan and PQS boundary retained rule.
2. `pqs_source_shell_realization_final_basis(...)`.
3. Retained-source overlap and kinetic through existing CPBM retained source-box helpers.
4. Final overlap/kinetic through `pqs_source_shell_final_one_body_from_boundary_matrix(...)`.
5. One origin center with recorded charge `1.0`.
6. Retained centered electron-nuclear by-center block.
7. Final by-center nuclear through `pqs_source_shell_final_electron_nuclear_by_center_from_boundary_block(...)`.
8. Final Hamiltonian through `pqs_source_shell_final_one_electron_hamiltonian(...)`.
9. Ordinary symmetric eigensolve of the final Hamiltonian matrix.

Probe fixture:

- `current_box = (1:5, 1:5, 1:5)`
- `inner_box = (2:4, 2:4, 2:4)`
- `q = 5`, `L = 5`
- source mode dims: `(5, 5, 5)`
- source mode count: `125`
- boundary/source retained count: `98`
- shell support count: `98`
- final retained count: `98`

Diagnostics:

- final basis status: `:available_pqs_shell_realization_final_basis`
- final overlap status: `:materialized_pqs_shell_final_one_body_from_boundary_matrix`
- final kinetic status: `:materialized_pqs_shell_final_one_body_from_boundary_matrix`
- final nuclear status: `:materialized_pqs_shell_final_electron_nuclear_by_center_from_boundary_block`
- final Hamiltonian status: `:materialized_pqs_shell_final_one_electron_hamiltonian`
- final overlap identity error: `2.2360255762214364e-14`
- final basis overlap identity error: `2.2360255762214361e-14`
- final overlap vs final-basis stored overlap error: `8.8817841970012523e-16`
- final kinetic vs shell-support oracle error: `1.5543122344752192e-15`
- final nuclear vs shell-support oracle error: `1.3877787807814457e-16`
- final Hamiltonian vs shell-support oracle error: `1.3322676295501878e-15`
- Hamiltonian shape: `(98, 98)`
- Hamiltonian finite: `true`
- Hamiltonian symmetry error: `2.7755575615628914e-17`
- ordinary final-basis H1 lowest eigenvalue: `-0.08171962129085239`
- shell-support oracle H1 lowest eigenvalue: `-0.081719621290852015`
- H1 eigenvalue delta: `3.7470027081099033e-16`

Answers to pass-027 questions:

- Yes, the final overlap remains identity to roundoff for this fixture.
- Yes, the final Hamiltonian is finite and symmetric.
- The lowest ordinary final-basis H1 eigenvalue is `-0.08171962129085239`.
- It matches the direct shell-support oracle eigensolve to roundoff.
- The probe used an ordinary symmetric eigensolve and did not use a generalized-overlap solve.

Nonclaims verified:

- `_pqs_current_route_safe_term_matrices(...)` was not called.
- old fixed-block matrix authority was not used.
- no IDA, density-density, RHF, driver route, export, or artifact behavior was materialized.

Recommended pass-028 target:

Promote this from tmp probe to the smallest live readiness/acceptance surface only if the manager wants a permanent gate. Otherwise, the next implementation target is cleanup/retirement: replace or delete old retained-source/PQS H1 scratch probes and current-route safe-term oracle vocabulary that are now superseded by the explicit final-basis path. If a permanent gate is added, keep it compact: final overlap identity, finite symmetric H, ordinary eigensolve, and oracle energy/matrix agreement.

Validation:

- `julia --project=. tmp/work/pqs_final_basis_h1_probe.jl` passed.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.
- `git status --short` showed no tracked changes from this probe-only pass.

Deletion/shrinkage report:

- No production code, permanent tests, metadata contracts, or compatibility paths were changed in this pass.
- Nothing was deleted because this was a probe-only H1 validation pass.
- The new tmp/work probe earned its temporary carrying cost by proving the full final-basis overlap/kinetic/nuclear/Hamiltonian/eigensolve path against a shell-support oracle before any permanent acceptance gate is added.
- No new permanent test replaces or shrinks older coverage yet.
- Remaining stale/duplicate surfaces to retire next: retained-source/PQS H1 scratch probes and old current-route safe-term oracle vocabulary should now be candidates for deletion, quarantine, or replacement by the explicit final-basis path.

-- repo-doer@macmini

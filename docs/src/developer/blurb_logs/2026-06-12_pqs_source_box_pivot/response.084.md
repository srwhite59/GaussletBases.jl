Pass 084 complete: audited the PQS H1/J density seam and recorded the current convention boundary.

Files changed:

- `docs/src/developer/pqs_near_term_final_basis_realization_plan.md`

Helper/convention map:

- `CartesianFinalBasisRealization.pqs_complete_core_shell_final_ida_weights(final_basis, support_weights)` projects support-row weights through the final coefficients. It is useful diagnostic/final-weight data, but it is not the current owner of the H1/J density-interaction boundary.
- `CartesianFinalBasisRealization.pqs_complete_core_shell_pre_final_density_interaction(final_basis, support_pair_raw_numerator, support_weights)` owns the corrected H1/J density convention:
  - input pair data is the raw support-pair numerator;
  - support weights are positive support-row IDA weights in `core_then_shell` order;
  - the density interaction is built in the localized pre-final positive-weight gauge;
  - weight division happens at that pre-final density-interaction boundary;
  - it explicitly records `signed_final_weight_division_used = false` and `raw_no_division_used = false`.
- `CartesianFinalBasisRealization.pqs_complete_core_shell_pre_final_orbital_self_coulomb(density_interaction, final_orbital_coefficients)` owns the current final-orbital consumption rule for J: map the final H1 orbital into the pre-final density gauge with the combined Lowdin cleanup, then evaluate the one-orbital restricted direct-minus-exchange diagnostic.
- `pqs_multilayer_complete_core_shell_h1_payload(...)` now owns the H1 one-body assembly seam, but it intentionally does not return support-space density inputs or H1 eigenvectors.

Exact missing inputs for a route-owned PQS H1/J payload:

- Support weights over the complete core/shell support rows in the same `core_then_shell` order as `final_basis.core_support_states` followed by `final_basis.shell_support_states`. Existing probes compute these from axis PGDG weights as `wx[ix] * wy[iy] * wz[iz]`.
- The support raw pair numerator matrix over the same complete support ordering. Existing probes build this from `pgdg_intermediate.pair_factor_terms_raw` and Coulomb expansion coefficients by summing product-axis raw numerator factors.
- The complete final basis, including `pre_final_coefficients`, `final_coefficients`, and `combined_lowdin_cleanup`; this is already available from `pqs_multilayer_complete_core_shell_final_basis(...)`.
- Lowest H1 orbital coefficients in the final basis. Current `pqs_complete_core_shell_final_h1_solve(...)` reports eigenvalues/lowest energy but not the eigenvector, so the side-13 H1/J probes recompute the eigendecomposition locally from the final Hamiltonian matrix.

What remains probe-local or missing:

- `_support_weights(...)` and `_support_pair_raw_numerator_matrix(...)` are still probe-local in `tmp/work` scripts.
- The side-13 H1/J and RHF probes still perform a local eigendecomposition to get the lowest H1 orbital.
- There is no route-owned compact H1/J payload that consumes the H1 payload plus final basis, materializes the density interaction, and reports H1/J together.
- No production producer currently packages the support raw numerator and support weights as a complete core/shell density-input object.

Smallest next implementation pass:

Add a narrow route-owned H1/J payload that consumes:

- `pqs_multilayer_complete_core_shell_h1_payload(...)`;
- the complete core/shell final basis;
- route-owned support weights;
- route-owned support raw pair numerator.

It should compute or expose the lowest H1 final-basis orbital coefficients, call `pqs_complete_core_shell_pre_final_density_interaction(...)`, then call `pqs_complete_core_shell_pre_final_orbital_self_coulomb(...)`. If the support density-input producer is not ready, the first implementation should stop at that exact blocker rather than preserve the probe-local helpers as a production contract.

Probe/test glue that would shrink if that lands:

- `tmp/work/pqs_multilayer_shell_side13_j_probe.jl` and the RHF/core-ladder probes can stop carrying local support-weight and raw-pair-numerator construction for the H1/J stage.
- The probes can stop recomputing the lowest H1 eigensystem solely to recover the orbital coefficients if the route-owned payload exposes them.
- Future compact checks can assert one route-owned H1/J payload result instead of preserving separate probe-local density vocabulary.

Docs update:

- Updated `docs/src/developer/pqs_near_term_final_basis_realization_plan.md` to replace the stale “future H1 payload” wording with the implemented H1 payload seam.
- Added the reviewed H1/J density convention language: raw numerator first, then localized pre-final positive-weight density boundary; not signed-final-weight division and not raw no-division.

Validation:

- `git diff --check` passed.

Deletion/shrinkage report:

- No source code, tests, or production metadata were added.
- No old code was deleted because this pass was an audit/design pass and only corrected docs.
- The docs now shrink conceptual drift by no longer describing the H1 payload as future work and by naming the exact H1/J density seam still missing.
- No new test was added; no older coverage was replaced in this pass.
- Remaining stale or duplicate surfaces to retire next: probe-local `_support_weights`, `_support_pair_raw_numerator_matrix`, and local H1 orbital eigendecomposition in side-13 H1/J/RHF probes, once a route-owned H1/J payload lands.

-- repo-doer@macmini

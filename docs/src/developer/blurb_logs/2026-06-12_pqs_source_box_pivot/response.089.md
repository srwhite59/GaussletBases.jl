Pass 089 complete: added the route-owned PQS multilayer H1/J payload helper.

Files changed:
- `src/pqs_multilayer_complete_core_shell_h1.jl`
- `tmp/work/pqs_complete_core_shell_pre_final_density_consumption_probe.jl` (ignored developer probe)

Function added:
- `pqs_multilayer_complete_core_shell_h1_j_payload(plan; final_basis, h1_payload, axis_weights, raw_pair_factor_terms, coulomb_expansion, metadata=(;))`

What it does:
- validates the multilayer shell source plan, complete core/shell final basis, and H1 payload;
- builds support weights with `pqs_multilayer_support_weights(...)`;
- builds the support raw pair numerator with `pqs_multilayer_support_pair_raw_numerator_matrix(...)`;
- materializes the reviewed pre-final positive-weight density interaction;
- extracts the lowest final H1 orbital from `h1_payload.final_hamiltonian`;
- materializes the self-Coulomb record through `pqs_complete_core_shell_pre_final_orbital_self_coulomb(...)`;
- returns a compact payload with status, blocker, density interaction, self-Coulomb record, summary, and metadata.

Convention guardrails preserved:
- signed-final-weight division is not implemented or accepted;
- raw-no-division final density is not implemented or accepted;
- density-normalized `pair_factor_terms` are not used as authority;
- RHF, GTO, driver wiring, exports, artifacts, side-13 reruns, q ladders, and fixture tuning were not added.

Probe result:
- `h1_j_payload_status = materialized_pqs_multilayer_complete_core_shell_h1_j_payload`
- `h1_j_payload_orbital_source = h1_payload_final_hamiltonian_lowest_eigenvector`
- `h1_j_payload_energy_reconstruction_error = 9.4368957093138306e-16`
- `H1 = -0.48047934800387126`
- `J = 0.6397851751855723`
- fixed-oracle `J = 0.6397857768997106`
- `J - fixed_oracle = -6.017141382930191e-7`
- exact hydrogen 1s reference used by the probe: `0.625`

The accepted H1/J result matches the previous accepted summary within the existing probe tolerance. Rejected signed-final-weight and raw-no-division comparisons remain probe-local diagnostics only.

Validation:
- `julia --project=. tmp/work/pqs_complete_core_shell_pre_final_density_consumption_probe.jl` passed.
- `julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl` passed: 46 tests.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.
- Whitespace scan on edited source/probe files found no trailing whitespace.

Deletion/shrinkage report:
- The accepted tmp/work H1/J path no longer manually wires support weights, support raw numerator, density interaction, lowest-orbital extraction, and self-Coulomb as separate accepted-route steps.
- Support-density inputs were already route-owned in pass 088; this pass moved the accepted H1 orbital extraction and H1/J payload construction behind the route-owned helper.
- The remaining probe-local duplication is limited to rejected-convention comparisons and their local orbital extraction, so those diagnostics do not become production concepts.
- No permanent test was added. The new helper is covered by the tmp/work diagnostic probe plus the existing focused H1 gate; adding another test would preserve probe vocabulary rather than protect a new endpoint.
- Remaining before RHF or fixture-rule work: this is still a compact H1/J diagnostic seam only. A future RHF path still needs an explicit route-owned density/HF contract rather than extending this helper into an acceptance driver.

-- repo-doer@macmini

Pass 088 complete: added route-owned complete core/shell support-density input helpers and updated the tmp/work H1/J probe to consume them.

Files changed:

- `src/GaussletBases.jl`
- `src/pqs_multilayer_support_density.jl`
- `tmp/work/pqs_complete_core_shell_pre_final_density_consumption_probe.jl` (ignored developer probe)

Functions added:

- `pqs_multilayer_support_weights(plan; axis_weights)`
  - consumes an available `pqs_multilayer_shell_source_plan`;
  - uses `core_support_states` followed by `shell_support_states`;
  - computes product support weights as `wx[ix] * wy[iy] * wz[iz]`;
  - accepts a common 1D axis-weight vector or x/y/z axis weights.
- `pqs_multilayer_support_pair_raw_numerator_matrix(plan; raw_pair_factor_terms, coulomb_expansion)`
  - consumes the same complete support ordering;
  - contracts raw `pair_factor_terms_raw` with positive Coulomb expansion coefficients;
  - accepts a common term-first raw pair-factor array or x/y/z axis-specific raw terms;
  - does not divide by weights;
  - does not consume density-normalized `pair_factor_terms` as authority;
  - does not materialize final density interaction.

Probe update:

- `tmp/work/pqs_complete_core_shell_pre_final_density_consumption_probe.jl` now builds a matching explicit-box `pqs_multilayer_shell_source_plan(...)` and uses:
  - `pqs_multilayer_support_weights(...)`;
  - `pqs_multilayer_support_pair_raw_numerator_matrix(...)`.
- The probe-local normalized-pair comparison remains local because it is only a rejected-convention diagnostic.
- H1 orbital coefficient extraction remains probe-local.

H1/J probe result:

- `support_density_input_source = pqs_multilayer_support_density_helpers`
- `source_plan_status = available_pqs_multilayer_shell_source_plan`
- H1 energy: `-0.48047934800387126`
- density interaction status: `materialized_pqs_complete_core_shell_pre_final_density_interaction`
- density gauge: `pre_final_localized_positive_weight`
- pre-final weights positive / negative / near-zero: `223 / 0 / 0`
- pre-final pair matrix symmetry error: `9.436895709313831e-16`
- H1 self-Coulomb J: `0.6397851751855723`
- fixed oracle J: `0.6397857768997106`
- J delta vs fixed oracle: `-6.017141382930191e-7`
- exact hydrogen 1s J reference: `0.625`
- signed-final-weight division remains rejected: `18.132786333403647`
- raw no-division final density remains rejected: `0.5624216201593974`

The H1/J values match the previous accepted summary for this compact complete core/shell density-consumption probe.

Validation:

- `julia --project=. tmp/work/pqs_complete_core_shell_pre_final_density_consumption_probe.jl`
  - passed and wrote the summary above.
- `julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl`
  - passed: 46 tests.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed.
- `git diff --check`
  - passed.
- Whitespace scan on touched/new source and probe files found no trailing whitespace.

Deletion/shrinkage report:

- Tmp/work probe-local duplication became unnecessary for support weights and raw support-pair numerator construction in the updated density-consumption probe.
- No tracked test was added because this is input production for the reviewed H1/J seam, and the requested validation remains a developer probe plus the existing H1 gate.
- No old tracked coverage was shrunk in this pass.
- Remaining probe-local pieces:
  - H1 orbital coefficient extraction still recomputes the eigendecomposition locally;
  - rejected normalized-pair and signed-final-weight comparison paths remain probe-only diagnostics;
  - no route-owned H1/J payload yet ties H1 payload, final basis, support density inputs, and self-Coulomb into one compact result.

-- repo-doer@macmini

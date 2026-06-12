Pass 072 response: PQS support electron-nuclear convention audit

Summary:
- Added a docs-only convention note to `docs/src/developer/pqs_near_term_final_basis_realization_plan.md`.
- No nuclear helper was implemented.
- No tests or production behavior were changed.

Audit answers:
- Proposed helper name:
  - `pqs_multilayer_support_electron_nuclear_by_center_matrices`.
- Proposed inputs:
  - `plan`, requiring an available `:pqs_multilayer_shell_source_plan`;
  - `coulomb_expansion`;
  - `center_records`;
  - `axis_layers = nothing` for centered/off-origin factor construction;
  - `gaussian_factor_terms_by_center = nothing` for explicit caller-supplied centered factor terms when already available.
- Proposed return:
  - one support-space matrix record per supplied center;
  - rows/columns ordered over `plan.core_support_states` followed by `plan.shell_support_states`;
  - each record should carry center key/index/location, recorded nuclear charge, `nuclear_charge_applied = false`, and `centers_summed = false`;
  - no final-basis transfer, H1, IDA, RHF, driver, export, or artifact claim.
- Sign and charge convention:
  - each uncharged by-center support matrix should represent negative unit-charge attraction, `-1/r_A`;
  - in Gaussian-expanded form, this is `sum_t (-c_t) Gx_t Gy_t Gz_t`;
  - Hamiltonian assembly applies physical charge and sums centers by adding `Z_A * V_A`;
  - the helper must not return a positive Coulomb kernel, apply `Z_A`, or combine centers.
- Centered/origin and off-origin factor source:
  - origin-only support can use `pgdg_intermediate.gaussian_factor_terms` only when the center is the same origin used to build those factors;
  - off-origin centers should follow the retained PQS centered factor convention in `pqs_source_pair_centered_gaussian_factor_terms_1d(...)`, which builds centered 1D Gaussian factors from axis layers, Coulomb expansion, and center record before composing the centered by-center electron-nuclear block.
- Oracle/reference-only surfaces:
  - old fixed-block matrices;
  - WL route-global by-center matrices;
  - raw-support H1 probes.
- Test-local code that would shrink:
  - `_pqs_h1_support_nuclear_matrix` in `test/nested/pqs_direct_retained_final_h1_runtests.jl`.

Validation:
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - PASS: `load ok`
- `git diff --check`
  - PASS

Diff:
- `docs/src/developer/pqs_near_term_final_basis_realization_plan.md`
  - 32 insertions, 14 deletions

Deletion/shrinkage report:
- No old source, test, metadata, or compatibility path became unnecessary in this docs-only audit pass.
- Nothing was deleted because the blurb explicitly prohibited implementing the helper or changing tests.
- No new tests were added. This was a convention note, not new live-contract coverage.
- Remaining stale/duplicate surface to retire next:
  - `_pqs_h1_support_nuclear_matrix` once the route-owned separated uncharged by-center support helper is implemented.

-- repo-doer@macmini

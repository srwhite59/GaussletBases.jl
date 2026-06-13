Pass 181 complete.

Files changed:
- `test/nested/pqs_direct_retained_final_h1_runtests.jl`

Files deleted:
- `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`

What changed:
- Added the H1-J / pre-final density-interaction diagnostic smoke to the
  corrected one-center fixed-q He PQS fixture.
- Used the existing route-owned/internal density input convention:
  - `GaussletBases.CartesianContractedParentMetrics._pqs_source_box_ida_factor_provenance(...)`
  - `provenance.axis_weights`
  - `provenance.raw_axis_pair_factor_terms`
- Added compact assertions for:
  - H1-J payload materialization
  - final dimension `419`
  - H1 energy carried through the H1-J summary
  - small H1 energy reconstruction error
  - materialized pre-final density interaction
  - density gauge `:pre_final_localized_positive_weight`
  - finite pre-final pair matrix
  - all-positive pre-final weights
  - finite positive self-Coulomb
  - RHF/GTO/driver/export/artifact flags remain false

Observed values:
- final dimension: `419`
- H1 energy: `-1.9866819751748936`
- H1 orbital energy reconstructed inside H1-J: `-1.9866819751749145`
- H1 energy reconstruction error: `2.0872192862952943e-14`
- self-Coulomb: `1.2261626003119184`
- density gauge: `pre_final_localized_positive_weight`
- density-interaction status:
  `materialized_pqs_complete_core_shell_pre_final_density_interaction`
- pre-final pair matrix finite: `true`
- pre-final weights all positive: `true`

Raw provenance source:
- Provenance helper:
  `GaussletBases.CartesianContractedParentMetrics._pqs_source_box_ida_factor_provenance`
- Provenance interaction path: `ida_gausslet_source_box`
- H1-J consumed `provenance.axis_weights` and
  `provenance.raw_axis_pair_factor_terms`.
- It did not consume retained diagnostic weights, density-normalized pair terms
  as authority, signed-final-weight division, raw-no-division, or fixed-block
  oracle data.

Deleted Be2 scaffold:
- Deleted `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`.
- Reason: Be2/Cr2 artifact comparison work is paused; this was a standalone
  route/artifact fingerprint scaffold, not the current atom-first physics
  target. A caller check found no active `src`/`test` caller or default-runner
  include, only historical handoff/log references.

Validation:
- `julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl`
  - passed, `51/51`, `6.3s`
- `julia --project=. -e 'include("test/nested/pqs_direct_retained_final_h1_runtests.jl"); fixture = _pqs_h1_complete_fixture(); center = (; center_key = :origin, center_index = 1, location = (0.0, 0.0, 0.0), charge = 2.0); h1_payload = GaussletBases.pqs_multilayer_complete_core_shell_h1_payload(fixture.plan; final_basis = fixture.final_basis, coulomb_expansion = fixture.expansion, center_records = (center,), gaussian_factor_terms_by_center = fixture.bundle.pgdg_intermediate.gaussian_factor_terms, metadata = (; fixture = :pqs_fixed_q_he_h1j_gate_probe)); provenance = GaussletBases.CartesianContractedParentMetrics._pqs_source_box_ida_factor_provenance(fixture.plan.bundles; expected_term_count = length(fixture.expansion.coefficients)); h1_j_payload = GaussletBases.pqs_multilayer_complete_core_shell_h1_j_payload(fixture.plan; final_basis = fixture.final_basis, h1_payload, axis_weights = provenance.axis_weights, raw_pair_factor_terms = provenance.raw_axis_pair_factor_terms, coulomb_expansion = fixture.expansion, metadata = (; fixture = :pqs_fixed_q_he_h1j_gate_probe)); s = h1_j_payload.summary; d = h1_j_payload.density_interaction; println("final_dimension=", s.final_dimension); println("h1_energy=", s.h1_energy); println("h1_orbital_energy=", s.h1_orbital_energy); println("h1_energy_reconstruction_error=", s.h1_energy_reconstruction_error); println("self_coulomb=", s.self_coulomb); println("density_gauge=", d.density_gauge); println("density_interaction_status=", d.status); println("pre_final_pair_matrix_finite=", d.pre_final_pair_matrix_finite); println("pre_final_weights_all_positive=", d.pre_final_weights_all_positive); println("provenance_interaction_path=", provenance.diagnostics.interaction_path)'`
  - passed; printed the observed values above
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed; printed `load ok`
- `git diff --check`
  - passed
- `git diff --numstat -- src test`
  - `36  0    test/nested/pqs_direct_retained_final_h1_runtests.jl`
  - `0   164  test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
  - source/test net: `+36/-164`, net `-128`

Git status:

```text
## main...origin/main
 M test/nested/pqs_direct_retained_final_h1_runtests.jl
 D test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl
```

Deletion/shrinkage report:
- deleted: standalone Be2/Cr2 route-artifact fingerprint scaffold
  `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
- simplified: the corrected fixed-q He test now carries the current H1 plus
  H1-J diagnostic seam in one focused atom-first gate
- quarantined: none
- not deleted because: the focused He H1/H1-J gate remains the active guard for
  the corrected `419`-dimensional one-center PQS basis
- exact remaining caller/blocker: no active source/test caller or default-runner
  include found for the deleted Be2 scaffold; only historical handoff/log
  references remain

-- repo-doer@macmini

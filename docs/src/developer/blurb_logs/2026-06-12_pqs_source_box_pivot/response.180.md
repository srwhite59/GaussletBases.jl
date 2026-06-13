Pass 180 complete.

Files changed:
- `test/nested/pqs_direct_retained_final_h1_runtests.jl`

Files deleted:
- `test/nested/pqs_source_box_route_driver_complete_core_shell_ham_payload_runtests.jl`

What changed:
- Restored a focused one-center He H1 smoke on the corrected fixed-q PQS fixture.
- Kept the pass-179 source/final-basis inventory assertions:
  - parent side `11`
  - core support `125`
  - three shell layers
  - source dims `(5, 5, 5)` for each shell layer
  - retained per shell `98`
  - final retained dimension `419`
- Added H1 payload construction through
  `GaussletBases.pqs_multilayer_complete_core_shell_h1_payload(...)` with the
  one-center `Z = 2.0` He nuclear center.
- Asserted only the live H1 smoke contract:
  - payload materialized
  - final H1 matrix finite
  - final H1 matrix symmetric
  - ordinary symmetric solve
  - finite negative lowest H1 energy
  - final dimension `419`
  - no generalized-overlap solve
  - no density-density, RHF, export, or artifact materialization

Observed H1 values:
- final dimension: `419`
- lowest H1 energy: `-1.9866819751748936`
- matrix finite: `true`
- symmetry max abs error: `2.3092638912203256e-14`
- solve kind: `ordinary_symmetric`
- generalized-overlap solve materialized: `false`

Old compact 223 Ham scaffold:
- Deleted `test/nested/pqs_source_box_route_driver_complete_core_shell_ham_payload_runtests.jl`.
- Reason: it was a standalone old compact `223` complete-core/shell Ham
  payload scaffold; `rg` found no live `src`/`test` caller or default-runner
  include, only historical handoff/log references.

Validation:
- `julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl`
  - passed, `35/35`, `5.2s`
- `julia --project=. -e 'include("test/nested/pqs_direct_retained_final_h1_runtests.jl"); fixture = _pqs_h1_complete_fixture(); center = (; center_key = :origin, center_index = 1, location = (0.0, 0.0, 0.0), charge = 2.0); h1_payload = GaussletBases.pqs_multilayer_complete_core_shell_h1_payload(fixture.plan; final_basis = fixture.final_basis, coulomb_expansion = fixture.expansion, center_records = (center,), gaussian_factor_terms_by_center = fixture.bundle.pgdg_intermediate.gaussian_factor_terms, metadata = (; fixture = :pqs_fixed_q_he_h1_gate_probe)); h = h1_payload.final_hamiltonian.hamiltonian_matrix; println("final_dimension=", h1_payload.h1.final_dimension); println("h1_lowest_energy=", h1_payload.h1.lowest_energy); println("h1_matrix_finite=", all(isfinite, h)); println("h1_symmetry_error=", maximum(abs.(h - transpose(h)))); println("solve_kind=", h1_payload.h1.solve_kind); println("generalized_overlap_solve_materialized=", h1_payload.h1.generalized_overlap_solve_materialized)'`
  - passed; printed the observed values above
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed; printed `load ok`
- `git diff --check`
  - passed
- `git diff --numstat -- src test`
  - `35  0    test/nested/pqs_direct_retained_final_h1_runtests.jl`
  - `0   115  test/nested/pqs_source_box_route_driver_complete_core_shell_ham_payload_runtests.jl`
  - source/test net: `+35/-115`, net `-80`

Git status:

```text
## main...origin/main
 M test/nested/pqs_direct_retained_final_h1_runtests.jl
 D test/nested/pqs_source_box_route_driver_complete_core_shell_ham_payload_runtests.jl
```

Deletion/shrinkage report:
- deleted: stale compact `223` Ham payload scaffold
  `test/nested/pqs_source_box_route_driver_complete_core_shell_ham_payload_runtests.jl`
- simplified: `test/nested/pqs_direct_retained_final_h1_runtests.jl` is now the
  single focused fixed-q `419` inventory-plus-H1 smoke for this one-center He
  target
- quarantined: none
- not deleted because: the focused H1 test remains the live guard for the
  corrected fixed-q PQS basis
- exact remaining caller/blocker: no remaining active caller found for the
  deleted 223 scaffold; historical handoff/log references remain only as record

-- repo-doer@macmini

# Cartesian Driver Endpoint Manifest

This manifest records current driver-owned endpoint roles. It is a compact
taxonomy aid, not a route framework.

| Driver input | Test file | Route role | Dimension | Endpoint ready | Current blocker | Default runner |
| --- | --- | --- | ---: | --- | --- | --- |
| `test/driver_inputs/he_pqs_q5_wlmap.jl` | `test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_runtests.jl` | He PQS q5 WL-map physical atom endpoint | 419 | yes | none | no |
| `test/driver_inputs/he_pqs_q5_wlmap.jl` | `test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_rhf_runtests.jl` | He PQS q5 WL-map private RHF diagnostic/endpoint check | 419 | yes | none | no |
| `test/driver_inputs/h2_pqs_q5_independent_source_box_r4.jl` | none | independent H2 PQS q5 source-box target/readiness input; physics flags default off | 471 expected | no | readiness-only input; use stage variants for final-basis/H1/H1-J diagnostics | no |
| `test/driver_inputs/h2_pqs_q5_independent_source_box_r4_final_basis.jl` | none | independent H2 PQS q5 final-basis diagnostic | 471 | no | `:missing_physical_gausslet_h1_builder` | no |
| `test/driver_inputs/h2_pqs_q5_independent_source_box_r4_h1.jl` | none | independent H2 PQS q5 H1 diagnostic | 471 | no | `:missing_physical_gausslet_h1_j_builder` | no |
| `test/driver_inputs/h2_pqs_q5_independent_source_box_r4_h1_j.jl` | none | independent H2 PQS q5 H1-J density diagnostic; not solver-ready | 471 | no | `:missing_physical_gausslet_rhf_or_solver_contract` | no |
| `test/driver_inputs/h2_pqs_q5_independent_source_box_r4_private_rhf.jl` | none | independent H2 PQS q5 private RHF diagnostic; not public solver/export-ready | 471 | no | `:private_rhf_diagnostic_not_public_solver_contract` | no |
| `test/driver_inputs/h2_fake_pqs_q5_wl_source_backed_r4.jl` | `test/nested/cartesian_ham_builder_h2_fake_pqs_wl_source_backed_r4_runtests.jl` | fake-PQS source-backed WL/QW 463 reproduction; not independent PQS | 463 | no, fake reproduction only | `:fake_pqs_source_backed_wl_reproduction_not_independent_pqs` | no |

# Cartesian Driver Endpoint Manifest

This manifest records current driver-owned endpoint roles. It is a compact
taxonomy aid, not a route framework.

| Driver input | Test file | Route role | Dimension | Endpoint ready | Current blocker | Default runner |
| --- | --- | --- | ---: | --- | --- | --- |
| `test/driver_inputs/he_pqs_q5_wlmap.jl` | `test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_runtests.jl` | He PQS q5 WL-map physical atom endpoint | 419 | yes | none | no |
| `test/driver_inputs/he_pqs_q5_wlmap.jl` | `test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_rhf_runtests.jl` | He PQS q5 WL-map private RHF diagnostic/endpoint check | 419 | yes | none | no |
| `test/driver_inputs/h2_pqs_q5_source_box_diagnostic_r4.jl` | `test/nested/cartesian_ham_builder_h2_pqs_q5_source_box_diagnostic_r4_runtests.jl` | H2 PQS q5 source-box diagnostic, not a physics endpoint | 221 | no | `:retained_atom_core_interiors_missing` | no |
| `test/driver_inputs/h2_fake_pqs_q5_wl_source_backed_r4.jl` | `test/nested/cartesian_ham_builder_h2_fake_pqs_wl_source_backed_r4_runtests.jl` | fake-PQS source-backed WL/QW 463 reproduction; not independent PQS | 463 | yes, as reproduction only | none | no |
| planned | planned | independent H2 PQS q5 physical gausslet endpoint | TBD | no | `:missing_independent_pqs_shared_shell_source_plan` | no |

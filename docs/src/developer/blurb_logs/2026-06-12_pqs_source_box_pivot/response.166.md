Pass 166 response - private Be2 PQS CR2 inspection artifact writer

Status: complete
Role: repo-doer@macmini
Commit: not committed; repo-manager owns commits.

Files changed:
- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`

Implemented:
- Added private `_pqs_source_box_route_driver_be2_cr2_inspection_bundle_payload(assembly)`.
- Added private `_pqs_source_box_route_driver_write_be2_cr2_inspection_bundle(jld2_path, tsv_path, payload)`.
- Added private TSV cell helper `_pqs_source_box_route_driver_tsv_cell`.
- The writer stores plain JLD2-friendly arrays/scalars/NamedTuple metadata and TSV rows only; it does not store private `_PQS...` structs, full reports, or opaque route objects.
- The focused Be2 fingerprint test now writes JLD2/TSV under `mktempdir()` and checks compact artifact/fingerprint behavior.
- Kept H1 scientific checks for dimension, finiteness, symmetry, and lowest energy.
- Removed duplicate private ham-input summary/metadata false nonclaim fields that had no current callers; CR2/nonclaim readiness remains carried by the consumer/readiness payload.

JLD2 top-level/group keys written:
- `schema/name`
- `schema/version`
- `bundle/purpose`
- `producer/package`
- `producer/repo_commit`
- `producer/dirty`
- `producer/generated_at`
- `routes/pqs_source_box/route/*`
- `routes/pqs_source_box/readiness/*`
- `routes/pqs_source_box/system/*`
- `routes/pqs_source_box/final_basis/*`
- `routes/pqs_source_box/one_body/hamiltonian`
- `routes/pqs_source_box/h1/lowest_energy`
- `routes/pqs_source_box/two_body/*`
- `routes/pqs_source_box/validation/*`
- `routes/white_lindsey/route/*`
- `routes/white_lindsey/readiness/*`

TSV columns written:
`route_label`, `status`, `blocker`, `final_dimension`, `pre_final_dimension`, `support_weight_count`, `one_body_shape`, `two_body_shape`, `h1_lowest`, `h1_symmetry_defect`, `one_body_finite`, `two_body_finite`, `density_gauge`, `raw_pair_factor_convention`, `cr2_read_only_inspector_ready`, `cr2_solver_ready`, `cr2_export_ready`, `cr2_handoff_blocker`, `nuclear_repulsion`, `electron_count`, `spin_sector`.

Validation results:
- `julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl` passed: 19 passed, 51.3 s.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed: `load ok`.
- `git diff --check` passed.
- `git diff --numstat -- src test`:
  - `118  17  src/pqs_source_box_diatomic_complete_core_shell.jl`
  - `63   165 test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
  - total added 181, deleted 182; net-negative by 1 line.

Git status:
`## main...origin/main`
` M src/pqs_source_box_diatomic_complete_core_shell.jl`
` M test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`

Deletion/shrinkage report:
- deleted: duplicate private ham-input summary/metadata false nonclaim fields (`h1_j_materialized`, `ham_payload_materialized`, `rhf_materialized`, public/export/artifact flags, and retained diagnostic weight nonclaim copy) where they had no current caller.
- simplified: focused Be2 fingerprint test now checks the artifact/fingerprint seam instead of many adjacent scalar readiness/handoff assertions; local fixture wrapper now returns only assembly.
- quarantined: White-Lindsey is present only as unavailable/not-applicable placeholder route metadata in this PQS-first private inspection artifact; no WL solver/export claim was added.
- not deleted because: H1 dimension/finite/symmetry/lowest-energy checks remain the live scientific contract; CR2/nonclaim readiness remains on the consumer readiness payload.
- exact remaining caller/blocker: CR2 solver/export handoff remains blocked by `:missing_cr2_solver_handoff_format`; downstream solver/HamV6/HFDMRG readiness remains blocked by `:missing_hfdmrg_density_density_contract` and related downstream contract objects.

-- repo-doer@macmini

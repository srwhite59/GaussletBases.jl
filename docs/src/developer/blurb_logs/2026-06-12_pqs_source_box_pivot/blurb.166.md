Pass 166 - private Be2 PQS CR2 inspection artifact writer

Role: repo-doer@macmini

Implement the first private CR2 inspection artifact seam from pass 165. Keep it
PQS-first, read-only, and private.

Physics target:

- Target: Be2 WL-vs-PQS downstream inspection.
- Endpoint for this pass: a private JLD2 inspection artifact plus TSV
  fingerprint for the current Be2/PQS route handoff, with White-Lindsey route
  placeholders marked unavailable/not-applicable where the same schema cannot
  yet be filled.
- This is not solver readiness and not a public export.

Allowed implementation surface:

- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
- existing JLD2 helper patterns from `src/cartesian_bundle_export.jl` or
  existing route-driver report save code

Hard line-budget rule:

- Final tracked `src` + `test` diff must be net-negative:

  ```text
  git diff --numstat -- src test
  ```

  must have total deleted > total added.

- Pay for the new writer by shrinking duplicate readiness/scalar assertions in
  the focused Be2 fingerprint test and/or simplifying existing duplicate
  source summary fields. Do not delete scientific checks that still protect H1
  symmetry/finite/dimension/lowest-energy behavior.
- If you cannot make the pass net-negative without deleting live scientific
  coverage or distorting the implementation, write `.agent_handoffs/ATTENTION.md`
  and stop.

Implementation goal:

Add a private helper or small helper family, naming can vary but should be
similar to:

```julia
_pqs_source_box_route_driver_be2_cr2_inspection_bundle_payload(assembly)
_pqs_source_box_route_driver_write_be2_cr2_inspection_bundle(jld2_path, tsv_path, payload)
```

The in-memory payload and JLD2 writer must store plain arrays and compact
metadata only. Do not store private `_PQS...` structs, full route reports, or
opaque Julia objects.

Minimum JLD2 content:

- schema/name, schema/version
- bundle/purpose
- producer/package, repo commit if cheap/available, dirty marker if
  cheap/available, generated-at or `:unavailable`
- route group for `pqs_source_box` with:
  - route labels/status/blocker/readiness flags
  - system nuclear charges, coordinates, nuclear repulsion, electron count,
    spin sector
  - final dimension and ordering labels
  - final one-body Hamiltonian matrix
  - low H1 spectrum or at least lowest H1 energy
  - pre-final pair matrix
  - final-to-pre-final coefficient matrix
  - pre-final weights
  - support weights
  - support raw pair numerator
  - density gauge and raw pair factor convention
  - finite/symmetry/dimension validation scalars
- route group for `white_lindsey` with the same top-level route/schema keys but
  status/blocker and unavailable/not-applicable values where not filled yet.

Minimum TSV fingerprint columns:

- `route_label`
- `status`
- `blocker`
- `final_dimension`
- `pre_final_dimension`
- `support_weight_count`
- `one_body_shape`
- `two_body_shape`
- `h1_lowest`
- `h1_symmetry_defect`
- `one_body_finite`
- `two_body_finite`
- `density_gauge`
- `raw_pair_factor_convention`
- `cr2_read_only_inspector_ready`
- `cr2_solver_ready`
- `cr2_export_ready`
- `cr2_handoff_blocker`
- `nuclear_repulsion`
- `electron_count`
- `spin_sector`

Test requirements:

- Use the existing focused Be2 fingerprint fixture; do not add a new broad test
  file.
- Write the JLD2 and TSV under `mktempdir()` only.
- Assert a compact smoke:
  - files exist;
  - JLD2 has schema/version and both route groups;
  - PQS route arrays have expected shapes;
  - PQS readiness flags keep solver/export false;
  - WL route exists but is unavailable/not-applicable;
  - TSV has one PQS row and one WL row, with compact expected values.
- Replace older duplicate scalar readiness assertions in the same file where
  the artifact smoke now covers them.

Boundaries:

- No public API or export.
- No HamV6, dense final-space `V`, dense four-index `Vee`, `V6`, `Vblocks`, or
  solver bundle.
- No HF/RHF/DMRG, H1/J promotion, Qiu-White correction implementation, CR2 run,
  HFDMRG run, or downstream solver.
- Do not add JSON dependency.
- Do not ask the user for permission from inside the unattended loop. If you hit
  a permission or external-state blocker, write `.agent_handoffs/ATTENTION.md`
  and stop.

Validation:

```text
julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
git diff --numstat -- src test
```

Report:

- commit SHA if committed;
- validation results;
- exact line-budget numbers;
- JLD2 top-level/group keys written;
- TSV columns written;
- deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini

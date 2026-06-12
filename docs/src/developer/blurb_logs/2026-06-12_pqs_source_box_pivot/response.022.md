Ran the PQS shell-support one-body source-boundary audit probe.

Artifacts:

- `tmp/work/pqs_shell_support_operator_boundary_audit_probe.jl`
- `tmp/work/pqs_shell_support_operator_boundary_audit_probe_output.txt`

What the probe did:

- Built the real cubic projected-q-shell fixture from pass 020:
  - `current_box = (1:5, 1:5, 1:5)`
  - `inner_box = (2:4, 2:4, 2:4)`
  - `q = 5`, `L = 5`
- Used old `_pqs_shell_realization_plan(...)` and `_pqs_product_box_support_overlap_matrix(...)` only as oracle/input sources.
- Built a CPBM/CRPS raw source plan and PQS boundary retained rule.
- Built existing CPBM retained-source overlap and kinetic blocks.
- Built shell-support overlap and kinetic oracle matrices directly from support states and axis metrics.
- Compared retained-source blocks against:

```text
O_boundary = P' * O_shell_support * P
```

- Fed the shell-support kinetic operator through `pqs_source_shell_projected_one_body_matrix(...)`.
- Did not call `_pqs_current_route_safe_term_matrices(...)`.
- Did not materialize H1, charge summing, electron-nuclear, IDA, density-density, RHF, driver routes, exports, artifacts, full-parent fallback, direct Cartesian fallback, or old fixed-block matrix authority.

Key diagnostics:

- `descriptor_kind = :projected_q_shell`
- `source_mode_dims = (5, 5, 5)`
- `source_mode_count = 125`
- `retained_count = 98`
- `shell_support_count = 98`
- `final_retained_count = 98`

Overlap:

- `retained_source_overlap_vs_shell_boundary_error = 8.8817841970012523e-16`
- `shell_boundary_overlap_final_error = 2.2204460492503131e-16`
- `retained_source_overlap_lowdin_final_error = 8.8817841970012523e-16`
- `retained_source_overlap_identity_error = 1.3322676295501878e-15`
- `shell_boundary_overlap_identity_error = 1.3322676295501878e-15`

Kinetic:

- `retained_source_kinetic_vs_shell_boundary_error = 1.3322676295501878e-15`
- `kinetic_projection_status = :materialized_pqs_shell_projected_one_body_matrix`
- `kinetic_projection_blocker = nothing`
- `kinetic_lowdin_boundary_crosscheck_error = 8.8817841970012523e-16`
- `kinetic_final_operator_shape = (98, 98)`

Answers:

- Retained-source overlap does equal the shell-projected boundary overlap for this real cubic PQS fixture, to roundoff.
- Lowdin-transforming that boundary overlap reproduces the final identity overlap, to roundoff.
- The same statement is demonstrated for kinetic, using a direct shell-support kinetic oracle built from support states and axis metrics.
- This means the existing CPBM retained-source overlap/kinetic blocks can be treated as boundary operators for this fixture; they are not merely full raw-source blocks that would be invalid under Lowdin cleanup.

Exact next blocker:

```text
:missing_pqs_final_one_body_from_boundary_operator
```

Recommended pass-023 implementation target:

Add the narrow CPBM helper that consumes an existing retained-source boundary operator block and the PQS final-basis object, then computes:

```text
O_final = L' * O_boundary * L
```

Start with overlap and kinetic only. Keep H1 blocked after that on the by-center electron-nuclear boundary/source term:

```text
:missing_pqs_shell_boundary_electron_nuclear_operator_source
```

Do not generate Hamiltonians or run H1 in that pass.

Validation:

- `julia --project=. tmp/work/pqs_shell_support_operator_boundary_audit_probe.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Deletion/shrinkage report:

- No tracked source, docs, or tests changed in this audit-only pass.
- Nothing was deleted yet because this only establishes that retained-source overlap/kinetic are valid boundary operators for the cubic fixture.
- The ignored `tmp/work` probe earned its temporary carrying cost by preventing a wrong implementation path: it shows Lowdin-transforming retained-source overlap/kinetic is valid here because those retained-source blocks equal `P' O_shell P` to roundoff, not because Lowdin alone is generally a raw-to-final transform.
- Remaining stale surfaces to retire next: once final overlap/kinetic-from-boundary materialization exists, bridge/readiness tests that only preserve `:shell_realization_not_materialized` can shrink further; retained-source H1 probes should remain diagnostic until the electron-nuclear boundary operator is available.

-- repo-doer@macmini

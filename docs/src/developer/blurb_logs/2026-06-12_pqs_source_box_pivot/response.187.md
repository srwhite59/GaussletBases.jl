Pass 187 complete: no-edit driver-owned He 419 PQS seam audit.

Summary:

- The visible driver spine in `bin/cartesian_ham_builder.jl` should be preserved as-is: defaults, optional include/override, then explicit calls to `cartesian_system`, `cartesian_recipe`, `cartesian_parent`, `cartesian_shells`, `cartesian_units`, `cartesian_transforms`, `cartesian_pair_terms`, `cartesian_assembly`, `cartesian_report`, `cartesian_materialization`, and `cartesian_save`.
- The current driver can partially describe the He 419 PQS case, but it cannot yet express the exact WL parent mapping used by the focused He gate.
- The complete-core/shell H1/H1-J route-owned construction already belongs in `cartesian_assembly`; the missing piece is a stable, compact physics artifact emitted by `cartesian_save`.

Can current driver variables express the He 419 PQS case?

Partially, not exactly.

Current usable driver names:

- `atom_symbols = ("He",)`
- `nuclear_charges = (2,)`
- `atom_locations = ((0.0, 0.0, 0.0),)`
- `route_family = :pqs_source_box`
- `route_kind = ...` is accepted as metadata, but no special `:one_center_fixed_q_complete_core_shell` behavior is keyed on it yet.
- `parent_axis_counts = (x = 11, y = 11, z = 11)` is the current spelling, not `parent_axis_count`.
- `parent_axis_probe_family = :G10` is the current axis-family input, not `parent_axis_family`.
- `q = 5`, `n_s = 5`, `reference_spacing = 1.0`, `tail_spacing = 10.0`.

Missing or mismatched names/seams:

- No current `system_kind` input; one-center classification is inferred from one center.
- No current `mapping_rule`, `Z`, or `d` input is consumed by the one-center parent-axis builder.
- `_cartesian_one_center_parent_basis_object` hard-codes `IdentityMapping()` and records `mapping = :IdentityMapping`.
- Existing `white_lindsey_Z` is only a materializer/export input for the White-Lindsey path; it does not drive PQS parent-axis construction.
- `run_h1`, `run_h1_j`, and `run_private_rhf` are not current driver controls. The PQS complete-core/shell H1/H1-J diagnostic payload is built unconditionally in `cartesian_assembly` for `route_family = :pqs_source_box`.

Smallest naming extension:

- Add parent mapping inputs to the driver defaults and `parent_inputs`, using the current naming style:
  - `parent_axis_family = :G10` as a clearer alias for `parent_axis_probe_family`, while preserving the existing name.
  - `parent_mapping_rule = :identity_mapping` or `:white_lindsey_atomic_mapping`.
  - `parent_mapping_Z = 2.0`
  - `parent_mapping_d = 0.3`
  - `parent_mapping_tail_spacing = tail_spacing`
- Teach `_cartesian_one_center_parent_basis_object` to choose `white_lindsey_atomic_mapping(Z = parent_mapping_Z, d = parent_mapping_d, tail_spacing = parent_mapping_tail_spacing)` for that rule.

Where H1/H1-J should be materialized:

- Keep construction in `cartesian_assembly`.
- Current assembly already builds `_PQSCompleteCoreShellDiagnosticRoutePayload` through:
  - `_pqs_source_box_route_driver_complete_core_shell_source_plan_payload`
  - `_pqs_source_box_route_driver_complete_core_shell_h1_payload`
  - `_pqs_source_box_route_driver_complete_core_shell_density_inputs`
  - `_pqs_source_box_route_driver_complete_core_shell_h1_j_diagnostic_payload`
  - `_pqs_source_box_route_driver_complete_core_shell_ham_payload`
- Do not move this into an opaque driver wrapper.
- Add a compact plain report/artifact payload derived from `assembly.complete_core_shell_diagnostic_route_payload.complete_core_shell_ham_payload`. `cartesian_report` should expose the small artifact payload; `cartesian_materialization` should not reconstruct route internals from scratch.

Current save/materialization state:

- `cartesian_save` delegates to `_pqs_source_box_route_driver_save`.
- `save_artifact`, `save_tsv`, `outfile`, and `tsvfile` already write a durable report/materialization JLD2 and TSV.
- `save_basis_artifact`, `save_ham_artifact`, `basisfile`, and `hamfile` are existing basis/Ham bundle knobs for current materializer paths.
- For `route_family = :pqs_source_box`, `_pqs_source_box_route_driver_materialization` still returns `:pending_source_box_retained_route`; it does not currently write a PQS physics artifact.
- The report currently exposes H1-J summary aliases:
  - `complete_core_shell_h1_j_final_dimension`
  - `complete_core_shell_h1_j_h1_energy`
  - `complete_core_shell_h1_j_self_coulomb`
  - `complete_core_shell_h1_j_density_gauge`
  - `complete_core_shell_h1_j_driver_route_materialized`
- Missing from stable artifact/report aliases: final-basis inventory details, retained-per-shell tuple, final overlap identity error, WL reference constants/deltas, and driver provenance fields.

Recommended first artifact layout:

Use `outfile` with `save_artifact = true` and add stable scalar/group keys beside the durable report, rather than asking tests to deserialize private route payloads.

Suggested minimal keys:

```text
config/input_path
config/route_family
config/route_kind
config/parent_mapping_rule
config/parent_mapping_Z
config/parent_mapping_d
config/tail_spacing
config/q
config/n_s

provenance/repo_commit
provenance/repo_dirty

basis/final_dimension
basis/core_support_count
basis/shell_support_count
basis/shell_layer_count
basis/retained_per_shell
basis/shell_final_retained_count
basis/final_overlap_identity_error

physics/h1_lowest
physics/h1_j_self_coulomb

density_interaction/status
density_interaction/density_gauge
density_interaction/raw_pair_factor_convention

comparison/reference_label
comparison/wl_h1_lowest
comparison/wl_h1_self_coulomb
comparison/delta_h1
comparison/delta_h1_j
```

Do not make final-basis self-overlap downstream working data. Keep only `basis/final_overlap_identity_error` as an identity-noise diagnostic.

Recommended driver test shape:

- Add checked-in input file `test/driver_inputs/he_pqs_q5_wlmap.jl`.
- Add a thin test that invokes `bin/cartesian_ham_builder.jl` with that input and temporary `outfile`/`tsvfile` overrides.
- Read `outfile` with JLD2 and assert:
  - final dimension `419`
  - shell layer count `3`
  - retained per shell `(98, 98, 98)`
  - final overlap identity error `< 1.0e-10`
  - H1 lowest `-1.991334820314074`
  - H1-J self-Coulomb `1.2420423900074902`
  - density gauge `:pre_final_localized_positive_weight`
  - WL H1 delta `9.649649361120893e-6`
  - WL J delta `-4.997485057112172e-6`

What the driver test should replace/delete:

- Once the driver artifact test exists, the hand-built construction in `test/nested/pqs_direct_retained_final_h1_runtests.jl` becomes obsolete:
  - `_pqs_h1_test_bundle`
  - `_pqs_h1_complete_fixture`
  - direct `shellify`/`lower_terminal_regions`/region-plan/source-plan/final-basis construction
  - direct H1/H1-J payload calls
  - direct WL delta assertions
- Because `test/nested/runtests.jl` currently includes `pqs_direct_retained_final_h1_runtests.jl`, the replacement pass should either delete that include and file, or shrink the file to a very small module-contract check only if a real non-driver contract remains uncovered.
- My recommendation: after the driver test covers the same physics and inventory artifact keys, delete the direct file and replace its include with the driver-owned test. Keep small internal tests only for isolated module contracts elsewhere.

Specific first implementation pass:

Files to edit:

- `bin/cartesian_ham_builder.jl`
- `src/pqs_source_box_route_driver_helpers.jl`
- `src/pqs_source_box_route_driver_reporting.jl`
- `test/nested/runtests.jl`

Input file to add:

- `test/driver_inputs/he_pqs_q5_wlmap.jl`

Test to add or replace:

- Add `test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_runtests.jl`.
- Remove `include("pqs_direct_retained_final_h1_runtests.jl")` from `test/nested/runtests.jl` and include the new driver-owned test.
- Delete `test/nested/pqs_direct_retained_final_h1_runtests.jl` once the new driver test passes.

Artifact keys to write:

- The minimal `config/*`, `provenance/*`, `basis/*`, `physics/*`, `density_interaction/*`, and `comparison/*` keys listed above.

Deletion/shrinkage target:

- Delete the direct He PQS test file after replacing it with the driver-owned artifact test. That should more than pay for the new small input file and thin test, keeping source/test/generator net-negative.

Validation command for the first coding pass:

```text
julia --project=. test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_runtests.jl
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

No tracked source/test/bin/generator changes were made in this audit.

Validation run:

```text
git status --short --branch
```

returned:

```text
## main...origin/main
```

```text
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

returned no output.

Deletion/shrinkage report:

- deleted: none in this no-edit audit.
- simplified: none in this no-edit audit.
- quarantined: none.
- not deleted because: pass 187 was explicitly no-edit; the direct He PQS test remains until a driver-owned artifact test can replace it.
- exact remaining caller/blocker: `test/nested/runtests.jl` includes `test/nested/pqs_direct_retained_final_h1_runtests.jl`; replacement is blocked on adding WL parent mapping inputs and a stable PQS physics artifact writer.

-- repo-doer@macmini

Pass 188 complete.

Files changed:

- `bin/cartesian_ham_builder.jl`
- `src/pqs_source_box_route_driver_helpers.jl`
- `src/pqs_source_box_route_driver_reporting.jl`
- `test/driver_inputs/he_pqs_q5_wlmap.jl`
- `test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_runtests.jl`
- `test/nested/runtests.jl`
- deleted `test/nested/pqs_direct_retained_final_h1_runtests.jl`

Driver input names added/used:

- `parent_axis_family`
- `parent_mapping_rule`
- `parent_mapping_Z`
- `parent_mapping_d`
- `parent_mapping_tail_spacing`
- `comparison_reference_label`
- `wl_h1_lowest`
- `wl_h1_self_coulomb`
- `input_path` in `save_inputs`, populated from the included driver input path

Parent mapping consumption:

- `bin/cartesian_ham_builder.jl` now carries the parent mapping inputs into `parent_inputs`.
- `_cartesian_one_center_parent_basis_object` consumes them through `_cartesian_one_center_parent_mapping`.
- Default remains `parent_mapping_rule = :identity_mapping`.
- `parent_mapping_rule = :white_lindsey_atomic_mapping` builds the parent axis with:

```julia
white_lindsey_atomic_mapping(
    Z = Float64(parent_mapping_Z),
    d = Float64(parent_mapping_d),
    tail_spacing = Float64(parent_mapping_tail_spacing),
)
```

Driver input added:

```text
test/driver_inputs/he_pqs_q5_wlmap.jl
```

It sets the He one-center PQS q=5/n_s=5 route, WL parent mapping parameters, artifact/TSV saving, and WL reference constants.

Artifact keys written:

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

Driver test:

- Added `test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_runtests.jl`.
- It runs `bin/cartesian_ham_builder.jl` with `test/driver_inputs/he_pqs_q5_wlmap.jl`, writes to a temp artifact, and checks the required artifact values.
- Replaced the old default nested include with the new driver-owned test.
- Deleted the direct hand-built H1/H1-J test.

Validation:

```text
julia --project=. test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_runtests.jl
```

passed:

```text
cartesian_ham_builder He PQS q5 WL-map artifact | 10 pass, 10 total, 4m07.4s
```

Runtime note: the focused driver test exceeded 60 seconds because this source edit triggered package precompilation and first-call compilation. The run printed about 58 seconds of precompilation; transform, pair terms, and assembly reported 12.4s, 25.2s, and 28.8s respectively, all mostly compilation time.

```text
julia --project=. -e 'using GaussletBases; println("load ok")'
```

passed with `load ok`.

```text
git diff --check
```

passed with no whitespace errors.

```text
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

reported:

```text
17	4	bin/cartesian_ham_builder.jl
65	6	src/pqs_source_box_route_driver_helpers.jl
36	4	src/pqs_source_box_route_driver_reporting.jl
22	0	test/driver_inputs/he_pqs_q5_wlmap.jl
35	0	test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_runtests.jl
0	163	test/nested/pqs_direct_retained_final_h1_runtests.jl
1	1	test/nested/runtests.jl
```

Line budget:

- added: `176`
- deleted: `178`
- net: `-2`

```text
git status --short --branch
```

returned:

```text
## main...origin/main
 M bin/cartesian_ham_builder.jl
 M src/pqs_source_box_route_driver_helpers.jl
 M src/pqs_source_box_route_driver_reporting.jl
 A test/driver_inputs/he_pqs_q5_wlmap.jl
 A test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_runtests.jl
 D test/nested/pqs_direct_retained_final_h1_runtests.jl
 M test/nested/runtests.jl
```

Deletion/shrinkage report:

- deleted: `test/nested/pqs_direct_retained_final_h1_runtests.jl`, including its hand-built parent/shellification/lowering/source-plan/final-basis/H1/H1-J fixture.
- simplified: the He q=5/n_s=5 WL-mapped PQS gate is now driver-owned: input file -> visible staged driver -> JLD2 artifact -> thin readback checks.
- quarantined: none.
- not deleted because: the new compact artifact writer remains needed so tests do not deserialize private route payload structs; driver defaults and mapping aliases remain needed for old-style include/override usage.
- exact remaining caller/blocker: no remaining tracked caller of `pqs_direct_retained_final_h1_runtests.jl`; `test/nested/runtests.jl` now includes the driver-owned He PQS artifact test.

-- repo-doer@macmini

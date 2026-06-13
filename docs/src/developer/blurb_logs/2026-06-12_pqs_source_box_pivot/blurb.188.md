Pass 188 - implement the driver-owned He 419 PQS artifact gate.

Purpose:

Move the accepted He q=5/n_s=5 PQS H1/H1-J physics gate out of the hand-built
test and into the visible driver path. The target is:

```text
driver input -> bin/cartesian_ham_builder.jl -> saved artifact -> read/check artifact
```

not:

```text
test reconstructs parent/shellification/source-plan/final-basis/H1/H1-J
```

Keep the old `nestpgg.jl` laboratory-script style. The modern driver should
remain visibly staged:

```text
cartesian_system
cartesian_recipe
cartesian_parent
cartesian_shells
cartesian_units
cartesian_transforms
cartesian_pair_terms
cartesian_assembly
cartesian_report
cartesian_materialization
cartesian_save
```

Do not replace that with an opaque `run_driver(config)` frontend.

Physics target:

```text
He, Z = 2
white_lindsey_atomic_mapping(Z=2,d=0.3,tail_spacing=10)
q = n_s = 5
fixed 5^3 core + three fixed-q shell sectors
final dimension = 419
retained per shell = (98, 98, 98)
H1 lowest = -1.991334820314074
H1-J self-Coulomb = 1.2420423900074902
density gauge = :pre_final_localized_positive_weight
WL H1 delta = +9.649649361120893e-6
WL J delta = -4.997485057112172e-6
```

Current state:

- The accepted focused physics gate is currently in:

  ```text
  test/nested/pqs_direct_retained_final_h1_runtests.jl
  ```

  That file still directly builds the parent basis, shellification/lowering,
  source plan, final basis, H1, and H1-J.

- Pass 187 found that the current driver can describe much of the case, but the
  one-center parent builder still hard-codes `IdentityMapping()`.

- The relevant driver/source surfaces are:

  ```text
  bin/cartesian_ham_builder.jl
  src/pqs_source_box_route_driver_helpers.jl
  src/pqs_source_box_route_driver_reporting.jl
  test/nested/runtests.jl
  test/nested/pqs_direct_retained_final_h1_runtests.jl
  ```

Exact task:

1. Add the smallest driver input support needed for the WL parent mapping.

   Use current driver naming style, but add clearer aliases where useful:

   ```julia
   parent_axis_family = :G10
   parent_axis_probe_family = :G10        # existing spelling remains accepted
   parent_axis_counts = (x = 11, y = 11, z = 11)

   parent_mapping_rule = :white_lindsey_atomic_mapping
   parent_mapping_Z = 2.0
   parent_mapping_d = 0.3
   parent_mapping_tail_spacing = tail_spacing
   ```

   Teach `_cartesian_one_center_parent_basis_object` or its immediate input
   path to use `white_lindsey_atomic_mapping(...)` for
   `parent_mapping_rule = :white_lindsey_atomic_mapping`.

   Preserve the old identity behavior as the default:

   ```julia
   parent_mapping_rule = :identity_mapping
   ```

2. Add a small checked-in driver input file:

   ```text
   test/driver_inputs/he_pqs_q5_wlmap.jl
   ```

   It should be readable and old-driver-like, not a giant config object. Use
   plain variable assignments. It should set:

   ```julia
   atom_symbols = ("He",)
   nuclear_charges = (2,)
   atom_locations = ((0.0, 0.0, 0.0),)

   route_family = :pqs_source_box
   route_kind = :one_center_fixed_q_complete_core_shell

   parent_axis_family = :G10
   parent_axis_counts = (x = 11, y = 11, z = 11)
   parent_mapping_rule = :white_lindsey_atomic_mapping
   parent_mapping_Z = 2.0
   parent_mapping_d = 0.3
   tail_spacing = 10.0
   parent_mapping_tail_spacing = tail_spacing

   q = 5
   n_s = 5
   fixed_source_mode_shape = true

   save_artifact = true
   save_tsv = true
   ```

   If a name is not consumed by the live driver and cannot be made meaningful
   without broad refactoring, stop and report the exact name/blocker instead of
   adding dead config fields.

3. Add a compact PQS physics artifact section to the existing save/report path.

   Use `outfile` with `save_artifact = true`. Do not require tests to deserialize
   private route payload structs.

   Minimal JLD2 keys to write:

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

   Provenance fields such as commit/dirty are useful if there is already a
   small local helper. Do not add a broad shelling-out provenance subsystem in
   this pass.

   Final-basis self-overlap policy:

   - write only `basis/final_overlap_identity_error` as a diagnostic;
   - do not pass or preserve a final-basis `S` matrix as downstream working
     data;
   - downstream final working bases are assumed orthonormal except for numerical
     identity noise.

4. Add the thin driver-owned test:

   ```text
   test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_runtests.jl
   ```

   The test should:

   - create a temporary output directory;
   - run `bin/cartesian_ham_builder.jl` with
     `test/driver_inputs/he_pqs_q5_wlmap.jl`;
   - pass temporary `outfile` and `tsvfile` overrides;
   - read the JLD2 artifact;
   - assert only the live physics/driver contract.

   Required checks:

   ```text
   basis/final_dimension == 419
   basis/shell_layer_count == 3
   basis/retained_per_shell == (98, 98, 98)
   basis/final_overlap_identity_error < 1.0e-10
   physics/h1_lowest ≈ -1.991334820314074
   physics/h1_j_self_coulomb ≈ 1.2420423900074902
   density_interaction/density_gauge == :pre_final_localized_positive_weight
   comparison/delta_h1 ≈ 9.649649361120893e-6
   comparison/delta_h1_j ≈ -4.997485057112172e-6
   ```

   Keep this test small. Do not reassert every intermediate payload field,
   nonclaim flag, blocker tuple, or route-shadow vocabulary item.

5. Replace the old direct test.

   Once the new driver-owned test passes:

   - remove `include("pqs_direct_retained_final_h1_runtests.jl")` from
     `test/nested/runtests.jl`;
   - add the new driver-owned test include;
   - delete `test/nested/pqs_direct_retained_final_h1_runtests.jl`.

   If you find a specific module-contract check in the direct file that remains
   live and is not covered by the driver artifact test, report the exact
   contract and caller. Do not preserve the whole hand-built route test by
   default.

Trust boundary:

- Do not add H2, Be2, Cr2, HFDMRG, DMRG, exports, public API, or public solver
  behavior.
- Do not add or wire private RHF driver behavior in this pass.
- Do not use supplemented WL 447 as the comparison baseline.
- Do not make the old one-shell 223/561/1115/1933 q-ladder the target.
- Do not create an opaque `run_driver(config)` wrapper.
- Do not add a second large fingerprint test.
- Do not request interactive approval or sandbox escalation during the loop. If
  approval is genuinely required, write `.agent_handoffs/ATTENTION.md` and stop.

Line-budget rule:

This pass edits implementation/test/driver surfaces, so it must be net-negative
for tracked source/test/bin/generator lines.

Measure exactly:

```text
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Acceptance condition:

```text
sum(deleted) > sum(added)
```

The intended deletion that should pay for the new driver input/test is:

```text
test/nested/pqs_direct_retained_final_h1_runtests.jl
```

If the pass cannot remain line-negative without deleting live scientific checks,
write `.agent_handoffs/ATTENTION.md` with the exact line-budget blocker and
stop. Do not preserve the old hand-built test and add the new driver test
beside it.

Validation:

Run focused validation only:

```text
julia --project=. test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_runtests.jl
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
git status --short --branch
```

If the driver test is expected to exceed 60 seconds, say so in the response
with the observed or estimated phase that dominates. Do not run broad nested
or integration suites for this pass.

Report:

- exact driver input names added/used;
- where the parent mapping rule is consumed;
- artifact keys written;
- driver test runtime;
- whether the direct hand-built test was deleted;
- line-budget arithmetic: added, deleted, net;
- validation results;
- deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

Write the result to `.agent_handoffs/response.188.md` and copy it to:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.188.md
```

-- repo-manager@macmini

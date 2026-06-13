Pass 187 - audit the driver-owned He 419 PQS handoff seam, no code yet.

Purpose:

We are pivoting the He q=5/n_s=5 PQS gate toward the driver style the user
wants. The model is the old laboratory-driver style in:

```text
/Users/srw/Dropbox/bin/nestpgg.jl
```

That script has the useful shape:

```text
defaults at top
include parameter file
command-line overrides
visible staged construction
@show / timing / diagnostics can be inserted anywhere
optional HF/ED/save blocks
stable JLD2 output
```

The current modern driver already has the same visible staged spine:

```text
bin/cartesian_ham_builder.jl

system = cartesian_system(...)
recipe = cartesian_recipe(...)
parent = cartesian_parent(...)
shells = cartesian_shells(...)
units = cartesian_units(...)
transforms = cartesian_transforms(...)
pairs = cartesian_pair_terms(...)
assembly = cartesian_assembly(...)
report = cartesian_report(...)
materialization = cartesian_materialization(...)
cartesian_save(...)
```

Do not replace that with an opaque `run_driver(config)` wrapper. The driver
must remain a visible laboratory script.

Current tracked He gate:

```text
test/nested/pqs_direct_retained_final_h1_runtests.jl
```

It still constructs the parent basis, shellification, lowering, source plan,
final basis, H1, and H1-J directly inside the test. That was right while
debugging the seam. The next target is to make the same case driver-owned:

```text
He, Z = 2
white_lindsey_atomic_mapping(Z=2,d=0.3,tail_spacing=10)
q = n_s = 5
fixed 5^3 core + 3 fixed-q shells
final dimension = 419
H1 lowest = -1.991334820314074
H1-J self-Coulomb = 1.2420423900074902
WL H1 delta = +9.649649361120893e-6
WL J delta = -4.997485057112172e-6
density gauge = :pre_final_localized_positive_weight
```

Task type:

No-edit audit/design seam. Do not implement yet.

Exact task:

Audit how to make the current He 419 PQS gate run through
`bin/cartesian_ham_builder.jl` and a small input file, while keeping the visible
driver stages intact.

Read/inspect these surfaces:

```text
bin/cartesian_ham_builder.jl
src/pqs_source_box_route_driver_helpers.jl
test/nested/pqs_direct_retained_final_h1_runtests.jl
/Users/srw/Dropbox/bin/nestpgg.jl
```

Focus on these questions:

1. Can the current driver input variables express a one-center He fixed-q
   complete-core/shell route with the WL parent mapping?

   Proposed small input file shape:

   ```text
   test/driver_inputs/he_pqs_q5_wlmap.jl
   ```

   Intended parameters:

   ```julia
   system_kind = :one_center_atom
   atom_symbols = ("He",)
   nuclear_charges = (2,)
   atom_locations = ((0.0, 0.0, 0.0),)

   route_family = :pqs_source_box
   route_kind = :one_center_fixed_q_complete_core_shell

   parent_axis_family = :G10
   parent_axis_count = 11
   mapping_rule = :white_lindsey_atomic_mapping
   Z = 2.0
   d = 0.3
   tail_spacing = 10.0

   q = 5
   n_s = 5
   fixed_source_mode_shape = true

   run_h1 = true
   run_h1_j = true
   run_private_rhf = false

   save_artifact = true
   save_tsv = true
   ```

   If the current driver uses different parameter names, report the exact
   current names and the smallest naming extension needed.

2. Where should H1/H1-J be materialized in the driver path?

   Current direct test calls:

   ```text
   pqs_multilayer_shell_source_plan
   pqs_multilayer_complete_core_shell_final_basis
   pqs_multilayer_complete_core_shell_h1_payload
   pqs_multilayer_complete_core_shell_h1_j_payload
   ```

   Identify whether `cartesian_assembly`, `cartesian_materialization`, or a
   small materialization subpayload should own this for the driver.

3. What should `cartesian_save` write for this first physics gate?

   Be precise about existing output concepts:

   ```text
   outfile
   tsvfile
   basisfile
   hamfile
   save_artifact
   save_tsv
   save_basis_artifact
   save_ham_artifact
   ```

   Recommend the smallest artifact layout that a thin test can read. It should
   contain enough for:

   ```text
   final dimension
   shell retained counts
   final overlap identity error
   H1 lowest
   H1-J self-Coulomb
   density gauge
   WL H1/J reference constants and deltas
   provenance: driver input path, commit, dirty flag if already available
   ```

   Do not ask for final-basis self-overlap as downstream working data. It can
   remain a diagnostic identity-noise check only.

4. What should the eventual driver test replace/delete?

   The desired end state is:

   ```text
   driver input -> driver output artifact -> read/check artifact
   ```

   not:

   ```text
   test reconstructs route internals
   ```

   Identify which parts of `pqs_direct_retained_final_h1_runtests.jl` become
   obsolete once the driver artifact test exists, and whether the whole file can
   eventually be deleted or whether a smaller module-contract fragment should
   remain.

5. What is the likely first implementation pass?

   Report a concrete next coding pass, but do not code it in pass 187. Include:

   ```text
   files to edit
   input file to add
   artifact keys to write
   test to add or replace
   deletion/shrinkage target to keep source/test/generator net-negative
   validation command
   ```

Trust boundary:

- Do not edit `src`, `test`, `bin`, docs, tracked driver inputs, or generators
  in this pass.
- Do not create a new opaque `run_driver(config)` function.
- Do not add H2, Be2, Cr2, HFDMRG, DMRG, exports, or artifacts.
- Do not add private RHF driver behavior yet.
- Do not make supplemented WL 447 part of this gate.
- Do not run long driver/materializer jobs. If you need a quick probe, keep it
  ignored under `tmp/work` and report why; otherwise stay read-only.

Line-budget/diff rule:

This is no-edit. These should be empty:

```text
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
git status --short --branch
```

Validation:

```text
git status --short --branch
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Report:

- whether current driver variables can express the He 419 PQS case;
- exact missing driver/materialization/save seams;
- recommended artifact key layout;
- recommended test shape;
- specific deletion/shrinkage target for the first coding pass;
- explicit confirmation that the visible staged driver spine is preserved;
- no tracked source/test/bin/generator changes.

Write the result to `.agent_handoffs/response.187.md` and copy it to:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.187.md
```

-- repo-manager@macmini

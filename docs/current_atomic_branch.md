# Current Atomic Branch

This page is the shortest current status read for the atomic line.

## 1. What the atomic branch is today

The present atomic story is:

1. the mature numerical radial workflow
2. an explicit one-electron `(l,m)` layer
3. static He / IDA-style interacting ingredients
4. direct / exchange / Fock helpers
5. a minimal UHF kernel

This is already a coherent small atomic line.

It is **not** yet a broad general atomic HF workflow.

## 2. Start here for the atomic line

Read these first:

- [`docs/first_radial_workflow.md`](first_radial_workflow.md)
- [`docs/recommended_atomic_setup.md`](recommended_atomic_setup.md)
- [`docs/atomic_ylm_layer.md`](atomic_ylm_layer.md)
- [`docs/atomic_ida_layer.md`](atomic_ida_layer.md)
- [`docs/atomic_ida_uhf.md`](atomic_ida_uhf.md)
- [`docs/hamiltonian_export_fullida_dense.md`](hamiltonian_export_fullida_dense.md)
- [`docs/hamiltonian_export_sliced_blocks.md`](hamiltonian_export_sliced_blocks.md)

Then run:

1. `examples/04_hydrogen_ground_state.jl`
2. `examples/15_atomic_hydrogen_ylm.jl`
3. `examples/16_atomic_ida_ingredients.jl`
4. `examples/22_atomic_ida_uhf.jl`
5. `examples/31_atomic_fullida_dense_export.jl`
6. `examples/32_atomic_sliced_export.jl`

## 3. What counts as current workflow documentation

These are the current-status pages for the atomic line:

- [`docs/atomic_ylm_layer.md`](atomic_ylm_layer.md)
- [`docs/atomic_ida_layer.md`](atomic_ida_layer.md)
- [`docs/atomic_ida_uhf.md`](atomic_ida_uhf.md)
- [`docs/hamiltonian_export_fullida_dense.md`](hamiltonian_export_fullida_dense.md)
- [`docs/hamiltonian_export_sliced_blocks.md`](hamiltonian_export_sliced_blocks.md)

Those five pages give the shortest current interpretation of the atomic
branch.

## 4. Supporting notes for the atomic line

The supporting-note stack is now grouped here:

- [`docs/atomic_mean_field_supporting_notes.md`](atomic_mean_field_supporting_notes.md)

Two other narrower atomic supporting notes remain useful:

- [`docs/gaunt_backend_note.md`](gaunt_backend_note.md)
- [`docs/atomic_angular_sectorization.md`](atomic_angular_sectorization.md)

These notes explain how the current atomic story was assembled, but they
should not be read as equally primary with the current workflow pages above.

## 5. Current interpretation

The current wording discipline for the atomic line is:

- radial is the mature numerical substrate
- `(l,m)` is the first explicit atomic angular layer
- IDA ingredients are explicit static data, not yet a broad solver framework
- direct / exchange / Fock / UHF are narrow current-model layers

If you want the broader package context after that, go back to:

- [`docs/index.md`](index.md)
- [`docs/example_guide.md`](example_guide.md)

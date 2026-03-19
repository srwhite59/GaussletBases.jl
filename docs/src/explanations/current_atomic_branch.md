# Current atomic branch

This page is the shortest current status read for the atomic line.

## What the atomic branch is today

The present atomic story is:

1. the mature numerical radial workflow
2. an explicit one-electron `(l,m)` layer
3. static He / IDA-style interacting ingredients
4. direct / exchange / Fock helpers
5. a minimal UHF kernel
6. dense and sliced Hamiltonian export for the current density-density model

This is already a coherent small atomic line. It is not yet a broad general
atomic HF framework.

## Start here

Within the current Documenter site, the best entry path is:

- [First radial workflow](../tutorials/first_radial_workflow.md)
- [Recommended atomic setup](../howto/recommended_atomic_setup.md)
- [Example guide](../howto/example_guide.md)

## Notes that are not yet migrated into the first site

The following current-workflow atomic notes still live in the flat `docs/`
tree and are not yet part of this first curated site:

- `atomic_ylm_layer.md`
- `atomic_ida_layer.md`
- `atomic_ida_uhf.md`
- `hamiltonian_export_fullida_dense.md`
- `hamiltonian_export_sliced_blocks.md`

The supporting-note chain for the atomic line also still lives in the flat
tree for now.

## Current interpretation

The wording discipline for the atomic line remains:

- the radial branch is the mature numerical substrate
- the `(l,m)` layer is the first explicit atomic angular layer
- the present interaction model is density-density / IDA, not a four-index
  Coulomb Hamiltonian
- direct / exchange / Fock / UHF are narrow current-model layers
- solver-facing export is already supported for the atomic line

That is the present atomic scope the library can stand behind.

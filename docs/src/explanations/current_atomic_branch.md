# Current atomic branch

This page is the shortest current user-facing status read for the mature
radial/atomic line.

It is the clearest current realization in the package of the main gausslet
idea from the papers: an atom-centered basis with orthonormality, locality,
variable resolution, and a compressed two-index or diagonal-style Coulomb
representation.

The most recent posted paper on this line is the radial gausslet manuscript on
arXiv:

- <https://doi.org/10.48550/arXiv.2603.22646>

## What the atomic branch is today

The present atomic story is:

1. the mature numerical radial workflow
2. an explicit one-electron `(l,m)` layer
3. static He / IDA-style interacting ingredients
4. direct / exchange / Fock helpers
5. a minimal UHF kernel
6. dense and sliced Hamiltonian export for the current density-density model

Scientifically, this branch adapts gausslets to atomic coordinates while
keeping the present electron-electron structure in a compact IDA form rather
than a four-index Coulomb tensor. The one-body side remains variational within
the basis, while the compression story is centered on the interaction.

This is already a coherent small atomic line. It is not yet a broad general
atomic HF framework.

## Angular research track

The next active scientific direction is the manuscript-facing angular gausslet
line. That branch is now tracked explicitly, but only as an experimental
research scaffold:

- [Angular research track](angular_research_track.md)

That page is the right place to read the current boundary:

- shell-local injected angular basis construction has **not** yet been ported
- Hooke is deferred for a later dedicated workflow line
- the near-term target is an atomic angular benchmark ladder:
  HF, small ED, and one DMRG-facing bridge

## Start here

Within the manual, the best entry path is:

- [First radial workflow](../tutorials/first_radial_workflow.md)
- [Recommended atomic setup](../howto/recommended_atomic_setup.md)
- [Example guide](../howto/example_guide.md)

## Notes that are not yet migrated into the first site

The following current-workflow atomic notes still live in the flat `docs/`
tree and are not yet part of the user-facing manual surface:

- `atomic_ylm_layer.md`
- `atomic_ida_layer.md`
- `atomic_ida_uhf.md`
- `hamiltonian_export_fullida_dense.md`
- `hamiltonian_export_sliced_blocks.md`

The supporting-note chain for the atomic line also still lives in the flat
tree for now. It should be treated as developer/supporting material rather
than as the main package manual.

## Current interpretation

The wording discipline for the atomic line remains:

- the radial branch is the mature numerical substrate
- the `(l,m)` layer is the first explicit atomic angular layer
- the present interaction model is density-density / IDA, not a four-index
  Coulomb Hamiltonian
- direct / exchange / Fock / UHF are narrow current-model layers
- solver-facing export is already supported for the atomic line

That is the present atomic scope the library can stand behind.

If you want the lower-priority architecture and supporting-note context after
that, continue with:

- [Developer Notes](../developer/index.md)

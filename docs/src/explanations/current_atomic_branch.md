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
line. That branch is now tracked explicitly, but it remains experimental and
manuscript-facing rather than part of the mature atomic user workflow:

- [Angular research track](angular_research_track.md)

That page is the right place to read the current boundary:

- shell-local injected angular basis construction and shell-to-atom assembly
  are now present in the repo
- the angular line already has one-electron, HF-style, and small-ED benchmark
  paths plus a direct in-memory HFDMRG payload handshake
- the exact common low-`l` reference already fits the current HamIO / HamV6
  consumer language cleanly, but the full mixed basis does not yet
- Hooke is deferred for a later dedicated workflow line
- the angular line should still be read as experimental, not as a mature
  public atomic workflow branch

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
- the newer injected-angular benchmark line is real, but it belongs to the
  separate experimental angular research track rather than to this mature
  atomic branch page

That is the present atomic scope the library can stand behind.

If you want the lower-priority architecture and supporting-note context after
that, continue with:

- [Developer Notes](../developer/index.md)

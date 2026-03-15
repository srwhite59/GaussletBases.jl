> **Note for new users:** this is a narrow research note about the current 1D contraction / hierarchy direction.  
> It is **not** the best place to start if you are new to GaussletBases.  
> Start instead with `README.md`, `docs/first_radial_workflow.md`, and `docs/recommended_atomic_setup.md`.

# Leaf-local 1D PGDG note

This note defines a very small PGDG-oriented Stage-2 slice.

The scope is intentionally narrow:

- 1D only
- hierarchy-driven local basis generation
- no atom or molecule geometry logic
- no historical nested-driver port
- no solver work

## Goal

The goal is to show that the new hierarchy layer can drive basis generation
cleanly.

This is *not* yet meant to reproduce the full historical PGDG or nested
workflows. It is only meant to prove that:

- a `HierarchicalBasisPartition1D` can serve as geometric input
- each leaf box can generate its own local Gaussian family
- the leaf-local families can be combined into one global construction object
- the existing matrix and representation layers still work on top of it

## Minimal generated object

The first generated object is:

- `LeafLocalPGDG1D`

It stores:

- the input hierarchy
- the leaf boxes used for generation
- a leaf-to-primitive index map
- the global `PrimitiveSet1D`
- an identity contraction matrix

In this first slice, the generated basis functions are just the generated
primitives themselves, so the contraction map is the identity.

That is enough to feed:

- `primitive_set(...)`
- `stencil_matrix(...)`
- `basis_metadata(...)`
- `basis_representation(...)`
- primitive overlap / position / kinetic builders

## Simplifying rules

The first version uses a deliberately simple local rule per leaf:

- choose a fixed number of local primitives per leaf
- place their centers evenly inside the leaf interval
- choose one local Gaussian width from the leaf size and local spacing

This means:

- refining a box increases local basis resolution because it creates more leaves
- untouched leaves generate exactly the same local structure as before

That is the main structural property we want to validate now.

## What this does not do yet

This slice does not yet include:

- adaptive refinement based on physics or errors
- Gaussian or GTO augmentation
- geometry-aware atom or molecule grouping
- historical nested basis-generation drivers

Those are later forks, not part of this first hierarchy-driven generation proof.

## Narrow Gaussian augmentation layer

The next small extension on top of `LeafLocalPGDG1D` is optional user-supplied
Gaussian augmentation.

The supported idea is still deliberately simple:

- start from an existing `LeafLocalPGDG1D`
- add extra Gaussian primitives to selected leaves, or to every leaf by a
  simple repeated rule
- keep everything Gaussian so the analytic-friendly overlap / position /
  kinetic route remains available

The first augmentation API is:

- `LeafGaussianSpec1D`
- `augment_leaf_pgdg(generator; by_leaf=..., every_leaf=...)`

This is still 1D only, still hierarchy-driven, and still not a named basis-set
library.

The generated object records which primitives came from:

- the original leaf-local generator
- the augmentation layer

so a downstream consumer can distinguish generated and augmented structure
without needing a permanent file format.

This still should be read as a useful prototype path for generation,
provenance, and hierarchy plumbing, not as the corrected conceptual path to the
historical nested construction. The corrected direction is described separately
in `docs/global_map_local_contraction.md`.

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

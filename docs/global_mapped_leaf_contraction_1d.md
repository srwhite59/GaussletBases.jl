# Global mapped layer and leaf contraction note

This note defines a narrow corrected implementation target.

The scope is intentionally small:

- 1D only
- one globally mapped common primitive layer over a region
- local contraction on leaf boxes only
- no parent-shell contraction yet
- no atom or molecule geometry logic
- no historical driver port
- no solver work

## Two-layer picture

The corrected path has two distinct layers.

### 1. Global mapped primitive layer

This is one common basis over the whole region, defined in a single global
coordinate map.

It should be useful in two ways:

- as the shared substrate for later local contraction
- as a simple basis or intermediate object that can already be used directly
  and may already be a reasonable endpoint in some calculations

So this layer is not just scaffolding.

### 2. Local leaf contraction layer

This is an optional refinement built on top of the global layer.

For each leaf box in an existing hierarchy:

- restrict the global layer to the primitives assigned to that leaf
- build local submatrices from the global primitive matrices
- do one local contraction there
- retain the resulting local functions

The local contraction does not define its own coordinate map. It only acts on a
submatrix of the globally defined primitive layer.

This first contraction layer is leaf-only. Parent-shell contraction is
deliberately deferred to a later step.

## Minimal object choices

The first narrow implementation can use objects in the spirit of:

- `GlobalMappedPrimitiveLayer1D`
- `LeafBoxContraction1D`
- `LeafBoxContractionLayer1D`

The global object should store:

- the global coordinate map
- the common `PrimitiveSet1D`
- the identity contraction map for the uncontracted basis
- basic metadata and representation support

The local contraction object should store:

- the leaf box it belongs to
- the primitive indices used locally
- the local contraction matrix
- the retained local centers or labels

The global contracted layer should then collect those local contractions into
one global contraction matrix over the shared primitive layer.

## First contraction rule

The first contraction rule should be intentionally simple.

A good first choice is:

- local Lowdin orthonormalization inside each leaf
- local position diagonalization on the orthonormalized subspace
- retain a small fixed number of local functions per leaf

This keeps the first implementation understandable while still showing the
correct architecture:

- one global common layer
- one optional local contraction layer

## What this is meant to prove

This slice is not meant to solve the full nested problem. It is meant to prove
that:

- the corrected global-first design fits the current primitive and hierarchy
  stack
- the uncontracted global mapped basis is usable directly
- local contraction can be added as a refinement layer without introducing
  separate per-leaf maps
- this direction is closer to the historical nested idea than
  `LeafLocalPGDG1D`, which is better viewed as a useful prototype

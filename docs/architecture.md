# Package architecture and current direction

This note gives the shortest useful map of what the package currently is.

It is not a replacement for the README or the example guide. Its purpose is to answer a different question:

**how do the main ideas in GaussletBases fit together?**

## 1. Broad foundation: ordinary gausslets

The broadest foundation of the package is the ordinary one-dimensional gausslet construction.

At that level, the main ideas are:

- a gausslet is built from an explicit Gaussian primitive layer
- mappings can redistribute resolution in physical space
- basis functions can be understood through their underlying primitive expansions
- matrices can be built either directly at the basis level or by building them on the primitive layer and contracting upward

This is the most general part of the present package. The radial and hierarchy lines both sit on top of it.

## 2. Mature public-facing path: radial work

The most mature public-facing path in the package is the radial line.

That is why the README and quickstart emphasize:

- radial basis construction
- explicit radial quadrature
- basis diagnostics
- radial one-body operators
- hydrogen as the first scientific validation step

This is not because radial gausslets are the whole point of the package. It is because the radial path is the clearest mature scientific workflow today.

## 3. Shared primitive layer

The key advanced structural idea is the shared primitive layer.

In practical terms, the package can expose:

- the common primitive set behind a basis
- the coefficient matrix that contracts primitives into basis functions
- primitive-level matrices
- basis-level matrices obtained by contraction

This is the bridge between “basis functions as final objects” and “basis functions as controlled contractions of a common Gaussian layer.”

That bridge matters scientifically because it makes matrix construction, locality, and alternative contraction strategies explicit instead of hidden.

## 4. Representation, partitions, and hierarchy

On top of the primitive layer, the package now has:

- a compact in-memory basis representation
- partitions of basis functions into interval boxes
- a simple parent-child hierarchy of those boxes

These layers are not yet full nested-gausslet machinery. They are the organizational scaffolding needed to study locality and local contraction cleanly.

## 5. Prototype line versus corrected direction

There are now two distinct advanced 1D lines in the repository.

### Prototype line

The prototype line is built around:

- `LeafLocalPGDG1D`
- leaf-local Gaussian generation
- leaf-local Gaussian augmentation

This line is useful because it proved that hierarchy-driven local generation, provenance tracking, and local enrichment can all be implemented cleanly.

But it is best understood as a prototype line, not the main conceptual direction.

### Corrected direction

The corrected direction is:

1. one global mapped primitive layer
2. optional local contraction on boxes or leaves of that common layer

This is the line represented most clearly by:

- `GlobalMappedPrimitiveLayer1D`
- `LeafBoxContractionLayer1D`
- `examples/13_global_leaf_contraction.jl`

This direction is closer to the historical nested idea because it keeps one shared global primitive basis and performs local contraction inside that common layer.

## 6. Why the global mapped layer matters

An important conceptual point is that the globally mapped uncontracted primitive/product basis is not just scaffolding.

It is useful in two ways:

- as the common substrate for later local contraction
- as a simple basis or intermediate object that may already be useful directly

So the package should not be read as “everything is waiting for contraction before it becomes meaningful.” The global mapped layer is already a scientifically meaningful object.

## 7. What is mature, what is advanced, what is still open

### Mature

- ordinary 1D gausslet construction
- radial basis construction and diagnostics
- radial one-body operator workflow
- hydrogen validation example

### Advanced but structurally real

- primitive-layer matrix construction
- basis contraction through a visible primitive layer
- basis representation
- partitions and hierarchy
- global mapped layer plus leaf contraction

### Still open

- the best local contraction criterion for a genuinely useful nested step
- parent-shell contraction beyond leaves
- geometry-aware grouping for atoms and molecules
- named Gaussian chemistry basis support
- higher-dimensional nested workflows

## 8. How to read the repository now

If you are new:

1. read the README
2. read `docs/first_radial_workflow.md`
3. run the first four examples
4. read `docs/intermediate_primitive_layer.md`
5. then read `docs/example_guide.md`

If your interest is specifically the current nested/contraction direction:

1. read `docs/intermediate_primitive_layer.md`
2. read `docs/global_map_local_contraction.md`
3. read `docs/global_mapped_leaf_contraction_1d.md`
4. run `examples/13_global_leaf_contraction.jl`
5. treat `11` and `12` as prototype side studies

## 9. Bottom line

The package should be understood this way:

- ordinary gausslets are the broad foundation
- radial gausslets are the current mature public-facing workflow
- primitive layers, contraction, partitions, and hierarchy are the structural bridge to future nested work
- the corrected current research direction is global mapped layer plus local contraction, not separate local basis worlds in each leaf

# The shared primitive layer

This is the first advanced note in the repository.

You do **not** need this page in order to run the radial examples or understand the hydrogen example. But if you want to understand how the package is built, or if you are interested in contraction and hierarchy, this is the right next note.

The basic idea is simple:

**the final basis is built from a shared Gaussian-type primitive layer, and the package makes that layer visible**

That matters because it turns the basis from a black box into something you can inspect, manipulate, and contract in a controlled way.

## 1. Why this layer exists

A gausslet basis function is not defined out of thin air. It is built from simpler Gaussian-type pieces.

For one function, that exact expansion is available through:

```julia
stencil(f)
```

For a whole basis, the package can expose:

- the common primitive layer behind the basis
- the coefficient matrix that combines those primitives into the final basis functions

That is the key structural idea behind the advanced side of GaussletBases.

## 2. The primitive set

The structured object for the common primitive layer is:

- `PrimitiveSet1D`

For a basis `b`, the preferred high-level call is:

```julia
P = primitive_set(b)
```

If you want the raw list of primitive function objects, you can also inspect:

```julia
primitives(b)
```

You do not need this for everyday radial work. It becomes useful when you want to understand how matrices are built, or when you want to study contraction and localization directly.

## 3. From primitives to basis functions

For a basis `b`, the coefficient matrix that contracts primitives into basis functions is:

```julia
C = stencil_matrix(b)
```

Mathematically, that means:

- each column of `C` gives one basis function in the shared primitive layer
- the rows follow the primitive ordering
- the columns follow the basis-function ordering

That is why the package can move cleanly between primitive-level and basis-level matrices.

## 4. Why contraction matters

Suppose you form a matrix at the primitive level, such as an overlap or kinetic matrix. If `C` is the contraction matrix, then the basis-level matrix is just the contracted version of that primitive matrix.

In code:

```julia
P = primitive_set(b)
C = stencil_matrix(b)

Smu = overlap_matrix(P)
Sb = contract_primitive_matrix(b, Smu)
```

This is one of the most important advanced ideas in the package:

**the basis is not a black box; it is a contraction built on a visible primitive layer**

That makes it easier to reason about representation, locality, and alternative contraction strategies.

## 5. Basis representations

The package also has a compact advanced object for bundling this information:

- `BasisRepresentation1D`

built by:

```julia
rep = basis_representation(b)
```

You can think of this as a convenient in-memory bundle containing:

- basis metadata
- the primitive layer
- the contraction matrix
- selected primitive-level matrices
- the corresponding contracted basis-level matrices

This is mainly for advanced users, method developers, and downstream experiments. It is not the first thing a new user needs.

## 6. Why this matters scientifically

This primitive-layer view is useful for at least three reasons.

### First: it explains where the basis really comes from

A basis function is not magical. It is a structured combination of simpler pieces.

### Second: it makes matrix construction transparent

You can see where a matrix is built and how it is carried upward by contraction.

### Third: it opens the door to local contraction and hierarchy

Once a common primitive layer is explicit, you can ask meaningful questions about:

- partitioning basis functions in space
- local retained spaces
- local versus global contraction
- hierarchy built on a common substrate

That is exactly where the newer 1D research line in the repository begins.

## 7. Partitions and hierarchy

The package includes a small 1D partition and hierarchy layer:

- `BasisPartition1D`
- `HierarchicalBasisPartition1D`

These objects are there to help study locality and box structure in a way that is explicit and easy to inspect.

The important point is that these are **advanced** tools, not required for ordinary radial calculations.

## 8. The current experimental direction

The current research direction in the repository is not “every leaf gets its own completely separate basis world.”

The cleaner present direction is:

1. one global mapped primitive layer
2. local contraction inside boxes or leaves of that common layer

That is why the package now includes objects such as:

- `GlobalMappedPrimitiveLayer1D`
- `LeafBoxContractionLayer1D`

This is a more faithful and more scientifically useful direction than treating each leaf as an unrelated local basis construction problem.

## 9. Where the prototype line fits

The package also contains a prototype line built around:

- `LeafLocalPGDG1D`

This is useful because it shows that hierarchy-driven local generation can work cleanly in code. But it should still be treated as a prototype, not as the main conceptual story of the package.

For most users, the right way to think about the current state is:

- radial calculations are the mature public-facing path
- primitive-layer contraction is the first advanced architectural idea
- leaf-local generation is a useful prototype line
- global mapped layer plus local contraction is the more faithful current research direction

## 10. What to read next

After this page, the most useful next reads are:

- [`example_guide.md`](example_guide.md) for the order of the advanced examples
- [`global_map_local_contraction.md`](global_map_local_contraction.md) for the correction away from per-leaf maps
- [`global_mapped_leaf_contraction_1d.md`](global_mapped_leaf_contraction_1d.md) for the current experimental target

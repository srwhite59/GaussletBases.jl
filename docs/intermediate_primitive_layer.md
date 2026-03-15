# Intermediate primitive layer note

This note defines the small Stage-2 slice added after the first radial release.

The goal is not to port nested gausslets yet. The goal is to validate the
shared intermediate layer that ordinary, radial, and later nested work can all
use.

## Minimal new objects

The first new public intermediate object is:

- `PrimitiveSet1D`

This is an explicit ordered set of lowest-level primitives. Those primitives may
be plain objects such as `Gaussian(...)`, or explicit mapped objects such as
`Distorted(Gaussian(...), mapping)`.

The first small export-oriented basis object is:

- `BasisMetadata1D`

This stores the basis kind, family name, mapping, centers, integral weights,
the shared `PrimitiveSet1D`, and the coefficient matrix that combines the
primitive layer into the final basis.

That is enough to support downstream export or reconstruction of matrix
workflows without committing to a much larger future schema.

## Minimal public API in this slice

This Stage-2 proof slice adds:

- `PrimitiveSet1D(primitives; name=nothing, labels=nothing)`
- `primitive_set(basis)`
- `overlap_matrix(set::PrimitiveSet1D)`
- `kinetic_matrix(set::PrimitiveSet1D)`
- `basis_metadata(basis)`

The matrix builders return ordinary dense matrices, just like the current radial
operator builders.

## Analytic and numerical backends

The public interface is intentionally *not* organized around backend choice.

The public call is simply:

```julia
S = overlap_matrix(set)
T = kinetic_matrix(set)
```

Under the hood, the implementation chooses between:

- an analytic path, when the primitive content supports it
- a numerical path, when explicit mapped or otherwise unsupported primitive
  content is present

For this first slice:

- plain full-line `Gaussian` primitive sets use an analytic path
- distorted sets such as `Distorted(Gaussian(...), mapping)` use numerical
  quadrature

That keeps the ordinary-gausslet analytic route alive while still supporting
explicitly distorted primitive sets through the same public calls.

## Proof-of-concept example

The first proof target is deliberately small:

1. start from the Gaussian primitive layer of an ordinary gausslet stencil
2. build its primitive overlap and kinetic matrices through the public API
3. build a distorted version of that primitive set
4. show that the same public matrix builders work there too

The comparison test for this slice checks that, for a plain Gaussian primitive
set, the analytic and numerical overlap builders agree.

## Feeding the primitive layer upward into a basis

For an existing basis, the shared primitive layer is provided by:

- `primitive_set(basis)`

The contraction map from that primitive layer to the final basis functions is
already the basis stencil matrix:

- `stencil_matrix(basis)`

So the basis-level contraction story is:

```julia
P = primitive_set(basis)
C = stencil_matrix(basis)

Smu = overlap_matrix(P)
Tmu = kinetic_matrix(P)

Sb = contract_primitive_matrix(basis, Smu)
Tb = contract_primitive_matrix(basis, Tmu)
```

That is the minimal basis-wide consumer for this layer. It shows that the same
primitive matrix object can feed an existing basis workflow without introducing
new shell or nested abstractions.

## Metadata versus live computation

`BasisMetadata1D` stores the structural data needed by downstream consumers:

- basis kind and family
- mapping
- centers and integral weights
- the shared primitive set
- the basis coefficient matrix

What it does *not* store are derived operator matrices such as overlap or
kinetic matrices. Those remain live computations, because they may come from
either analytic or numerical backends and may evolve as those backends improve.

## What this does not try to settle

This note does not define the final nested-gausslet API.

It also does not define a final export format for every future workflow. The
purpose of this slice is only to prove that:

- the package can represent an explicit primitive layer cleanly
- matrix builders can sit above that layer
- analytic and numerical evaluation can share one public interface

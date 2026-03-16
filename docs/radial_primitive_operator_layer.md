# Radial primitive operators on a shared quadrature grid

This note describes the next architectural step for the radial line.

The central idea is simple:

**the radial basis already has a common primitive layer, so radial one-body
operators should be allowed to live there too**

The goal of this step is not to derive analytic formulas for mapped radial
primitive integrals. The goal is to apply the same primitive-layer and
contraction architecture that now exists elsewhere in the repository to the
radial basis path.

## What the radial basis already has

For a radial basis `rb`, the package already exposes:

- a common primitive layer through `primitive_set(rb)`
- the contraction matrix from that primitive layer to the working basis through
  `stencil_matrix(rb)`

That means the structural ingredients are already present.

What has been missing is the corresponding primitive-space radial operator
layer.

## The right next move

The right next move is:

1. take the primitive layer behind the radial basis
2. evaluate primitive-space one-body matrices numerically on an explicit radial
   quadrature grid
3. contract those primitive matrices upward to the working radial basis

That gives the radial line the same clean story that now exists on the
ordinary-gausslet side:

- one common primitive layer
- one visible contraction map
- operators that can be understood first in primitive space and then in basis
  space

## Why not analytic mapped radial primitive formulas right now?

There are at least three reasons not to go in that direction yet.

### 1. The mapped radial primitives are already complicated enough

Once the primitives are explicitly mapped, the formulas become messy quickly.
That is even more true when injected `XGaussian`s are present near the origin.

### 2. The package already has a good explicit quadrature object

The radial line already uses an explicit `RadialQuadratureGrid`.

That grid is the natural numerical object on which to evaluate primitive-space
radial operators. Reusing it keeps the operator story consistent with the rest
of the radial package.

### 3. The architectural gain matters more right now than closed forms

At this stage, the main value is that the radial path will sit on the same
visible primitive/contraction architecture as the rest of the repository.

That is more important than winning a difficult analytic-integral campaign
before the larger radial and angular structure is settled.

## What operators belong in the first slice

The first slice should be the radial one-body operators only:

- overlap
- kinetic
- nuclear
- centrifugal

These are exactly the operators that already define the present one-electron
radial workflow.

The first implementation should therefore allow a user to do:

```julia
P = primitive_set(rb)

Smu = overlap_matrix(P, grid)
Tmu = kinetic_matrix(P, grid)
Vmu = nuclear_matrix(P, grid; Z = Z)
Cmu = centrifugal_matrix(P, grid; l = l)
```

and then contract upward with the existing basis contraction matrix.

## What this should prove

For a fixed radial basis and quadrature grid, the contracted primitive-space
matrices should reproduce the existing direct basis-space radial operators to
numerical precision.

That comparison is essential. If it works, then the radial line is no longer a
separate architectural story. It becomes another instance of the same general
primitive-layer pattern.

## Why this prepares the way for Ylm

Ylm and angular machinery should sit on top of a clean radial substrate.

The right substrate is:

- a common radial primitive layer
- quadrature-backed primitive-space one-body operators
- contraction to the working radial basis

That gives a clearer place to attach later angular structure. It also keeps the
radial line transparent: one can still see the shared primitive layer under the
basis functions and under the radial operators.

## Scope of this step

This step is intentionally narrow.

It does **not** include:

- analytic mapped radial primitive formulas
- Ylm
- two-body solve machinery
- a broader atom package layer
- solver infrastructure

It is simply the architectural move that makes the radial line sit cleanly on
the same primitive/contraction foundation as the rest of the repository.

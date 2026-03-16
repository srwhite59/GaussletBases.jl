# Current Status

GaussletBases now has one clearly established user-facing path and one newer advanced/research path.

## Established public-facing path: radial calculations

The most mature use of the package today is atom-centered radial work.

This includes:

- ordinary 1D, half-line, and radial gausslet bases
- explicit radial quadrature, separate from the basis
- radial diagnostics
- radial overlap, kinetic, nuclear, and centrifugal matrices
- a first explicit one-electron `(l,m)` atomic layer built on those radial matrices
- the current two-index IDA-style radial multipole matrices
- `RadialAtomicOperators`

For a new user, this is the clearest place to begin.

## First atomic angular step

The package now also has a first explicit atomic angular layer built from channels labeled by `(l,m)`.

This is still a narrow one-electron path, but it is an important step because it makes the hydrogen problem look like an actual radial-plus-angular atomic basis rather than only a radial test problem.

The present objects there are:

- `YlmChannel`
- `YlmChannelSet`
- `AtomicOneBodyOperators`

This is the right layer to add before the later interacting atomic path.

## Advanced but still settling: primitive layers and contraction

The package also has a real one-dimensional structure for studying:

- primitive Gaussian layers
- contraction from primitives to basis functions
- basis representations
- partitions and simple hierarchy
- a globally mapped common layer plus local contraction

This line is already useful for method development, but it is not yet as settled or beginner-friendly as the radial path.

## Prototype line

`LeafLocalPGDG1D` and its augmentation should still be viewed as useful prototypes.

They are scientifically helpful because they show that hierarchy-driven local generation can work cleanly. But they are not yet the clearest public conceptual path for the package.

The more canonical current research direction is:

- one global mapped primitive layer
- local contraction inside boxes or shells of that shared layer

## What is not in the package yet

The current repository does not yet provide:

- exact non-diagonal electron-electron operators
- the fuller interacting spherical-angular atomic layer beyond the present one-electron `(l,m)` path
- full nested 2D or 3D workflows
- named chemistry basis-set libraries
- solver layers such as HF, DMRG, or related workflows
- Python / Fortran interop layers

That is intentional. The package is still clarifying its scientific structure before it grows into those directions.

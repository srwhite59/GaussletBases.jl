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

## First static interacting atomic ingredients

The package now also has a narrow static interacting atomic layer built in the current IDA style.

This layer bundles:

- the one-body atomic blocks
- the radial multipole tables
- the angular Gaunt and M-summed kernel data
- explicit orbital indexing

It is an assembly layer only. It does not yet solve the many-electron problem.

## First minimal atomic mean-field line

The repository now also has the first small Hartree/Fock-style layer built on
top of the current atomic IDA ingredients.

That line now includes:

- the direct/Hartree term
- the exchange term
- the algebraic Fock-style combination
- a minimal UHF fixed-point kernel

This is already enough to run small He-like mean-field tests in the present
atomic IDA model.

But it should not be overdescribed. This is **not** yet a broad general HF
workflow. It is a narrow current-model mean-field layer built on:

- the radial basis
- the `(l,m)` atomic channel layer
- the current static IDA ingredients

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
- a broad general HF/RHF/UHF workflow beyond the present minimal atomic IDA kernel
- a larger interacting atomic workflow beyond the present tiny exact and minimal UHF layers
- full nested 2D or 3D workflows
- named chemistry basis-set libraries
- larger solver layers such as DMRG or related workflows
- Python / Fortran interop layers

That is intentional. The package is still clarifying its scientific structure before it grows into those directions.

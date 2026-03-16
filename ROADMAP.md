# Roadmap

This roadmap is meant as a short scientific guide to the next questions for GaussletBases.

It is not a schedule. It is a statement of which directions now look most valuable, especially if the package is meant to become useful to a slice of the electronic-structure community rather than only as an internal research notebook.

## Where the package now stands

Two directions are now real in the code:

1. a mature radial basis / quadrature / operator line
2. a first explicit one-electron atomic `(l,m)` layer on top of that radial line
3. a newer one-dimensional primitive-layer / contraction / hierarchy line

Those directions are both scientifically interesting, but they are not equally mature or equally valuable to outside users.

## Highest-value next additions for outside users

If the package is meant to become useful to method developers in atomic and related electronic-structure work, the next highest-value additions are likely these.

### 1. An exact non-diagonal radial electron-electron layer

The current package has the two-index IDA-style radial multipole matrices.

A natural next scientific step is the exact non-diagonal radial electron-electron object. That would make it possible to study more carefully where the present radial approximation is strong, where it is weak, and what the true cost/accuracy tradeoffs look like.

### 2. The first actual He / IDA-style solve on top of the present static atomic ingredients

The package now has:

- the explicit one-electron angular `(l,m)` layer
- the first static interacting IDA ingredients

The next atomic question is how to turn those explicit ingredients into a first useful interacting atomic calculation:

- radial one-body operators
- explicit `(l,m)` channels
- radial multipole data
- later Gaunt or related angular-coupling factors

That is the natural path toward helium and related atom-centered work.

### 3. Export and interoperability helpers

Once the radial and atomic pieces are clearer, export helpers become much more valuable.

The likely first targets are simple, explicit formats that another code can consume without needing to reimplement the basis construction.

## Important research questions inside the contraction line

The contraction/hierarchy line remains scientifically important, but it is currently more of a research track than a public-facing track.

The main questions there are:

### 1. What should define the first genuinely useful local retained space?

The immediate next question is still:

- what local contraction criterion should be used inside a box or shell?

The present simple retained spaces are good enough to demonstrate the architecture, but not yet enough to settle the scientific choice.

### 2. How should leaf contraction grow into parent-shell contraction?

The present corrected path is leaf-only.

If that line continues, the next structural extension is likely:

- parent-shell contraction on top of the existing leaf structure

### 3. When should geometry-aware grouping enter?

Eventually, one-dimensional interval boxes will no longer be enough.

The question is when it becomes scientifically worthwhile to move from:

- simple 1D locality and hierarchy

to:

- atom- or geometry-aware grouping

## What the roadmap is not promising yet

This roadmap is not a commitment to:

- a complete solver workflow
- immediate molecule-scale infrastructure
- a permanent exchange format
- a fully settled nested or PGDG public surface
- immediate Python/Fortran bindings

Those are possible later directions, but they are not the present center of gravity.

## Practical interpretation

If the goal is:

### Better public usefulness soon

then the package should prioritize:

1. exact radial electron-electron structure
2. the first interacting He / IDA-style atomic calculation on top of the present static ingredients
3. export/interoperability

If the goal is:

### Deeper basis/contraction research first

then the package should prioritize:

1. better local retained spaces
2. parent-shell contraction
3. geometry-aware grouping

Right now, the code contains seeds of both futures. The next major choice is which one should lead the public story.

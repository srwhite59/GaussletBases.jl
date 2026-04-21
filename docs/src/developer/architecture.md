# Architecture and current direction

This page gives the shortest useful map of how the main ideas in
`GaussletBases.jl` fit together.

If you are trying to use the package rather than understand its internal
structure, start in the [Manual](../manual/index.md) instead.

## Broad foundation

The broadest foundation of the package is the ordinary one-dimensional
gausslet construction:

- explicit Gaussian primitive layers
- coordinate mappings
- basis functions built from visible primitive expansions
- matrix construction either directly at the basis level or on the primitive
  layer and contracted upward

## Mature public-facing path

The most mature public-facing path in the package is still the radial line.
That is why the main onboarding path emphasizes:

- radial basis construction
- explicit radial quadrature
- basis diagnostics
- radial one-body operators
- hydrogen as the first scientific validation step

## Atomic line

On top of the radial substrate, the package now has:

- an explicit one-electron `(l,m)` layer
- a static interacting atomic IDA layer
- direct / exchange / Fock helpers
- a minimal UHF kernel
- dense and sliced export for downstream solver consumers

That is a real small atomic line, but not yet a broad atomic workflow package.

## Ordinary line

The ordinary Cartesian branch is the experimental line for mapped and hybrid
ordinary gausslets. Its current interpretation is:

- Coulomb-expansion first
- mild and hybrid regimes as the practical target
- numerical validation route plus experimental PGDG-style analytic backend

Within that branch, the pure bond-aligned routes now normalize through one
internal Cartesian build context before backend gating and operator assembly.
That seam now covers:

- pure bond-aligned direct-product ordinary/QW routes
- pure prebuilt nested fixed-block bond-aligned routes
- bond-aligned diatomic molecular supplement direct-product routes
- bond-aligned diatomic molecular supplement nested fixed-block routes
- atomic supplement direct-product routes
- atomic supplement nested fixed-block routes

Above that operator-side seam, the source-backed nested front-door layer now
normalizes through one internal nested-source context before route-specific
geometry source assembly and path packaging. That seam now covers:

- public bond-aligned diatomic nested source / fixed-block / diagnostics front
  doors
- experimental bond-aligned homonuclear chain nested wrappers
- experimental axis-aligned homonuclear square-lattice nested wrappers

Across those source-backed nested lines, the common glass-box subset is now:
`fixed_dimension`, `contract_audit`, `shared_shell_dimensions`,
`shared_shell_provenance`, and `leaf_count`, while family-specific topology and
split-policy diagnostics remain separate.

Experimental supplement-bearing routes still keep their separate orchestration
paths on purpose.

## Primitive, contraction, and hierarchy work

The package also has an advanced structural line built around:

- visible primitive layers
- explicit contraction matrices
- basis representations
- partitions and hierarchy
- global mapped primitive layers plus local contraction

That is important for the package’s long-term research direction, even though
it is not the first page a new user should start from.

## Current bottom line

The shortest package-shape summary is:

- ordinary gausslets are the broad foundation
- radial gausslets are the mature current workflow
- the atomic line sits on top of the radial substrate
- the ordinary mapped/hybrid line is promising but still experimental
- primitive layers and contraction are the structural bridge to later work

# Radial Extent Policy Note

This note records the intended separation between the three radial extents now
used in the repo.

## The three extents

- `rmax`
  - public, user-facing
  - the center of the last retained radial gausslet
  - this is the main scientific extent a user chooses
- `build_umax`
  - internal mapped-coordinate setup extent
  - used while building, orthonormalizing, and localizing the radial basis
  - allowed to include build padding
- `quadrature_umax`
  - internal mapped-coordinate active operator extent
  - used to choose the default physical quadrature extent for operators and diagnostics
  - should come from retained-basis support, not broader build padding

## Intended policy

The library should not teach users to tune internal quadrature extent as part
of the normal workflow. The front-door choice is `rmax`, together with the
mapping and other basis controls. The library then owns the internal setup and
quadrature extents.

That means:

- `radial_quadrature(rb)` is the normal public path
- `quadrature_rmax` remains only as an expert compatibility override
- large physical quadrature extents are not automatically bugs if they are
  genuinely implied by retained support under the mapping

## Current implementation rule

For the current radial path:

- `build_umax` is stored from the internal construction grid extent
- `quadrature_umax` is derived from retained radial centers together with the
  runtime family support width
- the default physical quadrature cutoff is `xofu(mapping(rb), quadrature_umax)`

So the active quadrature extent is no longer determined by broader primitive or
setup padding. These two internal extents answer different questions, so there
is no guarantee that `build_umax` must be numerically larger or smaller than
`quadrature_umax` on every map/family combination.

# Current ordinary branch

This page is the shortest current user-facing status read for the newer
ordinary / Cartesian line.

## What the ordinary branch is today

The current ordinary story is:

- Coulomb-expansion first, not 3D grid first
- one-dimensional mapped ordinary bases on each Cartesian axis
- exact Cartesian workflow primitives:
  - `cross_overlap`
  - `basis_projector`
  - `transfer_orbitals`
- an explicit backend split between validation and experimental PGDG-style
  analytic construction
- a separate paper-faithful Qiu-White residual-Gaussian route
- real one-center nested Cartesian fixed-block / consumer support
- real bond-aligned diatomic nested source / fixed-block / diagnostics /
  geometry payload support
- the old COMX-cleaned 1D hybrid route quarantined as legacy/internal code,
  not current workflow

This line is still experimental, but it is now a real second workflow surface
in the package rather than a hidden side note.

## Start here

Within the rendered docs, the best supporting pages are:

- [Example guide](../howto/example_guide.md)
- [Atomic and ordinary workflows reference](../reference/atomic_and_ordinary.md)
- [Qiu-White residual-Gaussian route](../algorithms/qiu_white_residual_gaussian_route.md)
- [Cartesian nested atomic nonrecursive route](../algorithms/cartesian_nested_atomic_nonrecursive_route.md)
- [Cartesian nested diatomic box policy](../algorithms/cartesian_nested_diatomic_box_policy.md)

The radial tutorial remains useful context because the package’s operator and
quadrature story was established there first, but it is no longer the whole
ordinary-branch story.

## Documentation authority for this branch

For the current ordinary / Cartesian story, use these pages in this order:

- this page for branch status and scope
- [Reference](../reference/index.md) for exported API surface
- [Algorithms](../algorithms/index.md) for path-specific construction recipes
- [Developer Notes](../developer/index.md) for lower-priority architecture and
  supporting material

Older flat `docs/*.md` notes remain available as supporting or historical
material, but they should not override this page as the current ordinary branch
summary.

## Current interpretation

The wording discipline for the ordinary line remains:

- `:numerical_reference` is the validation route
- the PGDG-style analytic route is good enough on the mapped ordinary backbone
- on that PGDG lane, distorted 1D primitives are replaced by plain Gaussian
  proxies before matrix assembly and numerical primitive quadrature is outside
  the intended production contract
- the public mapped ordinary backbone is the PGDG-capable surface today:
  `mapped_ordinary_one_body_operators` and the mapped Cartesian IDA path may
  use `:pgdg_experimental` or `:pgdg_localized_experimental`
- current public Qiu-White and nested Qiu-White routes remain
  `gausslet_backend = :numerical_reference` unless a route is explicitly
  documented otherwise; that includes the current diatomic and experimental
  chain/square nested QW routes
- exact overlap / projector / transfer are now first-class workflow primitives
- Gaussian-supplement comparisons should go through the separate paper-faithful
  3D Qiu-White route, not the old 1D COMX-cleaned surrogate path
- the public one-center White-Lindsey-style mapping language is now the
  `d`-driven helper documented in [Bases and mappings](../reference/bases_and_mappings.md),
  not the older fixed-`a` / `count -> s` historical note language
- hard pure mapped small-`c` cases remain stress tests
- `AsinhMapping` is the current working map, not final truth

The ordinary branch also includes narrow experimental producer-side chain and
square-lattice lines. Those are real code, but they should still be read as
experimental rather than as settled broad workflows.

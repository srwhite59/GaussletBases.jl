# Recommended atomic setup

This note is for users who want a sensible starting point for atom-centered radial calculations with GaussletBases.

The package exposes several method-level choices on purpose. That is useful for research, but it also means a new user needs one clear recipe to start from. This page gives that recipe.

## Short version

For an atom of nuclear charge `Z`, start here:

* `s = 0.2`
* `c = s / (2Z)`
* `tail_spacing = 10.0`
* `tails = 6`
* `odd_even_kmax = 6`
* `rmax = 30.0` bohr for a first-row atom
* `reference_spacing = 1.0`

In code:

```julia
using GaussletBases

Z = 2.0
s = 0.2

map = AsinhMapping(c = s / (2Z), s = s)

spec = RadialBasisSpec(:G10;
    rmax = 30.0,
    mapping = map,
    reference_spacing = 1.0,
    tails = 6,
    odd_even_kmax = 6,
)

rb = build_basis(spec)
grid = radial_quadrature(rb)
```

This is the recommended package starting point, not a hard-wired code default.

## What the main parameters mean

### `s`

`s` controls the overall radial spacing through the coordinate mapping.

* smaller `s` means finer radial resolution and usually more basis functions
* larger `s` means a coarser basis and fewer functions

For the package documentation, use `s = 0.2` as the standard starting value.

This is not a claim that `0.2` is uniquely best. In practice, larger values such as `0.3` to `0.5` can still give useful results, especially for exploratory work or smaller basis sizes. But if you are unsure, start at `0.2`.

### `c = s / (2Z)`

`c` sets the near-origin scale in a way that depends naturally on the nuclear charge.

The simple rule

```julia
c = s / (2Z)
```

puts more resolution near the nucleus as `Z` grows.

That is the first rule to try unless you have a specific reason to tune the mapping.

### `tail_spacing = 10.0`

The `AsinhMapping` used here includes a built-in linear tail term. Its job is to keep the spacing at large radius from growing without bound.

In practice, the default

```julia
tail_spacing = 10.0
```

is the standard choice in this package documentation.

### `tails = 6`

This controls the extra half-line completion functions used near the boundary.

For serious radial work, use:

```julia
tails = 6
```

You may see smaller values in tiny demos or smoke tests, but `6` is the normal recommendation.

### `odd_even_kmax = 6`

This controls the size of the small near-origin even block used in the radial construction.

Again, for serious work, use:

```julia
odd_even_kmax = 6
```

Smaller values are useful only for very small demos.

### `rmax`

`rmax` is the center of the last retained radial gausslet, expressed in bohr.

A good first guess is:

* `rmax = 30.0` for first-row atoms
* smaller values for compact test problems like hydrogen
* larger values for diffuse states or if you want safer retained-basis coverage

Unlike the internal setup and quadrature extents, `rmax` is a genuinely
physical modeling choice. It depends on the system you want to describe.

## What about x-gaussians?

The package allows extra near-origin functions through the `xgaussians` keyword.

These are useful because the radial construction has to satisfy the boundary condition at `r = 0`, and that can weaken some of the moment properties that make the diagonal interaction approximation work well. Adding one or two x-gaussians can improve that near-origin behavior.

### Simple recommendation

Use these rules of thumb:

* for a first one-electron test such as hydrogen:

  * start with `xgaussians = XGaussian[]`
* for more serious work where the IDA-style radial multipole matrices matter:

  * try one or two x-gaussians if the diagnostics suggest the near-origin behavior needs help

A simple one-function choice for experiments is:

```julia
xgaussians = [XGaussian(alpha = 0.2)]
```

A more serious two-function starting point, and now the standard paper-parity
preset, is:

```julia
xgaussians = [XGaussian(alpha = 0.0936), XGaussian(alpha = 0.0236)]
```

These are the published optimized widths for the `K = 6` radial-gausslet paper
construction. If you override them, check the basis diagnostics rather than
assuming another choice is automatically better.

## How to tell whether the basis is behaving well

After building a basis and quadrature grid, inspect the diagnostics:

```julia
diag = basis_diagnostics(rb)

diag.overlap_error
diag.D
```

Here `basis_diagnostics(rb)` chooses its own conservative integration grid internally. That makes it a good first check before you start tuning explicit quadrature settings.

Useful quick checks are:

* `diag.overlap_error` should be small
* `diag.D` should not be unusually large for the basis you are using

Here `D` is the aggregate mismatch between the nominal basis centers and the first-moment centers computed on the quadrature grid. Near the origin, that mismatch is one of the main signs that you may want to refine the setup, use a finer quadrature grid, or try x-gaussians.

## Construction grid versus quadrature grid

There are two separate numerical ideas here, and it helps to keep them distinct.

### 1. Basis construction controls

The internal construction grid is controlled through `build_basis`, for example:

```julia
rb = build_basis(spec; grid_h = 0.01, refine_grid_h = false)
```

Most users should leave this alone at first and use the built-in refinement logic.

### 2. Public quadrature grid

The public radial quadrature grid is created separately:

```julia
grid = radial_quadrature(rb)
```

This grid is the one used for one-body operators, moment centers, diagnostics, and the current radial multipole matrices.

For most users, this no-keyword call is the right place to start. The package
chooses the internal quadrature extent automatically from the retained radial
basis support and uses the default `accuracy = :high` quadrature profile
internally. `accuracy = :medium` is the cheaper exploratory option, and
`accuracy = :veryhigh` pushes the quadrature harder when you want additional
safety. The keywords `quadrature_rmax` and `refine` are expert overrides
underneath those profiles and are not part of the normal front-door workflow.

The basis and the quadrature grid are intentionally separate. That is one of the main ideas behind the package.

## Package recommendation versus paper settings

If you are trying to reproduce exact numbers from a manuscript or an older internal code, check the corresponding source carefully.

The package documentation uses:

* `s = 0.2`
* `c = s / (2Z)`

as the standard beginner-friendly recommendation.

Some paper calculations used a somewhat smaller `s`. So if your goal is reproduction of published numbers rather than a good package starting point, use the settings stated in the paper you are matching.

## A practical checklist

If you are starting a new atom-centered calculation, this is a good checklist:

1. start with `s = 0.2`
2. set `c = s / (2Z)`
3. use `tail_spacing = 10.0`
4. use `tails = 6`
5. use `odd_even_kmax = 6`
6. start with `rmax = 30.0` bohr for a first-row atom
7. build the default quadrature grid with `radial_quadrature(rb)`
8. check `basis_diagnostics(rb)`
9. if needed, raise the quadrature accuracy, increase the basis `rmax`, refine the quadrature manually, or try one or two x-gaussians

That is enough to get most first calculations off the ground without immediately disappearing into tuning details.

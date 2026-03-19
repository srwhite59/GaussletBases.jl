# Recommended atomic setup

This page gives one practical starting recipe for atom-centered radial
calculations.

If you want the API details for the constructors used here, also see:

- [Bases and mappings](../reference/bases_and_mappings.md)
- [Operators and diagnostics](../reference/operators_and_diagnostics.md)

## Short version

For an atom of nuclear charge `Z`, start with:

- `s = 0.2`
- `c = s / (2Z)`
- `tail_spacing = 10.0`
- `tails = 6`
- `odd_even_kmax = 6`
- `rmax = 30.0` bohr for a first-row atom
- `reference_spacing = 1.0`

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
    xgaussians = XGaussian[],
)

rb = build_basis(spec)
grid = radial_quadrature(rb)
```

This is the recommended package starting point, not a hard-wired default.

## What the main knobs mean

- `s` controls the overall radial spacing.
- `c = s / (2Z)` sets the near-origin scale in a way that naturally tightens
  as `Z` grows.
- `tail_spacing = 10.0` is the current standard linear-tail choice used with
  `AsinhMapping`.
- `tails = 6` and `odd_even_kmax = 6` are the current serious-work starting
  values rather than smoke-test values.
- `reference_spacing = 1.0` is the ordinary current starting spacing on the
  reference grid.

## What about x-gaussians?

For a first one-electron check such as hydrogen, it is reasonable to start
with:

```julia
xgaussians = XGaussian[]
```

If the near-origin diagnostics need help for the current radial interaction
approximation, try one or two `xgaussians` later.

## What to check

After building the basis, inspect:

```julia
diag = basis_diagnostics(rb)
diag.overlap_error
diag.D
```

Those are the quickest signs that the basis and quadrature are behaving
sensibly.

## What this page is for

This page is the practical tuning page. If you want the workflow taught more
gently before the advanced knobs appear, start with:

- [First radial workflow](../tutorials/first_radial_workflow.md)
- [Example guide](example_guide.md)

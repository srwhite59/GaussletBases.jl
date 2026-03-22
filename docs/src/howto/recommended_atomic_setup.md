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
- `rmax = 30.0` bohr for a first-row atom

In code:

```julia
using GaussletBases

Z = 2.0
s = 0.2

map = AsinhMapping(c = s / (2Z), s = s)

spec = RadialBasisSpec(:G10;
    rmax = 30.0,
    mapping = map,
)

rb = build_basis(spec)
grid = radial_quadrature(rb)
```

This is the recommended package starting point. It deliberately leaves the
deeper construction controls at their standard defaults.

## What the main knobs mean

- `s` controls the overall radial spacing.
- `c = s / (2Z)` sets the near-origin scale in a way that naturally tightens
  as `Z` grows.
- `rmax` is the physical outer radius of the basis.

## What about x-gaussians?

There are two practical regimes here.

For a first one-electron check such as hydrogen, it is reasonable to start
with no extra near-origin supplement:

```julia
xgaussians = XGaussian[]
```

For more serious atomic / IDA work, a simple two-function starting point is:

```julia
xgaussians = [XGaussian(alpha = 0.1), XGaussian(alpha = 0.025)]
```

Those extra functions are meant to improve the near-origin radial behavior
without changing the outer mapping family. If the diagnostics remain clean
without them for the problem you care about, it is also fine to leave them out.

So the front-door interpretation is:

- the default front-door behavior is the paper-style two-function supplement
- `xgaussian_count = 0` turns that supplement off
- `xgaussian_count = 1` keeps only the first preset function
- explicit `xgaussians = [...]` is the advanced override when you really want
  to pick the alphas yourself

## Advanced construction controls

The first-read code example intentionally hides the following controls:

- `reference_spacing = 1.0`
- `tails = 6`
- `odd_even_kmax = 6`
- `tail_spacing = 10.0`
- `xgaussian_count = 2`

Why:

- `reference_spacing = 1.0` is the standard default
- `tails` is construction padding, not a scientific modeling knob
- `odd_even_kmax` is a real accuracy knob, but not something a new user should
  choose before the basic workflow is working
- `tail_spacing = 10.0` is the standard `AsinhMapping` default
- `xgaussian_count = 2` is the current recommended radial default

So the practical rule is:

- do not expose these in the beginner examples
- revisit them only when tuning accuracy or debugging basis construction

## What to check

After building the basis, inspect:

```julia
diag = basis_diagnostics(rb)
diag.overlap_error
diag.D
```

Those are the quickest signs that the basis and quadrature are behaving
sensibly. `diag.overlap_error` is the first orthonormality check. `diag.D` is
the aggregate center/moment mismatch and is often the next thing to watch when
you are deciding whether the near-origin setup needs tightening.

## What this page is for

This page is the practical tuning page. If you want the workflow taught more
gently before the advanced knobs appear, start with:

- [First radial workflow](../tutorials/first_radial_workflow.md)
- [Example guide](example_guide.md)

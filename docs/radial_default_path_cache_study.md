# Radial Default-Path Cache Study

This note records a narrow timing study of the current default public radial
workflow after the radial extent split and the interval-sampled setup-grid
rewrite.

## Scope

The target path was the current recommended atomic recipe:

- family `:G10`
- `s = 0.2`
- `c = s / (2Z)`
- `Z = 2`
- `rmax = 30.0`
- `tails = 6`
- `odd_even_kmax = 6`
- two `xgaussians`
- normal `radial_quadrature(rb)` default path

The study script is:

- `tmp/work/radial_cache_study.jl`

## Main timings

Measured in one warmed Julia session after the interval-sampled setup-grid
change:

- `build_basis(spec)` first: about `6.65 s`
- `build_basis(spec)` second: about `6.50 s`
- `basis_diagnostics(rb)`: about `0.57` to `0.72 s`
- `radial_quadrature(rb)`: about `0.44 s`
- `basis_diagnostics(rb, grid)`: about `0.12 s`
- `atomic_operators(rb, grid; Z = 2, lmax = 2)`: about `1.60` to `1.67 s`
- full front-door sequence
  - `build basis -> diagnostics -> quadrature -> atomic_operators`: about `9.43 s`

Representative basis / grid size:

- radial basis length: `39`
- `quadrature_umax ≈ 60`
- default physical quadrature cutoff: about `225.18 bohr`
- default quadrature grid length: `3841`

For historical comparison, the same build step had previously been about
`171` to `173 s` before the interval-sampled setup-grid layer landed.

## Where time is going now

The standard public path is no longer dominated by basis construction in the
old pathological way.

The basis build is still the largest single step, but it is now in the same
single-digit-seconds regime as the rest of the front-door workflow rather than
being two orders of magnitude larger than everything else.

The duplicated work in the normal user sequence

```julia
diag = basis_diagnostics(rb)
grid = radial_quadrature(rb)
```

is still real, because `basis_diagnostics(rb)` internally builds a default
quadrature grid and then `radial_quadrature(rb)` builds it again. But the cost
is modest:

- `basis_diagnostics(rb)` minus `basis_diagnostics(rb, grid)` is about
  `0.45` to `0.60 s`

So the remaining default-path startup cost is spread across several moderate
steps rather than being dominated by one catastrophic construction phase.

## Cache conclusion

Caching is no longer the first thing the radial path needs.

If caching is considered later, the narrowest sensible cache target is still a
cache of the finalized `RadialBasis` construction result, keyed by the full
radial recipe:

- family name
- mapping type and parameters
- `rmax` or `count`
- `reference_spacing`
- `tails`
- `odd_even_kmax`
- `xgaussian` alphas
- package / table version marker

The useful cached object would still be the finalized basis data, not just a
quadrature grid:

- `primitive_data`
- `coefficient_matrix`
- `reference_center_data`
- `center_data`
- `integral_weight_data`
- `build_umax`

But after the interval-sampled setup-grid rewrite, that cache is now an
optional usability enhancement rather than an urgent repair for the default
public path.

## First run versus repeated runs

The repeated `build_basis(spec)` timing is still close to the first one, so
the remaining build cost is mostly real numerical work rather than pure
first-call compilation noise.

That means:

- a lazy on-demand cache would mainly help repeated runs
- it would not materially help the first ever build of a new spec
- a shipped/prebuilt cache would only make sense if the package wanted a small
  number of exact public radial recipes to feel nearly instant on first use

## Practical outcome

The interval-sampled build rewrite materially changed the cache decision.

The default radial path is now already practical enough to use directly, and a
cache study can stay secondary until there is a clearer need for repeated-run
convenience or shipped standard-recipe artifacts.

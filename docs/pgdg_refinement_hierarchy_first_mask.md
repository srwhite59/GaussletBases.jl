# PGDG Refinement Hierarchy First Mask

This note records the first clean implementation step for the proposed
one-dimensional distorted-gausslet PGDG refinement hierarchy.

The current practical choice is now settled:

- the first stored hierarchy mask is the analytic ternary `1 -> 1/3`
  Gaussian quasi-refinement mask
- the first practical default is `rho = 1.2`

This pass adds:

- one small internal mask object
- one-step application of that mask to a uniform Gaussian coefficient line
- repeated application of the same local ternary step

It does **not** add the full hierarchy machinery yet.

## Interpretation

The intended ladder remains:

- base PGDG spacing `1/3`: no extra refinement
- refined proxy spacing `1/9`: one application of the stored ternary mask
- refined proxy spacing `1/27`: two applications
- and so on

So the new mask is the repeated local transfer step for the hierarchy.

## Exact convention used

For the centered ternary step:

- coarse spacing `H = 1`
- fine spacing `h = 1/3`
- coarse width `Sigma = rho`
- fine width `sigma = rho / 3`
- `tau^2 = 8 rho^2 / 9`

The stored coefficients are the direct analytic quasi-refinement values

```math
c_k
=
\phi_\tau(k/3)
=
\frac{3}{4\sqrt{\pi}\rho}\exp\!\left(-\frac{k^2}{16\rho^2}\right),
```

with finite truncation only.

For the first practical default:

- `rho = 1.2`
- half-window `±8` coarse units
- support radius `24` fine-grid sites

This is the same truncation convention that matched the successful
analytic-versus-fitted mask comparison.

## Current internal API

This pass adds an internal mask object and narrow helpers in
`src/ordinary_pgdg_refinement_masks.jl`:

- `_TernaryGaussianRefinementMask1D`
- `_analytic_ternary_gaussian_refinement_mask(; rho = 1.2, half_window = nothing)`
- `_default_ternary_gaussian_refinement_mask()`
- `_apply_gaussian_refinement_mask(mask, coeffs; offset = 0)`
- `_apply_gaussian_refinement_mask_repeated(mask, coeffs; levels, offset = 0)`

The helper outputs are still simple coefficient-line objects:

- integer lattice offset
- coefficient vector

That is enough for the first hierarchy plumbing.

## Small note-level demo

```julia
mask = GaussletBases._default_ternary_gaussian_refinement_mask()

one_level = GaussletBases._apply_gaussian_refinement_mask(mask, [1.0])
two_levels = GaussletBases._apply_gaussian_refinement_mask_repeated(mask, [1.0]; levels = 2)
```

Operational reading:

- `one_level` is the centered `1 -> 1/3` local refinement mask
- `two_levels` is the repeated local `1 -> 1/9` coefficient line

## What This Does Not Settle Yet

This pass deliberately stops before deeper ordinary-branch integration.

It does **not** yet decide:

- how the hierarchy is threaded through all ordinary/Qiu–White paths
- whether the coarse `1/3` PGDG layer is used mainly as a basis-realization
  layer
- whether the finer `1/9`, `1/27`, ... levels are used mainly as auxiliary
  integral-evaluation layers

Those uses are still separate future questions.

The main point here is simpler:

- the first practical analytic mask is now real repo infrastructure
- the hierarchy no longer depends on an external study script to define its
  first local refinement step

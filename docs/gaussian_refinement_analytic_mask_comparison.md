# Gaussian Refinement Analytic Mask Comparison

This note records the next narrow comparison for the proposed Gaussian proxy
refinement hierarchy.

The repo now has a local analytic-mask note:

- `docs/gaussian_refinement_analytic_mask_note.md`

The present task is to compare that explicit analytic Gaussian
quasi-refinement mask against the fitted numerical mask already generated in
the recent refinement-mask study.

The practical question is:

- can the analytic formula be used directly as the default source of the first
  stored hierarchy mask?

## Setup

The same temporary study script was extended:

- `tmp/work/gaussian_refinement_feasibility.jl`

The comparison is kept narrow:

- ternary `1 -> 1/3` refinement
- main target `rho = 1.2`
- additional checks at `rho = 4/3` and `rho = sqrt(2)`

The fitted-mask side uses the same exact continuous `L2` fit as before. The
analytic-mask side uses the explicit Gaussian quasi-refinement formula.

## Analytic coefficient formula

Using the normalized Gaussian convention from the analytic note,

```math
\phi_\sigma(x)=\frac{1}{\sqrt{2\pi}\sigma}\exp\!\left(-\frac{x^2}{2\sigma^2}\right),
```

the `1 -> 1/3` centered refinement step has:

- coarse spacing `H = 1`
- fine spacing `h = 1/3`
- coarse width `\Sigma = \rho`
- fine width `\sigma = \rho/3`
- `\tau^2 = \Sigma^2 - \sigma^2 = 8\rho^2/9`

For the unnormalized Gaussian proxy convention used in the numerical study,
the natural analytic mask is

```math
c_k^{\mathrm{analytic}}
=
\phi_\tau(k/3)
=
\frac{1}{\sqrt{2\pi}\tau}\exp\!\left(-\frac{(k/3)^2}{2\tau^2}\right)
=
\frac{3}{4\sqrt{\pi}\rho}\exp\!\left(-\frac{k^2}{16\rho^2}\right).
```

This was used directly:

- same finite window as the fitted-mask study
- no extra renormalization

That is deliberate. The point here is to compare the direct analytic formula
to the fitted mask, not to refit or post-correct the analytic one.

## Positivity and symmetry

Both the fitted and analytic masks remain:

- symmetric
- strictly positive

So there is no qualitative disagreement at the mask level.

## Coefficient-by-coefficient comparison

### `rho = 1.2`

- max absolute coefficient difference: `5.12e-12`
- relative `L1` difference: `3.83e-11`
- relative `L2` difference: `2.30e-11`

Representative half-mask comparison:

```text
offset   fitted                 analytic               difference
0        3.526184897178e-01     3.526184897173e-01    -4.99e-13
1        3.376412458658e-01     3.376412458663e-01     5.06e-13
2        2.964194747957e-01     2.964194747952e-01    -5.29e-13
3        2.385936049246e-01     2.385936049252e-01     5.68e-13
4        1.760806735342e-01     1.760806735335e-01    -6.23e-13
...
20       1.017861433115e-08     1.017349850807e-08    -5.12e-12
```

### `rho = 4/3`

- max absolute coefficient difference: `1.54e-14`
- relative `L1` difference: `1.45e-13`
- relative `L2` difference: `8.17e-14`

So here the fitted and analytic masks are essentially identical at displayed
precision.

### `rho = sqrt(2)`

- max absolute coefficient difference: `1.31e-12`
- relative `L1` difference: `1.38e-11`
- relative `L2` difference: `7.62e-12`

The masks are still extremely close, though the fitted correction is more
visible than at `4/3`.

## Practical-behavior comparison

### Residue-sum / constant-reproduction ripple

#### `rho = 1.2`

- fitted: `2.0463e-11`
- analytic: `2.1708e-11`

#### `rho = 4/3`

- fitted: `5.3889e-14`
- analytic: `5.7264e-14`

#### `rho = sqrt(2)`

- fitted: `1.7135e-13`
- analytic: `2.1150e-14`

So the constant-reproduction behavior is essentially the same order throughout.

### Single-level coarse-Gaussian refinement error

#### `rho = 1.2`

- fitted: relative `L2` error `1.50510e-11`, max error `2.12863e-11`
- analytic: relative `L2` error `1.50596e-11`, max error `2.12845e-11`

#### `rho = 4/3`

- fitted: relative `L2` error `4.01462e-14`, max error `5.67764e-14`
- analytic: relative `L2` error `4.01507e-14`, max error `5.67752e-14`

#### `rho = sqrt(2)`

- fitted: relative `L2` error `3.29320e-15`, max error `1.21435e-15`
- analytic: relative `L2` error `2.52283e-14`, max error `1.15027e-15`

So the analytic mask is practically indistinguishable at `rho = 1.2` and
`rho = 4/3`. At `rho = sqrt(2)`, it is still very good, but the fitted mask
gains about one order of magnitude in the `L2` scalar.

### Repeated-refinement behavior

The same `1/3` mask was composed for `2` and `3` levels.

#### `rho = 1.2`

- level `2`
  - fitted: `1.50795e-11`
  - analytic: `1.50805e-11`
- level `3`
  - fitted: `1.50942e-11`
  - analytic: `1.51043e-11`

#### `rho = 4/3`

- level `2`
  - fitted: `4.01662e-14`
  - analytic: `4.01676e-14`
- level `3`
  - fitted: `4.01747e-14`
  - analytic: `4.01807e-14`

#### `rho = sqrt(2)`

- level `2`
  - fitted: `1.72113e-15`
  - analytic: `1.52371e-14`
- level `3`
  - fitted: `2.61432e-15`
  - analytic: `2.86025e-14`

Again, the main practical split is:

- `rho = 1.2` and `4/3`: analytic and fitted are the same for practical use
- `rho = sqrt(2)`: analytic is still very good, but the fitted correction is
  no longer completely negligible if one is chasing the last digit

## Practical conclusion

The main conclusion is:

- yes, the explicit analytic mask is essentially the same object as the fitted
  tabulated mask for the practical first hierarchy choice `rho = 1.2`

More specifically:

- for `rho = 1.2`, the coefficient differences are around `1e-12` to `1e-11`
- the practical refinement behavior is indistinguishable at that scale
- for `rho = 4/3`, the agreement is even tighter

So for the first practical stored hierarchy mask, the analytic formula can be
used directly as the default source.

## Recommendation

My recommendation is:

- use the analytic Gaussian quasi-refinement mask as the default source of the
  first stored `1 -> 1/3` hierarchy mask

with:

- `rho ≈ 1.2` as the first practical default

Reason:

- it matches the fitted numerical mask to about `1e-11`
- it keeps the same positivity/locality story
- it avoids needing a separate fitted table-generation step for the default
  mask

A small fitted correction still looks optional rather than necessary.

The only regime where a fitted correction may still be worth keeping in mind
is:

- more aggressive high-accuracy choices such as `rho = sqrt(2)`

where the analytic mask is still very good, but no longer matches the fitted
mask to the very last displayed digits under repeated refinement.

## Bottom line

For the first hierarchy level:

- the analytic mask is good enough to use directly
- the fitted table is not buying anything important at `rho ≈ 1.2`
- a fitted correction can remain a later optional refinement, not the default

So the first stored hierarchy mask should probably just be the explicit
analytic Gaussian quasi-refinement mask.

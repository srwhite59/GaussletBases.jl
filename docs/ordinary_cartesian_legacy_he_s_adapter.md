# Ordinary Cartesian Legacy He `s`-Gaussian Adapter Pass

## Status

Supporting note. Current relevance: active ordinary-branch basis-design note
for the first serious real-Gaussian supplement pass on top of the hybrid
ordinary branch.

## Goal

The calibrated residual-Gaussian interaction checks showed that, in the
paper-like core-spacing regime, the main remaining question is no longer the
bookkeeping of a tiny hand-picked Gaussian channel.

At that point, the hand-picked widths `[0.2, 0.6]` are no longer enough as a
serious hybrid test. The next step is to replace that toy supplement by a real
legacy-informed He `s`-Gaussian set.

This pass keeps the scope deliberately small:

- He only
- `s` functions only
- two standard basis choices:
  - `cc-pVTZ`
  - `cc-pVQZ`
- a tiny modern adapter for the legacy `BasisSets` / `ReadBasis.jl` format
- no broad Gaussian-basis subsystem

## Tiny adapter

The new narrow public entry point is:

```julia
legacy_s_gaussian_data("He", "cc-pVTZ")
```

It:

1. reads the legacy `BasisSets` file in the style of `ReadBasis.jl`
2. extracts only the `l = 0` shells
3. uses the primitive exponents directly
4. converts each legacy primitive

```text
exp(-zeta * r^2)
```

into the current one-dimensional Gaussian width convention

```text
exp(-0.5 * (x / width)^2),   width = 1 / sqrt(2 zeta)
```

## What is used in this first pass

This pass uses:

- all primitive `s` exponents from the chosen He basis
- uncontracted primitives only
- no diffuse filtering

So the old contraction coefficients are read but not used in the present
hybrid ordinary tests.

That is intentional. The first serious real-basis pass should avoid adding a
second subjective trimming rule at the same time.

## He `s` primitives used

For `cc-pVTZ`, the adapter uses all `s` primitives:

```text
zetas  = [234.0, 35.16, 7.989, 2.212, 0.6669, 0.2089]
widths = [0.046225016352102424, 0.11925059893763726, 0.2501720524494329,
          0.4754364132056024, 0.8658738891102185, 1.547090723905439]
```

For `cc-pVQZ`, the adapter uses all `s` primitives:

```text
zetas  = [528.5, 79.31, 18.05, 5.085, 1.609, 0.5363, 0.1833]
widths = [0.03075831259604325, 0.07940009594713392, 0.16643566632465154,
          0.31357362279453244, 0.5574513610066167, 0.9655640855770945,
          1.6515957995876271]
```

## Friendly-regime check

The main comparison was run in a paper-like hybrid regime:

- `count = 11`
- `s = 0.6`
- `xmax = 6`
- `dx_core ≈ 0.5166 bohr`
- `backend = :pgdg_localized_experimental`
- `Z = 2`

Results:

```text
pure ordinary:
  E1    = -1.9695420394204697
  <Vee> = 1.2392638967254135

hybrid with toy widths [0.2, 0.6]:
  E1    = -1.9897130930213345
  <Vee> = 1.2550226065837053

hybrid with He cc-pVTZ s primitives:
  E1    = -1.9997703172972017
  <Vee> = 1.2422720836593708

hybrid with He cc-pVQZ s primitives:
  E1    = -1.9999236696956322
  <Vee> = 1.2242858163500934
```

For the slightly looser `count = 11, s = 0.5` case:

```text
pure ordinary:
  E1    = -1.940025253727119
  <Vee> = 1.2356045484884919

hybrid toy:
  E1    = -1.9890163572242665
  <Vee> = 1.2703039947095944

hybrid cc-pVTZ s:
  E1    = -1.9997712332870128
  <Vee> = 1.3023008123593707

hybrid cc-pVQZ s:
  E1    = -1.9999235334654863
  <Vee> = 1.290552581728559
```

## Interpretation

This changes the practical picture in a useful way:

- the toy Gaussian supplement was not enough to judge the serious hybrid
  ordinary branch
- a real He `s` supplement is now easy to inject through the current hybrid
  interface
- with that real supplement, the one-body side becomes essentially hydrogenic
  for these small tests
- the remaining mixed behavior is on the interaction / basis-design side:
  realistic Gaussian support can improve `E1` strongly without automatically
  improving the same `1s^2` scalar

So the practical conclusion is:

- the hybrid ordinary branch is now using a realistic enough Gaussian
  supplement to be scientifically informative
- the next limiter is not the lack of a real Gaussian loader
- the next limiter is the choice of supplement and its interaction with the
  current density-density / IDA treatment

## Most natural next step

The next small step after this pass is not a He solver yet. It is likely one
of:

- a narrow filter/selection rule for the real He `s` supplement
- a comparison of toy vs legacy-informed supplements under the present hybrid
  interaction treatment
- only then a decision on whether the ordinary branch is ready for the first
  real He-style solver layer

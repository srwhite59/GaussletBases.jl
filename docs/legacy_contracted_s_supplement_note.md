# Contracted Legacy `s` Supplement Note

## Why this pass

The current legacy He supplement had still been using primitive `s` shells
directly. That kept the analytic primitive-integral route, but it also made
the added Gaussian channel more redundant than intended.

This pass switches the active legacy-informed supplement path to contracted
centered `s` functions while keeping the underlying primitive overlap /
kinetic / position / `x^2` / Gaussian-factor / pair-factor blocks analytic.
Contraction is applied only on the added-Gaussian side.

## Scope

This remains deliberately narrow:

- He only
- centered `s` supplement only
- legacy `BasisSets` / `ReadBasis.jl` style input
- active physical path now contracted
- primitive / uncontracted route retained only as a diagnostic mode

## Representation

The active object is still `LegacySGaussianData`, but it now carries both:

- primitive exponents / widths / centered primitive Gaussians
- a shell contraction matrix for the active supplement functions
- representative contracted widths / Gaussians for metadata and reporting

For the current He checks:

- `cc-pVTZ`: `6` primitives contracted to `3` centered `s` supplement functions
- `cc-pVQZ`: `7` primitives contracted to `4` centered `s` supplement functions

The contraction pattern follows the legacy `makeallcontractions(...)` style:
shell-by-shell contraction, with analytic primitive matrices underneath.

## Active consumers

The active paths now using the contracted supplement are:

- the ordinary hybrid surrogate path
- the base QW-PGDG reference path

The QW-PGDG path is the physically meaningful validation route here. The
hybrid combined-basis surrogate remains useful as a structural/diagnostic
path, but its raw `1s^2` scalar should not be treated as the reference target
for this contraction pass.

## Base QW-PGDG convergence check

Using:

- `basis = MappedUniformBasisSpec(:G10; s = 0.8, xmax = 6, reference_spacing = 1)`
- He `cc-pVTZ` contracted `s` supplement
- `Z = 2`
- `ordinary_cartesian_qiu_white_operators(...; interaction_treatment = :ggt_nearest)`

Results:

```text
count = 9:
  contracted supplement functions = 3
  kept residual modes             = 3
  E1                              = -1.9653729039662613
  <Vee>                           = 1.2571992738145692
  Δ(<Vee> - 1.25)                 = +0.007199273814569196
  constructor time                = 3.600558206 s

count = 11:
  contracted supplement functions = 3
  kept residual modes             = 2
  E1                              = -1.9936548803801342
  <Vee>                           = 1.2462429132727513
  Δ(<Vee> - 1.25)                 = -0.0037570867272487263
  constructor time                = 0.869931379 s

count = 13:
  contracted supplement functions = 3
  kept residual modes             = 2
  E1                              = -1.9971817407657657
  <Vee>                           = 1.247306574897531
  Δ(<Vee> - 1.25)                 = -0.0026934251024690603
  constructor time                = 2.08453247 s
```

## Interpretation

This is the useful result:

- the supplement channel is less redundant structurally
- the kept residual space drops from the earlier `3` modes to `2` modes at
  `count = 11` and `13`
- the `1s^2` scalar remains in the right physical regime and is already close
  to the `1.25` target
- the convergence pattern still looks clean rather than blocked by runtime or
  residual-selection pathology

So the contraction correction looks worthwhile. The active physical path is now
closer to the intended legacy model: analytic primitive integrals underneath,
but contracted `s` supplement functions at the added-channel level.

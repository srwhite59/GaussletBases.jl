# Ordinary Cartesian Residual-Gaussian Interaction Pass

## Status

Supporting note. Current relevance: active ordinary-branch validation note for
the hybrid/core-supported interaction treatment, still below the threshold for
opening a full ordinary He solver.

## Goal

The next ordinary-branch milestone after the first hybrid `1s^2` scalar checks
is not a He solve. It is a narrower interaction pass:

- keep the current hybrid basis interface
- keep the current one-body path
- separate the interaction contribution associated with the added core
  Gaussians from the current combined hybrid-basis treatment
- recheck the same doubly occupied noninteracting `1s^2` scalar for `Z = 2`

The scalar reference remains the hydrogenic value

```text
<Vee> = (5/8) Z = 1.25 Eh  for Z = 2.
```

## Why this comes before a He solve

The previous hybrid pass already showed that the core-supported route moves the
ordinary branch in the right direction:

- the lowest one-body orbital moved close to `-2`
- the `1s^2` IDA interaction moved much closer to `1.25`

But that earlier check was done through the current combined hybrid-basis
treatment. The next missing question was whether the interaction treatment of
the added Gaussian channel itself should be separated before attempting any He
solve.

## What is different in this pass

The current combined hybrid treatment builds the full hybrid basis first and
then applies the ordinary Cartesian density-density / IDA construction to that
combined localized basis.

The experimental residual-Gaussian treatment added here keeps the same hybrid
basis and one-body path, but changes the interaction side:

1. build the mapped ordinary backbone on the chosen backend
2. orthogonalize the added core Gaussians against that backbone to obtain
   residual-Gaussian directions
3. assign those residual-Gaussian directions to the nearest ordinary backbone
   centers for the purpose of the Gaussianized density-density / IDA pair
   factors
4. transfer that seed interaction back into the final localized hybrid basis

This is still a narrow validation ansatz, not a claim that the ordinary branch
has its final hybrid interaction model.

## Calibration lesson

The first hybrid scalar check at `count = 7` turned out to be too harsh a
regime to judge the residual-Gaussian approximation. In that case the pure
ordinary backbone itself was still too incomplete:

```text
count = 7
pure ordinary E1 = -1.4392418455757099
pure ordinary <Vee> = 0.9550147298987698
```

That is not yet the regime where the added Gaussian channel should be expected
to act only as a small correction. In the friendlier regime, the ordinary
backbone has to be closer to complete first.

## Calibrated friendly-regime checks

The first residual-Gaussian interaction tests used `s = 0.2` on small counts.
That turned out to be too coarse a calibration for the hybrid regime. In the
current code, `s` is a mapping-strength parameter rather than the actual
near-origin gausslet spacing, and the hybrid paper operates in a friendlier
core-spacing window closer to about `0.5`–`0.6` bohr.

The more relevant checks therefore use:

- `xmax = 6`
- explicit centered core Gaussians with widths `[0.2, 0.6]`
- `backend = :pgdg_localized_experimental`
- `Z = 2`
- mild maps chosen so the actual near-origin spacing `dx_core` is in the
  paper-like regime

Results:

```text
count = 9, s = 0.8, dx_core = 0.6507555285151014
pure ordinary:
  E1    = -1.9532986984681018
  <Vee> = 1.252447800305343

hybrid, combined basis treatment:
  E1    = -1.9886787450118553
  <Vee> = 1.3167775105956792

hybrid, residual-Gaussian interaction treatment:
  E1    = -1.9886787450118553
  <Vee> = 1.362965002197757
```

```text
count = 11, s = 0.5, dx_core = 0.6523335991349724
pure ordinary:
  E1    = -1.940025253727119
  <Vee> = 1.2356045484884919

hybrid, combined basis treatment:
  E1    = -1.9890163572242665
  <Vee> = 1.2703039947095944

hybrid, residual-Gaussian interaction treatment:
  E1    = -1.9890163572242665
  <Vee> = 1.3335876315874013
```

```text
count = 11, s = 0.6, dx_core = 0.5166101990265577
pure ordinary:
  E1    = -1.9695420394204697
  <Vee> = 1.2392638967254135

hybrid, combined basis treatment:
  E1    = -1.9897130930213345
  <Vee> = 1.2550226065837053

hybrid, residual-Gaussian interaction treatment:
  E1    = -1.9897130930213345
  <Vee> = 1.3133063073068527
```

## Interpretation

These calibrated checks change the interpretation of the residual-Gaussian
interaction pass:

- the earlier very small-`s` tests were useful mainly as a warning that
  underresolved hybrid IDA can behave badly
- once the ordinary backbone is in the paper-like `0.5`–`0.6` bohr core
  spacing regime, the pure ordinary branch is already quite good on the same
  scalar
- in that better-calibrated regime, the current combined hybrid-basis
  interaction treatment is actually closer to the hydrogenic `1.25` target
  than the present residual-Gaussian nearest-center ansatz
- the residual-Gaussian ansatz therefore remains useful as a diagnostic
  comparison, but it is not the preferred interaction treatment in the current
  friendly regime

That is still scientifically useful. It says that, once the backbone is
complete enough, the remaining limitation is more likely basis/orbital quality
than Gaussian interaction bookkeeping.

## Most natural next step

The next practical question is still not a He solve. It is:

- whether the ordinary hybrid branch should next adopt a small legacy-informed
  Gaussian-set loader for the friendly regime
- and then whether the resulting hybrid basis design is stable enough for the
  first real ordinary He-style solver layer

That loader question should still come after the interaction treatment is
understood, not before.

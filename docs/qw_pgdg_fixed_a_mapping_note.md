# QW-PGDG Fixed-`a` Mapping Note

## Why change families

The earlier fixed-`s` scans were useful as diagnostics, but they are not the
right convergence family for the Cartesian He QW-PGDG tests.

If `s` is held fixed while `count` increases, then the near-origin mapping
scale does not tighten in the natural way. For a more physical convergence
family, the near-origin asinh parameter `a` should be held fixed and `s`
should decrease as `count` increases so that the outer centers still land at
the chosen `±xmax`.

That makes:

- larger `count`  -> smaller `s`
- smaller `c = a s`
- smaller near-origin physical spacing

which is the right qualitative direction.

## First practical family

For He (`Z = 2`), the first practical Cartesian test family is taken from the
radial-style near-origin rule

```text
a = 1 / (2 Z) = 1 / 4.
```

For fixed `a`, odd `count`, and chosen `xmax`, the endpoint condition is still

```text
u(xmax) = (count - 1) / 2
```

with

```text
u(x) = x / tail_spacing + asinh(x / a) / s.
```

So `s` is determined directly by

```text
s = asinh(xmax / a) / ( ((count - 1) / 2) - xmax / tail_spacing ).
```

The implied Cartesian near-origin scale is then

```text
c = a s.
```

For these scans, `tail_spacing = 10` and `xmax = 10`.

## Base QW-PGDG fixed-`a` scan

Using:

- base PGDG intermediate only
- He `cc-pVTZ` contracted `s` supplement
- `interaction_treatment = :ggt_nearest`
- `a = 1/(2Z) = 1/4`
- `xmax = 10`

results are:

```text
count = 11:
  solved s         = 1.095545712016
  implied c = a s  = 0.273886428004
  dx_core          = 0.317784852947
  kept RG modes    = 1
  E1               = -1.981038528584
  <Vee>            = 1.231025632122
  Δ(<Vee> - 1.25)  = -0.018974367878
  constructor time = 4.208 s

count = 13:
  solved s         = 0.876436569613
  implied c = a s  = 0.219109142403
  dx_core          = 0.240874372377
  kept RG modes    = 2
  E1               = -1.994592509794
  <Vee>            = 1.247633256131
  Δ(<Vee> - 1.25)  = -0.002366743869
  constructor time = 2.183 s

count = 15:
  solved s         = 0.730363808011
  implied c = a s  = 0.182590952003
  dx_core          = 0.194735774804
  kept RG modes    = 1
  E1               = -1.997177156085
  <Vee>            = 1.247886543106
  Δ(<Vee> - 1.25)  = -0.002113456894
  constructor time = 4.487 s

count = 17:
  solved s         = 0.626026121152
  implied c = a s  = 0.156506530288
  dx_core          = 0.163856561294
  kept RG modes    = 2
  E1               = -1.998588440029
  <Vee>            = 1.248965557065
  Δ(<Vee> - 1.25)  = -0.001034442935
  constructor time = 9.263 s
```

## Interpretation

This family is cleaner than the earlier fixed-`s` scan in the ways that
matter most:

- the mapping parameters now move in the physically intended direction
- the near-origin spacing decreases systematically with count
- the higher-count `⟨Vee⟩` sequence (`13/15/17`) trends more cleanly toward
  `1.25`
- `E1` also moves toward `-2` as count increases

`count = 11` is still a bit coarse for this `xmax = 10`, but from `13`
upward the family looks like a sensible pre-nesting convergence line.

So this fixed-`a` family should replace the old fixed-`s` line as the default
Cartesian pre-nesting convergence family.

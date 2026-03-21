# Base PGDG GA One-Body Cross Correction

This note records the next correction on the base
`refinement_levels = 0` PGDG intermediate layer.

The previous validation/timing pass made two things clear:

- the base PGDG path is back in a physically meaningful `1s^2` regime
- the dominant runtime cost is still the active one-body-side `ga`
  cross-block construction

That same active midpoint path was also the most plausible source of the
remaining too-low one-body energy.

So the next correction is narrow:

- keep the shared PGDG intermediate layer intact
- keep the already-fast analytic/proxy `pair_ga` route
- move the remaining one-body-side `ga` blocks onto the same contracted
  analytic/proxy raw-space route

For the base PGDG layer, that means the active construction of

- `overlap_ga`
- `kinetic_ga`
- `position_ga`
- `x2_ga`
- `factor_ga`

should no longer rely on midpoint quadrature.

Instead, the intended structure is:

1. start from the base PGDG proxy/raw Gaussian line
2. build analytic primitive Gaussian cross blocks against the added Gaussian
   supplement
3. contract only on the gausslet side
4. feed those contracted `ga` blocks into the existing Qiu-White raw-space
   assembly

This keeps the clean layer separation and targets the actual remaining hot path
without reopening the refinement-hierarchy question yet.

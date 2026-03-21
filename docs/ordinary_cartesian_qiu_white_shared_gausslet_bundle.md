# Qiu-White Shared Gausslet Bundle Note

This note records the next refactor after the Qiu-White split-block correction.

The remaining bottleneck is now the gausslet-side one-dimensional block
construction itself:

- overlap `gg`
- kinetic `gg`
- position `gg`
- `x^2` `gg`
- Gaussian-factor `gg`
- pair-factor `gg`

That work is not specifically Qiu-White logic. It is generic mapped ordinary
infrastructure.

So the next correction in this pass is to pull that gausslet-side one-
dimensional bundle out of the Qiu-White constructor and make it shared.

The intended split remains:

- `gg` from shared mapped ordinary infrastructure
- `aa` analytically
- `ga` from the dedicated cross-block route
- then raw-space block matrices assembled from those pieces

After this refactor, the Qiu-White reference path should be responsible mainly
for:

- Gaussian-Gaussian analytic blocks
- gausslet-Gaussian cross blocks
- residual-space construction
- raw-to-final transforms
- nearest / GGT and MWG interaction handling

It should no longer privately rebuild the generic mapped-gausslet 1D bundle.

## Updated timing interpretation

On the same light paper-like He case

- `count = 9`
- `s = 0.8`
- `xmax = 6`
- `cc-pVTZ`
- `interaction_treatment = :mwg`

the refactored Qiu-White constructor still did not finish phase 1 within a
two-minute observation window.

Direct timings of the shared gausslet-side bundle internals show why:

- `basis_representation(...)`: about `15.12 s`
- `_x2_matrix(basis)`: about `18.97 s`
- `gaussian_factor_matrices(...)`: about `5.08 s`
- `_pair_gaussian_factor_matrices(...)`: still running after `60 s`

So the dominant remaining cost is not yet the later Qiu-White-specific dense
3D assembly or the split `ga` cross-block route. It is still the shared
mapped-gausslet pair-factor construction.

That means the refactor is still worthwhile:

- the Qiu-White path no longer owns or hides that generic cost
- the next optimization target is now clearly the shared gausslet-side
  pair-factor path
- only after improving that shared bundle will it be meaningful to judge the
  remaining dense 3D reference cost

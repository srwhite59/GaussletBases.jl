# Base PGDG Layer Validation And Timing Note

This note records the next narrow step after separating the base
`refinement_levels = 0` PGDG intermediate layer.

We are deliberately holding off on `1/9` and higher refinement levels for now.
The priority is the base PGDG auxiliary line itself:

- recover a clean physical `1s^2` validation scalar
- identify the real runtime bottleneck in the current full pipeline

The base-PGDG layer is now explicit shared infrastructure, so the immediate
question is no longer "how should the hierarchy work?" It is:

- does the base `1/3` PGDG auxiliary line support a sensible He-like `1s^2`
  check?
- and if runtime is still poor, is the remaining bottleneck now the true
  low-level Gaussian kernel or some higher numerical driver?

## Clean base-layer validation target

For the light He reference case

- basis family: `MappedUniformBasisSpec(:G10, count = 9, s = 0.8, xmax = 6)`
- Gaussian supplement: He `cc-pVTZ` `s` primitives
- Coulomb expansion: `doacc = false`
- nuclear charge: `Z = 2`
- PGDG layer: base `refinement_levels = 0`

the cleanest current base-PGDG test is the paper-faithful Qiu-White nearest /
GGT route built on the shared base PGDG intermediate line.

That run gives:

- lowest one-body orbital energy `E1 = -2.5258145386216424`
- `⟨Vee⟩ = 1.282969072925221`
- difference from the hydrogenic `1.25` target: `+0.032969072925221`
- final overlap error: about `4.0e-14`

So the base-PGDG route is now physically meaningful again on the `1s^2`
validation scalar. It is not exact, and the one-body energy is still too low,
but the `Vee` check is in the right regime.

## Timing diagnosis on the full base-PGDG pipeline

On the same light case, the constructor timing hook

- `ordinary_cartesian_qiu_white_operators(...; timing = true)`

reports approximately:

- shared gausslet-side 1D bundle: `16.26 s`
- split gausslet-Gaussian raw-block assembly: `93.82 s`
- 3D gausslet overlap assembly: `0.038 s`
- residual-space construction: `0.116 s`
- 3D gausslet one-body assembly: `0.676 s`
- raw one-body transform: `0.047 s`
- raw moment-matrix assembly: `0.151 s`
- 3D gausslet interaction assembly: `0.023 s`
- RG interaction assembly: `0.002 s`
- total: `111.14 s`

So the dominant cost is no longer the broad shared-bundle structure or the 3D
driver. It is overwhelmingly the split `ga` raw-block phase.

## Is the true hotspot now the low-level integral kernel?

For the current base PGDG layer, the answer is **no**.

The analytic low-level PGDG proxy kernels are already cheap on the same case:

- proxy scalar bundle total: about `1.89 s`
- primitive Gaussian-factor build: about `0.082 s`
- Gaussian-factor contraction: about `0.056 s`
- primitive pair-factor build: about `0.112 s`
- pair-factor contraction: about `0.045 s`

By contrast, directly timing the active Qiu-White `ga` cross build gives:

- `_qwrg_cross_1d_blocks_midpoint(...; include_pair = false)`: about `94.61 s`
- analytic `ga` pair primitive build: about `0.082 s`
- left contraction of those pair blocks: about `0.041 s`

That means the dominant remaining cost is **not** the analytic three-Gaussian /
Gaussian-factor / pair-factor primitive kernel layer.

It is the still-active midpoint-grid cross-quadrature path for the one-body-side
`ga` blocks:

- overlap `ga`
- kinetic `ga`
- position `ga`
- `x^2` `ga`
- Gaussian-factor `ga`

inside:

- `_qwrg_cross_1d_blocks_midpoint(...)`

## Bottom line

The base `refinement_levels = 0` PGDG intermediate layer is now the right
place to stand:

- the `1s^2` scalar is physically sensible again
- the shared PGDG auxiliary line is no longer the main architectural mistake
- the real hotspot is now clearly the midpoint `ga` cross-block driver, not
  the low-level analytic PGDG proxy kernel

So the next optimization target should be:

- replace the midpoint one-body-side `ga` cross-block construction with a
  contracted analytic / proxy raw-space route

not:

- more hierarchy rollout
- or premature optimization of the already-cheap low-level PGDG proxy kernels

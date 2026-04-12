# Gaussian Refinement Shape-Ratio Follow-up

This note records the next narrow numerical study for the proposed
one-dimensional Gaussian-refinement hierarchy.

The first two studies established that the refinement idea is numerically
feasible, but they only tested the baseline choice

- `rho = width / spacing = 1`

That is no longer the main question.

The next practical question is:

- how do accuracy, locality, and conditioning change when the finer Gaussians
  are made modestly wider than their spacing?

## Legacy clue

There is a clear legacy hint that `rho > 1` should be tested seriously:

- `InterpGau.jl` in `GaussletModules` contains
  `const widthratio = sqrt(2.0)`
- `testrt2.jl` uses the `cfrt2` / `Gaussletrt2` line with
  `sqrt(2.0)`

This pass does **not** adopt `sqrt(2)` blindly. It uses that older choice as a
numerical hint and compares it against smaller increases such as `1.1`, `1.2`,
and `4/3`.

## Study setup

The same temporary script was extended:

- `tmp/work/gaussian_refinement_feasibility.jl`

Main arithmetic:

- `BigFloat`, `256` bits

Studied shape ratios:

- `1.0`
- `1.1`
- `1.2`
- `4/3`
- `sqrt(2)`
- `1.5`

Studied refinement spacings for the coarse-Gaussian fit:

- `1/3`
- `1/9`

The script measures two things.

### A. Constant-lattice ripple

For a uniform infinite Gaussian lattice with spacing `h` and width `rho * h`,
the mean value over one cell is exactly

`sqrt(2π) * rho`

so the normalized lattice sum is sampled over one cell and the maximum ripple
from `1` is reported.

This quantity is scale-invariant under the `width = rho * spacing` rule, so it
depends only on `rho`, not on the absolute spacing.

### B. Coarse-Gaussian refinement quality

The unit-width centered coarse Gaussian is fit by the finer Gaussian lattice
using the same exact continuous `L2` normal equations as before.

For each `(spacing, rho)` pair, the script reports:

- relative `L2` error
- max pointwise error
- condition estimate
- solve residual
- coefficient `2`-norm
- coefficient RMS radius
- edge coefficient ratio
- smallest practical saturated window from a BigFloat sweep

## Constant-ripple results vs `rho`

Results:

| `rho` | max relative ripple | RMS ripple |
| --- | ---: | ---: |
| `1.0` | `5.35e-9` | `3.78e-9` |
| `1.1` | `8.48e-11` | `5.99e-11` |
| `1.2` | `9.05e-13` | `6.40e-13` |
| `4/3` | `1.15e-15` | `8.13e-16` |
| `sqrt(2)` | `1.43e-17` | `1.01e-17` |
| `1.5` | `1.03e-19` | `7.28e-20` |

So even a modest increase above `1` helps a lot:

- `1.1` improves the constant-lattice ripple by about two orders of magnitude
- `1.2` improves it by about four orders
- `4/3` and `sqrt(2)` make the ripple essentially negligible at displayed
  precision

## Coarse-Gaussian fit results vs `rho`

### Spacing `1/3`

| `rho` | chosen window | relative `L2` error | max pointwise error | cond | `||c||₂` | coeff RMS radius | edge ratio |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `1.0` | `±6` | `3.39e-8` | `4.80e-8` | `8.30e3` | `9.47e-1` | `6.67e-1` | `1.94e-9` |
| `1.1` | `±6` | `1.49e-9` | `2.10e-9` | `6.07e4` | `8.67e-1` | `6.58e-1` | `1.15e-9` |
| `1.2` | `±6` | `6.05e-11` | `8.55e-11` | `5.22e5` | `8.01e-1` | `6.48e-1` | `6.28e-10` |
| `4/3` | `±7` | `8.33e-13` | `1.18e-12` | `1.36e7` | `7.29e-1` | `6.33e-1` | `6.57e-14` |
| `sqrt(2)` | `±7` | `6.54e-14` | `9.24e-14` | `1.05e8` | `6.93e-1` | `6.24e-1` | `2.50e-14` |
| `1.5` | `±7` | `4.83e-15` | `6.83e-15` | `1.01e9` | `6.59e-1` | `6.12e-1` | `7.80e-15` |

### Spacing `1/9`

| `rho` | chosen window | relative `L2` error | max pointwise error | cond | `||c||₂` | coeff RMS radius | edge ratio |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `1.0` | `±6` | `4.97e-9` | `6.83e-9` | `9.50e3` | `1.60e0` | `7.03e-1` | `2.98e-8` |
| `1.1` | `±7` | `8.56e-11` | `1.21e-10` | `7.54e4` | `1.45e0` | `7.02e-1` | `4.20e-11` |
| `1.2` | `±8` | `1.06e-12` | `1.50e-12` | `7.28e5` | `1.33e0` | `7.01e-1` | `2.03e-14` |
| `4/3` | `±8` | `1.77e-15` | `2.48e-15` | `2.02e7` | `1.20e0` | `6.99e-1` | `2.25e-14` |
| `sqrt(2)` | `±9` | `2.68e-17` | `3.79e-17` | `1.81e8` | `1.13e0` | `6.98e-1` | `3.41e-18` |
| `1.5` | `±9` | `2.50e-19` | `3.54e-19` | `2.11e9` | `1.07e0` | `6.97e-1` | `3.50e-18` |

## Locality / conditioning tradeoff

The numerical tradeoff is now very clear.

Benefits of increasing `rho`:

- constant reproduction improves extremely quickly
- coarse-Gaussian refinement error drops by orders of magnitude
- coefficient norms decrease modestly
- coefficient RMS radii improve slightly rather than worsening

Costs of increasing `rho`:

- the overlap matrix condition number grows very quickly

The conditioning trend is the main practical limiter:

- around `8e3` to `1e4` at `rho = 1`
- around `6e4` to `8e4` at `rho = 1.1`
- around `5e5` to `7e5` at `rho = 1.2`
- around `1e7` to `2e7` at `4/3`
- around `1e8` at `sqrt(2)`
- around `1e9` at `1.5`

So the accuracy continues to improve, but past about `1.2` the condition
number starts rising very fast.

## Recommendation

For a first practical refined-PGDG hierarchy choice, the best compromise looks
to be:

- `rho ≈ 1.2`

Reason:

- it already improves constant-lattice ripple from `5e-9` to about `1e-12`
- it improves the coarse-Gaussian fit by roughly:
  - `3.39e-8 -> 6.05e-11` at spacing `1/3`
  - `4.97e-9 -> 1.06e-12` at spacing `1/9`
- it does this before the condition number reaches the `1e7` to `1e8` regime

`4/3` and `sqrt(2)` are numerically impressive, and the legacy `sqrt(2)` hint
clearly points in the right direction, but they probably overshoot the best
first practical compromise if locality and conditioning both matter.

So the practical ranking from this study is:

1. `rho ≈ 1.2` as the best first practical choice
2. `rho = 4/3` as a stronger-accuracy but more ill-conditioned option
3. `rho = sqrt(2)` as a high-precision legacy-style option when one is willing
   to tolerate much worse conditioning

## Bottom line

Spacing refinement alone is not the whole story.

The shape ratio matters a lot, and the numerical study says:

- yes, `rho > 1` is materially better
- yes, the legacy `sqrt(2)` clue was meaningful
- but for a first practical hierarchy level, a milder increase such as
  `rho ≈ 1.2` looks like the better starting point

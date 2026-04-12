# Gaussian Refinement Precision Follow-up

This note records the second, higher-precision pass on the one-dimensional
Gaussian-refinement study.

The first study already established the main feasibility result:

- `1/3` looks very good
- the coefficients are local and smooth
- the refinement hierarchy looks numerically promising

The remaining question was narrower:

- how much of the observed error is real approximation/truncation error
- and how much is just Float64 or solve sensitivity

## Purpose of this pass

This pass separates two effects explicitly:

- precision / solve error
- finite-window truncation error

The fitting formulation is unchanged from the first note:

- Gaussian convention:
  `exp(-0.5 * ((x - center) / width)^2)`
- target Gaussian:
  centered, width `1`
- fine fitting Gaussians:
  width = spacing = refinement level
- coefficients from the exact continuous `L2` normal equations
  `S c = b`

## Precision used

Main comparisons:

- `Float64`
- `BigFloat` at `256` bits

Confirmation:

- `BigFloat` at `512` bits for the largest centered `1/9` case

Condition numbers are reported as a stable estimate from the Float64 copy of
the same SPD matrix. The more direct precision diagnostic is the linear-system
solve residual, which drops to about `1e-77` at `256` bits and `1e-154` at
`512` bits.

## Fixed-window comparison at `±5`

This comparison keeps the finite window fixed and only changes arithmetic.

### `1/3`

- Float64:
  - relative `L2` error `4.614844855560e-08`
  - max pointwise error `4.756906712977e-08`
  - solve residual `6.82e-17`
- BigFloat 256:
  - relative `L2` error `3.854784902547e-08`
  - max pointwise error `4.756906850596e-08`
  - solve residual `2.37e-77`
- Float64 vs BigFloat coefficient difference:
  - relative `2`-norm difference `1.40e-13`

### `1/5`

- Float64:
  - relative `L2` error `1.397960639515e-07`
  - max pointwise error `2.182953457287e-07`
  - solve residual `1.52e-16`
- BigFloat 256:
  - relative `L2` error `1.401675367661e-07`
  - max pointwise error `2.182953457531e-07`
  - solve residual `1.87e-77`
- Float64 vs BigFloat coefficient difference:
  - relative `2`-norm difference `2.09e-13`

### `1/9`

- Float64:
  - relative `L2` error `4.005965543109e-07`
  - max pointwise error `7.633958258338e-07`
  - solve residual `1.32e-16`
- BigFloat 256:
  - relative `L2` error `4.005092013170e-07`
  - max pointwise error `7.633958258350e-07`
  - solve residual `2.41e-77`
- Float64 vs BigFloat coefficient difference:
  - relative `2`-norm difference `1.95e-13`

## What the fixed-window comparison says

The fixed-window non-monotone trend does **not** disappear in higher
precision.

That is an important result:

- the coefficients are essentially the same in Float64 and BigFloat
- the pointwise errors are essentially the same
- the solve residuals collapse to numerical zero in BigFloat

So the main limitation at fixed `±5` is not Float64 arithmetic and not the
normal-equation solve. It is the finite fitting window.

## BigFloat window sweep

Using `256`-bit arithmetic, the window sweeps were:

### `1/3`

- `±4`: relative `L2` error `4.71e-06`
- `±5`: relative `L2` error `3.85e-08`
- `±6`: relative `L2` error `3.39e-08`
- `±7`: relative `L2` error `3.39e-08`
- `±8`: relative `L2` error `3.39e-08`

So `1/3` is already saturated by about `±6`.

### `1/5`

- `±4`: relative `L2` error `2.19e-05`
- `±5`: relative `L2` error `1.40e-07`
- `±6`: relative `L2` error `8.34e-09`
- `±7`: relative `L2` error `8.33e-09`
- `±8`: relative `L2` error `8.33e-09`

So `1/5` needs one more step beyond `±5`, then saturates near `8.33e-09`.

### `1/9`

- `±4`: relative `L2` error `5.01e-05`
- `±5`: relative `L2` error `4.01e-07`
- `±6`: relative `L2` error `4.97e-09`
- `±7`: relative `L2` error `4.83e-09`
- `±8`: relative `L2` error `4.83e-09`

So `1/9` also needs a slightly larger window, then saturates near
`4.83e-09`.

## 512-bit confirmation

For the largest centered case tested here

- refinement `1/9`
- half-window `±8`

the `512`-bit result agrees with the `256`-bit result to the printed digits:

- relative `L2` error `4.827479847916e-09`
- max pointwise error `6.827087435065e-09`

That confirms the remaining digits there are not being set by the BigFloat
solver either.

## Interpretation

In sufficiently high precision and with a large enough window, the monotone
trend does reappear:

- `1/5` is genuinely better than `1/3`
- `1/9` is genuinely better than `1/5`

So the earlier non-monotone centered errors were mainly a finite-window
artifact, not a Float64 artifact.

## Practical digit count

A practical reading from the saturated BigFloat runs is:

- `1/3`: about `7` to `8` trustworthy digits
- `1/5`: about `8` trustworthy digits
- `1/9`: about `8` to almost `9` trustworthy digits

Those are approximation/truncation-limited digits for this building block, not
solve-limited digits.

## Bottom line

Yes, `1/3` still looks good enough for a practical first hierarchy level.

More specifically:

- `1/3` is already excellent as a coarse practical refinement level
- `1/5` and `1/9` do improve the approximation when the fitting window is
  large enough
- the present digits are mainly window-limited, not Float64-limited

So the hierarchy still looks numerically promising enough to pursue further,
and the next narrow question should be about how these Gaussian-scale fits map
onto the actual distorted-gausslet proxy line, not about whether the Gaussian
building block itself is numerically viable.

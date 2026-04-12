# Gaussian Refinement Mask-Quality Follow-up

This note records the next narrow numerical study for the proposed
Gaussian-refinement hierarchy.

The earlier studies already established the main feasibility point:

- refinement by a finer uniform Gaussian lattice works very accurately
- `rho > 1` materially improves the approximation quality

Once the refinement masks are tabulated offline in high precision, the linear
system condition number matters mainly for table generation. The more relevant
practical question then becomes:

- are the stored refinement masks themselves numerically healthy, or do they
  work by large positive/negative cancellation?

So this pass focuses on:

- coefficient locality
- sign structure
- cancellation
- constant-reproduction quality of the stored mask
- robustness under coefficient rounding

## Setup

The same temporary study script was extended:

- `tmp/work/gaussian_refinement_feasibility.jl`

Main arithmetic:

- `BigFloat`, `256` bits

The study uses the first hierarchy mask only:

- refinement spacing `1/3`

This is the right practical object to tabulate first. Repeated application of
that same mask then probes what happens at the `1/9` and `1/27` levels.

The Gaussian convention is unchanged:

- `exp(-0.5 * ((x-center)/width)^2)`

For a chosen shape ratio `rho`, the fitted coarse Gaussian has:

- coarse spacing `1`
- width `rho`

and the refining fine lattice has:

- spacing `1/3`
- width `rho/3`

Studied mask choices:

- `rho = 1.0` as a baseline
- `rho = 1.2`
- `rho = 4/3`
- `rho = sqrt(2)`

## Mask-quality diagnostics

For each mask, the script reports:

- the actual coefficient vector
- `max(abs(c))`
- `||c||₁`
- `||c||₂`
- signed sum
- absolute sum
- cancellation ratio `sum(abs(c)) / abs(sum(c))`
- number of sign changes
- edge coefficient size
- RMS radius

For constant reproduction, the relevant mask-level quantity is not the earlier
uniform-lattice ripple by itself. Instead, for the `1/3` refinement mask with
offsets `k`, define the three residue sums

- `a_p = sum(c_k for k ≡ p mod 3)`

If these three values are equal, a coarse constant refined by the stored mask
produces a uniform fine coefficient field. The reported constant-reproduction
ripple is therefore

- `max_p |a_p / mean(a) - 1|`

This is the right discrete constant-reproduction diagnostic for the stored
mask itself.

## Main mask-health result

The key conclusion is simple:

- all of the useful masks are symmetric and strictly positive
- the cancellation ratio is essentially `1`
- the masks are not working by delicate positive/negative cancellation

So once the masks are tabulated offline, the main remaining tradeoff is:

- accuracy versus support width/locality

not:

- cancellation fragility

## Summary table

| `rho` | window | `n` | `max|c|` | `||c||₁` | `||c||₂` | cancellation | sign changes | edge | RMS radius | residue ripple |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `1.0` | `±6` | `37` | `4.23e-1` | `3.0` | `9.47e-1` | `1.0` | `0` | `8.21e-10` | `6.67e-1` | `4.81e-8` |
| `1.2` | `±8` | `49` | `3.53e-1` | `3.0` | `8.65e-1` | `1.0` | `0` | `6.80e-12` | `8.00e-1` | `2.05e-11` |
| `4/3` | `±10` | `61` | `3.17e-1` | `3.0` | `8.20e-1` | `1.0` | `0` | `8.58e-15` | `8.89e-1` | `5.39e-14` |
| `sqrt(2)` | `±10` | `61` | `2.99e-1` | `3.0` | `7.97e-1` | `1.0` | `0` | `3.21e-13` | `9.43e-1` | `1.71e-13` |

Interpretation:

- all masks are positive and normalized cleanly
- increasing `rho` broadens the mask in a smooth way
- the discrete constant-reproduction ripple improves strongly with `rho`
- there is no sign-oscillation penalty hiding behind the higher-`rho` masks

## Representative coefficient vectors

By symmetry, it is enough to list the nonnegative offsets.

### `rho = 1.2`

Centered half-mask:

```text
[0.3526184897178, 0.3376412458658, 0.2964194747957, 0.2385936049246,
 0.1760806735342, 0.1191421433670, 0.07391285076557, 0.04204122126587,
 0.02192459199761, 0.01048307138134, 0.004595644142561, 0.001847162346556,
 0.0006807138239668, ...]
```

### `rho = 4/3`

Centered half-mask:

```text
[0.3173566407456, 0.3063934134781, 0.2757242276848, 0.2312778822818,
 0.1808243632112, 0.1317781650643, 0.08951455082834, 0.05667708254755,
 0.03344914384415, 0.01840035469173, 0.009434764244314, 0.004509200843061,
 0.002008777224947, ...]
```

### `rho = sqrt(2)`

Centered half-mask:

```text
[0.2992067103013, 0.2900010876019, 0.2640489950735, 0.2258530741158,
 0.1814780433896, 0.1369868140415, 0.09713819674973, 0.06470798911954,
 0.04049322488526, 0.02380473887634, 0.01314622537063, 0.006820171875694,
 0.003323886309508, ...]
```

All three are:

- smooth
- positive
- rapidly decaying

The main change as `rho` increases is simply a wider positive tail.

## Rounding robustness

Each high-precision reference mask was rounded to:

- `Float64`
- `16` decimal digits
- `12` decimal digits
- `8` decimal digits

The rounded masks were then tested on:

- mask constant reproduction via the residue ripple
- one-level coarse-Gaussian refinement error

### One-level rounding results

#### `rho = 1.2`

- reference: residue ripple `2.05e-11`, relative `L2` error `1.51e-11`
- `Float64`: essentially identical
- `12` digits: relative `L2` error `1.51e-11`
- `8` digits: relative `L2` error degrades to `1.01e-8`

#### `rho = 4/3`

- reference: residue ripple `5.39e-14`, relative `L2` error `4.01e-14`
- `Float64`: essentially identical
- `12` digits: relative `L2` error `4.99e-13`
- `8` digits: relative `L2` error `4.37e-9`

#### `rho = sqrt(2)`

- reference: residue ripple `1.71e-13`, relative `L2` error `3.29e-15`
- `Float64`: essentially identical
- `12` digits: relative `L2` error `3.90e-13`
- `8` digits: relative `L2` error `7.32e-9`

So:

- `Float64` storage is already effectively perfect for these masks
- `12` decimal digits still preserve the masks extremely well
- `8` digits are visibly worse, though still not catastrophic

## Repeated-refinement check

One simple accumulation test was added by repeatedly composing the same
`1/3` mask for:

- `2` levels (`1 -> 1/9`)
- `3` levels (`1 -> 1/27`)

The resulting refined Gaussian was compared on a dense grid against the
original coarse Gaussian.

### `rho = 1.2`

- level `3`, reference: relative RMS grid error `1.51e-11`
- level `3`, `Float64`: essentially identical
- level `3`, `12` digits: `1.52e-11`
- level `3`, `8` digits: `2.47e-8`

### `rho = 4/3`

- level `3`, reference: relative RMS grid error `4.02e-14`
- level `3`, `Float64`: essentially identical
- level `3`, `12` digits: `5.13e-13`
- level `3`, `8` digits: `5.64e-9`

### `rho = sqrt(2)`

- level `3`, reference: relative RMS grid error `2.61e-15`
- level `3`, `Float64`: essentially identical
- level `3`, `12` digits: `1.14e-12`
- level `3`, `8` digits: `9.95e-9`

So repeated refinement does **not** show bad accumulation when the masks are
stored at `Float64` or `12`-digit quality. The first clearly problematic
storage level is closer to `8` digits.

## Practical recommendation

The main practical conclusion is now stronger than before.

The first hierarchy mask should be chosen mainly by:

- locality
- support width
- achieved accuracy

not by fear of cancellation.

My recommendation for the first tabulated hierarchy mask is still:

- `rho ≈ 1.2`

Reason:

- it is already strictly positive and cancellation-free
- it improves the discrete constant reproduction from about `5e-8` to `2e-11`
- it improves the single-level refinement error to about `1.5e-11`
- it stays materially more local than `4/3` or `sqrt(2)`
- it remains robust under `Float64` or `12`-digit storage, even after `3`
  refinement levels

What changes after this mask-health study is the interpretation of the higher
`rho` choices:

- `4/3` becomes entirely acceptable if one wants a stronger-accuracy stored
  mask and is willing to accept a wider support
- `sqrt(2)` is also numerically healthy once tabulated offline, so the legacy
  clue was real
- but neither `4/3` nor `sqrt(2)` looks like the best **first** practical
  default, because the extra accuracy comes mainly by widening the positive
  mask rather than by fixing some cancellation pathology

## Bottom line

For the proposed refinement hierarchy:

- the masks are numerically healthy
- they do not rely on cancellation
- `Float64` tabulation is already good enough
- `12` decimal digits also look safe
- `rho ≈ 1.2` still looks like the best first practical default

So the next hierarchy work can treat mask health as a solved first-order
question and move on to how these tabulated masks interact with the distorted
gausslet proxy construction.

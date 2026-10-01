# Radial Multipole Stabilization Milestone

This note records the narrow operator-side robustness milestone introduced by
commit `346c4c6` (`Stabilize radial multipoles and record prototype
milestone`).

This is not another prototype milestone. The prototype milestone was
`53a16a4`, which introduced the named cached paper-parity radial prototype.
This follow-on milestone is about the downstream radial multipole / `Vee`
builder in `src/radial/operators.jl`.

## Problem before the patch

Before this patch, the integral-diagonal radial multipole builder used raw
power factors directly inside the prefix/suffix accumulation:

- `r^L`
- `r^-(L+1)`

That behavior was numerically fragile for larger multipole order `L` or wide
radial ranges. In those regimes, the raw-power path could lose stability and
produce non-finite intermediate or final data, allowing `Inf` or `NaN` values
to leak into downstream radial multipoles.

## Fix in `src/radial/operators.jl`

The current builder keeps the same prefix/suffix integral structure, but the
power handling is now stabilized.

The new path:

- works from `log(r)`
- builds scaled representations of:
  - `r^L`
  - `r^-(L+1)`
- accumulates the prefix/suffix data in scaled form
- recovers final matrix entries through an explicit guarded recovery helper
- folds in the end-of-builder normalization by basis integral weights through
  log-magnitude plus explicit sign handling

The old raw-power implementation is retained only as an internal comparison
helper for tests.

The stabilized builder also throws explicitly if scaled accumulation or
recovery still goes non-finite, rather than letting `Inf` or `NaN` propagate
silently downstream.

## What stayed the same

On ordinary benign cases, the stabilized builder agrees with the old raw
builder to numerical roundoff.

So this patch is not intended as a behavior change for well-conditioned
regimes. It is a robustness fix for the wide-range / high-`L` failure mode.

## Regression coverage

The new focused regression coverage in `test/runtests.jl` checks two things:

- benign-case parity:
  - the stabilized builder agrees with the old raw builder in a well-behaved
    regime
- risky synthetic case:
  - a wide-range, high-`L` case where the raw builder goes non-finite
  - while the stabilized builder remains finite and symmetric

That is the main trust story for this patch: unchanged behavior where the old
path was healthy, and materially improved robustness where the old path could
break.

## Why this matters downstream

This is not just an isolated helper cleanup.

The stabilized path builds radial multipole matrices, and those feed the
current two-index IDA / downstream `Vee` operator story. That makes this part
of the trusted operator path for current radial/angular atomic work, not just a
private numerical tweak.

## Validation status

Validation for the stabilization patch was already reported as passing through
the radial test group before this note was written.

The relevant command was:

- `env GAUSSLETBASES_TEST_GROUPS=radial JULIA_DEPOT_PATH=/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/tmp/julia_depot julia --project=. test/runtests.jl`

This note itself did not rerun tests in the current turn. So the trust
statement here relies on that already-reported passing radial validation for
the stabilization patch, rather than on a fresh rerun tied to writing this
note.

## Relation to the prototype milestone

- `53a16a4` was the named cached paper-parity prototype milestone
- `346c4c6` is a follow-on operator-side stabilization milestone affecting
  radial multipole / downstream `Vee` robustness

The prototype milestone settled the manuscript radial object itself. This
stabilization milestone hardens a later operator-building layer that consumes
radial bases and quadrature data.

## Follow-up: high-L underflow (2026-10)

The scaled path above removed the overflow (`Inf` / `NaN`) failure but not the
matching underflow. It scaled `r^L` by its global maximum (`r_max^L`) and
`r^-(L+1)` by its global maximum (`r_min^-(L+1)`, with `r_min` clamped to
`eps()`), so every inner product carried a factor of order
`(r_min / r_max)^L`. Once `L * log(r_max / r_min)` exceeded about 700 the
scaled sums fell below `floatmin`. `_recover_scaled_kernel_value` then returned
`0.0` silently. On the Hooke erf grids (`r` from about 1e-23 to 80) this happened
at `L = 16-17`, and every element was zero from `L = 18` on. The synthetic "risky"
regression above only checked finiteness and symmetry, so it did not catch this.

The kernel now uses a running normalization. With `x_p = W_p chi(r_p)` and
`rho_p = (r_{p-1}/r_p)^L` it accumulates

- `A_p = rho_p A_{p-1} + x_p`, the prefix sum normalized by `r_p^L`;
- `B_p = rho_{p+1} (B_{p+1} + x_{p+1}/r_{p+1})`, the suffix sum normalized by `r_p^-L`;

and forms `inner_p = A_p / r_p + B_p`. This avoids artificial global-scale
underflow; genuine Float64 range and cancellation limitations remain.
The recurrence requires finite positive sorted points; it does not sort them.
The global shifts and log-magnitude recovery helpers are gone.
The tests in `test/radial/runtests.jl` compare the new kernel with a BigFloat
evaluation of the same quadrature formula for `L` up to 132.

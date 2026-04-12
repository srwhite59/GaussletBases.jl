# Radial Multipole Stabilization Milestone

This note records the narrow operator-side robustness milestone introduced by
commit `346c4c6` (`Stabilize radial multipoles and record prototype
milestone`).

This is not another prototype milestone. The prototype milestone was
`53a16a4`, which introduced the named cached paper-parity radial prototype.
This follow-on milestone is about the downstream radial multipole / `Vee`
builder in `src/operators.jl`.

## Problem before the patch

Before this patch, the integral-diagonal radial multipole builder used raw
power factors directly inside the prefix/suffix accumulation:

- `r^L`
- `r^-(L+1)`

That behavior was numerically fragile for larger multipole order `L` or wide
radial ranges. In those regimes, the raw-power path could lose stability and
produce non-finite intermediate or final data, allowing `Inf` or `NaN` values
to leak into downstream radial multipoles.

## Fix in `src/operators.jl`

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

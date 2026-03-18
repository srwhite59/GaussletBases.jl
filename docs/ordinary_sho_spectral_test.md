# Ordinary PGDG SHO Spectral Test

This note records the next validation step for the ordinary mapped branch.

## Why raw `||ΔT||` is not the whole story

The earlier one-body fidelity passes showed a real residual mismatch on the
localized ordinary PGDG route, and the largest part of that mismatch was on
the kinetic matrix.

But a different kinetic matrix is not automatically a worse basis.

If two localized basis constructions span nearly the same low-energy smooth
subspace, they can have visibly different matrix entries while still giving
almost the same low-energy physics.

So the right next question is not:

- are the raw one-body matrices numerically close entry-by-entry?

It is:

- does the experimental PGDG route have comparably good low-energy smooth
  completeness in the regime where it is supposed to work well?

## Why a low-energy SHO test is a better probe

The one-dimensional harmonic oscillator is a clean probe of that question.

It uses only:

- the kinetic operator
- a smooth low-order polynomial potential

So it tests exactly the kind of low-momentum / low-order smooth completeness
that should be well captured if the two bases have nearly the same practical
span.

This is therefore a better validation step than another raw matrix-norm pass.

## How this connects to the White–Lindsey standard

The historical standard here is not exact basis identity.

It is closer to:

- nearly identical span/subspace
- almost identical fitting behavior
- very similar low-energy results in the intended regime

That is exactly what the SHO test is meant to probe.

If the numerical mapped route and the experimental PGDG route give nearly the
same low-lying SHO spectrum in the friendly mild/core-supported regime, then
the remaining raw kinetic mismatch is mostly a representation difference, not a
serious basis-quality failure.

## What this pass should compare

The narrow comparison should use:

- the numerical mapped route as the reference
- the localized experimental PGDG route as the candidate
- the friendlier hybrid/core-supported regime as the main target
- one harder pure mapped case only as a stress-test comparison

For each case, the useful outputs are:

- the lowest few SHO eigenvalues
- the exact SHO values
- backend-to-backend spectral differences
- backend-to-exact errors

That is the right decision point before spending effort on a more expensive
kinetic-side refinement.

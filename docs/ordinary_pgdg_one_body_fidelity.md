# Ordinary PGDG One-Body Fidelity

This note records the current one-body decision point on the mapped ordinary
branch.

The previous localized-backend pass established two things clearly:

- overlap is effectively solved on the localized PGDG route
- the static ordinary Cartesian `Vee` assembly is already close to the
  numerical reference in the mild mapped regime

So the remaining one-body issue is not overlap, and it is not the separable
interaction assembly. The remaining issue is `H1`.

## What the last diagnostic pass proved

The mild mapped comparison is done in the same regime that motivated the
ordinary-backend pivot:

- one global Asinh map
- `count = 5`
- fixed outer centers at `x = ±6`
- mild distortion strength `s = 0.5`

In that regime, a numerically localized reference shows:

- overlap error at roundoff
- Gaussianized nuclear-factor errors already only at the `10^-6` to `10^-5`
  level
- a much larger kinetic discrepancy
- an assembled `H1` gap almost entirely explained by the kinetic side

So the diagnostic conclusion is clear:

**the remaining one-body discrepancy is kinetic-dominated.**

## Why the aligned-kinetic route is useful but not the implementation target

The previous pass also showed that if the numerically localized kinetic matrix
is aligned to the localized experimental basis, the `H1` gap narrows
substantially.

That is valuable because it tells us:

- the kinetic operator is really the piece that still matters
- the possible remaining improvement is not huge or mysterious
- overlap cleanup and localization were not the missing ingredients

But that aligned-kinetic route should remain only a diagnostic oracle.

It is **not** the right implementation target for the ordinary PGDG backend,
because it pulls the numerically localized mapped basis and its kinetic matrix
back into the backend construction itself.

That would blur the intended architectural split:

- numerical mapped ordinary path = trusted validation/reference route
- PGDG-style mapped ordinary path = analytic experimental implementation route

So the aligned-kinetic construction is kept only as an internal benchmark.

## The present analytic task

The next implementation target is therefore narrower and cleaner:

- keep the localized PGDG route analytic
- improve only its kinetic fidelity
- leave overlap solved
- leave the already-good `Vee` side alone

The refinement used here is a more derivative-aware one-Gaussian proxy for each
distorted primitive.

The pre-COMX log-fit proxy was already good for span and projector agreement.
But kinetic accuracy is more sensitive to derivative information than overlap
or Gaussianized Coulomb factors are.

So the new localized experimental backend uses:

1. the same global mapped ordinary basis scaffold
2. one Gaussian proxy per distorted primitive
3. a short local fit to the logarithmic derivative of the distorted primitive
4. the usual cleanup / orthogonalization / COMX-style localization stage

This stays within the same analytic primitive/contraction framework. It is not
a new broad subsystem and not a numerical-reference pullback.

## How success should be judged

The standard remains the White-Lindsey one:

- nearly identical span/subspace
- almost identical fitting behavior
- good one-body accuracy in the mild-to-moderate distortion regime

Exact basis identity is still not the target.

The practical comparison for this pass is therefore:

- numerical reference
- pure analytic localized PGDG backend before the new kinetic-focused proxy
- internal aligned-kinetic oracle
- improved pure analytic localized PGDG backend

with the main metrics:

- one-dimensional kinetic error
- assembled 3D `H1` error
- overlap error
- `Vee` agreement

## Interpretation

If the improved pure analytic localized backend narrows the kinetic and `H1`
gaps materially while keeping overlap solved and `Vee` close, then the
ordinary one-body branch is much closer to being solver-ready.

If it does not, then the next step should still be another narrow analytic
kinetic refinement, not a He solver.

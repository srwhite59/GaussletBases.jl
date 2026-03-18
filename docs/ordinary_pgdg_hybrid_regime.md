# Ordinary PGDG Hybrid Regime

This note records a directional choice for the ordinary PGDG branch.

## Why the radial branch should stay numerical for now

The radial branch is not the natural place to force a PGDG-style analytic
construction right now.

Its analytic structure is genuinely harder:

- the half-line / `\Theta(r)` support is part of the basis definition
- the near-origin behavior is more singular
- the mapped primitive integrals are less forgiving than on the full line

So the current radial recommendation stays the same:

- keep the radial branch numerical
- keep the explicit primitive/contraction architecture there
- do not try to force a radial PGDG route yet

## Why the ordinary Cartesian branch is the right analytic PGDG target

The full-line ordinary Cartesian branch is where PGDG-style analytic work makes
sense:

- the primitive layer can stay Gaussian
- the Coulomb expansion is separable
- the one-body and IDA pieces can be assembled from one-dimensional factors
- the analytic primitive/contraction route is realistic there

That is the branch where the current experimental PGDG backend should keep
evolving.

## Why the small-`c`, no-core pure mapped tests are harsher than the intended regime

The earlier pure mapped ordinary tests were useful, but they are harsher than
the practical White-Lindsey operating regime:

- they push toward smaller `c` / tighter near-origin distortion
- they use no explicit core Gaussian support
- they ask the mapped ordinary backbone to do all of the near-nuclear work by
  itself

That is better understood as a stress-test regime.

It is not the right regime for deciding whether the ordinary PGDG route is
practically viable.

## The next practical route

The friendlier and more historically faithful direction is:

- mild-to-moderate full-line mapping
- explicit centered core Gaussian augmentation
- then the same overlap cleanup / orthogonalization / COMX localization on the
  combined one-dimensional space

In other words:

- mapped ordinary gausslet backbone
- plus a few explicit core Gaussians
- with the numerical mapped route retained as the reference
- and the localized analytic PGDG route treated as the candidate

## What this pass is testing

The goal here is not a large basis-library system.

It is a small hybrid construction that can answer one narrower question:

**does the ordinary PGDG analytic route behave more like the practical
White-Lindsey picture when the basis is given a mild map and explicit core
Gaussian support?**

That is the right next benchmark before any ordinary-branch He-style solver
layer is opened.

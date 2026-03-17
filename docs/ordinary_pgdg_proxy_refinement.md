# Ordinary Mapped PGDG: Refining the Primitive Proxy

This note records the next step after the experimental COMX follow-on pass.

## 1. What the COMX pass established

The previous pass answered an important structural question.

It showed that:

- overlap cleanup and COMX/localization fix the representation-level overlap and
  position structure
- but they do **not** materially reduce the hydrogen energy gap

So COMX was worth testing, but it is not the main missing ingredient in the
ordinary mapped analytic path.

## 2. What the next task therefore is

The next task is to improve the **analytic primitive proxy itself**.

The current local-linear proxy is useful, but it is still only a first-order
physical-space approximation to the explicitly distorted primitive layer.

That means the main question is no longer:

- should we add more localization machinery?

It is now:

- can we choose a better plain-Gaussian proxy for the distorted primitives,
  while keeping the one-body operator algebra analytic?

## 3. What “success” should mean here

The historical White–Lindsey standard is not exact basis-function identity.

The stronger and more relevant claim is closer to:

- the DG and PGDG functions are very similar
- the fitting behavior is almost identical
- the spans or retained subspaces are nearly identical

So the right comparison target is:

- nearly identical span / projector behavior
- almost identical fitting quality

not:

- pointwise identity of the localized basis vectors

## 4. The refinement used in this pass

The present refinement is a **short weighted log-quadratic Gaussian fit**.

For each distorted primitive, instead of choosing the physical-space Gaussian
only from the first-derivative scaling of the map, the refined proxy fits

```text
log f(x) ≈ a0 + a1 (x - x0) + a2 (x - x0)^2
```

over a short physical-space window around the mapped primitive center `x0`,
using weights proportional to `f(x)^2`.

That gives:

- a shifted Gaussian center
- a corrected Gaussian width
- a fitted amplitude

This is still a one-Gaussian-per-primitive proxy. It does **not** yet add the
larger historical PGDG machinery.

This is still deliberately narrow:

- one-dimensional
- one-body only
- still experimental/internal
- still not the full historical PGDG driver path

## 5. What should be compared now

The right three-way comparison is still:

1. current mapped numerical ordinary path
2. current analytic proxy
3. refined analytic proxy

And the right diagnostics are:

- direct subspace/projector agreement
- fitting/projection agreement for plain Gaussians and odd `x`-Gaussian test
  functions as secondary probes
- hydrogen energy only as an end-to-end check

That is the right standard for deciding whether the ordinary branch should now
pivot more decisively toward the PGDG-style analytic one-body route.

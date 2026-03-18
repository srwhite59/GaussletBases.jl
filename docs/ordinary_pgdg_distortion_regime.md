> **Status:** supporting development note. For the current ordinary branch,
> read `docs/current_ordinary_branch.md` first. For this note chain, use
> `docs/ordinary_pgdg_supporting_notes.md`.

# Ordinary Mapped PGDG: Distortion Regimes

This note records the next decision point after the first refined analytic
primitive proxy.

## 1. Why distortion strength matters

The right DG-versus-PGDG comparison is not independent of the map.

If the coordinate distortion is mild, then a plain-Gaussian proxy for the
distorted primitive layer should have a better chance of recovering nearly the
same span after contraction. If the distortion becomes very strong, then one
plain Gaussian per distorted primitive is a much harsher approximation, and the
agreement should be expected to degrade.

So the relevant scientific question is not:

- does the proxy remain almost identical under any distortion at all?

It is:

- does the proxy behave in a White–Lindsey-like way in the mild-to-moderate
  regime that is actually useful for mapped ordinary bases?

## 2. What the comparison target should be

The historical target is not exact basis-vector identity.

The more faithful target is:

- exact restoration of overlap structure at the PGDG stage where that should
  happen
- approximate restoration of the center/position structure
- very similar functions
- almost identical fitting behavior
- nearly identical spans or retained subspaces

That means the most important diagnostics here are:

- principal-angle / subspace metrics
- projector differences
- fitting/projection errors for simple Gaussian and odd `x`-Gaussian probes
- hydrogen energy only as an end-to-end check

## 3. Why strong distortion should be treated differently

A deliberately strong map is useful as a stress test.

But it is not the right main success criterion for the present experimental
ordinary PGDG line. If the proxy works very well in a mild or moderate mapped
regime and then degrades under a much stronger distortion, that is scientifically
useful information, not a failure of the whole direction.

So this pass should be interpreted in three bands:

- mild / physically relevant
- moderate
- strong stress-test

Only after that separation is clear does it make sense to decide whether the
ordinary one-body branch is ready to lean more heavily on the analytic PGDG-style
route.

# Atomic Residual Injection Optional Refinement

This note records one possible future refinement for the atomic Cartesian
QW residual route. It is not the current required fix.

## Current Status

The current repo-side blocker on the large atomic nested lane was resolved by
hardening residual-block stabilization after the `:near_null_only` keep step.
That fix should remain the active contract unless a later case shows a real
need for something stronger.

So:

- current contract: keep residual modes with retained residual-overlap
  eigenvalue `> 1e-8`
- current fix: stabilize the kept residual block numerically before downstream
  center extraction and transforms
- current result: large full-rank retained residual blocks can build cleanly
  without aggressive pruning

## Optional Future Refinement

If later work shows that weaker-but-useful residual modes should be retained
without pushing the ordinary added-residual channel into instability, consider
an **injection** split.

The important contract is:

- injection is **replace**, not add
- the removed mode is one direction in the existing approximate span, usually a
  linear combination of many gausslet-like channels, not one single basis
  function
- the injected exact mode enters through the full post-contraction working
  span, not through a one-orbital swap heuristic

## Intended Atomic Translation

For a future atomic residual split, the retained post-null-cutoff residual
space could be divided into:

- added residual channels
- injected residual channels

with a policy like:

- `lambda <= tau_null`: discard as numerical null
- `tau_null < lambda < tau_add`: inject by replacement
- `lambda >= tau_add`: add as an ordinary residual channel

The injected channels would then:

- use exact post-contraction linear algebra for overlap and one-body matrices
- avoid the fragile residual-center / residual-width / effective-Gaussian path
- be treated in the approximate two-particle surrogate spirit already used by
  the residual route, rather than as fully explicit new Gaussian residual
  channels

## Why This Is Deferred

This design is attractive, but it is not currently required to make the large
atomic nested route work. The present hardening pass showed that the failing
large retained blocks were still full-rank and could be stabilized without
dropping or reclassifying modes.

So the future rule is:

- do not reach for injection first
- use it only if a later case shows that stabilization plus the near-null
  cutoff is no longer enough, or if physically useful weaker modes need a more
  careful treatment than ordinary residual addition

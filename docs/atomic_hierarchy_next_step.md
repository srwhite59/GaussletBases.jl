# Atomic Hierarchy Next-step Note

This pass chooses the first hierarchy-design step after the trusted atomic
nonrecursive anchors.

## Starting Point

Active atomic footing:

- `legacy_atomic_gaussian_supplement("He", "cc-pVTZ"; lmax = 0)`

Trusted nonrecursive atomic states:

- unnested hybrid QW reference
- shell-plus-core hybrid as the conservative nested anchor
- corrected complete-shell hybrid as the first trusted reduced anchor

This is enough to stop doing anchor-establishment work and decide what the
first post-nonrecursive hierarchy step should be.

## Three Plausible Directions

### 1. Stop at the current corrected complete-shell anchor for now

Pros:

- zero new source-language risk
- the current atomic comparison state is already trustworthy

Cons:

- it does not advance the hierarchy
- it leaves the remaining direct core block untouched
- it does not answer whether the current complete-shell source language can be
  continued hierarchically inside the retained core

Conclusion:

- not the right next implementation step
- only acceptable as a pause if the next atomic hierarchy move were still
  structurally unclear

### 2. Add one more flat nonrecursive shell-sequence refinement

This means extending the current complete-shell sequence by one more shell-like
replacement step without yet refining the retained core as its own hierarchical
piece.

Why this is not the right continuation:

- the trusted corrected complete-shell anchor is already a full four-layer
  shell decomposition of the `13^3` working cube down to a direct `5^3` core
- each current complete shell layer contributes the same source-language
  content:
  - faces: `6 * (4 * 3) = 72`
  - edges: `12 * 3 = 36`
  - corners: `8`
  - total per shell layer: `116` columns
- the current trusted reduced anchor is therefore
  - `589 = 4 * 116 + 125`
  - where the remaining direct block is exactly the `5^3 = 125` core

So a further flat shell-sequence move would not introduce a new useful
structural ingredient. It would only delay the real question:

- how should the remaining direct `5^3` core itself be nested?

It is also the more awkward design direction because it forces new shell policy
decisions at the smallest box immediately, but without yet making that box an
explicit hierarchical object in its own right.

Conclusion:

- not the best next step
- the marginal nonrecursive extension is less natural than making the core
  itself the first hierarchical target

### 3. Introduce the first hierarchical core-only refinement inside the trusted complete-shell working box

This means:

- keep the current trusted corrected complete-shell outer four-shell object
  unchanged
- keep the fixed-block consumer unchanged
- replace only the remaining direct `5^3` core by one local hierarchical
  object built from the original parent-space functions assigned to that core
  region
- do not re-coarsen already-renormalized outer shell functions

Why this is the right continuation:

- it targets the only remaining unreduced direct block in the trusted reduced
  anchor
- it reuses the existing complete shell language instead of inventing a new
  flat sequence policy
- it matches the practical W&L-style picture better than strong recursion:
  choose a nested box hierarchy, compress the assigned parent-space pieces
  level by level, and stop at a small core
- it fits the current code structure directly:
  - `_nested_shell_sequence_from_core_block(...)` already separates
    shell-layer assembly from the core block source
  - the current fixed-block adapter and QW nearest/GGT consumer do not need a
    rewrite
- it keeps the outer trusted anchor frozen while testing hierarchy locally in
  the smallest remaining block

This is the first point where a hierarchical core refinement becomes
structurally justified rather than premature.

## Recommendation

The next atomic hierarchy step should be:

- the first hierarchical core-only refinement inside the trusted corrected complete-shell anchor

In practical terms:

- keep the trusted corrected complete-shell outer sequence as-is
- act only on the retained direct `5^3` core
- introduce one local hierarchical replacement for that core, built from the
  parent-space functions assigned to the `5^3` region
- then compare against the current trusted corrected complete-shell anchor

This is better than one more flat nonrecursive shell step because it attacks
the actual remaining direct block instead of expanding the outer sequence again.

## What The Next Implementation Pass Should Actually Do

Keep it narrow:

1. Construct one local core-only nested object inside the current direct `5^3`
   core from the original parent-space functions assigned to that region.
2. Feed that new core object through the existing
   `_nested_shell_sequence_from_core_block(...)` route.
3. Keep the outer complete-shell layers unchanged.
4. Compare only on the established atomic He fixed-`a` family with the same
   shared `lmax = 0` supplement.

The diagnostic gates should be the same trusted-anchor gates:

- overlap quality
- finite positive transformed fixed-block weights
- projected fixed-only interaction transfer
- parent low-energy one-body capture
- final nearest/GGT `E1`
- final nearest/GGT `⟨Vee⟩`
- fixed dimension and runtime relative to the current `589`-dimensional
  corrected complete-shell anchor

## What Is Still Too Speculative

Not yet justified:

- a broad new flat shell-sequence policy beyond the trusted complete-shell line
- full multi-level recursion in the strong algebraic sense over the whole
  working cube
- diatomic placement/policy
- true physical `lmax = 1` supplement claims

Those should wait until the first hierarchical core-only refinement is either
trusted or clearly shown not to buy enough.

## When To Move On To Diatomics

The roadmap should stay atomic until one of these happens:

- the first hierarchical core-only atomic refinement is trusted and
  numerically useful
- or it is cleanly ruled out as not worth the added source-language complexity

That is the right point to stop atomic hierarchy work and move to diatomics on
top of a settled atomic comparison basis.

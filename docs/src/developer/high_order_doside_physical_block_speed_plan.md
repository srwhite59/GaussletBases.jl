# High-Order Doside Physical-Block Speed Plan

Date: 2026-05-01

Status: active experimental design note, not a production merge plan

## Purpose

This note records the current design discussion for making the experimental
high-order doside line practical on larger parent lattices.

It is not another validation note for the already-built experimental lane.
Instead, it is a forward plan for a faster construction path that preserves the
main scientific intent while avoiding the slow dense parent-lattice route that
currently limits the application study.

This note should be read together with:

- [high_order_doside_experimental_plan.md](high_order_doside_experimental_plan.md)
- [high_order_doside_distorted_parent_followup.md](high_order_doside_distorted_parent_followup.md)

## Long-Term Integration Goal

Although the current higher-order doside work is still an experimental line,
the design should be chosen so that it can fit naturally into the main repo if
it succeeds.

That means the target is not a permanently separate experimental mechanism.
The target is:

- a new higher-order nested-basis capability that could eventually live
  alongside the existing mainstream routes

So even while the implementation remains quarantined for now, the design should
prefer:

- concepts that match the repo's existing construction style
- operator flows that can share the main infrastructure naturally
- basis objects that could later coexist with ordinary routes without awkward
  one-off contracts

It should avoid, as much as practical:

- special-purpose ideas that only make sense inside a permanent side lane
- experimental-only abstractions that would become a burden if the route moves
  into normal use

## Current PGDG Contract Clarification

The recent timing and accuracy diagnostics force one important contract point to
be stated explicitly.

For this line of work:

- the established PGDG plus COMX construction is part of the repo baseline
  contract
- it should not be treated as an ad hoc replaceable stand-in for the distorted
  gausslet basis
- it is only expected to match the distorted gausslet basis approximately

So when the higher-order route is tested on a PGDG backend, the right question
is not:

- can one invent a different analytic local Gaussian fit that agrees better
  with the distorted gausslet basis?

The right question is:

- is the current higher-order route using the established PGDG plus COMX
  construction correctly and consistently?

This matters because some recent diagnostic experiments compared the current
baseline PGDG construction against alternative local Gaussian fits such as
log-fit and derivative-fit variants. Those experiments were useful only as
off-contract sensitivity probes. They do **not** justify changing the baseline
PGDG contract casually.

The working rule for this note is therefore:

- do not redesign PGDG itself as part of the higher-order speed work unless
  that decision is made explicitly and separately

## Current PGDG Diagnostic Interpretation

The present diagnostic picture is now fairly specific.

What appears to be working:

- primitive analytic Gaussian kinetic integrals are correct
- contraction of those primitive kinetic matrices into the PGDG basis is
  correct
- numerical-reference quadrature is stable on the ordinary Coulomb-Gaussian
  terms
- the new higher-order reduced path can use the PGDG-localized backend quickly

What is still problematic:

- on the small distorted higher-order test, the current PGDG-localized
  one-electron result still differs too much from the numerical distorted
  reference
- that gap shows up mainly through the kinetic side
- replacing only the PGDG-localized kinetic by the oracle/corrected kinetic
  almost removes the He+ energy miss on that test

What that means for the current plan:

- same-backend PGDG reduced-versus-direct agreement is the required algebra
  check for the high-order route
- numerical-reference differences should be interpreted as backend-fidelity
  diagnostics, not as reduced-route algebra failures
- the remaining problem should be treated first as a fidelity or usage issue in
  the higher-order consumption of the baseline PGDG plus COMX object
- it should **not** be treated first as a reason to redesign PGDG itself

The off-contract log-fit and derivative-fit experiments showed that the higher-
order result is sensitive to the local Gaussian representation, but they should
be read narrowly:

- they are evidence about sensitivity
- they are not evidence that the baseline PGDG contract should be replaced

So the next diagnostic priority is:

- audit the higher-order route against the intended PGDG plus COMX contract
- check whether the higher-order route is consuming the same basis and operator
  conventions as the established ordinary route
- keep alternative local-fit constructions only as diagnostics unless a
  separate decision is made to broaden the PGDG contract

## Current Stage-4.5 Timing Finding

The first coarse reduced-one-body timing probe has now been run.

The important result is that the transformed-shell mechanics are **not** the
current bottleneck.

On the reduced one-body path:

- building the transformed shell block is cheap
- contracting the one-dimensional operators is cheap
- assembling the reduced three-dimensional one-body operator is cheap

The dominant cost on the distorted path is instead:

- building the distorted one-dimensional parent one-body factors

The probe also confirmed that the routine distorted path no longer spends time
on the old dense parent-space reference comparison.

In short:

- the dense parent comparison was successfully removed from the routine
  distorted validation path
- but the present reduced route is still slowed mainly by rebuilding expensive
  distorted one-dimensional operator data

That means the next speed target is not shell selection and not reduced
three-dimensional assembly.

The next speed target is:

- reuse or restructure the distorted one-dimensional one-body/operator data so
  they are not rebuilt expensively inside the reduced higher-order path

The first reuse pass has now confirmed that caching this distorted
one-dimensional parent one-body package is very effective for repeated work on
the same parent:

- first distorted parent one-body build on the `5`-point test parent:
  about `148 s`
- second build on the same prepared axis data and same exponents:
  about `0.015 s`

So the expensive part is not "using the one-dimensional parent package" after
it exists. The expensive part is the first construction of that package.

This means the next fast route should be organized around:

- build the distorted one-dimensional parent package once per parent and
  exponent set
- reuse it across all sides, shells, and local higher-order reductions on that
  parent

In particular, the next optimization target remains:

- one-dimensional parent one-body construction and reuse

and not:

- shell-selection mechanics
- reduced three-dimensional assembly
- or redesign of the baseline PGDG contract

## Current Interpretation

The main scientific signal from the recent distorted-parent study is:

- the higher-order structured basis itself appears to be the main win over the
  current lower-order nesting route
- residual correction layers can help, but they do not appear to be the main
  source of the advantage
- the present practical limit is speed, not the basic idea

The fake-occupancy study against Gaussian target spaces was especially useful
for this interpretation.

In plain language:

- if one asks how much chemically relevant target-space content is still
  missing, the structured higher-order basis already improves materially over
  the lower-order route
- the extra refinement from residual functions looks more like cleanup than the
  primary source of the gain

That means the next major question is not "can the current dense
parent-projection path be pushed a little further?"

It is:

- how can the higher-order structured basis be built in the same kind of
  efficient reduced representation that makes the conventional nesting line
  practical?

## Main Speed Diagnosis

The current experimental line is slowed mainly by large dense parent-space work.

The important point is not just that the parent lattice is large. It is that
the present route repeatedly falls back to a representation where basis
construction and operator assembly are expressed directly in the full parent
lattice.

The conventional lower-order route in this repo is faster because it uses a
different pattern:

1. build local one-dimensional contraction data
2. combine them into three-dimensional tensor-product blocks
3. keep the basis in that reduced block representation as long as possible
4. only then apply smaller changes of basis or contractions on top

That is the intended speed model for the next higher-order implementation.

The key practical goal is therefore:

- escape the parent lattice early
- do most work in reduced block spaces
- apply any later mixing or rotation only after the basis dimension is already
  much smaller

## Core Reframing

The right primitive object is not "the shell" by itself.

The right primitive object is:

- a transformed three-dimensional cube or rectangular block

and that transformed block is built as:

- the tensor product of three transformed one-dimensional blocks

This gives a cleaner picture:

- the full transformed block is the primary local object
- a shell basis is obtained only by omitting inner transformed block functions
- the shell is therefore a selection inside the full transformed-block
  framework, not a separate construction rule

That point matters for both correctness and speed.

## Physical-Space Contract

The transform should be defined in physical space, not on the undistorted
reference lattice.

Reason:

- the distorted parent functions already encode strong nuclear-region behavior
- that physical scrunching is exactly what the transformed local basis is meant
  to capture
- building the transform first in the undistorted reference coordinate and then
  trying to repair the physical-space distortion afterward is not the intended
  contract

So the working contract should be:

- on each local one-dimensional distorted block, fit the distorted parent
  functions to ordinary physical-coordinate polynomials
- then form the three-dimensional block transform as the product of those three
  one-dimensional transforms

For modest local side lengths, the one-dimensional construction should simply
follow the same polynomial idea already used in the ordinary gausslet line:

- fit to `1, x, x^2, x^3, ...`

The note does not assume that this is the final answer for arbitrarily large
local side lengths. For larger `n_s`, there may be conditioning or shape issues
that require refinement. But for the present target regime, the ordinary
polynomial-fit view is the intended starting point.

## Primary Objects

This plan uses the following conceptual objects.

### Parent distorted basis

The original distorted Cartesian product basis on the working lattice.

This is the expensive object that the new route should leave as early as
possible.

### One-dimensional physical block transform

For a local interval on one coordinate axis, build a contraction matrix whose
columns reproduce ordinary physical-coordinate polynomials over that interval.

This is the one-dimensional higher-order primitive.

### Three-dimensional transformed block

For a local cube or rectangular block, form the tensor product of the three
one-dimensional transforms.

This gives the full local transformed block basis.

### Structured shell subset

Select only the outer transformed block functions.

This is the structured retained basis that plays the role of the higher-order
shell basis.

### Residual Gaussian supplement

Keep the existing idea that Gaussian-derived correction functions may be added
on top of the structured higher-order basis.

This note treats residual Gaussians as the main likely completeness-correction
layer. It does not assume that a separate residual-gausslet program is worth
carrying at the same time.

## Why Residual Gaussians Matter To The Design

This point changes the priorities.

If the final working route already intends to add Gaussian correction functions,
then the higher-order gausslet line does not need to repair every last local
completeness defect by itself.

The higher-order structured basis should aim first to do the main job well:

- provide a better reduced core basis than the current lower-order route

The Gaussian supplement can then handle part of the remaining chemistry-oriented
completeness correction.

This suggests a simpler architecture:

1. structured higher-order transformed blocks
2. shell selection as needed
3. Gaussian supplement and residual Gaussian correction on top

The current working assumption is therefore:

- do not spend early design effort on a large residual-gausslet correction
  layer unless later evidence shows that the Gaussian supplement is not enough

## Why The Later Rotations Are Not The Main Concern

The discussion so far suggests that later local rotations or cleanup transforms
are not the main performance risk.

Reason:

- once the basis has already been reduced from the parent lattice to block-size
  local spaces, those later transformations act on much smaller matrices
- even if they are only sparse or block-sparse rather than perfectly trivial,
  they should still be much cheaper than repeated dense parent-space projection

So the design priority is:

- optimize the early contraction into reduced block spaces

not:

- optimize every later local mixing step first

## Operator Strategy

The intended operator strategy is:

1. assemble the local contraction matrices first
2. use them to move operators from the parent representation into reduced block
   representations
3. combine or transform operators further only after the parent dimension has
   already been removed from the problem

This means the main algebra should look like:

- parent operators exist only as the source representation
- contraction matrices are the first-class data
- operator application in the higher-order line should mostly happen after
  contraction, not before

In plain terms:

- the construction should be "build contraction matrices first, then apply them
  to operators"
- not "build finished higher-order basis functions in the full parent space and
  project everything densely afterward"

That is the same general reason the conventional nesting route remains fast.

## Nesting Expectations

The desired nesting standard is:

- exact in the zero-distortion limit
- clean and well-conditioned for the distorted cases that actually matter
- any residual incompleteness small enough that Gaussian correction functions
  can repair it

It is already understood that distortion may spoil exact local span identities
that hold in the zero-distortion case.

That is acceptable if all of the following remain true:

- the structured higher-order basis still captures the main gain
- the remaining defects are not large or pathological
- the Gaussian supplement can recover the important missing chemistry

This is a more realistic target than demanding perfect distorted-case local
shell theorems at every step.

## Main Risks And Tricky Aspects

These are the points that should be kept in mind during implementation.

### 1. Conditioning of one-dimensional polynomial fits

For moderate local side lengths, ordinary polynomial fitting is the intended
construction.

But conditioning may worsen as local `n_s` grows.

Questions to watch:

- when do the one-dimensional moment fits become numerically fragile?
- does the distortion make some intervals much worse than others?
- is local re-centering or scaling needed before solving the fit?

### 2. Clean merging of neighboring distorted blocks

The speed plan assumes that neighboring transformed blocks can be assembled into
a usable global structured basis without excessive global cleanup.

Questions to watch:

- how much overlap-metric cleanup is actually needed?
- can cleanup remain local or shell-local?
- does distortion create long-range overlap pollution that breaks the cheap
  picture?

### 3. Reuse across side lengths

The practical route needs clean nesting from smaller to larger local side
lengths.

Questions to watch:

- can one-dimensional transformed blocks for one `n_s` be reused naturally
  inside larger `n_s` constructions?
- can shell increments be expressed by small local updates?
- does distortion destroy the simple transfer story?

### 4. Operator contraction boundary

The main win only materializes if operators are moved into reduced block spaces
early enough.

Questions to watch:

- which one-body pieces contract cleanly immediately?
- which two-body or interaction pieces still force parent-space work?
- can mixed gausslet/Gaussian couplings use the same reduced-block logic, or do
  they remain a separate expensive path?

### 5. Even and odd families

The previous experiments already showed that even and odd family behavior can
diverge in subtle ways.

Questions to watch:

- should the first fast implementation target only one parity family?
- do even-side blocks need special care in block-center conventions?
- is the same shell-selection rule numerically safe in both families?

### 6. Scope creep from residual corrections

There is a risk of blurring the main structured-basis problem together with
residual-correction design.

This note recommends against that.

The first fast implementation should answer:

- can the structured higher-order block-contraction basis itself be built
  efficiently and work well?

Residual Gaussians can then be layered on afterward.

## Recommended Working Contract

The current best working contract is:

1. Use physical-coordinate polynomial fitting on local one-dimensional distorted
   blocks.
2. Build full three-dimensional transformed blocks as tensor products of those
   one-dimensional transforms.
3. Treat the shell basis as a selected subset of those transformed block
   functions.
4. Escape the parent lattice as early as possible.
5. Carry operators through contraction matrices rather than by repeated dense
   parent projection.
6. Use residual Gaussian supplements as the main remaining completeness
   correction layer.

## Recommended Implementation Order

The next work should stay incremental.

### Phase 1: One-dimensional physical-block transform

Build and validate:

- one distorted one-dimensional interval
- one polynomial-fit contraction matrix
- diagnostics for overlap, conditioning, and polynomial reproduction

Questions to answer:

- does the one-dimensional fit behave well for the target local side lengths?
- how sensitive is it to distortion strength?

### Phase 2: Single three-dimensional full block

Build:

- one transformed three-dimensional full block from the tensor product of the
  three one-dimensional transforms

Validate:

- dimension
- overlap quality
- comparison against the current higher-order experimental full-block span on a
  small case

### Phase 3: Shell selection inside the transformed block

Implement:

- shell selection by omitting inner transformed block functions

Validate:

- shell count
- shell-block overlap quality
- relationship to the current experimental structured basis on a small case

### Phase 4: Reduced operator contraction

Implement:

- contraction matrices first
- one-body operator contraction into the reduced block basis
- the smallest practical interaction path that avoids parent-space fallback
  where possible

Measure:

- build time split in plain language
- parent-space work that remains
- memory footprint

### Phase 5: Structured higher-order plus Gaussian supplement

Layer on:

- the Gaussian supplement / residual Gaussian logic already used in the hybrid
  routes

Validate:

- He+ behavior on a small distorted-parent case
- fake-occupancy capture against Gaussian target bases
- timing comparison against the current lower-order route

### Phase 6: Larger parent lattices

Only after the reduced-block path is working should the larger parent study be
repeated.

The purpose of the speed work is exactly to make that larger-parent comparison
practical.

## How To Judge Success

The fast higher-order line should be considered successful if it achieves most
of the following:

- retains the main structured-basis advantage seen in the distorted-parent
  studies
- reduces reliance on dense parent-space projection substantially
- makes larger parent lattices practical enough for meaningful He+ comparison
- keeps the mathematical contract understandable
- works naturally with the Gaussian supplement layer

It does not need, in the first pass, to:

- eliminate every small distorted-case incompleteness
- provide a residual-gausslet correction story
- become a production public API

## What To Keep In Mind During Tricky Local Decisions

When detailed implementation choices become distracting, return to these main
points:

- the full transformed block is the primitive object
- the shell is only a selection inside that object
- the transform is defined in physical space
- the real speed win comes from escaping the parent lattice early
- later local rotations are secondary
- the structured higher-order basis is the main gain
- Gaussian residual functions are the expected cleanup layer

## Bottom Line

The present dense parent-projection route has already done its job as a proof
and diagnostic lane.

The next meaningful step is not another small patch to that route.

It is a new higher-order implementation shaped like the conventional fast
nested line:

- one-dimensional local contractions first
- three-dimensional tensor-product blocks second
- shell selection inside those transformed blocks
- operator contraction after the basis has already been reduced
- Gaussian residual correction on top if needed

That is the most coherent path currently visible for making the higher-order
line scientifically useful on the larger distorted-parent lattices that still
matter, while still keeping open the path for it to become a mainstream repo
capability later.

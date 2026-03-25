# Diatomic Endcap Slab Design Note

This note records the current design direction for the next diatomic
box/refinement question after the midpoint-slab correction and the move of the
odd-on-symmetric rule to the actual `doside` boundary.

It is a design note for repo-manager, not yet an algorithm page.

## Current Landed Context

What is already landed and should be treated as current baseline:

- bond-aligned diatomic distortion is real in the ordinary QW line
- midpoint-slab preservation for odd-length homonuclear bond-axis splits is
  real in code
- the odd-on-symmetric retained-count rule now lives at the actual 1D
  `doside` / `COMX` boundary
- the ordinary and nested diatomic H2 lines are physically active, with a real
  molecular `s/p` supplement and nonzero residual-Gaussian sector
- the current lightweight visualization/debug path is good enough to inspect
  source geometry, compressed centers, and `doside` traces directly

Relevant current policy/code references:

- [cartesian_nested_diatomic_box_policy.md](/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/docs/src/algorithms/cartesian_nested_diatomic_box_policy.md)
- [diatomic_h2_completion_state.md](/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/docs/diatomic_h2_completion_state.md)
- [cartesian_nested_faces.jl](/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/src/cartesian_nested_faces.jl)

## Structural Clarification Already Made

One important reframing from the recent work should be recorded explicitly:

- the symmetry-preserving odd-retain rule is not best understood as a
  shared-shell special case
- it is better understood as a rule on the actual 1D contraction boundary:
  whenever a local interval passed to `doside` is symmetric about zero, an
  even provisional retained count should be reduced by one before the `COMX`
  localization

That rule now matches the structure of the code:

- `_nested_doside_retained_count(...)`
- `_nested_doside_1d(...)`

This means:

- shared-shell symmetric tangential contractions are covered
- symmetric child-shell tangential contractions are covered
- symmetric edge contractions are covered when they pass through the same 1D
  `doside` path

What is still *not* implied:

- this is not yet a full adaptive local retain-count policy
- it does not settle heteronuclear behavior
- it does not settle chain behavior
- it does not decide when a whole slab or region should exist as a separate
  object in the first place

## Visualization/Interpretation Clarification

Recent inspection work also clarified what should and should not be used as
human-facing geometry language.

For human-facing summaries and plots:

- avoid relying on `parent` / `child` terminology
- that language is only locally meaningful for one split and does not scale
  well to multisplit hierarchies
- avoid relying on detailed region labels as the main interpretation surface
  if the picture already shows the geometry

Preferred user-facing summary language:

- `final gausslet count / initial gausslet count`
- then separately `+ N residual Gaussians` when present

This is more stable than:

- `parent`
- `child`
- or an over-specific region vocabulary that may not survive deeper hierarchy

The main role of the visualization layer is therefore:

- show the geometry by eye
- show whether the source geometry and compressed hybrid geometry look
  sensible
- support count-based summaries, not hierarchy jargon

## New Design Question: Endcap Slabs

The next candidate refinement idea is not another midpoint policy change.

It is a distinct idea:

- allow an outer `nx × ny × 1` slab near the end of a child box to be
  contracted if it is visibly over-resolved relative to its local physical
  position

Example current motivating picture:

- after the midpoint-slab correction, the central symmetry issue is fixed
- in the larger `H2` debug case, some outer slab-like regions still look finer
  than needed
- one natural example is a `7 × 7 × 1` slab that might plausibly become
  `5 × 5 × 1`
- the remaining core behind it might then be something like `7 × 7 × 4`

This is not the same as the midpoint slab:

- the midpoint slab is required by homonuclear symmetry for odd-length split
  intervals
- an endcap slab would be optional
- its purpose is basis reduction subject to local completeness

So the governing question is:

- when is a slab fine enough relative to its position that it should be
  contracted?

## Why A Pure Geometric Rule Is Not Enough

The current conclusion from the recent discussions is:

- a raw box-size rule such as `7 -> 5` is not enough by itself
- the correct decision should depend on the local physical spacing implied by
  the distortion
- the right criterion should therefore be phrased as a local angular or
  spacing-resolution test

That leads to the next refinement:

- do not define the target resolution from `1 / nside` analytically
- define it from the actual localized `doside` geometry produced by the
  retained basis

## Empirical Angular Calibration For `nside`

The key idea is:

- `nside = 5` should correspond to the angular resolution actually produced by
  a representative `doside` contraction that keeps `5` localized points
- this is not something known exactly in advance from a simple formula

For atoms, a reference contraction such as:

- a symmetric `7 -> 5` face-side contraction

defines an empirical angular scale for `nside = 5`.

That scale can then be measured as follows:

1. Choose a representative symmetric local interval.
2. Run the tentative 1D `doside` contraction with retained count `5`.
3. Measure the actual localized center spacing returned by the `COMX`
   localization.
4. Convert that spacing to an angular spacing at the relevant radius.

So the target angular resolution is not:

- `1/5`
- or `asin(1/5)`
- or any other purely analytic placeholder

Instead it is:

- the empirical angular resolution produced by the actual retained/localized
  basis for the chosen reference case

## Candidate Endcap-Slab Acceptance Rule

The current best design shape is:

1. Identify a candidate endcap slab.
2. Propose a tentative contraction, for example `7 × 7 × 1 -> 5 × 5 × 1`.
3. For each tangential direction, make a tentative `doside` call using the
   proposed retained count.
4. Read off the actual localized center spacings from the trial result.
5. Convert those spacings to a local angular spacing using the slab position
   relative to the nearest nucleus.
6. Accept the slab contraction only if the realized angular spacing is not
   coarser than the calibrated `nside = 5` reference scale, up to a chosen
   tolerance band.

In short:

- tentative `doside`
- measure actual spacing
- reject if too coarse

This is better than a fixed geometry rule because it:

- respects the actual distortion
- uses the real localized basis geometry
- naturally inherits the odd-on-symmetric rule already built into `doside`

## Multiple Endcaps

Another point that should be recorded explicitly:

- the goal is not to refine exactly one endcap
- the goal is to make the basis as small as possible without losing
  completeness

So an endcap-slab rule should be iterative in principle:

1. test the outermost candidate slab
2. contract it if admissible
3. examine the next slab
4. repeat until the next trial would become too coarse

That means the intended policy is greedy and local:

- refine slabs when you can
- stop when the tentative localized spacing crosses the admissible bound

This is a much better long-term framing than:

- a one-off hand-tuned slab exception
- or a hard-coded “always peel one slab” rule

## Continuity Guard

One additional guard from the recent discussion is worth recording:

- the slab should not be contracted below the effective transverse resolution
  already used on the adjacent retained structure

So the acceptance test should likely compare against both:

- an empirical `nside`-based angular target
- a continuity requirement with the adjacent core/face structure

This avoids a case where the slab passes a crude absolute threshold but makes
an abrupt resolution jump relative to its neighbor.

## What The Current H2 Evidence Suggests

The recent H2 debug sequence supports the following interpretation:

- the earlier missing-center problem was real and was correctly fixed by moving
  the odd-on-symmetric rule to the 1D `doside` boundary
- after that fix, the remaining outer-visual sparsity is no longer a midpoint
  symmetry problem
- it is now reasonable to ask whether some slab-like regions are simply
  over-resolved and should be contracted

That is exactly why endcap slabs should now be considered.

But the evidence does **not** yet justify:

- a final slab-refinement rule
- a general adaptive retain-count policy
- a heteronuclear extension
- a chain extension

Those still need a separate diagnostic/calibration pass.

## Recommended Next Diagnostic Pass

Before implementing endcap slabs, the next pass should be diagnostic-only.

It should do two things.

### 1. Atomic `nside = 5` calibration

Measure the empirical angular resolution corresponding to `nside = 5` for a
small family of representative symmetric atomic contractions.

Suggested outputs:

- local interval
- localized centers
- central spacing
- radius from nucleus
- inferred angular spacing

The point is to determine whether `nside = 5` corresponds to a reasonably
stable angular-resolution band across several representative cases.

### 2. Candidate H2 endcap-slab trials

For the current bond-aligned `H2` line:

- identify one or more outer `7 × 7 × 1` slab candidates
- tentatively contract them to `5 × 5 × 1`
- measure the realized localized spacings and inferred angles
- compare those against the atomic `nside = 5` reference band
- decide whether one slab, or more than one slab, is admissible

The output of that pass should be sufficient to answer:

- whether endcap slabs are justified at all
- whether more than one slab is admissible
- whether the next implementation should be a narrow slab-refinement code pass

## Algorithm-Page Decision

This should stay a design note for now.

Recommended documentation decision:

- do **not** create a new algorithm page yet
- do **not** expand the diatomic box-policy page with a detailed endcap-slab
  algorithm before the rule is calibrated and landed

Why:

- the midpoint-slab rule is already stable enough to live in the algorithm page
- the endcap-slab rule is not yet stable enough
- the correct acceptance criterion still depends on a calibration study

What should happen later:

- once the endcap-slab rule is implemented and validated, fold the final rule
  into the existing diatomic box-policy algorithm page rather than creating a
  separate standalone page unless the slab logic becomes large enough to justify
  one

## Repo-Manager Recommendation

The next repo-manager-level decision should be:

- approve one narrow diagnostic/calibration pass for empirical `nside`
  angular resolution and tentative endcap-slab trials
- defer algorithm-page promotion until the rule is actually calibrated and
  landed

That keeps the current diatomic line disciplined:

- midpoint symmetry already corrected
- odd-on-symmetric `doside` rule already moved to the right structural
  boundary
- endcap slabs treated as the next design question, not prematurely as a fixed
  policy

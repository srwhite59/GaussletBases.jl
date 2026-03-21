# PGDG Refinement Bundle Next Step

This note records the next narrow integration step after landing the first
analytic ternary PGDG refinement mask.

The first mask is now real repo infrastructure:

- analytic ternary `1 -> 1/3` mask
- default `rho = 1.2`
- local one-step and repeated application helpers

The next step should **not** be a full hierarchy rollout. The right next move
is narrower:

- let one shared one-dimensional PGDG proxy bundle optionally consume that
  ternary mask as a refined auxiliary line

## Intended role

The immediate target is:

- refined auxiliary integral-evaluation lines such as `1/9` and `1/27`

while leaving the current coarse `1/3` PGDG basis-realization story separate.

So this next step should not reinterpret the base `1/3` path. It should add an
optional refined proxy line that can be used for one-dimensional Coulomb data
construction while the working basis remains what it already is.

## Proposed integration point

The cleanest shared integration point is the existing one-dimensional mapped
bundle layer in:

- `src/ordinary_mapped_backends.jl`

Specifically, the likely target is the shared internal constructor around:

- `_mapped_ordinary_gausslet_1d_bundle(...)`

That is the right place because it already feeds:

- ordinary Cartesian IDA
- Qiu–White reference work
- and related mapped ordinary constructions

So one narrow extension there can benefit multiple ordinary consumers without
wiring the full hierarchy everywhere.

## Proposed narrow API shape

The next step should stay internal. A reasonable first shape would be one
optional refinement setting on the shared bundle builder, for example:

- `refinement_levels = 0` for the current base `1/3` auxiliary line
- `refinement_levels = 1` for an auxiliary `1/9` line
- `refinement_levels = 2` for an auxiliary `1/27` line

or an equivalent internal object carrying:

- base proxy line
- stored ternary mask
- number of repeated local refinement applications

The main point is:

- the refined line should be derived from the stored analytic mask
- not rebuilt as a separate unrelated proxy construction

## What the next step should do

For the first narrow integration pass, the shared bundle should be able to:

1. start from the current base `1/3` PGDG Gaussian proxy line
2. apply the stored ternary mask zero, one, or two times
3. build the refined auxiliary Gaussian coefficient line
4. use that refined line in the one-dimensional Coulomb-data construction

The intended near-term use is the one-dimensional auxiliary integral layer:

- overlap / kinetic / position / `x^2`
- Gaussian-factor terms
- pair-factor terms

not a full rebasing of every ordinary path.

## What should stay separate for now

This next step should still leave open:

- whether coarse `1/3` PGDG is used as a basis-realization layer
- whether finer `1/9`, `1/27`, ... levels are used only as auxiliary
  integral-evaluation representations
- how far the hierarchy eventually propagates into all ordinary/Qiu–White
  routes

Those are later integration questions.

## Bottom line

The next functional step should be:

- one shared 1D PGDG proxy bundle optionally consuming the analytic ternary
  refinement mask to produce refined auxiliary lines

That is narrow enough to keep the hierarchy rollout disciplined, but useful
enough to expose the first real `1/9` / `1/27` integration point.

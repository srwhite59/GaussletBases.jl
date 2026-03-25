# Angular Research Track

This page records the smallest current repo-side statement of the angular
gausslet line.

It is an active research track, not a frozen public workflow branch.

## Current interpretation

The present angular direction in `GaussletBases.jl` should be read as:

- manuscript-facing and experimental
- centered on shell-local injected angular basis work that has not yet been
  imported into the package
- intentionally narrower than the mature radial/atomic producer line

The repo now contains only the first controlled scaffold needed to start that
import:

- the existing atomic `(l,m)` and Gaunt/sectorized IDA foundation
- a read-only curated sphere-point-set access layer

It does **not** yet contain:

- shell-local injected angular basis construction
- the full `sphgatomps`-style workflow
- optimizer/sweep drivers for sphere-point generation
- manuscript figure scripts

## Near-term target

The near-term scientific target is an **atomic angular benchmark ladder**:

1. HF
2. small ED
3. one DMRG-facing bridge

That is the first clean way to prove the angular line through the same
producer/consumer boundary now established for the current atomic IDA branch.

## Curated point-set scaffold

The current repo-side angular scaffold exposes a small curated subset of
optimized sphere-point sets through:

- `curated_sphere_point_set_orders()`
- `curated_sphere_point_set(order)`

These helpers are experimental and read-only. They are there so later
shell-local basis construction can start from stable curated point-set data
without importing the optimization workflow itself.

The vendored subset currently lives in:

- `data/angular/curated_sphere_points.toml`

and carries:

- point-set cardinality/order
- Cartesian coordinates on `S^2`
- nearest-neighbor spread ratio `nn_ratio`
- source tag / source project / source artifact note

## What is deferred

The following items remain explicitly deferred:

- Hooke as its own later dedicated line
- heteronuclear/angular coupling beyond the first atomic benchmark ladder
- any claim that the present angular data API is stable enough to freeze

Hooke remains important, but it should come later as its own workflow/paper
line after the first atomic angular benchmark ladder is real.

## Next import boundary

If this scaffold stays sound, the next real scientific import should be:

- shell-local injected angular basis construction

That is the first step that turns the current angular placeholder into a live
scientific branch inside `GaussletBases`.

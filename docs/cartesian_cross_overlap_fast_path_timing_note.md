# Cartesian Cross-Overlap Fast-Path Timing Note

This note records the first concrete timing milestone for the public Cartesian
cross-overlap / projector / bundle path.

The main point is narrow:

- the public exact `cross_overlap(...)` path is now much faster on the
  supported atomic hybrid lane
- the mathematical contract did not change
- the fast path now survives bundle export/import cleanly enough to remain fast
  on disk

This is a runtime/producer improvement, not a baseline or numerical-contract
change.

## Scope

The optimization covered the public Cartesian overlap stack:

- representation-level exact `cross_overlap(...)`
- projector / orbital transfer built on that overlap
- bundle export/import carrying the overlap sidecars needed to preserve the
  cheap path on disk

The relevant implementation lives in:

- `src/cartesian_basis_representation.jl`
- `src/cartesian_bundle_export.jl`
- `src/cartesian_bundle_io.jl`

The important policy point is unchanged:

- final working-basis transfer uses exact cross overlap `S_BA`
- final-basis self-overlaps are diagnostic only
- dense-reference tests still check the exact result rather than accepting a
  looser approximation

## Benchmark Case

The pinned case is:

- atomic hybrid He
- `ns = 5`
- `parent_count = 13`
- supplement `cc-pVTZ`
- `lmax = 1`
- source final dimension `526`
- target final dimension `428`
- one occupied column transferred

This is not meant to be the largest production case. It is a stable small
reference point that exercises the exact hybrid overlap path honestly.

## Before / After Summary

The decisive measured changes were:

- in-memory exact `S_BA` build:
  - before: about `28.32 s`
  - after: about `0.195 s`
- disk exact `S_BA` build:
  - before: about `18.92 s`
  - after: about `0.212 s`

The final multiply was never the bottleneck:

- `S_BA * C_occ` after the cleanup is about `3e-4 s`

So the practical workflow cost of:

- reading two bundles
- building exact `S_BA`
- transferring one occupied column

is now on the order of a few tenths of a second on this lane, not tens of
seconds.

## What Was Actually Fixed

The large slowdown was not mainly in the 3D factorized algebra itself.

The key missed fast path was in the 1D overlap layer:

- `_basis_cross_overlap_1d(...)` was missing the stored-overlap shortcut
  because semantically identical axes were not matching the old shortcut test
- that forced expensive numerical rebuilds of distorted
  `Distorted{Gaussian,AsinhMapping}` primitive cross overlaps
- the atomic lane was also redundantly rebuilding identical `x` and `y` axis
  work
- bundle import discarded the stored overlap sidecars, so the disk path still
  fell back to the expensive numerical route even after the in-memory fix

The fix set therefore did three things:

- let the 1D same-basis path hit the stored overlap route on this lane
- reuse identical axis work where the atomic geometry really permits it
- round-trip the relevant 1D overlap sidecars through the bundle layer

## Exactness Was Preserved

This timing milestone did not relax the trusted answer.

The tests still treat the dense reference as the referee. The reading is:

- same exact answer
- much faster path
- same public API surface

Representative validation on the pinned lane stayed clean:

- projector agreement: exact to roundoff on the tested path
- transfer agreement against explicit `S_BA * C_A`: exact to roundoff
- disk projector / transfer agreement against memory: exact to roundoff on the
  tested path

## Where The Time Went Before

Before the 1D fix, the measured in-memory `S_BA` build on the pinned case was
roughly:

- total: `28.32 s`
- `x` axis `_basis_cross_overlap_1d(...)`: `9.34 s`
- `y` axis `_basis_cross_overlap_1d(...)`: `9.36 s`
- `z` axis `_basis_cross_overlap_1d(...)`: `9.41 s`
- surrounding dense contractions: microseconds

So the old bottleneck was unequivocally the 1D overlap rebuild, not the final
matrix multiply and not the outer 3D contractions.

## Where The Time Goes Now

After the 1D fast-path and bundle-sidecar fixes, the same in-memory case is
roughly:

- total overlap build in transfer benchmark: `0.195 s`
- axis-table stage: about `3.45e-4 s`
- `x` axis basis-cross: about `5.5e-5 s`
- `y` axis: reused from `x`
- `z` axis basis-cross: about `3.6e-5 s`

The disk path is now similar:

- bundle read: about `0.036 s`
- overlap build: about `0.212 s`

So the 1D overlap layer is no longer the dominant cost on this lane.

## Next Hotspot

After this milestone, the next measured hotspot moved to the 3D
Cartesian-Cartesian block construction, at roughly:

- `0.18 - 0.20 s`

On the pinned He `ns = 5` lane, that remaining cost is small enough that it is
probably not worth immediate further optimization by itself.

If later larger production cases show that block becoming annoying again, that
is the next honest optimization target.

## Reading Of The Milestone

The correct reading of this work is:

- the public Cartesian overlap / projector / bundle stack reached a real
  runtime milestone
- the same exact answer is now produced much faster
- the persistence layer now preserves the fast path instead of defeating it
- the next bottleneck is identified clearly rather than guessed

That is the point this note is meant to pin.

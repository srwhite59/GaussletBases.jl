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

The repo now contains the first controlled scaffold needed to start that
import, plus the first shell-local experimental construction on top of it:

- the existing atomic `(l,m)` and Gaunt/sectorized IDA foundation
- a read-only curated sphere-point-set access layer
- a shell-local injected angular basis constructor on one curated sphere-point
  set
- a first shell-to-atom angular assembly layer over shell radii
- a first atom-side one-electron angular benchmark path built on that assembly

It does **not** yet contain:

- the full `sphgatomps`-style workflow
- optimizer/sweep drivers for sphere-point generation
- manuscript figure scripts
- the later HF / small-ED / DMRG-facing ladder stages

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
without importing the optimization workflow itself. The first shell-local basis
layer is now available through:

- `build_shell_local_injected_angular_basis(order; ...)`
- `build_shell_local_injected_angular_basis(point_set; ...)`
- `shell_local_injected_angular_diagnostics(shell)`

This remains an experimental research-track surface, not a frozen public API.

The first shell-to-atom assembly layer is now available through:

- `assign_atomic_angular_shell_orders(shell_radii; ...)`
- `build_atomic_shell_local_angular_assembly(shell_radii; ...)`
- `atomic_shell_local_angular_diagnostics(assembly)`

This assembly layer is still angular-only. It stops before Coulomb assembly,
HF/ED workflow wiring, or end-to-end atomic angular drivers.

The first narrow atom-side benchmark path is now available through:

- `build_atomic_injected_angular_one_body_benchmark(radial_ops; ...)`
- `atomic_injected_angular_one_body_diagnostics(benchmark)`

This benchmark stays on the one-electron central-potential line. It is there to
prove that the shell-local angular assembly supports a real atom-side Galerkin
calculation before the later benchmark ladder stages are imported.

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

- the next atomic angular benchmark ladder step on top of the first one-electron
  assembly benchmark

That is the next step that turns the current one-electron benchmark branch into
the later HF / small-ED / bridge-facing ladder inside `GaussletBases`.

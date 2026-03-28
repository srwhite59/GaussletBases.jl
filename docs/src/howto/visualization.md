# Visualization Utilities

`GaussletBases.jl` includes a small `viz/` subtree for inspecting emitted
point-cloud and path datasets.

These utilities are still lightweight and format-driven. They are not a full
interactive visualization framework. But they are useful for quickly checking:

- bond-aligned 2D projections
- bond-aligned 3D point clouds
- 3D ordering paths or traversal overlays

## What lives in `viz/`

The current user-facing entry points are:

- `viz/showpoints2d.jl`
- `viz/showpoints3d.jl`
- `viz/showpath3d.jl`

The 2D path uses `CairoMakie` and writes static figures.
The 3D viewers use `GLMakie` for interactive inspection.

All three scripts support the same simple activation pattern from the repo
checkout:

```bash
julia --project=viz -e 'using Pkg; Pkg.instantiate()'
julia --project=viz viz/showpoints2d.jl INPUT.dat OUTPUT.png "Optional title"
julia --project=viz viz/showpoints3d.jl INPUT.dat "Optional title"
julia --project=viz viz/showpath3d.jl INPUT.dat "Optional title"
```

For a metadata-only read without opening a viewer window:

```bash
julia --project=viz viz/showpoints3d.jl --describe INPUT.dat
julia --project=viz viz/showpath3d.jl --describe INPUT.dat
```

## Expected input style

The current viewers are intentionally simple and operate on emitted dataset
files rather than on internal package objects.

Typical inputs include:

- 2D projected point datasets with `# dataset ...` headers
- 3D point datasets with `# box ...` and optional `# path ...` blocks
- concatenated path blocks for ordering or traversal inspection

These are the same lightweight file styles used by current ordinary/nested
geometry helpers and paper-side emitters.

## When to use them

The viewers are most useful when:

- checking geometric locality or extent after basis emission
- inspecting nested or Qiu-White point arrangements
- looking at bridge/order traversals without writing custom plotting code

They are intentionally not required for the main package workflows, but they
are a useful inspection layer for users who want to see the emitted geometry
rather than only read diagnostics.

## Notes

- The 3D viewers depend on the `viz/` environment, which includes `GLMakie`.
- The 3D viewers require a working GUI/OpenGL path on the local machine.
- The `showpoints2d.jl` path remains the easiest route for producing a static
  figure that can be dropped into notes or reports.
- The viewers are still rudimentary, so they are linked lightly from the docs
  rather than treated as a primary workflow surface.

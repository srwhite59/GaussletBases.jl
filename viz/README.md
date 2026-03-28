# Visualization Utilities

This directory contains lightweight viewers for emitted point-cloud and path
datasets.

Current entry points:

- `showpoints2d.jl`
- `showpoints3d.jl`
- `showpath3d.jl`

From the repository root:

```bash
julia --project=viz -e 'using Pkg; Pkg.instantiate()'
julia --project=viz viz/showpoints2d.jl INPUT.dat OUTPUT.png "Optional title"
julia --project=viz viz/showpoints3d.jl INPUT.dat "Optional title"
julia --project=viz viz/showpath3d.jl INPUT.dat "Optional title"
```

For metadata-only inspection without opening a viewer:

```bash
julia --project=viz viz/showpoints3d.jl --describe INPUT.dat
julia --project=viz viz/showpath3d.jl --describe INPUT.dat
```

Notes:

- `showpoints2d.jl` uses `CairoMakie` and writes a static image.
- `showpoints3d.jl` and `showpath3d.jl` use `GLMakie` for interactive 3D
  viewing.
- The input formats are simple emitted dataset files rather than package
  objects. See the script sources for the supported header conventions.

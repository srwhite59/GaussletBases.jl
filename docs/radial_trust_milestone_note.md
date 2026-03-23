# Radial Extent And Interval-sampled Build Milestone

This note records commit `ae7ed92` (`Land radial extent and interval-sampled build milestone`).

Milestone name, in one sentence:

- this milestone lands the radial line in a coherent and trusted state after
  cleaning up its extent semantics and replacing the old dense setup-grid build
  with a much faster interval-sampled construction path

Old public/runtime problem, in one sentence:

- the radial public path mixed the user-facing meaning of `rmax` with internal
  numerical extents, while the old dense setup-grid construction made
  `build_basis(spec)` so slow that the standard radial workflow remained a live
  stabilization problem

New extent contract, in one sentence:

- `rmax` is the public user-facing last-center extent, `build_umax` and
  `quadrature_umax` are internal library-owned extents, and `quadrature_rmax`
  remains only an expert override

Construction-side change, in one sentence:

- the runtime family tables are now trimmed to machine-significant tails while
  the full high-precision tables are preserved internally, and the old dense
  setup-grid build has been replaced by interval-sampled construction without
  changing the public `RadialBasis` form

The supporting notes for this milestone are:

- [radial_extent_policy_note.md](/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/docs/radial_extent_policy_note.md)
- [radial_default_path_cache_study.md](/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/docs/radial_default_path_cache_study.md)

The trust sweep script is:

- [radial_trust_sweep.jl](/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/tmp/work/radial_trust_sweep.jl)

## Numerical Trust Sweep

### Standard Front-door Hydrogen

| Quantity | Value |
| --- | ---: |
| basis length | `35` |
| `diag.overlap_error` | `8.797053967946056e-6` |
| `diag.D` | `1.3554300216193975e-7` |
| `E0` | `-0.4999999972288096` |

Recipe:

- `:G10`
- `Z = 1`
- `s = 0.2`
- `c = s / (2Z)`
- `rmax = 30`
- default two-`xgaussian` front-door supplement

### `(l,m)` One-body Hydrogen, `lmax = 2`

| Lowest energies | Value |
| --- | ---: |
| lowest energy | `-0.4999999945176998` |
| lowest `l = 1` energy | `-0.12499999988503978` |
| lowest `l = 2` energy | `-0.05555046143346199` |

This keeps the expected hydrogenic pattern through the current atomic one-body
route.

### Smaller Nondefault Radial Recipe

| Quantity | Value |
| --- | ---: |
| basis length | `21` |
| `diag.overlap_error` | `4.0596245811904994e-7` |
| `diag.D` | `2.4278240367355398e-4` |
| `E0` | `-0.49999997523710266` |

Recipe:

- `:G10`
- `Z = 1`
- `s = 0.3`
- `c = s / (2Z)`
- `rmax = 12`
- no `xgaussians`

### Runtime Table Trim Versus Preserved Internal Tables

| Quantity | Value |
| --- | ---: |
| trimmed runtime `:G10` radius | `75` |
| preserved internal high-precision `:G10` radius | `132` |
| max sampled value difference after trimming | `4.440892098500626e-16` |

So the shorter runtime tables are materially cleaner for runtime support logic
without changing representative basis values at a scientifically meaningful
level.

## Timing

Standard timing recipe:

- `:G10`
- `Z = 2`
- `s = 0.2`
- `c = s / (2Z)`
- `rmax = 30`
- `tails = 6`
- `odd_even_kmax = 6`
- two `xgaussians`

### Old Versus New Build Time

| Step | Time |
| --- | ---: |
| old `build_basis(spec)` | about `173 s` |
| new `build_basis(spec)` | about `6.5 s` |

### Current Full Front-door Sequence

| Step | Time |
| --- | ---: |
| `build_basis(spec)` | about `6.5 s` |
| `basis_diagnostics(rb)` | about `0.6 s` |
| `radial_quadrature(rb)` | about `0.44 s` |
| `atomic_operators(rb, grid; Z = 2, lmax = 2)` | about `1.6 s` |
| full front-door sequence | about `9.4 s` |

So the standard radial path is no longer dominated by a construction cost large
enough to block normal public use.

## Validation

Validation run for this milestone:

- `julia --project=. test/runtests.jl`
- `julia --project=docs docs/make.jl`

## Current Trust Status

- the radial line is now coherent and trusted
- the public story is now consistent: users choose `rmax`, and the library
  owns the internal build and quadrature extents
- the interval-sampled construction path is trusted for the current radial
  line, not just as a provisional optimization

## Roadmap Consequence

- radial is no longer the main live stabilization front
- basis caching is no longer the first urgent fix
- the next likely repo priorities are:
  - promote the best current nonrecursive nesting state
  - or do the first serious radial-vs-Cartesian comparison

In short, this milestone changes the repo priorities because the radial branch
has moved from “still stabilizing” to “scientifically usable and operationally
coherent.”

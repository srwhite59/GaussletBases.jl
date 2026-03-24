# Radial Interval-sampled Build And Extents

## Pseudocode

1. Start from one `RadialBasisSpec`.
   The front-door radial recipe chooses:
   - family
   - one public extent: `rmax` or `count`
   - mapping
   - reference spacing
   - seed-window parameters
   - optional `xgaussian` supplement
   Code: `src/bases.jl`

2. Separate the public extent from the internal extents.
   The current contract is:
   - `rmax`: public last-center extent
   - `build_umax`: internal setup/build extent in mapped `u`
   - `quadrature_umax`: internal active quadrature extent in mapped `u`
   Code: `src/bases.jl`, `src/quadrature.jl`

3. Determine the internal build target in reference/mapped space.
   The construction logic derives:
   - odd-count seed window
   - internal build/setup extent
   - setup grid extent large enough for localization and cleanup
   Code: `src/bases.jl`

4. Build the setup-grid sampled layer with interval support, not dense full
   columns.
   For each shifted seed gausslet and each injected `xgaussian`, construct a
   setup-grid object only on the interval where it is numerically active.
   Code: `src/bases.jl`

5. Use the runtime machine-significant family support for that setup layer.
   The runtime family tables are trimmed after the last coefficient with
   `abs(c) >= 1e-17`, while the full high-precision family tables remain in the
   repo as internal reference data.
   Code: `src/families.jl`, `src/internal/families_high_prec.jl`

6. Assemble the seed overlap and position data by overlap-range weighted dot
   products.
   The interval-sampled layer supplies:
   - seed Gram matrices
   - seed position matrices
   - seed-to-`xgaussian` cross terms
   without storing long dense zero-tail columns.
   Code: `src/bases.jl`

7. Keep the later dense seed-space linear algebra unchanged.
   The present radial build still does:
   - odd/even cleanup
   - delta-direction removal
   - COMX localization
   - final QR/sign-fix flow
   - primitive contraction into the public `RadialBasis`
   Code: `src/bases.jl`

8. Build the final public basis in the same form as before.
   The interval-sampled objects are construction-time helpers only. Once the
   final localized basis is formed, they are discarded and the public result is
   the ordinary `RadialBasis`.
   Code: `src/bases.jl`

9. Derive the default quadrature extent from retained support, not setup
   padding.
   The current automatic rule takes:
   - retained radial centers
   - runtime family support width
   - retained `xgaussian` support
   and builds `quadrature_umax` from those active retained objects.
   Code: `src/bases.jl`, `src/quadrature.jl`

10. Build the physical quadrature grid from that retained-support extent.
    The automatic default path maps `quadrature_umax` to physical space and
    refines until the public-quality quadrature checks are satisfied.
    `quadrature_rmax` remains only as an expert override.
    Code: `src/quadrature.jl`, `src/diagnostics.jl`

## References

- Supporting notes:
  - `docs/radial_extent_policy_note.md`
  - `docs/radial_trust_milestone_note.md`
  - `docs/radial_default_path_cache_study.md`

## What This Builds

This page records the current trusted radial construction route:

- fast interval-sampled basis construction
- public `rmax` semantics separated from internal numerical extents
- automatic quadrature extent derived from retained support

It is the current source-of-truth algorithm page for the radial front door.

## Current Repo Status

The repo now has the radial line in the state described on this page:

- runtime family tables use machine-significant tails
- full high-precision family tables are preserved internally for reference
- basis construction uses interval-sampled setup-grid objects
- the standard radial path is no longer dominated by `build_basis`
- `rmax`, `build_umax`, and `quadrature_umax` are now distinct concepts in both
  code and docs

The normal public story is therefore:

- users choose `rmax`
- the library owns `build_umax`
- the library owns `quadrature_umax`

## Why This Matters

This algorithm page replaces two older confusions:

- treating `rmax` as if it directly set every internal extent
- treating the old dense setup-grid build as if it were the defining radial
  algorithm rather than an implementation bottleneck

The present route is both faster and conceptually cleaner.

## Implementation Notes

Recommended code-comment style:

```julia
# Alg Radial-Build step 4: Sample shifted seeds and xgaussians only on their
# active setup-grid intervals.
# See docs/src/algorithms/radial_interval_sampled_build_and_extents.md.
```

Guidelines:

- keep step numbers aligned with this page
- use this page for the trusted radial build/extents contract
- keep cache-policy or milestone numerics in supporting notes rather than
  growing this page into a status log

# Cartesian Nested Diatomic Coordinate Distortion

## Pseudocode

1. Restrict the first molecular distortion policy to one distinguished-axis
   geometry family.
   The first supported family is:
   - bond-aligned homonuclear diatomics
   - immediate extension to heteronuclear diatomics on the same axis
   - immediate extension to linear chains on that same axis
   Arbitrary non-linear geometries are deliberately deferred.
   Companion geometry page:
   [Cartesian nested diatomic box policy](cartesian_nested_diatomic_box_policy.md)

2. Keep the coordinate-slicing assumption.
   The molecular 3D basis still comes from three one-dimensional mappings:
   - one mapping on the distinguished molecular axis
   - one mapping on the first transverse axis
   - one mapping on the second transverse axis
   The rectangular tensor-product grid is retained; this is the point of
   contact with later box shrinkage and splitting.
   Historical support:
   - `QiuWhite_source.tex`
   - `PureGaussianGausslet.PGGbackbone3D(...)`

3. Define the target local core spacing as the primary physical constraint.
   The mapping policy should be written in terms of the desired local physical
   spacing at each nucleus, not only in terms of internal map parameters.
   For the first homonuclear diatomic pass:
   - use one common target local spacing `d_core`
   For later heteronuclear extension:
   - allow a per-atom target local spacing `d_core,I`
   - for the first concrete heteronuclear start, this means explicit targets
     such as `d_core,He` and `d_core,H`
   Historical support:
   - the paper states that the spacing at each nucleus should stay at the
     intended local core spacing
   - the legacy `PGGbackbone3D(...)` path uses `coresp` as that target spacing

4. On the bond axis, use a combined multi-center one-dimensional mapping.
   The bond-axis mapping density is the sum of atomic contributions, with one
   contribution per nucleus on that axis, plus an asymptotic constant tail
   density.
   The mapping must satisfy:
   - monotone increasing `u(z)`
   - asymptotic tail spacing set by the chosen far-field density / width
   - local spacing at each nucleus equal to the requested `d_core` or
     `d_core,I`
   - for the first homonuclear diatomic case, symmetry about the bond midpoint
     with the same target local spacing at both nuclei
   - for the first heteronuclear diatomic case, the same combined bond-axis
     family but with unequal per-atom local spacing targets carried honestly
     into the fit
   Historical support:
   - the paper explicitly says the molecular mapping density is the sum of the
     atomic ones, with a slight modification
   - `CoordinateMapping.getmapping(...)` builds exactly such a combined
     multi-center mapping from per-center contributions

5. Record the legacy implementation-facing rule on the bond axis.
   In the old PGDG line:
   - `getNGgaussletonly(...)` gathers the unique nuclear coordinates on each
     Cartesian axis
   - `PGGbackbone3D(...)` builds one backbone per axis
   - on a multi-center axis it calls `getmapping(...)`
   - `getmapping(...)` solves for a combined inverse-sqrt-density mapping whose
     local density at each nucleus matches the requested local spacing
   - `alignAtoms!(...)` then rescales the mapping so the nuclei sit close to
     integer `u`-grid locations
   Historical support:
   - `CoordinateMapping.jl`
   - `PureGaussianGausslet.jl`

6. On the transverse axes of a bond-aligned diatomic, use a single-center
   atomic-style mapping at the shared transverse coordinate.
   For a bond aligned with `z`, both nuclei share the same `x = 0` and `y = 0`
   projection. This shared projection is exactly why the first repo policy does
   **not** build a separate multi-center transverse mapping.
   Instead:
   - use one atomic-style transverse map centered at the shared transverse
     coordinate
   - choose its local target spacing from the same atomic core-spacing policy
   - for the first heteronuclear pass, if the requested spacings differ, use
     the finer of the two coincident projected spacings on that transverse
     axis
   - phrase that as the first practical rule for unequal atoms, not as a final
     general theorem
   - if the nuclei do not share the same transverse projection, that geometry
     is out of scope for this first page and should not be forced into this
     policy
   Historical support:
   - the paper says that for a linear chain only the chain direction needs the
     multi-center adjustment

7. Extend the same bond-axis rule to linear chains.
   For a chain aligned with the distinguished axis:
   - build one combined multi-center mapping on that axis from all nuclei
   - keep the two transverse directions on the same single-center
     atomic-style mapping family
   - only the distinguished axis needs the multi-center spacing adjustment
   Historical support:
   - `QiuWhite_source.tex`

8. Keep the first repo implementation close to the legacy combined
   inverse-sqrt-density style on the bond axis.
   The first modern implementation should not invent a qualitatively different
   multi-center family if the goal is to reproduce the old policy faithfully.
   The preferred first target is therefore:
   - a combined smooth one-dimensional molecular mapping on the bond axis
   - with the same conceptual ingredients as the old `getmapping(...)`
     inverse-sqrt construction
   This remains the intended first heteronuclear bond-axis family too.
   A fitted surrogate family is acceptable only if it preserves the same
   implementation-facing constraints:
   - local spacing at each nucleus
   - monotonicity
   - asymptotic tail spacing
   - smoothness good enough for the current mapped Cartesian basis machinery

9. Treat `alignAtoms!` as historical guidance, not a mandatory public API.
   The old line used `alignAtoms!(...)` as a further normalization/alignment
   step after building the combined bond-axis map. For the new repo, the key
   requirement is the resulting grid geometry:
   - nuclei should land close to stable parent-grid locations
   - the mapping should be reproducible and explicit
   The new code may realize that through a direct construction or a later
   alignment pass, but it should preserve the same practical effect.

10. Stop before arbitrary non-linear molecular distortion policies.
    This page does not settle:
    - bent molecules
    - axis-dependent splits without one distinguished molecular axis
    - general multi-center 3D distortion beyond coordinate-slicing
    Those need a different geometry language and should not be folded into this
    first diatomic policy.

## References

- Companion geometry page:
  [Cartesian nested diatomic box policy](cartesian_nested_diatomic_box_policy.md)
- Historical provenance:
  - `/Users/srw/Library/CloudStorage/Dropbox/chatarchive/references/papers/canonical/QiuWhite_source.tex`
  - `/Users/srw/Dropbox/GaussletModules/CoordinateMapping.jl`
  - `/Users/srw/Dropbox/GaussletModules/PureGaussianGausslet.jl`
  - `/Users/srw/Library/CloudStorage/Dropbox/chatarchive/reports/software_reviews/nestpgg3d_family_map_2026-03-15.md`
  - `/Users/srw/Library/CloudStorage/Dropbox/chatarchive/reports/software_reviews/white_lindsey_run_provenance_2026-03-15.md`

## What This Frames

This page records the first molecular coordinate-distortion policy for the
Cartesian nested line:

- one combined molecular mapping on the distinguished bond/chain axis
- two single-center atomic-style mappings on the transverse axes
- explicit target local spacing at each nucleus
- rectangular coordinate-slicing preserved so the box policy can operate on the
  resulting grid

It is the mapping half of the first diatomic geometry policy.

## Historical Support vs Modernized Policy

What comes directly from the paper and old code:

- molecular coordinate-slicing still uses one 1D mapping per Cartesian axis
- the molecular mapping density on a multi-center axis is the sum of atomic
  contributions, with a slight modification to keep the local spacing at each
  nucleus fixed
- for linear chains, only the chain direction needs that multi-center
  adjustment
- the old PGDG line builds axiswise backbones with `PGGbackbone3D(...)`,
  calls `getmapping(...)` on a multi-center axis, and then applies
  `alignAtoms!(...)`

What is a modernized repo policy choice:

- start with bond-aligned homonuclear diatomics
- extend next only to heteronuclear diatomics and linear chains on the same
  distinguished axis
- define the target in terms of requested local physical spacing `d_core` or
  `d_core,I`
- require the first homonuclear bond-axis map to be symmetric about the bond
  midpoint
- keep the first heteronuclear bond-axis map on the same combined
  inverse-sqrt-density family, but with unequal per-atom local spacing targets
- use a single-center atomic-style map on the transverse axes
- for unequal atoms, choose the transverse target spacing from the tighter /
  heavier side by taking the finer of the two requested local spacings
- rely on that transverse policy only when both nuclei project to the same
  transverse coordinate
- keep the first bond-axis implementation close to the legacy combined
  inverse-sqrt-density style rather than introducing a qualitatively different
  family immediately
- treat the legacy alignment step as an implementation device, not the public
  conceptual API

## Current Recommendation

This policy is precise enough for the first implementation pass if the code
stays within the intended scope:

- one distinguished bond axis
- one combined multi-center bond-axis map
- one atomic-style transverse map family
- one explicit local-spacing target per nucleus
- for the first homonuclear case, bond-axis symmetry about the midpoint
- shared transverse projection at the two nuclei

The first implementation should therefore target:

- bond-aligned ordinary-QW `HeH+` first
- heteronuclear nested fixed-block `HeH+` second
- one combined bond-axis inverse-sqrt-density map with explicit per-atom local
  spacing targets, such as `d_core,He` and `d_core,H`
- one shared transverse map on the common transverse projection, controlled by
  the finer of the two requested local spacings
- one validation that the resulting mapped centers support the split/no-split
  box policy on the companion diatomic box page
- the existing visualization/debug path kept in the loop before broader
  promotion: representative `xz` views and 3D center/source inspection through
  the same regenerate-first workflow already used on the `H2` line

It should not yet target arbitrary molecules or a general molecular mapping
API.

## Implementation Notes

Recommended code-comment style once implementation starts:

```julia
# Alg Nested-Diatomic-Map step 4: Build the bond-axis combined multi-center
# mapping so the local spacing at each nucleus matches the requested core
# spacing.
# See docs/src/algorithms/cartesian_nested_diatomic_coordinate_distortion.md.
```

Guidelines:

- keep this page focused on the one-dimensional distortion policy
- keep box shrink/split rules on the diatomic box-policy page
- keep atomic shell details on the atomic nonrecursive page
- do not extend this page to arbitrary non-linear geometries until that policy
  is genuinely settled

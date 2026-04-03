# Cartesian Nested Diatomic Box Policy

## Pseudocode

1. Restrict the first molecular nesting policy to one distinguished-axis
   geometry family.
   The first supported family is:
   - bond-aligned diatomics
   - linear chains aligned with the same distinguished axis
   Arbitrary non-linear geometries are deliberately deferred.
   References:
   - `QiuWhite_source.tex`
   - `GaussletModules/Boxes.jl`

2. Choose one Cartesian axis as the molecular axis and align the atoms to it.
   For the first pass, treat the bond or chain axis as the distinguished
   direction. The two transverse directions stay equivalent and are handled by
   the existing rectangular shell language.
   Historical support:
   - the paper uses coordinate-slicing on a rectangular grid
   - for linear chains, only the mapping parameters along the chain direction
     need special adjustment
   Reference: `QiuWhite_source.tex`

3. Start from one large shared rectangular parent box around all atoms.
   Before any split, the molecular fixed line is one shared box on the parent
   Cartesian grid. Shrink that box shell-by-shell at large radius using the
   same local shell language already established for the atomic case.
   Primitive shell language:
   [Cartesian nested face construction](cartesian_nested_face_construction.md)

4. Delay splitting until the parent box remains long enough along the
   distinguished axis to support two honest child boxes.
   The first split policy is:
   - split only along the distinguished bond/chain axis
   - require more than `2 * nside` raw sites along that long direction
   - require the split to leave each child with enough raw sites to continue as
     an atomic-style box rather than a thin sliver
   This is the first implementation-level meaning of:
   `N_parallel > 2 * nside`.
   Historical support:
   - `Boxes.jl` only splits when the long direction is sufficiently larger than
     the transverse directions and large enough relative to `doside`

5. Choose split planes by nearest physical midpoint in mapped-grid index
   space, and treat the homonuclear midpoint slab as a symmetry special case.
   For a diatomic there is one midpoint. For a chain, order the atoms along the
   distinguished axis and use the physical midpoint between each neighboring
   pair. The split plane should stay attached to the actual mapped grid rather
   than to a continuous idealized coordinate.
   The first heteronuclear rule is deliberately narrow:
   - split only along the distinguished bond axis
   - choose the split index as the parent-grid index nearest the physical
     midpoint between the two nuclei
   - do not reserve a shared midpoint slab by default
   - if the discrete midpoint choice is tied, assign the extra row to the
     tighter/heavier side rather than inventing a heteronuclear slab rule
   The first explicit homonuclear correction is narrower and symmetry-driven:
   - for a bond-aligned homonuclear diatomic whose bond-axis working interval
     has odd length, reserve the midpoint index as a shared `nx × ny × 1` slab
   - split the remaining bond-axis rows into equal left/right child boxes
   - treat the midpoint slab as a direct shared region between the two child
     subtrees rather than assigning it to either child
   So the midpoint slab is a homonuclear symmetry special case, not the
   default diatomic rule.
   Historical support:
   - `Boxes.jl` computes physical midpoints between neighboring atoms and picks
     the nearest grid index

6. Require child boxes to be roughly cubic in physical extent.
   A permitted split must not create long thin slivers. After mapping the index
   box to physical coordinates:
   - the bond-axis physical width of each child should stay comparable to the
     transverse physical widths
   - if the split would make the child much narrower than the transverse box,
     do not split yet and continue shell shrinkage on the shared parent box
   Historical support:
   - `Boxes.jl` rejects narrow splits and adds extra bars when needed to keep
     shells more square

7. Keep shell resolution roughly matched across directions.
   The current fixed retain tuples are only provisional. The guiding policy is:
   - local side contractions should aim to preserve roughly matched physical
     spacing across shell directions
   - if a fixed global retain tuple visibly violates that constant-resolution
     idea, it should be treated as provisional rather than as the intended
     long-term rule
   The current landed homonuclear correction is intentionally narrower than a
   full adaptive rule:
   - at the 1D `doside` / `COMX` boundary, if the local interval is symmetric
     about zero and the provisional retained count is even, reduce it by one
     so the localized side space keeps an odd center-bearing pattern
   - this is a structural symmetry correction at the actual contraction
     boundary, not yet a general adaptive retain-count formula
   - this structural `doside` correction does not by itself settle the broader
     heteronuclear or chain box policy
   This page does not yet freeze an adaptive local-side-count formula. That
   should wait for explicit 1D `doside` / `COMX` diagnostics.
   Endcap slabs are also still deferred by default for the first
   heteronuclear start.

8. After the split, treat each child as an atomic-style subtree.
   Once the shared parent box has split:
   - each child box inherits the same shell language as the atomic route
   - each child shrinks shell-by-shell independently
   - each child remains defined on original parent-space rows assigned to that
     box region, not by re-coarsening already-renormalized functions
   If a homonuclear midpoint slab is present, it stays as a direct shared
   region between the child subtrees rather than belonging to either one.
   Landed atomic route:
   [Cartesian nested atomic nonrecursive route](cartesian_nested_atomic_nonrecursive_route.md)

9. Extend to linear chains by repeated midpoint splitting on the same axis.
   The immediate extension beyond a diatomic is:
   - one shared parent box for the full chain
   - repeated midpoint splits only along the distinguished chain axis
   - atomic-style shell subtrees on the resulting children
   This keeps the policy inside one geometry family:
   one distinguished axis plus rectangular shells.

10. Stop before arbitrary non-linear molecular box policies.
   The first diatomic/chain page does not settle:
   - bent geometries
   - multiple distinguished axes
   - general Voronoi-like or adaptive 3D cell policies
   Those require a different geometry language and should not be smuggled into
   this first implementation.

## References

- Primitive shell language:
  [Cartesian nested face construction](cartesian_nested_face_construction.md)
- Landed atomic subtree route:
  [Cartesian nested atomic nonrecursive route](cartesian_nested_atomic_nonrecursive_route.md)
- Historical provenance:
  - `/Users/srw/Library/CloudStorage/Dropbox/chatarchive/reports/software_reviews/white_lindsey_run_provenance_2026-03-15.md`
  - `/Users/srw/Library/CloudStorage/Dropbox/chatarchive/reports/software_reviews/nestpgg3d_family_map_2026-03-15.md`
  - `/Users/srw/Dropbox/GaussletModules/Boxes.jl`
  - `/Users/srw/Library/CloudStorage/Dropbox/chatarchive/references/papers/canonical/QiuWhite_source.tex`

## What This Frames

This page records the intended first box policy for molecular Cartesian
nesting:

- one large shared rectangular box around a bond-aligned diatomic or linear
  chain
- shell shrinkage at large radius before any split
- midpoint-based splits only along the distinguished molecular axis, with the
  shared midpoint slab reserved only for odd-length homonuclear working
  intervals
- child boxes that then continue as atomic-style shell subtrees

It is a geometry-policy page, not an implementation page.

## Historical Support vs Modernized Policy

What is historically supported by the paper/provenance:

- molecular coordinate-slicing still lives on a distorted rectangular grid
- linear chains are a special one-axis family
- the legacy `Boxes.jl` logic shrinks boxes shell-by-shell and uses midpoint
  splits on the long direction with anti-sliver safeguards
- the paper-era driver family clearly targeted diatomics and linear chains

What is a modernized repo policy choice:

- start with bond-aligned diatomics only
- extend next to heteronuclear diatomics and linear chains on the same
  distinguished axis
- make the split rule explicit as nearest-physical-midpoint index selection
  plus a homonuclear-only midpoint-slab correction
- require `N_parallel > 2 * nside` before splitting
- require child boxes to stay roughly cubic in physical extent
- record constant-resolution shell matching as a policy guideline, with one
  landed narrow odd-on-symmetric `doside` correction but without yet fixing
  the general adaptive formula
- keep endcap slabs deferred by default for the first heteronuclear pass
- defer arbitrary non-linear geometries until a different policy is written

## Current Recommendation

This policy is now precise enough for the first heteronuclear implementation
pass if the code stays within the intended scope:

- one distinguished molecular axis
- one shared parent box first
- midpoint splits only along that axis
- explicit shared midpoint slab for odd-length homonuclear working intervals
- atomic-style shell language reused after the split

The first narrow heteronuclear implementation now follows this target shape:

- bond-aligned ordinary-QW `HeH+` first
- heteronuclear nested fixed-block `HeH+` second
- the same split-eligibility and anti-sliver checks already used on the
  homonuclear line
- nearest-physical-midpoint split placement without a default heteronuclear
  midpoint slab
- endcap slabs still deferred by default
- the existing visualization/debug path used to check the first heteronuclear
  geometry before broader promotion:
  ordinary and nested representative `xz` views, plus 3D center/source
  inspection through the same regenerate-first workflow already used on `H2`

It should not yet target arbitrary molecules or a final adaptive local-side
formula.

## Implementation Notes

Recommended code-comment style while this narrow heteronuclear line remains in
checkpoint scope:

```julia
# Alg Nested-Diatomic step 5: For an odd-length homonuclear working interval,
# reserve the midpoint row as a shared slab and split the remaining bond-axis
# rows into equal left/right child boxes.
# See docs/src/algorithms/cartesian_nested_diatomic_box_policy.md.
```

Guidelines:

- keep this page focused on box-policy decisions
- keep atomic shell details on the atomic nonrecursive page
- keep residual-Gaussian completion on the QW residual-Gaussian page
- do not extend this page to arbitrary non-linear geometries until that policy
  is genuinely settled

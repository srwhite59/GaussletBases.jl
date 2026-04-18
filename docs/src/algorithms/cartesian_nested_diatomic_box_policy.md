# Cartesian Nested Diatomic Box Policy

## Pseudocode

1. Restrict this policy to the landed bond-aligned homonuclear diatomic nested
   source builder.
   The geometry family is:
   - one distinguished bond axis
   - two equivalent transverse axes
   - one shared mapped Cartesian parent box before any split

   This page does not define the final policy for arbitrary molecules, and it
   does not replace the separate chain or lattice policies.

2. Start from one shared working box and peel outer shells first.
   The source builder still uses the established rectangular-shell language on
   the parent mapped grid:
   - shrink the working box shell-by-shell at large radius
   - build shell layers directly from the real mapped basis lines
   - defer any bond-axis split until the current unsplit box is genuinely
     elongated

   Reference:
   [Cartesian nested face construction](cartesian_nested_face_construction.md)

3. Decide whether splitting is even allowed from the unsplit box shape.
   The landed split policy uses a dedicated unsplit-box aspect guard:

   - `min_unsplit_parallel_to_transverse_ratio_for_split = 3.0`

   Let `W_parallel` be the physical width of the current working box along the
   bond axis, and `W_short` the shorter of the two transverse physical widths.
   Splitting is not considered unless:

   - `W_parallel > 3.0 * W_short`

   This is a separate gate from the child anti-sliver test below. It answers:
   “is the unsplit box elongated enough that a split is worth considering at
   all?”

4. Keep the existing child anti-sliver guard as a second, separate test.
   After a candidate split plane is chosen, the resulting child boxes must
   still be roughly cubic in physical extent. The landed child guard remains:

   - `min_parallel_to_transverse_ratio = 0.4`

   For each child, the bond-axis physical width must satisfy:

   - `W_child_parallel >= 0.4 * max(W_child_x, W_child_y)`

   This guard answers a different question from step 3:
   “if a split is attempted, do the resulting children remain honest 3D boxes
   instead of slivers?”

5. Choose the split plane on the real mapped bond axis.
   When the split guards are satisfied:

   - split only along the distinguished bond axis
   - choose the split index from the actual mapped bond-axis line
   - use the nearest physical midpoint index

   For the bond-aligned homonuclear route, an odd bond-axis working interval
   may reserve the midpoint row as a shared slab:

   - `nx × ny × 1`
   - shared between the left and right child subtrees

   This midpoint slab is a homonuclear symmetry special case. It is not the
   general rule that decides whether splitting happens.

6. Build an ideal cube-style angular reference from the short-side scale.
   The landed retain-count policy no longer treats the distorted current box as
   its own standard. Instead it constructs an ideal local reference:

   - take the shortest physical side of the current box
   - build an ideal cube from that scale
   - evaluate a target angular band on that ideal cube
   - widen the acceptable band by:
     - `reference_fudge_factor = 1.2`

   This reference provides the comparison target for bond-axis retain-count
   decisions on both shells and the inner core.

7. Choose shell bond-axis retain counts adaptively against that angular band.
   The outer shell stage is still the ordinary complete-shell stage, but the
   bond-axis retain count is no longer hard-wired to one fixed long-side
   contraction.

   For each shell segment on the long direction:

   - generate candidate retained counts
   - compare each candidate against the ideal-reference angular band
   - choose the smallest acceptable retained count
   - if the direct parent line is already too coarse, treat the case as
     parent-limited rather than blaming the contraction

   So the landed shell policy is:
   - ordinary shell language on the outside
   - adaptive bond-axis retain count on the long direction
   - no fixed “always keep the long side at count X” rule

8. Use a nonuniform inner-core policy when the remaining core is still
   elongated.
   The important landed change is inside the unsplit core:

   - do not force one uniform retained count for the whole inner block
   - treat the bond-axis contraction nonuniformly across transverse rows
   - evaluate each row against the same ideal-reference angular band

   In other words, the core is trimmed row-by-row in the transverse plane
   rather than assigning one global bond-axis retain count to the whole inner
   region.

9. Protect near-nucleus rows from inner-core trimming.
   The landed auto protection rule is:

   - `protect_rows = div(nside, 2) - 1`

   Examples:
   - `ns = 5 -> 1`
   - `ns = 7 -> 2`
   - `ns = 9 -> 3`

   This protection applies only to the nonuniform inner-core trimming policy.
   Rows near the nuclear axis are forced direct there instead of being trimmed.

10. Build the source from the real basis, not a sandbox surrogate.
    The real implementation uses:
    - the actual mapped basis lines carried by the repo basis object
    - the actual QW / L&W `d = core_spacing` used to build that basis

    The algorithm should therefore be read as a policy on the real mapped
    basis, not as a rule on an idealized uniform grid.

## Code Pointers

- Split-eligibility geometry:
  - `src/cartesian_nested_diatomic.jl:_nested_bond_aligned_diatomic_split_geometry`
- Ideal-reference angular band:
  - `src/cartesian_nested_diatomic.jl:_nested_diatomic_reference_band`
- Adaptive shell retain-count choice:
  - `src/cartesian_nested_diatomic.jl:_nested_diatomic_adaptive_shell_retention`
- Nonuniform inner-core construction:
  - `src/cartesian_nested_diatomic.jl:_nested_bond_aligned_diatomic_nonuniform_core_block`
  - `src/cartesian_nested_diatomic.jl:_nested_bond_aligned_diatomic_sequence_for_box`
- Auto near-nucleus protection:
  - `src/cartesian_nested_diatomic.jl:_nested_diatomic_resolve_core_near_nucleus_protect_rows`
- Real diatomic source entry point:
  - `src/ordinary_qw_nested_frontends.jl:bond_aligned_diatomic_nested_fixed_source`
- Real diagnostics entry point:
  - `src/ordinary_qw_nested_frontends.jl:bond_aligned_diatomic_nested_geometry_diagnostics`

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

This page now describes the landed real-repo box policy for the bond-aligned
homonuclear diatomic nested source builder:

- one shared parent box first
- shell shrinkage before splitting
- a two-guard split rule
- adaptive shell retain counts against an ideal cube-style angular reference
- nonuniform inner-core trimming on the bond axis
- near-nucleus row protection during that inner-core trimming

It is still a geometry-policy page. It does not redefine supplement selection,
residual orthogonalization, overlap assembly, or bundle storage.

## Real Implementation Notes

The real repo path is close to the 2D sandbox policy that motivated this
rewrite, but it is not identical. The important recorded correction is:

- in the real path, `ns = 5`, `R = 1.0` remains unsplit

So the public algorithm page should not simply copy the sandbox split table.
The real source builder must remain the referee.

The representative real He2 scan used the spacing rule:

- `d = 1.2 / (Z * (ns - 3))`

For He:
- `ns = 5 -> d = 0.3`
- `ns = 7 -> d = 0.15`
- `ns = 9 -> d = 0.10`

Qualitatively, the landed real path gives:
- `ns = 5`: split at `R = 0.8, 1.4, 2.0, 4.0, 6.0`, but unsplit at `R = 1.0`
- `ns = 7`: unsplit only at `R = 0.8`
- `ns = 9`: unsplit only at `R = 0.8`

That is close to the sandbox trend, but the real code path is the one that
should be documented.

## Current Recommendation and Support Boundary

The landed support boundary is:

- the updated source geometry policy is active on the real bond-aligned
  homonuclear diatomic nested source path
- the same helper machinery now also affects some diatomic split diagnostics
  shared with heteronuclear fixtures
- this page does not claim that chain, lattice, or arbitrary molecular routes
  already use the same full policy

So the current public recommendation is:

- use this page as the algorithm reference for the bond-aligned homonuclear
  diatomic source builder
- treat heteronuclear, chain, and broader molecular routes as separate support
  boundaries unless their pages are updated explicitly

This keeps the documentation aligned with the landed code rather than
overclaiming unlanded generalizations.

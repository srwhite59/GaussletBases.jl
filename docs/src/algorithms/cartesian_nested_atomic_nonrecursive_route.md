# Cartesian Nested Atomic Nonrecursive Route

## Pseudocode

1. Start from one finalized Cartesian parent fixed block on the atomic QW line.
   In the current repo this is the finalized Cartesian fixed block produced
   before residual-Gaussian completion.
   Code: `src/ordinary_qw_operator_assembly.jl`,
   `src/ordinary_mapped_backends.jl`

2. Choose one centered rectangular working cube in the parent index lattice.
   The nonrecursive atomic route works on that parent cube only. It does not
   re-coarsen already-renormalized shell functions.
   Code: `src/cartesian_nested_faces.jl`

3. Decompose the outer shell annuli of that cube into disjoint geometric
   strata.
   For each shell layer:
   - faces: exactly one boundary coordinate
   - edges: exactly two boundary coordinates
   - corners: exactly three boundary coordinates
   - retained core: no boundary coordinates
   Code: `src/cartesian_nested_faces.jl`

4. Build the local contraction objects on the original parent-space rows
   assigned to each stratum.
   The primitive rules are:
   - face piece: tangential `doside × doside`
   - edge piece: one free-interval `doside`
   - corner piece: direct corner column in the first pass
   - retained core: direct parent block at the present nonrecursive level
   Code: `src/cartesian_nested_faces.jl`
   Primitive algorithm page:
   [Cartesian nested face construction](cartesian_nested_face_construction.md)

5. Enforce both disjointness and completeness on the working cube.
   Every parent fixed-space row in the working cube must belong to exactly one
   object in the source language:
   - shell face piece
   - shell edge piece
   - shell corner piece
   - retained core block
   Code: `src/cartesian_nested_faces.jl`

6. Assemble the nonrecursive nested fixed source from those shell layers plus
   the retained core.
   The current landed atomic line has two important nonrecursive states:
   - shell-plus-core
   - corrected complete-shell
   The corrected complete-shell line is the first trusted reduced atomic
   anchor beyond shell-plus-core.
   Code: `src/cartesian_nested_faces.jl`

7. Propagate the fixed-fixed operator packet through the same contraction maps.
   The carried packet includes:
   - overlap
   - kinetic
   - position
   - `x^2`
   - Gaussian-factor terms
   - pair-factor terms
   The IDA pair transfer must be weight-aware:

   ```math
   M = D(w)\,V\,D(w), \qquad
   M' = C^T M C, \qquad
   w' = C^T w, \qquad
   V' = D(1/w')\,M'\,D(1/w').
   ```

   Code: `src/cartesian_nested_faces.jl`

8. Build the nested fixed-block adapter without changing the downstream
   consumer algebra.
   The adapter carries:
   - the combined coefficient matrix `C`
   - transformed fixed-block weights
   - the transformed fixed-fixed packet
   - fixed-Gaussian cross blocks contracted from parent raw blocks through `C`
   Gaussian-Gaussian blocks stay on the existing analytic route.
   Code: `src/cartesian_nested_faces.jl`, `src/ordinary_qw_raw_blocks.jl`,
   `src/ordinary_qw_operator_assembly.jl`

9. Hand the adapted fixed block to the existing QW residual-Gaussian route.
   This page stops at the nested fixed-block / transferred-packet stage. The
   residual-Gaussian complement, orthogonalization, and final hybrid completion
   remain the responsibility of:
   [Qiu-White residual-Gaussian route](qiu_white_residual_gaussian_route.md)

10. Judge the nonrecursive atomic route by fixed-block fidelity and reduced
    hybrid quality, not by shell geometry alone.
    The trusted diagnostics are:
    - overlap quality
    - finite positive transformed weights
    - projected fixed-only interaction transfer
    - parent low-energy one-body capture
    - final nearest/GGT `E1`
    - final nearest/GGT `⟨Vee⟩`
    - dimension/runtime against the unnested and shell-plus-core anchors

## References

- Supporting notes:
  - `docs/cartesian_nested_shell_plus_core.md`
  - `docs/cartesian_nested_complete_shell_layer.md`
  - `docs/cartesian_nested_sequence_coverage_fix.md`
  - `docs/cartesian_nested_ida_weight_transfer_fix.md`
  - `docs/atomic_hybrid_anchor_comparison.md`
- Primitive algorithm page:
  [Cartesian nested face construction](cartesian_nested_face_construction.md)
- Fixed-block consumer handoff:
  [Qiu-White residual-Gaussian route](qiu_white_residual_gaussian_route.md)

## What This Builds

This page records the landed nonrecursive atomic nesting route on the current
Cartesian QW line:

- one parent working cube in the original parent-space basis
- one complete shell-language decomposition of that cube
- one reduced fixed block formed from shell pieces plus retained core
- one transferred fixed-block operator packet in the same consumer language the
  QW residual-Gaussian route already reads

It does not describe:

- strong algebraic recursion
- diatomic box policy
- residual-Gaussian completion

Those are deliberately separate questions.

## Current Repo Status

The repo now has the full nonrecursive atomic fixed-block route described on
this page:

- primitive `doside` and face-product pieces
- complete shell layers including faces, edges, and corners
- coverage enforcement on the working cube
- shell-plus-core and corrected complete-shell fixed sources
- nested fixed-block adapters consumed by the current nearest/GGT QW path
- weight-aware fixed-block IDA transfer

The current trusted atomic nonrecursive states are:

- unnested hybrid QW reference
- shell-plus-core hybrid as the conservative nested anchor
- corrected complete-shell hybrid as the trusted reduced atomic anchor

The next atomic step, if atomic hierarchy work is ever reopened, should be
judged against that corrected complete-shell anchor rather than against earlier
primitive shell milestones.

## Relation To Other Pages

- [Cartesian nested face construction](cartesian_nested_face_construction.md)
  records the primitive local contraction language:
  `doside`, face products, and the first shell packet.
- This page records the landed nonrecursive atomic fixed-block route built from
  those primitives.
- [Qiu-White residual-Gaussian route](qiu_white_residual_gaussian_route.md)
  then takes over once the nested fixed block and its transferred packet are in
  hand.

## Code Map

The current landed implementation is concentrated in
`src/cartesian_nested_faces.jl`:

- `_nested_doside_1d(...)`
  primitive local `doside` construction for step `4`
- `_nested_complete_rectangular_shell(...)`
  complete shell-layer construction with faces, edges, and corners for steps
  `3`-`6`
- `_nested_assert_sequence_coverage(...)`
  explicit coverage/disjointness enforcement for step `5`
- `_nested_shell_plus_core(...)`
  shell-plus-core source construction for step `6`
- `_nested_shell_sequence_from_core_block(...)`
  general nonrecursive shell-sequence assembly for step `6`
- `_nested_fixed_block(...)`
  fixed-block adapter construction for step `8`

The downstream handoff point lives in `src/ordinary_qw_operator_assembly.jl`:

- `ordinary_cartesian_qiu_white_operators(fixed_block::_NestedFixedBlock3D, ...)`
  consumes the adapted fixed block and then follows the QW residual-Gaussian
  route from step `9` onward

The code comments on this line are present but still a bit uneven; this page is
the place to tighten that mapping over time.

## Implementation Notes

Recommended code-comment style:

```julia
# Alg Nested-Atomic step 7: Transfer the fixed-block IDA packet through the
# nonrecursive shell map with weight-aware reweighting.
# See docs/src/algorithms/cartesian_nested_atomic_nonrecursive_route.md.
```

Guidelines:

- keep step numbers aligned with this page
- use this page for the landed nonrecursive atomic route
- use the face-construction page only for primitive local contraction logic
- keep residual-Gaussian completion comments on the QW residual-Gaussian page

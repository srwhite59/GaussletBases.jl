# Round 005 Reconciliation

Status: reconciled into Slice A freeze candidate

Repo-manager reconciled the recursive-projection freeze review into
`docs/src/developer/cartesian_hamiltonian_producer_design.md`. The design
remains unapproved and is not implementation authority.

## Accepted Changes

- Clarified that `HP-FN-00` performs projection plus shell-local Lowdin only and
  does not own final sign canonicalization.
- Assigned support-weight derivation, direct-weight validation, sign
  canonicalization, and `CartesianTerminalBasisBlock` construction to
  `HP-FN-01`.
- Added the projection threshold rule: if a previous-block residual is already
  below `projection_atol`, do not subtract and do not enlarge effective support.
- Replaced "factorized or incremental" wording with the explicit block-action
  contract `C_left' * S_lr * C_right`.
- Stated that at most one terminal support-pair workspace may be live, with
  tiling required when it exceeds the cap.
- Removed per-object/helper line targets; Slice A now has one budget target and
  redesign threshold.
- Required the reviewed Cr2 fixture to produce a real terminal basis. Distorted
  COMX rejection is allowed only when the typed inventory actually contains
  distorted COMX.

## Deferred Or Rejected Changes

- Did not bind Slice A IDs.
- Did not bind `AGENTS.md`.
- Did not touch `src`, `test`, `bin`, or `tools`.
- Did not run the final recursive-projection spike.

## Freeze Implication

The design is ready for one final uncommitted recursive-projection spike. If the
spike passes, Slice A can be frozen without another broad review round.

## Implementation Freeze Status

Not frozen. Not implementation authority. No source implementation should begin
until Slice A approval is explicit.

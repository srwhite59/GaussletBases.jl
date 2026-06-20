# Round 002 Reconciliation

Status: reconciled into Slice A freeze candidate

Repo-manager reconciled the ChatGPT-Pro milestone review into
`docs/src/developer/cartesian_hamiltonian_producer_design.md`. The design
remains unapproved and is not implementation authority.

## Accepted Changes

- Status changed to draft v3 Slice A freeze candidate.
- Scope now says only Slice A should freeze first; B/C/D remain future
  candidates until Slice A reports empirical support and memory facts.
- Added effective-support semantics for terminal basis blocks.
- Moved sign canonicalization into Slice A terminal basis finalization.
- Added direct-sector validations for identity overlap and positive finite IDA
  weights.
- Rejected `HP-RES-01` as a persistent terminal-basis result wrapper.
- Added `HP-FN-00`, an explicit previous-block projected shell realization
  operation.
- Rejected/deferred `HP-CHANGE-01` as insufficient for the corrected
  projection problem.
- Removed `parent_dims` from `CartesianTerminalBasisRealization`.
- Marked `HP-FN-03`, `HP-FN-04`, and `HP-FN-05` as future candidates.
- Added Slice A target/redesign-threshold/net-deletion budget.
- Added an allowed uncommitted numerical spike with required measurements.
- Updated the mechanical gate expected Slice A design ID set.

## Deferred Or Rejected Changes

- Did not bind `AGENTS.md`.
- Did not mark any Slice A ID as approved implementation authority.
- Did not touch `src`, `test`, `bin`, or `tools`.
- Did not settle B/C operator ownership beyond leaving those surfaces future
  candidates.
- Did not split `HP-FN-04`; this is deferred until after Slice A data exists.

## Design V3 Summary

Design v3 is a Slice A freeze candidate. The only surfaces intended for near
approval are:

- `HP-OBJ-01`
- `HP-OBJ-02`
- `HP-FILE-01`
- `HP-FN-00`
- `HP-FN-01`
- `HP-FN-02`

The intended first implementation slice is terminal basis realization only. It
must project each PQS shell against all previous blocks, return effective
support plus sign-canonicalized coefficients, audit cross-block overlaps, and
delete the current terminal preflight path in the same merge.

## Remaining Questions Before Freezing Slice A

1. Should the listed Slice A freeze candidates now be marked approved, or should
   they receive one final manager/user signoff first?
2. Is `HP-FN-00` specific enough for implementation, or should the numerical
   spike run before approval?
3. Is the 225 added-source-line redesign threshold acceptable for Slice A?

## Implementation Freeze Status

Not frozen. Not implementation authority. No source implementation should begin
until Slice A approval is explicit.

# Round 003 Reconciliation

Status: reconciled into Slice A freeze candidate

Repo-manager reconciled the atomic/diatomic unification recommendation into
`docs/src/developer/cartesian_hamiltonian_producer_design.md`. The design
remains unapproved and is not implementation authority.

## Accepted Changes

- Added a binding atomic/diatomic unification invariant.
- Stated that one-center and bond-aligned diatomic geometry may differ only
  before terminal support, retained, and transform records exist.
- Required one-center atomic, contact-core H2, and separated Cr2 terminal plans
  to use the same terminal-basis realization entry point in Slice A.
- Added `HP-WIRE-01`, generic terminal-basis stage integration.
- Forbade Slice A dispatch on system classification, atom count, route kind,
  bond axis, and terminal role vocabulary.
- Updated Slice A validation and the allowed numerical spike to include
  one-center terminal records.
- Recorded that the old atomic multilayer/common-H1 path is migration-only and
  should be deleted after the generic one-body route reproduces the atomic
  endpoint.
- Added future-slice constraints that one-body, IDA, and final producer code
  consume `CartesianTerminalBasisRealization` without atomic/diatomic branches.

## Deferred Or Rejected Changes

- Did not require the full atomic Hamiltonian in Slice A.
- Did not add a committed atomic smoke test.
- Did not bind `AGENTS.md`.
- Did not touch `src`, `test`, `bin`, or `tools`.

## Design V3 Update Summary

Slice A now freezes a generic terminal-basis boundary, not merely a generic
diatomic boundary. The intended first implementation slice must realize the
terminal basis for one-center atomic, H2, and Cr2 records through the same
`pqs_terminal_basis_realization(...)` path and return the same
`CartesianTerminalBasisRealization` object.

## Implementation Freeze Status

Not frozen. Not implementation authority. No source implementation should begin
until Slice A approval is explicit.

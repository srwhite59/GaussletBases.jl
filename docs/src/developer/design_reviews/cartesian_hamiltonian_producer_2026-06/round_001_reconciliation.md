# Round 001 Reconciliation

Status: reconciled into design v2 candidate

Repo-manager reconciled `round_001_consolidated_review.md` into
`docs/src/developer/cartesian_hamiltonian_producer_design.md`. The design
remains unapproved and is not implementation authority.

## Inputs

- `round_001_prompt.md`
- `round_001_consolidated_review.md`
- `docs/src/developer/cartesian_hamiltonian_producer_design.md`

## Accepted Changes

- Changed the design status to a draft v2 candidate for milestone review.
- Demoted positive HP-* registry entries from approved to candidate. The draft
  now says only a later frozen design revision can mark items approved for
  implementation.
- Added the explicit rule that each PQS shell candidate must be projected out
  of every already-retained earlier terminal block in terminal order before
  shell-local Lowdin.
- Stated that `:terminal_pqs_cross_block_projection_required` is a construction
  failure, not a mergeable endpoint.
- Added a localized one-basis IDA convention:
  `electron_electron_ida[a,b]` is consumed by
  `(pq|rs)_IDA = sum_ab X[a,p]X[a,q] V[a,b] X[b,r]X[b,s]`.
- Added the sign-canonicalized positive localized IDA weight requirement.
- Added the `support_states` derivation invariant for
  `CartesianTerminalBasisBlock`.
- Added the `basis === nothing` iff `blocker !== nothing` invariant for the
  two-field terminal-basis result.
- Updated `HP-FN-03` to be an accumulating helper with `scale::Float64 = 1.0`.
- Marked owners for `HP-FN-03` and `HP-FN-04` as unresolved before approval.
- Added an IDA block-pair tiling rule when a full support-pair workspace would
  exceed the reviewed Cr2 memory budget.
- Changed implementation line budgets from hard caps to manager-reviewed
  targets when correctness or validation would otherwise be at risk.

## Rejected Or Deferred Changes

- Did not approve any HP-* item as final implementation authority. The whole
  registry remains candidate/rejected until a later freeze.
- Did not split `HP-FN-04` into multiple helper IDs yet. The design now records
  it as candidate and underspecified; ChatGPT-Pro should advise whether to
  split it or keep one helper with subrequirements.
- Did not choose the exact owner file for one-body and IDA assembly. The design
  names candidate ownership and asks for review.
- Did not add source, test, tool, driver, or `AGENTS.md` changes.
- Did not settle whether the current projected-shell descriptor can implement
  accumulated previous-block projection directly or whether one small explicit
  projection operation must be registered.

## Design V2 Summary

Design v2 is still a compact authority candidate, not an approved authority. It
now has a clearer numerical contract:

- terminal PQS shells are projected against all earlier retained terminal
  blocks before shell-local Lowdin;
- cross-block overlap failure stops construction;
- direct identity sectors remain implicit;
- one-body operators and IDA are assembled in the same canonicalized localized
  final basis;
- `CartesianIDAHamiltonian` receives `K`, separated unit-center `U_A`, and a
  one-basis localized `electron_electron_ida`;
- IDA construction must tile large support-pair workspaces when required by the
  Cr2 memory budget.

The registry is now a candidate registry. Implementation may not begin from it
until manager/user approval freezes specific IDs.

## Remaining Questions For ChatGPT-Pro

1. Is the Section 4.3 one-basis localized IDA convention the right public
   mathematical contract for `CartesianIDAHamiltonian`?
2. Should `HP-FN-04` remain one candidate helper with subrequirements, or should
   sign-gauge, raw pair-factor tiling, and final contraction be separate design
   IDs?
3. Can previous-block projection be safely expressed through the existing
   projected-shell descriptor path, or should the design approve one explicit
   projection operation?
4. Which existing file should own blockwise one-body and IDA assembly without
   reviving the retired pair-materialization framework?
5. Is the 1.2 GiB Cr2 peak-memory target appropriate, or should v2 define a
   lower target plus an absolute hard cap?
6. Are any candidate objects still unnecessary before Slice A?

## Implementation Freeze Status

Not frozen. Not implementation authority. Do not bind `AGENTS.md` or start
source implementation until the milestone review is complete and a later
revision explicitly approves the production surfaces.

# Round 004 Reconciliation

Status: reconciled into Slice A freeze candidate

Repo-manager reconciled the uncommitted terminal projection spike into
`docs/src/developer/cartesian_hamiltonian_producer_design.md`. The design
remains unapproved and is not implementation authority.

## Accepted Changes

- Clarified that shell projection, shell overlap, and Lowdin cleanup are not
  terminal contract fields today and must be constructed inside the terminal
  realizer, not added as metadata or staged summaries.
- Required production Slice A to use recursively projected coefficients and
  effective supports from previous PQS blocks.
- Required later direct records to be cross-checked against all earlier direct
  and PQS blocks.
- Clarified that one-center terminal records must be connected to typed
  terminal records if the current one-center path still needs route-skeleton
  shape input.
- Added a Cr2 performance guardrail: projection and overlap audits should use
  factorized or incremental support-overlap construction rather than dense
  scratch matrices as the production strategy.
- Updated Slice A validation and spike reporting to distinguish true recursive
  projection from the spike's shell-local approximation.

## Deferred Or Rejected Changes

- Did not freeze Slice A IDs yet.
- Did not bind `AGENTS.md`.
- Did not touch `src`, `test`, `bin`, or `tools`.
- Did not require full one-body or IDA assembly in Slice A.

## Freeze Implication

The spike supports the Slice A direction: terminal records are reachable for
one-center, H2, and Cr2; ranks remained stable; and projection appeared
well-conditioned. It also prevents immediate freeze because the production
design must explicitly require recursive projected coefficients and avoid the
dense scratch pattern used by the spike.

## Implementation Freeze Status

Not frozen. Not implementation authority. No source implementation should begin
until Slice A approval is explicit.

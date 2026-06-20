# Cartesian Hamiltonian Producer Design Reviews

This directory tracks review rounds for
`docs/src/developer/cartesian_hamiltonian_producer_design.md`.

The design document is still a draft. Review files in this directory are
discussion artifacts, not implementation authority. Production source work may
start only after a reviewed design revision is explicitly approved and bound by
repo policy.

## Operating Model

Repo-manager owns this lane.

Each round should produce:

- a reviewer prompt;
- one consolidated review, or separate numerical/performance/minimality reviews
  when explicitly requested;
- a reconciliation note;
- a revised design document, if changes are accepted.

Reviewers may propose document edits only. They must not edit `src`, `test`,
`bin`, or `tools`.

ChatGPT-Pro should be used as a milestone reviewer after local review and
manager reconciliation, not as a per-pass checker.

## Current Milestone

Round 001 should produce a reviewed v2 candidate of
`cartesian_hamiltonian_producer_design.md` plus a reconciliation note ready for
ChatGPT-Pro review.

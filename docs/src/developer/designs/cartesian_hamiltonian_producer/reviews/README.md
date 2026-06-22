# Cartesian Hamiltonian Producer Design Reviews

This directory preserves review rounds for the Cartesian Hamiltonian producer
design.

The compact current authority now lives one directory up:

- [../current.md](../current.md)
- [../registry.md](../registry.md)
- [../invariants.md](../invariants.md)
- [../implementation_slices.md](../implementation_slices.md)

The full June 2026 design document is preserved at:

- [../history/cartesian_hamiltonian_producer_design_2026-06_full.md](../history/cartesian_hamiltonian_producer_design_2026-06_full.md)

Review files are historical discussion and reconciliation artifacts. They are
not normal startup reading and do not override the compact current authority.

Historical recursive/previous-block projection material in these reviews is
stale. The live Slice A authority rejects previous-block projection and
effective-support growth for PQS terminal shell realization; use the compact
`registry.md`, `invariants.md`, and `implementation_slices.md` contract
instead.

## Historical Operating Model

Repo-manager owned this design-review lane.

Each round produced, when needed:

- a reviewer prompt;
- one consolidated review, or separate numerical/performance/minimality reviews;
- a reconciliation note;
- a revised design document when changes were accepted.

Reviewers proposed document edits only. They did not edit `src`, `test`, `bin`,
or `tools`.

ChatGPT-Pro was used as a milestone reviewer after local review and manager
reconciliation, not as a per-pass checker.

The original Round 001 milestone was to produce a reviewed v2 candidate plus a
reconciliation note. That milestone is historical; the compact current
authority now records the implemented A/B/C/D base handoff.

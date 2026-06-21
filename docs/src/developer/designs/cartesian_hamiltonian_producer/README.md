# Cartesian Hamiltonian Producer

This directory is the compact current authority for the Cartesian/PQS
Hamiltonian producer.

Status: Slice A, Slice B, Slice C1, Slice C2, and Slice D base handoff are
implemented for the internal base PQS Hamiltonian path. R1 public base
producer implementation is approved only for the narrow H/H2 scope recorded in
`r1_public_base_producer.md` and `registry.md`. Broad driver UX, Cr2
stress/performance validation, and non-base or supplement Hamiltonians remain
deferred.

Agents should read first:

- [current.md](current.md)
- [registry.md](registry.md)
- [invariants.md](invariants.md)
- [implementation_slices.md](implementation_slices.md)
- [Algorithm implementation index](../../algorithm_implementation_index.md)

Approved amendments:

- [R1 public base producer](r1_public_base_producer.md) defines the approved
  minimal public base Hamiltonian producer surface for first origin-centered H
  and z-axis H2 implementation.

Historical material:

- [history/cartesian_hamiltonian_producer_design_2026-06_full.md](history/cartesian_hamiltonian_producer_design_2026-06_full.md)
  preserves the full June 2026 design document.
- [reviews/README.md](reviews/README.md) preserves design review rounds and reconciliations.

Historical files are useful for context and audit trails, but they are not
normal startup reading and should not override the compact current authority.

Strategic context:

- [Cartesian long-range roadmap](../../roadmaps/cartesian_long_range_roadmap.md)
  sequences future public-producer, unification, supplement, high-order, and
  Cr2 validation work. It is not implementation authority.

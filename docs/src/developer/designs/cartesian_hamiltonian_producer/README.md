# Cartesian Hamiltonian Producer

This directory is the compact current authority for the Cartesian/PQS
Hamiltonian producer.

Status: Slice A, Slice B, Slice C1, Slice C2, and Slice D base handoff are
implemented for the internal base PQS Hamiltonian path. R1 public base
producer implementation is approved only for the narrow H/H2 scope recorded in
`r1_public_base_producer.md` and `registry.md`. R3-A/B/C residual-GTO/MWG
augmentation is implemented for the narrow H2 endpoint and compact artifact.
The R3 usability lane approves only a non-exported supported facade for H2 and
internal/performance-supported Be2 artifacts. The owner-local
residual-selection source correction is approved with H2 self-Coulomb
`0.4574265214362075`; broad public API, driver UX, Cr2 validation, ECP, EGOI,
RHF, and solver work remain deferred.

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
- [R3 residual-GTO/MWG augmentation](r3_residual_gto_mwg_augmentation.md)
  records the implemented R3-A/B/C residual basis, exact augmented one-body,
  MWG/IDA Hamiltonian, and compact artifact provenance path for the first H2
  endpoint.
- [R3 usability supplemented workflow](r3_usability_supplemented_workflow.md)
  approves only a non-exported supported facade that constructs H2 and
  internal/performance-supported Be2 supplemented artifacts from system, base
  basis, and supplement specs.

Candidate amendments:

- Cr2-readiness measurement, public supplemented workflow/export, and
  basis/supplement-realism lanes remain candidate-only until separately
  approved.

Historical material:

- [history/cartesian_hamiltonian_producer_design_2026-06_full.md](history/cartesian_hamiltonian_producer_design_2026-06_full.md)
  preserves the full June 2026 design document.
- [reviews/README.md](reviews/README.md) preserves design review rounds and reconciliations.

Historical files are useful for context and audit trails, but they are not
normal startup reading and should not override the compact current authority.
Recursive/previous-block projection review, spike, and performance material is
stale: terminal supports are owned local rows, cross-block overlap is
structurally zero, and post-`d2bf139c` Be2 measurements supersede earlier
recursive-projection Be2 measurements.

Strategic context:

- [Cartesian long-range roadmap](../../roadmaps/cartesian_long_range_roadmap.md)
  sequences future public-producer, unification, supplement, high-order, and
  Cr2 validation work. It is not implementation authority.

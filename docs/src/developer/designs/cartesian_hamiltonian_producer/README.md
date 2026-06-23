# Cartesian Hamiltonian Producer

This directory is the compact current authority for the Cartesian/PQS
Hamiltonian producer.

Status: Slice A, Slice B, Slice C1, Slice C2, and Slice D base handoff are
implemented for the internal base PQS Hamiltonian path. R1 public base
producer implementation is approved only for the narrow H/H2 scope recorded in
`r1_public_base_producer.md` and `registry.md`. Residual Gaussian basis,
exact augmented operators, and residual MWG/IDA interaction now belong to the
internal `CartesianResidualGaussians` module. The R3 usability lane approves
only a non-exported supported facade for H2 and internal/performance-supported
Be2 artifacts. A neutral Cartesian Gaussian raw-block owner is approved for the
uncharged by-center nuclear slice and the narrow non-nuclear
overlap/kinetic/moment slice. A narrow R3/RG terminal `G-G` product-matrix
optimization lane is approved separately. Broad public API, driver UX, Cr2
validation, ECP, EGOI, RHF, and solver work remain deferred.

Agents should read first:

- [current.md](current.md)
- [registry.md](registry.md)
- [invariants.md](invariants.md)
- [implementation_slices.md](implementation_slices.md)
- [Residual Gaussian domain module](residual_gaussian_domain_module.md) for
  current RG algorithm authority
- [Cartesian Gaussian raw blocks - nuclear slice](cartesian_gaussian_raw_blocks_nuclear.md)
  for the neutral uncharged nuclear raw-block owner
- [Cartesian Gaussian raw blocks - non-nuclear slice](cartesian_gaussian_raw_blocks_non_nuclear.md)
  for the neutral overlap/kinetic/moment raw-block owner
- [R3 terminal G-G product matrices](r3_terminal_gg_product_matrices.md)
  for the narrow final-basis `G-G` product-matrix optimization lane
- [Algorithm implementation index](../../algorithm_implementation_index.md)

Approved amendments:

- [R1 public base producer](r1_public_base_producer.md) defines the approved
  minimal public base Hamiltonian producer surface for first origin-centered H
  and z-axis H2 implementation.
- [R3 residual-GTO/MWG augmentation](r3_residual_gto_mwg_augmentation.md)
  records implementation history and compact artifact provenance for the first
  H2 endpoint. Current residual Gaussian algorithm authority is the domain
  module page below.
- [R3 usability supplemented workflow](r3_usability_supplemented_workflow.md)
  approves only a non-exported supported facade that constructs H2 and
  internal/performance-supported Be2 supplemented artifacts from system, base
  basis, and supplement specs.
- [Residual Gaussian domain module](residual_gaussian_domain_module.md)
  is the canonical current RG algorithm contract for residual basis selection,
  exact augmented operators, matched-width Gaussian descriptors, and residual
  IDA interaction blocks.
- [Cartesian Gaussian raw blocks - nuclear slice](cartesian_gaussian_raw_blocks_nuclear.md)
  approves only the neutral owner for exact uncharged by-center Gaussian
  nuclear `G-A`/`A-A` raw blocks shared by Residual Gaussian and Qiu-White
  consumers.
- [Cartesian Gaussian raw blocks - non-nuclear slice](cartesian_gaussian_raw_blocks_non_nuclear.md)
  approves only the neutral owner for exact non-nuclear Gaussian
  overlap/kinetic/coordinate-moment/second-moment `G-A`/`A-A` raw blocks.
- [R3 terminal G-G product matrices](r3_terminal_gg_product_matrices.md)
  approves only the terminal final-basis `G-G` product matrices used by
  residual-Gaussian exact augmented operators.

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

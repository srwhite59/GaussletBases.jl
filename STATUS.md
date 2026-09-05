# Current Status

This is the current capability/trust summary, not a changelog or a roadmap.
A supported interface specifies what may be constructed; it does not promise
convergence, arbitrary geometry, or a complete electronic-structure solution.

## Recommended starting point

The radial/atomic line remains the mature beginner workflow: ordinary 1D,
half-line and radial bases, quadrature and diagnostics, radial one-body
operators, atomic `(l,m)` channels, static IDA ingredients, direct/exchange/Fock
helpers, a minimal UHF kernel, and dense/sliced atomic export.
Start with README.md and the [Manual](docs/src/manual/index.md).

## Released bounded interfaces

- [Projected q-shells (PQS)](docs/src/manual/projected_q_shells.md) is the
  standard current Cartesian construction and default base route.
  `cartesian_base_hamiltonian` constructs an unsupplemented, uncorrected
  all-electron localized IDA Hamiltonian for origin-centered positive-charge
  atoms or equal-label/equal-charge homonuclear z-axis diatomics.
  Both `q` and `ns` are supported: PQS uses `q = ns`; the White-Lindsey
  alternative uses `q = ns - 2`. Supplying both requires consistency.
- [Reference-density Hartree screening](docs/src/manual/reference_density_hartree_screening.md)
  assembles a correction from a supplied field in the same orthonormal basis
  and order, with a separate energy scalar. It does not generate or fit the
  field, correct exchange, or run SCF, and is distinct from PQS construction.
- [External Cartesian GTO transfer](docs/src/manual/external_cartesian_gto_transfer.md)
  reads validated version-1 packets and preserves raw projection, including
  capture loss. Closest-determinant cleanup is separate and caller-thresholded.
  The checkpoint-only PySCF/NumPy exporter is optional; it supplies neither a
  Hamiltonian nor a complete restart and does not broaden destination geometry.

Exact Cartesian `cross_overlap`, `basis_projector`, and `transfer_orbitals`
remain supported on their working representation families. The residual-GTO
system constructor supports same-construction overlap and external import.

## Expert and research workflows

Mapped ordinary and Qiu-White residual-Gaussian routes retain their documented
backend and supplement restrictions. One-center nested fixed-block and
bond-aligned diatomic source/fixed-block/geometry workflows are supported
within their individual contracts; they are not general molecular solvers.
See the [Current ordinary branch](docs/src/explanations/current_ordinary_branch.md).

Angular injection, contraction/hierarchy, and chain/square nested producers
remain experimental. The sliced hydrogen-chain producer is a separate expert
Hamiltonian constructor, not a solver. `cartesian_base_working_basis` exposes
expert/unstable staged state; PRF/injection mechanics remain private.

## Legacy and exclusions

The old 1D COMX-cleaned hybrid route is legacy/internal, not the current
ordinary workflow. Flat supporting-note history does not override current
manuals and contracts. Broad exact four-index interaction workflows, general
molecular SCF/post-HF, and native DMRG remain outside the supported package
workflow; the optional PySCF exporter is not a general interoperability layer.

# Current Cartesian Hamiltonian Producer Authority

This page is the compact live status page for future agents. It is not the
algorithm manual. Read it first to understand what is implemented and where the
current authority lives.

Normal startup reading:

- `README.md` for orientation;
- `current.md` for live status;
- `registry.md` for approved IDs, ownership, files, and function surfaces;
- `invariants.md` for architecture-wide guardrails;
- `residual_gaussian_domain_module.md` for the canonical Residual Gaussian
  algorithm contract;
- `cartesian_gaussian_raw_blocks_nuclear.md` for the neutral uncharged nuclear
  raw-block owner;
- `docs/src/developer/algorithm_implementation_index.md` for existing kernels
  and donor paths.

Historical design and review material remains under `history/`, `reviews/`,
and the R3 amendment pages. Those files are evidence and rationale, not normal
startup authority when they conflict with the compact current files.

## Live Status

The internal base PQS Hamiltonian lane is implemented for origin-centered H and
Cartesian z-axis H2. This is internal numerical producer authority, not broad
public API polish.

Implemented base path:

- Slice A terminal basis realization over owned terminal supports;
- Slice B final-basis one-body assembly;
- Slice C localized IDA matrix assembly and existing
  `CartesianIDAHamiltonian` construction;
- Slice D base PQS materialization handoff;
- R1 public base facade for the approved H/H2 scope and fixed
  `producer_provenance/` artifact group.

Implemented Residual Gaussian path:

- Residual Gaussian basis construction, exact augmented one-body/moment
  operators, and residual MWG/IDA interaction live in
  `src/cartesian_residual_gaussians/`;
- compatibility entry points and artifact/facade hooks remain in
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`;
- the H2 owner-local residual-GTO/MWG endpoint has augmented dimension `489`
  and lowest-orbital IDA self-Coulomb `0.4574265214362075` within `1.0e-10`;
- Be2 remains an internal/performance-supported usability proxy, not a
  committed public gate;
- Cr2 residual diagnostics are not facade support and do not authorize full Cr2
  supplemented Hamiltonian or artifact production.

Approved neutral raw-block nuclear owner:

- `src/cartesian_gaussian_raw_blocks/` is approved for exact uncharged
  by-center Cartesian Gaussian nuclear `G-A` and `A-A` raw blocks only;
- `src/cartesian_gaussian_axis_integrals.jl` is approved only for the
  `HP-CGAI-FN-01` in-place axis integral table fill helper consumed by that
  neutral nuclear owner;
- the owner may be consumed by Residual Gaussian and Qiu-White code after
  behavior-preserving parity;
- it does not own overlap, kinetic, moments, terminal projection, residual
  Gaussian transforms, Qiu-White route objects, caches, reports, artifacts, or
  public API.

## Residual Gaussian Authority

`residual_gaussian_domain_module.md` is the canonical current RG algorithm
contract. Do not duplicate the full residual-selection, exact-transform, or
MWG-interaction algorithm in this file.

Essential live guardrails:

- residual Gaussian directions are selected separately on each atom and then
  merged once;
- residual occupation is a physical retained-direction measure, not numerical
  rank;
- exact augmented one-body/moment transforms are not the MWG approximation;
- MWG descriptors are rotation-dependent and must be computed from the final
  merged residual basis;
- RG does not own artifact writing, artifact provenance, basis loading, facade
  parsing, parent lattice construction, terminal topology, raw Gaussian
  nuclear block formula ownership, or public exports.

Exact uncharged by-center Gaussian nuclear raw blocks are a separate neutral
kernel authority under `cartesian_gaussian_raw_blocks_nuclear.md`. RG consumes
those blocks as exact operator inputs; RG does not own their raw construction.

The old R3-A/R3-B/R3-C labels are implementation-history labels. The current
domain names are:

- `build_residual_gaussian_basis`;
- `transform_augmented_operator`;
- `moment_matched_gaussians`;
- `assemble_residual_ida_interaction`.

## Current Boundaries

Base pair/assembly role:

- Cartesian pair terms remain a structural precursor for the base route;
- `cartesian_assembly` remains topology/report context only;
- production Hamiltonian construction must use realized terminal basis and
  final-basis matrix kernels, not source-plan/report reconstruction.

Artifact ownership:

- base Hamiltonian files use the existing `CartesianIDAHamiltonian` artifact
  shape plus the fixed `producer_provenance/` group;
- supplemented R3/RG writes may add only the approved compact
  `supplement_provenance/` group through the existing helper outside the RG
  module;
- `read_cartesian_ida_hamiltonian` ignores provenance and reads matrices
  normally.

Deferred lanes:

- public-driver polish and canonical user-facing workflow cleanup;
- Cr2 stress/performance and Cr2 support decision;
- non-base/supplement public workflow;
- ECP/EGOI/RHF/solver/HamV6 work;
- artifact/public API decisions beyond the approved compact provenance groups.

Before implementation, confirm the approved ID and source surface in
`registry.md`; if the needed surface is absent, do a docs-only amendment first.

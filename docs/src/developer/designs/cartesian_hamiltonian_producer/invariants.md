# Cartesian Hamiltonian Producer Invariants

These are architecture-wide guardrails. Detailed algorithms live in the module
contract pages, especially `residual_gaussian_domain_module.md` for Residual
Gaussian work.

## Authority

- `registry.md` owns approved IDs, files, functions, and source ownership.
- `current.md` owns live status and startup direction.
- `residual_gaussian_domain_module.md` owns the canonical Residual Gaussian
  math/algorithm contract.
- Historical design and review files are rationale unless explicitly promoted
  by the compact current authority.

## Terminal Basis

- Terminal basis blocks are supported only on their owned terminal support
  regions.
- Previous-block projection, recursive projection, and effective-support growth
  onto earlier terminal regions are forbidden.
- Cross-block overlap is zero by structure because parent rows are orthonormal
  and terminal supports are disjoint. A nonzero structural overlap means
  duplicated support rows, wrong row restriction, wrong ownership, or indexing
  error.
- Cross-block kinetic, nuclear, one-body, and IDA terms may be nonzero and must
  be assembled blockwise.
- Do not form a global parent overlap, global support operator, global dense
  coefficient matrix, global final overlap, or global Lowdin path.

## Lowdin

- Shell-local and residual merge Lowdin operations use
  `inv(sqrt(Symmetric(overlap)))`.
- Do not hand-roll eigendecomposition transforms as the production Lowdin
  implementation.
- Do not floor small eigenvalues to preserve directions that the approved
  residual or merge policy rejects.

## One-Body And IDA Assembly

- Final-basis one-body operators are assembled from terminal block pairs using
  realized block supports and coefficients.
- Nuclear attraction uses term-first Gaussian contraction and reusable 1D factor
  data where possible.
- Neutral Cartesian Gaussian raw-block nuclear kernels may own only exact
  uncharged by-center `G-A` and `A-A` raw matrices, `U_A = -1/r_A`. They do not
  apply physical nuclear charges, perform terminal projection, transform into
  residual bases, or own overlap/kinetic/moment blocks.
- Neutral Cartesian Gaussian raw-block non-nuclear kernels may own only exact
  overlap, kinetic, coordinate moment, and second-moment `G-A` and `A-A` raw
  matrices. They do not perform terminal projection, optimize final-basis
  `G-G` product matrices, transform into residual bases, or own nuclear/MWG
  interaction blocks.
- R3/RG terminal `G-G` product-matrix optimization may own only final-basis
  product matrices for kinetic, coordinate moments, and second moments used by
  residual-Gaussian exact augmented operators. It must not change terminal
  basis realization, `G-A`/`A-A` raw blocks, unit-nuclear Gaussian sums, IDA,
  MWG, residual selection, or route/public/artifact surfaces.
- Slice C/base IDA assembly produces a real final-basis matrix and then uses the
  existing `CartesianIDAHamiltonian`; no Hamiltonian wrapper or result payload is
  approved.

## Residual Gaussian Guardrails

The canonical RG algorithm contract is
`residual_gaussian_domain_module.md`. Keep only these non-negotiable rules here:

- residual basis directions are selected separately on each physical owner atom
  and then merged once;
- residual occupation is not numerical rank, an integral weight, or a tolerance
  repair knob;
- exact augmented one-body/moment transformation is not the MWG approximation;
- raw Cartesian Gaussian block construction is not Residual Gaussian basis
  selection or augmented-operator transformation;
- MWG centers and widths are computed from the final merged residual basis and
  are not invariant under arbitrary residual rotations;
- weight-aware final-basis density normalization is required for PQS-shell
  `V_GM` blocks;
- RG does not own artifact writing, artifact provenance, basis loading, facade
  parsing, public exports, driver workflow, or route-stage/report fields.

## Provenance And Metadata

- Metadata may contain provenance only. It must not carry transforms,
  coefficients, matrices, source plans, runtime inventories, or numerical data
  needed by algorithms.
- The base `producer_provenance/` group and supplemented
  `supplement_provenance/` group are compact artifact provenance, not algorithm
  inputs.
- Existing Hamiltonian readers must continue to read matrices normally and
  ignore provenance unless a later approved reader contract says otherwise.

## Public Boundaries

- R1 public base scope is origin-centered H and Cartesian z-axis H2 only.
- R3/RG supplemented usability remains internal unless a later public/export
  amendment approves it.
- Be2 is an internal performance/usability proxy, not a committed public gate.
- Cr2 diagnostics do not authorize Cr2 facade support, full augmented
  Hamiltonian construction, artifacts, or stress gates.

## Carrying Cost

- New source surfaces require an approved ID and owner in `registry.md`.
- Replacement work must switch live callers and delete obsolete helper pressure
  in the same lane unless an exact live caller blocks deletion.
- Tests should validate physics endpoints or small active module contracts, not
  stale helper names, blocker vocabulary, route metadata, or transition-only
  wrappers.

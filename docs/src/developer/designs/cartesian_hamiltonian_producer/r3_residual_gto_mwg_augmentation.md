# R3 Residual-GTO/MWG Augmentation

Status: R3-A, R3-B, and R3-C are implemented for the narrow H2
residual-GTO/MWG endpoint. R3-A provides deterministic residual-GTO basis
construction, exact augmented one-body matrices, and exact moments. R3-B
provides the same-construction in-memory MWG/IDA Hamiltonian path. R3-C
provides compact supplemented artifact provenance in the existing Hamiltonian
artifact shape. The residual selection authority is corrected and approved for
source implementation: owner-local residual occupation selection plus final
merge Lowdin replaces the implemented global candidate-order selection.
Future residual Gaussian source ownership is the
`CartesianResidualGaussians` domain module approved in
`residual_gaussian_domain_module.md`; the R3-A/B/C labels and
`pqs_terminal_residual_gto.jl` helper names are implementation history.

This document records the R3 path for generic residual-GTO/MWG augmentation of
the Cartesian base Hamiltonian producer. It approves only the R3-A, R3-B, and
R3-C IDs listed below for deterministic residual-basis construction, exact
augmented one-body and moment matrices, residual MWG/IDA interaction blocks,
the first in-memory supplemented `CartesianIDAHamiltonian{Float64}`, and
compact supplemented artifact provenance. It does not approve public API
expansion, driver or tool workflow, broad payload/status/report objects, ECP,
EGOI, Be2/Cr2 validation, or broad residual-basis serialization.

The first R3 review did not approve the IDs unchanged. This revision tightened
the residual basis contract, split the implementation lane, kept artifact
provenance compact after the in-memory numerical endpoint was accepted, and now
records the R3 closeout status before public workflow or Cr2 expansion.

## Approved IDs And Split

R3-A approved scope: residual basis plus exact one-body and moments.

- `HP-R3-OBJ-01` - residual-GTO augmentation object.
- `HP-R3-FN-01` - deterministic residual-basis construction.
- `HP-R3-FN-02` - exact augmented one-body and moment assembly.
- `HP-R3-TEST-01` - first H2 augmented one-body endpoint validation.

Approved R3-A source owner/path:

```text
Owner module: CartesianFinalBasisRealization
Source file: src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
```

`HP-R3-OBJ-01`, `HP-R3-FN-01`, and `HP-R3-FN-02` may be implemented only in
that source file. Existing `CartesianFinalBasisRealization` module include
plumbing may include the file only to expose the approved internal R3-A
surfaces; it does not approve public API or export expansion.

R3-B approved scope: MWG interaction and in-memory Hamiltonian.

- `HP-R3-FN-03` - residual MWG/IDA interaction assembly and existing
  `CartesianIDAHamiltonian{Float64}` construction.

R3-C approved scope: compact supplemented artifact provenance.

- `HP-R3-ART-01` - compact supplemented artifact provenance in the existing
  Cartesian IDA Hamiltonian artifact.

No R3 approved ID authorizes a monolithic public producer, broad
payload bundle, driver workflow expansion, status/result object, pair/assembly
workflow expansion, solver work, ECP, EGOI, Be2 validation, or Cr2 stress
gate. R3-C approves only the compact `supplement_provenance/` group described
below.

## Goal

R3 restores supplement augmentation as a generic final-basis augmentation, not
as an H2-specific residual-GTO sidecar and not as a route-stage
preflight/status graph.

The durable long-range boundary is:

```text
validated base construction specification
-> base terminal final basis G and base G-G operator blocks
+ Gaussian supplement candidates A
-> residual block R = orthogonal complement of A against fixed G
-> augmented basis B = [G, R]
-> exact augmented K, U_A, and moment matrices
-> residual MWG interaction blocks
-> CartesianIDAHamiltonian
-> optional supplemented artifact / consumer handoff
```

The base Hamiltonian and terminal basis must be produced from the same
validated construction specification. An arbitrary
`CartesianIDAHamiltonian{Float64}` of matching dimension is not sufficient
provenance for augmentation. Existing base G-G matrix blocks may be reused, but
augmentation is not a post-hoc wrapper around an opaque Hamiltonian.

Parent-axis numerical ownership invariant: once the Cartesian parent lattice
and axis bundle are realized, that construction is the authority for reusable
parent-only one-dimensional numerical data, including overlap, kinetic,
coordinate, second moment, integral weights, Gaussian factor terms, raw
pair-factor terms, and exponent ordering. Parent-supplement cross tables are
derived augmentation work data because they additionally depend on the
validated supplement, Gaussian expansion, and physical centers. R3 may derive
those tables construction-locally from the parent-axis source and reuse them
across residual overlap, exact kinetic, coordinate moments, second moments, and
by-center nuclear attraction. They are not metadata, report fields, route-stage
fields, artifacts, public API, or global mutable caches.

## Donor Inventory

| Surface | Classification | R3 use |
| --- | --- | --- |
| `CartesianGaussianShellSupplementRepresentation3D` in `src/cartesian_basis_representation.jl` | active reusable representation | Use as the explicit Gaussian candidate sector `A`. The orbital representation carries label, angular powers, center, exponents, coefficients, and primitive normalization, but no owner-nucleus index. |
| `legacy_atomic_gaussian_supplement`, `legacy_bond_aligned_diatomic_gaussian_supplement`, `legacy_bond_aligned_heteronuclear_gaussian_supplement` in `src/legacy_basis_adapter.jl` | active donor with migration naming | Use initially for named-basis loading and deterministic Cartesian shell generation. Do not freeze `legacy_*` names into new R3 public API. |
| `cpb_mixed_gto_overlap_block`, `cpb_mixed_gto_position_operator_block`, `cpb_mixed_gto_x2_operator_block`, `cpb_mixed_gto_kinetic_operator_block` in `src/CartesianCPBBlockProviders.jl` | active reusable kernels / H2-valid reference path | Use as local mixed gausslet/GTO block kernels and oracle/comparison route. The Be2 audit below shows repeated CPB-per-terminal-block construction is not the preferred R3-A scaling path. R3 should consume individual numerical blocks, not pass a vague persistent provider bundle. |
| `cpb_gto_overlap_operator_block`, `cpb_gto_position_operator_block`, `cpb_gto_x2_operator_block`, `cpb_gto_kinetic_operator_block` in `src/CartesianCPBBlockProviders.jl` | active reusable kernels / H2-valid reference path | Use for supplement self overlap, moment, and kinetic blocks when local CPB construction is the right scale; do not require this route for R3-A optimization. |
| `cpb_mixed_gto_nuclear_by_center_block`, `cpb_gto_nuclear_by_center_block`, `cpb_gto_supplement_local_operator_bundle` in `src/CartesianCPBBlockProviders.jl` | active reusable kernels / H2-valid reference path | Valid local mixed/self by-center nuclear blocks and exact moment ingredients. Not the preferred Be2-scale R3-A provider seam after the donor-kernel audit. |
| `_cartesian_supplement_cross_overlap`, `_cartesian_factorized_basis_supplement_cross`, `gto_overlap_matrix` in `src/cartesian_cross_overlap.jl` and `src/cartesian_gto_probes.jl` | active probes / reusable overlap kernels | Use for validation and possibly for small residual-space construction. Avoid making probe-only block-index and occupancy APIs the R3 production boundary. |
| `fit_radial_ylm_to_solid_harmonic_gto`, `radial_ylm_fit_cartesian_gto_adapter`, `project_cartesian_gto_to_supplement_subspace` in `src/radial_ylm_gto_bridge.jl` | active conversion/projection helpers | Useful for supplement candidate generation and small subspace projection diagnostics; not the Hamiltonian assembly boundary. |
| `src/ordinary_qw_residuals.jl` | oracle/reference plus algorithm donor | Reuse residual keep-policy, stabilization, and diagnostics ideas. Do not carry ordinary-QW route objects or thresholds into R3 authority automatically. |
| `src/ordinary_qw_raw_blocks.jl` | active donor pattern for R3-A exact block organization | Use the analytic 1D-table cross/self organization as an implementation pattern inside the approved R3 owner file. Do not edit this file, expose a new shared helper API, copy QW route objects, or pass broad provider payloads. |
| `src/ordinary_qw_operator_assembly.jl` | oracle/reference plus MWG semantic donor | Reference for moment-matched Gaussian semantics and residual-containing interaction blocks. Do not copy QW payload/status structure. |
| `src/gaussian_coulomb_reference.jl` | oracle/reference only | Use for tiny dense Gaussian Coulomb checks, not production scaling. |
| `gto_overlap_matrix`, `gto_occupancy_matrix` public probes in `src/cartesian_gto_probes.jl` | diagnostics / oracle only | Keep for user diagnostics and validation; do not route production through broad probe APIs. |
| `test/docs/cartesian_ham_builder_policy_runtests.jl` | superseded negative policy guard | The canonical-driver policy is now amended by `cartesian_driver_usability_workflow.md`. Update or delete this test during driver implementation if it conflicts with the compact functional driver, while preserving the ban on route-provider diagnostics in the canonical driver. |

Current-symbol audit, June 2026:

- no live `src`, `test`, `tools`, or `bin` definitions remain for
  `_pqs_source_box_route_driver_diatomic_physical_gausslet_supplement_*`,
  `supplement_preflight`, `:missing_provider_gto_supplement_blocks`,
  `:missing_residual_mwg_representation`,
  `:missing_combined_density_density_readiness`,
  `residual_gto_provider_blocks`, or `pqs_h2_residual_gto_handoff`;
- remaining `physical_gausslet` names in live code are route/RHF vocabulary,
  not R3 supplement-preflight wrappers, and are not R3 deletion targets unless
  a later audit finds a direct R3 caller;
- history, blurb logs, and the manager running log intentionally retain old
  supplement-preflight evidence as historical material only.

## R3-A Residual Basis Contract

R3-A constructs an augmented basis `[G, R]` and exact one-body/moment matrices.
It does not construct MWG interaction blocks, final Hamiltonians, artifacts, or
cleanup commits.

### Candidate Ownership Rule

`CartesianGaussianShellOrbitalRepresentation3D` does not carry an explicit
owner-nucleus index. For first R3-A scope, every supplement candidate center
must match exactly one physical nucleus after validated input normalization.
Owner indices are derived from that center match only. Parsing owner identity
from labels such as `"a_"` or `"b_"` is forbidden.

A later design may amend the typed representation to carry explicit owner
indices. R3-A does not require that representation change.

### Residualization Algebra

Corrected R3 residual selection is owner-local. The global raw-candidate
symmetric Lowdin and global raw-column pivoted-Cholesky selection are
superseded and are not the R3 residual algorithm for future work.

Let `G` be the fixed, orthonormal base terminal final basis and `A` the ordered
Gaussian supplement candidates. Assemble the mixed overlap one terminal block
at a time from exact mixed overlaps; do not form a global parent-space
coefficient matrix.

```text
X = G' S A
S_AA = A' S A
```

Partition candidate indices by physical owner center. For owner `a`, let
`A_a` be the candidate subblock and `X_a = G' S A_a`. The projected
owner-local residual candidates are

```text
B_a = (I - P_G) A_a
```

and their owner-local residual Gram matrix is

```text
M_a = B_a' S B_a = S_AaAa - X_a' X_a
```

Interpret `M_a` as the nonzero spectrum of a unit-candidate-occupation
residual density for that owner. Its eigenvalues are residual occupations:
they measure how much atom-local GTO content remains outside the fixed
gausslet span. They are not merely numerical-rank diagnostics.

Do not insert `S_AA^-1` to make a proper nonorthogonal occupied density. That
would turn the construction into a projector and destroy the residual
occupation information.

### Owner-Local Selection And Merge

For each owner:

```text
M_a = U_a Lambda_a U_a'
K_a = {mu : lambda_a_mu > eta_RG}
R_a0 = B_a U_a[:, K_a] Lambda_a[K_a, K_a]^(-1/2)
```

The cutoff `eta_RG` is a physical residual-occupation cutoff. It is separate
from numerical negative-eigenvalue and stabilization tolerances. Owner-local
modes below `eta_RG` are discarded because they are already represented by the
gausslet basis. Eigenvalue flooring must not be used to retain low-occupation
residual modes.

Each retained owner-local residual block satisfies:

```text
G' S R_a0 = 0
R_a0' S R_a0 = I
```

The owner-local orientation inside a retained subspace may use deterministic
natural residual modes, owner-local symmetric orthogonalization of a selected
GTO-like span, or a later approved localization rotation. The selection
criterion remains the owner-local residual occupation spectrum regardless of
that orientation choice.

Approved first orientation rule:

- diagonalize `M_a` for owner-local residual-occupation selection;
- if all owner-local modes are retained, preserve the owner candidate order and
  use owner-local symmetric Lowdin `B_a M_a^(-1/2)` so the residual block
  retains the donor-style atom-centered GTO orientation;
- if rank loss occurs, use deterministic owner-local natural residual modes
  from the retained `M_a` eigenspace, ordered by decreasing residual
  occupation with deterministic sign canonicalization;
- after either owner-local orientation, apply the final inter-owner merge
  Lowdin described below.

Concatenate owner-local residual blocks:

```text
R0 = [R_1^0, R_2^0, ...]
S_merge = R0' S R0
```

The diagonal owner blocks of `S_merge` are identity; only inter-owner overlaps
remain. Perform one final symmetric Lowdin over `S_merge`:

```text
R = R0 S_merge^(-1/2)
```

This final merge Lowdin is legitimate. It merges already selected and
orthonormal owner-local residual sectors. If `S_merge` is nearly singular
after owner-local low-occupation modes are removed, diagnose duplicated
interatomic residual directions or an inadequate occupation policy; do not use
flooring to preserve the directions.

The final raw-basis representation remains:

```text
R = G T_G + A T_A
```

`T_G` and `T_A` after the final merge are the numerical authority consumed by
exact one-body, moment, and MWG construction.

Approved cutoff and tolerance policy:

```text
eta_RG = 1.0e-8
tau_neg_abs = 1.0e-12
tau_neg_rel = 1.0e-12
tau_merge_abs = 1.0e-12
tau_merge_rel = 1.0e-12
final_orthogonality_tol = 1.0e-10
```

For each owner-local `M_a`, use
`tau_neg = max(tau_neg_abs, tau_neg_rel * max(lambda_max, 1.0))`.
Eigenvalues below `-tau_neg` are construction errors. Eigenvalues between
`-tau_neg` and `eta_RG` are discarded. Eigenvalues greater than `eta_RG` are
retained residual occupations.

For the final merge, let
`tau_merge = max(tau_merge_abs, tau_merge_rel * max(lambda_max(S_merge), 1.0))`.
Any merge eigenvalue below `-tau_merge` is a construction error. Any merge
eigenvalue `<= tau_merge` is a near-singular merge error. Do not floor merge
eigenvalues to preserve directions. After final merge, require `G' S R` and
`R' S R - I` errors below `final_orthogonality_tol`.

Measurement evidence for this approval:

- H2, Be2, Cr2 q4, and Cr2 q5 pass owner-local selection with final
  orthogonality below `1.0e-10`;
- no rank loss was observed under trial residual-occupation cutoffs
  `1.0e-8` or `1.0e-7`;
- the corrected H2 lowest-orbital IDA self-Coulomb scalar is
  `0.4574265214362075`.

### HP-R3-OBJ-01 Approved Fields

The R3 residual object is a numerical object, not a status/result payload. It
does not store matrices in metadata.

Approved exact fields:

| Field | Contract |
| --- | --- |
| `base_dimension::Int` | number of base final functions `n_G` |
| `candidate_count::Int` | number of raw supplement candidates `n_A` |
| `residual_dimension::Int` | retained residual rank `n_R` |
| `candidate_labels::Vector{String}` | deterministic raw candidate labels |
| `candidate_owner_indices::Vector{Int}` | derived from exact center-to-nucleus match |
| `candidate_centers::Vector{NTuple{3,Float64}}` | raw candidate centers in candidate order |
| `residual_source_owner_indices::Vector{Int}` | owner center for each final residual function |
| `residual_occupations::Vector{Float64}` | retained owner-local residual occupations associated with final residual functions before final merge |
| `owner_retained_counts::Vector{Int}` | retained residual count by owner |
| `residual_labels::Vector{String}` | deterministic residual labels derived from owner and local retained mode order |
| `T_G::Matrix{Float64}` | shape `n_G x n_R`; base-side residual transform |
| `T_A::Matrix{Float64}` | shape `n_A x n_R`; supplement-side residual transform |
| `occupation_cutoff::Float64` | approved residual-occupation cutoff `eta_RG = 1.0e-8` |
| `tau_neg_abs::Float64` | absolute negative-eigenvalue error tolerance |
| `tau_neg_rel::Float64` | relative negative-eigenvalue error tolerance |
| `tau_merge_abs::Float64` | absolute final-merge near-singular threshold |
| `tau_merge_rel::Float64` | relative final-merge near-singular threshold |
| `selection_rule::Symbol` | owner-local residual occupation selection rule |
| `orientation::Symbol` | owner-local orientation plus final merge convention |
| `sign_rule::Symbol` | `:largest_T_A_entry_positive` |

Invariants:

- `size(T_G) == (base_dimension, residual_dimension)`;
- `size(T_A) == (candidate_count, residual_dimension)`;
- `length(candidate_labels) == candidate_count`;
- `length(candidate_owner_indices) == candidate_count`;
- `length(candidate_centers) == candidate_count`;
- `length(residual_source_owner_indices) == residual_dimension`;
- `length(residual_occupations) == residual_dimension`;
- `sum(owner_retained_counts) == residual_dimension`;
- retained `residual_occupations` are greater than `occupation_cutoff`;
- `G' S R` is zero within tolerance;
- `R' S R` is identity within tolerance.

`retained_candidate_indices` is not part of the corrected persistent object
contract because a retained residual mode may be a linear combination of
several atom-local candidates. Full discarded owner spectra belong in ignored
diagnostics unless a later consumer-driven artifact amendment promotes them.

## R3-A Exact One-Body And Moment Assembly

R3-A must assemble exact augmented operators by transforming raw `[G, A]`
blocks. For any exact one-body operator `O`, including kinetic, every
uncharged by-center nuclear attraction, Cartesian position, and Cartesian
second moment:

```text
O_raw = [O_GG  O_GA
         O_AG  O_AA]
```

Using `T_G = -X T` and `T_A = T`:

```text
O_GR = O_GG T_G + O_GA T_A
```

```text
O_RR =
    T_G' O_GG T_G
  + T_G' O_GA T_A
  + T_A' O_AG T_G
  + T_A' O_AA T_A
```

The augmented operator is:

```text
O_aug = [O_GG  O_GR
         O_GR' O_RR]
```

This applies to:

- `K`;
- every uncharged `U_A`;
- `x`, `y`, `z`;
- `x^2`, `y^2`, `z^2`.

The moment matrices are needed for MWG construction and are not contained in
`CartesianIDAHamiltonian`. R3-A must therefore keep exact moment matrices in
the R3-A construction boundary, not recover them from a base Hamiltonian.

Be2 donor-kernel audit clarification: `HP-R3-FN-02` does not require CPB
providers as the only exact-block implementation. The ignored
`tmp/work/be2_r3a_qw_cross_donor_audit.jl` comparison showed that the QW
analytic 1D-table cross/self organization matches the current CPB-provider
`G-A` and `A-A` blocks at roundoff on the Be2 proxy while avoiding repeated
CPB-per-terminal-block work. R3-A may therefore build full
parent-by-supplement `G-A` matrices once for overlap, kinetic, position,
second moments, and by-center nuclear attraction, project the parent rows
through terminal blocks, reuse the once-built overlap `G-A` block for
`X = G' S A`, and use the once-built `G-A`/`A-A` blocks in augmented operator
assembly.

This allowance was the tactical R3-local bridge for the first exact R3-A
operator implementation. A full parent-by-supplement cross matrix is allowed
because it is rectangular cross data, not a forbidden parent-by-parent global
operator or dense global pair matrix. For uncharged by-center nuclear
`G-A`/`A-A` raw blocks, the later `HP-CGRB-*` authority supersedes the
R3-local ownership story and approves a neutral raw-block owner plus
behavior-preserving Residual Gaussian and Qiu-White rewiring. That later
authority still does not approve a broad QW API, persistent provider bundles,
public API, artifact/provenance, driver/tool workflow, R3-C expansion, Be2 as a
committed validation gate, or Cr2 facade/artifact support.

The production target is to derive parent-supplement 1D cross tables once per
augmentation construction and reuse them. The immediate R3-local bridge may
tolerate one duplicate overlap construction between residualization and full
exact-operator assembly to avoid adding a persistent raw-block bundle, but it
must not rebuild cross tables per terminal block or per operator.

R3-A validation gates:

- candidate ownership is derived by exact center-to-nucleus match;
- base H2 endpoint facts still hold before augmentation;
- candidate count/order is exactly the frozen H2 fixture below;
- residual rank is measured and reported, not assumed;
- `G' S R` is near zero;
- `R' S R` is near identity;
- base `G-G` blocks equal the base construction blocks within tolerance;
- augmented `K`, every `U_A`, and all moment matrices are finite and symmetric;
- first one-body eigenvalue is variational:
  `E1_aug <= E1_base + epsilon`.
- no Be2 or Cr2 validation is part of the first R3-A gate.

## R3-B MWG/IDA And In-Memory Hamiltonian

R3-B starts only after R3-A is accepted. It constructs residual MWG descriptors,
residual-containing IDA blocks, and the in-memory augmented
`CartesianIDAHamiltonian{Float64}`.

R3-B was first reapproved with a corrected weight-aware compact-path scalar for
the then-current global candidate-order residual basis. That scalar is
superseded for corrected R3 work because the residual basis construction
itself is superseded. Under the approved owner-local residual-selection
convention, the H2 lowest-orbital IDA self-Coulomb scalar is
`0.4574265214362075`.

Previously superseded scalar history:

- `0.457435475059184` came from a retired private
  `[pre_final_pqs, residual_gto]` density-gauge diagnostic;
- `0.4574331709135599` came from direct parent-density-normalized `G-M`
  insertion into the final basis;
- `0.4574256036192161` came from the global candidate-order R3-A residual
  basis with the corrected weight-aware `G-M` contraction.

The owner-local residual-selection correction may change the MWG interaction
because MWG is not invariant under rotations of the residual subspace. Exact
augmented one-body eigenvalues should remain invariant when the retained span
is unchanged. Do not add a residual width scale factor and do not relax
tolerance to preserve any stale scalar.

Approved R3-B source owner/path/function:

```text
Owner module: CartesianFinalBasisRealization
Source file: src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
Function: pqs_terminal_residual_gto_augmented_hamiltonian
```

`HP-R3-FN-03` may be implemented only in that file. The function is internal
module surface, not a public API or export. The approved R3-B construction
shape is same-construction: callers supply the base Hamiltonian, terminal
realization, bundle/axis source, supplement, atom locations, and nuclear
charges from the same base construction; the function constructs the accepted
R3-A residual object, exact augmented one-body/moment matrices, residual MWG
descriptors, `V_GM`, `V_MM`, and the existing
`CartesianIDAHamiltonian{Float64}` inside one call.

Conceptual boundary:

```text
pqs_terminal_residual_gto_augmented_hamiltonian(
    base_hamiltonian,
    basis::CartesianTerminalBasisRealization,
    bundles,
    supplement,
    atom_locations,
    nuclear_charges;
    expansion = nothing,
)::CartesianIDAHamiltonian{Float64}
```

The exact Julia signature may be adjusted to match local types and naming. The
contract is that callers no longer pass independently constructed residual or
augmented-operator objects into the internal R3-B boundary. Existing
lower-level R3-A and R3-B helpers may remain and may be reused; this amendment
does not require deleting or replacing them. File-local helpers are allowed
only for residual construction, exact augmented operator construction,
residual MWG descriptor extraction, `G-M` block assembly, `M-M` block
assembly, final `V_aug` construction, and endpoint validation.

The same-construction path may recompute or locally reuse the one-shot
parent-by-supplement exact-block family. It must not add a persistent
raw-block bundle/cache object. One duplicate overlap build between
residualization and full exact-operator assembly remains acceptable unless it
can be removed locally without a new persistent shape.

### MWG Moment Convention

For normalized residual `r` and Cartesian axis `alpha`:

```text
c_ralpha = <r | alpha | r>
v_ralpha = <r | alpha^2 | r> - c_ralpha^2
sigma_ralpha = sqrt(2 * v_ralpha)
```

The factor of two is required by the repo's Gaussian width convention.

Rules:

- moments are computed from the actual final residual functions using the
  exact transformed R3-A moment matrices;
- the effective MWG function is a separable Gaussian with center `c_r` and
  axis widths `sigma_r`;
- off-diagonal covariance is deliberately omitted as part of the separable MWG
  approximation;
- `v_ralpha <= 0`, nonfinite moments, or nonfinite/nonpositive widths are
  construction errors;
- residual unsquared integrals are not used as base-PQS IDA weights.

### Interaction Convention

Let `M` denote the moment-matched effective Gaussians for residual functions,
not the original supplement candidates and not exact residual-GTO Coulomb
integrals.

```text
V_aug = [V_GG_base  V_GM
         V_GM'      V_MM]
```

Normalization convention: R3-B consumes density-normalized donor factors for
`M-M` pair blocks directly. For `G-M`, donor support-to-M factors are
parent-density normalized and must be transformed to final-basis density
normalization at each terminal block:

```text
support_weights = wx .* wy .* wz
final_weights = C' * support_weights
C_density = C .* support_weights ./ final_weights'
V_GM_block = C_density' * V_support_M
```

Direct blocks use identity/final weights consistently and therefore agree with
direct insertion. PQS shell blocks must use the weight-aware contraction above.
After this contraction, the final augmented `V` assembly performs no additional
division by base final weights or effective-Gaussian weights.

`V_GG_base` is reused from the accepted base localized IDA path. `V_GM` and
`V_MM` are MWG approximations, not exact residual-GTO Coulomb integrals.

The in-memory Hamiltonian is constructed as:

```text
CartesianIDAHamiltonian(
    K_aug,
    U_aug_by_center,
    V_aug,
    base_hamiltonian.nup,
    base_hamiltonian.ndn;
    nuclear_charges = base_hamiltonian.nuclear_charges,
    nuclear_positions = base_hamiltonian.nuclear_positions,
)
```

`K_aug` and `U_aug_by_center` must come from the accepted R3-A augmented
one-body construction. The `G-G` blocks of `K_aug`, every `U_A`, and `V_aug`
must match the supplied base Hamiltonian. The base Hamiltonian is allowed only
as the same-construction base object for the validated H2 fixture; accepting an
arbitrary dimension-compatible Hamiltonian is not an R3-B contract.

The same-construction function must reject inconsistent inputs with clear
errors rather than status/result payloads. At minimum, `base_hamiltonian`,
`basis`, and `bundles` must correspond to the same base construction, and
`atom_locations` / `nuclear_charges` must match the physical centers used by
that construction. The design does not require a new provenance object or
parent-stage field to prove this; local numerical and dimensional consistency
checks are allowed.

### Fast Coulomb Organization

R3-B must:

- build axis pair-factor tensors once outside terminal-block loops;
- retain the Coulomb Gaussian expansion term as the short inner reduction;
- contract `G-M` blockwise through terminal basis blocks;
- build `M-M` directly as an `n_R x n_R` term-first contraction;
- avoid materializing one full 3D matrix per Gaussian term;
- avoid persistent `nterms x support x residual` workspaces;
- preserve the existing bounded-workspace policy.

R3-B validation gates:

- residual centers and widths are finite and widths are positive;
- `V_GM` and `V_MM` are finite;
- `V_GM` matches an independent weight-aware final-basis density-normalized
  check built from support weights, final weights, and support-to-M donor
  values; matching only the final self-Coulomb scalar is not sufficient;
- `V_aug` is symmetric;
- base `V_GG` block is unchanged;
- returned object is the existing `CartesianIDAHamiltonian{Float64}`;
- augmented dimension is `489` for the first H2 fixture;
- the lowest augmented one-body orbital has IDA self-Coulomb
  `0.4574265214362075` within `1.0e-10`;
- no Hamiltonian wrapper/result/status payload is introduced.

The existing standalone H2 R3 endpoint gate remains the validation surface. It
must exercise the same-construction path, keep the returned augmented
dimension `489`, keep the self-Coulomb value above, and keep an independent
weight-aware `V_GM` check either in that standalone test or in ignored
validation. This amendment does not approve a new test file.

The validation scalar above comes from the weight-aware compact-path R3-B
equations after owner-local residual selection and final inter-owner merge
Lowdin: current R3-A exact residual moments, `sigma = sqrt(2v)`, final-basis
density-normalized `G-M` contraction, direct density-normalized `M-M`
contraction, and no final-width or tolerance fitting. Old private provider
payloads, artifact fields, driver diagnostics, the retired pre-final density
gauge, the direct parent-density `G-M` scalar, and the global candidate-order
R3-A scalar are donor history only, not production authority.

The full R3-A moment matrices `x`/`y`/`z`/`x^2`/`y^2`/`z^2` remain approved
outputs of the H2 construction boundary. This hardening note does not replace
them with residual-diagonal-only moments.

## R3-C Artifact And Provenance

`HP-R3-ART-01` is approved. The existing
`CartesianIDAHamiltonian{Float64}` is sufficient as the in-memory numeric
Hamiltonian object when the augmented matrices satisfy the existing contract.
The artifact remains the existing `write_cartesian_ida_hamiltonian` file shape
plus a compact `supplement_provenance/` group.

The supplemented artifact schema is compact. It must derive
provenance from the validated R3 construction specification, not recover
`base_route` or input data from an in-memory Hamiltonian.

Approved compact provenance keys:

| Key | Value |
| --- | --- |
| `supplement_provenance/provenance_version` | `1` |
| `supplement_provenance/producer` | `:cartesian_residual_gto_mwg_augmentation` |
| `supplement_provenance/supplement_policy` | `:mwg_residual_gto` |
| `supplement_provenance/basis_by_center` | center-indexed basis labels |
| `supplement_provenance/lmax` | `Int` |
| `supplement_provenance/uncontracted` | `Bool` |
| `supplement_provenance/width_filtering` | `nothing` or explicit width policy |
| `supplement_provenance/candidate_count` | `Int` |
| `supplement_provenance/owner_counts` | center-indexed candidate counts |
| `supplement_provenance/base_dimension` | `Int` |
| `supplement_provenance/residual_dimension` | retained residual rank |
| `supplement_provenance/augmented_dimension` | final Hamiltonian dimension |
| `supplement_provenance/augmented_basis_order` | `:base_then_residual` |
| `supplement_provenance/residual_basis_convention` | `:owner_local_residual_occupation_final_merge_lowdin` |
| `supplement_provenance/rank_rule` | owner-local residual occupation selector |
| `supplement_provenance/occupation_cutoff` | `1.0e-8` |
| `supplement_provenance/tau_neg_abs` | `Float64` |
| `supplement_provenance/tau_neg_rel` | `Float64` |
| `supplement_provenance/tau_merge_abs` | `Float64` |
| `supplement_provenance/tau_merge_rel` | `Float64` |
| `supplement_provenance/mwg_convention_version` | `1` |
| `supplement_provenance/mwg_convention` | `:separable_moment_matched_density_normalized` |
| `supplement_provenance/one_body_source` | `:exact_transformed_raw_blocks` |
| `supplement_provenance/interaction_source` | `:weight_aware_residual_mwg_ida_blocks` |
| `supplement_provenance/validation_check_labels` | compact symbols naming accepted checks, including `:h2_lowest_augmented_one_body_orbital_ida_self_coulomb` when writing the H2 validation fixture |
| `supplement_provenance/h2_self_coulomb_reference` | `0.4574265214362075` for the validation fixture, otherwise `nothing` |

The compact schema must not store:

- the full residual eigenvalue vector;
- MWG center matrices;
- MWG width matrices;
- `T_G` or `T_A`;
- candidate labels;
- dense moment matrices;
- full construction inputs.

If those arrays become consumer-critical, a later explicit residual-basis
artifact group must promote them together: residual transforms, candidate
labels/order, owner indices, retained candidate indices, eigenvalue
diagnostics, centers, and widths. Calling only centers and widths "provenance"
would be misleading.

Approved implementation surface:

- `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` owns the
  R3-C supplemented artifact write path.
- It may call the existing `write_cartesian_ida_hamiltonian` and then add
  `supplement_provenance/` keys to the same JLD2 file, following the R1
  provenance pattern.
- No edit to `src/cartesian_ida_hamiltonian.jl` is approved by this amendment.
  If a future implementation cannot add the provenance group from the R3 owner
  file, the exact writer seam must return for a separate docs-only amendment.

Validation gate:

- H2 supplemented artifact write/readback;
- readback Hamiltonian `K`, every uncharged `U_A`, and `V` match the returned
  in-memory Hamiltonian exactly or within tight roundoff deltas;
- provenance keys match the validated R3 construction specification;
- R3-B lowest-orbital IDA self-Coulomb matches
  `0.4574265214362075` within `1.0e-10`.

## First Validation Fixture

First implementation proxy: z-axis H2 with the existing public/base H2
geometry and the frozen contracted two-center H/cc-pVTZ supplement fixture.

Fixture:

```text
system: z-axis H2
placement: same H/cc-pVTZ atomic source on both physical H centers
lmax: 1
uncontracted: false
max_width filtering: none
candidate ordering: center-major, then source-shell/contraction order,
                    then Cartesian component order
candidate count: 18 total
owner counts: 9 on center 1 and 9 on center 2
owner rule: exact candidate-center match to one physical nucleus
```

`lmax = 1` alone does not imply 18 candidates. The exact count is a property
of this contracted cc-pVTZ fixture and must be asserted directly. The retained
residual rank must be measured.

R3-A first endpoint is H2 augmented one-body and moments, not a
residual-basis-only scaffold commit. R3-B adds MWG interaction and the existing
in-memory Hamiltonian for the same H2 fixture, using the corrected
weight-aware compact-path self-Coulomb baseline. R3-C adds compact artifact
provenance for the same supplemented Hamiltonian shape.

Cr2 remains a later stress/consumer-readiness milestone, not the first R3
correctness gate.

## R3 Closeout Status

Implemented narrow R3 scope:

- R3-A: residual-GTO augmentation object, deterministic residual construction,
  deterministic rank-loss handling, exact augmented `K`, uncharged `U_A`, and
  exact moment matrices `x`/`y`/`z`/`x^2`/`y^2`/`z^2`;
- R3-B: same-construction internal augmented Hamiltonian path, residual
  moment-matched Gaussian descriptors, weight-aware final-basis
  density-normalized `V_GM`, direct density-normalized `V_MM`, unchanged base
  `V_GG`, and existing `CartesianIDAHamiltonian{Float64}`;
- R3-C: compact `supplement_provenance/` group added to the existing
  Hamiltonian artifact shape, with validation-only readback.

Legacy H2 closeout facts before the owner-local selection correction:

- augmented dimension `489`;
- lowest augmented one-body orbital IDA self-Coulomb
  `0.4574256036192161` within `1.0e-10` for the global candidate-order
  residual basis;
- tracked standalone endpoint exercises the same-construction path and
  independent weight-aware `V_GM` check;
- H2 supplemented artifact write/readback preserves Hamiltonian matrices within
  tight deltas and stores compact provenance from the validated construction
  specification.

These facts remain useful regression evidence for the old compact path, but
the old self-Coulomb scalar and residual-basis provenance must not be
preserved as targets after owner-local residual selection and final merge
Lowdin are implemented. The corrected H2 self-Coulomb scalar is
`0.4574265214362075`.

Be2 measurement conclusions:

- The R3-A exact-operator bottleneck was repeated CPB-per-terminal-block
  mixed/self block construction, not residualization math. On the Be2 proxy it
  measured about `43.2 s` and `35.4 GiB`.
- The one-shot parent-by-supplement analytic block organization measured about
  `1.94 s` and `2.1 GiB` on the same proxy, with roundoff agreement for tested
  overlap, kinetic, coordinate, second-moment, and by-center nuclear blocks.
- Be2 R3-B MWG/IDA storage and runtime were modest at residual rank `26`.
  Bounded/streamed MWG term storage remains a high-rank and Cr2 guardrail, but
  is not urgent before choosing the next lane at this proxy scale.

Remaining deferred hardening:

- Cr2-readiness measurement only: candidate count, retained rank, memory, and
  time forecast before any full Cr2 Hamiltonian stress run;
- high-rank residual MWG streaming if rank growth makes the current dense
  residual term storage costly;
- nonallocating or bounded-allocation validation reductions for large dense
  matrices;
- public or supported internal supplemented workflow for H2/Be2 artifacts;
- basis and supplement realism beyond the first H2 contracted cc-pVTZ fixture.

Follow-on lane: `r3_usability_supplemented_workflow.md` approves a
non-exported supported facade for H2 and internal/performance-supported Be2
supplemented artifacts so consumers can request GTO/MWG artifacts without
composing private R3 calls. A Cr2-readiness lane should stay measurement-only
until that workflow and its input policy are exercised.

## Owner-Local Selection Measurement Pass

The owner-local measurement pass is closed for source-approval purposes. It
reported H2, Be2, Cr2 q4, and Cr2 q5 final orthogonality below `1.0e-10`, no
rank loss at trial residual-occupation cutoffs `1.0e-8` or `1.0e-7`, and H2
self-Coulomb `0.4574265214362075`.

Future diagnostics may live under ignored `tmp/work`; they must not edit
source, tests, tools, drivers, artifacts, or public APIs.

Use the existing `X`, `S_AA`, and `candidate_owner_indices` data. For each
owner:

```text
M_a = S_AA[I_a, I_a] - X[:, I_a]' X[:, I_a]
diagonalize M_a
report residual occupation spectrum
retain trial modes above eta_RG
normalize retained owner-local modes
```

Then concatenate owner blocks, form `S_merge`, report its eigenvalue range and
condition, apply the final symmetric Lowdin, and report final orthogonality.

Diagnostics should report:

- atom-local residual occupation spectra separately;
- retained counts per owner under each trial cutoff;
- owner-local orthogonality errors;
- `S_merge` eigenvalue range and condition;
- final merged `G' S R` and `R' S R` errors.

Forbidden in the measurement pass:

- global raw-candidate Lowdin as residual-content selection;
- global raw-column pivoted-Cholesky selection as residual-content selection;
- eigenvalue flooring to retain tiny residual occupations;
- width filtering as conditioning repair;
- Cr2 facade expansion or full Hamiltonian construction;
- operator, MWG, artifact, public API, driver, bin, or tool work.

The existing `_qwrg_residual_space_by_owner` owner-local organization is the
primary donor. Do not create a new residual framework unless measurement shows
the donor cannot express the required owner-local selection and final merge.

## R3-A Approval Evidence

Manager log Pass 048 records the old measurement-only R3-A residual-spectrum
spike for the superseded global candidate-order construction. The owner-local
approval above supersedes it for future source work.

Evidence:

- fixture: public/base z-axis H2;
- supplement: contracted H/cc-pVTZ on both physical H centers;
- `lmax = 1`;
- `uncontracted = false`;
- no width filtering;
- 18 supplement candidates total;
- 9 candidates per H center;
- deterministic center-major/source-shell/Cartesian-component order;
- base final dimension `471`;
- full residual rank `18`;
- residual metric eigenvalue range approximately `3.05e-4` to `1.35e-2`;
- condition estimate approximately `44.3`;
- `inv(sqrt(Symmetric(S_R)))` in candidate order succeeded for the superseded
  global residual construction;
- measured `R' S R` identity error approximately `4.4e-15`.

This evidence explains why the first H2 global construction was numerically
well behaved. It is no longer residual-selection authority for Cr-like
supplements or future R3 implementation.

## Deletion And Retirement Targets

Current live-code audit found no active `src`, `test`, `tools`, or `bin`
definitions for the old supplement-preflight/provider blocker family named in
earlier plans. Do not preserve those names as R3 authority.

After R3-A:

- delete any reintroduced supplement-preflight fixture or wrapper whose only
  purpose is to announce missing overlap, one-body, or moment provider blocks;
- keep named-basis loading, Gaussian representation, and CPB provider kernels;
- update or delete `test/docs/cartesian_ham_builder_policy_runtests.jl` under
  the driver usability lane if it conflicts with the amended canonical-driver
  policy.

After R3-B:

- retire duplicate H2-specific residual/MWG construction paths if parity shows
  they no longer have a named consumer;
- retain `src/ordinary_qw_residuals.jl`,
  `src/ordinary_qw_raw_blocks.jl`, and
  `src/ordinary_qw_operator_assembly.jl` as donor/oracle code until active
  callers and parity obligations are reviewed;
- retain `src/ordinary_cartesian_ida.jl` hybrid residual MWG primitives until a
  focused audit classifies them.

After R3-C:

- remove blocked artifact/status fields superseded by the real supplemented
  artifact, if any remain in live code at that time;
- remove historical docs only when explicitly doing documentation
  reorganization; history and blurb logs are not deletion targets by default.

## Implemented / Corrected First Target

The implemented first target is R3-A plus the narrow R3-B in-memory
interaction continuation and R3-C compact artifact provenance. The corrected
source target replaces the legacy residual block with owner-local
residual-occupation selection:

```text
base z-axis H2
-> contracted two-center H/cc-pVTZ lmax-1 candidates
-> owner-local residual-occupation selection and final merge Lowdin
-> exact augmented K, U_A, x/y/z, and x2/y2/z2
-> H2 augmented one-body endpoint
-> residual moment-matched Gaussians with sigma = sqrt(2v)
-> weight-aware V_GM and direct V_MM MWG/IDA blocks
-> in-memory supplemented CartesianIDAHamiltonian
-> H2 lowest-orbital IDA self-Coulomb endpoint
-> existing Hamiltonian artifact plus compact supplement_provenance/
```

The first committed/standalone endpoint gate is approved by `HP-R3-TEST-01`
in `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`. R3-A checks
the H2 augmented one-body/moment endpoint. R3-B may extend that same standalone
file, despite its R3-A name, only with the first in-memory supplemented
Hamiltonian checks described above. It must check `G' S R`, `R' S R`, base G-G
block equality, finite/symmetric augmented `K`, uncharged `U_A`, and moment
matrices, `E1_aug <= E1_base + epsilon`, finite/symmetric `V_aug`, unchanged
base `V_GG`, independent weight-aware `V_GM`, returned
`CartesianIDAHamiltonian{Float64}`, augmented dimension `489`, and the
owner-local lowest-orbital IDA self-Coulomb `0.4574265214362075` within
`1.0e-10`. It must not assert private pair/assembly/report/status behavior,
add the standalone file to `test/runtests.jl`, or run Be2 or Cr2.

Do not start Cr2 stress tests from this R3 closeout. The immediate follow-on
usability lane is recorded in `r3_usability_supplemented_workflow.md`; public
export, Cr2-readiness, and basis/supplement-realism lanes still require
separate amendments.

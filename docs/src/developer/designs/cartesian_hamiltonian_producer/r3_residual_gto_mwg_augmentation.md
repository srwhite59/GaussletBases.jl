# R3 Residual-GTO/MWG Augmentation

Status: R3-A approved; R3-B and R3-C remain candidate-only.

This document records the R3 path for generic residual-GTO/MWG augmentation of
the Cartesian base Hamiltonian producer. It approves only the R3-A IDs listed
below for deterministic residual-basis construction plus exact augmented
one-body and moment matrices. It does not approve R3-B MWG/IDA, supplemented
`CartesianIDAHamiltonian` construction, artifacts, public API expansion, driver
or tool workflow, broad payload/status/report objects, ECP, EGOI, Be2/Cr2
first-gate validation, or implementation of the R3-B/R3-C candidate IDs.

The first R3 review did not approve the IDs unchanged. This revision tightens
the residual basis contract, splits the implementation lane, and defers
artifact provenance until the in-memory numerical endpoint exists.

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

R3-B candidate-only scope: MWG interaction and in-memory Hamiltonian.

- `HP-R3-FN-03` - residual MWG/IDA interaction assembly and existing
  `CartesianIDAHamiltonian{Float64}` construction.

R3-C candidate-only scope: artifact provenance and cleanup.

- `HP-R3-ART-01` - compact supplemented artifact provenance and retirement
  cleanup after R3-A/B numerical acceptance.

No R3-A approved ID authorizes a monolithic public producer, broad payload
bundle, driver workflow expansion, status/result object, pair/assembly workflow
expansion, solver work, ECP, EGOI, MWG/IDA, supplemented Hamiltonian
construction, artifact provenance, Be2 first-gate validation, or Cr2 stress
gate.

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

## Donor Inventory

| Surface | Classification | R3 use |
| --- | --- | --- |
| `CartesianGaussianShellSupplementRepresentation3D` in `src/cartesian_basis_representation.jl` | active reusable representation | Use as the explicit Gaussian candidate sector `A`. The orbital representation carries label, angular powers, center, exponents, coefficients, and primitive normalization, but no owner-nucleus index. |
| `legacy_atomic_gaussian_supplement`, `legacy_bond_aligned_diatomic_gaussian_supplement`, `legacy_bond_aligned_heteronuclear_gaussian_supplement` in `src/legacy_basis_adapter.jl` | active donor with migration naming | Use initially for named-basis loading and deterministic Cartesian shell generation. Do not freeze `legacy_*` names into new R3 public API. |
| `cpb_mixed_gto_overlap_block`, `cpb_mixed_gto_position_operator_block`, `cpb_mixed_gto_x2_operator_block`, `cpb_mixed_gto_kinetic_operator_block` in `src/CartesianCPBBlockProviders.jl` | active reusable kernels | Use as local mixed gausslet/GTO block kernels. R3 should consume individual numerical blocks, not pass a vague persistent provider bundle. |
| `cpb_gto_overlap_operator_block`, `cpb_gto_position_operator_block`, `cpb_gto_x2_operator_block`, `cpb_gto_kinetic_operator_block` in `src/CartesianCPBBlockProviders.jl` | active reusable kernels | Use for supplement self overlap, moment, and kinetic blocks. |
| `cpb_mixed_gto_nuclear_by_center_block`, `cpb_gto_nuclear_by_center_block`, `cpb_gto_supplement_local_operator_bundle` in `src/CartesianCPBBlockProviders.jl` | active reusable kernels | Preferred first provider seam for local mixed/self by-center nuclear blocks and exact moment ingredients. |
| `_cartesian_supplement_cross_overlap`, `_cartesian_factorized_basis_supplement_cross`, `gto_overlap_matrix` in `src/cartesian_cross_overlap.jl` and `src/cartesian_gto_probes.jl` | active probes / reusable overlap kernels | Use for validation and possibly for small residual-space construction. Avoid making probe-only block-index and occupancy APIs the R3 production boundary. |
| `fit_radial_ylm_to_solid_harmonic_gto`, `radial_ylm_fit_cartesian_gto_adapter`, `project_cartesian_gto_to_supplement_subspace` in `src/radial_ylm_gto_bridge.jl` | active conversion/projection helpers | Useful for supplement candidate generation and small subspace projection diagnostics; not the Hamiltonian assembly boundary. |
| `src/ordinary_qw_residuals.jl` | oracle/reference plus algorithm donor | Reuse residual keep-policy, stabilization, and diagnostics ideas. Do not carry ordinary-QW route objects or thresholds into R3 authority automatically. |
| `src/ordinary_qw_raw_blocks.jl` | oracle/reference | Reference for raw GTO/PQS block organization and supplement plumbing. Prefer current CPB provider kernels for R3 implementation. |
| `src/ordinary_qw_operator_assembly.jl` | oracle/reference plus MWG semantic donor | Reference for moment-matched Gaussian semantics and residual-containing interaction blocks. Do not copy QW payload/status structure. |
| `src/gaussian_coulomb_reference.jl` | oracle/reference only | Use for tiny dense Gaussian Coulomb checks, not production scaling. |
| `gto_overlap_matrix`, `gto_occupancy_matrix` public probes in `src/cartesian_gto_probes.jl` | diagnostics / oracle only | Keep for user diagnostics and validation; do not route production through broad probe APIs. |
| `test/docs/cartesian_ham_builder_policy_runtests.jl` | negative policy guard | Currently asserts the canonical driver does not expose `residual_gto_provider_blocks`. Retain until a later public-driver design changes that policy. |

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

Let `G` be the fixed, orthonormal base terminal final basis and `A` the ordered
Gaussian supplement candidates. Assemble the mixed overlap one terminal block
at a time from exact mixed overlaps; do not form a global parent-space
coefficient matrix.

```text
X = G' S A
S_AA = A' S A
S_R = S_AA - X'X
A_perp = A - G X
```

For a supplement-side residual transform `T`,

```text
R = (A - G X) T
```

and the raw-basis representation is:

```text
T_A = T
T_G = -X T
R = G T_G + A T_A
```

The transform must satisfy:

```text
T' S_R T = I
```

### Orientation And Rank Policy

Residual basis orientation matters because each final residual function is
later replaced by an effective moment-matched Gaussian. R3 must not use raw
eigenvectors of `S_R` as the residual functions.

Deterministic candidate convention:

- preserve candidate ordering from the validated supplement fixture;
- project all candidates against `G`;
- use residual eigenvalues as rank diagnostics, not as final orientation;
- if all candidate directions are retained, use symmetric Lowdin in candidate
  order: `T = S_R^(-1/2)`;
- if rank loss occurs, select independent candidate columns with a
  deterministic rank-revealing rule before Lowdin, then apply symmetric Lowdin
  to the selected projected candidates;
- tie breaks in any rank-revealing step use lower original candidate index;
- selected candidates are returned to ascending original candidate order before
  the Lowdin step;
- final residual column signs are canonicalized by making the largest-magnitude
  entry of the corresponding `T_A` column positive, with lower candidate index
  breaking ties.

Approved cutoff rule:

```text
lambda_max = max(maximum(eigenvalues(Symmetric(S_R))), 0.0)
tau_keep = max(tau_abs, tau_rel * lambda_max)
tau_neg = max(tau_neg_abs, tau_neg_rel * max(lambda_max, 1.0))
```

Policy:

- eigenvalues `< -tau_neg` are construction errors;
- eigenvalues in `[-tau_neg, tau_keep]` are roundoff-negative or near-null and
  are discarded;
- eigenvalues `> tau_keep` define the residual numerical rank.

Frozen R3-A threshold values:

```text
tau_abs = 1.0e-10
tau_rel = 1.0e-10
tau_neg_abs = 1.0e-12
tau_neg_rel = 1.0e-12
```

The values are part of the persistent residual object and later artifact
provenance; ordinary-QW values remain donor evidence, not automatic authority.

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
| `retained_candidate_indices::Vector{Int}` | selected raw candidate columns in ascending candidate order |
| `residual_labels::Vector{String}` | deterministic residual labels derived from retained candidate labels and owner counts |
| `T_G::Matrix{Float64}` | shape `n_G x n_R`; base-side residual transform |
| `T_A::Matrix{Float64}` | shape `n_A x n_R`; supplement-side residual transform |
| `residual_metric_eigenvalues::Vector{Float64}` | diagnostic eigenvalues of `S_R`, not final orientation |
| `tau_abs::Float64` | absolute keep threshold |
| `tau_rel::Float64` | relative keep threshold |
| `tau_neg_abs::Float64` | absolute negative-eigenvalue error threshold |
| `tau_neg_rel::Float64` | relative negative-eigenvalue error threshold |
| `rank_rule::Symbol` | deterministic rank selector, e.g. `:pivoted_cholesky_candidate_order` |
| `orientation::Symbol` | `:selected_candidate_order_symmetric_lowdin` |
| `sign_rule::Symbol` | `:largest_T_A_entry_positive` |

Invariants:

- `size(T_G) == (base_dimension, residual_dimension)`;
- `size(T_A) == (candidate_count, residual_dimension)`;
- `length(candidate_labels) == candidate_count`;
- `length(candidate_owner_indices) == candidate_count`;
- `length(candidate_centers) == candidate_count`;
- `length(retained_candidate_indices) == residual_dimension`;
- `G' S R` is zero within tolerance;
- `R' S R` is identity within tolerance.

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
It remains candidate-only and is not approved by R3-A.

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
`G-M` and `M-M` pair blocks, matching the ordinary-QW donor convention. The
final augmented `V` assembly inserts those density-normalized blocks directly;
no second division by base final weights or effective-Gaussian weights occurs.

`V_GG_base` is reused from the accepted base localized IDA path. `V_GM` and
`V_MM` are MWG approximations, not exact residual-GTO Coulomb integrals.

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
- `V_aug` is symmetric;
- base `V_GG` block is unchanged;
- returned object is the existing `CartesianIDAHamiltonian{Float64}`;
- no Hamiltonian wrapper/result/status payload is introduced.

## R3-C Artifact And Provenance

`HP-R3-ART-01` is deferred until after R3-A and R3-B numerical endpoints are
accepted. The existing `CartesianIDAHamiltonian{Float64}` is sufficient as the
in-memory numeric Hamiltonian object when the augmented matrices satisfy the
existing contract.
It remains candidate-only and is not approved by R3-A.

The first supplemented artifact schema should be compact. It should derive
provenance from the validated R3 construction specification, not recover
`base_route` or input data from an in-memory Hamiltonian.

Candidate compact provenance keys:

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
| `supplement_provenance/residual_basis_convention` | `:selected_candidate_order_symmetric_lowdin` |
| `supplement_provenance/rank_rule` | deterministic rank selector |
| `supplement_provenance/tau_abs` | `Float64` |
| `supplement_provenance/tau_rel` | `Float64` |
| `supplement_provenance/tau_neg_abs` | `Float64` |
| `supplement_provenance/tau_neg_rel` | `Float64` |
| `supplement_provenance/mwg_convention_version` | `1` |
| `supplement_provenance/mwg_convention` | `:separable_moment_matched_density_normalized` |
| `supplement_provenance/one_body_source` | `:exact_transformed_raw_blocks` |
| `supplement_provenance/interaction_source` | `:density_normalized_mwg_blocks` |

The first compact schema should not store:

- the full residual eigenvalue vector;
- MWG center matrices;
- MWG width matrices;
- `T_G` or `T_A`;
- candidate labels.

If those arrays become consumer-critical, a later explicit residual-basis
artifact group must promote them together: residual transforms, candidate
labels/order, owner indices, retained candidate indices, eigenvalue
diagnostics, centers, and widths. Calling only centers and widths "provenance"
would be misleading.

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

R3-A first endpoint should be H2 augmented one-body and moments, not a
residual-basis-only scaffold commit. R3-B should then add MWG interaction and
the existing in-memory Hamiltonian. Be2 follows as the first
performance/realism proxy before Cr2.

Cr2 remains a later stress/consumer-readiness milestone, not the first R3
correctness gate.

## R3-A Approval Evidence

Manager log Pass 048 records the measurement-only R3-A residual-spectrum spike
used for this approval.

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
- `inv(sqrt(Symmetric(S_R)))` in candidate order succeeded;
- measured `R' S R` identity error approximately `4.4e-15`.

The spectrum is far from the frozen threshold band, so symmetric Lowdin
residualization in candidate order is not marginal for the first H2 fixture.

## Deletion And Retirement Targets

Current live-code audit found no active `src`, `test`, `tools`, or `bin`
definitions for the old supplement-preflight/provider blocker family named in
earlier plans. Do not preserve those names as R3 authority.

After R3-A:

- delete any reintroduced supplement-preflight fixture or wrapper whose only
  purpose is to announce missing overlap, one-body, or moment provider blocks;
- keep named-basis loading, Gaussian representation, and CPB provider kernels;
- keep `test/docs/cartesian_ham_builder_policy_runtests.jl` unless R3 later
  amends canonical-driver policy.

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

## Approved First Implementation Target

The approved first implementation target is R3-A only:

```text
base z-axis H2
-> contracted two-center H/cc-pVTZ lmax-1 candidates
-> deterministic residual-GTO block
-> exact augmented K, U_A, x/y/z, and x2/y2/z2
-> H2 augmented one-body endpoint
```

The first committed/standalone endpoint gate is approved by `HP-R3-TEST-01`
for the H2 augmented one-body/moment endpoint only. It must check `G' S R`,
`R' S R`, base G-G block equality, finite/symmetric augmented `K`, uncharged
`U_A`, and moment matrices, and `E1_aug <= E1_base + epsilon`. It must not
assert private pair/assembly/report/status behavior and must not run Be2 or
Cr2.

Do not approve R3-B or R3-C until R3-A produces the exact one-body/moment
endpoint. Do not start Cr2 stress tests before the H2 and Be2 augmented paths
are numerically and organizationally accepted.

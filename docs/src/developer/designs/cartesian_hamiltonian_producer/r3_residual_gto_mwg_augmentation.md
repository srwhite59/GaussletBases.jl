# R3 Residual-GTO/MWG Augmentation Candidate

Status: candidate design amendment, not implementation authority.

This document sketches the R3 path for generic residual-GTO/MWG augmentation of
the Cartesian base Hamiltonian producer. It does not approve source work, tests,
tools, driver edits, public API expansion, artifact changes, or implementation
of the candidate IDs below.

## Candidate IDs

- `HP-R3-OBJ-01` - augmented basis/residual block construction object.
- `HP-R3-FN-01` - residual-GTO/MWG augmentation boundary.
- `HP-R3-FN-02` - augmented one-body operator assembly.
- `HP-R3-FN-03` - augmented IDA/MWG interaction assembly.
- `HP-R3-ART-01` - supplemented Hamiltonian artifact provenance schema.
- `HP-R3-TEST-01` - first supplemented endpoint/proxy validation.

All R3 IDs are candidate-only until explicitly approved in `registry.md`.

## Goal

R3 restores supplement augmentation as a generic final-basis augmentation, not
as an H2-specific residual-GTO sidecar and not as a route-stage
preflight/status graph.

The common construction boundary should be:

```text
base terminal final basis G
+ Gaussian supplement candidates A
-> residual block R = orthogonal complement of A against fixed G
-> augmented basis B = [G, R]
-> augmented K / U_A / V
-> CartesianIDAHamiltonian
-> supplemented artifact / consumer handoff
```

The base `CartesianIDAHamiltonian` may provide reusable `G-G` matrix blocks and
physical center data, but augmentation is not a post-hoc patch around an opaque
completed base Hamiltonian. The supplemented Hamiltonian is assembled in the
augmented basis. Mixed `G-R` and residual `R-R` one-body blocks, residual MWG
descriptors, and residual-containing IDA terms must be constructed before the
final augmented Hamiltonian is handed to consumers.

## Donor Inventory

| Surface | Classification | R3 use |
| --- | --- | --- |
| `CartesianGaussianShellSupplementRepresentation3D` in `src/cartesian_basis_representation.jl` | active reusable representation | Use as the explicit Gaussian candidate sector `A`. It carries orbitals plus light source metadata. |
| `legacy_atomic_gaussian_supplement`, `legacy_bond_aligned_diatomic_gaussian_supplement`, `legacy_bond_aligned_heteronuclear_gaussian_supplement` in `src/legacy_basis_adapter.jl` | active donor with migration naming | Use initially for named-basis loading and deterministic Cartesian shell generation. Do not freeze `legacy_*` names into new R3 public API. |
| `cpb_mixed_gto_overlap_block`, `cpb_mixed_gto_position_operator_block`, `cpb_mixed_gto_x2_operator_block`, `cpb_mixed_gto_kinetic_operator_block` in `src/CartesianCPBBlockProviders.jl` | active reusable kernels | Use for local mixed gausslet/GTO one-body provider blocks once route support tiling/row ownership is explicit. |
| `cpb_gto_overlap_operator_block`, `cpb_gto_position_operator_block`, `cpb_gto_x2_operator_block`, `cpb_gto_kinetic_operator_block` in `src/CartesianCPBBlockProviders.jl` | active reusable kernels | Use for supplement self one-body blocks. |
| `cpb_mixed_gto_nuclear_by_center_block`, `cpb_gto_nuclear_by_center_block`, `cpb_gto_supplement_local_operator_bundle` in `src/CartesianCPBBlockProviders.jl` | active reusable kernels | Preferred first provider seam for local mixed/self overlap, moments, kinetic, and by-center nuclear blocks. |
| `_cartesian_supplement_cross_overlap`, `_cartesian_factorized_basis_supplement_cross`, `gto_overlap_matrix` in `src/cartesian_cross_overlap.jl` and `src/cartesian_gto_probes.jl` | active probes / reusable overlap kernels | Use for validation and possibly for small residual-space construction. Avoid making probe-only block-index and occupancy APIs the R3 production boundary. |
| `fit_radial_ylm_to_solid_harmonic_gto`, `radial_ylm_fit_cartesian_gto_adapter`, `project_cartesian_gto_to_supplement_subspace` in `src/radial_ylm_gto_bridge.jl` | active conversion/projection helpers | Useful for supplement candidate generation and small subspace projection diagnostics; not the Hamiltonian assembly boundary. |
| `src/ordinary_qw_residuals.jl` | oracle/reference plus possible algorithm donor | Reuse the residual-space math, keep policy, and stabilization ideas. Do not carry ordinary-QW route objects into Cartesian/PQS production. |
| `src/ordinary_qw_raw_blocks.jl` | oracle/reference | Reference for raw GTO/PQS block organization and supplement plumbing. Prefer current CPB provider kernels for R3 implementation. |
| `src/ordinary_qw_operator_assembly.jl` | oracle/reference | Reference for MWG residual moment matching and residual-containing interaction semantics. Do not copy QW payload/status structure. |
| `src/gaussian_coulomb_reference.jl` | oracle/reference only | Use for tiny dense Gaussian Coulomb checks, not production scaling. |
| `gto_overlap_matrix`, `gto_occupancy_matrix` public probes in `src/cartesian_gto_probes.jl` | diagnostics / oracle only | Keep for user diagnostics and validation; do not route production through broad probe APIs. |
| `_pqs_source_box_route_driver_diatomic_physical_gausslet_supplement_request_payload`, `_representation_payload`, `_preflight_payload` | migration-only wrappers | Historical staging surfaces. They can guide input facts but should be retired or renamed away from physical-gausslet/fake-PQS vocabulary once R3 has a consumer-owned boundary. |
| `independent_h2_pqs_supplement_preflight` artifact/status groups | migration-only evidence | Useful as an audit trail that H2/H/cc-pVTZ lmax-1 had 18 candidate orbitals and was blocked at provider blocks. Do not promote preflight/status graph shape. |
| Fake-PQS/WL supplement reproduction paths and scalar comparisons | deletion candidates | Must not be used as independent-PQS R3 evidence. Delete or quarantine after the generic endpoint replaces them. |

## Construction Choice

R3 should extend the localized final basis before final supplemented K/U/V
assembly:

1. Keep the base terminal final basis `G` fixed. Do not rotate or re-Lowdin the
   base PQS basis.
2. Build Gaussian supplement candidates `A` from an explicit supplement
   representation.
3. Build cross and self overlaps:

   ```text
   X = G' S A
   S_AA = A' S A
   S_R = S_AA - X'X
   ```

   when `G` is orthonormal. Use the general metric form only if a later design
   proves it is required.
4. Diagonalize `S_R`, retain positive residual eigenvalues above cutoff, and
   construct residual block `R = (A - G X) L_R`.
5. Assemble the augmented one-body and interaction matrices in basis `[G, R]`.
6. Return the existing `CartesianIDAHamiltonian{Float64}` with augmented
   matrices.

Adding supplement blocks around an already-complete base Hamiltonian is not the
primary construction. It is allowed only as an optimization when it is exactly
equivalent to assembling the augmented matrices in `[G, R]`, with the base
`G-G` blocks reused unchanged and all mixed/residual blocks explicitly
constructed.

## Basis And Gauge Differences From Base PQS

Residual-GTO/MWG is not base PQS shell realization:

- Base PQS terminal blocks use positive final-integral sign canonicalization and
  require positive final IDA weights.
- Residual Gaussian candidates are projected out of the fixed base final basis.
  The fixed base basis is not changed.
- Residual overlap eigenvalues near the cutoff are expected; discarded
  near-null directions are normal.
- Retained residual metric eigenvalues must be finite and positive above the
  selected cutoff, but residual coefficients, one-body matrix elements, and
  moment-derived intermediate quantities must not be forced through the base
  PQS positive-weight gauge.
- Near-zero or negative residual integrals and sign-changing contributions are
  allowed where the residual construction or physical operator permits them.
  Invalid MWG widths, nonfinite moments, or negative retained metric
  eigenvalues are errors, not gauge-fix opportunities.
- No metadata-carried matrices, report payload clouds, readiness/status graphs,
  or route-stage field expansion is part of the R3 production boundary.

## Candidate Boundary

Candidate conceptual boundary:

```julia
augment_cartesian_ida_hamiltonian_with_residual_gto_mwg(
    terminal_basis_realization::CartesianTerminalBasisRealization,
    base_hamiltonian::CartesianIDAHamiltonian{Float64},
    supplement::CartesianGaussianShellSupplementRepresentation3D,
    provider_blocks,
    residual_options,
)::CartesianIDAHamiltonian{Float64}
```

This is conceptual only. The final approved function name, owner file, and
argument types must be frozen by a later approval pass.

The boundary must own:

- supplement candidate validation;
- cross/self overlap construction against the fixed base final basis;
- residual metric diagonalization and cutoff;
- residual transform and diagnostics;
- exact one-body mixed/self residual blocks;
- MWG residual descriptors for residual-containing IDA terms;
- final augmented `K`, by-center uncharged `U_A`, and `V`;
- existing `CartesianIDAHamiltonian` construction.

The boundary must not own:

- solver work;
- ECP or EGOI corrections;
- public driver route-stage vocabulary;
- status-only preflight payloads;
- route-global matrices stored in metadata;
- Cr2 stress/performance gates.

## Artifact And Provenance

The existing `CartesianIDAHamiltonian{Float64}` is sufficient as the in-memory
numeric Hamiltonian object if the augmented matrices satisfy the existing
contract: dense finite/symmetric `K`, by-center uncharged `U_A`, localized IDA
`V`, charges, positions, and electron counts in one orthonormal working basis.

A supplemented artifact needs additional provenance. Candidate `HP-R3-ART-01`
extends the final Hamiltonian file with fixed keys under:

```text
supplement_provenance/
```

Candidate schema:

| Key | Value |
| --- | --- |
| `supplement_provenance/provenance_version` | `1` |
| `supplement_provenance/producer` | `:cartesian_residual_gto_mwg_augmentation` |
| `supplement_provenance/base_producer` | `:cartesian_base_hamiltonian` |
| `supplement_provenance/base_route` | copied from base `producer_provenance/route` |
| `supplement_provenance/supplement_policy` | `:mwg_residual_gto` |
| `supplement_provenance/supplement_kind` | `:cartesian_gaussian_shell_supplement` |
| `supplement_provenance/supplement_basis` | basis label, e.g. `"H/cc-pVTZ"` or `"Be/..."` |
| `supplement_provenance/supplement_lmax` | `Int` |
| `supplement_provenance/supplement_orbital_count` | `Int` |
| `supplement_provenance/base_final_dimension` | `Int` |
| `supplement_provenance/residual_dimension` | retained residual rank |
| `supplement_provenance/augmented_dimension` | final Hamiltonian dimension |
| `supplement_provenance/residual_metric_cutoff` | `Float64` |
| `supplement_provenance/residual_metric_eigenvalues` | `Vector{Float64}` of retained values |
| `supplement_provenance/residual_cross_overlap_max` | `Float64` |
| `supplement_provenance/residual_self_overlap_error` | `Float64` |
| `supplement_provenance/mwg_kind` | `:matched_width_gaussian` |
| `supplement_provenance/mwg_centers` | `residual_dimension x 3` `Matrix{Float64}` |
| `supplement_provenance/mwg_widths` | `residual_dimension x 3` `Matrix{Float64}` |
| `supplement_provenance/one_body_source` | `:exact_mixed_and_residual_blocks` |
| `supplement_provenance/interaction_source` | `:base_ida_plus_residual_mwg_ida` |

This candidate does not approve a new Hamiltonian wrapper, separate manifest,
public provenance reader, HamV6 export, solver-ready export, or consumer API.

## Validation Proxy

Recommended first implementation proxy: z-axis H2 with the existing public/base
H2 geometry and a small H-centered Cartesian GTO supplement such as H/cc-pVTZ
with `lmax = 1`.

Reasons:

- prior supplement audits already established a fake-free independent H2 PQS
  preflight path and an 18-orbital H/cc-pVTZ lmax-1 supplement representation;
- H2 is small enough to debug residualization, mixed/self blocks, and artifact
  provenance without Cr2-scale noise;
- it exercises a molecular supplement placement and mixed gausslet/GTO blocks.

First validation should check:

- base H2 endpoint still matches the approved base facts before augmentation;
- Gaussian candidate count and labels are deterministic;
- residual metric rank and eigenvalues are finite and cutoff-controlled;
- `G' S R` is small and `R' S R` is near identity;
- augmented `K`, every `U_A`, and `V` are finite and symmetric;
- the base `G-G` block matches the base Hamiltonian within tolerance;
- `CartesianIDAHamiltonian{Float64}` is returned with augmented dimension;
- artifact readback of the Hamiltonian still works;
- `supplement_provenance/` keys match the candidate schema;
- no private pair/assembly/report/status field is asserted as public behavior.

Recommended second proxy: Be2 with the corrected base PQS path, after H2 proves
the boundary. Be2 is the first useful performance/realism proxy for Cr2-facing
work, but it should not be the first place to debug the residual construction
contract.

Cr2 remains a later stress/consumer-readiness milestone, not the first R3
correctness gate.

## Deletion And Retirement Targets

After the first approved R3 implementation slice lands and has a real
supplemented endpoint:

- delete or quarantine fake-PQS/WL supplement reproduction artifacts and tests
  that imply independent-PQS supplement readiness;
- retire `physical_gausslet`-named supplement request/representation/preflight
  wrappers from the production path, or keep them only as historical migration
  compatibility with no public authority;
- remove status-only supplement preflight/report payload clouds once the
  consumer-owned augmentation boundary exists;
- delete blocked artifact fields whose only purpose was to report
  `:missing_provider_gto_supplement_blocks`,
  `:missing_residual_mwg_representation`, or
  `:missing_combined_density_density_readiness`;
- keep CPB provider kernels and their focused tests while they are active
  construction kernels;
- keep ordinary-QW residual and MWG paths as oracle/reference until R3 parity is
  reviewed, then classify any duplicated donor code for deletion.

## Recommendation

First implementation should be H2, not Be2 or Cr2:

```text
base z-axis H2
-> H/cc-pVTZ lmax-1 supplement candidates
-> residual-GTO block
-> MWG residual descriptors
-> augmented CartesianIDAHamiltonian
-> supplemented artifact provenance
```

Then run Be2 as the first performance/realism proxy before Cr2.

This R3 candidate deliberately does not approve implementation. A follow-up
review should decide whether the first approved slice starts with provider
blocks only, residualization only, or a complete tiny H2 augmented endpoint.

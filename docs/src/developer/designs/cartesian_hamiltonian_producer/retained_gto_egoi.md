# Retained-GTO Local-Product EGOI

Lifecycle:

- `HP-RG-PROTECT-EGOI-AUDIT-01` is a completed historical measurement, not
  active implementation authority;
- `HP-RG-PROTECT-EGOI-FN-01` is approved pending internal source authority and
  is not implemented in committed source;
- `HP-RG-PROTECT-EGOI-TEST-01` is the approved pending validation contract.

This page is the canonical contract for the protected-localized
retained-original-GTO EGOI lane. The registry owns lifecycle, permission, and
file surfaces. This page owns target selection, product and mask semantics,
diagnostics, failure behavior, and exclusions.

The generic matrix-level EGOI routines are already implemented in
`src/hamiltonian_corrections.jl`. The protected retained-GTO helper described
here is not. Uncommitted worktree additions in that file do not change this
lifecycle and are not accepted by this documentation contract.

## Purpose

The pending helper would construct a small, local density-density correction
for interaction errors seen by compact retained original supplement GTOs in a
protected-localized basis. It is not a general correction framework and does
not choose the protected basis itself.

The approved first physical target is:

```text
retained original supplement GTOs
mapped from compact retained source indices
owner-balanced s1+s2 only
represented in native protected-localized L coordinates
local symmetric products as first-class products
M2-local interaction variables
```

The correction remains an in-memory `DeltaV`. It does not define a corrected
artifact or solver workflow.

## Dependencies

The lane depends on:

- the [protected-localized basis convention](protected_localized_basis.md) for
  compact retained-source provenance, native `L`, and inherited-site `Vee_L`;
- the [protected-localized artifact contract](protected_localized_artifact.md)
  only when a measurement consumes a persisted protected member;
- committed matrix-level routines:
  - `egoi_target_product_matrix`;
  - `egoi_target_coulomb_matrix`;
  - `egoi_density_density_correction`;
  - `egoi_stationary_hamiltonian_correction`.

Those routines remain the generic EGOI algebra. The protected helper may
compose them or add only the narrow symmetric/local-product constrained layer
needed by this contract. It must not duplicate or fork the matrix-level law.

## Retained Original Targets

Target orbitals come from original supplement GTO columns whose source indices
survive compact residual selection. They are not inferred from final rows,
localized sectors, broad protected directions, or orbital labels detached from
source provenance.

For the initial homonuclear convention, select the retained original `s1+s2`
pair on each owner. Selection must be owner-balanced: the same approved target
class must exist on both atoms. Missing, duplicated, or ambiguous source-index
mapping is a hard failure.

The following are not targets under this authority:

- broad protected-`Z` directions;
- residualized RG functions;
- final protected-localized basis rows;
- atom-HF or packet occupied orbitals;
- `s3`, `p`, or `d` retained channels.

## Native-L Representation

Let `Qtarget` contain the retained original GTO targets represented in the
actual protected-localized basis:

```text
size(Qtarget) = (size(Vee_L, 1), target_count)
```

Rows follow native `L` order. Do not z-sort `Qtarget`, reinterpret artifact
sector ranges as localized target labels, or rebuild targets from row-center
metadata. The representation must use the same protected geometry and
transform that produced `H1_L` and `Vee_L`.

Report each target's owner, source index, channel, norm, Gram contribution,
and projection loss. Material projection loss, rank loss caused by incorrect
source mapping, or disagreement with protected-source provenance blocks the
helper. Do not hide it with normalization or a broader target class.

## First-Class Products And Acceptance Blocks

For each atom, form symmetric products of its local retained targets. With
owner `A`, the first-class product set is:

```text
AA = { q_a .* q_b : a <= b, owner(a) = owner(b) = A }
```

and analogously for owner `B`. The exact-target and residual acceptance metric
must include:

- `AA-AA`: Coulomb among owner-A local products;
- `BB-BB`: Coulomb among owner-B local products;
- `AA-BB` and its transpose: inter-atom Coulomb between the two local-product
  sets.

`AA-BB` is not an AB overlap-product target. Products of one A orbital and one
B orbital, `q_a .* q_b` with different owners, are excluded as first-class
products. The completed audit found no justification for promoting those AB
overlap products.

The acceptance diagnostics must keep the three local-product blocks separate.
For `AA-BB`, also report diagonal-diagonal, diagonal-offdiagonal, and
offdiagonal-offdiagonal contributions.

## M2 Interaction Mask

The only approved initial interaction mask is `M2`. In native `L` order,
define the row scale

```text
ell_i = max(nearest-neighbor scale_i, core_spacing)
```

and permit an interaction correction only when

```text
r_ij <= 1.75 * max(ell_i, ell_j)
```

where `r_ij` is the distance between the two native-row centers. The mask must
be symmetric and derived from the actual protected-localized row geometry.

Every long-range or otherwise disallowed entry must remain exactly zero:

```text
max_disallowed_delta_v == 0
```

Do not weaken this to a tolerance, post hoc truncation, or reporting-only
check. A nonzero disallowed entry is a construction failure.

## Pending Helper Result

The approved source lane may return an internal in-memory result containing:

- corrected interaction `Vee_L + DeltaV`;
- symmetric `DeltaV` in native order;
- compact target, product, mask, residual, and correction diagnostics.

It may not persist a new result shape or artifact field under this authority.
The correction must not mutate the inherited-site meaning of the input
`Vee_L` or be described as direct `C' * V * C` interaction rotation.

## Required Diagnostics

The helper and validation must report:

- target labels, source indices, owners, and channels;
- target norms, Gram matrix, projection loss, and represented rank;
- first-class product count, symmetric product rank, singular values, and
  numerical rank threshold;
- residual Frobenius and maximum norms before and after correction for
  `AA-AA`, `BB-BB`, and `AA-BB`;
- the `AA-BB` diagonal class split;
- `DeltaV/V` Frobenius norm, maximum, `p95`, and median by allowed-variable
  class;
- saturated-variable counts by class;
- `max_disallowed_delta_v`, which must be exactly zero;
- corrected interaction finiteness and symmetry;
- low-Fock spectral shift diagnostics;
- parity with the accepted measurement convention when replay data are used.

Regularization or rank deficiency must remain visible. A fit that relies on
large saturated corrections, loses the accepted residual behavior, or creates
a bad low-Fock mode does not pass validation merely because a matrix was
returned.

## Validation Contract

`HP-RG-PROTECT-EGOI-TEST-01` is pending with the source helper. No committed
protected retained-GTO helper test currently exists.

The approved validation shape is:

- package load and `git diff --check`;
- bounded H, Be, and Be2 retained-original-GTO smokes for `s1` and `s1+s2`;
- projection, rank, residual, correction-size, saturation, symmetry, and
  low-Fock checks from this contract;
- an ignored Cr2 replay of the accepted `s1+s2`/`M2` measurement;
- exact `max_disallowed_delta_v = 0`;
- no production Cr2 HF and no committed large Cr2 fixture.

The exact focused test file should be assigned when implementation resumes;
this extraction does not create one or accept the current worktree WIP.

## Failure Behavior

Stop rather than broadening policy when any of these occurs:

- retained source indices cannot be mapped uniquely or owner balance fails;
- `Qtarget` is materially underrepresented or inconsistent with native `L`;
- local-product or exact-target dimensions disagree;
- mask geometry is missing, nonfinite, or nonsymmetric;
- any disallowed `DeltaV` entry is nonzero;
- corrected `Vee` is nonfinite or nonsymmetric;
- required block, rank, saturation, or low-Fock diagnostics are unavailable.

Poor residual reduction, large corrections, widespread saturation, or bad
low-Fock behavior fail the validation gate. They do not authorize adding AB
overlap products, promoting channels, weakening `M2`, or writing an artifact.

## Ownership

Approved pending primary source surface:

- `src/hamiltonian_corrections.jl`.

Optional narrow mapping support, only if required to obtain retained source
provenance or transform-ready `Qtarget`:

- `src/cartesian_residual_gaussians/residual_basis.jl`;
- `src/cartesian_residual_gaussians/augmented_operators.jl`.

No public export, new source file, artifact writer, or workflow caller is
approved.

## Explicit Non-Goals

This contract does not approve:

- corrected protected-localized artifact variants or provenance fields;
- ladder transfer, bundles, or restart sidecars;
- occupied-first or additive-reference construction;
- rho0 or screened-Hartree changes;
- solver, HF, MP2-NO, or public-driver integration;
- Cr2 production energy claims;
- RG/injection selection-policy changes;
- broad protected-`Z`, atom-HF, packet, final-row, or residualized-RG targets;
- AB overlap products as default first-class targets;
- `s3`, `p`, or `d` target promotion.

## Historical Evidence

The completed measurement sequence and numerical values remain in manager
running-log Passes 302-308. Pass 309 records a transient uncommitted source
WIP; it is not evidence that the helper reached committed implementation.
Probe paths and machine-local output directories remain historical evidence,
not this canonical contract.

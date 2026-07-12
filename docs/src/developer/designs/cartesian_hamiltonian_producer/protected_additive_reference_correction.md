# Protected Additive Atomic Reference Correction

Status: implemented narrow internal, opt-in facility under
`HP-RG-PROTECT-ADDREF-FN-01` and
`HP-RG-PROTECT-ADDREF-TEST-01`. The implementation landed in commit
`0b778a676`.

This protected replacement path remains implemented. The separate
[numerical-complete residual basis](numerical_complete_residual_basis.md) is an
implemented alternative that preserves `G` and appends `R_num`; it does
not alter this contract or reuse protected artifact semantics.

This lane builds one protected-localized molecular member whose basis contains
the complete span of all placed atomic packet occupied spaces, then constructs
the existing screened direct-Hartree correction in that member's native `L`
order. It is not public workflow, artifact, solver, or endpoint authority.

## Physics Target

For a two-center molecule such as Cr2, place one converged atomic reference
packet on each nucleus. The protected basis must contain both packet occupied
spaces. The reference density remains the additive atomic density:

```text
P0 = P0_A + P0_B
q0 = diag(P0)
J0 = J_A + J_B
E0 = E_AA + E_BB + 2 E_AB
```

`E0` is the no-half Coulomb quantity used by the existing correction helper.
For more than two packets, the general rule is:

```text
E0 = sum_a E_aa + 2 * sum_{a<b} E_ab
```

The orthonormalized union used to protect the basis is not the reference
density. Each original packet occupied block and its occupations define its
own additive contribution to `P0`.

## Current Live Seams

The implementation in `src/cartesian_protected_ladder_bundle.jl` builds
`R_compact`, embeds each packet occupied block, and passes those blocks with
the exact residual object to
`staged_protected_original_injection_geometry(...)`. The separate
`occupied_first_injection_geometry(...)`, documented in
[occupied-first injection geometry](occupied_first_injection.md), is tested but
is not called by this builder. It is not a direct replacement because the
protected path operates over `M = [G, R_compact]` under the existing
[protected-localized basis convention](protected_localized_basis.md).

Existing reusable owners already provide:

- raw exact reference Hartree `GG`, `GA`, and `AA` blocks in
  `src/cartesian_gaussian_raw_blocks/mixed_hartree_blocks.jl`;
- protected fixed-sector and localized one-body transformation from the
  protected-localized basis owner;
- the existing `ScreenedHartreeCorrection` algebra in
  `src/cartesian_reference_density/screened_hartree_correction.jl`.

The implemented path combines staged protected-original geometry with
mandatory occupied blocks, fast placed fitted-potential `GG/GA/AA`, additive
packet reference assembly, and the final `F -> L` handoff while reusing the
owners above.

This is an explicit internal opt-in path. A protected member built without
placed reference packets must preserve the current geometry, `H1_L`, `Vee_L`,
ordering, and artifact behavior.

The implemented composition seam is a private sibling of the current member
builder:

```text
_plb_build_additive_reference_member(recipe, stages, placements)
    -> (member, correction, reference)
```

Each normalized placement is a small stable record containing `packet`,
`owner_index`, `center`, and explicit `supplement_indices`. The sibling shares
the existing member-build core rather than duplicating member construction.
`_plb_build_member(recipe, stages)` remains the unchanged no-reference path.

## One Compact Residual Construction

The protected path builds `R_compact` once. The staged protected geometry
consumes that exact `CartesianResidualGaussianBasis`; it does not rerun ordered
compact-first selection to reconstruct a nominally equivalent `M`.

The residual object carries the native source fact required to protect
originals corresponding to compact residuals in one vector-backed internal
field:

```text
compact_source_candidate_indices::Union{Nothing,Vector{Int}}
```

For ordered compact-first MGS, it records the sorted unique set of accepted
source candidates whose originals are protected. It is not a one-to-one label
for final residual columns, which may rotate during final merge cleanup. For
selection rules without native accepted-source semantics, it is `nothing`.
Do not parse residual labels to recover this fact. No artifact field or public
result is approved.

The old duplicate compact-selection block inside staged protected geometry is
deleted; the live caller passes the already-built residual object.

## Mandatory Occupied Span

For each placed packet `a`, embed its original occupied coefficients into the
combined molecular supplement coordinates:

```text
Y_a in A coordinates
n_a = packet occupations
```

Packet self-integrity must first be established by exactly recomputing the
stored packet overlap fingerprint. Mapping must then use exact packet atom and
basis identity, function count, owner-local supplement indices, placement
center, ordered labels, angular powers, and packet-to-molecular column order.
Label-only picks are forbidden.

The translated or reconstructed owner-local molecular overlap block is an
equivalent numerical representation, not packet storage bytes. After the exact
structural checks, require

```text
norm(S_AA[indices, indices] - packet.overlap, Inf) <= overlap_atol
```

with the unchanged default `overlap_atol = 1e-10`. A different mapped-block
SHA-256 fingerprint is diagnostic and is not by itself a failure. Return one
nested internal overlap-mapping summary with the stored, recomputed, and mapped
fingerprints; mapped exact-match boolean; maximum absolute and infinity-norm
differences; and tolerance. Do not flatten these diagnostics into staged or
artifact fields.

Validate each packet block separately:

```text
Y_a' S_AA Y_a = I
sum(n_a) = packet electron count
packet RHF converged
packet self-fingerprint is exact
packet identity/order/owner/placement maps exactly
owner-local overlap agrees numerically within overlap_atol
```

Form a full-rank `S_AA`-orthonormal basis for the union:

```text
Y_all   = [Y_1, Y_2, ...]
Y_union = orth_SAA(Y_all)
```

The union Gram spectrum, retained rank, and each packet block's recovery from
`Y_union` are required diagnostics. A true Gram-null direction may be removed
from the union basis, but every original packet occupied block must still be
recovered at roundoff.

Build mandatory protected content in this order:

1. `Y_union`;
2. current compact-original protected directions, orthogonalized against
   `Y_union` and rank-cleaned;
3. optional supplement directions only after the complete mandatory block is
   fixed.

The current staged representability (`s_cut`) and fake-RDM (`occ_cut`) policy
for optional broad directions remains unchanged and applies only in the
`S_AA` complement of the mandatory block. Mandatory occupied directions are
never discarded by a cutoff. If their projection into `M` is rank deficient
or fails the active staged representability threshold, stop and report
insufficient compact-main-basis support.

The final geometry remains replacement, not append:

```text
Z = [Z_mandatory, Z_optional]
F = [Z, M Q_perp]
L = F W
```

Required recovery checks cover every original packet block in `Z`, `F`, and
`L`, not only the orthonormalized union.

## Additive Reference Density In L

Keep the original embedded `Y_a` blocks separate from `Y_union`. Given native
localized raw coefficient maps `G_L` and `A_L`, represent each packet occupied
block by the final-basis cross overlap:

```text
C_aL = G_L' * X_GA * Y_a + A_L' * S_AA * Y_a
```

The final basis is orthonormal. Do not introduce generalized final-basis
overlap logic.

Validate each block independently:

```text
C_aL' C_aL = I
Tr(C_aL * Diagonal(n_a) * C_aL') = sum(n_a)
```

Then assemble:

```text
P0_L = sum_a C_aL * Diagonal(n_a) * C_aL'
q0_L = diag(P0_L)
```

Do not globally orthogonalize `[C_1L, C_2L, ...]` when constructing `P0_L`.
Inter-packet occupied overlap is a reported physical/reference diagnostic, not
a reason to rotate the additive atomic densities.

## Placed Fitted-Potential Hartree Field

The practical packet path uses each packet's fitted potential to evaluate the
same fitted density cloud field efficiently. The neutral raw-block owner in
`src/cartesian_gaussian_raw_blocks/mixed_hartree_blocks.jl` applies an explicit
placed spherical Gaussian potential term list through the existing `GG`, `GA`,
and `AA` axis factors and block assembly. The packet and ladder owners do not
duplicate those analytic Gaussian loops.

For every packet placement, build raw blocks:

```text
J_a_raw = (GG_a, GA_a, AA_a)
```

Sum them before the protected transform:

```text
J0_raw = sum_a J_a_raw
J0_F   = transform_protected_original_fixed_sector_exact_hartree(J0_raw, geometry)
J0_L   = W' * J0_F * W
```

`J0_L` must remain in native protected-localized order. The density-fit exact
mixed-Hartree path remains the small-system oracle; it is not the practical
Cr2 construction.

Every placement must use a packet produced by the ordinary determinant ->
density-fit -> radial-potential-fit pipeline. Packets carrying retired
`potential_fit/moment_polish/*` provenance are rejected and must be
regenerated. The additive lane must not re-fit packet potentials, adjust their
coefficients for molecular separations, or add a scalar anchor patch.

## Fitted-Potential Consistency Reporting

The density fits define `E_aa` and `E_ab`; the ordinary fitted potentials
provide approximate fields `J_a`. In native `L` order report:

```text
epsilon_aa = Tr(P_a * J_a) - E_aa
epsilon_ab = Tr(P_a * J_b) + Tr(P_b * J_a) - 2 E_ab
epsilon_total = sum_a epsilon_aa + sum_{a<b} epsilon_ab
```

The assembled total must also satisfy

```text
epsilon_total = Tr(P0_L * J0_L) - E0
```

up to numerical assembly tolerance. Report the total, every available self
term, and every available pair cross term. Their algebraic decomposition is a
strict check; their physical magnitude is a fitted-potential approximation
diagnostic and is not required to be below `1e-8 Ha`. Exact/density-fit oracle
fields retain strict energy identities.

## Additive Reference Self-Energy

Packet density fits, not potential-fit Gaussians, define reference Coulomb
self and cross energies. Use the packet-local compact Coulomb role explicitly.

For each placement pair, evaluate:

```text
E_aa = (rho_a | rho_a)
E_ab = (rho_a | rho_b)
```

The cross helper belongs with atomic packet density-cloud evaluation in
`atomic_hf_reference_packets.jl`. Validate `E_ab = E_ba` and compare the total
against a combined-cloud compact oracle on the bounded Be2 gate. Do not omit
the factor of two on cross terms.

## Existing Correction Algebra

The internal additive entry point in the existing same-basis screened-Hartree
owner validates each coefficient block separately rather than requiring the
concatenated occupied columns to be globally orthonormal.

With native-order `Vee_L`, return the existing correction object:

```text
Delta_J0 = J0_L - Diagonal(Vee_L * q0_L)
C        = 0.5 * q0_L' * Vee_L * q0_L - 0.5 * E0
```

The helper reuses the existing energy and derivative anchor core. Do not copy
the correction formula into the ladder owner.

`H1_L` and `Vee_L` remain unchanged. The correction is returned separately in
memory for an explicit off/on comparison. No corrected artifact is approved.

## Implemented Ownership

- `src/cartesian_residual_gaussians/residual_basis.jl` owns compact residual
  source indices, the mandatory occupied union, and staged protected geometry.
- `src/cartesian_residual_gaussians/augmented_operators.jl` owns native-`L`
  packet representation and protected fixed/localized Hartree transforms.
- `src/cartesian_gaussian_raw_blocks/mixed_hartree_blocks.jl` owns placed
  fitted-potential `GG/GA/AA` assembly.
- `src/cartesian_reference_density/atomic_hf_reference_packets.jl` owns packet
  identity and embedding checks, fitted fields, and compact density-cloud
  self/cross energies.
- `src/cartesian_reference_density/screened_hartree_correction.jl` owns
  additive `P0_L/q0_L`, consistency diagnostics, and correction algebra.
- `src/cartesian_protected_ladder_bundle.jl` owns the private composition
  sibling and returns the member, correction, and reference diagnostics.

The registry owns exact source authority. This lane adds no public export,
artifact field, corrected artifact, or persistent workflow object. Diagnostics
remain compact module-owned records, and the ordinary no-reference member path
must remain numerically unchanged.

## Validation And Evidence

Commit `0b778a676` is the implementation evidence. Its committed checks cover:

- mandatory occupied-union rank and recovery in `test/misc/runtests.jl`;
- packet embedding and identity failures in the atomic-packet and
  screened-Hartree nested tests;
- separate packet-block traces, additive `P0_L/q0_L`, self/cross density-cloud
  energy, raw fitted-potential blocks, and correction-anchor algebra in
  `test/nested/cartesian_screened_hartree_correction_runtests.jl`.

Exact same-packet and numerically equivalent mapped overlaps pass only after
packet identity checks; overlap error above `1e-10`, corrupted packet storage,
or changed owner/order metadata fails. The accepted bounded Be2 smoke used the
now-retired determinant-moment polish, so only its structural, recovery, and
derivative evidence remains relevant. Its forced sub-`1e-8 Ha` fitted-field
result is historical false-start evidence, not validation of the ordinary fit
or an endpoint claim. Current consumption of those legacy references requires
regenerated ordinary packets.

## Failure Rules

Stop and report if:

- stored packet overlap fingerprint cannot be reproduced exactly;
- packet atom/basis, function count, owner, placement, labels, angular powers,
  or packet-to-molecular order does not map exactly;
- mapped owner-local overlap differs from packet overlap by more than
  `overlap_atol` in matrix infinity norm;
- a packet occupied block is not separately recovered at roundoff;
- the mandatory union is not stably representable by `M`;
- fitted-potential `GA/AA` requires duplicated Gaussian kernels rather than the
  neutral raw-block owner;
- cross self-energy cannot be evaluated in the documented compact density
  convention;
- a packet carries retired moment-polish provenance;
- the reported total fitted-potential consistency error does not agree with
  its assembled self/cross decomposition within numerical tolerance;
- the additive reference can be formed only by globally orthogonalizing packet
  occupied blocks;
- the correction requires changing `H1_L`, rotating `Vee_L`, or adding an
  artifact/public/solver path.

Do not turn a failed mandatory occupied direction into an RG, silently drop it,
weaken the active representability threshold, or reject an otherwise valid
ordinary fitted potential solely because its reported consistency magnitude
exceeds `1e-8 Ha`.

## Explicitly Separate Later Lanes

This authority does not approve:

- protected one-center atom compactness; the current two-owner compactness
  assumption remains unchanged;
- counterpoise artifacts or retention of separated kinetic/unit-nuclear
  matrices in protected bundle artifacts;
- compact-to-high transfer machinery;
- public driver/API/export/default changes;
- corrected protected-localized artifacts;
- solver, HF, MP2-NO, or production CR2 workflow;
- `Vee` transformation, `C' V C`, or four-index interactions;
- EGOI, exchange, rho0 row-gauge, residual-selection, injection-threshold, or
  mapping-default changes;
- endpoint or publication claims.

Before any future compact/high transfer authority, run a same-commit audit that
compares `X`, `S_AA`, compact residual geometry, protected `G_L/A_L`, native
ordering, and exact final-basis cross overlap. Coulomb accuracy changes PGDG
factor tables but does not by itself establish a basis change. Historical row
count differences are not transfer evidence.

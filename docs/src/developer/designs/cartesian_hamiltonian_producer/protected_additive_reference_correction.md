# Protected Additive Atomic Reference Correction

Status: implemented narrow internal, opt-in facility under
`HP-RG-PROTECT-ADDREF-FN-01` and
`HP-RG-PROTECT-ADDREF-TEST-01`.

This protected replacement path remains implemented. The separate
[numerical-complete residual basis](numerical_complete_residual_basis.md) is an
approved-pending alternative that preserves `G` and appends `R_num`; it does
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

The current protected builder in
`src/cartesian_protected_ladder_bundle.jl` calls
`staged_protected_original_injection_geometry(...)`. The implemented
[occupied-first injection geometry](occupied_first_injection.md) is tested but
not wired into that builder. It is not a direct replacement because the
protected path operates over `M = [G, R_compact]` under the existing
[protected-localized basis convention](protected_localized_basis.md).

Existing reusable owners already provide:

- raw exact reference Hartree `GG`, `GA`, and `AA` blocks in
  `src/cartesian_gaussian_raw_blocks/mixed_hartree_blocks.jl`;
- protected fixed-sector and localized one-body transformation from the
  protected-localized basis owner;
- the existing `ScreenedHartreeCorrection` algebra in
  `src/cartesian_reference_density/screened_hartree_correction.jl`.

The implemented path combines mandatory occupied geometry, fast placed
fitted-potential `GG/GA/AA`, additive packet reference assembly, and the final
`F -> L` handoff while reusing the owners above.

This is an explicit internal opt-in path. A protected member built without
placed reference packets must preserve the current geometry, `H1_L`, `Vee_L`,
ordering, and artifact behavior.

The implemented composition seam is a private sibling of the current member
builder:

```text
_plb_build_additive_reference_member(recipe, stages, placements)
    -> (member, correction)
```

Each normalized placement is a small stable record containing `packet`,
`owner_index`, `center`, and explicit `supplement_indices`. The sibling must
share the existing member-build core; it must not duplicate member
construction. `_plb_build_member(recipe, stages)` remains the unchanged
no-reference path.

## One Compact Residual Construction

The protected path must build `R_compact` once. The staged protected geometry
must consume that exact `CartesianResidualGaussianBasis`; it must not rerun
ordered compact-first selection to reconstruct a nominally equivalent `M`.

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
same fitted density cloud field efficiently. Add a neutral internal helper in
`src/cartesian_gaussian_raw_blocks/mixed_hartree_blocks.jl` that applies an
explicit placed spherical Gaussian potential term list to the existing `GG`,
`GA`, and `AA` machinery. It must reuse the current axis factors and mixed/self
block assembly; do not duplicate analytic Gaussian loops in the packet or
ladder owner.

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

The helper must reuse the existing energy and derivative anchor core. Do not
copy the correction formula into the ladder owner.

`H1_L` and `Vee_L` remain unchanged. The correction is returned separately in
memory for an explicit off/on comparison. No corrected artifact is approved.

## Approved Source Surface

- `src/cartesian_residual_gaussians/residual_basis.jl`
  - carry native compact residual source indices if needed;
  - consume the already-built compact residual in staged geometry;
  - build the mandatory occupied union before optional selection;
  - delete/delegate duplicate compact selection.
- `src/cartesian_residual_gaussians/augmented_operators.jl`
  - reuse the protected fixed-sector raw operator transform;
  - represent original packet occupied blocks in native `L` order;
  - apply the existing localized `W` transform to `J0_F`.
- `src/cartesian_gaussian_raw_blocks/mixed_hartree_blocks.jl`
  - add one neutral explicit fitted-potential `GG/GA/AA` entry point reusing
    existing factor/block kernels.
- `src/cartesian_reference_density/atomic_hf_reference_packets.jl`
  - validate packet self-integrity and exact owner-local mapping;
  - accept only numerically equivalent mapped overlap blocks and return one
    nested overlap-mapping summary;
  - validate/place packet occupied blocks;
  - expose placed fitted-potential raw blocks through the neutral owner;
  - evaluate compact density-cloud cross energies.
- `src/cartesian_reference_density/screened_hartree_correction.jl`
  - assemble additive `P0_L/q0_L` from separately validated blocks;
  - call the existing correction/anchor core.
- `src/cartesian_protected_ladder_bundle.jl`
  - a private additive-reference sibling that reuses the existing member-build
    core and returns the member plus existing correction object;
  - no new public recipe input and no artifact field.

No new source file, public export, artifact schema, or persistent workflow
object is approved. Related diagnostics must be nested in compact module-owned
records rather than added as a flat staged field cloud.

The embedding-equivalence follow-on itself may edit only
`src/cartesian_reference_density/atomic_hf_reference_packets.jl` and, if
diagnostic forwarding is directly required, the existing private additive-
reference caller. The amendment does not reopen the other implemented source
surfaces above.

The moment-polish retirement follow-on may use the already-approved packet and
screened-Hartree owners to remove the retired fit and energy rejection. This
additive lane may edit only its existing private caller when needed to report
the total/self/cross consistency decomposition or reject retired packets. It
must not change protected geometry, `H1_L`, `Vee_L`, placement algebra, or the
ordinary no-reference path.

The no-reference protected member path must remain numerically unchanged.

Target source budget is at most 350 added lines across the approved files,
with deletion of the duplicate compact-selection construction counted
separately. If the numerical work needs substantially more or needs a second
correction model, stop and report the missing abstraction.

## Validation

Committed correctness coverage may extend only:

- `test/misc/runtests.jl` for a tiny mandatory-union/additive-density contract;
- `test/nested/cartesian_screened_hartree_correction_runtests.jl` for additive
  block validation and anchor algebra.

Do not add a new committed test file or binary fixture in this lane.

Focused embedding coverage must prove that exact same-packet mapping passes;
a numerically equivalent translated/reconstructed block with a different hash
passes; an infinity-norm difference above `1e-10` fails; and a corrupted stored
packet fingerprint fails. Reordered labels or powers, wrong owner, wrong
center, or changed packet-to-molecular order remain hard failures. Rerun Cr2
preflight only after these focused source tests and existing packet/additive-
reference tests pass.

The first end-to-end acceptance gate is an ignored, source-backed, physically
padded Be2 construction using two converged Be core `2e` cc-pV5Z, `lmax = 1`
packets. Use driver-style padding of at least `10` bohr and inspect terminal
due diligence before interpreting the result.

The Be2 gate may use one ignored `tmp/work/*.jl` probe and durable text/TSV
outputs under `/Users/srw/dmrgtmp`. Do not commit its generated packet,
Hamiltonian, or matrix fixtures.

Required Be2 evidence:

- omitted-reference protected member parity with the current path;
- one compact residual basis is shared by geometry and operators;
- mandatory occupied union rank/Gram spectrum and roundoff recovery for each
  packet before and after localization;
- per-packet represented trace `2`, total `P0/q0` charge `4`;
- explicit `E_AA`, `E_BB`, `E_AB`, reversed-cross parity, and total
  `E0 = E_AA + E_BB + 2E_AB`;
- finite/symmetric placed raw blocks, `J0_F`, `J0_L`, and `Delta_J0`;
- strict derivative/algebra checks;
- fitted-potential total consistency error and its self/cross decomposition,
  without a `1e-8 Ha` magnitude gate;
- optional staged capture/fake-RDM counts after mandatory inclusion;
- exact confirmation that the unscreened `H1_L` and `Vee_L` inputs were not
  mutated;
- terminal parent bounds, axis counts, padding/radius, final dimension,
  retained counts, shell/slab topology, and warning flags;
- phase timings and carrying-cost report.

No endpoint energy or SCF assertion is required. The earlier padded Be2 run
used the now-retired determinant-moment polish; its structural, recovery, and
derivative evidence remains useful, but its forced sub-`1e-8 Ha` energy result
is historical false-start evidence. Regenerate ordinary Be/Ne/Cr packets and
rerun the bounded construction before further consumption. The additive
implementation itself remains accepted; this is not a production claim or
repo test.

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

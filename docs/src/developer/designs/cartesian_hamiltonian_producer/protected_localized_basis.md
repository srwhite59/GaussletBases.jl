# Protected-Localized Basis Convention

Status: implemented internal, default-off subsystem. The compact-main design
rationale is recorded under `HP-RG-PROTECT-INJECT-DESIGN-01`; staged geometry
and exact one-body source ownership are implemented under
`HP-RG-PROTECT-INJECT-FN-01` / `HP-RG-PROTECT-INJECT-TEST-01` and
`HP-RG-PROTECT-ONEBODY-FN-01` / `HP-RG-PROTECT-ONEBODY-TEST-01`. The durable
interaction decision comes from completed audit
`HP-RG-PROTECT-VEE-AUDIT-01`.

This page is the canonical numerical contract for protected-original
replacement, its localized one-particle basis, and the inherited-site
two-index interaction convention. The registry owns permission, lifecycle,
and exact implementation surfaces. Separate subsystem pages own persistence,
[row-locality metadata and artifact behavior](protected_localized_artifact.md),
[retained-GTO EGOI](retained_gto_egoi.md),
[ladder transfer/bundles](protected_localized_ladder.md), occupied-first
reference geometry, and additive screened-Hartree assembly.

## Purpose

The construction adds selected original contracted Gaussian directions to a
compact gausslet-plus-residual main space without appending a second,
independent basis sector. Original directions replace a represented subspace
of the main basis. The resulting basis is then localized back toward the main
site order so exact one-body improvement can coexist with the existing
two-index IDA/MWG interaction.

The durable sequence is:

```text
build compact residuals first
form M = [G, R_compact]
select protected and represented original Gaussians Z
replace their represented M subspace: F = [Z, M Q_perp]
localize F toward M site order: L = F W
transform exact one-body operators into L
inherit the pre-injection site-order Vee_M as Vee_L
```

## Compact-First Main Space

`G` is the orthonormal terminal gausslet/final-PQS basis. `R_compact` is the
already-built compact/narrow residual Gaussian basis selected by the ordered
compact-first policy. Define:

```text
M = [G, R_compact]
```

The geometry and exact operator paths must consume the same residual object.
Do not rebuild compact selection inside protected geometry, infer source
indices from labels, or construct a nominally equivalent second `R_compact`.

The residual object supplies its `G` and supplement coefficient blocks
`T_G/T_A` and the accepted compact source indices. Those indices identify the
original supplement Gaussians whose span is protected first.

## Protected Originals And Gaussian Cleanup

Let `A` be the original supplement and `S_AA` its overlap matrix. Original
directions corresponding to accepted compact residual sources form the initial
protected block. Orthonormalize that block in `S_AA` and discard only genuine
Gaussian metric linear dependence using a rank tolerance of the established
absolute/relative form.

Remaining original candidates are first projected out of the protected block
in the Gaussian metric and then Gram-cleaned in their own complement. This is
Gaussian cleanup only. It answers whether candidate Gaussian columns are
linearly independent; it does not answer whether the compact main space can
represent them.

The protected span must survive later cleanup. Report its self-overlap,
cross-overlap with accepted broad directions, and minimum/maximum protected
span singular values.

## Broad Candidate Representability

For a cleaned broad candidate subspace `W`, form its overlap with the compact
main space:

```text
B_W = M' S W
```

Diagonalize `B_W' B_W` and retain only the subspace whose representability
singular values meet the active `s_cut`. In the implemented staged variant,
apply the fake-RDM eigenspace filter only after this representability filter,
retaining occupations that meet `occ_cut`. Shape/localization classification
may remain diagnostic; it is not an alternative to the subspace gates.

The selected broad block is `Z_broad`. With the protected block
`Z_protected`, define:

```text
Z = [Z_protected, Z_broad]
B = M' S Z
```

Gaussian Gram cleanup and representability are different contracts:

- Gram cleanup removes dependent Gaussian directions.
- Representability asks whether an independent physical Gaussian direction is
  supported by `M`.

A good-norm broad direction that fails representability is an
insufficient-main-basis diagnostic. Report its owner/channel information and
reject it. It must not become an MWG residual channel, and the cutoff must not
be weakened merely to retain it.

## Replacement Fixed Sector

Construct an orthonormal complement `Q_perp` to the columns of `B` in the
coordinate space of `M`:

```text
B' Q_perp = 0
Q_perp' Q_perp = I
F = [Z, M Q_perp]
```

`F` has the same dimension as `M`; injection is replacement, not append. The
source implementation obtains `Q_perp` from a complete QR factorization of
`B` and carries the complement in raw `G/A` coefficients:

```text
G_perp = Q_G + T_G Q_R
A_perp =       T_A Q_R
```

Required geometry checks include:

- full-rank, well-conditioned `B` and its singular spectrum;
- `Z' S_AA Z` near identity;
- `Z' S M Q_perp` near zero;
- sampled/block `F' S F` near identity;
- compact-main `G-R` and `R-R` identity errors;
- protected-span preservation;
- retained/dropped representability and fake-RDM summaries;
- explicit confirmation that rejected directions did not become MWG
  residuals.

## Exact Fixed-Sector One-Body Operators

For an exact one-body operator supplied as Gaussian/gausslet blocks
`O_GG`, `O_GA`, and `O_AA`, form its dense fixed-sector matrix from the raw
coefficients of `[Z, M Q_perp]`:

```text
O_ZZ = Z' O_AA Z
O_ZQ = Z' (O_AG G_perp + O_AA A_perp)
O_QQ = [G_perp,A_perp]' O [G_perp,A_perp]
O_F  = [O_ZZ O_ZQ; O_ZQ' O_QQ]
```

This is an exact one-body congruence in the represented `G/A` space. The
implemented helper applies it to kinetic energy, every uncharged unit nuclear
matrix, and their charge-weighted assembled one-body Hamiltonian:

```text
H1_F = K_F + sum_A Z_A U_A,F
```

The output remains dense and in memory. Validate dimensions, finite values,
symmetry, trace/low-spectrum diagnostics where used, and consistency with the
same geometry that produced `F`.

## Protected-Localized Basis

The fixed sector `F` contains the desired original Gaussian span but is not
itself the site ordering required by the two-index IDA/MWG interaction. Let
`C = B` and use the same `Q_perp` to build the overlap/reconstruction matrix:

```text
A_loc = [C'; Q_perp']
Sbar  = A_loc' A_loc
W     = A_loc Sbar^(-1/2)
L     = F W
```

`Sbar` must be positive definite. `W'W` must be identity to tolerance, and

```text
M' S L = [C, Q_perp] W
```

must remain close to the identity/site correspondence. Report diagonal
deviation, maximum off-diagonal magnitude, and normalized Frobenius deviation
for this localization map.

The exact localized one-body Hamiltonian is:

```text
H1_L = W' H1_F W
```

## Inherited-Site Interaction

The two-index density-density interaction is not a one-body operator. The
accepted protected-localized convention therefore uses:

```text
Vee_L = sym(Vee_M)
```

in the inherited pre-injection `M` site order. The symmetrization is roundoff
cleanup, not a basis rotation. Validate matching dimensions, finite values,
symmetry, and the localization/site diagnostics above. This is the
angular-gausslet interpretation: protected injection
improves the localized one-particle basis while IDA/MWG remains attached to
the localized gausslet-like sites.

A direct interaction congruence such as:

```text
C' V C
```

is rejected, not an alternative implementation. The completed audit showed
that rotating a two-index IDA/MWG density-density matrix this way fails
null/projected many-electron energy invariance even when matrix
back-transforms look accurate. Do not revive it, hide it behind a diagnostic
option, or use it for fixed-density transfer.

Required interaction checks include:

- `H1_L` and `Vee_L` dimensions agree;
- both matrices are finite and symmetric;
- `L` is orthonormal;
- `M' S L` remains close to site identity;
- low one-body modes do not acquire anomalous broad-sector weight;
- bounded physics probes show no residual/broad occupation incentive from the
  inherited-site convention.

## Failure Behavior

The source constructors fail immediately when candidate metadata or dimensions
do not match the residual object, compact source indices are unavailable,
Gaussian cleanup loses all required rank, no broad subspace survives an active
required gate, `Sbar` is not positive definite, or one-body/interaction
matrices are nonfinite or dimensionally inconsistent.

The geometry/localization helpers also return acceptance diagnostics rather
than embedding every physics threshold. A caller must stop acceptance or
endpoint interpretation when the protected span is not preserved, `B` is
rank deficient or insufficiently represented, fixed-sector orthogonality
fails, localized site correspondence is poor, or symmetry diagnostics exceed
their reviewed tolerance.

Do not repair failure by rebuilding `R_compact`, parsing labels, converting a
rejected broad direction into an MWG residual, weakening representability,
dropping a protected direction, rotating `Vee`, or changing artifact/solver
workflow.

## Current Consumer Boundary

The geometry and operator helpers are internal. Their ordinary no-reference
consumer is explicit and default-off; ordinary PQS/RG construction remains
unchanged. The current protected member path in
`src/cartesian_protected_ladder_bundle.jl` composes the staged geometry, exact
fixed-sector one-body matrix, and protected-localized inherited-site
Hamiltonian. Persistence and bundle behavior remain under their separate
authorities and are not specified here.

[Occupied-first injection geometry](occupied_first_injection.md) is a
dependency for identified mandatory HF occupied spaces, but it is not a direct
call-site replacement for this staged `M = [G, R_compact]` construction.
[Protected additive atomic reference correction](protected_additive_reference_correction.md)
owns the separate composition that adds placed packet occupied protection and
screened-Hartree reference fields.

## Ownership And Validation Surfaces

Core source ownership:

- `src/cartesian_residual_gaussians/residual_basis.jl` owns staged protected
  geometry, Gaussian cleanup, representability/fake-RDM selection, and
  geometry diagnostics.
- `src/cartesian_residual_gaussians/augmented_operators.jl` owns exact
  fixed-sector one-body transformation and the localized/inherited-site
  in-memory convention.

Current internal consumer, under separate authority:

- `src/cartesian_protected_ladder_bundle.jl`.

The source/test IDs were validated through bounded default-path smokes and
ignored source-backed replays rather than a dedicated committed unit-test
file. Existing validation/evidence surfaces are:

- `tmp/work/cr2_source_backed_staged_protected_geometry_probe.jl`;
- `tmp/work/cr2_protected_onebody_dense_source_replay.jl`;
- `tmp/work/cr2_protected_vee_algebra_debug.jl`;
- `tmp/work/protected_localized_injection_vee_probe.jl`;
- `docs/src/developer/reports/cr2_staged_subspace_filter_870498b54/README.md`;
- `docs/src/developer/reports/cr2_protected_onebody_audit_eaf05a38c/README.md`;
- manager running-log Passes 235, 253, 259, 269, and 270.

## Explicit Non-Goals

This contract does not own or approve:

- occupied-first selection beyond the dependency link above;
- protected additive-reference construction beyond its dependency link;
- artifact schema, readback, convention-version fields, or row-locality
  metadata, which belong to the
  [protected-localized artifact contract](protected_localized_artifact.md);
- EGOI targets or corrections, which belong to
  [retained-GTO local-product EGOI](retained_gto_egoi.md);
- ladder transfer, cross overlaps, bundles, or restart sidecars, which belong
  to [protected-localized ladder bundles](protected_localized_ladder.md);
- rho0/reference-density history or screened-Hartree formula changes;
- public driver/API/default changes, solver workflow, or endpoint claims;
- new residual policy, fake-RDM hierarchy, mapping defaults, or Cr2-specific
  production behavior.

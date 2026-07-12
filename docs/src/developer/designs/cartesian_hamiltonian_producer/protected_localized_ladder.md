# Protected-Localized Ladder Bundles

Lifecycle:

- `HP-RG-PROTECT-LADDER-XFER-AUDIT-01` is a completed historical
  measurement;
- `HP-RG-PROTECT-LADDER-BUNDLE-FN-01` is an implemented internal opt-in
  facility;
- `HP-RG-PROTECT-LADDER-BUNDLE-TEST-01` is completed acceptance evidence with
  no continuing permission.

This page is the canonical contract for protected-localized ladder transfer
and bundle persistence. The registry owns lifecycle, permission, and source
surfaces. This page owns bundle identity and layout, comparability gates,
transfer direction, ordering, diagnostics, readback, and current limitations.

Protected member numerics and member artifact persistence remain owned by the
[protected-localized basis convention](protected_localized_basis.md) and
[protected-localized artifact contract](protected_localized_artifact.md).
This page does not redefine either contract.

## Purpose

The facility builds a small ladder of comparable protected-localized
Hamiltonians, normally differing in `ns`, and packages the basis-to-basis data
needed to inspect convergence. It replaces repeated consumer-side scripts with
one internal, opt-in directory bundle.

The facility can:

- build protected-localized inherited-site Hamiltonian members;
- prove that adjacent members share the same physical parent lattice and
  supplement;
- compute exact final-basis cross overlaps;
- optionally transfer saved occupied orbitals;
- evaluate the transferred fixed density in the target Hamiltonian;
- write machine-readable sidecars and bounded human-readable summaries.

It does not run UHF or continue a solver.

## Bundle Identity

The implemented bundle identity is:

```text
artifact_kind       = :protected_localized_ladder_bundle
format_version      = 1
convention_id       = :protected_localized_inherited_site_ladder_v1
convention_version  = 1
```

No alternate bundle convention is implemented. The manifest reader rejects a
wrong `artifact_kind` and unsupported `format_version`, and returns the stored
`convention_id` and `convention_version`. Current source does not separately
reject mismatched convention fields; this is a readback-hardening limitation,
not permission to interpret another convention as compatible.

## Directory Layout

The implemented directory shape is:

```text
bundle_root/
  manifest.jld2
  members/nsN/protected_localized_hamiltonian.jld2
  transfers/S_nsB_nsA.jld2
  restarts/nsB_from_nsA_occupied_orbitals.jld2   optional
  summaries/ladder_members.tsv
  summaries/parent_lattice.tsv
  summaries/transfers.tsv
  summaries/restarts.tsv                          optional
  summaries/stages.tsv
```

The manifest records bundle identity, source/current commit provenance,
requested and materialized `ns` values, member/transfer/restart paths,
geometry, electron counts, basis controls, parent-lattice proof, summary
location, and one Hamiltonian-wide Coulomb-expansion summary.

Each member is an ordinary
`:protected_localized_inherited_site_ida_hamiltonian` artifact. Its identity,
matrices, native ordering, row locality, sectors, compatibility, and reader
behavior are defined only by the
[protected-localized artifact contract](protected_localized_artifact.md).

## Comparability Requirements

For every adjacent source/target pair, require all of the following before
constructing a transfer sidecar:

1. The actual parent-axis center arrays have identical counts and agree on
   every axis within the internal parent proof tolerance.
2. Supplement labels and centers are identical in the same order.
3. Geometry, basis controls, and supplement construction come from the shared
   ladder recipe.
4. Every member resolves to the same Coulomb-expansion summary.

A parent-lattice, supplement, or Coulomb-expansion mismatch is a hard failure.
Matching final dimensions or row indices is not a substitute for these proofs.

The current implementation computes transfers only for adjacent entries in
the caller-supplied `ns_values` sequence.

## Exact Final-Basis Transfer

For source final basis `L_A` and target final basis `L_B`, the sidecar stores

```text
S_BA = <L_B | L_A>
```

with target rows and source columns. Orbital transfer is exactly

```text
C_B = S_BA * C_A
```

The cross overlap is assembled from the reconstructed protected in-memory
representations, including terminal gausslet cross action, mixed
gausslet-supplement overlaps, supplement cross overlap, and the source/target
raw protected-localized coefficient blocks.

Final working bases are orthonormal. Their self-overlaps are construction or
debug diagnostics only and are not normal transfer inputs. Do not introduce a
generalized-overlap solve, inverse self metric, or generalized final-basis
workflow.

## Hamiltonian Evaluation Boundary

Transfer moves coefficients only. It does not move a Hamiltonian or
interaction:

- do not transform source `H1_L` into the target basis;
- do not transform source `Vee_L`;
- do not use `C' * V * C` or any interaction rotation.

When saved orbitals are supplied, fixed-density evaluation uses the target
member's native-order `H1_L` and inherited-site `Vee_L`. The source state
contributes only the transferred occupied coefficients. The reported total
also uses the target geometry's nuclear repulsion.

Transfer trace and orthogonality diagnostics determine whether that fixed
density is physically interpretable. The bundle helper records the values; it
does not impose a solver-quality acceptance threshold or run continuation.

## Ordering And Restarts

Member matrices, cross-overlap rows/columns, transferred coefficients, and
new restart sidecars use native protected-localized order.

An input orbital file may declare:

- `:native` or `:inherited_M_site_order`; or
- `:z_order`, `:z`, or `:bond_axis_z`.

For z-ordered source coefficients, convert to native source order with the
source member's row-locality inverse:

```text
C_native[i, :] = C_z[native_to_z_order[i], :]
```

Unknown order labels are rejected. Restart sidecars written by the facility
always record `site_order_kind = "native"` and contain target-native `psiup`
and `psidn`.

The z permutations remain metadata. They never reorder canonical member
matrices or change the transfer definition.

## Sidecar Contracts

Each cross-overlap sidecar uses:

```text
artifact_kind   = :protected_localized_ladder_cross_overlap
format_version  = 1
convention_id   = :protected_localized_inherited_site_ladder_v1
source_ns
target_ns
S_BA
```

Each optional restart sidecar uses:

```text
artifact_kind   = :protected_localized_ladder_transferred_orbitals
format_version  = 1
convention_id   = :protected_localized_inherited_site_ladder_v1
source_ns
target_ns
source_orbitals
site_order_kind = "native"
psiup
psidn
diagnostics/*
fixed_density_energy/*
```

These are ladder-owned sidecars, not fields added to protected member
artifacts. The separately governed external-GTO representation sidecar in
[External GTO orbital import](external_gto_orbital_import.md) is not a ladder
transfer/restart object and does not change this layout or identity.

## Required Diagnostics

Member summaries report:

- `ns`, member path, base/compact/final dimensions, and protected/broad counts;
- `B_min`, median, maximum, and threshold counts;
- `H1_L` / `Vee_L` finiteness and symmetry.

Parent proof reports, for each axis:

- source and target counts;
- equal-count result;
- maximum center difference;
- pass/fail status.

Transfer summaries report:

- source and target `ns`;
- cross-overlap dimensions;
- identical-supplement result;
- minimum, median, and maximum singular values;
- counts below the documented singular-value thresholds.

When orbitals are transferred, also report:

- source and target alpha/beta traces;
- alpha/beta trace loss;
- target occupied-block orthogonality errors;
- target fixed-density one-body, Hartree, exchange, electronic, nuclear, and
  total energies;
- source-orbital and restart paths.

The bundle records stage time, allocation, and GC summaries plus source
artifact, source commit, current commit, geometry, basis controls, electron
counts, and Coulomb-expansion provenance.

## Readback

`read_protected_localized_ladder_bundle(...)` accepts either the bundle
directory or its manifest path. It always returns manifest metadata and may
optionally load:

- protected member artifacts through the protected member reader;
- cross-overlap matrices with source/target `ns`;
- restart alpha/beta coefficient matrices with source/target `ns`.

The current lightweight sidecar readers expose their payloads but do not
repeat every manifest identity, dimension, or convention check. Conforming
consumers use the manifest and canonical v1 identities together. This current
readback boundary is not authority for malformed or alternate-convention
sidecars.

## Current Representation Limitation

Protected member artifacts persist `H1_L`, `Vee_L`, locality, sectors, and
provenance, but not the raw protected-localized `G_L` / `A_L` coefficients
needed to assemble a new exact cross overlap.

Therefore, cross-overlap construction currently requires reconstructing both
protected members in memory from the recipe. A bundle cannot derive a new
transfer from member artifacts alone. Raw-`L` coefficient fields remain
forbidden. `HP-REP-XGTO-PROTECT-SIDECAR-*` separately implements saving an exact
external `S_LG` once while one in-memory member exists; it does not let the
bundle reconstruct other overlaps or add member fields.

## Failure Behavior

Stop rather than writing a transfer when:

- parent lattice proof fails;
- supplement labels, centers, or ordering differ;
- member Coulomb-expansion summaries differ;
- exact in-memory protected representations are unavailable;
- source orbital dimensions do not match the source final basis;
- source site order is unknown;
- member artifact readback rejects a member;
- manifest kind or format is unsupported.

Large transfer loss, poor singular values, or bad occupied orthogonality make
the transferred state unsuitable for physical interpretation. They do not
authorize generalized overlap, Hamiltonian transformation, or interaction
rotation.

## Ownership And Evidence

Primary implemented owner:

- `src/cartesian_protected_ladder_bundle.jl`.

Related existing owners:

- `src/cartesian_representation_transfer.jl` for the general final-basis
  cross-overlap transfer contract;
- `src/cartesian_ida_hamiltonian.jl` for protected member artifacts.

Implementation evidence:

- source commit `3eaa812a9`;
- manager running-log Pass 313.

The accepted H2 smoke and ignored Cr2 replay implement
`HP-RG-PROTECT-LADDER-BUNDLE-TEST-01`. No committed large Cr2 fixture or
production HF requirement belongs to this validation contract.

## Explicit Non-Goals

This contract does not approve:

- changes to protected member artifact identity or schema;
- EGOI targets, corrections, or artifacts;
- additive-reference or screened-Hartree behavior;
- rho0/reference-density changes;
- solver or UHF continuation;
- source-Hamiltonian or source-`Vee` transformation;
- generalized final-basis overlap transfer;
- ladder-owned representation sidecars or raw-`L` member fields; the separate
  external-GTO sidecar remains outside this contract;
- package exports, public API, or driver defaults;
- Cr2 production claims.

## Historical Audit

`HP-RG-PROTECT-LADDER-XFER-AUDIT-01` is completed historical evidence. Its
same-parent transfer measurements established the cross-overlap-only and
target-Hamiltonian evaluation conventions later implemented by the bundle.
Numerical paths and detailed values remain in manager running-log history, not
this canonical contract.

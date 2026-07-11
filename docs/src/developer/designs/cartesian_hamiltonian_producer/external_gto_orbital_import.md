# External GTO Orbital Import

Status: the general import facility and narrow protected-localized composition
are implemented under `HP-REP-XGTO-IMPORT-FN-01` and
`HP-REP-XGTO-IMPORT-TEST-01`. The durable representation sidecar is implemented
under
`HP-REP-XGTO-PROTECT-SIDECAR-FN-01` and
`HP-REP-XGTO-PROTECT-SIDECAR-TEST-01`.

This is representation-transfer infrastructure. It is not a physics endpoint,
Hamiltonian correction, or solver workflow.

## Purpose

Build exact overlaps between a GaussletBases final working basis and an
explicit external Gaussian AO basis, then import external HF orbitals by the
normal final-basis transfer rule:

```text
S_FG = <final basis F | external GTO basis G>
C_F  = S_FG * C_G
```

`C_G` is the external AO coefficient matrix, for example from PySCF. The final
working basis is intended to be orthonormal, so import uses only the cross
overlap. This is the same final-basis transfer convention recorded in the repo
overlap policy.

The external AO self-overlap is validation data only:

```text
C_G' * S_GG * C_G ~= I
```

Do not use `S_GG` to build a generalized final-basis solve or a generalized
downstream final-basis transfer.

## Approved IDs

- `HP-REP-XGTO-IMPORT-FN-01`
- `HP-REP-XGTO-IMPORT-TEST-01`
- `HP-REP-XGTO-PROTECT-SIDECAR-FN-01`
- `HP-REP-XGTO-PROTECT-SIDECAR-TEST-01`

## Source Ownership

Approved source surface:

- `src/cartesian_external_gto_import.jl`
- `src/GaussletBases.jl` for include/export wiring
- `src/cartesian_representation_transfer.jl` only for shared transfer
  diagnostics, if needed
- `src/cartesian_gto_probes.jl` only for narrow reuse around
  `gto_overlap_matrix(...)` and
  `_cartesian_final_gto_cross_overlap_handoff(...)`, without changing their
  numerical contracts
- `src/cartesian_protected_ladder_bundle.jl` only if a small internal accessor
  is required to expose an already-built member's `base`, supplement, `G_L`,
  and `A_L`; no ladder manifest or sidecar contract may change

The implementation should be a small wrapper/facility around existing exact
overlap kernels, especially:

```text
gto_overlap_matrix(working, probes)
```

Do not put this under screened-Hartree, EGOI, residual-GTO selection, or a
Cr2-specific workflow.

Allowed compact API/result shapes:

- a read/write helper for explicit external GTO orbital packets;
- a compact packet/result/diagnostic object for imported orbitals;
- a direct import helper such as `import_external_gto_orbitals(working, packet)`;
- an internal protected-member composition helper such as
  `protected_localized_external_gto_import(member, packet)`;
- optional spin-resolved alpha/beta output when the packet supplies both spin
  channels.

These shapes must stay focused on orbital import and capture diagnostics. They
must not become Hamiltonian, solver, or interaction-transform payloads.

## External Packet Contract

The repo reader consumes an explicit external basis/state packet. It must not
depend on PySCF at test time and must not infer AO order from label strings
when explicit basis data is available.

Required packet contents:

- centers;
- resolved angular powers or an already-resolved Cartesian/spherical
  convention;
- exponents;
- contraction coefficients;
- AO labels;
- AO ordering/fingerprint;
- `S_GG` and its fingerprint;
- spin-resolved MO coefficient matrices where applicable;
- occupations;
- provenance, including source code/tool, source basis, molecule/atom
  geometry, and ordering convention.

If the source calculation used spherical harmonics, a CR2/PySCF-side writer may
resolve the external basis into explicit Cartesian/GTO coefficient data before
the packet reaches this repo. The repo reader should validate the resolved
data; it should not guess spherical ordering from strings as the accepted
construction rule.

## Import Diagnostics

The import facility must define/report:

- `S_FG = <F|G_external>`;
- `S_GG = <G_external|G_external>` for validation only;
- `C_F_alpha = S_FG * C_G_alpha`;
- `C_F_beta = S_FG * C_G_beta`, when beta orbitals are present;
- `C_G' * S_GG * C_G` orthogonality errors;
- `C_F' * C_F` capture/norm diagnostics;
- density trace capture/loss by spin;
- worst orbital capture;
- labels/fingerprint/order checks;
- packet provenance and basis dimensions.

Capture less than identity is a projection/capture diagnostic. It is not by
itself a request to transform the source Hamiltonian or to introduce a
generalized final-basis overlap workflow.

## Protected-Localized Composition

For one already-built protected-localized member, form

```text
raw_to_L = [member.raw.G_L; member.raw.A_L]
S_LG     = <L | G_external>
C_L      = S_LG * C_G
```

by calling the existing exact
`_cartesian_final_gto_cross_overlap_handoff(...)` with the member's fixed/base
working basis, supplement, `raw_to_L`, and packet probes. `S_LG` has native
protected-localized rows and external AO columns. The helper must verify the
handoff orientation `:final_by_gto`, finiteness, final dimension, external AO
dimension, and packet identity before importing either spin block.

The returned alpha/beta `C_L` matrices are the direct projected coefficients.
They are deliberately not orthonormalized or converted into solver restart
orbitals. A consumer may use the same projected state for screened and
unscreened calculations, but solver-ready orthonormalization and all HF state
management remain consumer responsibilities.

This composition requires the in-memory member because protected Hamiltonian
artifacts intentionally omit `G_L/A_L`. It must not reconstruct `S_LG` from
`H1_L`, `Vee_L`, row centers, labels, or sector metadata.

## Metric-Aware Principal Capture

The external AO basis is nonorthogonal. Full external-subspace capture must
therefore use `S_GG`, not raw AO coefficient squares or `svdvals(S_LG)` alone.
For the diagnostic eigendecomposition

```text
S_GG = U * Diagonal(s) * U'
tau  = max(1e-12, 1e-10 * maximum(s))
Q_G  = U_keep * Diagonal(s_keep^(-1/2))
T_LG = S_LG * Q_G
```

report the retained source-metric rank, discarded/tiny eigenvalue count,
source-metric eigenvalue range, principal singular values `svdvals(T_LG)`, and
projected Gram eigenvalues `eigvals(T_LG' * T_LG)`. Eigenvalues below `-tau`
or capture outside `[-1e-8, 1 + 1e-8]` fail. Tiny roundoff may be clamped for
reporting only. Tiny positive metric modes below `tau` are excluded only from
this diagnostic and are reported; they do not change `C_L = S_LG*C_G`.

For each occupied spin block, the packet has already established
`C_G' * S_GG * C_G = I`. Its metric-aware capture is therefore directly

```text
K_occ = C_L' * C_L
```

with per-orbital diagonal captures, occupied principal singular values,
projected Gram eigenvalues, and occupation-weighted trace loss. Do not define
shell/angular capture from raw squared AO coefficients. A later angular
decomposition must be separately designed in the external AO metric.

## Protected Representation Sidecar

The sidecar is a standalone representation artifact. It is not a field of the
protected Hamiltonian artifact and not a ladder transfer or restart sidecar.
Its fixed identity is:

```text
artifact_kind       = :protected_localized_external_gto_representation
format_version      = 1
convention_id       = :protected_localized_external_gto_native_v1
convention_version  = 1
site_order_kind     = :native
orientation         = :final_by_external
```

The v1 JLD2 layout uses these required groups:

```text
artifact_kind
format_version
convention_id
convention_version
site_order_kind
orientation

cross_overlap/S_LG
cross_overlap/fingerprint_sha256
cross_overlap/final_dimension
cross_overlap/external_dimension

external/ao_count
external/ao_labels
external/ordering_fingerprint_sha256
external/S_GG_fingerprint_sha256
external/alpha_coefficients_fingerprint_sha256
external/beta_coefficients_fingerprint_sha256       optional with beta
external/provenance/*

imported/alpha/coefficients
imported/alpha/occupations
imported/beta/coefficients                           optional with beta
imported/beta/occupations                            optional with beta

protected/final_dimension
protected/artifact_kind
protected/convention_id
protected/recipe_fingerprint_sha256
protected/H1_L_fingerprint_sha256
protected/Vee_L_fingerprint_sha256
protected/source_artifact
protected/member_artifact                           optional
protected/source_commit
protected/current_commit
protected/basis_controls/*
protected/geometry_inputs/*
coulomb_expansion/*

diagnostics/cross_overlap/*
diagnostics/source_metric/*
diagnostics/alpha/*
diagnostics/beta/*                                  optional with beta
```

Matrix fingerprints use the same deterministic Float64 byte convention as the
existing external overlap fingerprint. Optional beta groups are all-or-none.
Unknown extra keys do not change v1 semantics; missing required keys fail.
The ordering and `S_GG` fingerprints define external-basis compatibility for
later imports. Source coefficient fingerprints bind the stored alpha/beta
imports to the saved packet state; they do not prevent applying the same
`S_LG` to a different coefficient set with the same validated AO basis and
fresh capture diagnostics.

Required numerical payload:

- `S_LG`, stored explicitly with final rows and external AO columns;
- direct imported alpha coefficients and occupations;
- direct imported beta coefficients and occupations when present;
- no orthonormalized solver orbitals, Hamiltonian, or interaction matrix.

Required external identity:

- AO count and ordered AO labels;
- packet ordering fingerprint and `S_GG` fingerprint;
- source alpha/beta coefficient fingerprints;
- bounded packet/source provenance.

Required protected-member identity:

- final dimension and native-order convention;
- protected member convention identity and paired artifact path when
  available;
- recipe/source/current commit provenance and a compact recipe fingerprint;
- basis and geometry controls;
- Hamiltonian-wide Coulomb policy summary;
- `H1_L` and `Vee_L` fingerprints used only to bind the sidecar to the member,
  not to transform either matrix.

Required diagnostics:

- `S_LG` dimensions, finiteness, maximum magnitude, and fingerprint;
- packet ordering, `S_GG`, and source-orthogonality validation;
- per-spin capture matrices, per-orbital captures, occupied principal values,
  source/captured density traces, and trace losses;
- full external-metric rank/eigenvalue and principal-capture diagnostics;
- native-order and protected-member identity checks.

Readback must reject unrecognized identity/version/order/orientation,
nonfinite or dimensionally inconsistent matrices, changed `S_LG` fingerprint,
malformed spin blocks, or inconsistent diagnostic lengths. When a packet or
protected member/artifact is supplied for validation, its fingerprints,
dimensions, native ordering, recipe/member identity, and Coulomb policy must
match. When the original external packet is supplied, recompute and validate
`C_L = S_LG*C_G` for every stored spin block. A later packet/state with the
same AO ordering and `S_GG` may instead produce a new direct import from the
saved `S_LG`; it does not overwrite or relabel the stored import. Missing raw
`G_L/A_L` in a protected artifact is not a readback error; the saved `S_LG` is
precisely the durable representation object.

At the current Cr2 dimensions, `6945 x 448` Float64 `S_LG` storage is about
25 MB and is intentionally retained. Do not replace it with reconstruction
metadata or lossy compression in this lane.

## Consumer Boundary

The repo owns exact representation transfer, packet/member identity, capture
diagnostics, and sidecar persistence. Consumers such as CR2 own
orthonormalization into solver-ready orbitals, screened/unscreened HF starts,
per-sweep occupied-subspace overlap logging, energies, spin diagnostics,
owner-local moments, compact-residual occupation, and physical interpretation.
None of those consumer responsibilities belongs in this sidecar helper.

## Tests

Approved test surface:

- `test/nested/cartesian_external_gto_import_runtests.jl`

Tests should be correctness-only and small. They should not require PySCF at
test time. Synthetic or repo-built explicit external packets are acceptable.

Required test coverage:

- exact overlap shape, finiteness, and validation of `S_GG`;
- `C_G' * S_GG * C_G` validation;
- `C_F = S_FG * C_G` import;
- capture diagnostics;
- invariance under occupied-space rotations;
- mismatch/fingerprint/order failure;
- spin-resolved alpha/beta handling if present.

The protected extension and sidecar tests must additionally cover:

- exact parity with a direct
  `_cartesian_final_gto_cross_overlap_handoff(...)` call;
- native `S_LG` orientation and `C_L = S_LG*C_G` for both spins;
- metric-aware full/occupied capture diagnostics and occupied-space rotation
  invariance;
- temporary-file sidecar roundtrip;
- rejection of wrong kind/version/order/orientation, malformed dimensions,
  changed fingerprints, and mismatched packet/member identity;
- absence of `G_L`, `A_L`, `H1_L`, and `Vee_L` payloads.

No energy assertions, SCF convergence claims, molecule-specific physics gates,
Cr2-sized committed fixtures, or Cr2 production gates are approved.

## Explicit Exclusions

Forbidden in this lane:

- Hamiltonian transforms;
- `C' V C`;
- `Vee` or source-interaction transforms;
- generalized final-basis overlap workflow;
- raw-coefficient `s/p/d/f+` capture accounting;
- solver workflow;
- screened-Hartree changes;
- EGOI changes;
- residual selection or injection policy changes;
- production Cr2 claims;
- PySCF dependency in repo tests;
- protected Hamiltonian artifact fields for `G_L/A_L` or `S_LG`;
- ladder manifest, transfer-sidecar, or restart-sidecar changes;
- reuse of ladder sidecar identities for external-GTO representation.

The implementation remains in the existing source and test files and owns one
compact internal protected-import/readback result shape. It adds no module,
source file, public export, driver input, protected/ladder artifact schema, or
fixture. Persistent raw protected coefficients, a generalized final metric,
or changes outside the approved source surface remain forbidden.

## Decision Rule

If the import can be implemented as a bounded packet reader plus
`gto_overlap_matrix`-based transfer with clear capture diagnostics, proceed
under this authority. Protected composition may additionally consume one
already-built member and persist the exact native `S_LG` sidecar above. If
correct import requires Hamiltonian transformation, interaction rotation,
generalized final-basis metric logic, reconstruction from a protected artifact,
or PySCF-dependent repo tests, stop and request a new design amendment.

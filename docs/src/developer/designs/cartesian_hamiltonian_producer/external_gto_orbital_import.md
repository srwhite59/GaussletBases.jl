# External GTO Orbital Import

Status: approved source/design authority under `HP-REP-XGTO-IMPORT-FN-01`
and `HP-REP-XGTO-IMPORT-TEST-01`.

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

## Source Ownership

Approved source surface:

- `src/cartesian_external_gto_import.jl`
- `src/GaussletBases.jl` for include/export wiring
- `src/cartesian_representation_transfer.jl` only for shared transfer
  diagnostics, if needed
- `src/cartesian_gto_probes.jl` only for narrow reuse around
  `gto_overlap_matrix(...)`, without changing its numerical contract

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

No energy assertions, SCF convergence claims, molecule-specific physics gates,
or Cr2 production gates are approved.

## Explicit Exclusions

Forbidden in this lane:

- Hamiltonian transforms;
- `C' V C`;
- `Vee` or source-interaction transforms;
- generalized final-basis overlap workflow;
- solver workflow;
- screened-Hartree changes;
- EGOI changes;
- residual selection or injection policy changes;
- production Cr2 claims;
- PySCF dependency in repo tests.

## Decision Rule

If the import can be implemented as a bounded packet reader plus
`gto_overlap_matrix`-based transfer with clear capture diagnostics, proceed
under this authority. If correct import requires Hamiltonian transformation,
interaction rotation, generalized final-basis metric logic, or PySCF-dependent
repo tests, stop and request a new design amendment.

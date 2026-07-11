# External GTO Orbital Import

The general external-GTO importer and the protected-localized representation
sidecar are implemented representation infrastructure. They do not own a
Hamiltonian correction, solver state, or physics endpoint.

## Lifecycle

`HP-REP-XGTO-IMPORT-FN-01` / `TEST-01` and
`HP-REP-XGTO-PROTECT-SIDECAR-FN-01` / `TEST-01` are implemented. The general
importer landed in `0452ab581`; protected composition and the sidecar landed in
`702aa4a62`. Manager-log Passes 337, 363, and 365 record acceptance.

## Transfer Convention

For an orthonormal final working basis `F` and explicit external Gaussian AO
basis `G`, import uses only the cross overlap:

```text
S_FG = <F|G_external>
C_F  = S_FG * C_G
```

The external self-overlap validates the packet and source orbitals:

```text
C_G' * S_GG * C_G ~= I
```

`S_GG` is not a final-basis metric. It is never used to construct a
generalized final-basis solve, source-Hamiltonian transform, or interaction
transform. Capture below one reports representation loss.

## Implemented Ownership And Entry Points

`src/cartesian_external_gto_import.jl` owns all packet, result, protected
composition, sidecar, and saved-representation logic. `src/GaussletBases.jl`
exports the general packet/import API and includes the owner file.

Exported entry points are `ExternalGTOOrbitalSpinBlock`,
`ExternalGTOOrbitalPacket`, `ExternalGTOOrbitalSpinImport`,
`ExternalGTOOrbitalImportResult`, `external_gto_overlap_fingerprint`,
`external_gto_ordering_fingerprint`, and `import_external_gto_orbitals`.
Protected entry points remain module-qualified and unexported:
`protected_localized_external_gto_import`,
`write_protected_localized_external_gto_representation`,
`read_protected_localized_external_gto_representation`, and
`external_gto_import_from_saved_representation`.

`src/cartesian_gto_probes.jl` owns `gto_overlap_matrix(...)` and the reused
exact `_cartesian_final_gto_cross_overlap_handoff(...)`.
`src/cartesian_representation_transfer.jl` records the related cross-overlap
final-basis convention but is not called by this implementation.
`src/cartesian_protected_ladder_bundle.jl` supplies the existing in-memory
protected member shape and commit helper; no ladder accessor or schema was
added.

There is no external packet file reader, PySCF parser, or external-file format
in this facility. Callers construct the explicit packet in memory.

## External Packet

`ExternalGTOOrbitalSpinBlock` stores exactly `spin`, `coefficients`, and
`occupations`.

`spin` is `:restricted`, `:alpha`, or `:beta`. A packet without beta accepts a
restricted or alpha block. A spin-resolved packet requires alpha and beta roles
exactly. Coefficient rows must match AO count; orbital columns must match the
occupation vector; all values must be finite.

`ExternalGTOOrbitalPacket` stores exactly:

```text
probes
centers
angular_powers
exponents
contraction_coefficients
primitive_normalizations
ao_labels
ordering_fingerprint
S_GG
S_GG_fingerprint
alpha
beta
provenance
```

The constructor receives an explicit resolved Cartesian supplement
representation, `S_GG`, and spin blocks. It derives center, angular-power,
primitive, contraction, and normalization arrays from `probes`. `provenance`
is a caller-owned compact value; the repo does not infer an external AO order
or spherical convention from label strings.

The ordering fingerprint covers ordered labels, powers, centers, exponents,
coefficients, and primitive normalization. The overlap fingerprint is over the
Float64 matrix bytes. Import recomputes both fingerprints by default, checks
`S_GG` finiteness and symmetry, and recomputes `<probes|probes>` numerically.
The existing `validate_fingerprints=false` option reports fingerprint booleans
without making them fatal; it does not bypass the physical `S_GG` checks.
Protected composition always uses strict fingerprint validation.

## General Import Result

`import_external_gto_orbitals(working, packet)` validates the packet before
calling `gto_overlap_matrix(working, packet.probes)`. It returns
`ExternalGTOOrbitalImportResult` with:

```text
cross_overlap_size
cross_overlap_finite
ordering_fingerprint_valid
S_GG_fingerprint_valid
S_GG_symmetry_error
S_GG_expected_error
alpha
beta
```

The general result does not retain `S_FG`. Each
`ExternalGTOOrbitalSpinImport` stores:

```text
spin
imported_coefficients
source_orthogonality_error
capture_matrix
orbital_captures
density_trace_source
density_trace_capture
density_trace_loss
worst_orbital_capture
```

Here `imported_coefficients = S_FG*C_G`,
`capture_matrix = C_F'*C_F`, and orbital captures are its diagonal. Density
trace values use the packet occupations. Source MOs that fail the configured
`C_G'*S_GG*C_G` orthogonality tolerance are rejected. The returned coefficients
are direct projections and are not solver-orthonormalized.

## Protected-Localized Composition

For an already-built protected member, the internal composition forms:

```text
raw_to_L = [member.raw.G_L; member.raw.A_L]
S_LG     = <L|G_external>
C_L      = S_LG * C_G
```

The exact handoff must report `:final_by_gto`; `S_LG` rows must equal the
member's native protected dimension and columns must equal external AO count.
The returned `ProtectedLocalizedExternalGTORepresentation` stores exactly:

```text
cross_overlap
imported
source_metric
occupied_capture
identity
```

It does not orthonormalize alpha or beta imports. Construction requires the
in-memory member because protected Hamiltonian artifacts intentionally omit
`G_L/A_L`. It never reconstructs `S_LG` from matrices, centers, labels, or
sector metadata.

### Metric-aware capture

The full nonorthogonal external span is whitened only for diagnostics:

```text
S_GG = U * Diagonal(s) * U'
tau  = max(1e-12, 1e-10 * maximum(s))
Q_G  = U_keep * Diagonal(s_keep^(-1/2))
T_LG = S_LG * Q_G
```

`source_metric` reports retained rank, discarded count, `tau`, metric
eigenvalue extrema, `svdvals(T_LG)`, and eigenvalues of `T_LG'*T_LG`.
Materially negative metric eigenvalues or capture outside
`[-1e-8, 1+1e-8]` fail. Tiny modes are excluded only from this diagnostic; they
do not alter `C_L = S_LG*C_G`.

For each source-metric-orthonormal spin block, `occupied_capture` reports the
principal singular values and projected-Gram eigenvalues of `C_L` directly.
Raw squared AO coefficients do not define angular-channel capture.

## Protected Representation Sidecar

The standalone sidecar identity is fixed:

```text
artifact_kind      = :protected_localized_external_gto_representation
format_version     = 1
convention_id      = :protected_localized_external_gto_native_v1
convention_version = 1
site_order_kind    = :native
orientation        = :final_by_external
```

`S_LG` is stored with final rows and external AO columns. The complete live v1
key families are:

```text
cross_overlap/{S_LG,fingerprint_sha256,final_dimension,external_dimension}

external/{ao_count,ao_labels,ordering_fingerprint_sha256,
          S_GG_fingerprint_sha256,alpha_coefficients_fingerprint_sha256}
external/beta_coefficients_fingerprint_sha256              optional
external/provenance/repr

imported/alpha/{coefficients,occupations}
imported/beta/{coefficients,occupations}                    optional

protected/{final_dimension,artifact_kind,convention_id,
           recipe_fingerprint_sha256,H1_L_fingerprint_sha256,
           Vee_L_fingerprint_sha256,source_artifact,
           source_commit,current_commit}
protected/member_artifact                                  optional
protected/basis_controls/{nesting,ns,core_spacing,s_factor,basisname,lmax}
protected/geometry_inputs/{atom_symbols,nuclear_charges,atom_locations,nup,ndn}

coulomb_expansion/{policy,doacc,term_count,del,s,c,maxu}

diagnostics/cross_overlap/{final_dimension,external_dimension,finite,max_abs,
                           fingerprint_sha256}
diagnostics/source_metric/{retained_rank,discarded_count,tau,
                           metric_eigenvalue_min,metric_eigenvalue_max,
                           projected_gram_eigenvalues,principal_singular_values,
                           S_GG_symmetry_error,S_GG_expected_error}
diagnostics/alpha/{spin,source_orthogonality_error,capture_matrix,
                   orbital_captures,density_trace_source,density_trace_capture,
                   density_trace_loss,worst_orbital_capture,
                   projected_gram_eigenvalues,principal_singular_values}
diagnostics/beta/{same fields as alpha}                     optional
```

The Coulomb group delegates to the current centralized summary serializer; it
does not define Coulomb policy. Optional beta keys are all-or-none. The sidecar
stores no external `S_GG` matrix, `G_L`, `A_L`, `H1_L`, or `Vee_L` payload.
Matrix fingerprints use the deterministic Float64 byte convention from
`external_gto_overlap_fingerprint`; H1/Vee fingerprints bind member identity
only.

## Read, Write, And Reimport

The writer verifies that the in-memory `S_LG` still matches its fingerprint.
Readback requires every applicable v1 key and rejects wrong identity, version,
ordering, orientation, dimensions, nonfinite values, changed fingerprints,
malformed spin roles, partial beta groups, incomplete Coulomb provenance, or
inconsistent capture diagnostics.

Optional readback validation may supply:

- the original packet, which also verifies source coefficient/occupation
  identity and recomputes every stored import;
- the in-memory protected member;
- a protected Hamiltonian artifact, validated by content so relocation alone
  does not fail.

`external_gto_import_from_saved_representation(value, packet)` accepts a new
orbital state only when ordered labels, ordering fingerprint, and `S_GG`
fingerprint match the saved external AO basis. It then returns fresh direct
imports and capture diagnostics without overwriting the sidecar's original
state. Unknown extra keys are not interpreted as v1 semantics.

This sidecar is not a protected Hamiltonian field and not a ladder manifest,
transfer sidecar, or restart sidecar. Those identities and schemas remain
unchanged.

## Validation And Failure Boundary

The committed test owner is:

```text
test/nested/cartesian_external_gto_import_runtests.jl
```

It covers generic restricted/spin-resolved imports, fingerprint/order/metric
failures, occupied rotations, exact protected handoff parity, rectangular
capture dimensions, alpha/beta direct imports, sidecar key inventory and
roundtrip, saved-overlap reimport, packet/member/artifact validation, relocated
artifact acceptance, and payload tampering. It uses synthetic repo-owned
packets and no PySCF dependency or molecular endpoint.

The facility must not add or perform:

- packet-file or PySCF readers;
- solver-ready orbital orthonormalization or HF state management;
- Hamiltonian, source-Hamiltonian, or `Vee` transforms;
- direct `C' V C` or a generalized final-basis metric workflow;
- raw-coefficient angular capture claims;
- protected artifact or ladder-sidecar schema changes;
- screened-Hartree, EGOI, residual-selection, injection, or Cr2 physics claims.

Consumers own solver initialization, orthonormalization, sweep diagnostics,
energies, spin interpretation, and endpoint conclusions.

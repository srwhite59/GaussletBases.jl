# External Cartesian GTO Interchange

This contract owns a small public interchange between an external real
Cartesian Gaussian AO calculation and the existing GaussletBases
`ExternalGTOOrbitalPacket` importer. It removes recurring AO-order and
normalization reconciliation from consumer scripts without adding another
Gaussian overlap implementation or a quantum-chemistry driver.

## Lifecycle And Ownership

The implementation is approved and pending under:

- `HP-REP-XGTO-INTERCHANGE-FN-01` for the versioned bundle reader and
  convention reconciliation;
- `HP-REP-XGTO-PYSCF-EXPORT-FN-01` for one optional thin PySCF exporter;
- `HP-REP-XGTO-CLOSESTDET-FN-01` for the explicit closest-determinant
  operation; and
- `HP-REP-XGTO-INTERCHANGE-TEST-01` for the frozen fixture and bounded public
  validation.

These are new records. They do not broaden maintenance authority for
`HP-REP-XGTO-IMPORT-*`, which continues to own the existing in-memory packet
and import operation.

## User Workflow

An external user exports one already-computed PySCF orbital state:

```text
python bin/export_pyscf_cartesian_gto.py \
    --checkpoint calculation.chk \
    --output external_gto_bundle.toml
```

Julia then reads a validated packet and uses the existing importer:

```julia
packet = read_external_cartesian_gto_packet("external_gto_bundle.toml")
imported = import_external_gto_orbitals(working, packet)
```

The direct import remains the raw projection

```text
S_FG = <F|G_external>
C_F  = S_FG * C_G
```

and is not silently orthonormalized. A determinant consumer may request a
separate same-construction operation:

```julia
start = closest_external_gto_determinant(
    working,
    packet;
    minimum_gram_eigenvalue = 0.999,
)
```

The caller chooses the minimum acceptable occupied-subspace capture. No
default capture threshold is provided.

## Public Julia Surface

Implementation may add exactly these two root exports:

```julia
read_external_cartesian_gto_packet(path; overlap_atol=1e-10,
    overlap_rtol=1e-10)::ExternalGTOOrbitalPacket

closest_external_gto_determinant(working, packet;
    minimum_gram_eigenvalue)
```

The reader returns the existing packet type. The closest-determinant operation
returns a documented, fixed small `NamedTuple`; it does not justify a new
public type. At most one new public type may be proposed only by stopping and
showing that the required validation diagnostics cannot be represented
clearly with the existing packet/import types and this final result record.

The closest result contains the unchanged direct import plus per-spin records
with:

```text
coefficients
gram_eigenvalues
minimum_gram_eigenvalue
principal_angles_radians
maximum_principal_angle_radians
orthonormality_error
```

The top-level result contains `imported`, `alpha`, and optional `beta`. It
does not retain `S_FG`, create solver state, or persist a Hamiltonian.

## Version 1 Bundle

The v1 interchange is one manifest path and one same-stem binary sibling:

```text
external_gto_bundle.toml
external_gto_bundle.f64
```

The manifest records the exact payload basename. Absolute payload paths,
directory traversal, alternate extensions, and multiple payloads are rejected.
Other files near the bundle are irrelevant and are never scanned.

`manifest.toml` is strict, versioned, and reader-independent. It records:

- format kind and version;
- all length and exponent units;
- scalar type, byte order, and array storage order;
- producer and exporter names and versions;
- source-state kind and original spherical/Cartesian representation;
- declared all-electron nuclear-Coulomb Hamiltonian kind, molecular charge,
  and spin sector;
- nuclear charges and positions;
- ordered explicit AO records;
- restricted or alpha/beta orbital-block identities and occupations;
- array shapes, byte offsets, element counts, and SHA-256 hashes;
- whole-payload SHA-256;
- source-checkpoint hash when a checkpoint supplied the state.

Every AO record is present in the exact external Cartesian order and contains:

```text
ao_index_1based
label
atom_index_1based
shell_index_1based
contraction_index_1based
angular_powers = [lx, ly, lz]
center_bohr
exponents_bohr_inverse_square
contraction_coefficients
```

The shell, contraction, and angular fields are descriptive checks on the
ordered record. They are not keys for sorting or reconstructing a permutation.
Basis names are provenance only and never authorize lookup or replacement of
the explicit primitive data.

The `.f64` payload concatenates only the named dense arrays in manifest order:

```text
source_overlap
alpha_coefficients
alpha_occupations
beta_coefficients       # spin-resolved bundles only
beta_occupations        # spin-resolved bundles only
```

Every array is IEEE-754 Float64, little-endian, and column-major. The reader
rejects missing, duplicate, overlapping, out-of-range, nonfinite, wrongly
shaped, wrongly hashed, or trailing payload data. The bundle carries no
kinetic, nuclear-attraction, Fock, ERI, Hamiltonian, density-fit, or solver
array.

Version 1 uses Bohr for centers and inverse-square Bohr for exponents. The
reader requires those exact unit identifiers rather than becoming a general
unit-conversion layer.

## Bounded QCSchema Fit Decision

QCSchema was checked before defining this format. The current official
[wavefunction schema](https://molssi-qc-schema.readthedocs.io/en/latest/auto_wf.html)
states that AO functions use CCA/libint order. Its
[basis definition at the checked revision](https://github.com/MolSSI/QCSchema/blob/5390e6f11d21847e4e7ca2ad14a97594f957cb2d/qcschema/dev/definitions.py)
stores shell exponents, coefficients, and a spherical/Cartesian flag, but does
not specify an exact per-AO normalization/order record. The current
wavefunction fields provide orbital coefficients and occupations but no
standard source-overlap matrix field.

Using QCSchema here would therefore require a nonstandard extension for exact
PySCF order, normalization, and `S_GG`, plus a material QCElemental/QCSchema
dependency or an ambiguous partial parser. The v1 interchange instead uses
the minimal TOML plus deterministic Float64 payload above. This decision does
not authorize general QCSchema reading, writing, or translation.

## PySCF Export Convention

`bin/export_pyscf_cartesian_gto.py` is an optional exporter, not a PySCF
calculation driver. It may load one existing molecular checkpoint and emit the
v1 bundle. PySCF and NumPy remain dependencies of that external command's
environment; neither becomes a Julia package dependency.

For an already-Cartesian PySCF state, the exporter preserves the coefficient
matrix and exact PySCF Cartesian AO order. For a real spherical state it must
use the documented PySCF transformation:

```python
X = mol.cart2sph_coeff(normalized="sp")
C_cart = X @ C_sph
```

PySCF documents `normalized="sp"` as the libcint convention and verifies
`X.T @ S_cart @ X == S_sph` in its
[`cart2sph_coeff` API](https://pyscf.org/_modules/pyscf/gto/mole.html).
The exporter must not use the paper-local global solve
`solve(S_cart, S_cart_sph)` as the production conversion.

The exporter obtains explicit exponents and input contraction coefficients
from the loaded molecule. In particular, it uses `mol.bas_ctr_coeff(ib)` or
the exactly equivalent PySCF operation that removes libcint's radial primitive
normalization before the Julia-side axiswise primitive construction. A raw
libcint environment coefficient is not an equivalent interchange value. The
exporter emits Cartesian labels and powers in the exact live AO sequence and
computes `S_GG` from the same Cartesian realization of the loaded molecule.
It verifies every exported spin block with

```text
C_cart' * S_GG * C_cart ~= I
```

before writing. It exports only columns with strictly positive source
occupation and records their original MO indices and occupations. This is not
an occupation cutoff or determinant claim; fractional occupations remain
representable by the packet but are rejected by the closest-determinant
operation.

The first exporter supports finite, nonperiodic, real RHF and collinear UHF
checkpoint states for an ordinary all-electron nuclear-Coulomb Hamiltonian. It
rejects GHF, spinor, complex, periodic, ECP/pseudopotential, externally modified
one-body, and claimed ROHF-initial-determinant semantics. The manifest records
that supported Hamiltonian kind; it does not store a Hamiltonian matrix. The
exporter invokes no SCF calculation and accepts no basis name as a substitute
for explicit checkpoint basis data.

## Julia Convention Reconciliation

The Julia reader consumes AO records in manifest order. It constructs raw
`CartesianGaussianShellOrbitalRepresentation3D` probes using the existing
axiswise-normalized Cartesian primitive convention and one
`CartesianGaussianShellSupplementRepresentation3D` with
`supplement_kind = :external_cartesian_gto`.

Let `S_raw` be the self overlap from the existing Gaussian overlap kernel and
`S_src` the exported source overlap. For each AO, the reader derives only the
positive diagonal factor

```text
d_i = sqrt(S_src[i,i] / S_raw[i,i])
```

and scales that probe's contraction coefficients by `d_i`. It then recomputes
the complete probe overlap and requires

```text
norm(S_probe - S_src, Inf)
    <= overlap_atol + overlap_rtol * max(norm(S_src, Inf), 1)
```

Both diagonal inputs and every factor must be finite and strictly positive.
This is an AO normalization reconciliation, not a metric solve. The existing
MO coefficients remain unchanged because the scaled probes now represent the
actual source AOs. Failure of full-matrix parity is a convention mismatch and
must stop construction.

The reader then constructs `ExternalGTOOrbitalSpinBlock` values and the
existing `ExternalGTOOrbitalPacket`. Existing ordering and `S_GG`
fingerprints remain authoritative for later import. There is no runtime
permutation ledger, label-based reorder, basis lookup, or second overlap
implementation.

## Closest Orthonormal Determinant

The closest operation calls the existing importer on the supplied
same-construction `working` object. For each imported occupied block it forms

```text
G  = C_F' * C_F
C0 = C_F * G^(-1/2)
```

using the repository Lowdin convention

```text
C0 = C_F * inv(sqrt(Symmetric(G)))
```

It separately evaluates the symmetric Gram eigenvalues for diagnostics and
reports all eigenvalues and principal angles

```text
theta_i = acos(sqrt(clamp(lambda_i, 0, 1)))
```

and fails unless `minimum_gram_eigenvalue` is finite and in `(0,1]`, every
eigenvalue is finite and in `[-1e-10,1+1e-10]`, and the minimum eigenvalue is
at least the caller's threshold. The clamp above is diagnostic-only after that
bound check. It also requires `norm(C0' * C0 - I, Inf) <= 1e-10`. It does not
floor, drop, or append a direction.

The direct `imported_coefficients` stored in the import result are not mutated.
Within a fixed `1e-10` representation tolerance, restricted determinant blocks
require occupation two and alpha-only or alpha/beta determinant blocks require
occupation one. Empty, zero, negative, fractional, or mixed restricted
occupations fail this operation even though the generic packet/import
representation may report them.

The v1 closest operation accepts only a same-construction object exposing the
Hamiltonian used with its working basis. It checks final dimension, exact
nuclear charges, nuclear positions within `1e-10` Bohr, the declared
all-electron nuclear-Coulomb kind, and the Hamiltonian electron sector against
the determinant occupations. These are the v1 physical-Hamiltonian identity
checks; no source Hamiltonian matrix is present or inferred. It accepts no
detached `C_F` and no separately supplied Hamiltonian. A custom external
one-body operator or a working basis without this identity boundary remains
eligible only for raw import, not for the public closest-determinant operation.

## Supported And Rejected Cases

Supported in v1:

- explicit real molecular Cartesian contracted GTO records;
- source states exported directly as Cartesian or transformed from PySCF real
  spherical orbitals by `cart2sph_coeff(normalized="sp")`;
- restricted, alpha-only, and collinear alpha/beta packet imports;
- arbitrary source basis names when all AO primitives are explicit;
- raw projected coefficients and separately requested full-rank determinant
  cleanup.

Rejected or out of scope:

- runtime AO permutation tables or structural-key sorting;
- PySCF calculation setup, SCF execution, or basis-name lookup;
- Julia checkpoint, HDF5, NPZ, Molden, or QCSchema readers;
- kinetic, nuclear, Fock, ERI, Hamiltonian, or HFDMRG bundle data;
- GHF, spinor, complex, periodic, ECP, and ROHF-initial-determinant claims;
- generalized final-basis metrics or a second Gaussian overlap path;
- silent rank loss, fractional determinant semantics, Hamiltonian mismatch, or
  automatic acceptance of low captured norm;
- paper request identifiers in any public name, type, field, or error.

## Validation Contract

The committed validation uses one small two-center PySCF H2/cc-pVTZ fixture,
which contains a Cartesian `d` shell. CI reads its frozen same-stem `.toml` and
`.f64` files without Python or PySCF installed. The existing public residual-
GTO owner must reuse its already-built H2 working basis and verify:

- payload and ordered-AO integrity;
- full source-overlap parity after probe scaling;
- source MO metric identity;
- exact reconstruction of the existing `S_FG*C_G` import;
- occupied-subspace capture diagnostics;
- closest-determinant orthonormality and unchanged raw projection;
- malformed units, dimensions, hashes, AO order, source overlap,
  normalization, occupations, geometry, sector, and capture threshold;
- one synthetic multi-column occupied rotation-covariance identity.

The public test may use only root-public interfaces and ordinary Julia/stdlib
operations. It must not call `_cartesian_supplement_cross_overlap` or another
private overlap helper to establish expected results.

After the bounded fixture passes, acceptance requires one read-only replay of
the existing C2/aug-cc-pV6Z paper artifacts. It must reproduce the accepted
occupied subspace using the direct PySCF Cartesian transformation and the new
reader, without copying the large packet, NWChem snapshot, or historical
permutation ledger into the repository. If this direct path disagrees, stop
and report the exact AO order or normalization mismatch.

Optional scheduled PySCF regeneration is a later lifecycle decision. The
initial implementation adds no workflow.

## Implementation Surfaces And Budgets

Approved implementation is limited to:

```text
src/cartesian_external_gto_interchange.jl
src/GaussletBases.jl
bin/export_pyscf_cartesian_gto.py
test/driver_public/cartesian_residual_gto_mwg_system_runtests.jl
test/driver_public/external_cartesian_gto_h2_ccpvtz_v1.toml
test/driver_public/external_cartesian_gto_h2_ccpvtz_v1.f64
```

Limits are stop-and-report bounds, not compression targets:

- Julia source plus root wiring: `340` preferred, `400` hard added lines;
- Python exporter: `220` preferred, `280` hard lines;
- existing public-test extension: `180` preferred, `240` hard added lines;
- frozen two-file fixture: `256 KiB` hard total;
- exactly one new Julia source file, one Python exporter, one logical frozen
  fixture, and one existing public test owner;
- exactly two planned root exports and no planned public type;
- zero `Project.toml`, dependency, workflow, artifact, release, or existing
  packet/import source changes unless a demonstrated blocker returns here.

Do not split helpers or tests into extra files to evade a limit. If a readable,
validated implementation exceeds a hard bound, make no implementation commit
and report the smallest justified adjustment.

## Release Boundary

This public facility is not present in immutable `v0.2.0-rc1`. Including it in
v0.2 requires a separately authorized new release candidate, with its own
clean replay and reader-facing review. This design grants no version change,
tag, GitHub release, registration, citation, or final-release action.

## Failure Rule

Make no implementation commit if the path requires a runtime permutation
ledger, duplicates Gaussian overlap mathematics, adds a mandatory Julia
dependency, weakens full source-overlap parity, silently repairs rank loss, or
cannot reproduce the accepted C2 occupied subspace through PySCF's direct
Cartesian transformation. Report the exact convention mismatch or missing
existing kernel to the repo-design-manager.

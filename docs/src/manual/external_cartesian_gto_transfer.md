# External Cartesian GTO Transfer

GaussletBases can transfer occupied orbitals from an ordinary molecular PySCF
Hartree-Fock checkpoint into a supported orthonormal Cartesian working basis.
The transfer preserves the raw projection and its capture diagnostics. A
separate optional step can prepare the closest orthonormal determinant inside
the captured subspace.

## Export A Checkpoint

First run a real molecular RHF or collinear UHF calculation in PySCF and save
its checkpoint. From a GaussletBases checkout, convert that existing checkpoint
to a same-stem version-1 `.toml`/`.f64` bundle:

```bash
python bin/export_pyscf_cartesian_gto.py \
    --checkpoint state.chk \
    --output state.toml
```

This is a checkpoint-only exporter: it runs no SCF calculation. PySCF and
NumPy are optional dependencies of this external Python command, not Julia
package dependencies.

Running the exporter is also a producer attestation that the checkpoint came
from the ordinary all-electron nuclear-Coulomb Hartree-Fock Hamiltonian. The
checkpoint fields alone cannot prove that DFT or a custom one-body operator was
absent.

## Construct The Working Object

Construct the GaussletBases object with the same nuclei, positions, and
electron sector as the checkpoint. For example, a two-center residual-GTO/MWG
working object can be built entirely through public names:

```julia
using GaussletBases

nuclei = [(0.0, 0.0, -2.0), (0.0, 0.0, 2.0)]
system = (;
    atom_symbols = ["H", "H"],
    nuclear_charges = [1.0, 1.0],
    atom_locations = nuclei,
    nup = 1,
    ndn = 1,
)
basis = (;
    q = 5,
    core_spacing = 0.5,
    xmax_parallel = 6.0,
    xmax_transverse = 4.0,
)
supplement = (;
    basis_by_center = ["cc-pVTZ", "cc-pVTZ"],
    lmax = 1,
    uncontracted = false,
    width_filtering = nothing,
)
working = cartesian_residual_gto_mwg_system(
    system; basis, supplement)
```

The returned working object deliberately has a narrow public surface. Its
Hamiltonian is `working.hamiltonian`; the transfer operations consume the
same-construction overlap data without exposing a general basis carrier.

## Read, Inspect, And Project

Read the bundle, inspect the final-by-source overlap, and retain the unchanged
raw projection:

```julia
packet = read_external_cartesian_gto_packet("state.toml")
S_FG = gto_overlap_matrix(working, packet.probes)
raw = import_external_gto_orbitals(working, packet)

size(S_FG) == raw.cross_overlap_size
raw.alpha.orbital_captures
raw.alpha.density_trace_loss
raw.alpha.worst_orbital_capture
```

For UHF packets, inspect `raw.beta` as well. The import is exactly
`C_F = S_FG * C_G`; it does not repair, renormalize, or rotate the projected
columns. Raw projection can lose norm when the working basis does not contain
the complete occupied space.

## Optional Determinant Preparation

If a downstream solver requires an orthonormal determinant, choose and record
an explicit capture requirement for that calculation:

```julia
minimum_capture = 0.99  # example choice; replace with your justified criterion
start = closest_external_gto_determinant(
    working,
    packet;
    minimum_gram_eigenvalue = minimum_capture,
)

start.alpha.gram_eigenvalues
start.alpha.principal_angles_radians
start.alpha.orthonormality_error
start.imported.alpha.imported_coefficients == raw.alpha.imported_coefficients
```

The minimum Gram eigenvalue is the invariant occupied-subspace capture
criterion. GaussletBases supplies no universal threshold. The cleanup
orthonormalizes only the captured occupied subspace and cannot recover a
missing direction. `start.imported` preserves the unchanged raw projection
beside the cleaned `start.alpha` and optional `start.beta` blocks.

Closest-determinant preparation requires occupations of `2` for each
restricted orbital, or `1` for each alpha and beta orbital. Other positive
occupations may be represented by a packet and inspected through the raw
projection, but they do not define this determinant cleanup.

## Supported Boundary

Version 1 supports real finite-molecule RHF, alpha-only, and collinear UHF
checkpoints. Basis-only export, export from a live mean-field object, GHF,
spinors, complex orbitals, periodic systems, ECP or pseudopotential states,
ROHF determinant semantics, and arbitrary Hamiltonian transfer are unsupported.

The result is a set of transferred orbitals and capture diagnostics. It is not
an HFDMRG state, another solver's restart, or a transformed Hamiltonian.

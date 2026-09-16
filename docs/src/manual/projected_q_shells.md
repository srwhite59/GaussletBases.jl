# Projected q-shells (PQS)

Projected q-shells (PQS) are a Cartesian terminal-basis construction for
localized gausslet Hamiltonians. New users should still begin with the
[radial workflow](../tutorials/first_radial_workflow.md); this page is the
smallest current entrance to the bounded molecular PQS route.

## PQS and White-Lindsey shells

Both methods start from the same mapped three-dimensional parent basis and the
same physical assignment of parent rows to cores, shells, and slabs. They
differ in how each shell's retained functions are formed.

White-Lindsey shells use face, edge, and corner products of one-dimensional
contractions. PQS instead forms a filled source box, selects its boundary
product modes, restricts those modes to the rows owned by the physical shell,
and applies a symmetric Lowdin orthogonalization only within that shell.

The detailed constructions are documented in:

- [Cartesian PQS and IDA overview](../algorithms/cartesian_ida_overview.md)
- [PQS shell construction](../algorithms/pqs_shell_construction.md)
- [IDA Hamiltonian and counterpoise](../algorithms/ida_hamiltonian_and_counterpoise.md)

## Supported public construction

The public base-Hamiltonian facade currently accepts positive-charge,
origin-centered one-center atoms and equal-label/equal-charge homonuclear
diatomics with distinct centers on the Cartesian z axis. It constructs an
unsupplemented, uncorrected localized-IDA Hamiltonian; it does not run a solver.

This PQS-only H2+ construction uses the scientific local order `q = 4`:

```julia
using GaussletBases

system = (;
    atom_symbols = ["H", "H"],
    nuclear_charges = [1.0, 1.0],
    atom_locations = [(0.0, 0.0, -1.0), (0.0, 0.0, 1.0)],
    nup = 1,
    ndn = 0,
)
basis = (;
    q = 4,
    core_spacing = 0.6,
    xmax_parallel = 3.0,
    xmax_transverse = 2.0,
    tail_spacing = 2.8,
    nesting = :pqs,
)
ham = cartesian_base_hamiltonian(system; basis = basis)
h1 = one_body_hamiltonian(ham)
Vee = ham.electron_electron_ida
E_nn = nuclear_repulsion(ham)
```

Here `Vee` is the matrix called ``V_{ee}`` on the algorithm pages. The snippet
only constructs operator data; `h1`, `Vee`, and `E_nn` are not solver results.
For PQS, `q = ns`, so this is the PQS half of Example 39's matched `ns = 4`
fixture with identical physical and extent inputs. White-Lindsey instead uses
`q = ns - 2`.

[`examples/39_pqs_h2plus.jl`](https://github.com/srwhite59/GaussletBases.jl/blob/main/examples/39_pqs_h2plus.jl)
builds a small H2+ case with nuclei at `z = -1,+1` bohr. The PQS and
White-Lindsey calculations share every physical input and differ only in the
`nesting` choice. The example checks the stored geometry and one-electron
sector, the common 293-function dimension, finite symmetric one-body and IDA
interaction matrices, nuclear repulsion, and the lowest-eigenpair residual.

Run it from a checkout after instantiating the project:

```bash
julia --project=. examples/39_pqs_h2plus.jl
```

The printed one-body energies are a bounded construction smoke, not basis
convergence or publication evidence. SCF and correlated solvers, Gaussian
supplements, parent residual functions, screening, and paper-scale campaigns
remain consumer or external workflows.

## Expert finite collinear systems

`cartesian_collinear_working_basis(z, Z; ...)` accepts strictly ordered z-axis
positions and positive nuclear charges, subject to the existing mapping's
validity limits. Electrons remain three-dimensional; boundaries are finite/open.
Without `q`, all construction controls remain explicit; for example:

```julia
expansion = coulomb_gaussian_expansion(doacc=false)
z, Z = [-2.4, 0.0, 2.4], ones(3)
working = cartesian_collinear_working_basis(z, Z;
    core_spacing=0.6, transverse_spacing=0.45,
    padding_parallel=3.0, padding_transverse=6.0, core_side=3,
    angular_reference_count=5, outer_face_count=9,
    tail_spacing=2.8, angular_resolution_scale=1.4, expansion)
operators = cartesian_collinear_operators(working, z, Z; expansion)
```

This small construction returns complete `one_body` and `electron_electron_ida`
matrices plus `nuclear_repulsion`, not a reweightable/artifact Hamiltonian.
The opaque working handle supports existing GTO overlap and raw orbital import.
Converge parent padding/spacing and retention independently. Angular reference
count is not a constant source q; outer count applies to both slab-face axes
and must fit each face. Truncation can introduce transverse asymmetry. These
controls are a tested fixture, not chemical accuracy or long-chain scaling
claims. Dense operators remain quadratic; no solver or periodic extension is
included. The compact sliced-chain approximation remains a separate capability.

For hydrogen, the scientific-q entrance resolves the standard source recipe:

```julia
working = cartesian_collinear_working_basis([-1.8, 0.0, 1.8], ones(3);
    q=5, expansion=coulomb_gaussian_expansion(doacc=true))
```

Integer `q >= 3` sets both spacings to `1.2/(q-1)` bohr, odd core width to
`q` (odd) or `q+1` (even), angular reference to `q`, and angular scale to 1.4.
Shell transverse orders are exactly `(q,q)`; longitudinal order remains the
existing all-nucleus selector's result, including its lower-band-limited
fallbacks. This is not a guarantee that both angular-band criteria succeed.
Tail spacing defaults to 2.8 and both paddings to 10 bohr. Outer-face count
starts at `q` but is a separate convergence control for both in-plane axes;
it must fit the source intervals. Padding, tail spacing and outer count may
be supplied explicitly. Neither these defaults nor successful construction
certifies diffuse convergence.
Supplying **both** spacings selects a fixed-parent comparison, not a standard
scaled-q ladder. Explicit core/reference/scale values must match the recipe.
Scientific q currently requires unit charges; mixed positive charges retain
the fully explicit expert entrance. Mapping rejection is unchanged.

### Gaussian-supplemented collinear systems

Use explicit Cartesian Gaussian candidates with the same working basis:

```julia
packet = read_external_cartesian_gto_packet("state.toml")
system = cartesian_residual_gto_mwg_system(working, z, Z;
    supplement=packet.probes, expansion)
ham = system.hamiltonian
raw = import_external_gto_orbitals(system, packet)
```

The supplied ordered positions and positive charges define the potential and
candidate ownership independently of the construction geometry. Each candidate
must lie exactly on one supplied nucleus; contracted candidates are supported.
The result retains the actual terminal-plus-residual basis for cross overlaps
and raw import. `ham` contains complete `one_body`, `electron_electron_ida` and
`nuclear_repulsion`, not electron-sector, artifact or reweighting operations.
The existing atom/diatomic constructor keeps its `CartesianIDAHamiltonian` return.

Residual selection and orthonormalization use the existing cutoff and merge.
Nuclear contributions are accumulated without per-center final matrices.
Base IDA is unchanged; residual-containing interactions use integral-normalized
moment-matched Gaussians (MWG), not exact four-index Coulomb integrals. Raw import
does not clean up a determinant or repair capture loss. Small primitive tests
do not establish diffuse-space convergence or H-chain accuracy; converge the
supplement and parent/retained spaces for scientific calculations.

## Matched H2+ comparison

[`examples/41_pqs_h2plus_table1.jl`](https://github.com/srwhite59/GaussletBases.jl/blob/main/examples/41_pqs_h2plus_table1.jl)
is the slower fixed public comparison intended for release validation. It
reconstructs one `21 x 21 x 29` parent and matched `1285`-function PQS and
White-Lindsey terminal spaces at the declared Table I parameters. The returned
rows report parent-state capture, electronic and total energies, contraction
error, and error relative to a caller-supplied independent total energy.

Run it explicitly when validating a candidate release:

```bash
julia --project=. examples/41_pqs_h2plus_table1.jl /tmp/pqs_h2plus_table1.tsv
```

The comparison fixes its construction parameters and exposes no scan or staged
producer objects. Its cross-platform acceptance uses declared physical
tolerances rather than raw-byte equality.

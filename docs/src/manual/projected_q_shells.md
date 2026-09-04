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

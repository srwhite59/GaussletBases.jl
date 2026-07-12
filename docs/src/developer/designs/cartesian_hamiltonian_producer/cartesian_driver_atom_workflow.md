# Cartesian Driver Atom Workflow

Status: implemented canonical-driver source contract for explicit
origin-centered one-center atoms under `HP-DRV-ATOM-FN-01` and
`HP-DRV-ATOM-WIRE-01`. The paired test ID is completed evidence with no
continuing permission. Producer-side atom semantics remain owned by
[R1 one-center base atoms](r1_one_center_base_atoms.md).

## Live Selection

The driver selects a base atom with:

```julia
Natom = 1
basisname = nothing
```

There is no live `mode = :base` input. `basisname === nothing` is the base
selection used by the script; a non-`nothing` basis label selects the
separately governed supplemented-atom composition path.

The checked-in defaults describe H2 (`Natom = 2`, `nup = 1`, `ndn = 1`). An
atom run must supply electron counts consistent with its explicit charge. For
example, changing only `Natom` to `1` does not produce a valid neutral H input;
neutral H also requires one total electron.

## Atom Contract

For `Natom = 1`, the driver constructs:

```julia
system = (
    atom_symbols = [String(atom)],
    nuclear_charges = [Float64(Z)],
    atom_locations = [(0.0, 0.0, 0.0)],
    nup = nup,
    ndn = ndn,
)
```

The producer requires finite positive integer-valued `Z`, nonnegative integer
spin-sector counts, and `nup + ndn == round(Int, Z)`. The `atom` value is a
provenance label only; it never supplies charge, electron count, spin, basis,
or ECP behavior. Translated atoms are rejected.

The driver constructs the atom basis from:

```julia
(
    ns = ns,
    core_spacing = core_spacing,
    s_factor = s_factor,
    parent_axis_family = gausslet_family,
    nesting = nesting,
    source_span = source_span,
    radius = padding,
)
```

Thus `padding`, not `ns`, is the physical atom-box radius. `ns` controls the
source/nesting size, and route-local `q` is derived from `nesting`. The driver
does not pass legacy `d`; `core_spacing` is the sole public near-nucleus scale.
It also does not expose atom-specific `reference_spacing`, `tail_spacing`, or
`coulomb_accuracy`, so the producer defaults apply.

## Execution And Artifact

The exported base facade remains:

```julia
cartesian_base_hamiltonian(system; basis, hamfile)
```

The driver does not call that facade as one opaque operation. It invokes the
same non-exported working-basis, product, unit-nuclear, localized-IDA, and
assembly stages directly for timing and due-diligence presentation. Base
assembly writes the ordinary Cartesian IDA artifact to the nonempty `hamfile`.
With `check_file = true`, the driver reads it through
`read_cartesian_ida_hamiltonian` and checks its dimension.

`print_contract = true` prints the explicit atom inputs and the standard
terminal due-diligence report before an energy or residual result is
interpreted. `expected_dimension`, when supplied, must match exactly.

## Source Ownership

The `HP-DRV-ATOM-*` driver contract owns atom input normalization and wiring
only in `bin/cartesian_ham_builder.jl`. Public atom validation and construction
in `src/cartesian_base_hamiltonian.jl` remain under `HP-R1-ATOM-*`, while
supplemented atom composition remains under `HP-COMP-SUPPATOM-*`. This driver
contract owns no package source, test, tool, or committed input fixture.

## Failure Behavior

Unsupported or malformed atom inputs throw, normally with `ArgumentError`.
This includes nonpositive or noninteger charge, nonneutral electron count,
invalid `ns`, spacing, nesting, source span, or `s_factor`, an empty artifact
name, and any translated or multicenter shape presented as an atom. Numerical,
writer, and readback failures propagate.

## Non-Goals

This driver contract does not broaden the public facade, infer element data,
accept translated atoms, add ECP/pseudopotential or solver/RHF behavior, alter
Residual Gaussian or IDA/MWG conventions, create an atom-only Hamiltonian
builder, add artifact fields, or expose route diagnostics. Supplemented atoms
reuse the shared composition contract and are not redefined here.

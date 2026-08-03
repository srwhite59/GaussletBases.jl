# R1 Public Base Producer

Status: implemented exported public facade for unsupplemented, uncorrected,
all-electron Cartesian IDA Hamiltonians. This page is the canonical facade,
input, composition, and base-producer provenance contract. Explicit charged
electron sectors are approved under `HP-R1-ESECTOR-*` and remain pending
source implementation.

## Owned IDs

- `HP-R1-FILE-01` - implemented public producer source owner;
- `HP-R1-FN-01` - implemented exported facade;
- `HP-R1-CORE-FN-01` - implemented public near-nucleus spacing convention;
- `HP-R1-WIRE-01` - implemented report-free staged composition;
- `HP-R1-ART-01` - implemented base-producer provenance sidecar;
- `HP-R1-TEST-01` - implemented standalone public endpoint gate.
- `HP-R1-ESECTOR-FN-01` - approved explicit charged-sector relaxation;
- `HP-R1-ESECTOR-TEST-01` - approved sector-independence validation.

Exact lifecycle, permission, source, test, and dependency metadata is recorded
in [the registry](registry.md).

## Meaning Of Base

`base` means the Hamiltonian before Gaussian supplements, screened-Hartree or
other corrections, fragment/counterpoise changes, and solver processing. The
output uses the terminal Cartesian basis, exact assembled one-body matrices,
and the producer's localized IDA electron-electron matrix.

Supplements, corrected Hamiltonians, ECPs, arbitrary molecular geometry, and
solver workflows are separate contracts. They must not be added by redefining
this facade.

## Public Interface

The function is exported from `GaussletBases`:

```julia
cartesian_base_hamiltonian(
    system::NamedTuple;
    basis::NamedTuple,
    hamfile::Union{Nothing,AbstractString} = nothing,
)::CartesianIDAHamiltonian{Float64}
```

Implementation owner:

```text
src/cartesian_base_hamiltonian.jl
```

Export/include owner:

```text
src/GaussletBases.jl
```

The return is the existing `CartesianIDAHamiltonian{Float64}` directly. The
facade never returns a wrapper, status object, report, payload, or partial
result.

## System Input

`system` must contain exactly:

- `atom_symbols::AbstractVector`;
- `nuclear_charges::AbstractVector`;
- `atom_locations::AbstractVector`;
- `nup`;
- `ndn`.

Center collections must have equal lengths. Each location is a fixed
three-element tuple with finite coordinates. Nuclear charges are converted to
`Float64` and must be finite and positive. `nup` and `ndn` must be nonnegative
integers and may not be `Bool`. Unknown or missing keys fail.

### One-Center Atom

The approved atom scope is:

- exactly one center at `(0.0, 0.0, 0.0)`;
- finite positive nuclear charge;
- explicit nonnegative `nup` and `ndn` valid for the realized orbital
  dimension, with positive total electron count under the current
  compatibility container;
- atom symbol used as an explicit provenance label, not as charge, electron,
  basis, or ECP authority.

The broader atom rationale and exclusions are canonical in
[R1 one-center base atoms](r1_one_center_base_atoms.md).

### Homonuclear Z-Axis Diatomic

The approved molecular scope is:

- exactly two centers;
- equal atom-symbol labels and equal finite positive charges;
- both centers on the Cartesian z axis with distinct finite z coordinates;
- explicit nonnegative `nup` and `ndn` valid for the realized orbital
  dimension, with positive total electron count under the current
  compatibility container.

There is no translation or rotation step. Heteronuclear, x/y-aligned,
shifted-parallel, generally oriented, or ECP systems fail validation.

## Electron-Sector Independence

For fixed nuclear charges, nuclear positions, basis input, and numerical
policy, the bare Cartesian construction is independent of `nup` and `ndn`.
The terminal basis, kinetic matrix, unit-nuclear matrices, assembled physical
one-body matrix, IDA interaction, and nuclear repulsion must therefore be
identical when only the electron sector changes. The net molecular charge is
derived rather than independently supplied:

```text
Q_net = sum(nuclear_charges) - (nup + ndn).
```

The compatibility facade continues to require explicit `nup` and `ndn` and
continues to return `CartesianIDAHamiltonian`, which stores those counts as
sector metadata. This pass does not add an optional neutral default,
`net_charge` input, derived-charge field, public accessor, artifact key, or
cache-key convention. Existing neutral calls remain unchanged.

The same invariant composes through the already-implemented supplemented
atom/diatomic path: changing only the sector must not change supplement
loading, residual selection, augmented one-body matrices, or MWG interaction.
This authorizes no supplement algorithm change.

A future architecture may separate a particle-number-independent Cartesian
operator from an electron-sector problem object. That is accepted design
direction only: `cartesian_operator`, `with_electron_sector`, replacement
types, artifact migration, and generalized correction/cache machinery require
separate authority. Sector- or reference-dependent corrections must remain
explicit objects and must not alter the meaning of the bare producer.

## Basis Input

`basis` is a plain `NamedTuple`. It must contain `core_spacing` and at least one
of `ns` or legacy-compatible `q`. Size values must be positive integers and
may not be `Bool`.

Geometry-specific required fields:

| System | Required extent fields |
| --- | --- |
| one-center atom | `radius` |
| z-axis diatomic | `xmax_parallel`, `xmax_transverse` |

All spacings and extents must be finite and positive.

Implemented optional fields and defaults:

| Field | Default / rule |
| --- | --- |
| `parent_axis_family` | `:G10`; no other family is accepted |
| `reference_spacing` | `1.0` |
| `tail_spacing` | `10.0` |
| `nesting` | `:pqs`; accepted values are `:pqs` and `:wl` |
| `source_span` | `:ordinary`; `:mapped_comx` is PQS-only |
| `s_factor` | `1.0`, finite and positive |
| `coulomb_accuracy` | `:compact`; current source also accepts `:high` |

`nesting`, `source_span`, and `coulomb_accuracy` may be supplied as symbols or
strings and are normalized to symbols before construction.

Public `ns` and route-local `q` semantics belong to
[nesting and supplement composition](nesting_supplement_composition_plan.md).
Durable calls use `ns`; `q` remains a compatibility input. When both are
present, `q` must match the value derived from `ns` and `nesting`.

`core_spacing` is the one public near-nucleus physical scale. It is not
`reference_spacing`, box extent, or a hidden element default. One-center
mapping uses the resolved `core_spacing` as its internal mapping `d`.

Public `d` is deprecated compatibility input for one-center atoms only. If
present, it must be finite, positive, and exactly equal to `core_spacing`.
Diatomics reject `d`; public `parent_mapping_d` is unsupported. The expert
mapping-strength exception is separately canonical in
[PQS mapping `s_factor`](pqs_mapping_s_factor.md).

The complete three-tier `:compact | :standard | :high` design belongs to
[Coulomb accuracy policy](coulomb_accuracy_policy.md). At the Pass 378 source
baseline, this facade still validates only `:compact | :high`; `:standard` is
approved there but has not yet landed in `src/cartesian_base_hamiltonian.jl`.
This R1 reconciliation does not change that implementation lifecycle.

Unknown basis keys, invalid types, incompatible geometry-specific fields, and
unsupported policy combinations fail rather than being ignored.

## Report-Free Composition

The facade uses the current staged producer, without route reports or old
pair/assembly materialization wrappers:

```text
cartesian_base_working_basis
-> cartesian_base_products
-> cartesian_base_unit_nuclear
-> cartesian_base_vee
-> cartesian_base_hamiltonian_assembly
-> optional base artifact write
```

One resolved input and Coulomb expansion are carried through the parent,
terminal basis, one-body, unit-nuclear, and IDA interaction construction. The
assembly returns the existing Hamiltonian with explicit charges, positions,
and electron counts.

Numerical terminal realization, blockwise one-body construction, localized
IDA, and Hamiltonian assembly are canonical in
[Terminal basis and base assembly](terminal_basis_and_base_assembly.md). This
facade does not duplicate those algorithms or expose their internal objects.

## Artifact Behavior

If `hamfile === nothing`, no file is written. A nonempty path writes the
existing version-1 `:cartesian_ida_hamiltonian` artifact and the facade still
returns the same in-memory Hamiltonian. Empty paths fail; file-system errors
propagate normally. Production does not require readback.

The minimal matrix payload and `read_cartesian_ida_hamiltonian(...)` are owned
by `src/cartesian_ida_hamiltonian.jl`. The reader reconstructs the Hamiltonian
from its standard matrix and physical metadata keys and ignores additional
provenance groups.

### `producer_provenance/`

`HP-R1-ART-01` owns one fixed base-producer sidecar in the same file. Current
keys are:

```text
provenance_version        producer
nesting                   route
ns                        q
q_rule                    ns_source
core_spacing              s_factor
reference_spacing         tail_spacing
parent_axis_family        parent_axis_counts
mapping_kind              mapping_d
mapping_s_factor          mapping_s_standard
mapping_s_effective       radius
xmax_parallel             xmax_transverse
atom_symbols              nuclear_charges
atom_locations            nup
ndn                       final_dimension
```

`producer` is `:cartesian_base_hamiltonian`. `route` is derived truthfully
from system kind and nesting. `ns`, `q`, `q_rule`, and `ns_source` preserve the
public-size normalization. Atom artifacts record
`:white_lindsey_atomic_mapping`, resolved `mapping_d`, and `radius`; diatomic
artifacts record `mapping_d = nothing` and their explicit extents.

Mapping-strength provenance records the requested factor, standard value, and
effective value. The sidecar is consumer provenance only; construction stages
must not read it back as numerical input.

The broader `hamiltonian_manifest/`, `recipe_provenance/`, and
`coulomb_expansion/` groups are governed by
[Cartesian Hamiltonian artifact manifest](cartesian_hamiltonian_artifact_manifest.md)
and the Coulomb policy. They are not duplicated R1 schemas.

## Library Facade Versus Driver

The public library facade is the function documented above. The human-facing
canonical driver is the separate script:

```text
bin/cartesian_ham_builder.jl
```

The driver owns editable defaults, trusted input-file/override handling,
terminal due-diligence presentation, coarse stage timing, artifact checks, and
optional supplemented workflow. It currently calls the same staged producer
functions directly so timings remain visible. Driver variables are not facade
keywords, and driver defaults are not hidden library defaults.

Driver authority is canonical in
[Cartesian driver usability workflow](cartesian_driver_usability_workflow.md).

## Failure Behavior

Malformed public requests throw, normally with `ArgumentError`, before
expensive construction where practical. This includes:

- missing or unknown keys;
- non-vector center collections or non-tuple locations;
- mismatched center counts;
- nonfinite/nonpositive charges, spacings, or extents;
- invalid `ns`, `q`, electron counts, geometry, nesting, source span, mapping,
  or Coulomb policy;
- translated atoms, unsupported diatomics, mismatched atom `d`, any diatomic
  `d`, or empty `hamfile`.

Numerical construction and file-system failures propagate. The facade never
converts failure into `nothing`, readiness flags, blocker symbols, or partial
artifacts.

## Validation

The committed standalone gate is:

```text
test/driver_public/cartesian_base_hamiltonian_runtests.jl
```

It validates the exported facade for one-center H and z-axis H2, direct
Hamiltonian return, finite/symmetric matrices, endpoint values, compact/high
Coulomb behavior, atom and molecular geometry rejection, key/type/spacing
validation, deprecated `d` behavior, artifact write/readback, and base
producer provenance. The current H2 fixture has final dimension `487`.

`HP-R1-ESECTOR-TEST-01` additionally approves exact neutral/charged He and
He2 operator parity in this existing file, supplemented He2 parity in
`test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`, and charged
artifact/readback sector preservation. It adds no new test or fixture file.

This is a focused public integration gate, not default tiny-unit coverage. It
must not assert private route-stage inventories or teach the human-facing
driver as the facade implementation.

## Non-Goals

This contract does not authorize:

- supplement-algorithm changes, numerical-complete/protected bases, or
  corrections;
- solver, RHF/UHF, fragment, counterpoise, ECP, or Cr2 workflow;
- translated atoms, arbitrary orientation, or heteronuclear molecules;
- new public symbols, input objects, result wrappers, reports, or status
  payloads;
- new artifact kinds, matrix keys, public provenance readers, or schema
  expansion;
- mapping, nesting, Coulomb, terminal-basis, or driver policy changes.

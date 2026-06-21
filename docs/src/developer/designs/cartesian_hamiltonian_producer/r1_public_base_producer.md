# R1 Public Base Producer Candidate

Status: candidate design amendment, not implementation authority.

This document defines the proposed minimal public base Cartesian Hamiltonian
producer for first H/H2 use. It does not approve source work, tests, tools,
driver edits, new report fields, new artifact formats, or implementation of
the candidate IDs below.

## Candidate IDs

- `HP-R1-FILE-01` - public base producer source file.
- `HP-R1-FN-01` - public base Hamiltonian producer facade.
- `HP-R1-WIRE-01` - report-free base producer wiring from the facade to the
  existing terminal-basis and Hamiltonian construction path.
- `HP-R1-TEST-01` - small committed public endpoint test/example for the
  facade.

Both IDs are candidate-only until explicitly approved in `registry.md`.

## Ownership And Export

Candidate owner: top-level `GaussletBases` public API.

Candidate source file:

```text
src/cartesian_base_hamiltonian.jl
```

If approved, `cartesian_base_hamiltonian` should be exported from
`src/GaussletBases.jl`. No other public symbol is proposed by this candidate.

The source file may use existing internal driver helpers, but it must not make
those helpers public API and must not expose route-stage vocabulary in the
public call shape.

## Public Call Shape

The proposed public entry point is a single function:

```julia
cartesian_base_hamiltonian(
    system::NamedTuple;
    basis::NamedTuple,
    method::Symbol = :pqs_source_box,
    route::Symbol = :auto,
    output::NamedTuple = (;),
)::CartesianIDAHamiltonian{Float64}
```

The function returns the existing `CartesianIDAHamiltonian{Float64}` directly.
It must not return a wrapper, status object, materialization payload, report
mirror, or `(value, status)` pair.

The first implementation scope is base H and bond-aligned base H2. Broader
atoms, general molecules, WL/QW unification, supplements, corrections, solver
handoff, and Cr2-scale performance remain later roadmap lanes unless separately
approved.

## Input Shape

All public inputs are plain `NamedTuple` groups. This candidate does not
introduce public input structs, builder objects, config payloads, or keyword
field clouds attached to intermediate stages.

`system` is a `NamedTuple`.

Required fields for first H/H2 scope:

- `atom_symbols`
- `nuclear_charges`
- `atom_locations`
- `nup`
- `ndn`
- `radius`

Optional fields:

- `bond_axis`
- `bond_length`
- `map_backend`

`basis` is a `NamedTuple`.

Required fields for first H/H2 scope:

- `q`
- `n_s`
- `core_spacing`
- `xmax_parallel`
- `xmax_transverse`
- `parent_axis_family`

Optional fields with R1 defaults:

- `reference_spacing = 1.0`
- `tail_spacing = 10.0`
- `q_to_core_spacing_rule = :standard_pqs_ns_equals_q`
- `parent_axis_bundle_backend = :pgdg_localized_experimental`
- `parent_mapping_rule = :identity_mapping`
- `parent_mapping_Z = nothing`
- `parent_mapping_d = nothing`
- `parent_axis_counts = nothing`

`method` / `route` choose the construction method:

- first candidate method: `:pqs_source_box`
- first candidate route: automatic one-center atomic or bond-aligned diatomic
  base PQS route selection from `system`
- allowed route overrides, if needed: `:one_center_atomic_base` and
  `:bond_aligned_diatomic_base`
- private route-kind symbols are not part of the public API

`output` is a `NamedTuple`.

Allowed fields:

- `save_ham_artifact::Bool = false`
- `hamfile`

No other output fields are part of the minimal base R1 contract:

- no basis artifact is part of this base R1 contract
- no TSV/report artifact is part of this base R1 contract

Grouped inputs are call-boundary convenience only. They must not become
persistent payload structs, status summaries, metadata field clouds, or report
mirrors.

## First Public Examples

One-center H:

```julia
h_system = (;
    atom_symbols = ("H",),
    nuclear_charges = (1,),
    atom_locations = ((0.0, 0.0, 0.0),),
    nup = 1,
    ndn = 0,
    radius = 4.0,
)

h_basis = (;
    q = 5,
    n_s = 5,
    parent_axis_family = :G10,
    core_spacing = 0.5,
    xmax_parallel = 4.0,
    xmax_transverse = 4.0,
)

h_ham = cartesian_base_hamiltonian(h_system; basis = h_basis)
```

Bond-aligned H2:

```julia
h2_system = (;
    atom_symbols = ("H", "H"),
    nuclear_charges = (1, 1),
    atom_locations = ((0.0, 0.0, -2.0), (0.0, 0.0, 2.0)),
    nup = 1,
    ndn = 1,
    bond_axis = :z,
    bond_length = 4.0,
    radius = 4.0,
)

h2_basis = (;
    q = 5,
    n_s = 5,
    parent_axis_family = :G10,
    core_spacing = 0.5,
    xmax_parallel = 6.0,
    xmax_transverse = 4.0,
)

h2_ham = cartesian_base_hamiltonian(
    h2_system;
    basis = h2_basis,
    output = (;
        save_ham_artifact = true,
        hamfile = "h2_cartesian_ida_hamiltonian.jld2",
    ),
)
```

The examples intentionally do not call `cartesian_pair_terms`,
`cartesian_assembly`, `cartesian_report`, or route-internal fixture controls.

## Output And Artifact Contract

The output is the existing `CartesianIDAHamiltonian{Float64}`.

If `save_ham_artifact = true`, the implementation must use the existing
`write_cartesian_ida_hamiltonian` writer and existing readback path. If
`save_ham_artifact = true` and `hamfile` is missing, `nothing`, or an empty
string, the facade must throw `ArgumentError` before writing. It must not invent
a default file name for the public API. Parent-directory and file-system errors
from the existing writer should propagate normally.

This candidate does not approve:

- a new artifact file shape;
- a new artifact manifest;
- a new basis/provenance artifact;
- a new materialization wrapper;
- durable status/result-kind fields.

## Error Behavior

Unsupported or malformed public requests must throw clear exceptions, normally
`ArgumentError` for input validation failures:

- missing required `system`, `basis`, or `output` fields;
- unsupported `method`;
- unsupported `route`;
- unsupported atom count, route geometry, or non-H/H2 system in the first R1
  scope;
- inconsistent field lengths for symbols, charges, positions, or electron
  counts;
- requested artifact write without `hamfile`.

The facade must not return status objects, readiness summaries, blocker
symbols, `nothing`, or partial payloads for unsupported public requests.
Internal numerical construction errors may propagate after validation.

## Public Workflow

The recommended public workflow is:

```text
system / specification
-> parent and route geometry
-> terminal basis realization
-> Hamiltonian production
-> optional existing Hamiltonian artifact
```

`cartesian_pair_terms` and `cartesian_assembly` are absent from the recommended
public base workflow. They may remain temporarily for legacy script and report
compatibility, but R1 public documentation and examples must not teach users to
call them for base Hamiltonian construction.

## Mapping From Current Private Driver

The current private driver and harness group inputs like this:

| Current private group | R1 public input |
| --- | --- |
| `system_inputs` | `system` |
| `spacing_inputs` and `parent_inputs` | `basis` |
| `route_inputs` | `method` / `route` |
| `materialization_inputs` | `output` |
| `save_inputs`, TSV/report options | not part of the minimal base producer |

The current private execution spine is:

```text
cartesian_system
-> cartesian_recipe
-> cartesian_parent
-> cartesian_shells
-> cartesian_units
-> cartesian_transforms
-> cartesian_pair_terms
-> cartesian_assembly
-> cartesian_report
-> cartesian_materialization
```

The R1 facade should instead map to:

```text
cartesian_system / recipe / parent / shells / units
-> cartesian_transforms
-> terminal_basis_realization
-> base Hamiltonian materialization
-> optional write_cartesian_ida_hamiltonian
```

The direct dependency to remove is report recovery through
`cartesian_assembly`: today `cartesian_report` recovers `route_skeleton` and a
low-order shellification summary from `cartesian_assembly`. R1 must not add a
new base consumer of `cartesian_pair_terms` or `cartesian_assembly` just to feed
that report dependency. The implementation must either use a report-free base
materialization boundary or narrow existing report construction without adding
new report field clouds.

## Validation Gates For First Implementation

First implementation must validate:

- package load;
- H base Hamiltonian construction through the public facade with a reviewed
  one-center H system;
- H2 base Hamiltonian construction through the public facade;
- returned type is `CartesianIDAHamiltonian{Float64}`;
- H2 final dimension remains `471`;
- H2 `one_body_hamiltonian(ham)` lowest value remains
  `-0.79460371733658908` within reviewed tolerance;
- H2 localized IDA self-Coulomb remains `0.4569117646737212` within reviewed
  tolerance;
- `K`, every unit `U_A`, and `V` are finite and symmetric within reviewed
  tolerance;
- invalid public requests throw clear `ArgumentError`s rather than returning
  status/blocker objects;
- `save_ham_artifact = true` without `hamfile` throws `ArgumentError`;
- optional artifact write/readback uses `write_cartesian_ida_hamiltonian` and
  has zero one-body readback delta;
- R0 warm/cold baseline is not materially regressed without explanation;
- no `cartesian_pair_terms` or `cartesian_assembly` call appears in the
  recommended public example path.

Validation policy decision: ignored validation is not enough for a new public
facade. R1 should include one small committed public endpoint test or executable
example after `HP-R1-TEST-01` is explicitly approved. The test/example should
exercise only the public facade, H/H2 endpoint facts, and existing
Hamiltonian-artifact readback. It must not assert route-internal stage fields,
status symbols, report mirrors, terminal role vocabulary, or pair inventories.

## Deletion And Shrinkage Targets

The R1 implementation should shrink or quarantine:

- public documentation that teaches pair/assembly stages as required base
  Hamiltonian workflow;
- use of `cartesian_pair_terms` and `cartesian_assembly` in the recommended
  public base path;
- report-only dependency on `route_skeleton` and low-order shellification
  summary when it exists only to support base materialization;
- driver/harness-only materialization options that are not part of the public
  base producer contract.

Do not delete useful local product-box, 1D factor, terminal-block, or
term-first contraction kernels merely because their historical file or helper
name contains "pair". Those kernels remain donor/oracle inventory for R2/R3
classification.

## Forbidden Additions

The R1 base producer candidate forbids:

- new payload/result/status wrapper around `CartesianIDAHamiltonian`;
- new status, blocker, readiness, materialized, or result-kind field clouds;
- new report fields duplicating terminal bases, matrices, factors, or raw pair
  tensors;
- metadata-carried numerical data;
- new artifact shapes;
- public solver controls;
- supplement, correction, fragment, counterpoise, or Cr2 stress functionality;
- broad compatibility adapters that preserve the old pair/assembly story.

If any of these become necessary, stop and write a separate docs-only design
amendment before source work.

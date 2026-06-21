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

All four R1 IDs are candidate-only until explicitly approved in `registry.md`.

## Ownership And Export

Candidate owner: top-level `GaussletBases` public API.

Candidate source file:

```text
src/cartesian_base_hamiltonian.jl
```

If approved, `cartesian_base_hamiltonian` should be exported from
`src/GaussletBases.jl`. No other public symbol is proposed by this candidate.

The name "base" is permanent terminology in this candidate: it means
unsupplemented, uncorrected, all-electron localized-IDA Hamiltonian
construction. Supplements, corrections, fragments, counterpoise, solver
handoff, and non-base Hamiltonians are separate roadmap lanes. If a later API
absorbs those directly, it should be a broader API such as
`cartesian_hamiltonian`, not a redefinition of this base facade.

The source file may use existing internal driver helpers, but it must not make
those helpers public API and must not expose route-stage vocabulary in the
public call shape.

## Public Call Shape

The proposed public entry point is a single function:

```julia
cartesian_base_hamiltonian(
    system::NamedTuple;
    basis::NamedTuple,
    hamfile::Union{Nothing,AbstractString} = nothing,
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

Unknown keys in public input groups must throw `ArgumentError`. The facade must
normalize accepted public inputs immediately into fixed internal data and must
not carry arbitrary public `NamedTuple` shapes through the staged
implementation.

`system` is a `NamedTuple`. Center-sized collections must be vectors or other
`AbstractVector` values, not variable-size tuples.

Required fields for first H/H2 scope:

- `atom_symbols::AbstractVector`
- `nuclear_charges::AbstractVector{<:Real}`
- `atom_locations::AbstractVector{<:NTuple{3}}`
- `nup`
- `ndn`

No optional `system` fields are part of the R1 public contract. Bond axis and
bond length are derived from `atom_locations`. `map_backend` is a construction
choice, not physical-system data, and is not public in R1.

`basis` is a `NamedTuple`.

Required common fields:

- `q`
- `core_spacing`

Conditional fields:

- one-center H requires `radius`
- bond-aligned H2 requires:
  - `xmax_parallel`
  - `xmax_transverse`

R1 derives `n_s = q` internally.

Optional public fields with R1 defaults:

- `parent_axis_family = :G10`
- `reference_spacing = 1.0`
- `tail_spacing = 10.0`

Fixed private choices in R1:

- `method = :pqs_source_box`
- `q_to_core_spacing_rule = :standard_pqs_ns_equals_q`
- `parent_axis_bundle_backend = :pgdg_localized_experimental`
- `parent_mapping_rule = :identity_mapping`

These are not public keywords or accepted `basis` fields in R1:

- `n_s`
- `bond_axis`
- `bond_length`
- `map_backend`
- `parent_axis_counts`
- `parent_axis_bundle_backend`
- `parent_mapping_rule`
- `parent_mapping_Z`
- `parent_mapping_d`

`xmax_parallel` and `xmax_transverse` are accepted only for bond-aligned H2.
They are unknown keys for one-center H.

Routing is implicit. One center selects the one-center base route. Two centers
must be bond-aligned and select the bond-aligned diatomic base route. R1 has no
public `method`, `route`, or output group selector. New method or route
keywords should wait until R2 or later provides a second real public method.

No basis artifact, TSV, report artifact, or output group is part of this base
R1 contract.

## First Public Examples

One-center H:

```julia
h_system = (;
    atom_symbols = ["H"],
    nuclear_charges = [1.0],
    atom_locations = [(0.0, 0.0, 0.0)],
    nup = 1,
    ndn = 0,
)

h_basis = (;
    q = 5,
    core_spacing = 0.5,
    radius = 4.0,
)

h_ham = cartesian_base_hamiltonian(h_system; basis = h_basis)
```

Bond-aligned H2:

```julia
h2_system = (;
    atom_symbols = ["H", "H"],
    nuclear_charges = [1.0, 1.0],
    atom_locations = [(0.0, 0.0, -2.0), (0.0, 0.0, 2.0)],
    nup = 1,
    ndn = 1,
)

h2_basis = (;
    q = 5,
    core_spacing = 0.5,
    xmax_parallel = 6.0,
    xmax_transverse = 4.0,
)

h2_ham = cartesian_base_hamiltonian(
    h2_system;
    basis = h2_basis,
    hamfile = "h2_cartesian_ida_hamiltonian.jld2",
)
```

The examples intentionally do not call `cartesian_pair_terms`,
`cartesian_assembly`, `cartesian_report`, or route-internal fixture controls.

## Output And Artifact Contract

The output is the existing `CartesianIDAHamiltonian{Float64}`.

If `hamfile !== nothing`, the implementation must use
`write_cartesian_ida_hamiltonian(hamfile, ham)` after constructing the
Hamiltonian and still return `ham`. If `hamfile` is `nothing`, no artifact is
written. If `hamfile` is an empty string, the facade must throw `ArgumentError`
before writing. Parent-directory and file-system errors from the existing
writer should propagate normally.

The production facade must not automatically read back the artifact after
writing. `read_cartesian_ida_hamiltonian` is for validation and users.

This candidate does not approve:

- a new artifact file shape;
- a new artifact manifest;
- a new basis/provenance artifact;
- a new materialization wrapper;
- durable status/result-kind fields.

## Error Behavior

Unsupported or malformed public requests must throw clear exceptions, normally
`ArgumentError` for input validation failures:

- missing required `system` or `basis` fields;
- unknown `system` or `basis` fields;
- unsupported atom count, route geometry, or non-H/H2 system in the first R1
  scope;
- non-bond-aligned two-center geometry;
- inconsistent field lengths for symbols, charges, positions, or electron
  counts;
- empty `hamfile`;
- center-sized tuple inventories for atom symbols, charges, or locations.

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
| `spacing_inputs` and `parent_inputs` | `basis` subset plus fixed R1 defaults/private choices |
| `route_inputs` | fixed internally to base PQS with geometry-inferred route |
| `materialization_inputs` | `hamfile` keyword only |
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
-> shared base Hamiltonian constructor
-> optional write_cartesian_ida_hamiltonian
```

The direct dependency to remove is report recovery through
`cartesian_assembly`: today `cartesian_report` recovers `route_skeleton` and a
low-order shellification summary from `cartesian_assembly`. R1 must not add a
new base consumer of `cartesian_pair_terms` or `cartesian_assembly` just to feed
that report dependency. The implementation must either use a report-free base
materialization boundary or narrow existing report construction without adding
new report field clouds.

## Report-Free Shared Constructor Seam

`HP-R1-WIRE-01` candidate wiring owns one private shared constructor seam:

```julia
_cartesian_base_ida_hamiltonian(
    terminal_basis_realization,
    parent_axis_bundle_object,
    atom_locations::Vector{NTuple{3,Float64}},
    nuclear_charges::Vector{Float64},
    nup::Int,
    ndn::Int,
)::CartesianIDAHamiltonian{Float64}
```

Candidate owner file:

```text
src/pqs_source_box_low_order_materialization.jl
```

Allowed caller files after approval:

- `src/pqs_source_box_low_order_materialization.jl`
- `src/cartesian_base_hamiltonian.jl`

The existing report-bound materialization path and the new public facade should
call this same private constructor. The implementation should move the current
K/unit-`U_A`/localized-IDA/`CartesianIDAHamiltonian` orchestration behind this
seam and delete the duplicated report-bound numerical orchestration it
replaces. It must not fabricate a private report object, duplicate the
Hamiltonian builder in the new public file, or leave two parallel base
Hamiltonian construction paths.

## Validation Gates For First Implementation

First implementation must validate:

- package load;
- H base Hamiltonian construction through the public facade with a reviewed
  one-center H system;
- H2 base Hamiltonian construction through the public facade;
- returned type is `CartesianIDAHamiltonian{Float64}`;
- H one-body baseline remains `-0.49855234726272035` within `1.0e-10`;
- H2 final dimension remains `471`;
- H2 `one_body_hamiltonian(ham)` lowest value remains
  `-0.79460371733658908` within reviewed tolerance;
- H2 localized IDA self-Coulomb remains `0.4569117646737212` within reviewed
  tolerance;
- `K`, every unit `U_A`, and `V` are finite and symmetric within reviewed
  tolerance;
- invalid public requests throw clear `ArgumentError`s rather than returning
  status/blocker objects;
- unknown public input keys throw `ArgumentError`;
- empty `hamfile` throws `ArgumentError`;
- non-`nothing` `hamfile` writes with `write_cartesian_ida_hamiltonian`;
- validation readback with `read_cartesian_ida_hamiltonian` has zero one-body
  readback delta;
- R0 warm/cold baseline is not materially regressed without explanation;
- no `cartesian_pair_terms` or `cartesian_assembly` call appears in the
  recommended public example path.

Validation policy decision: ignored validation is not enough for a new public
facade. R1 should include one small committed public endpoint test or executable
example after `HP-R1-TEST-01` is explicitly approved. The test/example should
exercise only the public facade, H/H2 endpoint facts, and existing
Hamiltonian-artifact readback. It must not assert route-internal stage fields,
status symbols, report mirrors, terminal role vocabulary, or pair inventories.

Candidate test path:

```text
test/driver_public/cartesian_base_hamiltonian_runtests.jl
```

Candidate cadence: integration/endpoint gate, not default tiny unit coverage.
Run it when touching the public facade, base producer wiring, Hamiltonian
artifact writer/reader, or terminal base-Hamiltonian path. Do not add it to the
ordinary fast per-edit test set without explicit repo-manager approval.

Test artifact behavior: use `mktempdir()` for `hamfile` output.

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

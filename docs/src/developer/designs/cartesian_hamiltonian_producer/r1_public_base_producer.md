# R1 Public Base Producer

Status: approved design amendment for first R1 implementation.

This document defines the minimal public base Cartesian Hamiltonian producer
for first H/H2 use. It approves only the R1 IDs listed below and only within
the scope recorded in `registry.md`. It does not approve tools, driver edits,
new report fields, new artifact formats, broad public-driver polish, or
implementation outside these IDs.

## Approved IDs

- `HP-R1-FILE-01` - public base producer source file.
- `HP-R1-FN-01` - public base Hamiltonian producer facade.
- `HP-R1-WIRE-01` - report-free base producer wiring from the facade to the
  existing terminal-basis and Hamiltonian construction path.
- `HP-R1-TEST-01` - small committed public endpoint test/example for the
  facade.

All four R1 IDs are approved for first implementation in `registry.md`.

## Ownership And Export

Approved owner: top-level `GaussletBases` public API.

Approved source file:

```text
src/cartesian_base_hamiltonian.jl
```

`cartesian_base_hamiltonian` should be exported from `src/GaussletBases.jl`.
No other public symbol is approved by this R1 design.

The name "base" is permanent terminology in this design: it means
unsupplemented, uncorrected, all-electron localized-IDA Hamiltonian
construction. Supplements, corrections, fragments, counterpoise, solver
handoff, and non-base Hamiltonians are separate roadmap lanes. If a later API
absorbs those directly, it should be a broader API such as
`cartesian_hamiltonian`, not a redefinition of this base facade.

The source file may use existing internal driver helpers, but it must not make
those helpers public API and must not expose route-stage vocabulary in the
public call shape.

## Public Call Shape

The approved public entry point is a single function:

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

The first implementation scope is origin-centered base H and Cartesian z-axis
aligned base H2. Broader atoms, x/y-aligned diatomics, generally oriented
molecules, WL/QW unification, supplements, corrections, solver handoff, and
Cr2-scale performance remain later roadmap lanes unless separately approved.

## Input Shape

All public inputs are plain `NamedTuple` groups. This design does not
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

No optional `system` fields are part of the R1 public contract. Bond length is
derived from `atom_locations`; the bond axis must be Cartesian `z` in R1.
`map_backend` is a construction choice, not physical-system data, and is not
public in R1.

Scalar and collection validation rules:

- `q` must be a positive integer;
- `core_spacing`, `reference_spacing`, `tail_spacing`, `radius`,
  `xmax_parallel`, and `xmax_transverse` must be finite and positive when
  present;
- all coordinates and nuclear charges must be finite;
- first scope supports only H and H2 with symbols and charges consistent with
  hydrogen nuclei;
- `nup` and `ndn` must be nonnegative integers and must give the supported
  total electron count for the requested H or H2 system;
- one-center H must be at `(0.0, 0.0, 0.0)`;
- z-axis H2 centers must have finite `x == 0`, finite `y == 0`, and distinct
  finite `z` coordinates.

`basis` is a `NamedTuple`.

Required common fields:

- `q`
- `core_spacing`

Conditional fields:

- one-center H requires `radius`
- z-axis H2 requires:
  - `xmax_parallel`
  - `xmax_transverse`

R1 derives `n_s = q` internally.

Optional public fields with R1 defaults:

- `parent_axis_family = :G10`
- `reference_spacing = 1.0`
- `tail_spacing = 10.0`

The reviewed R1 one-center H endpoint is not the omitted-`reference_spacing`
default. It must use explicit public `reference_spacing = 0.3` in the example
and endpoint test to reproduce the reviewed H baseline
`-0.49855234726272035`. Omitting `reference_spacing` leaves the general public
default at `1.0`, which is an allowed user choice but not the reviewed R1 H
baseline.

Fixed private choices in R1:

- `method = :pqs_source_box`
- `q_to_core_spacing_rule = :standard_pqs_ns_equals_q`
- `parent_axis_bundle_backend = :pgdg_localized_experimental`
- `parent_mapping_rule = :identity_mapping`

The current private H2 setup still needs an internal radius-like domain value.
The public z-axis H2 facade must derive that private value as
`max(xmax_parallel, xmax_transverse)`. This derived value is not a public
`basis` field.

For one-center H, public `reference_spacing` maps to the private one-center
mapping parameter historically named `parent_mapping_d`. The facade must not
accept public `d` or public `parent_mapping_d`, and no later stage may reset
that value to a hidden default. After initial lattice/parent construction, the
resolved value is provenance, not a staged algorithm input.

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

`xmax_parallel` and `xmax_transverse` are accepted only for z-axis H2.
They are unknown keys for one-center H.

Routing is implicit. One supported center selects the one-center base route.
Two supported centers must be z-axis H2 and select the z-axis diatomic base
route. R1 performs no translation or rotation. Cartesian x/y alignment and
arbitrary molecular orientation must throw `ArgumentError` and are deferred.
R1 has no public `method`, `route`, or output group selector. New method or
route keywords should wait until R2 or later provides a second real public
method.

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
    reference_spacing = 0.3,
)

h_ham = cartesian_base_hamiltonian(h_system; basis = h_basis)
```

Z-axis H2:

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

R1 amends the artifact contract narrowly: when `hamfile !== nothing`, the final
Hamiltonian file must also preserve the normalized public input provenance
needed by consumers to identify the run. This provenance records the accepted
public `system` and `basis` values after validation and default resolution,
including the resolved `reference_spacing` whether explicit or defaulted. It
must not be consumed by downstream construction stages after initial
lattice/parent setup, and it must not become a report, status object, wrapper
payload, separate manifest, or second artifact file.

This design does not approve:

- a new artifact file shape beyond the R1 normalized public input provenance
  stored in the final Hamiltonian file;
- a new artifact manifest;
- a separate basis/provenance artifact;
- a new materialization wrapper;
- durable status/result-kind fields.

## Error Behavior

Unsupported or malformed public requests must throw clear exceptions, normally
`ArgumentError` for input validation failures:

- missing required `system` or `basis` fields;
- unknown `system` or `basis` fields;
- unsupported atom count, route geometry, or non-H/H2 system in the first R1
  scope;
- non-z-axis H2 geometry, including x/y-aligned or generally oriented H2;
- inconsistent field lengths for symbols, charges, positions, or electron
  counts;
- nonpositive or nonfinite spacing, radius, extent, coordinate, or charge
  values;
- unsupported electron counts for H or H2;
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

`HP-R1-WIRE-01` owns one private shared constructor seam:

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

Approved owner file:

```text
src/pqs_source_box_low_order_materialization.jl
```

This PQS-named owner is acceptable for the explicitly PQS-only R1 migration.
It is not permanent method-neutral ownership for R2.

Allowed caller files:

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
  one-center H system and explicit `reference_spacing = 0.3`;
- H2 base Hamiltonian construction through the public facade;
- returned type is `CartesianIDAHamiltonian{Float64}`;
- H one-body baseline remains `-0.49855234726272035` within `1.0e-10` for the
  explicit-`reference_spacing = 0.3` public example;
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
- validation readback preserves normalized public input provenance, including
  the H example's `reference_spacing = 0.3`;
- R0 warm/cold baseline is not materially regressed without explanation;
- no `cartesian_pair_terms` or `cartesian_assembly` call appears in the
  recommended public example path.

Validation policy decision: ignored validation is not enough for a new public
facade. `HP-R1-TEST-01` approves one small committed public endpoint
test/example. The test/example should exercise only the public facade, H/H2
endpoint facts, and existing Hamiltonian-artifact readback. It must not assert
route-internal stage fields, status symbols, report mirrors, terminal role
vocabulary, or pair inventories.

Approved standalone integration gate:

```text
test/driver_public/cartesian_base_hamiltonian_runtests.jl
```

Invocation:

```text
julia --project=. test/driver_public/cartesian_base_hamiltonian_runtests.jl
```

Cadence: standalone integration/endpoint gate, not default tiny unit coverage.
Run it when touching the public facade, base producer wiring, Hamiltonian
artifact writer/reader, or terminal base-Hamiltonian path. `HP-R1-TEST-01` does
not approve adding this gate to `test/runtests.jl`; a runner edit would require
explicit repo-manager approval.

Test artifact behavior: use `mktempdir()` for `hamfile` output.

The H endpoint test must use the public example's explicit
`reference_spacing = 0.3`. It must not accept a hidden public `d` keyword or
`parent_mapping_d` field as a substitute.

The test must enforce the R1 geometry contract: x/y-aligned H2,
shifted-parallel H2, and generally oriented H2 fail with `ArgumentError`
before expensive construction, while the reviewed z-axis H2 endpoint succeeds.

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

The R1 base producer design forbids:

- new payload/result/status wrapper around `CartesianIDAHamiltonian`;
- new status, blocker, readiness, materialized, or result-kind field clouds;
- new report fields duplicating terminal bases, matrices, factors, or raw pair
  tensors;
- metadata-carried numerical data;
- new artifact shapes except the approved normalized public input provenance in
  the final Hamiltonian file;
- public solver controls;
- supplement, correction, fragment, counterpoise, or Cr2 stress functionality;
- broad compatibility adapters that preserve the old pair/assembly story.

If any of these become necessary, stop and write a separate docs-only design
amendment before source work.

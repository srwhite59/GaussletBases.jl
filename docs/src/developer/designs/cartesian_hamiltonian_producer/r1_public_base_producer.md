# R1 Public Base Producer

Status: approved design amendment for first R1 implementation. The later
`r1_one_center_base_atoms.md` amendment relaxes the one-center base atom scope
from H-only to explicit origin-centered all-electron atoms; this document
remains the baseline R1 H/H2 contract.

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
- `HP-R1-ART-01` - fixed producer-provenance schema in the final Hamiltonian
  file.
- `HP-R1-TEST-01` - small committed public endpoint test/example for the
  facade.

All five R1 IDs are approved for first implementation in `registry.md`.

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

The first implementation scope was origin-centered base H and Cartesian z-axis
aligned base H2. `HP-R1-ATOM-*` separately approves explicit origin-centered
all-electron one-center atoms through the same base facade and shared
atom/diatomic producer machinery. X/y-aligned diatomics, generally oriented
molecules, translated atoms, WL/QW unification, supplements, corrections,
solver handoff, and Cr2-scale performance remain later roadmap lanes unless
separately approved.

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
- first H/H2 scope supports only symbols and charges consistent with hydrogen
  nuclei; `HP-R1-ATOM-*` separately relaxes one-center atoms to explicit
  origin-centered all-electron inputs without element lookup;
- `nup` and `ndn` must be nonnegative integers and must give the supported
  total electron count for the requested H or H2 system;
- one-center H must be at `(0.0, 0.0, 0.0)`;
- z-axis H2 centers must have finite `x == 0`, finite `y == 0`, and distinct
  finite `z` coordinates.

`basis` is a `NamedTuple`.

Required common fields:

- `ns`
- `core_spacing`

Conditional fields:

- one-center H requires:
  - `radius`
- z-axis H2 requires:
  - `xmax_parallel`
  - `xmax_transverse`

`HP-COMP-NS-FN-01` amends the public size field. Durable public examples should
use `ns`, the requested cube/source/nesting size. Route-local `q` is derived
after selecting `nesting`: `q = ns` for `nesting = :pqs`, and `q = ns - 2` for
`nesting = :wl`. Legacy public `q` may remain temporarily only as
compatibility, with consistency checks when both `ns` and `q` are supplied.

`HP-COMP-ATOMBOX-FN-01` amends the one-center atom sizing contract: `radius`
is the public physical box extent authority for one-center atoms, and parent
axis counts must be derived from that extent plus `core_spacing` / the existing
spacing policy. `ns` remains source/nesting resolution metadata and must not be
interpreted as the direct parent side-count rule.

Optional public fields with R1 defaults:

- `parent_axis_family = :G10`
- `reference_spacing = 1.0`
- `tail_spacing = 10.0`

The reviewed R1 one-center H endpoint uses explicit public
`core_spacing = 0.3` plus the general `reference_spacing = 1.0` default.
Public `d` is deprecated and is not part of the durable producer contract.
If a temporary compatibility path accepts `d`, it must require
`d == resolved core_spacing`; mismatches must throw `ArgumentError`.

Fixed private choices in R1:

- `method = :pqs_source_box`
- `q_to_core_spacing_rule = :standard_pqs_ns_equals_q` for PQS after
  normalizing public `ns`
- `parent_axis_bundle_backend = :pgdg_localized_experimental`

The current private H2 setup still needs an internal radius-like domain value.
The public z-axis H2 facade must derive that private value as
`max(xmax_parallel, xmax_transverse)`. This derived value is not a public
`basis` field.

One-center H private wiring:

```text
spacing_inputs.reference_spacing = basis.reference_spacing
parent_inputs.parent_mapping_rule = :white_lindsey_atomic_mapping
parent_inputs.parent_mapping_d = resolved basis.core_spacing
```

`core_spacing` is the authoritative public near-core physical spacing. In the
White-Lindsey atom mapping, the internal `parent_mapping_d` is this same
resolved physical scale. The `Z` dependence is the mapping-shape rule
`core_range = sqrt(core_spacing / Z)` and
`mapping_strength = sqrt(core_spacing * Z)`, not a second public knob.
`reference_spacing`, `tail_spacing`, and box/domain controls remain separate
concepts. The facade must not accept public `parent_mapping_d`. Later
`HP-PQS-MAP-SFACTOR-*` approves only the expert scalar `s_factor` as a narrow
mapping-strength override; it does not revive public `d` or
`parent_mapping_d`.

These are not public keywords or accepted `basis` fields in R1:

- `n_s`
- `q` in durable public examples after `HP-COMP-NS-FN-01`
- `bond_axis`
- `bond_length`
- `map_backend`
- `parent_axis_counts`
- `parent_axis_bundle_backend`
- `parent_mapping_rule`
- `parent_mapping_Z`
- `parent_mapping_d`
- `d` in durable public examples

`xmax_parallel` and `xmax_transverse` are accepted only for z-axis H2.
They are unknown keys for one-center H. Public `d` is deprecated; if
temporarily accepted for one-center H, it must equal resolved `core_spacing`,
and it remains an unknown key for z-axis H2.

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
    ns = 5,
    core_spacing = 0.3,
    radius = 4.0,
    reference_spacing = 1.0,
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
    ns = 5,
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

## Artifact Provenance Schema

`HP-R1-ART-01` amends the artifact contract narrowly: when
`hamfile !== nothing`, the final Hamiltonian file must also preserve a fixed,
versioned producer-provenance record for consumers. The existing Hamiltonian
matrix payload remains the existing writer/readback payload; R1 provenance is
additional data in the same file, not a separate artifact or report.

Approved JLD2 provenance key prefix:

```text
producer_provenance/
```

Approved common keys:

| Key | Value |
| --- | --- |
| `producer_provenance/provenance_version` | `1` |
| `producer_provenance/producer` | `:cartesian_base_hamiltonian` |
| `producer_provenance/nesting` | public construction family, `:pqs` or `:wl` |
| `producer_provenance/route` | truthful base route label derived from `(input.kind, input.nesting)` |
| `producer_provenance/ns` | normalized public requested cube/source/nesting size |
| `producer_provenance/q` | derived route-local `q` consumed by route construction |
| `producer_provenance/q_rule` | `:pqs_ns_equals_q` or `:wl_ns_minus_2` |
| `producer_provenance/ns_source` | `:public_ns` or `:legacy_q_compatibility` |
| `producer_provenance/core_spacing` | public `basis.core_spacing::Float64` |
| `producer_provenance/reference_spacing` | resolved public `basis.reference_spacing::Float64` |
| `producer_provenance/tail_spacing` | resolved public `basis.tail_spacing::Float64` |
| `producer_provenance/parent_axis_family` | resolved public `basis.parent_axis_family::Symbol` |
| `producer_provenance/parent_axis_counts` | realized `NTuple{3,Int}` |
| `producer_provenance/mapping_kind` | route mapping symbol |
| `producer_provenance/mapping_d` | `Float64` for one-center atoms, `nothing` for z-axis diatomics |
| `producer_provenance/radius` | `Float64` for one-center atoms, `nothing` for z-axis diatomics |
| `producer_provenance/xmax_parallel` | `nothing` for one-center atoms, `Float64` for z-axis diatomics |
| `producer_provenance/xmax_transverse` | `nothing` for one-center atoms, `Float64` for z-axis diatomics |
| `producer_provenance/atom_symbols` | `Vector{String}` |
| `producer_provenance/nuclear_charges` | `Vector{Float64}` |
| `producer_provenance/atom_locations` | `Matrix{Float64}` with one row per center |
| `producer_provenance/nup` | `Int` |
| `producer_provenance/ndn` | `Int` |
| `producer_provenance/final_dimension` | `Int` |

Route-specific values:

- one-center H or explicit one-center atom with `nesting = :pqs`:
  - `producer_provenance/route = :one_center_pqs_base`
  - `producer_provenance/nesting = :pqs`
  - `producer_provenance/mapping_kind = :white_lindsey_atomic_mapping`
  - `producer_provenance/mapping_d = resolved basis.core_spacing`
  - `producer_provenance/radius = basis.radius`
  - `producer_provenance/xmax_parallel = nothing`
  - `producer_provenance/xmax_transverse = nothing`
- one-center H or explicit one-center atom with `nesting = :wl`:
  - `producer_provenance/route = :one_center_wl_base`
  - `producer_provenance/nesting = :wl`
  - `producer_provenance/mapping_kind = :white_lindsey_atomic_mapping`
  - `producer_provenance/mapping_d = resolved basis.core_spacing`
  - `producer_provenance/radius = basis.radius`
  - `producer_provenance/xmax_parallel = nothing`
  - `producer_provenance/xmax_transverse = nothing`
- z-axis H2, or explicit homonuclear z-axis all-electron diatomic after
  `HP-COMP-BASEDIAT-FN-01`, with `nesting = :pqs`:
  - `producer_provenance/route = :z_axis_diatomic_pqs_base`
  - `producer_provenance/nesting = :pqs`
  - `producer_provenance/mapping_kind = :multicenter_pqs_mapping`
  - `producer_provenance/mapping_d = nothing`
  - `producer_provenance/radius = nothing`
  - `producer_provenance/xmax_parallel = basis.xmax_parallel`
  - `producer_provenance/xmax_transverse = basis.xmax_transverse`
- z-axis H2, or explicit homonuclear z-axis all-electron diatomic after
  `HP-COMP-BASEDIAT-FN-01`, with `nesting = :wl`, after
  `HP-COMP-WLDIAT-FN-01` succeeds:
  - `producer_provenance/route = :z_axis_diatomic_wl_base`
  - `producer_provenance/nesting = :wl`
  - `producer_provenance/mapping_kind = resolved WL diatomic parent mapping symbol`
  - `producer_provenance/mapping_d = nothing`
  - `producer_provenance/radius = nothing`
  - `producer_provenance/xmax_parallel = basis.xmax_parallel`
  - `producer_provenance/xmax_transverse = basis.xmax_transverse`

The `:z_axis_diatomic_wl_base` value is a truthful route-label value under the
existing schema, not an artifact schema change. It may be written only after the
WL z-axis diatomic base path succeeds under `HP-COMP-WLDIAT-FN-01` /
`HP-COMP-WLDIAT-TEST-01`. Unsupported `(kind, nesting)` combinations must throw
before artifact writing rather than writing a PQS-oriented route label.

`format_version` for the existing Cartesian IDA Hamiltonian matrix payload is
not changed by R1. Existing `read_cartesian_ida_hamiltonian` must continue to
read the Hamiltonian matrices while ignoring the provenance keys. Validation
may inspect the `producer_provenance/` keys directly. A public provenance
reader, a separate manifest, or a new Hamiltonian wrapper requires another
design amendment.

The provenance is for consumer tracking. It must not be consumed by downstream
construction stages after initial lattice/parent setup, and it must not become
a report, status object, wrapper payload, separate manifest, or second artifact
file.

This design does not approve:

- a new artifact file shape beyond the `HP-R1-ART-01`
  `producer_provenance/` keys stored in the final Hamiltonian file;
- a new artifact manifest;
- a separate basis/provenance artifact;
- a new materialization wrapper;
- durable status/result-kind fields.

## Error Behavior

Unsupported or malformed public requests must throw clear exceptions, normally
`ArgumentError` for input validation failures:

- missing required `system` or `basis` fields;
- unknown `system` or `basis` fields;
- unsupported atom count, route geometry, or unsupported system in the active
  R1/R1-atom scope;
- non-z-axis diatomic geometry, including x/y-aligned or generally oriented
  diatomics;
- inconsistent field lengths for symbols, charges, positions, or electron
  counts;
- nonpositive or nonfinite spacing, radius, extent, coordinate, or charge
  values;
- unsupported electron counts, including non-neutral all-electron diatomics;
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

The old route-driver execution spine was:

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

`cartesian_materialization` and the paired route-driver print/save wrappers are
now approved for retirement under `HP-RETIRE-DRV-MAT-*`. Do not copy this old
spine into the canonical driver or future public workflow.

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

Later authority `HP-PQS-COULOMB-ACCURACY-FN-01` supersedes the implicit
expansion part of this signature. A focused caller scan should delete this
helper if it is no longer live. If it remains live, it must receive the
already-resolved `CoulombGaussianExpansion` explicitly and must not select
compact accuracy internally.

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
  one-center H system, explicit `core_spacing = 0.3`, and
  `reference_spacing = 1.0`;
- H2 base Hamiltonian construction through the public facade;
- returned type is `CartesianIDAHamiltonian{Float64}`;
- H one-body baseline remains `-0.49855234726272035` within `1.0e-10` for the
  `core_spacing = 0.3`, `reference_spacing = 1.0` public example;
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
- temporary `d` for one-center H, if accepted, must match resolved
  `core_spacing`; mismatches throw `ArgumentError`;
- `d` for z-axis H2 throws `ArgumentError`;
- empty `hamfile` throws `ArgumentError`;
- non-`nothing` `hamfile` writes with `write_cartesian_ida_hamiltonian`;
- validation readback with `read_cartesian_ida_hamiltonian` has zero one-body
  readback delta;
- artifact validation preserves the fixed `producer_provenance/` keys,
  including H `core_spacing = 0.3`, `reference_spacing = 1.0`,
  `mapping_kind = :white_lindsey_atomic_mapping`, and resolved
  `mapping_d = 0.3`;
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
`core_spacing = 0.3` with `reference_spacing = 1.0`. It must reject public
`parent_mapping_d` and any temporary `d` value that differs from resolved
`core_spacing`.

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
- new artifact shapes except the approved `HP-R1-ART-01`
  `producer_provenance/` keys in the final Hamiltonian file;
- public solver controls;
- supplement, correction, fragment, counterpoise, or Cr2 stress functionality;
- broad compatibility adapters that preserve the old pair/assembly story.

If any of these become necessary, stop and write a separate docs-only design
amendment before source work.

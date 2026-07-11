# R1 One-Center Base Atoms

Status: approved narrow source authority to relax the one-center base producer
from hydrogen-only validation to explicit origin-centered all-electron atoms.
This is not supplemented-atom, translated-atom, ECP, solver, or broad molecule
authority.

## Decision

Approve `cartesian_base_hamiltonian(...)` to accept explicit one-center base
atom systems beyond H when the atom is origin-centered and all physical data is
supplied by the caller.

This amendment extends the existing R1 public base facade input scope. It does
not create a new public function, new export, new driver mode, new artifact
schema, or new route vocabulary. The route remains the unsupplemented,
uncorrected, all-electron localized-IDA base Hamiltonian route.

The current committed endpoint remains origin-centered H. Non-H atoms such as
Be or Cr are allowed as explicit user-run or ignored validation/stress inputs
after H remains clean. A committed non-H reference gate requires a later
decision with an accepted reference value or endpoint criterion.

## Shared Workflow Constraint

Atoms and diatomics must share the same producer workflow except at the narrow
places where the physics actually differs. The approved public workflow remains:

```text
system / specification
-> parent and route geometry
-> terminal basis realization
-> Hamiltonian production
-> optional existing Hamiltonian artifact
```

One-center atoms may have atom-specific geometry and shellification inputs, and
future supplemented atoms may have atom-specific residual placement. Those
differences must still feed the same underlying terminal-basis, one-body, IDA,
Hamiltonian-construction, and artifact-writing machinery whenever possible.

This amendment does not approve:

- an atom-only Hamiltonian builder;
- a parallel atom materialization path;
- separate atom route-stage/report/status objects;
- atom-only one-body or IDA orchestration when the shared kernels already
  apply;
- provenance or metadata fields used to route algorithmic data around the
  shared workflow.

If implementation discovers that a non-H atom cannot use the existing shared
workflow after geometry/shellification normalization, it must stop and report
the exact missing shared seam rather than adding an atom-specific workaround.

## Approved IDs

- `HP-R1-ATOM-FN-01` - explicit one-center all-electron base atom facade scope.
- `HP-R1-ATOM-WIRE-01` - one-center atom normalization to the existing R1 base
  construction path and shared atom/diatomic producer machinery.
- `HP-R1-ATOM-TEST-01` - validation gates for the one-center base atom
  relaxation.

## Approved File

Approved source file:

```text
src/cartesian_base_hamiltonian.jl
```

No other `src`, `test`, `tools`, `bin`, or committed input-fixture file is
approved by this source lane.

## Public Input Scope

The public call shape remains unchanged:

```julia
cartesian_base_hamiltonian(
    system::NamedTuple;
    basis::NamedTuple,
    hamfile::Union{Nothing,AbstractString} = nothing,
)::CartesianIDAHamiltonian{Float64}
```

Approved one-center atom `system` rules:

- exactly one center;
- `atom_symbols`, `nuclear_charges`, and `atom_locations` are vectors or other
  `AbstractVector` values, not variable-size tuples;
- the single atom location is exactly `(0.0, 0.0, 0.0)` in the supported
  public contract;
- the nuclear charge is supplied explicitly, finite, positive, and
  integer-valued;
- `nup` and `ndn` are explicit nonnegative integers;
- neutral all-electron count is derived from the explicit charge:
  `nup + ndn == round(Int, only(nuclear_charges))`;
- the atom symbol is provenance/user labeling only and must not be used to infer
  charge, electron count, spin, basis, or ECP behavior.

Approved one-center atom `basis` rules:

- required durable public fields: `ns`, `core_spacing`, and `radius`;
- route-local `q` is derived from `ns` and `nesting` under
  `HP-COMP-NS-FN-01`;
- temporary legacy `q` compatibility, if kept, must normalize to `ns` and
  reject inconsistent `ns`/`q` pairs;
- optional fields retain existing R1 defaults:
  `reference_spacing = 1.0`, `tail_spacing = 10.0`, and
  `parent_axis_family = :G10`;
- public `d` is deprecated and is not part of the durable atom contract; if a
  temporary compatibility path accepts it, it must equal resolved
  `core_spacing`;
- public `parent_mapping_Z`, `parent_mapping_d`, `parent_mapping_rule`,
  backend controls, and parent axis counts remain unsupported.

Private one-center mapping remains:

```text
spacing_inputs.reference_spacing = basis.reference_spacing
parent_inputs.parent_mapping_rule = :white_lindsey_atomic_mapping
parent_inputs.parent_mapping_d = resolved basis.core_spacing
parent_inputs.parent_mapping_Z = only(system.nuclear_charges)
```

The public charge, not the atom symbol, supplies `Z`.
The White-Lindsey `Z` dependence remains the internal mapping-shape rule
`core_range = sqrt(core_spacing / Z)` and
`mapping_strength = sqrt(core_spacing * Z)`. Future automatic presets may
derive `core_spacing` from `Z`, for example through a fixed
`core_spacing * Z` family, but after resolution `core_spacing` is the single
authoritative near-nucleus spacing. `reference_spacing`, `tail_spacing`, and
box controls remain separate concepts. Later `HP-PQS-MAP-SFACTOR-*` approves
only the expert scalar `s_factor` as a narrow mapping-strength override; it
does not revive public `d`, public `parent_mapping_d`, or element defaults.

## Physical Parent Extent

`HP-COMP-ATOMBOX-FN-01` is implemented. Public `basis.radius` is the
one-center physical parent-extent authority. The producer maps that radius
through the existing White-Lindsey atomic mapping and spacing policy, rounds
to an odd axis count, and applies the public-`ns` direct-core side only as a
minimum:

```text
mapped_count = odd count covering basis.radius
direct_core_side = isodd(ns) ? ns : ns + 1
parent_side = max(mapped_count, direct_core_side)
```

Thus `ns` controls source/nesting resolution; it does not replace atom
radius or driver padding as physical box size. The shared `ns`/`q` rules
are canonical in
[Nesting/supplement composition](nesting_supplement_composition_plan.md), and
direct-core parity is canonical in
[Public ns direct-core side parity](public_ns_core_side_parity.md).

## Artifact Contract

No new artifact schema is approved. Existing `HP-R1-ART-01`
`producer_provenance/` keys are sufficient for one-center base atoms:

- `route = :one_center_pqs_base`;
- `mapping_kind = :white_lindsey_atomic_mapping`;
- `mapping_d = resolved basis.core_spacing`;
- `radius = basis.radius`;
- `atom_symbols`, `nuclear_charges`, `atom_locations`, `nup`, `ndn`, and
  `final_dimension` record the explicit public input and result.

Do not add `mapping_Z`, element-table fields, spin labels, ECP fields, status
fields, a separate manifest, or a provenance reader in this lane.

## Validation

`HP-R1-ATOM-TEST-01` approves validation only for this one-center base atom
relaxation:

- `git diff --check`;
- package load;
- existing origin-centered H public facade endpoint remains unchanged,
  including the reviewed `core_spacing = 0.3`, `reference_spacing = 1.0`
  baseline and internal `parent_mapping_d = core_spacing`;
- optional ignored/user-run Be or Cr one-center base atom artifact
  write/readback using explicit charge, spin sectors, origin geometry, and
  basis controls;
- finite/symmetric `K`, unit `U_A`, and IDA `V` for ignored/user-run non-H
  atom checks;
- clear `ArgumentError` for translated atom input, mismatched temporary `d`,
  noninteger or nonpositive charge, nonneutral electron count, or
  element-table/default requests where practical.

No new committed test file, committed non-H atom fixture, public reference
scalar, solver run, supplemented atom endpoint, ECP gate, or translated-atom
gate is approved by this ID.

Temporary validation inputs should live under ignored `tmp/work`.

## Forbidden

This amendment does not approve:

- supplemented atom Hamiltonians;
- Residual Gaussian, MWG/IDA, raw-block, or terminal-kernel convention changes;
- translated one-center atoms;
- multi-center atom/molecule broadening beyond existing R1 H2 and approved ZDI
  supplemented diatomics;
- element lookup tables, inferred nuclear charges, inferred spin, or
  element-specific defaults;
- ECP or pseudopotential support;
- solver/RHF workflow;
- public API redesign or new export;
- artifact schema changes or new provenance keys;
- route diagnostics, report/status/payload fields, metadata-carried numerical
  data, or new route objects;
- driver changes beyond already-approved `HP-DRV-ATOM-*` authority;
- committed atom fixture files or new committed tests.

## Line Budget And Failure Rule

Line budget:

- at most `80` added `src` lines;
- net simplification expected where H-specific validation is replaced by
  explicit one-center atom validation;
- no new committed test, tool, driver, or input-fixture file.

Failure rule: if implementation requires source edits outside
`src/cartesian_base_hamiltonian.jl`, changes to private materialization
owners, new artifact keys, translated atom support, supplemented atom support,
ECP behavior, solver workflow, element lookup/default tables, committed
fixtures/tests, route/report/status/payload expansion, or an atom-only
producer/materialization path, stop and request a separate docs-only amendment.

# Cartesian Driver Atom Workflow

Status: approved narrow canonical-driver workflow for explicit one-center base
atom inputs. This is driver workflow authority only. It does not broaden the
base producer, Residual Gaussian, artifact schema, solver, or public API
contracts.

## Decision

Approve a compact atom workflow in `bin/cartesian_ham_builder.jl` for
origin-centered one-center atoms in `mode = :base`.

This approval is for base atom driver output only. The driver may accept an
explicit one-center atom system, validate the visible workflow input, and call
the existing `cartesian_base_hamiltonian(...)` facade with `hamfile` when the
base facade already supports the requested atom. Current committed validation
remains the approved origin-centered H endpoint. The separate `HP-R1-ATOM-*`
amendment may broaden the base facade to explicit origin-centered
all-electron atoms; the driver may consume that support through the same
base-facade call without new driver authority.

Supplemented atom Hamiltonians are not approved in this pass. The Residual
Gaussian owner-local selection algorithm is compatible in spirit with a single
owner, but the current supported supplemented facade is scoped to diatomics.
One-center supplemented atoms need a separate amendment that names the facade
surface, basis-loading path, validation scalar or endpoint, and artifact
provenance behavior.

## Approved IDs

- `HP-DRV-ATOM-FN-01` - canonical-driver explicit base atom workflow.
- `HP-DRV-ATOM-WIRE-01` - driver-to-base-facade wiring for one-center atoms.
- `HP-DRV-ATOM-TEST-01` - validation gates for the base atom workflow.

## Approved File

Approved source file:

```text
bin/cartesian_ham_builder.jl
```

No `src`, `test`, `tools`, other `bin`, or committed input-fixture file is
approved by this atom-driver lane.

## System Scope

Approved system input shape:

- `atom_symbols::AbstractVector`;
- `nuclear_charges::AbstractVector`;
- `atom_locations::AbstractVector{<:NTuple{3}}`;
- `nup`;
- `ndn`.

Approved first workflow scope:

- exactly one center;
- origin-centered location `(0.0, 0.0, 0.0)`;
- finite positive nuclear charge supplied explicitly;
- `nup` and `ndn` are explicit nonnegative integers;
- neutral all-electron count uses the explicit charge:
  `nup + ndn == round(Int, only(nuclear_charges))`;
- no element lookup table, element-specific default, isotope/default-charge
  inference, ECP, pseudopotential, or solver behavior.

Translated atoms are not approved. One-center base construction and the
reviewed R1 endpoint use the origin-centered atom contract; supporting
translated atoms would need a separate design decision for mapping,
provenance, and endpoint validation.

## Basis Scope

The driver may pass the explicit base `basis` group to the existing base facade.
For the current approved H endpoint, the visible basis fields include:

- `q`;
- `core_spacing`;
- `radius`;
- explicit `d`;
- optional `reference_spacing`;
- optional `tail_spacing`;
- optional `parent_axis_family`.

`d` remains the public one-center mapping control and has no driver-created
hidden default. `reference_spacing` remains a separate reference-grid spacing.

## Wiring Contract

`HP-DRV-ATOM-FN-01` may add driver-local input normalization for one-center
base atom systems. It must keep the driver compact and copyable: visible
defaults, optional trusted input file, command-line overrides, compact summary,
coarse timing, artifact write, and optional readback.

`HP-DRV-ATOM-WIRE-01` may call:

```julia
cartesian_base_hamiltonian(system; basis, hamfile)
```

for base atom construction. The driver must not compose package-internal
route-stage helpers, terminal basis objects, raw-block providers, or
report/status/payload objects.

If the requested one-center atom is outside the existing base facade support,
the driver or facade must throw a clear `ArgumentError` before or at the
supported producer boundary. This driver approval does not itself authorize
changing `src/cartesian_base_hamiltonian.jl`; broader one-center base atom
facade support is governed by `HP-R1-ATOM-*`.

## Validation

`HP-DRV-ATOM-TEST-01` approves validation only for the base atom driver
workflow:

- `git diff --check`;
- package load;
- origin-centered H base driver artifact write/readback using explicit
  `atom_symbols`, `nuclear_charges`, `atom_locations`, `nup`, `ndn`, and
  one-center basis fields;
- optional ignored negative checks for non-origin atom input, nonneutral
  electron count, missing `d`, or unsupported atom input;
- no committed atom fixture, committed test file, solver run, supplemented
  atom endpoint, or translated-atom gate.

Temporary project input files for validation should live under ignored
`tmp/work`.

## Forbidden

This amendment does not approve:

- supplemented atom Hamiltonians;
- changes to Residual Gaussian selection, MWG/IDA conventions, raw blocks, or
  terminal kernels;
- changes to `src/cartesian_base_hamiltonian.jl` or any other source file
  outside the canonical driver, except under separate `HP-R1-ATOM-*`
  authority;
- general atom support beyond the existing base facade and `HP-R1-ATOM-*`;
- translated atoms;
- element-specific defaults beyond visible example inputs;
- ECP or pseudopotential support;
- solver/RHF workflow;
- new artifact schema or provenance keys;
- public API/export redesign;
- committed atom fixture files or committed tests;
- private route-stage controls, diagnostics, metadata/status/report fields, or
  payload objects.

## Line Budget And Failure Rule

Line budget:

- at most `80` added `bin` lines;
- net simplification preferred if existing H-only driver branches are replaced
  by explicit one-center input validation;
- no new committed test, tool, or input-fixture file.

Failure rule: if implementation requires source changes outside
`bin/cartesian_ham_builder.jl`, supplemented atom support, translated atoms,
ECP, solver workflow, artifact schema changes, route diagnostics,
metadata/status/report fields, committed fixtures/tests, or element
lookup/default tables, stop and request a separate docs-only amendment.

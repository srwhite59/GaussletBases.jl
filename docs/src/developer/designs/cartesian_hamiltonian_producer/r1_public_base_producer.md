# R1 Public Base Producer Candidate

Status: candidate design amendment, not implementation authority.

This document defines the proposed minimal public base Cartesian Hamiltonian
producer for first H/H2 use. It does not approve source work, tests, tools,
driver edits, new report fields, new artifact formats, or implementation of
the candidate IDs below.

## Candidate IDs

- `HP-R1-FN-01` - public base Hamiltonian producer facade.
- `HP-R1-WIRE-01` - report-free base producer wiring from the facade to the
  existing terminal-basis and Hamiltonian construction path.

Both IDs are candidate-only until explicitly approved in `registry.md`.

## Public Call Shape

The proposed public entry point is a single function:

```julia
cartesian_base_hamiltonian(system;
    basis,
    method = :pqs_source_box,
    route = :auto,
    output = (;),
)::CartesianIDAHamiltonian{Float64}
```

The function returns the existing `CartesianIDAHamiltonian{Float64}` directly.
It must not return a wrapper, status object, materialization payload, report
mirror, or `(value, status)` pair.

The first implementation scope is base H and bond-aligned base H2. Broader
atoms, general molecules, WL/QW unification, supplements, corrections, solver
handoff, and Cr2-scale performance remain later roadmap lanes unless separately
approved.

## Inputs

`system` describes the physical system:

- `atom_symbols`
- `nuclear_charges`
- `atom_locations`
- `nup`
- `ndn`
- optional geometry facts such as `bond_axis`, `bond_length`, and `radius`

`basis` describes basis and spacing choices:

- `q`
- `n_s`
- `reference_spacing`
- `tail_spacing`
- `core_spacing`
- `xmax_parallel`
- `xmax_transverse`
- parent-axis family and backend choices needed by the existing route

`method` / `route` choose the construction method:

- first candidate method: `:pqs_source_box`
- first candidate route: automatic one-center atomic or bond-aligned diatomic
  base PQS route selection from `system`
- explicit route override may be allowed only for reviewed H/H2 fixtures

`output` controls optional artifact behavior:

- `hamfile`
- `save_ham_artifact::Bool`
- no basis artifact is part of this base R1 contract
- no TSV/report artifact is part of this base R1 contract

Grouped inputs are call-boundary convenience only. They must not become
persistent payload structs, status summaries, metadata field clouds, or report
mirrors.

## Output And Artifact Contract

The output is the existing `CartesianIDAHamiltonian{Float64}`.

If `save_ham_artifact = true`, the implementation must use the existing
`write_cartesian_ida_hamiltonian` writer and existing readback path. This
candidate does not approve:

- a new artifact file shape;
- a new artifact manifest;
- a new basis/provenance artifact;
- a new materialization wrapper;
- durable status/result-kind fields.

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
- optional artifact write/readback uses `write_cartesian_ida_hamiltonian` and
  has zero one-body readback delta;
- R0 warm/cold baseline is not materially regressed without explanation;
- no `cartesian_pair_terms` or `cartesian_assembly` call appears in the
  recommended public example path.

No new committed smoke/test file is approved by this candidate. If a committed
public example or smoke is needed, add or approve a separate test/example ID
before implementation.

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

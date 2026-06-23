# Cartesian Driver Usability Workflow

Status: approved narrow source authority for making
`bin/cartesian_ham_builder.jl` a compact functional Hamiltonian producer
driver. This is workflow authority only; it is not new Hamiltonian algorithm,
raw-kernel, solver, or diagnostic authority.

## Decision

The previous canonical-driver freeze served its purpose: it kept private route
diagnostics, ladder controls, and provider switches from shaping the production
implementation. The base and residual-GTO/MWG producer paths are now coherent
enough that the canonical driver should prove those paths work together.

Approve `bin/cartesian_ham_builder.jl` as the standard compact, copyable,
human-facing Cartesian Hamiltonian producer driver. It should follow the
scientific-driver pattern:

```text
visible editable defaults
-> optional trusted project input file
-> command-line key=value overrides
-> compact normalized run summary
-> coarse timed user-facing phases
-> supported base or supplemented Hamiltonian construction
-> artifact write
-> optional readback check
```

Producing a Hamiltonian artifact is the full success endpoint. Staged
milestones are allowed only when they are user-facing workflow milestones, not
private route-stage diagnostics.

## Approved IDs

- `HP-DRV-FILE-01` - canonical Cartesian driver source file.
- `HP-DRV-FN-01` - compact functional driver workflow.
- `HP-DRV-TEST-01` - validation gates for the driver usability lane.

## Approved File

```text
bin/cartesian_ham_builder.jl
```

No other `bin`, `tools`, `src`, `test`, or committed input-fixture file is
approved by this lane. Project-specific copied drivers are acceptable consumer
artifacts, but they are not part of the canonical driver authority.

## Input Shape

The canonical invocation should support:

```text
julia --project=. bin/cartesian_ham_builder.jl [input.jl] [key=value ...]
```

Approved input policy:

- keep visible editable defaults near the top of the driver;
- optionally read one trusted local Julia input file for project-specific
  defaults;
- apply later command-line `key=value` overrides after the input file;
- treat the input file and overrides as trusted local scientific workflow
  inputs, not as a security boundary or safe public parser;
- keep the normalized internal configuration compact.

Approved user-facing configuration concepts:

- `mode = :base` or `mode = :supplemented`;
- `system` specification;
- base `basis` specification;
- optional `supplement` specification for the supported residual-GTO/MWG path;
- `hamfile`;
- coarse `timing`;
- compact `print_summary`;
- optional `readback_check`.

The implementation may choose a compact concrete representation, such as a
small `NamedTuple` or driver-local variables, but it must not grow route-stage
payloads, report mirrors, status objects, provider bundles, or metadata
carriers.

## Allowed Workflow

The driver may call only supported producer surfaces:

- the approved base producer facade for base H/H2 scope;
- the approved non-exported residual-GTO/MWG usability facade for supported
  supplemented H2 and internal/performance-supported Be2 scope;
- the approved Hamiltonian artifact writer and readback check.

Allowed user-facing phases:

- validate/normalize input;
- build base Hamiltonian;
- build supplemented Hamiltonian;
- write artifact;
- read back artifact for validation;
- print compact summary and coarse timing.

Driver printing and timing are allowed when they use user-facing labels and
small summaries. Examples include mode, system, basis label, final dimension,
artifact path, self-Coulomb/check scalar for known endpoints, elapsed phase
times, and readback deltas.

## Not Approved

This lane does not approve:

- private route-stage controls, stop-after internals, ladder probes, stage
  markers, fixture hacks, or diagnostic knobs;
- underscored package helper calls from the driver;
- raw-block provider switches;
- residual algorithm, MWG/IDA, raw-block, or terminal-kernel changes;
- route reports, status symbols, payload dumps, metadata field clouds, or
  report mirrors;
- allocation probes, benchmarking harnesses, or Cr2-specific stress controls;
- solver/RHF/ECP/EGOI/HamV6 workflow;
- public API/export changes;
- artifact schema changes beyond the existing approved provenance groups;
- committed test files or committed driver-input fixtures.

Diagnostics and ladder probing remain in `tools/` or ignored `tmp/work`
scripts, not in the canonical driver.

## Validation

`HP-DRV-TEST-01` approves validation only for the compact driver workflow:

- `git diff --check`;
- package load;
- H2 base driver run writes a `CartesianIDAHamiltonian` artifact and optional
  readback passes;
- H2 supplemented driver run writes a supplemented
  `CartesianIDAHamiltonian` artifact with the approved compact
  `supplement_provenance/` group and optional readback passes;
- optional ignored Be2 usability run if the implementation touches the
  supplemented mode;
- generic explicit homonuclear z-axis Cr2 stress is approved only by
  `HP-R3U-ZDI-WIRE-01` and remains ignored/user-run after H2/Be2 validation;
- no Cr2-specific driver run, solver run, committed test file, or committed
  input fixture is approved by this ID.

Temporary project input files for validation should live under ignored
`tmp/work`.

## Line Budget And Failure Rule

Line budget:

- at most `150` added `bin` lines;
- the canonical driver should remain compact and copyable;
- no new committed test or tool file.

Failure rule: if implementation requires a parser framework, source files
outside `bin/cartesian_ham_builder.jl`, committed input fixtures, route-stage
diagnostics, status/report/payload expansion, artifact schema changes, or
Cr2-specific workflow support, stop and request a new docs-only amendment.

## Deferred

Cr2 facade/artifact workflow, public API/export, solver/RHF/ECP/EGOI/HamV6,
diagnostic harnesses, route-stage cleanup, and broad driver feature polish
remain deferred.

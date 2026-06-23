# Cartesian Driver Usability Workflow

Status: approved narrow source authority for making
`bin/cartesian_ham_builder.jl` a compact functional Hamiltonian producer
driver, plus a narrow internal staged producer surface in
`src/cartesian_base_hamiltonian.jl` so the driver can execute visible
physics-level construction stages. This is workflow authority only; it is not
new Hamiltonian algorithm, raw-kernel, solver, or diagnostic authority.

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
-> visible public system / basis / optional supplement contract construction
-> compact normalized run summary
-> coarse timed user-facing physics-level construction stages
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
- `HP-DRV-STAGE-FN-01` - internal visible physics-stage producer surface.
- `HP-DRV-STAGE-WIRE-01` - canonical driver wiring to the staged surface.
- `HP-DRV-STAGE-TEST-01` - validation gates for staged driver execution.
- `HP-DRV-TEST-01` - validation gates for the driver usability lane.

## Approved File

```text
bin/cartesian_ham_builder.jl
src/cartesian_base_hamiltonian.jl
```

`src/cartesian_base_hamiltonian.jl` is approved only for the staged producer
surface described below. No other `bin`, `tools`, `src`, `test`, or committed
input-fixture file is approved by this lane. Project-specific copied drivers
are acceptable consumer artifacts, but they are not part of the canonical
driver authority.

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

- `basisname = nothing` for base mode, or `basisname !== nothing` for
  supported supplemented diatomic mode;
- `system` specification;
- base `basis` specification;
- optional `supplement` specification for the supported residual-GTO/MWG path;
- `hamfile`;
- `padding`;
- `check_file`;
- `print_contract`;
- `print_timing`;
- `expected_dimension`.

Compact summary printing and artifact readback checks may remain part of the
workflow, but they are not open-ended run-level hooks. They must not introduce
additional route, diagnostic, artifact-schema, or solver controls.

The implementation may choose a compact concrete representation, such as a
small `NamedTuple` or driver-local variables, but it must not grow route-stage
payloads, report mirrors, status objects, provider bundles, or metadata
carriers.

## Public Contract Construction

The driver must visibly construct the public contract objects before it calls
an approved producer facade:

```text
system
basis
optional supplement
hamfile
```

For base runs, `basisname = nothing` selects the base facade path. For
supplemented runs, `basisname !== nothing` selects the supported residual-GTO
/ MWG supplemented diatomic facade and becomes the visible supplement basis
label. `basisname !== nothing` must reject `Natom == 1`; supplemented atoms
remain unapproved until a separate facade/RG amendment approves them.

`padding` is a public physical box-padding control, not a route-stage field.
For one-center atoms, it means the box padding around the atom and maps to the
one-center base facade's `radius`. For z-axis diatomics, it means the padding
around the two nuclei and maps internally to the existing facade extents; under
the current origin-based z-axis contract this is equivalent to
`xmax_parallel = max(abs(z_i)) + padding` and
`xmax_transverse = padding`. If a future translated/general geometry needs
different box-centering semantics, that requires a separate amendment.

The driver may print the public contract when `print_contract = true`. The
printed contract is limited to the public `system`, `basis`, optional
`supplement`, `hamfile`, and run-level hooks. It must not print route-stage
objects, private helper inputs, raw-provider choices, report payloads,
metadata clouds, artifact schemas, or internal matrices.

## Public Hooks

The only approved compact run-level hooks are:

- `check_file`;
- `print_contract`;
- `print_timing`;
- `expected_dimension`.

These hooks exist for human expert review and Codex-controlled artifact checks.
`print_timing` may print coarse user-facing phase timings. `expected_dimension`
may validate the returned Hamiltonian dimension and throw a clear error on
mismatch. `check_file` may write a compact check record containing public
contract facts, artifact path, final dimension, expected-dimension result,
readback deltas, and coarse timing. It must not become an artifact schema,
provenance dump, status object, route report, allocation log, or solver output.

## Visible Physics-Level Stages

The canonical driver must not hide the supplemented construction entirely
inside one opaque facade call. It should execute and optionally time visible
physics-level stages while still avoiding route-internal choreography.

Approved visible stages are:

- construct public `system`, `basis`, and optional `supplement`;
- build the base working basis / terminal realization and base Hamiltonian;
- load or build the Gaussian supplement basis when `basisname !== nothing`;
- build residual Gaussian augmentation;
- build exact augmented operators;
- assemble the supplemented Hamiltonian;
- write and check the artifact.

These are workflow stages, not old route stages. They must not expose
`cartesian_parent`, `cartesian_shells`, `cartesian_units`,
`cartesian_pair_terms`, `cartesian_assembly`, reports, route skeletons,
retained-rule internals, raw-block providers, or diagnostic stop points.

`HP-DRV-STAGE-FN-01` approves a small non-exported, non-underscored,
driver-facing staged producer surface in `src/cartesian_base_hamiltonian.jl`.
The exact function names may follow local style, but their purpose must match
the approved construction stages above. Public contract construction remains a
driver responsibility, and artifact write/check remains an approved writer and
readback responsibility. The source owner may factor the existing base and
residual-GTO/MWG facade bodies into these stage functions so the canonical
driver can call named construction stages without calling underscored helpers
directly.

The staged surface may return existing domain objects and small fixed-key
ephemeral stage products needed by the next approved stage. It must not create
a broad payload/report/status object, route-stage object, persistent cache,
runtime-keyed field cloud, artifact schema, or new public API/export.

The existing one-call facades may remain as convenience wrappers around the
same staged implementation. The canonical driver should prefer the staged
surface when it needs visible execution and timing.

## Allowed Workflow

The driver may call only supported producer surfaces:

- the approved staged producer surface for base working-basis/base-Hamiltonian
  construction, including the origin-centered atom workflow approved by
  `HP-DRV-ATOM-*`;
- the approved non-exported residual-GTO/MWG usability facade for supported
  supplemented H2 and internal/performance-supported Be2 scope, or the staged
  producer surface that factors that facade into visible physics-level stages;
- the approved Hamiltonian artifact writer and readback check.

Allowed user-facing phases:

- validate/normalize input;
- construct public contract;
- build base working basis and base Hamiltonian;
- build Gaussian supplement;
- build residual augmentation;
- build exact augmented operators;
- assemble supplemented Hamiltonian;
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
- private helper contract construction instead of public `system` / `basis` /
  `supplement` objects;
- opaque one-call supplemented construction as the only canonical-driver
  execution shape;
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
- public contract printing/checking for at least one base run when the driver
  construction code changes;
- visible stage timing or summary for base construction and supplemented H2
  construction when staged wiring changes;
- H2 base driver run writes a `CartesianIDAHamiltonian` artifact and optional
  readback passes;
- origin-centered H atom driver run writes a `CartesianIDAHamiltonian` artifact
  and optional readback passes under `HP-DRV-ATOM-TEST-01`;
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

- at most `150` added `bin` lines and at most `150` added `src` lines;
- the canonical driver should remain compact and copyable;
- no new committed test or tool file.

Failure rule: if implementation requires a parser framework, source files
outside `bin/cartesian_ham_builder.jl` and `src/cartesian_base_hamiltonian.jl`,
committed input fixtures, route-stage diagnostics, status/report/payload
expansion, artifact schema changes, public API/export changes, or Cr2-specific
workflow support, stop and request a new docs-only amendment.

## Deferred

Cr2 facade/artifact workflow, public API/export, solver/RHF/ECP/EGOI/HamV6,
diagnostic harnesses, route-stage cleanup, and broad driver feature polish
remain deferred.

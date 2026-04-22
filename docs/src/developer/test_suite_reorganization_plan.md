# Test Suite Reorganization Plan

## Purpose

Record the current repo-manager judgment about how the test surface in
`GaussletBases` should be reorganized now that the main internal Cartesian
consolidation line has largely landed.

This note is intended to be durable. It should be updated as the test split
progresses instead of spawning disconnected flat planning notes.

See also:

- [Architecture and current direction](architecture.md)
- [Cartesian source/build unification plan](cartesian_source_build_unification_plan.md)

## Current judgment

`test/runtests.jl` is now large enough that file size itself is a maintenance
problem. The issue is not only line count. It is that too many ownership
domains are mixed into one file:

- radial
- core mapped/cartesian
- nested
- ordinary/QW
- diatomic
- angular
- IDA
- docs/examples/misc

The internal code architecture is now cleaner than the test architecture. That
gap should be reduced.

The right next move is:

- keep `test/runtests.jl` as the runner
- split tests into group-aligned includes
- move only genuinely shared harness helpers into `test/support/`
- trim repeated contract checks only after the structure is clearer

This is now worth doing because the code seams have stabilized:

- normalized operator-build contexts
- source-backed nested front-door context
- nested glass-box contract
- capability-driven validators

The structural runner split has now also landed far enough that the remaining
test-surface work is no longer “how should the file be split?” but “which
remaining semantic checks are stale, duplicated, or badly placed?”

## Current state

As of 2026-04-22:

- `test/runtests.jl` is down to roughly 3.1k lines
- the runner has explicit group gating through `GAUSSLETBASES_TEST_GROUPS`
- extracted group runners now exist for:
  - `radial`
  - `core`
  - `nested`
  - `ordinary`
  - `diatomic`
  - `angular`
  - `ida`
  - `misc`
  - `docs`
- the remaining named inline tail is `examples`
- recent consolidation work created opportunities to remove repeated
  contract-shape assertions and stale snapshot-style checks

## Target shape

The target shape is:

1. `test/runtests.jl` remains the single entry runner
2. domain tests move into included files
3. only real shared helpers move into `test/support/`
4. repeated contract checks are trimmed carefully after the split

Suggested directory shape:

- `test/runtests.jl`
- `test/support/`
- `test/radial/`
- `test/core/`
- `test/nested/`
- `test/ordinary/`
- `test/diatomic/`
- `test/angular/`
- `test/ida/`
- `test/docs/`
- `test/examples/`

This should be treated as ownership cleanup, not as a test-framework rewrite.

## Reorganization principles

- preserve current test-group behavior
- preserve current test names unless a rename adds real clarity
- do not hide important fixture logic behind too much abstraction
- prefer local helpers inside a domain file when they are not genuinely shared
- keep legacy/internal route tests clearly quarantined instead of mixed through
  active public-contract coverage
- do not combine unrelated domains just to reduce file count

## Planned phases

### Phase 1: runner split and first low-risk extraction

Status: done

Scope:

- keep `test/runtests.jl` as the single runner
- extract only the global harness pieces that really need to be shared
- introduce `include(...)`-based domain files
- move the lowest-risk stable domains first

Recommended first extraction target:

- `radial`
- `core`

Why these first:

- they are already relatively self-contained
- they are not the most entangled part of the recent Cartesian refactors
- they validate the include-based structure with lower scientific churn risk

### Phase 2: extract Cartesian contract-heavy domains

Status: done

Scope:

- `nested`
- `ordinary`
- `diatomic`

Goal:

- make the recent source/build/operator contract work easier to review and
  maintain
- keep fixtures and contract checks closer to the domains they protect

### Phase 3: extract angular / IDA / smoke surfaces

Status: mostly done

Scope:

- `angular`
- `ida`
- `docs`
- `examples`
- remaining misc/repl/export checks as appropriate

Goal:

- separate algebra/contract failures from smoke/export/docs failures

Current reading:

- `angular`, `ida`, `misc`, and `docs` are already extracted
- `examples` remains inline behind its own gate and can now be treated as
  optional final polish rather than a blocking architecture problem

### Phase 4: trim repeated contract checks

Status: active next step

Only after the split is stable:

- reduce repeated assertions that merely restate the same unified contract
- prefer one canonical helper where several blocks now test the same internal
  context/contract shape
- keep route-specific behavior checks where they still carry unique value
- remove stale snapshot-style checks that encode an older repo story rather
  than a durable current contract

## Next bounded chunk

The current bounded chunk should be:

- semantic cleanup of the extracted `docs` test contract
- narrow follow-on cleanup of known stale expectations such as
  `working_box`-based test assumptions where they no longer match the live
  object contract

Do not reopen the runner structure unless the optional `examples` extraction
later proves worth doing.

## Main risks

- over-abstracting test helpers so the tests become harder to read
- moving too many domains at once and making failures difficult to localize
- hiding legacy/internal tests rather than cleanly quarantining them
- mixing structural split work with semantic test rewrites
- trimming duplicated checks before the new file ownership is stable

## Non-goals

This reorganization effort does not mean:

- rewriting the scientific meaning of tests
- replacing the current group selection interface
- building a heavyweight custom test framework
- deleting legacy/internal checks just because they are legacy

## Current status

As of 2026-04-22:

- the structural runner split is effectively complete for the named domain
  groups
- the docs-consistency test has already been refreshed away from stale
  snapshot-style string inventories
- the main remaining work on this line is semantic cleanup and trimming, not
  file-structure surgery
- optional final polish remains possible for the still-inline `examples` tail

## Update rule

When a test reorganization chunk lands, update this note rather than opening a
new disconnected planning note unless the topic has clearly split into a
different line.

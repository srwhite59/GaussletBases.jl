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
- keep the fail-closed Cartesian authority checker as its own required CI job;
  docs tests should verify permanent command wiring and markers, not duplicate
  prose inventories or transition-era filenames
- keep legacy/internal route tests clearly quarantined instead of mixed through
  active public-contract coverage
- do not combine unrelated domains just to reduce file count

## Structured-state rule

Do not propagate a new route concept through the driver by adding many scalar
fields to every staged `NamedTuple`.

If a concept has internal structure, make it a compact object owned by the
appropriate module, carry that object through the stages, and expose only a
small summary or fingerprint. Final human-facing reports may expose a few
scalar aliases derived from the object, but compatibility aliases should be
temporary and minimal.

Bad pattern:

```julia
foo_available
foo_status
foo_count
foo_keys
foo_kinds
foo_kind_counts
foo_inventory
foo_summary
foo_materialized
```

especially when the same field cloud is repeated through units, transforms,
pairs, assembly, and reports.

Good pattern:

```julia
foo = FooModule.plan_or_summary(...)
foo_summary = FooModule.summary(foo)
```

Implementation rules:

- before adding more than three related fields to a staged object, stop and
  define a module-owned object or compact summary
- before copying the same field group across two stages, stop and carry the
  object instead
- repeated scalar pass-through fields are a code smell and should trigger
  refactoring

This applies broadly to shellification, lowering, selected lowering, CRC
sidecars, final retained units, pair inventories, operator plans, reports, and
future route concepts.

## Staged metadata assertion rule

Do not compare large staged metadata objects with `==` or `===`.

This is a hard testing rule for route metadata, CRC sidecars, staged summaries,
route inventories, and deeply nested `NamedTuple` objects. Tests must not assert
that a huge sidecar object is identical across stages, and they must not rely on
failure paths that print or type-infer the entire staged object.

Instead, test compact, stable summaries:

- `status`
- counts
- keys, roles, and kind tuples
- booleans
- materialization flags
- short missing-reason tuples

If a stage carries a sidecar forward, tests should prove the contract through
these small fields. Whole-object equality on staged metadata is brittle,
expensive, and can dominate runtime through specialization, deep traversal, or
failure rendering.

Tests for structured staged state should compare compact summaries or
fingerprints, not whole nested objects and not long rows of copied scalar
fields.

## Test runtime policy

Every test should have a runtime class, and long tests require justification.
The goal is to stop the pattern of repeatedly running full route construction
after each mechanical field-carry pass.

Runtime classes:

1. Tiny contract tests

   Goal: seconds.

   Use for every small pass when a pure helper, constructor, fingerprint, parser,
   or local contract is the thing being edited. These tests should not build the
   full driver route unless that route build is itself the contract.

2. Stage propagation tests

   Goal: under roughly 30 seconds after compilation.

   Use for staged field movement and route metadata propagation. These tests
   must compare compact fingerprints or stable summaries, not full staged
   objects.

3. Integration tests

   Goal: allowed to be slow, but explicitly marked.

   Use at baton boundaries, before merging, or when behavior crosses a real
   stage boundary. Do not use integration tests as the per-pass gate for
   mechanical metadata propagation.

4. Long tests

   Anything expected to exceed roughly 2 minutes must carry:

   - a named reason
   - the feature it validates
   - why a smaller test is insufficient
   - a suggested cadence: per-pass, baton-end, nightly, or manual

Hard Codex operating rule:

- before running any test expected to take more than 60 seconds, explain why it
  is necessary
- if a shorter contract, fingerprint, parse, or focused stage test would
  validate the edit, write or run that instead
- do not run full route integration tests repeatedly during mechanical
  field-carry passes
- do not run full integration tests with `--compiled-modules=no`; reserve that
  for parse/load diagnostics only
- time routine Julia tests and probes with Julia-level timing, for example:

  ```julia
  t = @elapsed include("tmp/work/script.jl")
  println("elapsed_s=", t)
  ```

  Avoid routine `/usr/bin/time` wrappers. If OS-level memory data such as
  maximum RSS is genuinely needed, use a stable wrapper script such as
  `tools/time_julia` and approve that wrapper prefix once rather than approving
  broad `/usr/bin/time` usage.

Recommended per-pass cadence:

- documentation-only: `git diff --check`
- syntax-only or dependency-blocked work: parse touched Julia files
- new pure helper/module contract: run the helper's direct test only
- stage field propagation: run a fingerprint/summary propagation test only
- driver behavior change: run one focused staged test
- numerical/materialization change: run the relevant integration test, with the
  expected runtime called out first

Repo-specific classification:

- `test/nested/cartesian_pair_stage_low_order_policy_runtests.jl` is an
  integration gate, not a per-pass gate

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

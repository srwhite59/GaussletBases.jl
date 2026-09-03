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

## Angular Optional-Consumer Control Flow

`HP-ANGULAR-TEST-CONTROL-FN-01` and
`HP-ANGULAR-TEST-CONTROL-TEST-01` authorize one bounded repair to the angular
runner's optional HFDMRG handshake. The current file-level early return is not
an acceptable absence policy: it can stop the included angular owner after
`507` of `61,930` assertions and silently omit the remaining `61,423`.

The repair must preserve these boundaries:

- a missing optional HFDMRG checkout produces one visible skip limited to the
  handshake testset;
- if the configured checkout path exists, loading or importing HFDMRG must
  succeed or fail the test normally; no catch-all may reinterpret an import or
  load failure as absence;
- the angular owner contains no file-level early return, and every test after
  the optional handshake continues to execute;
- `angular` remains an explicitly selectable group and continues to include
  `test/angular/runtests.jl`, but it is removed from `_FAST_TEST_GROUPS` because
  the complete research suite takes roughly 16 minutes;
- all numerical assertions, angular implementations, public APIs,
  dependencies, workflows, and other group definitions remain unchanged.

Implementation is limited to `test/runtests.jl` and
`test/angular/runtests.jl`. The direct control-flow repair should be about
`+6/-4` lines and may add at most `12` control-flow lines. If needed, at most
`12` focused regression lines may be added to the existing
`test/docs/runtests.jl`; no new test file, mock package, configuration or
environment interface, helper abstraction, or framework is authorized.

Acceptance requires an absent-checkout run with a visible handshake skip and
proof that a later angular test executed; a present-checkout run or equivalent
focused evidence showing that import failures propagate; one complete explicit
`angular` run on a supported Julia version; exact confirmation that `fast` no
longer selects `angular`; and package-load, authority, documentation, and diff
checks. The roughly 16-minute complete run is justified once at implementation
acceptance because only the full owner proves that the former early-return
boundary no longer truncates the suite. Fixed-radial extraction, documentation
of angular exports, `ShellLocalAngularProfileKey` de-export, CI ownership, and
all other angular classification remain separate.

Implementation commit `22d25e8c5061abafa5a37ea35d4a51adfa4b9a72`
closed this repair with `+6/-9` lines in the two runner files. Its
absent-checkout run completed in `929.68` seconds with `61,907` passes and one
visible skip, including tests after the handshake. Direct inspection confirms
that an existing checkout now loads without a catch-all, `angular` remains
explicitly selectable, and only its `fast` membership was removed. No focused
docs regression was needed or added.

## Fixed-Radial Angular Public CI Extraction

`HP-ANGULAR-PUBLIC-CI-FN-01` and `HP-ANGULAR-PUBLIC-CI-TEST-01`
authorize one separate extraction of the supported three-level fixed-radial
angular sequence into the Julia `1.10` Supported-floor gate. The new group is
named `angular_public`. It must be available through explicit selection and
the existing `all` selection, but neither `angular_public` nor the complete
`angular` research group belongs to the `fast` alias.

The test change must move and condense, not duplicate, the existing `Atomic
fixed-radial angular sequence` block from `test/angular/runtests.jl` into the
single planned owner
`test/driver_public/angular_fixed_radial_sequence_runtests.jl`. Its compact,
preferably table-driven checks must preserve representative coverage of:

- the `[10, 15, 32]` three-level construction and common radial-basis identity;
- shell identifiers, centers, dimensions, and angular-profile identities;
- both adjacent overlap sidecars and the direct `10 -> 32` sidecar;
- dense level and overlap payload semantics; and
- one representative serialization round trip with its sequence, level,
  profile, gauge, labels, and pair-kind identity metadata.

The focused owner should add `60-80` lines and must not exceed `95`; the whole
superseded block must be deleted in the same commit, and total tracked test LOC
must decrease. `test/runtests.jl` may add only the group symbol and one narrow
include. `.github/workflows/ci.yml` may replace only the Supported-floor group
list with
`core,ida,cartesian,examples,radial,misc,angular_public`.
`test/docs/runtests.jl` may add at most `6` focused lines requiring that exact
selection and rejecting bare `angular`. No source, dependency, fixture
framework, numerical policy, cache, example, workflow row/job, Julia version,
timeout, trigger, permission, or other group may change.

Acceptance requires the focused owner to pass on Julia `1.10` with reported
runtime (about nine seconds expected), and the remote Supported-floor job must
visibly execute it. One combined `angular,angular_public` run on one supported
Julia version is sufficient to prove that extraction did not remove or
duplicate coverage; do not repeat that expensive validation on another Julia
line or during docs-only closeout. The full angular research group remains
outside per-push CI. `all` must retain the new owner and `fast` must exclude
both angular groups. Package load, authority/self-test, docs tests, Documenter,
all three unchanged named CI gates, and diff checks must pass. Documentation of
the seven angular exports and de-export of `ShellLocalAngularProfileKey` remain
separate.

Implementation commit `db156966a4a3a7bf2f685fa0f89312afca7b4280`
closed the extraction. It replaced the `163`-line research-suite block with an
`80`-line, `83`-assertion table-driven public owner, for a total test delta of
`+88/-165` and net reduction of `77` lines. The owner is explicitly selectable,
included by `all`, excluded from `fast`, and is the only angular group added to
the Julia `1.10` Supported-floor selection. Its focused Julia `1.10` run passed
`83/83` in `11.17` seconds. The single authorized combined
`angular,angular_public` Julia `1.12` run passed in `930.34` seconds; it is
accepted evidence and is not repeated during closeout. Remote CI run
`33509253422` passed all three named gates and visibly ran `angular_public`
`83/83`. Full `angular` remains direct-run research coverage outside `fast` and
per-push CI. No source, API, numerical policy, fixture, dependency, job, or
release behavior changed.

## Scheduled Cartesian Internal Maintenance Gate

`HP-CARTESIAN-INTERNAL-MAINTENANCE-CI-FN-01` and
`HP-CARTESIAN-INTERNAL-MAINTENANCE-CI-TEST-01` authorize one bounded internal
numerical-maintenance gate. It is not ordinary CI and grants no public,
release, or scientific-policy status to the tests it executes.

The only approved workflow path is:

```text
.github/workflows/cartesian-internal-maintenance.yml
```

It must use Julia `1.12` and one sequential read-only job with job-level
`contents: read`. The only triggers are manual `workflow_dispatch` and the
weekly cron `17 10 * * 3`. Pushes, pull requests, `main`, and tags must not
trigger it. The job instantiates once, checks package load once, and then runs
exactly these existing files, in this order, as five separate Julia processes:

```text
test/nested/cartesian_atomic_hf_reference_packet_runtests.jl
test/nested/cartesian_occupied_first_injection_runtests.jl
test/nested/cartesian_external_gto_import_runtests.jl
test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl
test/nested/cartesian_screened_hartree_correction_runtests.jl
```

Separate processes are mandatory because the standalone files have overlapping
constants and helper names. Commit `1f550899f` implemented the exact `29`-line
workflow. Its first manual dispatch reached the initial `10`-minute job limit
only during the final suite. Commit `b8e89dae0` then changed only
`timeout-minutes: 10 -> 15`, with a `+1/-1` delta. Manual run `32657600781`
passed in `9m59s`: atomic packet `117/117`, protected sidecar `49/49`, R3A
`464/464` plus facade `64/64`, and screened Hartree `85/85`. The maintained
ceiling is therefore `15` minutes. Any further timeout change requires a new
evidence-backed amendment. Upload no artifact and add no helper script, test
framework, status payload, or coverage dependency.

Pass 592 independently classified the unchanged occupied-first owner as unique,
live physical coverage suitable for this lane. At the relocated Step 4B head it
passed `64/64` in `27.0` seconds, covering packet-driven Be/Ne PQS overlap,
mandatory occupied recovery, optional selection, and terminal due diligence;
the synthetic `misc` checks do not duplicate those endpoints. Commit
`9367afd5ac47b81afe4ad8e2170f274fd36f362f` added it only after the atomic
packet with a two-line workflow step and ten focused policy-test lines. Manual
run `33801669502` passed in `10m17s`: atomic packet `117/117`, occupied-first
`64/64`, protected sidecar `49/49`, R3A `464/464` plus facade `64/64`, and
screened Hartree `85/85`. CI run `33801664774` and Docs run `33801664806`
passed. No owner, command, timeout, trigger, permission, or public CI changed.

The workflow record owns only the existing workflow. The current validation
amendment may add focused policy assertions only in `test/docs/runtests.jl`.
Every numerical file is an unchanged execution input under its existing owner.
Do not edit or split the governance-heavy R3A file or any other numerical owner.

The original Pass 510 screening-deduplication transaction removed only two
demonstrated duplicate classes:

- generic supplied-field malformed-input assertions already protected by the
  `22`-check public owner
  `test/driver_public/screened_hartree_runtests.jl`;
- packet readback/schema failures already protected by the `117`-check atomic
  owner `test/nested/cartesian_atomic_hf_reference_packet_runtests.jl`.

That accepted deletion preserved packet correction, occupied embedding,
additive `P0/q0`, self/cross accounting, translated density fields,
`GG/GA/AA` checks, fitted-field validation, and consumer-specific rejection.
It is maintenance evidence, not authority to edit the screening owner again.

The represented molecular-Hartree suite remains excluded direct-run evidence,
blocked on its incomplete scaling-owner lifecycle. Do not run, reclassify, or
schedule it. Do not edit `test/runtests.jl`, the existing
three-row public CI workflow, production source, APIs, dependencies, manifests,
examples, versions, tags, releases, or immutable RC1 state.

The gate is accepted for weekly and manual maintenance. Each run must retain
the five independent processes, existing order, read-only permissions, and
Julia `1.12`. Contract maintenance requires authority check/self-test, docs
tests, package load, Documenter, manager-log bound, YAML inspection, diff
checks, and the relevant remote workflow evidence. The schedule does not
promote these internal suites into ordinary public CI or publication gates.

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

The five-process scheduled Cartesian internal-maintenance gate is in
maintenance. Represented molecular Hartree remains blocked on its scaling-owner
decision and requires a separate lifecycle assignment.

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

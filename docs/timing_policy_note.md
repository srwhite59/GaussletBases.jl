# Timing Policy Note

This note records the intended repo-wide timing policy for upcoming work.

The immediate trigger was the Cartesian hybrid overlap/transfer path, where a
large slowdown was noticed too late and the code did not already expose which
stage owned the cost.

The purpose of this policy is simple:

- timing should be easy to add
- timing should be on by default
- timing output should stay readable
- expensive phases should be broken down until the real cost center is obvious

## Main Policy

Timing should be treated as a normal observability layer of the repo, not as an
occasional after-the-fact debugging tool.

The policy is:

- collect timing by default
- make it easy to disable globally
- instrument semantic phases, not every tiny helper
- use nested timing scopes
- report only timings above a small threshold
- if a scope is slow and unexplained, add one more timing layer inside it

This is meant to catch large slowdowns early, before wrong performance stories
or unnecessary complexity accumulate around them.

## Global Timing Switch

The first implementation step should be one repo-global timing gate.

Requirements:

- default timing state: `on`
- one simple way to disable timing globally
- available both for normal runs and tests

A likely shape is:

- environment-variable control, for example `GAUSSLETBASES_TIMING=0`
- runtime setter, for example `set_timing!(false)`

The call sites should not need to thread an explicit on/off variable through
the code.

## Macro Surface

The preferred interface is one lightweight macro with `@time`-like ergonomics,
but connected to repo timing state rather than immediate ad hoc printing.

The intended style is:

```julia
@timeg "build S_BA" begin
    ...
end
```

or

```julia
@timeg "cartesian-supplement block" expr
```

The macro should:

- consult the global timing controller
- record a named scope when timing is enabled
- add minimal clutter at the call site

The macro should not require a per-call flag argument. That would make the code
messier immediately.

## Nested Timing

Nested timing should be supported.

The desired behavior is:

- an outer timing scope records total time for a semantic phase
- inner timing scopes record its major parts
- reports show the timing tree, pruned by thresholds

This is better than making each scope responsible only for itself with no
parent-child relation. The tree structure is what makes it easy to see where a
slow phase is spending time.

Implementation-wise, this can be done with:

- a small global timing configuration
- a task-local or thread-local stack of active timing nodes
- a function-level timing helper
- a tiny `@timeg` macro lowering to that helper

The complexity belongs in the timing subsystem, not in the algorithms that use
it.

## Threshold Policy

The first intended thresholds are:

- expansion threshold: `10.0` seconds
- drop threshold: `0.1` seconds

Meaning:

- timings below `0.1 s` should normally be suppressed from reports
- if a scope exceeds `10 s`, its children should be shown
- recursively continue showing children that exceed the drop threshold

This gives a practical recursive rule:

- large phases get broken down
- tiny noise does not dominate the report

The exact numerical thresholds may later change, but this is the right first
policy.

## Instrumentation Rule

The repo should not blindly time every tiny function and inner loop.

That would create:

- code bloat
- timing noise
- misleading reports

Instead, timing should be inserted at semantic boundaries, such as:

- basis construction
- factorized basis extraction
- axis-table construction
- Cartesian-Cartesian overlap block build
- Cartesian-supplement overlap block build
- final contraction
- bundle read/write
- HF iteration phases

If a timed scope is still too large to explain a slowdown, then that scope
should be refined by adding one more layer of timing inside it.

So the recursion is real, but it is driven by evidence rather than by
instrumenting everything from the start.

## Test And Benchmark Practice

Timing should be visible during testing and development, but there is an
important Julia-specific rule:

- do not confuse compilation/warmup time with runtime cost

So timing-oriented tests and benchmarks should usually:

- warm the path once
- then record timing on subsequent execution

The first use of the timing system should emphasize observability, not timing
assertions.

That means:

- timing should be collected and inspectable
- but tests should not initially fail on runtime thresholds

Only later, once stable workflows are well understood, should explicit
performance regression thresholds be added for selected cases.

## Temporary Fine-Grained Timing

During active optimization work, it is acceptable to add narrower temporary
timing scopes inside a suspicious phase.

That temporary instrumentation is useful when:

- a path is known to be slow
- the current timing tree does not identify the cost center

Once the cost center is understood and the path is stable, these finer scopes
may be:

- removed
- collapsed
- or left in place if they continue to pay their way

This is not considered bad code bloat. It is normal investigation machinery.

## Current Example

The immediate motivating example is the Cartesian hybrid overlap/transfer path.

The repo needed to know separately:

- bundle read / reconstruction
- exact cross-overlap construction
- overlap application to occupied columns
- any orthonormalization cleanup
- diagnostics

That path should have exposed its stage timings early.

The policy in this note is meant to make that kind of visibility normal rather
than exceptional.

## First-Pass Design Target

The first implementation should aim for:

- one repo-global timing enable/disable switch
- one lightweight `@timeg` macro
- nested timing tree collection
- thresholded timing report
- semantic instrumentation at major phases

That is enough to start using the policy productively.

It is better to land that cleanly than to overdesign the first version.

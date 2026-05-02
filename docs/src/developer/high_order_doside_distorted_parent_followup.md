# High-Order Doside Distorted-Parent Follow-Up

Date: 2026-04-25

Status: active experimental follow-up note, not a production merge plan

## Purpose

Record the broadened repo-side viewpoint that follows the initial
high-order-doside experimental plan in:

- `docs/src/developer/high_order_doside_experimental_plan.md`

The original plan was the right bounded way to prove the construction, but it
is no longer broad enough as the application plan. The new question is not
whether the undistorted construction is algebraically real. That question has
already been answered well enough to move on.

The new question is whether the lane remains useful in the distorted mapped
regime where one would actually want to use it.

## Why The Scope Must Broaden

The initial sandbox handoff was written from a clean undistorted setting. That
was appropriate for proving:

- full tensor-shell span correctness
- inner-to-outer overlap-metric orthogonalization
- projected one-electron correctness
- ee-only IDA action correctness on a low-noise control case

It was not broad enough to answer the application question, because the real
target route is not the undistorted parent basis. It is a distorted/localized
parent basis where:

- basis quality is no longer a purely clean algebra question
- one-electron performance has to be judged against the currently viable lower
  order route
- ee-only IDA error can become large enough to obscure basis-performance
  differences

That makes the interpretation change explicit:

- the undistorted lane is now the oracle/control lane
- the distorted lane is the application/decision lane

## What The Current Repo Work Already Established

The committed experimental lane already established, on the bounded uniform
undistorted parent-grid regime:

- the high-order full-tensor shell construction is real inside the repo
- the stack basis matches the bounded full-union reference on the tested
  numerical-reference cases
- projected one-electron He+ behavior is sane and variational over shell
  additions
- ee-only IDA He singlet behavior is internally consistent on the same clean
  cases
- spacing, parent-box, conditioning, and moment-risk diagnostics are all
  acceptable on that bounded control

Those are strong results, but they should now be read correctly:

- they validate the construction and its clean control behavior
- they do not yet validate the distorted application route

## The Two Questions That Matter Now

The next experimental phase should treat these as the primary questions.

### 1. He+ performance question

How much better is the new high-order treatment for He+ when measured as:

- accuracy versus retained function count

and compared against:

- the current lower-order viable route
- on the same distorted parent-basis family
- against the same higher-quality distorted-parent He+ reference

This is the real performance question. It should not be phrased first as
"energy quality on the undistorted control", because that is no longer the
decision regime.

### 2. He ee-only IDA interference question

Is the electron-electron IDA error on the distorted He problem small enough
that it does not mask the He+ performance gain of the new basis?

This needs to be stated as an interference criterion, not just a smoothness
criterion.

If:

- `Δ_perf` = the He+ performance gain of the new higher-order basis over the
  lower-order route on a common distorted parent family
- `Δ_IDA` = the distorted-case ee-only IDA error on He over a bounded matched
  comparison

then He can only serve as a meaningful discriminator if:

- `Δ_IDA` is comfortably smaller than `Δ_perf`

If `Δ_IDA` is of the same order as `Δ_perf`, then distorted He with ee-only IDA
is too confounded to decide much about the basis construction.

## Practical Interpretation Of The Current Evidence

The current committed undistorted results do support one useful claim:

- ee-only IDA is not obviously catastrophic on the clean bounded control lane

But they do not yet support the stronger application claim:

- ee-only IDA is small enough on the distorted lane not to interfere with the
  performance judgment

Likewise, the current He+ results show that the higher-order stack behaves
correctly on the clean control lane, but they do not yet provide the decision
metric that matters for the real route:

- He+ accuracy versus retained function count on the distorted parent basis,
  compared directly against the lower-order route

## Revised Experimental Priorities

The follow-up experimental priorities should therefore be:

1. Keep the current undistorted lane as a control/oracle lane.
2. Start a bounded distorted-parent He+ benchmark as the main performance
   study.
3. Run a reduced distorted-parent He ee-only IDA audit to estimate whether the
   He+ performance signal survives the ee approximation.
4. Only then use distorted-parent He as a meaningful basis discriminator.

## Recommended Distorted-Parent Benchmark Shape

The next experimental benchmark should stay bounded, but it should be framed as
an application study rather than another clean algebra study.

Recommended boundaries:

- still experimental only
- still 3D first
- still `doside = 5`
- still odd side ladders first
- same high-order shell construction idea
- same current lower-order route as the baseline comparator
- distorted parent basis allowed and expected
- no public API commitment

Primary benchmark output should be:

- He+ error versus retained function count

for:

- lower-order route
- new high-order doside route

on:

- the same distorted parent-basis family
- the same geometry/spacing family

against:

- a better He+ distorted-parent reference on that same family

## Recommended He ee Audit Shape

The He audit should stay reduced and instrumental.

Its job is not to prove the final physical quality of He on this lane. Its job
is to estimate whether ee-only IDA is accurate enough that He still says
something useful about the basis construction.

Recommended output:

- a bounded distorted-parent He comparison where a better ee reference or
  tighter benchmark is still practical
- an estimate of `Δ_IDA`
- a direct comparison of `Δ_IDA` against the measured He+ performance signal
  `Δ_perf`

Decision rule:

- if `Δ_IDA << Δ_perf`, continue using distorted-parent He as a comparative
  basis-quality probe
- if `Δ_IDA ~ Δ_perf`, then He is too confounded and He+ should remain the main
  performance discriminator for this lane

## What Should Remain Quarantined

Even with the broadened viewpoint, the following should still remain explicitly
experimental:

- all high-order doside lane types, helpers, and drivers
- all distorted-parent extensions of that lane
- all benchmark scripts and comparison harnesses
- any IDA-error interpretation specific to this lane

This still should not widen into:

- public API commitments
- production nested fixed-block/QW integration
- bundle/export contracts
- claims of support for arbitrary `doside`, arbitrary distortions, or general
  distorted nested families

## Recommended Immediate Next Chunk

The next repo chunk should be a bounded distorted-parent He+ performance
benchmark, not more undistorted-lane refinement.

That chunk should answer:

- on a common distorted parent family, does the high-order doside lane buy
  better He+ accuracy per retained function than the current lower-order route?

Only after that should the reduced distorted-parent He ee-only IDA interference
audit be used to decide whether He remains informative enough for this lane.

## Bottom Line

The original undistorted plan was the right first plan, but it was not broad
enough as an application plan.

The current repo evidence justifies a more informed framing:

- undistorted high-order doside is now the control/oracle lane
- distorted-parent high-order doside is now the real application lane
- He+ accuracy versus retained function count is the primary performance metric
- ee-only IDA on He must now be judged by whether it is small enough not to
  interfere with that performance metric

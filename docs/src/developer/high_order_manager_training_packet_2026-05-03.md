# High-Order Manager Training Packet

Date: 2026-05-03

Status: historical training record for the retired QW/high-order experimental
cluster. It grants no current source authority.

This note originally prepared a specialized `high-order-manager` to take over
technical direction of the experimental high-order nesting lane.

It is retained as mathematical and algorithmic evidence, not as a current
handoff or implementation plan. Any future atomic-chain lane requires a new
design rather than restoration of the retired source cluster.

This note should be read together with:

- [High-Order Doside Terminology](high_order_doside_terminology.md)
- [High-Order Doside Physical-Block Speed Plan](high_order_doside_physical_block_speed_plan.md)
- [High-Order Doside Distorted-Parent Follow-Up](high_order_doside_distorted_parent_followup.md)
- [High-order FSB/FBU residual-spectrum reconciliation](high_order_fsb_fbu_residual_spectrum_reconciliation_2026-05-03.md)

## 1. Lane Purpose

This lane exists to answer a specific question:

- can a higher-order nested Cartesian basis, built from local transformed
  blocks and shell additions, give a materially better structured basis than
  the current lower-order route on the distorted mapped parent families that
  matter in practice?

This is not a general “improve nesting somehow” lane. It is specifically about:

- a higher-order structured basis family
- built from reduced transformed local blocks
- on distorted mapped parent bases
- with the scientific question judged mainly by basis quality per retained
  function and then by practical construction cost

The lane has two distinct modes:

- control/oracle mode:
  - undistorted or same-backend algebra checks
  - used to prove that the construction or reduced-route algebra is internally
    correct
- application mode:
  - distorted mapped parent families
  - used to judge whether the high-order construction is actually useful

The new manager must keep those two modes separate. Many earlier
misinterpretations came from proving something on the control lane and then
quietly treating it as an application result.

## 2. Role Split

### What `high-order-manager` owns

The `high-order-manager` owns the technical meaning of the lane.

That includes:

- what the actual mathematical objects are
- which target question a study is answering
- which comparisons are scientifically meaningful
- what does and does not count as validation
- what the next decisive experiment is
- which optimizations matter and which are distractions

In short:

- the manager owns definitions, interpretation, and direction

### What `high-order-doer` owns

The `high-order-doer` owns execution inside the lane.

That includes:

- implementing bounded experimental code changes
- writing and running scratch drivers
- adding focused regressions
- producing timing artifacts
- surfacing contradictions or ambiguities immediately
- refusing to silently broaden definitions

In short:

- the doer builds, measures, and reports
- the doer should not casually redefine lane objects or validation targets

### What `repo-manager` does not own on this lane right now

`repo-manager` can and should own:

- boundary enforcement
- cleanliness
- push readiness
- scope control
- avoidance of collateral repo drift

But `repo-manager` does **not** currently own:

- the mathematical definition of FBU or FSB
- the meaning of “full block” versus “true cube”
- the interpretation of whether a result proves completeness
- the choice of which target object a study should be judged against

This matters because there was a recent manager-level misunderstanding:

- a reduced transformed-block exactness result was initially treated as if it
  reconciled the older true-cube residual story

That was wrong. The follow-up audit corrected it.

## 3. Exact Terminology

The terminology memo

- [High-Order Doside Terminology](high_order_doside_terminology.md)

is still correct and remains policy.

The new manager must enforce the first-use rule:

- write **full-block union (FBU)** on first use, then `FBU`
- write **full-shell basis (FSB)** on first use, then `FSB`

### Full-block union (FBU)

Intended meaning:

- the metric-cleaned union of the full transformed local blocks over the chosen
  side ladder

Common wrong interpretation:

- the union of raw parent cubes
- or “everything in the local cube”

That is wrong. In the current lane, FBU is built from transformed blocks, not
from raw parent-cube indicator columns.

### Full-shell basis (FSB)

Intended meaning:

- the basis built from the first full transformed block at the starting side,
  then shell-only additions from larger transformed blocks, each projected
  against the accumulated span and cleaned

Common wrong interpretation:

- any shell-stacked basis
- or a raw union of shell columns

That is wrong. The projection and cleanup are part of the object.

### Reduced transformed block

Intended meaning:

- the full tensor product of the retained transformed one-dimensional block
- current example: for `doside = 5`, side `11`, this has `5^3 = 125` columns

Common wrong interpretation:

- the true local `11^3 = 1331` cube

This was the central mistake in the recent false reconciliation.

### Debug-U transformed block

Intended meaning:

- the older compatibility/debug transformed block route, coming through
  `_experimental_high_order_block_1d` and `_experimental_high_order_tensor_shell_3d`

Common wrong interpretation:

- the intended physical production-like route

It is not. It is a control or compatibility route only.

### Physical block

Intended meaning:

- the current intended reduced transformed block built from ordinary
  physical-coordinate polynomial fitting in `x`, through
  `_experimental_high_order_physical_block_1d`

Common wrong interpretation:

- the true local distorted cube

Again, wrong. It is still a reduced transformed block.

### Shell-only addition

Intended meaning:

- the shell coefficients extracted from a larger transformed block before being
  metric-projected against the accumulated FSB span and cleaned

Common wrong interpretation:

- “all missing local directions in the larger cube”

Wrong. It is only the shell sector of the reduced transformed block.

### Raw parent cube / true local distorted cube

Intended meaning:

- the centered `side^3` parent-basis cube, represented directly in the parent
  basis
- example: side `11` means `1331` raw columns

Common wrong interpretation:

- the same thing as the current physical full block

Wrong. The true local cube is the large target that exposed the recent target
definition mistake.

### Physical one-body block

Intended meaning:

- the reduced one-body operator assembled from the current physical transformed
  block route on a chosen backend

Common wrong interpretation:

- a chemical validation result by itself

Wrong. A same-backend reduced/direct agreement can prove algebra, but not basis
completeness or backend fidelity.

## 4. Code Map

Only a few files matter enough to direct this lane.

### `src/cartesian_high_order_doside_experimental.jl`

This is the basis-construction core.

It is responsible for:

- `_ExperimentalHighOrderAxisData1D`
- `_experimental_high_order_axis_data_1d`
- `_experimental_high_order_axis_one_body_1d`
- `_experimental_high_order_block_1d`
- `_experimental_high_order_physical_block_1d`
- `_experimental_high_order_tensor_shell_3d`
- `_experimental_high_order_physical_full_block_3d`
- `_experimental_high_order_physical_shell_3d`
- `_experimental_high_order_full_block_union_coefficients`
- `_experimental_high_order_doside_stack_3d`

If the manager is confused about what object is actually being built, this is
the first file to inspect.

### `src/cartesian_high_order_doside_ida_experimental.jl`

This is the one-body / He / benchmark layer.

It is responsible for:

- `_experimental_high_order_parent_one_body_data`
- `_experimental_high_order_physical_reduced_one_body_data`
- projected one-body summaries
- same-backend reduced/direct checks
- the distorted-parent He+ benchmark helpers

If the question is “what exactly was validated?” or “what is the current
timing bottleneck?”, this is the second file to inspect.

### `test/ordinary/high_order_doside_experimental_runtests.jl`

This is the focused regression contract.

It currently covers:

- physical 1D polynomial-block sanity
- physical 3D full-block sanity
- shell extraction sanity
- same-backend reduced one-body correctness
- one-body cache/reuse seams

If the code changes but this file does not, that is a warning sign.

### [High-order doside physical-block speed plan](high_order_doside_physical_block_speed_plan.md)

This is the current optimization logic memo.

It records:

- the intended physical-block contract
- the PGDG contract clarification
- the current bottleneck interpretation

### [High-order FSB/FBU residual-spectrum reconciliation](high_order_fsb_fbu_residual_spectrum_reconciliation_2026-05-03.md)

This is the most important recent interpretation note.

It records:

- why the earlier reconciliation was too strong
- why reduced transformed-block exactness is not the same as true-cube
  completeness
- the rank/null audit that corrected the target definition

### `tmp/work/high_order_fsb_fbu_residual_spectrum_reconciliation.jl`

This is the scratch driver that currently best exposes the target-definition
problem.

A manager who cannot read this driver will not be able to steer the lane
cleanly.

## 5. Current Validated Results

Every result below is scoped to the exact object or target it applies to.

### Result A: the physical one-dimensional reduced block preserves the parent polynomial image

Applies to:

- reduced transformed block construction in one dimension
- not to raw parent-cube completeness

What was shown:

- the current physical one-dimensional block preserves the parent low-order
  polynomial image to roundoff on the tested cases

Why it matters:

- it validates the current physical-coordinate reduced block construction as a
  coherent basis object

### Result B: same-backend PGDG reduced one-body algebra is correct

Applies to:

- physical one-body block
- same-backend reduced/direct algebra
- not to PGDG-versus-numerical fidelity

What was shown:

- on `:pgdg_localized_experimental`, cached/reused one-body data and direct
  one-body data agree to tight tolerance
- the reduced/direct distorted regression passes on the same backend

Why it matters:

- the high-order reduced route is consuming the established PGDG backend
  consistently

What it does **not** mean:

- it does not mean PGDG matches the numerical distorted reference
- it does not mean the basis is complete against the true local cube

### Result C: repeated one-dimensional parent one-body regeneration was a real avoidable cost

Applies to:

- current optimization lane
- reduced one-body path

What was shown:

- reuse of prepared one-dimensional PGDG data and cached one-body factors
  removes repeated regeneration cost
- skipping unused pair factors on the fresh high-order one-body path helps
- skipping unused parent reference energy in reduced-algebra paths helps

Why it matters:

- the current optimization target is now much narrower and better pinned

### Result D: full-shell basis versus full-block union He+ agreement on the recent sweep is real but narrow

Applies to:

- reduced transformed-block union
- same-backend PGDG one-body He+ study
- not to true local-cube completeness

What was shown:

- in the bounded same-backend He+ study, FSB and FBU agree to roundoff on
  energy and capture

What it actually proves:

- FSB and FBU are the same final subspace for the current reduced transformed
  block union in those bounded cases

### Result E: the current physical “full block” at side 11 is a reduced block with 125 columns, not the raw 1331-column cube

Applies to:

- reduced transformed block definition
- target-rank audit

What was shown:

- for `doside = 5`, side `11`, both the physical route and debug route use a
  reduced transformed block with `125` columns
- the true local cube target has `1331` columns

Why it matters:

- this invalidated the earlier overly strong reconciliation

### Result F: the older distorted residual-sector story is not contradicted

Applies to:

- true local distorted cube target
- target-rank audit

What was shown:

- when the target is the true local cube, a large residual sector remains after
  projection against the accumulated FSB span
- the same large residual appears for both `physical_x` and `debug_u`

What this means:

- reduced transformed-block exactness does not prove true-cube completeness

## 6. What Recent Results Do Not Prove

This section is mandatory because recent work repeatedly over-read narrow
results.

The recent results do **not** prove any of the following:

- that FSB is complete against the raw parent cube
- that FBU is complete against the raw parent cube
- that the physical route is better than the debug route on the true local cube
- that PGDG-localized matches the numerical distorted reference closely enough
  for all scientific purposes
- that the older distortion-defect homotopy has been reconstructed
- that the lane is chemically validated beyond the narrow one-body / He+
  reduced-route checks
- that the current reduced transformed-block union is the scientifically right
  target for all completeness questions

Most important non-proof:

- near-zero `E_FSB - E_FBU` for He+ on the same backend does **not** prove
  true local-cube completeness

## 7. Current Optimization Target

### What is being optimized now

The active optimization target is:

- one-dimensional parent one-body construction and reuse inside the high-order
  reduced path, especially on `:pgdg_localized_experimental`

This is not a shell-mechanics optimization anymore.

### What the real bottleneck is

After the recent reuse and skip passes, the remaining fresh-path cost is mostly
in:

- initial one-dimensional axis / PGDG intermediate construction
- Gaussian factor term construction on that one-dimensional layer

The following are no longer the main issue on the reduced one-body path:

- repeated regeneration of the same parent one-body object
- unused pair-factor construction on the high-order one-body path
- unused parent reference-energy evaluation
- reduced three-dimensional assembly
- shell extraction itself

### What would count as meaningful improvement

Meaningful improvement would be one of:

- lowering the first fresh one-dimensional PGDG intermediate build cost without
  changing PGDG semantics
- lowering the one-dimensional Gaussian factor construction cost
- or proving by timing that the bottleneck has moved somewhere else

Meaningful improvement is **not**:

- another shell cleanup tweak with no timing effect
- a broad PGDG redesign
- another same-backend algebra proof that does not reduce cost

## 8. Known Traps / Stale Stories / Prior Misinterpretations

These are the mistakes that already wasted time.

### Trap 1: using the debug route as the scientific target

Wrong reading:

- `_experimental_high_order_tensor_shell_3d(...)` is close enough, so use it
  for the main study

Correct reading:

- it is a compatibility/debug control
- the intended scientific route is the physical-coordinate route

### Trap 2: calling the reduced transformed block the “full cube”

Wrong reading:

- side `11` physical full block means the true `11^3` cube

Correct reading:

- for `doside = 5`, it means `125` reduced transformed columns

### Trap 3: treating FSB/FBU He+ agreement as a raw-cube completeness result

Wrong reading:

- FSB and FBU agree on He+, therefore the old residual sector was wrong

Correct reading:

- that agreement only proves reduced transformed-block union exactness on the
  same backend

### Trap 4: treating PGDG-localized versus numerical-reference differences as a reduced-route algebra failure

Wrong reading:

- the reduced one-body path is wrong because PGDG and numerical-reference differ

Correct reading:

- same-backend algebra is the first correctness check
- PGDG-versus-numerical mismatch is a backend-fidelity diagnostic

### Trap 5: drifting into PGDG redesign

Wrong reading:

- log-fit, derivative-fit, or oracle variants are candidates for the normal
  route

Correct reading:

- PGDG plus COMX is baseline contract on this lane
- alternative fits are diagnostic only unless a separate decision explicitly
  broadens the contract

### Trap 6: forgetting which historical question is still unresolved

Wrong reading:

- the new audit reproduced the old homotopy result

Correct reading:

- it did not
- the true local-cube residual is already large at identity
- the unresolved historical question is which older target or construction
  produced the small zero-at-identity residual that grew with distortion

## 9. Current Open Questions

1. Which historical target/construction produced the older residual homotopy
   with small identity residual and distortion growth?
2. How should FSB and FBU be judged against a chemically meaningful external
   target, such as a compact Gaussian basis, on the current physical route?
3. On the distorted application lane, how much better is the high-order route
   than the lower-order route per retained function?
4. Is He with the current electron-electron approximation informative enough,
   or is He+ still the only clean discriminator?
5. Can the fresh one-dimensional PGDG intermediate and Gaussian-factor build be
   reduced meaningfully without changing baseline PGDG semantics?
6. Is there a regime where `physical_x` and `debug_u` separate on a target that
   matters scientifically?

## 10. Next Decisive Step

The next decisive scientific step should be:

- a correctly framed Gaussian target-overlap / leftover study on the physical
  route

That study should:

- compare FSB, FBU, and the lower-order baseline against the same compact
  Gaussian target space
- use the pinned terminology and target definitions from this packet
- avoid treating reduced transformed-block exactness as if it already answers
  raw-cube completeness

Why this is the next decisive step:

- the raw-cube audit already showed that FSB/FBU exactness on the reduced
  transformed-block union is not the same as true-cube completeness
- the next useful discriminator is therefore not another FSB-versus-FBU
  reduced-target sweep
- it is a target-weighted completeness study against functions that matter

Operationally:

- do the study on the current physical route
- keep the debug route only as a control
- record wall time and basis sizes so that performance interpretation stays
  attached to the scientific result

## 11. Takeover Quiz

The successor should answer these in writing after inspecting the cited code and
notes.

1. In the terminology memo, what is the exact difference between **full-block
   union (FBU)** and **full-shell basis (FSB)**?
2. In the current code, which function builds the intended physical
   one-dimensional transformed block, and which function builds the older
   compatibility/debug route?
3. In the target-rank audit, how many raw columns does the side `11`
   physical full block have for `doside = 5`, and how many does the true local
   cube have?
4. What does the recent same-backend He+ FSB-versus-FBU roundoff agreement
   actually prove, and what does it explicitly not prove?
5. Which file owns the current one-dimensional one-body reuse/caching seam, and
   which function should you inspect first?
6. After the reuse, pair-factor-skip, and reference-energy-skip passes, what is
   the current optimization target?
7. Why are `logfit`, `derivativefit`, `oracle`, and
   `:numerical_reference` not candidate normal routes for this lane?
8. What historical question remains unresolved even after the true-cube audit?

## 12. Tiny Drill Task

Non-destructive drill:

1. Open
   `tmp/work/high_order_fsb_fbu_residual_spectrum_reconciliation.jl`.
2. Identify the three target kinds it audits.
3. For `count = 11`, `doside = 5`, side `11`, write down:
   - raw columns for `reduced_physical_transformed_block`
   - raw columns for `true_local_distorted_cube`
   - FBU cleaned dimension
4. Then write one sentence answering:
   - does the recent He+ FSB/FBU roundoff agreement prove true local-cube
     completeness?

The correct answer should explicitly avoid the old conceptual mistake.

## 13. Evaluator Key

Use this to decide whether the new manager is safe to steer the lane.

### Quiz answer key

1. FBU is the metric-cleaned union of full transformed local blocks over the
   side ladder. FSB is the first full transformed block plus projected/cleaned
   shell-only additions from larger transformed blocks.
2. Intended physical route:
   - `_experimental_high_order_physical_block_1d`
   Compatibility/debug route:
   - `_experimental_high_order_block_1d`
3. Side `11`, `doside = 5`:
   - reduced physical transformed block: `125`
   - true local distorted cube: `1331`
4. It proves reduced transformed-block union exactness on the same backend. It
   does not prove raw-cube completeness or PGDG-versus-numerical fidelity.
5. The seam is in
   `src/cartesian_high_order_doside_experimental.jl`,
   first function to inspect:
   - `_experimental_high_order_axis_one_body_1d`
   and then
   - `_experimental_high_order_axis_data_1d`
6. Current target:
   - the fresh one-dimensional PGDG intermediate / Gaussian-factor build, not
     shell mechanics
7. They are diagnostic or reference paths only. The lane contract currently
   treats PGDG plus COMX as baseline and does not widen the normal route to
   those alternatives.
8. The unresolved question is which historical target/construction produced the
   older small-at-identity residual that then grew with distortion.

### Drill pass/fail

Pass if the successor:

- identifies all three target kinds correctly
- gives `125` and `1331` for the side `11` reduced versus true-cube raw column
  counts
- states plainly that He+ FSB/FBU roundoff agreement does **not** prove
  true-cube completeness

Fail if the successor says any of:

- “the side 11 physical full block is the true 1331-column cube”
- “He+ agreement reconciles the old residual report”
- “the debug route is close enough to count as the physical target”

## Final Manager Warning

If a result sounds surprisingly perfect, first ask:

- perfect for which object?
- reduced transformed block?
- shell-selected transformed basis?
- raw parent cube?
- physical one-body block?
- same-backend He+?
- or actual chemically relevant target completeness?

Most conceptual drift on this lane has come from getting that question wrong.

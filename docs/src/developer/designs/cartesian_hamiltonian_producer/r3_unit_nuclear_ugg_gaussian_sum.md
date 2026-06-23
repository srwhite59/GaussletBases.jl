# R3 Unit-Nuclear U_GG Gaussian-Sum Optimization

Status: approved narrow source authority for reducing R3/RG exact-operator
allocation in terminal final-basis unit-nuclear `U_GG` Gaussian-sum
construction. This is not `HP-R3GG-*`, `HP-CGRB-*`, or route/setup authority.

## Decision

The `HP-R3REM-AUDIT-01` measurement separated the remaining post-`954c86cd`
Cr2 q4 exact augmented-operator allocation. The exact-operator wrapper now
measures about `5.7739s / 4680.627 MiB`; the largest in-wrapper owner is
unit-nuclear `U_GG` factor lookup plus Gaussian-sum construction at about
`2.0447s / 1856.819 MiB`.

Crossed or non-target buckets:

- neutral non-nuclear raw blocks: `0.1894s / 860.736 MiB`;
- neutral nuclear raw blocks: `0.6316s / 15.765 MiB`;
- terminal `G-G` kinetic/moment products: `0.8352s / 733.701 MiB`;
- augmented nuclear transforms only: `0.0125s / 179.268 MiB`.

Approve a narrow source lane for terminal final-basis unit-nuclear `U_GG`
Gaussian-sum construction. Do not use this lane for route/stage setup,
raw-block setup, neutral raw-block kernels, residual Gaussian algorithms,
MWG/IDA, public workflow, or Cr2 artifact/facade work.

## Approved IDs

- `HP-R3UN-FN-01` - terminal final-basis unit-nuclear `U_GG` Gaussian-sum
  allocation reduction.
- `HP-R3UN-TEST-01` - validation gates for the narrow `U_GG` optimization.

## Scope

Approved owner:

```text
Owner module: CartesianFinalBasisRealization
```

Approved source files:

```text
src/cartesian_final_basis_realization/pqs_terminal_one_body.jl
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
```

The first source pass should prefer `pqs_terminal_one_body.jl`. Edits to
`pqs_terminal_residual_gto.jl` are approved only for narrow caller wiring needed
to pass function-local scratch or use the optimized helper from the R3 exact
augmented-operator path.

Approved target functions:

```text
_accumulate_terminal_gaussian_sum!
_terminal_gaussian_sum_action
```

Exact names may follow existing local code, but the approved surface is only
the existing terminal Gaussian-sum path used to construct unit-nuclear `U_GG`
matrices.

## Allowed Implementation Shapes

`HP-R3UN-FN-01` may:

- reuse function-local scratch/workspace across Gaussian-sum terms and center
  calls;
- accumulate terminal Gaussian-sum contributions in-place into the caller's
  destination;
- reduce avoidable allocation in factor lookup and terminal Gaussian-sum
  action construction;
- introduce small internal scratch arguments or file-local helpers if they
  remain inside `CartesianFinalBasisRealization` and do not create persistent
  state;
- simplify or delete obsolete allocation-heavy code inside the targeted
  Gaussian-sum path once parity is established.

The mathematical operator is unchanged: each unit-nuclear `U_GG` block remains
the exact final-basis uncharged nuclear attraction matrix assembled from the
approved term-first Gaussian-sum factors and terminal basis representation.

## Not Approved

This amendment does not approve:

- neutral nuclear or non-nuclear `G-A`/`A-A` raw-block changes;
- terminal `G-G` kinetic, coordinate-moment, or second-moment product changes;
- residual Gaussian selection, orientation, transforms, exact augmented
  transform semantics, MWG, or IDA changes;
- Qiu-White semantic changes or route objects;
- route/stage setup cleanup, raw-block setup cleanup, parent construction, or
  terminal basis realization changes;
- persistent caches, persistent workspace objects, broad Gaussian-sum
  frameworks, provider bundles, metadata/report/status/payload fields,
  artifacts, public API/export, committed tests, Cr2 facade support, or Cr2
  artifact workflow.

## Validation

`HP-R3UN-TEST-01` approves validation only for this narrow optimization:

- existing H2 Residual Gaussian endpoint unchanged;
- Be2 Residual Gaussian facade/readback unchanged except for allowed
  timing/allocation improvement;
- Cr2 q4 exact-operator audit reports before/after unit-nuclear `U_GG`
  allocation and total wrapper allocation;
- Cr2 q4 replay parity for unit-nuclear `U_GG` blocks and final exact
  augmented operators at roundoff;
- exact operators remain finite and symmetric.

Recommended commands for a source pass:

```text
git diff --check
julia --project=. -e 'using GaussletBases; println("load ok")'
julia --project=. test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl
julia --project=. tmp/work/be2_r3u_facade_measurement.jl
julia --project=. tmp/work/cr2_exact_operator_allocation_audit.jl
```

If the Cr2 audit does not directly compare `U_GG` blocks, add an ignored
`tmp/work` parity probe. Do not add a committed test file under this ID.

## Line Budget And Failure Rule

Line budget:

- at most `100` added `src` lines total;
- net simplification is expected through deletion/simplification of
  allocation-heavy Gaussian-sum helper code;
- no new committed test file.

Failure rule: if this needs a persistent cache/workspace object, a broad
Gaussian-sum framework, files outside the approved source files, source edits
outside the terminal unit-nuclear `U_GG` path, or more than the line budget,
stop and request a new docs-only amendment before coding.

## Deferred

Route/stage setup, raw-block setup, neutral raw-block kernels, residual
Gaussian algorithm changes, MWG/IDA, public supplemented workflow/export, Cr2
facade support, and Cr2 artifact workflow remain deferred.

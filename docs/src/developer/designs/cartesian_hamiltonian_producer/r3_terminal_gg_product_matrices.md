# R3 Terminal G-G Product Matrices

Status: approved narrow source authority for reducing R3/RG exact-operator
`G-G` product-matrix allocation. This is not `HP-CGRB-NN-*` raw-block
authority.

## Decision

The neutral Cartesian Gaussian raw-block lanes have crossed their measured
Cr2 q4 bottlenecks. After Pass 086B, neutral non-nuclear raw blocks measure
about `0.1873s / 860.736 MiB`, while the Cr2 exact-operator wrapper is about
`6.4184s / 9043.987 MiB`. The remaining exact-operator allocation is now
reported as dominated by terminal final-basis `G-G` product matrices and
unrelated route/stage setup, not neutral `G-A`/`A-A` raw blocks.

Approve a narrow terminal product-matrix optimization lane owned by
`CartesianFinalBasisRealization`.

## Approved IDs

- `HP-R3GG-FN-01` - R3/RG terminal `G-G` product-matrix construction
  optimization.
- `HP-R3GG-TEST-01` - validation gates for the narrow optimization.

## Scope

`HP-R3GG-FN-01` approves only the `G-G` product matrices used by
`pqs_terminal_residual_gto_augmented_operators(...)`:

- kinetic `K_GG`;
- coordinate moments `x_GG`, `y_GG`, `z_GG`;
- second moments `x2_GG`, `y2_GG`, `z2_GG`.

Approved source owner:

```text
Owner module: CartesianFinalBasisRealization
```

Approved source files:

```text
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
src/cartesian_final_basis_realization/pqs_terminal_one_body.jl
```

The first source pass should prefer `pqs_terminal_residual_gto.jl` only. Edit
`pqs_terminal_one_body.jl` only if a small internal terminal-product workspace
or multi-product helper is needed to avoid repeated function-local buffer
allocation in the existing terminal product assembly. Any helper added there
must remain internal to `CartesianFinalBasisRealization`.

## Allowed Implementation Shapes

Allowed changes:

- accumulate the three kinetic-axis product contributions into one `K_GG`
  destination instead of allocating and summing three full matrices;
- when a same-construction caller already has the base Hamiltonian kinetic
  matrix, reuse that exact `G-G` kinetic block rather than recomputing it,
  provided equality to the existing construction is validated;
- build coordinate and second-moment `G-G` product matrices one axis at a time
  and transform each immediately, instead of retaining larger intermediate
  tuples longer than needed;
- share function-local terminal-product scratch/workspace across consecutive
  product assemblies in the same exact-operator construction;
- delete or simplify `_r3a_product_matrix(...)` if it is replaced by narrower
  accumulation/fill helpers with no remaining live caller.

Any new helper should use domain-local names such as:

```text
_r3a_fill_product_matrix!
_r3a_accumulate_product_matrix!
_r3a_terminal_gg_product_blocks
```

Exact names may follow local style, but they must not introduce route-stage,
status, report, payload, cache, artifact, public API, or raw-block vocabulary.

## Not Approved

This amendment does not approve:

- `G-A` or `A-A` Gaussian raw-block changes;
- nuclear raw-block changes;
- unit-nuclear `U_A` Gaussian-sum changes;
- IDA or MWG interaction changes;
- terminal basis realization changes;
- residual Gaussian selection, orientation, or transform changes;
- Qiu-White semantic changes or Qiu-White route objects;
- parent construction or parent-stage fields;
- persistent caches, broad provider bundles, metadata, report/status/payload
  fields, artifact schema, public API, or exports;
- Cr2 facade support, full Cr2 Hamiltonian workflow, or Cr2 artifact workflow.

## Validation

`HP-R3GG-TEST-01` approves validation only for the narrow G-G optimization:

- existing H2 Residual Gaussian endpoint unchanged;
- Be2 Residual Gaussian usability/performance measurement unchanged except for
  allowed timing/allocation improvement;
- Cr2 q4 `K_GG`, coordinate moment `G-G`, and second-moment `G-G` product
  matrices match the current construction at roundoff, as ignored validation;
- augmented exact operators remain finite and symmetric;
- base `G-G` block equality checks in the existing H2 endpoint still pass;
- Cr2 q4 exact-operator allocation is remeasured after parity.

Recommended commands for a source pass:

```text
git diff --check
julia --project=. -e 'using GaussletBases; println("load ok")'
julia --project=. test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl
julia --project=. tmp/work/be2_r3u_facade_measurement.jl
julia --project=. tmp/work/cr2_exact_operator_allocation_audit.jl
```

If the Cr2 audit does not already compare `G-G` product matrices directly, add
an ignored `tmp/work` parity probe for `K_GG`, `x/y/z`, and `x2/y2/z2`.

## Line Budget And Shrinkage

Line budget:

- at most `100` added `src` lines total;
- no new committed test file in the first source pass;
- net source simplification is expected through deletion or simplification of
  `_r3a_product_matrix(...)` or equivalent repeated allocation paths.

Failure rule: if the optimization needs a persistent workspace/cache object,
a broad product-operator framework, files outside the approved source files, a
new public/internal payload, or more than the line budget, stop and request a
new docs-only amendment before coding further.

## Deferred

Route/stage setup allocation remains a separate measured cost. Do not use this
lane to optimize route setup, base producer stages, public workflow, artifacts,
or Cr2 facade readiness.

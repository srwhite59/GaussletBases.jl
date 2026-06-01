# PQS Source-Box Operator Framework

This note is the framework contract for the projected q-shell (PQS) operator
lane. It exists to prevent architecture drift across multi-pass work,
especially after handoffs, reviews, or context compaction.

This is not a public API contract and not a claim that the production route
already consumes this machinery. It is the current design frame for private
PQS operator work and retained-unit unification.

## Purpose

PQS should provide a compact, source-box-first way to build high-order
Cartesian retained spaces and operator blocks.

The practical goal is:

```text
small raw product boxes
-> one-dimensional operator factors
-> retained transforms / realization transforms
-> retained operator blocks
```

The lane should move PQS toward a usable route for molecular calculations
without turning shell-row support contractions into the algorithm. Shell-row
representations are useful for construction checks, current-route authority
comparisons, and debug oracles. They are not the desired computational
primitive for PQS operator assembly.

## Source Documents

Future work in this lane should read this page first, then the detailed policy
notes:

- `docs/src/developer/raw_product_source_retained_transform_policy.md`
- `docs/src/developer/projected_q_shell_policy.md`
- `docs/src/developer/performance_review_contracts.md`

If those documents and this one disagree, stop and update the framework or
policy explicitly before continuing implementation. Do not resolve the conflict
silently in code.

## Key Spaces

The lane uses four spaces that must remain distinct.

### Parent Lattice

The parent lattice is the full Cartesian gausslet grid. It owns the physical
axis points, parent indexing, PGDG operator data, and any parent-level
coordinate/order contract.

Parent rows are not the preferred loop space for PQS operator blocks when a
smaller raw product source exists.

### Raw Product Source Box

A raw product source box is a small local product domain:

```text
axis intervals + source_mode_dims = (nx, ny, nz)
-> 1D source transforms
-> source product modes
```

Once the axis intervals and total source-mode dimensions are chosen, the 1D
source transforms should be deterministic under the selected rule. For the
current PQS lane, these are COMX/source transforms over the full source box.

Operator factors belong here first. For kinetic:

```text
K_raw =
    Kx x Sy x Sz
  + Sx x Ky x Sz
  + Sx x Sy x Kz
```

where each factor is a 1D cross-axis factor between the left and right source
boxes.

### Retained Transform

A retained transform maps raw source modes to retained columns. In the compact
mode-selected PQS reference path, this can be boundary COMX-product mode
selection.

For product/doside units, this is the product/doside retained transform or an
identity selector over a source box.

The retained transform is the object used in the operator formula:

```text
O_retained = T_left' * O_raw_product * T_right
```

### Shell Realization

Shell realization is separate from raw product-box operator construction. For
the current shell-supported PQS fixture, realization is:

```text
selected product-box modes
-> project to shell rows
-> full-rank symmetric Lowdin cleanup
-> isometric shell-supported retained representation
```

If the active current-route object is shell-realized, the realization transform
must be represented as part of the retained transform from source-box space to
retained columns. It must not force the operator algorithm to loop over shell
rows except as an oracle.

## Intended Data Flow

The intended operator data flow is:

```text
left raw product source box
right raw product source box
-> 1D cross-axis factors
-> raw product-box pair operator
-> left/right retained or realization transforms
-> retained operator block
```

For any pair of retained units that have raw product source boxes:

```text
O_block = T_left' * O_raw(left_source, right_source) * T_right
```

This is the organizing principle for:

- PQS/PQS blocks;
- PQS/product blocks;
- product/product blocks;
- future local/Gaussian one-body source-box blocks;
- future GTO/source-box transfer helpers where applicable.

Support-local contraction may still be used to verify the result:

```text
support rows + explicit retained coefficients -> oracle block
```

but this support-row path must be labeled as oracle/debug/reference. It is not
the algorithmic target for PQS.

## Required Objects

Future implementation should keep these objects conceptually separate even if
the initial code uses private named tuples or adapters.

### Raw Product Box Plan

Owns:

- parent axis intervals;
- total source-mode dimensions;
- 1D source transforms and local coefficients;
- source-mode ordering;
- PGDG/no-fallback provenance diagnostics where available;
- no retained rule.

Must not own:

- shell-row support coefficients;
- Lowdin cleanup;
- packet construction;
- retained weight semantics.

### Retained Rule Or Transform

Owns:

- source-box identity;
- retained dimension;
- transform from source modes to retained columns;
- rule kind, such as boundary selector, product/doside transform, identity
  selector, or shell-realization transform;
- diagnostics describing whether the transform is intended algorithmic input
  or oracle-only.

Must not pretend retained PQS weights are positive quadrature weights.

### Retained Unit

Owns:

- column range;
- role and route provenance;
- support ownership metadata;
- retained transform metadata;
- capability labels for operator families.

The unit may carry support-local coefficients for debug, current-route
compatibility, or authority comparison. Those coefficients do not by
themselves define the intended PQS operator algorithm.

### Pair Plan

Owns:

- left/right retained units;
- left/right source boxes, when available;
- supported terms;
- 1D cross-factor construction plan;
- selected algorithm path;
- oracle path, if present;
- cost and provenance diagnostics.

The pair plan is the place to say whether a pair has a source-box algorithm
available, is oracle-only, or is unsupported.

## Invariants

These invariants are lane-wide.

- Operators for product-structured retained units should be built in raw
  source-box spaces first.
- Shell-row projection and Lowdin cleanup belong to realization, not raw-box
  operator construction.
- The support-local shell-row path is an oracle/debug/reference path unless a
  framework update explicitly says otherwise.
- Retained PQS weights are debug/reference metadata only. They must not become
  IDA division weights or positive quadrature-weight carriers.
- PGDG analytic provenance and "no numerical fallback invoked by this helper"
  are different claims and should be reported separately.
- Current fixed-block or packet authority remains authoritative until a private
  source-box path has exact checks on the relevant fixtures.
- Private source-box helpers do not imply public/default route adoption.
- A full retained matrix diagnostic is allowed for small fixtures, but it is
  not the intended scaling pattern.

## Must-Not-Become Patterns

The lane should stop if it starts drifting into any of these patterns.

- "Optimize support-local fallback" as the PQS operator strategy.
- Labeling shell-row contractions as the active PQS algorithm.
- Treating Lowdin cleanup alone as the full source-to-retained transform.
- Mixing shell-row realization into raw product-box operator construction.
- Claiming PQS/product or PQS/PQS production support from oracle-only evidence.
- Treating retained PQS weights as IDA-safe weights.
- Adding local/ECP/Gaussian/MWG/interaction work before the safe source-box
  operator path is structurally clean.
- Installing sidecars, packets, QW/Hamiltonian hooks, or public/default routes
  from a diagnostic helper.
- Keeping stale docs or helper names that encode the wrong contract.

## Diagnostic Versus Algorithmic Paths

Use these labels consistently.

Algorithmic source-box path:

- builds 1D source-box factors;
- forms or streams raw product-box pair operators;
- applies retained transforms;
- avoids shell-row loops except for verification.

Shell-realization transform path:

- starts from source-box modes;
- projects selected modes to shell rows;
- applies Lowdin cleanup;
- returns an isometric retained transform or equivalent retained columns;
- can be used in `T_left' * O_raw * T_right`.

Support-local oracle path:

- uses explicit support-row coefficients;
- contracts rows directly;
- validates current-route semantics;
- may be slow;
- must not be described as the route algorithm.

Current-route authority comparison:

- compares private helper output to fixed block or sequence packet fields where
  those fields exist;
- proves agreement for those fields and fixtures only;
- does not imply adoption.

The recent `_pqs_current_route_*` helpers should be read through this lens.
In particular, `_pqs_current_route_retained_pair_policy(...)`,
`_pqs_current_route_safe_term_matrices(...)`, and
`_pqs_current_route_safe_term_authority_comparison(...)` are validation and
oracle scaffolding for the current shell-realized route. They are not the
algorithmic policy for PQS shell pairs, and their support-local shell-row
contracts must not be copied into source-box implementation as the target
algorithm.

## Current Evidence

The current evidence supports the framework direction, but it is still private
and partial.

What is established:

- Raw product-box PQS self references exist for safe terms.
- PQS/product source-box mixed blocks exist for overlap, position, `x2`, and
  kinetic, including nontrivial product retained transforms.
- Product/product source-box shadow blocks exist, showing the vocabulary is not
  PQS-only.
- Cross-PQS source-box references exist for compatible raw-box fixtures.
- A private q4 current-route retained-unit inventory covers all retained
  columns and has a route-wide safe-term authority comparison.
- A Be2-like strict PQS q5 inventory-shape check pins the 8-unit, 36-pair,
  multi-shared-shell route shape.
- A blockwise timing probe shows that shell-realized PQS/PQS support-local
  oracle contraction is the dominant cost center.

What is not established:

- Production packet or fixed-block adoption of the source-box algorithm.
- A complete shell-realized PQS/PQS source-box realization-transform block for
  the current route.
- Local/Gaussian one-body terms in the source-box framework.
- MWG/IDA interaction support through retained PQS source-box transforms.
- CR2 validation or public/default route readiness.
- That dense retained full-matrix diagnostics scale to larger routes.

## Performance Expectations

The desired cost model is controlled by source-box dimensions and retained
dimensions, not by shell support-row counts.

For PQS/PQS safe terms, the expected algorithmic shape is:

```text
1D source-box factors
-> small source-box operator action
-> small retained transforms
-> retained block
```

The observed Be2-like timing shows why this matters. A sampled shell-realized
PQS/PQS support-local oracle block of shape `(98, 98)` took tens of seconds
for overlap/kinetic even though allocation was small. That is evidence that
the oracle path is CPU-expensive, not evidence that support-row contraction
should become the algorithm.

Future performance reports should separate:

- inventory construction time;
- source-box factor time;
- retained-transform application time;
- oracle comparison time;
- full retained matrix materialization time, when used;
- allocation and dense matrix materialization claims.

## Open Questions

These questions are intentionally unresolved.

- What is the exact private representation of a shell-realization transform
  suitable for `T_left' * O_raw * T_right`?
- Should the first shell-realized PQS/PQS implementation materialize the raw
  source-box pair operator or stream separable factors directly through the
  transforms?
- What is the smallest Be2-like sampled validation that proves same-shell and
  cross-shell PQS/PQS blocks without running a full route-wide oracle?
- How should support-dense atom-box fallback enter the broader retained-unit
  vocabulary without pretending to be product/doside?
- When should local/Gaussian one-body terms be introduced after safe terms?
- What evidence is required before any packet/fixed-block adoption discussion?

## Framework Update Rule

Update this document before implementation when any of these changes are
proposed:

- changing the meaning of PQS retained transforms;
- treating shell-row support-local contraction as more than an oracle;
- adding a new retained-unit kind to the PQS route vocabulary;
- adding a new operator family such as local/Gaussian/MWG/IDA;
- moving from diagnostic helpers to construction, packet, fixed-block, QW, or
  public/default adoption;
- changing weight semantics;
- changing performance expectations or validation gates.

Every PQS baton loop should cite this document in `RUN.md` and in each blurb.
Each architecture-sensitive PQS baton loop should also create or refresh a
loop-local `INVARIANTS.md` that quotes the governing invariants from this
framework. Before each pass and before each commit, reread this framework and
the loop-local invariants file. If implementation pressure suggests the
framework is wrong or incomplete, stop and make the framework update explicit.

## Next Intended Correction

The next source-box algorithmic correction should not optimize support-local
fallback. It should implement a private shell-realized PQS/PQS safe-term block
that:

1. builds overlap and kinetic in raw product-box space from 1D factors;
2. applies exact shell-realization transforms to retained columns;
3. compares against support-local shell-row contraction only as an oracle;
4. starts with same-shell fixtures;
5. then extends to cross-shell fixtures if the same-shell result is exact and
   inexpensive.

No packet adoption, fixed-block construction adoption, QW/Hamiltonian routing,
local/ECP/Gaussian/MWG/interaction work, public/default route, or CR2 claim is
part of that correction.

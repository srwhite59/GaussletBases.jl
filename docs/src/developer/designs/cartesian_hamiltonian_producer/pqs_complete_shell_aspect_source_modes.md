# Matched PQS/White-Lindsey Complete-Shell Aspect Source Modes

Status: implemented under `HP-PQS-ASPECTSHELL-FN-01`; committed validation is
completed under `HP-PQS-ASPECTSHELL-TEST-01`. Both records are in maintenance.

This page is the canonical contract for selecting one aspect-aware outer
source shape for a shared complete shell and consuming it in matched
bond-aligned diatomic PQS and White-Lindsey constructions. Machine authority
owns source/test permission and lifecycle; this page owns the source-shape
algorithm, family-specific retention, validation, and exclusions.

## Purpose

A rectangular physical shell must not silently receive a cubic source span
merely because the public transverse size is `ns`. It also must not use an
aspect-aware longitudinal size in PQS while White-Lindsey independently
retains the cubic `ns-2` count on every axis. The latter compares different
shell dimensions rather than two retained constructions of one source shell.

For an eligible shared z-axis complete shell, select one outer shape:

```text
source_mode_shape = (ns, ns, L)
```

where `L` comes from the existing angular-resolution selector rather than raw
index aspect. PQS and White-Lindsey consume that same shape but construct
different retained subspaces of the same dimension.

## Construction Boundary

The bond-axis source length cannot be selected honestly in the early terminal
region builder, which sees region metadata but not all parent/bundle and
retention facts. It also cannot be repaired downstream after retained and
support records have frozen a shape.

The owner is the terminal low-order route path after:

```text
shellification
+ parent/bundle construction
```

and before:

```text
lowering contract inventory
retained-unit plans
retained-unit transform contracts
support records
terminal retained-rule plans
due-diligence rows
final realization
```

`_pqs_source_box_route_driver_aspect_shell_lowering_plan(...)` currently
enriches PQS `:complete_shell` contracts. The approved correction generalizes
that same private operation to the eligible White-Lindsey complete-shell
contract. It must not create a second selector, shellifier, or route-global
plan.

Other region kinds, one-center construction, and nonshared shells retain their
existing policies. White-Lindsey's scalar route-local count remains a live
fallback for those ineligible cases.

## Angular-Resolution Selection

For each eligible complete shell, the route owner calls the established
diatomic source-box dimension planner with:

- the parent axis bundle and atomic locations;
- the shell outer and inner-exclusion boxes;
- selected transverse outer count `ns`;
- bond axis;
- current complete-shell retention policy and support count;
- the established shared-shell angular-resolution scale.

The selector returns:

```text
source_mode_dims
raw_source_dims
selected_q
raw_q
raw_L
axis_selector_retained_counts
```

For bond axis `z`, `raw_q` is the common transverse outer source size and
`raw_L` is the selected longitudinal size. The enriched lowering contract
records one authoritative shape and policy:

```text
source_mode_shape = source_mode_dims
source_mode_policy = :diatomic_shared_shell_adaptive_angular_source_box
```

A simple physical-aspect estimate may remain a due-diligence comparison. It is
not the construction rule and must not overwrite the angular-band result.

## Family-Specific Retention

Let the common outer shape be:

```text
D = (Dx, Dy, Dz) = (ns, ns, L)
```

PQS retains its existing boundary quotient of the outer product span:

```text
N_PQS = prod(D) - prod(D .- 2)
```

White-Lindsey splits the same shell support into native facets, edges, and
corners. Its one-dimensional contraction count on a free axis is the
corresponding inner count:

```text
K = D .- 2 = (ns - 2, ns - 2, L - 2)
```

The sum of all White-Lindsey stratum products must satisfy:

```text
N_WL = prod(D) - prod(K) = N_PQS
```

This is dimension parity, not subspace identity. PQS projection and
White-Lindsey product contraction remain the scientific objects being
compared. Public `ns` and route-local `q = ns-2` keep their existing meanings;
the WL route-local scalar is the transverse part of `K`, not a second outer
source shape.

## Downstream Shape Consistency

Every downstream representation of an eligible shell must consume the same
`source_mode_shape`:

- lowering contract metadata;
- retained-unit and support metadata;
- retained transform contracts;
- due-diligence actual-shape rows;
- final realization validation.

The current PQS terminal route carries noncubic `source_mode_shape` through
its lowering, support, retained-rule, and terminal-realization records. It
must require equal transverse dimensions, derive transverse `q` and bond-axis
`L`, and must not replace `L` with `q` after an authoritative shape exists.

White-Lindsey carries the existing `source_mode_shape` to its stratum children
and derives each free-axis count from the corresponding component minus two.
Do not add another shape field or encode the longitudinal count as a new
scalar. The existing `white_lindsey_retained_count_1d` remains only the
fallback for shells to which this matched policy does not apply.

This authorizes one narrow legacy-stage metadata read of the already-carried
`source_mode_shape` in the WL realizer. It does not authorize another metadata
key, a generic metadata-driven algorithm, or additional downstream metadata
reads.

Due diligence may repeat the common shape on native WL stratum rows, but
full-shell retained-count validation must aggregate all strata with the same
physical shell identity. It must not compare one facet or edge count with the
full-shell boundary count.

## Source Ownership

Current maintenance ownership is limited to:

- `src/pqs_source_box_route_driver_helpers.jl` for angular-dimension
  selection and matched lowering-contract enrichment;
- `src/cartesian_base_hamiltonian.jl` for due-diligence shape reporting and
  route-aware aggregate interpretation;
- `src/cartesian_final_basis_realization/white_lindsey_terminal_basis_realization.jl`
  for axis-specific WL one-dimensional counts.

The obsolete contracted-parent/multilayer planning route was retired under
`HP-RETIRE-CONTRACTED-PARENT-FN-01`. It is not an aspect-shell owner or
fallback implementation.

The older angular selector remains in its established diatomic owner. No
second selector or user-facing aspect parameter is introduced.

## Implementation Budget

- At most `60` added source lines across the approved existing files.
- At most `30` added test lines in
  `test/driver_public/cartesian_base_hamiltonian_runtests.jl`.
- No new source, helper, test, probe, fixture, module, driver, or tool file.
- No new struct, persistent object, metadata key, report field, artifact key,
  or compatibility adapter.
- `bin/pqs_paper_h2_driver.jl` is replayed unchanged.

If the correction cannot fit these limits while preserving the existing
family boundaries, make no source commit and report the missing reusable
operation.

## Evidence And Validation

The implemented PQS replay replaced hidden cubic `(5,5,5)` shapes with
angular-band selections such as:

```text
(5,5,8), (5,5,7), (5,5,6), (5,5,6)
```

The 2026-07-29 paper audit then found that the standard H2+/H2 comparison used
these PQS shapes but left White-Lindsey at `(3,3,3)`. At the accepted `R=2`
point, common core/slab counts were `275/50`, while complete-shell totals were
`960` for PQS and `784` for WL, producing bare dimensions `1285/1109`. Those
values remain valid for the two implemented bases but are provisional evidence
for a matched-method comparison.

Commit `6e2c97704` implemented matched consumption. Its clean unchanged-driver
replay preserved the parent and both PQS rows exactly, changed WL shared-shell
columns from `784` to `960`, and produced matched bare/supplemented dimensions
`1285/1303` for both families. The rebuilt WL fixed-target fingerprints are
`80b446467a5d0529d6e428e7fa54b134ed172f638cfb07bb7098efd71a0cb5d9`
and
`0d749a687a87f7b8beb8453095ae6301d8f3c74ea958aaac8db0ee9f89744898`.
These are production measurements; external paper-oracle interpretation
remains a consumer validation step.

The correction acceptance checks include:

- unchanged parent fingerprint, shell regions, owned supports, direct cores,
  slabs, public `ns`, and route-local `q`;
- exact common outer-shape agreement between PQS and WL for every eligible
  shared shell;
- PQS retained count and aggregated WL stratum count both equal to
  `prod(D)-prod(D.-2)`;
- equal direct-core, complete-shell, slab, and bare final dimensions for the
  matched construction;
- exact PQS dimension, coefficient, H1, and interaction parity with the
  pre-correction implementation;
- exact shape agreement across lowering, retained, support, due-diligence, and
  final realization records;
- block-local WL Gram identity and exact terminal column accounting;
- finite and symmetric base one-body and interaction matrices;
- no WL missing-shape or per-child/full-shell warning mismatch.

Validation uses one bounded committed H2 regression in the existing public
Cartesian base test and a clean replay of the existing private paper driver.
The paper replay must preserve the parent and PQS rows, rebuild bare and
supplemented WL operators and interactions from the corrected terminal basis,
and rerun any external same-density oracle against the new WL state
fingerprints. No old WL energy or fingerprint is an acceptance target.

## Failure Behavior

Stop if:

- the angular selector cannot produce finite, valid dimensions;
- any component of `source_mode_shape .- 2` is nonpositive;
- transverse dimensions disagree;
- PQS and eligible WL contracts carry different outer shapes;
- aggregated WL stratum count differs from the PQS boundary count;
- support ownership or shell-local projection would need to change;
- a noncubic WL count cannot reach final realization through the existing
  stratum products;
- implementation requires a new route object, metadata field, persistent
  carrier, driver branch, or broad payload.

Do not repair failure by restoring cubic shapes, guessing `L` from index
aspect, weakening support checks, changing retained counts after the contract
is frozen, or teaching downstream stages to patch conflicting shapes.

## Explicit Non-Goals

This contract does not approve:

- public driver/API inputs or automatic source-shape tuning;
- shell ownership, direct/core identity, thin-slab, or angular-z-extension
  changes;
- one-center or nonshared White-Lindsey retention changes;
- route-local `q`, supplement selection, residual cutoff, or interaction
  formula changes;
- artifact schema or provenance changes;
- injection, EGOI, screened-Hartree, solver, or Cr2 production work;
- route-global materialization revival or a report/payload framework.

Terminal due diligence remains advisory reporting. It exposes actual and
expected physical-aspect shapes but does not itself choose `L`.

The separately implemented
[semantic shell-q override](pqs_semantic_shell_q_overrides.md) may rerun this
selector for a private PQS-only refinement or coarsening. This amendment does
not add a matched White-Lindsey override or change that lane's rejection of WL
requests.

# PQS Complete-Shell Aspect-Aware Source Modes

Status: implemented internal construction policy under
`HP-PQS-ASPECTSHELL-FN-01` and
`HP-PQS-ASPECTSHELL-TEST-01`.

This page is the canonical contract for aspect-aware source dimensions of
shared complete shells in bond-aligned diatomic PQS construction. The registry
owns source/test permission and lifecycle; this page owns the source-shape
algorithm, stage boundary, validation, and exclusions.

## Purpose

A rectangular physical shell must not silently receive a cubic
`(q,q,q)` source span merely because `q` is the transverse PQS size. That
under-resolves bond-axis angular content and changes retained counts,
transforms, Hamiltonian matrices, and downstream energies.

Terminal due diligence first exposed this failure. The implemented policy now
uses:

```text
source_mode_shape = (q, q, L)
```

for a z-axis diatomic shared complete shell, with `L` selected by the
existing angular-resolution machinery rather than by raw index aspect.

## Construction Boundary

The bond-axis source length cannot be selected honestly in the early terminal
region contract builder, which sees region metadata but not all parent/bundle
and retention facts. It also cannot be repaired only in the multilayer source
plan after retained-rule and support records have already frozen a shape.

The implemented owner is the terminal low-order route path after:

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

`_pqs_source_box_route_driver_aspect_shell_lowering_plan(...)` enriches only
PQS `:complete_shell` contracts for shared molecular shells of bond-aligned
z-axis diatomics. Other region kinds, White-Lindsey lowering, one-center
construction, and nonshared shells retain their existing policies.

## Angular-Resolution Selection

For each eligible complete shell, the route owner calls the established
diatomic source-box dimension planner with:

- the parent axis bundle and atomic locations;
- the shell outer and inner-exclusion boxes;
- selected transverse `q`;
- bond axis;
- current complete-shell retention policy and support count;
- the established shared-shell angular-resolution scale.

The old angular-band machinery selects axis retained counts and returns:

```text
source_mode_dims
raw_source_dims
selected_q
raw_q
raw_L
axis_selector_retained_counts
```

For bond axis `z`, `raw_q` is the common transverse source size and
`raw_L` is the selected longitudinal size. The enriched lowering contract
records one authoritative shape and policy:

```text
source_mode_shape = source_mode_dims
source_mode_policy = :diatomic_shared_shell_adaptive_angular_source_box
```

A simple physical aspect estimate may remain a due-diligence comparison. It is
not the construction rule and must not overwrite the angular-band result.

## Downstream Shape Consistency

Every downstream representation of the shell must consume the same
`source_mode_shape`:

- lowering contract metadata;
- retained-unit and support metadata;
- raw-product retained rules;
- multilayer region and source plans;
- due-diligence actual-shape rows;
- final realization validation.

The multilayer source plan accepts noncubic `raw_source_dims`, requires both
transverse dimensions to agree, derives `q` and bond-axis `L`, and calls
`_nested_projected_q_shell_layer(...)` with explicit
`raw_source_dims/selected_q/q/L`. It must not replace `L` with `q` or
reconstruct a shape from the inner box after an authoritative contract shape
exists.

The retained complete-shell count follows the established raw-product boundary
rule, including:

```text
prod(source_mode_shape) - prod(source_mode_shape .- 2)
```

where that rule applies.

## Implemented Source Ownership

The implementation is in:

- `src/pqs_source_box_route_driver_helpers.jl` for angular-dimension
  selection and lowering-contract enrichment;
- `src/pqs_multilayer_shell_region_plan.jl` for carrying the contract shape;
- `src/pqs_multilayer_shell_source_plan.jl` for noncubic source realization;
- `src/cartesian_base_hamiltonian.jl` for corrected due-diligence expected
  shape reporting.

The older angular selector remains in its established diatomic owner. No
second selector or user-facing aspect parameter is introduced.

## Validation And Evidence

The accepted bounded H2/H2+ replay replaced hidden cubic `(5,5,5)` shapes
with angular-band selections:

```text
(5,5,8), (5,5,7), (5,5,6), (5,5,6)
```

with retained counts `146`, `130`, `114`, and `114`. The base final
dimension changed from the old cubic value `767` to `879`, as expected for
a basis-policy correction.

Acceptance checks include:

- exact shape agreement across lowering, retained, support, source-plan,
  due-diligence, and final realization records;
- noncubic multilayer source-plan acceptance;
- retained-count agreement with the active boundary rule;
- finite/symmetric base Hamiltonian matrices;
- artifact/readback matrix parity;
- disappearance of the stale
  `rectangular_physical_shell_cubic_source_modes` warning;
- continued reporting of physical-aspect differences where the advisory
  estimate and angular-band policy differ.

Old scalar dimensions or energies tied to cubic complete shells must be
remeasured, not preserved through compatibility logic.

## Failure Behavior

Stop if:

- the angular selector cannot produce finite, valid dimensions;
- transverse dimensions disagree;
- any downstream record carries a different source shape;
- support ownership or shell-local projection would need to change;
- a noncubic source shape cannot reach final realization without a new broad
  route or payload object.

Do not repair failure by restoring cubic shapes, guessing `L` from index
aspect, weakening support checks, or changing retained counts after the
contract is frozen.

## Explicit Non-Goals

This contract does not approve:

- public driver/API inputs or automatic source-shape tuning;
- White-Lindsey source-mode or retained-basis changes;
- shell ownership, direct/core identity, thin-slab, or angular-z-extension
  policy changes;
- artifact schema/provenance changes;
- residual/RG/MWG/IDA, injection, EGOI, or screened-Hartree changes;
- route-global materialization revival or a report/payload framework;
- solver workflow or Cr2 production claims.

Terminal due diligence remains advisory reporting. It exposes actual and
expected physical-aspect shapes but does not itself choose `L`.

The separately approved-pending
[semantic shell-q override](pqs_semantic_shell_q_overrides.md) may rerun this
same angular-band selector with a larger transverse source size for one matched
shared shell. It does not replace or modify the ordinary aspect policy.

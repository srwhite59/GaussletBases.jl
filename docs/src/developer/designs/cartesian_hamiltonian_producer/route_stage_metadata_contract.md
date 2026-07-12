# Route/Stage Metadata Contract

Status: the five route/stage source IDs are implemented. Their paired test IDs
are completed implementation-time evidence with no continuing permission.

This page is the canonical current contract for internal route/stage inventory,
plan, and carrier semantics established by those cleanup families. The
registry owns ID lifecycle, source/test surfaces, and maintenance permission;
this page does not independently authorize source work. Family-selective route
recipes and supported geometry/nesting/supplement composition are owned by
[Nesting/supplement composition](nesting_supplement_composition_plan.md).

## Boundary

Variable-size route facts are data, not type parameters. The inventories and
plan collections converted under the cleanup IDs use vectors or vector-backed
records. Dictionaries may provide label lookup, but dictionary iteration is
not ordering authority. Small fixed mathematical values such as three-axis
coordinates and dimensions may remain `NTuple{3,Int}`. The separately owned
`ShellificationPlan.terminal_regions` field remains tuple-backed in current
source and was not converted by these IDs. Accessor compatibility for the
converted inventories preserves facts, order, and iteration/indexing semantics
where consumed; it does not preserve their obsolete variable-length `Tuple`
return types.

"Metadata" here names the typed internal planning layer; it does not redefine
generic metadata or artifact provenance. These plans must not be persisted as
route dumps, become reader inputs, or create a parallel report/payload workflow.
A record marked metadata-only or not materialized is not numerical
construction. The existing validated source-axis transform-fact seam is a
separately governed exception; this contract neither expands nor reinterprets
it.

## Current Owners

- `src/pqs_source_box_route_driver_helpers.jl` owns ordered retained-unit and
  pair-family rows and the compact stage handoffs.
- `src/pqs_source_box_route_driver_skeletons.jl` retains ownership of upstream
  route-skeleton construction; the carrier cleanup did not retire it.
- `src/cartesian_terminal_shellification_geometry.jl` owns the private
  shellification-region inventory derived in scaffold order. It is descriptive
  shellification metadata, not a second lowering-plan owner.
- `src/cartesian_raw_product_sources/` owns `RawProductBoxPlan` source-mode
  inventories and their accessors.
- `src/cartesian_terminal_lowering/` owns available and selected terminal
  lowering contracts.
- `src/cartesian_retained_units/` owns the ordered retained-unit plan;
  `src/cartesian_retained_unit_transform_contracts/` owns one ordered transform
  contract per retained unit.
- `src/pqs_source_box_diatomic_complete_core_shell.jl` owns the active terminal
  topology support-region and terminal retained-rule plans consumed by
  `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`.

## Inventory And Plan Semantics

The retained-unit vector is the primary route inventory. Derived source-box,
source-dimension, retained-count, range, and pair-family count rows preserve
their producer order. The internal ordered-row carrier stores label and value
vectors plus a label-to-index lookup. Labels remain values; adding units or
families must not create a new concrete inventory type.

`RawProductBoxPlan` keeps fixed three-axis source dimensions and mode
coordinates while storing the variable source-mode and source-column sequences
as vectors. `source_mode_indices(plan)` follows the plan's declared ordering;
the current standard ordering is `:x_major_y_major_z_fast`. Column numbers must
remain aligned with that sequence, and retained-mode selection must preserve
the same mode/column association.

`TerminalLoweringPlan.available_contracts` and its selected `contracts` are
vectors. `available_contracts(plan)`, `selected_contracts(plan)`, and
`contracts(plan)` expose their respective order. Available alternatives and
selected route contracts remain distinct concepts even when the default plan
copies the selected sequence into both. The per-contract `source_cpbs` tuple
was outside this cleanup and is not made obsolete by it.

`RetainedUnitPlan.units` follows selected lowering-contract order.
`RetainedUnitTransformContractPlan.contracts` is a vector built one-for-one in
retained-unit order, and `transform_contracts(plan)` preserves that order.
Summaries may expose counts and kinds, but they are not alternate inventory
authorities.

## Stage Handoffs

The staged path remains:

```text
cartesian_parent
-> cartesian_shells  # constructs and consumes the route skeleton locally
-> cartesian_units
-> cartesian_transforms
-> cartesian_pair_terms
-> cartesian_assembly
-> cartesian_report
```

`cartesian_shells` constructs the route skeleton locally, extracts ordered
retained-unit and pair-entry vectors plus compact route/shellification facts,
and does not carry the complete skeleton onward. `cartesian_units` derives
ordered count/range rows and keeps
the retained-unit plan needed for downstream construction.
`cartesian_transforms` keeps the transform-contract plan, terminal retained-rule
plan, and terminal basis realization needed by consumers, without copying the
earlier shellification and lowering payloads through every later stage.

Later pair, assembly, and report stages carry compact pair and route summaries.
Their continued existence is outside these cleanup families: reduced metadata
is not evidence that the route skeleton, pair stage, assembly stage, report
stage, or any tool may be deleted or retired.

## Deterministic Order

The following associations must remain aligned:

1. shellification regions, private region records, and terminal lowering
   contracts;
2. selected lowering contracts, retained units, and transform contracts;
3. terminal support records, retained-rule records, and realized terminal
   blocks;
4. raw source modes, source columns, retained modes, and retained columns;
5. retained-unit labels, counts, ranges, and final retained dimension;
6. pair entries, pair keys, and pair-family counts.

Recomputation of a small summary is allowed only from the canonical ordered
object at that stage. A lookup dictionary, compact summary, status, or
fingerprint must not silently become a competing order source.

## Complete-Core-Shell Boundary

Complete-core-shell construction remains active for the supported one-center
atom and homonuclear z-axis diatomic producer paths. Common shellification may
emit `:complete_shell` regions; route-specific lowering then selects PQS filled
source-box or White-Lindsey boundary-stratum behavior. Terminal topology,
support coverage, retained counts, and realization continue through the plans
owned above.

This active construction is distinct from the deleted
`pqs_multilayer_complete_core_shell_rhf.jl` payload stack. That RHF machinery
was retired by `HP-RETIRE-CCS-RHF-*`; its deletion did not retire
complete-core-shell geometry, lowering, or terminal realization. Any surviving
compatibility/report field with RHF vocabulary is not executable RHF authority
and must not be used to restore the retired workflow.

## Preservation Guardrails

- Preserve route selection, shellification, lowering, retained rules, support
  coverage, terminal realization, numerical matrices, and deterministic order.
- Preserve the public producer and human-facing driver contracts, existing
  artifact schema/readback, and manifest behavior.
- Do not reintroduce runtime-keyed `NamedTuple` inventories, variable-size
  tuple carriers, duplicate lowering-contract mirrors, or adapters that hide
  those old shapes under new names.
- Do not infer source-cleanup authority, route-skeleton retirement, stage or
  tool deletion, diagnostic/report expansion, new payloads, or a public API
  from this metadata contract.
- Do not change raw Gaussian blocks, Residual Gaussian selection, MWG/IDA,
  Qiu-White semantics, numerical kernels, or Cr2 workflow in this lane.

Any later authorized maintenance should validate the affected order/alignment
contract and the smallest surviving base or supplemented endpoint. Terminal
realization changes additionally require its focused endpoint; metadata-only
maintenance does not justify a Cr2 gate.

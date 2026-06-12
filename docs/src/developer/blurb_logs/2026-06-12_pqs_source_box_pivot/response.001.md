PQS source-box-first restart audit completed. This pass was read-only except
for this response file and the tracked curated response copy.

Files/functions inspected:

- `AGENTS.md`
- `docs/src/developer/root_live_baton_loop.md`
- `.agent_handoffs/RUN.md`
- `.agent_handoffs/state.md`
- `.agent_handoffs/DOER_STARTUP.md`
- `docs/src/developer/pqs_source_box_operator_framework.md`
- `docs/src/developer/raw_product_source_retained_transform_policy.md`
- `docs/src/developer/projected_q_shell_policy.md`
- `docs/src/developer/cartesian_route_retirement_ledger.md`
- `src/cartesian_raw_product_sources/CartesianRawProductSources.jl`
- `src/cartesian_raw_product_sources/records.jl`
- `src/cartesian_raw_product_sources/source_mode_indices.jl`
- `src/cartesian_raw_product_sources/axis_transform_facts.jl`
- `src/cartesian_raw_product_sources/summaries.jl`
- `src/cartesian_retained_unit_transform_contracts/unit_contracts.jl`
- `src/cartesian_retained_unit_transform_contracts/summaries.jl`
- `src/cartesian_pair_operator_plans/plan_inventory.jl`
- `src/cartesian_pair_block_materialization/preflight.jl`
- `src/cartesian_pair_block_materialization/pqs_source_safe_terms.jl`
- `src/cartesian_pair_block_materialization/pqs_source_shell_bridge.jl`
- `src/cartesian_pair_block_materialization/pqs_source_final_readiness.jl`
- `src/cartesian_nested_faces.jl`
  - `_cartesian_raw_product_box_plan`
- `src/CartesianContractedParents.jl`
  - `_cartesian_raw_box_pqs_retained_rule_transform`
  - `_cartesian_raw_product_source_retained_transform`
  - `_cartesian_raw_product_source_pair_operator_packet`
  - `_cartesian_raw_product_source_pair_plan`
  - `_cartesian_raw_product_source_pair_plan_audit`
- `src/CartesianContractedParentMetrics.jl`
  - `_pqs_raw_product_box_plan`
  - `_pqs_product_source_box_pair_plan`
  - `_pqs_product_source_box_reference_blocks_from_pair_plan`
  - `_pqs_product_source_box_reference_blocks`
  - `_pqs_pqs_source_box_pair_plan`
  - `_pqs_pqs_source_box_reference_blocks_from_pair_plan`
  - `_pqs_pqs_source_box_reference_blocks`
  - `_pqs_pqs_source_box_*density_density*`
  - `_pqs_pqs_source_box_*local_gaussian*`
- `test/nested/cartesian_raw_product_sources_contract_runtests.jl`
- `test/nested/cartesian_retained_unit_transform_contracts_runtests.jl`
- `test/nested/cartesian_pair_operator_plans_contract_runtests.jl`
- `test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
- `test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl`
- `test/nested/pqs_raw_product_box_plan_probe_runtests.jl`
- `test/nested/pqs_source_box_route_driver_report_runtests.jl`
- `test/nested/pqs_route_driver_one_center_materializer_probe_runtests.jl`

Audit answers:

1. Raw product source-box ownership

The module-spine owner is `CartesianRawProductSources`.

- `RawProductBoxPlan` in `src/cartesian_raw_product_sources/records.jl`
  owns:
  - source CPB;
  - source intervals/shape;
  - total `source_mode_dims`;
  - `source_mode_count`;
  - deterministic `source_mode_indices`;
  - source-mode column indices;
  - source-mode ordering;
  - metadata-only axis transform facts.
- `raw_product_box_plan(source_box; source_mode_dims, source_key,
  source_mode_ordering, metadata)` constructs the module-owned source facts.
- `source_mode_indices(...)` currently supports the z-fast ordering
  `:x_major_y_major_z_fast`.
- `AxisSourceTransformFact` deliberately records transform facts as metadata
  only. Default facts have `coefficient_status = :not_materialized` and
  `coefficient_matrix = nothing`.

The module intentionally does not own retained rules, shell projection, Lowdin,
final retained units, pair blocks, Hamiltonians, artifacts, or IDA weights.

The current module-spine producer that attaches raw product source facts to PQS
units is `CartesianRetainedUnitTransformContracts`:

- `_raw_product_source_contract_metadata(unit)` in
  `src/cartesian_retained_unit_transform_contracts/unit_contracts.jl`
  recognizes `:pqs_shell_retained_unit`;
- `_pqs_source_mode_dims(unit)` uses `unit.metadata.source_mode_shape`, or
  falls back to `(q, q, q)`;
- `_pqs_source_cpb(unit)` requires exactly one filled/codimension-0 source CPB;
- it then calls `CartesianRawProductSources.raw_product_box_plan(...)` and
  stores `raw_product_source_plan`, `raw_product_source_summary`, and
  `raw_product_source_plan_status` in transform-contract metadata.

Older source-box helpers still exist outside the module spine:

- `_cartesian_raw_product_box_plan(...)` in `cartesian_nested_faces.jl` builds
  numerical local axis transforms and wraps a `RawProductBoxPlan`;
- `_pqs_raw_product_box_plan(...)` in `CartesianContractedParentMetrics.jl`
  is the old PQS-shaped raw-plan adapter with boundary selector and operator
  factors;
- `_cartesian_raw_product_source_retained_transform(...)` and
  `_cartesian_raw_product_source_pair_plan(...)` in
  `CartesianContractedParents.jl` are fixture/oracle raw-source retained
  transform and pair-plan surfaces.

Those older surfaces should be treated as oracle/reference or migration
material, not the new route authority.

2. PQS/PQS raw source-space one-body blocks

The module-spine materialization owner is
`CartesianPairBlockMaterialization`, especially
`src/cartesian_pair_block_materialization/pqs_source_safe_terms.jl`.

Live safe-term entry points are:

- `pqs_source_pair_overlap_block(record; overlap_1d)`
- `pqs_source_pair_position_block(record; axis, overlap_1d, position_1d)`
- `pqs_source_pair_x2_block(record; axis, overlap_1d, x2_1d)`
- `pqs_source_pair_kinetic_block(record; overlap_1d, kinetic_1d)`
- plural plan-level variants:
  - `pqs_source_pair_overlap_blocks(...)`
  - `pqs_source_pair_position_blocks(...)`
  - `pqs_source_pair_x2_blocks(...)`
  - `pqs_source_pair_kinetic_blocks(...)`
- selector entry points:
  - `pqs_source_pair_one_body_block(...)`
  - `pqs_source_pair_one_body_blocks(...)`

They consume ready `PairBlockMaterializationRecord`s whose
`materialization_path == :pqs_source_pair_preflight`, require source-mode dims
and ordering in metadata, and build dense raw source-space blocks from
caller-supplied 1D factors. Kinetic is assembled as:

```text
Kx Sy Sz + Sx Ky Sz + Sx Sy Kz
```

Their results explicitly report:

- `block_space = :raw_product_source_modes`;
- `source_operator_blocks_materialized = true`;
- `shell_realization_materialized = false`;
- `final_pair_blocks_materialized = false`;
- `operator_blocks_materialized = false`;
- no Hamiltonian/export/artifact materialization.

The preflight records feeding these helpers are made in:

- `pair_block_materialization_plan(...)` in `preflight.jl`;
- `_is_pqs_pqs_source_pair_preflight(record)`;
- `_pqs_source_pair_preflight_status(...)`;
- `_pair_block_materialization_record_metadata(...)`.

The source operator path originates in `CartesianPairOperatorPlans`:

- `_source_operator_path(pair)` returns
  `:pqs_source_cpb_1d_factor_path` when either side is
  `:pqs_shell_retained_unit`;
- `_final_block_path(...)` returns
  `:source_block_realization_then_final_block` when either realization path is
  `:shell_projection_lowdin_planned`.

3. Shell-realization/support-row surfaces and classification

Shell-realization surfaces exist, but should be oracle/adapter/debug only at
this point.

Current module-spine metadata bridge/readiness:

- `pqs_source_pair_shell_realization_bridge_summary(...)` in
  `pqs_source_shell_bridge.jl` records the planned handoff from source-space
  blocks to shell realization. It is metadata only.
- `pqs_source_pair_final_block_readiness_summary(...)` in
  `pqs_source_final_readiness.jl` currently remains blocked by
  `:shell_realization_not_materialized` for available source-space blocks.

Current retained-unit transform contract metadata:

- `:pqs_source_modes_boundary_selection_shell_realization_contract`
  is the transform path for PQS retained units.
- `:shell_projection_lowdin_planned` is the realization path.

Older shell/support-row oracle surfaces:

- `_cartesian_raw_product_source_retained_transform(...)` on fixed-block or
  sidecar payloads produces shell-realized support-local fixture transforms
  with `transform_kind = :boundary_projection_lowdin`.
- `_cartesian_raw_box_pqs_retained_rule_transform(...)` produces the raw-box
  boundary selector metadata transform with
  `transform_kind = :boundary_comx_product_mode_selection` and explicitly
  reports no shell projection/Lowdin.
- old PQS/product and PQS/PQS CCPM helpers compare against support-local or
  explicit raw-box references.

Classification:

- Raw source-space blocks are the intended route primitive.
- Shell projection plus Lowdin is a later realization stage.
- Support-row/shell-row contraction and old fixed-block sidecars should remain
  oracle, adapter, or debug surfaces until a reviewed shell-realization
  implementation consumes source-space blocks. They should not become the
  source-box operator algorithm.

4. Smallest first retained-rule object

The smallest useful object is a source-mode retained rule that lives with the
raw source facts and is independent of shell realization. It should represent:

```text
RawProductBoxPlan + PQSBoundaryProductModeSelector + retained transform
```

For the first pass it can be a narrow internal object, conceptually:

```text
PQSRetainedRule
  source_key
  source_mode_dims
  source_mode_ordering
  retained_rule_kind = :boundary_comx_product_mode_selection
  retained_mode_indices
  retained_column_indices
  retained_count
  transform_kind = :source_mode_column_selector
  shell_realization_materialized = false
  lowdin_cleanup_used = false
```

and a corresponding transform:

```text
PQSRetainedTransform
  raw_source_plan
  retained_rule
  transform_matrix or selector representation
  retained_dimension
  transform_stage = :raw_product_box_modes_to_retained_boundary_modes
```

For the first implementation, storing selector columns is enough. A dense
`125 x 98` selector matrix is acceptable for a tiny module-contract fixture if
it makes `T_left' * O_source * T_right` explicit, but production shape should
prefer selector metadata and only materialize a dense transform when a caller
asks for it.

This object should not include shell projection, Lowdin cleanup, support-local
coefficients, final IDA weights, global ranges, Hamiltonian facts, or driver
report aliases.

5. Boundary product selector feasibility

Yes. The first retained rule can be the explicit boundary product selector:

```text
source modes: 5 x 5 x 5 = 125
interior modes: 3 x 3 x 3 = 27
retained boundary modes: 125 - 27 = 98
```

This matches:

- `projected_q_shell_policy.md`, which states `PQS(5,5)` has
  `5^3 - 3^3 = 98`;
- the old integration checks in
  `test/nested/pqs_projected_q_shell_local_layer_integration_runtests.jl`;
- `_cartesian_raw_box_pqs_retained_rule_transform(...)`, which already records
  `:boundary_comx_product_mode_selection` as the raw-box retained rule.

The important implementation detail is that this is source-mode boundary
selection, not shell-row boundary support selection.

6. Exact next implementation seam

Recommended next blurb:

Implement a tiny internal source-mode retained-rule object in or directly next
to `CartesianRawProductSources`, then thread it into
`CartesianRetainedUnitTransformContracts` for PQS units.

Concrete first source surfaces:

- `src/cartesian_raw_product_sources/records.jl`
  - add a compact `PQSRetainedRule` or similarly named internal record;
  - add a constructor such as
    `pqs_boundary_product_mode_retained_rule(raw_plan)` or
    `pqs_boundary_product_mode_retained_rule(source_mode_dims;
    source_mode_ordering = :x_major_y_major_z_fast)`;
  - compute retained boundary mode indices/columns by selecting source-mode
    tuples with at least one local index equal to `1` or that axis dimension;
  - record `retained_count = 98` for `(5,5,5)`.
- `src/cartesian_raw_product_sources/summaries.jl`
  - expose a compact summary: status, dims, retained count, selector rule,
    no shell realization, no Lowdin.
- `src/cartesian_retained_unit_transform_contracts/unit_contracts.jl`
  - when `_raw_product_source_contract_metadata(unit)` builds an available raw
    product source plan for a `:pqs_shell_retained_unit`, also attach the new
    retained rule/summary to metadata.

Then, in a follow-up implementation after that object exists:

- update `CartesianPairBlockMaterialization` PQS source block materialization
  or a new retained-block helper to apply:

```text
O_retained = T_left' * O_source * T_right
```

using selector columns for the source-mode boundary retained rule.

Do not start with shell realization. Do not add Lowdin or support-row
contraction.

7. Likely oracle/removal targets after retained source-mode contraction works

Likely oracle-only or removable surfaces:

- CCPM `_pqs_pqs_source_box_reference_blocks*` and
  `_pqs_product_source_box_reference_blocks*` should become oracle/reference
  checks rather than the active source-box retained-block implementation.
- `_pqs_pqs_source_box_explicit_boundary_selection_reference(...)` should
  remain a compact oracle only, not production logic.
- Shell/support-row tests in
  `pqs_projected_q_shell_local_layer_integration_runtests.jl` that assert
  support-local fixture vocabulary should be shrunk once source-mode retained
  contraction is covered through module-level contracts.
- `pqs_source_pair_shell_realization_bridge_summary(...)` and
  `pqs_source_pair_final_block_readiness_summary(...)` should stay as
  metadata bridge/readiness surfaces until real shell realization exists, but
  they should not grow into another metadata framework.
- Driver/report compatibility in
  `pqs_source_box_route_driver_helpers.jl`, `pqs_source_box_route_driver_skeletons.jl`,
  and `pqs_source_box_route_driver_reporting.jl` remains compatibility
  scaffolding until module-owned source-box retained transforms can feed a
  real route state.

Deletion/shrinkage report:

- Nothing was deleted in this read-only audit.
- No existing surface is obsolete yet because the first retained source-mode
  rule/transform object has not been implemented in the module spine.
- Once retained source-mode contraction lands, old CCPM source-box reference
  helpers and integration tests that preserve support-row/shell-realization
  vocabulary should be reduced to oracle-only coverage or deleted where no
  longer called.
- Shell realization should stay adapter/oracle/debug only until a later blurb
  explicitly implements that stage.

Recommended next implementation blurb:

Add a compact PQS boundary source-mode retained rule under
`CartesianRawProductSources` and attach it to PQS retained-unit transform
contracts.

Scope for that pass:

- source modes only;
- `(5,5,5)` boundary selector contract with retained count `98`;
- selector columns/mode indices and compact summary;
- no shell projection;
- no Lowdin;
- no pair-block materialization yet except perhaps a tiny constructor contract
  test;
- no driver/report changes.

Validation run:

- `git status --short --branch`
- `git diff --check`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`

Validation result:

- branch status before response: `## main...origin/main`
- `git diff --check`: passed
- load check: `load ok`

-- repo-doer@macmini

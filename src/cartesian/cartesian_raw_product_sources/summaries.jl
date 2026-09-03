# Compact raw product source summaries.

"""
    summary(plan::RawProductBoxPlan)

Return a compact, stable summary of a raw product source plan.

The summary is intended for reports and tests that need source shape,
source-mode dimensions, source-mode count, ordering, axis-transform statuses,
and materialization flags. It deliberately does not expose large source-mode
inventories or numerical transform data.
"""
function summary(plan::RawProductBoxPlan)
    return (;
        object_kind = :raw_product_box_plan_summary,
        status = :available_raw_product_box_plan,
        source_key = plan.source_key,
        source_shape = plan.source_shape,
        source_mode_dims = plan.source_mode_dims,
        source_mode_count = plan.source_mode_count,
        source_mode_ordering = plan.source_mode_ordering,
        axis_transform_statuses =
            Tuple(fact.coefficient_status for fact in plan.axis_transform_facts),
        materialized = plan.materialized,
        retained_rule_materialized = false,
        shell_realization_materialized = false,
        pair_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
    )
end

"""
    summary(rule::PQSBoundaryProductModeRetainedRule)

Return a compact summary of a PQS retained source-mode boundary rule.

This summary exposes rule kind, retained count, ordering, and materialization
nonclaims. It does not duplicate retained mode inventories or numerical
transform data.
"""
function summary(rule::PQSBoundaryProductModeRetainedRule)
    return (;
        object_kind = :pqs_boundary_product_mode_retained_rule_summary,
        status = :available_pqs_boundary_product_mode_retained_rule,
        source_key = rule.source_key,
        source_mode_dims = rule.source_mode_dims,
        source_mode_ordering = rule.source_mode_ordering,
        retained_rule_kind = rule.retained_rule_kind,
        retained_count = rule.retained_count,
        transform_kind = rule.transform_kind,
        shell_realization_materialized = rule.shell_realization_materialized,
        lowdin_cleanup_used = rule.lowdin_cleanup_used,
        transforms_materialized = false,
        coefficient_maps_materialized = false,
        pair_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
    )
end

"""
    unavailable_summary(status, blocker = nothing)

Return the compact summary shape used when a raw product source plan is not
available.

This keeps downstream metadata paths explicit about why raw source facts are
absent, while preserving the same materialization flags as an available
metadata-only plan.
"""
function unavailable_summary(status::Symbol, blocker = nothing)
    return (;
        object_kind = :raw_product_box_plan_summary,
        status,
        blocker,
        source_shape = nothing,
        source_mode_dims = nothing,
        source_mode_count = 0,
        source_mode_ordering = nothing,
        axis_transform_statuses = (),
        materialized = false,
        retained_rule_materialized = false,
        shell_realization_materialized = false,
        pair_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
    )
end

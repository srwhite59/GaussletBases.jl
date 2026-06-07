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

# Compact raw product source summaries.

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

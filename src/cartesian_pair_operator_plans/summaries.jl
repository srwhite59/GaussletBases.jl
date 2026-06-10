# Compact pair-operator plan summaries.

function _count_by_value(values, field::Symbol)
    counts = Dict{Any,Int}()
    order = Any[]
    for value in values
        if !haskey(counts, value)
            counts[value] = 0
            push!(order, value)
        end
        counts[value] += 1
    end
    return Tuple((; field => value, count = counts[value]) for value in order)
end

function _pair_operator_plan_summary(
    policy::PairOperatorPlanPolicy,
    unit_pair_plan::CartesianUnitPairs.UnitPairPlan,
    records,
    route_core_plan,
)
    unit_pair_summary = CartesianUnitPairs.summary(unit_pair_plan)
    blockers = Tuple(record.blocker for record in records)
    blocked_count = count(blocker -> !isnothing(blocker), blockers)
    status =
        blocked_count == 0 ?
        :available_pair_operator_plan :
        :blocked_pair_operator_plan
    return (;
        object_kind = :cartesian_pair_operator_plan_summary,
        status,
        blocker =
            blocked_count == 0 ?
            nothing :
            first(blocker for blocker in blockers if !isnothing(blocker)),
        policy_kind = policy_kind(policy),
        retained_unit_count = unit_pair_summary.retained_unit_count,
        unit_pair_count = unit_pair_summary.pair_count,
        pair_operator_plan_count = length(records),
        expected_pair_operator_plan_count =
            unit_pair_summary.expected_upper_triangular_pair_count,
        pair_families = Tuple(record.pair_family for record in records),
        pair_family_counts =
            _count_by_value((record.pair_family for record in records), :pair_family),
        source_operator_path_counts =
            _count_by_value(
                (record.source_operator_path for record in records),
                :source_operator_path,
            ),
        transform_path_counts =
            _count_by_value((record.transform_path for record in records), :transform_path),
        realization_path_counts =
            _count_by_value(
                (record.realization_path for record in records),
                :realization_path,
            ),
        final_block_path_counts =
            _count_by_value((record.final_block_path for record in records), :final_block_path),
        materialization_status_counts =
            _count_by_value((record.status for record in records), :materialization_status),
        blocker_counts = _count_by_value(blockers, :blocker),
        blocked_pair_operator_plan_count = blocked_count,
        route_core_pair_operator_plan_inventory_available =
            !isnothing(route_core_plan.inventory),
        route_core_pair_operator_plan_inventory_status = route_core_plan.status,
        route_core_pair_operator_plan_inventory_blocker = route_core_plan.blocker,
        route_core_pair_operator_plan_count =
            isnothing(route_core_plan.inventory) ?
            0 :
            CartesianRouteCore.pair_operator_plan_count(route_core_plan.inventory),
        route_core_pair_operator_plan_blocked_count =
            isnothing(route_core_plan.inventory) ?
            0 :
            count(
                plan -> !isnothing(CartesianRouteCore.blocker(plan)),
                CartesianRouteCore.pair_operator_plans(route_core_plan.inventory),
            ),
        materialized = false,
        source_operator_blocks_materialized = false,
        final_pair_blocks_materialized = false,
        operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
    )
end

function unavailable_summary(status::Symbol, blocker = nothing)
    return (;
        object_kind = :cartesian_pair_operator_plan_summary,
        status,
        blocker,
        retained_unit_count = 0,
        unit_pair_count = 0,
        pair_operator_plan_count = 0,
        expected_pair_operator_plan_count = 0,
        pair_families = (),
        pair_family_counts = (),
        source_operator_path_counts = (),
        transform_path_counts = (),
        realization_path_counts = (),
        final_block_path_counts = (),
        materialization_status_counts = (),
        blocker_counts = (),
        blocked_pair_operator_plan_count = 0,
        route_core_pair_operator_plan_inventory_available = false,
        route_core_pair_operator_plan_inventory_status = :not_available,
        route_core_pair_operator_plan_inventory_blocker = blocker,
        route_core_pair_operator_plan_count = 0,
        route_core_pair_operator_plan_blocked_count = 0,
        materialized = false,
        source_operator_blocks_materialized = false,
        final_pair_blocks_materialized = false,
        operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
    )
end

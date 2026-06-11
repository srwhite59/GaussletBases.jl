# Compact retained-unit pair summaries.

function _pair_family_counts(pairs)
    counts = Dict{Symbol,Int}()
    order = Symbol[]
    for pair in pairs
        family = pair.pair_family
        if !haskey(counts, family)
            counts[family] = 0
            push!(order, family)
        end
        counts[family] += 1
    end
    return (
        Tuple(order),
        Tuple((; pair_family = family, count = counts[family]) for family in order),
    )
end

function _unit_pair_plan_summary(
    policy::UnitPairPolicy,
    retained_unit_plan::CartesianRetainedUnits.RetainedUnitPlan,
    pairs,
    route_core_inventory,
)
    pair_families, pair_family_counts = _pair_family_counts(pairs)
    retained_unit_count = length(CartesianRetainedUnits.units(retained_unit_plan))
    return (;
        object_kind = :cartesian_unit_pair_plan_summary,
        status = :available_unit_pair_plan,
        policy_kind = policy_kind(policy),
        retained_unit_count,
        pair_count = length(pairs),
        expected_upper_triangular_pair_count =
            div(retained_unit_count * (retained_unit_count + 1), 2),
        pair_families,
        pair_family_counts,
        pair_storage = pairs isa UnitPairIndexTable ? :unit_pair_index_table :
            :diagnostic_legacy_rich_unit_pair_records,
        rich_unit_pair_records_stored = pairs isa Vector{UnitPairRecord},
        route_core_pair_inventory_available =
            !isnothing(route_core_inventory.inventory),
        route_core_pair_inventory_status = route_core_inventory.status,
        route_core_pair_inventory_blocker = route_core_inventory.blocker,
        route_core_pair_count =
            isnothing(route_core_inventory.inventory) ?
            0 :
            length(CartesianRouteCore.pair_entries(route_core_inventory.inventory)),
        route_core_pair_missing_final_unit_indices =
            route_core_inventory.missing_indices,
        materialized = false,
        pair_inventory_materialized = false,
        source_operator_blocks_materialized = false,
        operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
    )
end

function unavailable_summary(status::Symbol, blocker = nothing)
    return (;
        object_kind = :cartesian_unit_pair_plan_summary,
        status,
        blocker,
        retained_unit_count = 0,
        pair_count = 0,
        expected_upper_triangular_pair_count = 0,
        pair_families = (),
        pair_family_counts = (),
        route_core_pair_inventory_available = false,
        route_core_pair_inventory_status = :not_available,
        route_core_pair_inventory_blocker = blocker,
        route_core_pair_count = 0,
        route_core_pair_missing_final_unit_indices = (),
        materialized = false,
        pair_inventory_materialized = false,
        source_operator_blocks_materialized = false,
        operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
    )
end

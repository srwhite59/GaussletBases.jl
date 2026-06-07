# Compact retained-unit pair summaries.

function _count_by_symbol(values, field::Symbol)
    counts = Dict{Symbol,Int}()
    order = Symbol[]
    for value in values
        if !haskey(counts, value)
            counts[value] = 0
            push!(order, value)
        end
        counts[value] += 1
    end
    return Tuple((; field => value, count = counts[value]) for value in order)
end

function _unit_pair_plan_summary(
    policy::UnitPairPolicy,
    retained_unit_plan::CRU.RetainedUnitPlan,
    pairs,
    route_core_inventory,
)
    pair_families = Tuple(pair.pair_family for pair in pairs)
    retained_unit_count = length(CRU.units(retained_unit_plan))
    return (;
        object_kind = :cartesian_unit_pair_plan_summary,
        status = :available_unit_pair_plan,
        policy_kind = policy_kind(policy),
        retained_unit_count,
        pair_count = length(pairs),
        expected_upper_triangular_pair_count =
            retained_unit_count * (retained_unit_count + 1) ÷ 2,
        pair_families,
        pair_family_counts = _count_by_symbol(pair_families, :pair_family),
        route_core_pair_inventory_available =
            !isnothing(route_core_inventory.inventory),
        route_core_pair_inventory_status = route_core_inventory.status,
        route_core_pair_inventory_blocker = route_core_inventory.blocker,
        route_core_pair_count =
            isnothing(route_core_inventory.inventory) ?
            0 :
            length(CRC.pair_entries(route_core_inventory.inventory)),
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

# Compact retained-unit plan summaries.

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

function _unique_symbols(values)
    seen = Set{Symbol}()
    ordered = Symbol[]
    for value in values
        if !(value in seen)
            push!(seen, value)
            push!(ordered, value)
        end
    end
    return Tuple(ordered)
end

function _retained_unit_plan_summary(
    policy::RetainedUnitPolicy,
    lowering_plan::CTL.TerminalLoweringPlan,
    planned_units,
)
    unit_kinds = Tuple(unit.unit_kind for unit in planned_units)
    lowering_kinds = Tuple(unit.lowering_kind for unit in planned_units)
    terminal_region_keys = _unique_symbols(
        Tuple(unit.terminal_region_key for unit in planned_units),
    )
    route_core_statuses =
        Tuple(unit.metadata.route_core_sidecar_status for unit in planned_units)

    return (;
        object_kind = :cartesian_retained_unit_plan_summary,
        status = :available_retained_unit_plan,
        policy_kind = policy_kind(policy),
        lowering_policy_kind = CTL.summary(lowering_plan).policy_kind,
        selected_contract_count = length(CTL.selected_contracts(lowering_plan)),
        retained_unit_count = length(planned_units),
        unit_kinds,
        unit_kind_counts = _count_by_symbol(unit_kinds, :unit_kind),
        lowering_kind_counts = _count_by_symbol(lowering_kinds, :lowering_kind),
        terminal_region_count = length(terminal_region_keys),
        route_core_final_unit_count = count(!isnothing(unit.route_core_final_unit) for unit in planned_units),
        route_core_final_unit_available_count = count(==(:available), route_core_statuses),
        route_core_final_unit_blocked_count = count(==(:blocked), route_core_statuses),
        materialized = false,
        transforms_materialized = false,
        coefficient_maps_materialized = false,
        pair_inventory_materialized = false,
        operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
    )
end

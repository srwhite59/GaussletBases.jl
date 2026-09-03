# Compact summaries for terminal lowering plans.

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

function _terminal_lowering_summary(
    policy::TerminalLoweringPolicy,
    regions,
    available,
    selected,
)
    selected_kinds = Tuple(contract.lowering_kind for contract in selected)
    available_kinds = Tuple(contract.lowering_kind for contract in available)
    selected_source_cpb_counts =
        Tuple(length(contract.source_cpbs) for contract in selected)

    return (;
        object_kind = :cartesian_terminal_lowering_plan_summary,
        status = :available_terminal_lowering_plan,
        policy_kind = policy_kind(policy),
        terminal_region_count = length(regions),
        available_contract_count = length(available),
        selected_contract_count = length(selected),
        available_contract_kinds = available_kinds,
        selected_contract_kinds = selected_kinds,
        available_contract_kind_counts =
            _count_by_symbol(available_kinds, :lowering_kind),
        selected_contract_kind_counts =
            _count_by_symbol(selected_kinds, :lowering_kind),
        selected_contract_source_cpb_counts = selected_source_cpb_counts,
        selected_source_cpb_count = sum(selected_source_cpb_counts; init = 0),
        all_terminal_regions_have_selected_contract =
            length(selected) == length(regions),
        materialized = false,
        coefficient_maps_materialized = false,
        transforms_materialized = false,
        retained_spaces_materialized = false,
        final_retained_units_materialized = false,
        pair_inventory_materialized = false,
        operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
    )
end

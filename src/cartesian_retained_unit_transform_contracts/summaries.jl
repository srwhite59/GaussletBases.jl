# Compact retained-unit transform-contract summaries.

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

function _retained_unit_transform_contract_plan_summary(
    policy::RetainedUnitTransformContractPolicy,
    retained_unit_plan::CRU.RetainedUnitPlan,
    contracts,
)
    transform_paths = Tuple(contract.transform_path for contract in contracts)
    realization_paths = Tuple(contract.realization_path for contract in contracts)
    blockers = Tuple(
        contract.blocker
        for contract in contracts
        if !isnothing(contract.blocker)
    )

    return (;
        object_kind = :cartesian_retained_unit_transform_contract_plan_summary,
        status =
            isempty(blockers) ?
            :available_retained_unit_transform_contract_plan :
            :blocked_retained_unit_transform_contract_plan,
        blocker = isempty(blockers) ? nothing : first(blockers),
        policy_kind = policy_kind(policy),
        retained_unit_count = length(CRU.units(retained_unit_plan)),
        transform_contract_count = length(contracts),
        transform_path_counts =
            _count_by_symbol(transform_paths, :transform_path),
        realization_path_counts =
            _count_by_symbol(realization_paths, :realization_path),
        blocked_contract_count = length(blockers),
        raw_product_source_plan_available_count =
            _raw_product_source_plan_available_count(contracts),
        raw_product_source_plan_blocked_count =
            _raw_product_source_plan_blocked_count(contracts),
        materialized = false,
        transforms_materialized = false,
        coefficient_maps_materialized = false,
        pair_inventory_materialized = false,
        operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
    )
end

function _raw_product_source_plan_status(contract)
    haskey(contract.metadata, :raw_product_source_plan_status) ||
        return nothing
    return contract.metadata.raw_product_source_plan_status
end

function _raw_product_source_plan_available_count(contracts)
    return count(
        contract -> _raw_product_source_plan_status(contract) ===
                    :available_raw_product_box_plan,
        contracts,
    )
end

function _raw_product_source_plan_blocked_count(contracts)
    return count(contracts) do contract
        status = _raw_product_source_plan_status(contract)
        return !isnothing(status) && status !== :available_raw_product_box_plan
    end
end

function unavailable_summary(status::Symbol, blocker = nothing)
    return (;
        object_kind = :cartesian_retained_unit_transform_contract_plan_summary,
        status,
        blocker,
        policy_kind = :unavailable_retained_unit_transform_contracts,
        retained_unit_count = 0,
        transform_contract_count = 0,
        transform_path_counts = (),
        realization_path_counts = (),
        blocked_contract_count = 0,
        raw_product_source_plan_available_count = 0,
        raw_product_source_plan_blocked_count = 0,
        materialized = false,
        transforms_materialized = false,
        coefficient_maps_materialized = false,
        pair_inventory_materialized = false,
        operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
    )
end

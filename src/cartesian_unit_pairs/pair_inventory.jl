# Retained-unit to upper-triangular pair inventory rules.

"""
    unit_pair_plan(retained_unit_plan; policy = MetadataOnlyUnitPairs())

Build metadata-only upper-triangular retained-unit pairs from a retained-unit
plan. Pair planning starts from retained units; this function does not inspect
shellification regions, CPBs, or lowering contracts directly.
"""
function unit_pair_plan(
    retained_unit_plan::CRU.RetainedUnitPlan;
    policy::UnitPairPolicy = MetadataOnlyUnitPairs(),
    metadata = (;),
)
    retained_units = CRU.units(retained_unit_plan)
    isempty(retained_units) &&
        throw(ArgumentError("unit_pair_plan requires at least one retained unit"))

    route_core_inventory = _route_core_pair_inventory_or_nothing(retained_units)
    route_core_pairs =
        isnothing(route_core_inventory.inventory) ?
        nothing :
        CRC.pair_entries(route_core_inventory.inventory)

    pair_records = UnitPairRecord[]
    pair_index = 0
    for left_index in eachindex(retained_units)
        for right_index in left_index:length(retained_units)
            pair_index += 1
            left = retained_units[left_index]
            right = retained_units[right_index]
            sidecar =
                isnothing(route_core_pairs) ? nothing : route_core_pairs[pair_index]
            push!(
                pair_records,
                _unit_pair_record(
                    pair_index,
                    left,
                    right,
                    left_index,
                    right_index,
                    sidecar,
                    route_core_inventory,
                ),
            )
        end
    end

    pairs = Tuple(pair_records)
    plan_summary =
        _unit_pair_plan_summary(policy, retained_unit_plan, pairs, route_core_inventory)
    return UnitPairPlan(
        policy,
        retained_unit_plan,
        pairs,
        route_core_inventory.inventory,
        plan_summary,
        NamedTuple(metadata),
    )
end

function _route_core_pair_inventory_or_nothing(retained_units)
    missing_indices = Tuple(
        index for (index, unit) in enumerate(retained_units)
        if isnothing(unit.route_core_final_unit)
    )
    if !isempty(missing_indices)
        return (;
            inventory = nothing,
            status = :blocked_missing_route_core_final_units,
            blocker = :missing_route_core_final_units,
            missing_indices,
        )
    end

    route_core_units = Tuple(unit.route_core_final_unit for unit in retained_units)
    try
        inventory = CRC.unit_pair_inventory(
            route_core_units;
            metadata = (; source = :cartesian_unit_pairs),
        )
        return (;
            inventory,
            status = :available_route_core_pair_inventory,
            blocker = nothing,
            missing_indices = (),
        )
    catch err
        return (;
            inventory = nothing,
            status = :blocked_route_core_pair_inventory_error,
            blocker = sprint(showerror, err),
            missing_indices = (),
        )
    end
end

function _unit_pair_record(
    pair_index::Int,
    left::CRU.RetainedUnitRecord,
    right::CRU.RetainedUnitRecord,
    left_index::Int,
    right_index::Int,
    route_core_sidecar,
    route_core_inventory,
)
    sidecar_status =
        isnothing(route_core_sidecar) ?
        route_core_inventory.status :
        :available_route_core_unit_pair
    return UnitPairRecord(
        (left.unit_key, right.unit_key),
        pair_index,
        _pair_family(left, right),
        left,
        right,
        left_index,
        right_index,
        left.unit_key,
        right.unit_key,
        left.unit_kind,
        right.unit_kind,
        route_core_sidecar,
        false,
        (;
            route_core_pair_sidecar_status = sidecar_status,
            route_core_pair_sidecar_blocker =
                isnothing(route_core_sidecar) ? route_core_inventory.blocker : nothing,
        ),
    )
end

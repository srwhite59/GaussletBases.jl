# Retained-unit to upper-triangular pair inventory rules.

"""
    unit_pair_plan(retained_unit_plan; policy = MetadataOnlyUnitPairs())

Build metadata-only upper-triangular retained-unit pairs from a retained-unit
plan. Pair planning starts from retained units; this function does not inspect
shellification regions, CPBs, or lowering contracts directly.
"""
function unit_pair_plan(
    retained_unit_plan::CartesianRetainedUnits.RetainedUnitPlan;
    policy::UnitPairPolicy = MetadataOnlyUnitPairs(),
    metadata = (;),
)
    retained_units = CartesianRetainedUnits.units(retained_unit_plan)
    isempty(retained_units) &&
        throw(ArgumentError("unit_pair_plan requires at least one retained unit"))

    route_core_inventory = _route_core_pair_inventory_or_nothing(retained_units)
    pairs = unit_pair_index_table(
        retained_units;
        metadata = (;
            retained_pair_source = :upper_triangular_unit_index_table,
            rich_unit_pair_record_storage = :not_stored,
            route_core_pair_sidecars_stored_on_pairs = false,
        ),
    )
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
        inventory = CartesianRouteCore.unit_pair_inventory(
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

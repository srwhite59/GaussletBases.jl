# Retained-unit pair planning records.

"""
    UnitPairPolicy

Policy marker for converting retained-unit plans into pair-inventory metadata.
"""
abstract type UnitPairPolicy end

"""
    MetadataOnlyUnitPairs()

Default pair policy. It builds upper-triangular retained-unit pair records and
RouteCore pair sidecars where all retained units have RouteCore final-unit
sidecars. It does not materialize operator blocks.
"""
struct MetadataOnlyUnitPairs <: UnitPairPolicy
end

policy_kind(::MetadataOnlyUnitPairs) = :metadata_only_unit_pairs

"""
    UnitPairRecord

Metadata-only upper-triangular pair of retained units. The optional RouteCore
sidecar is bookkeeping for later pair/operator planning, not a numerical
operator block.
"""
struct UnitPairRecord
    pair_key::Tuple{Symbol,Symbol}
    pair_index::Int
    pair_family::Symbol
    left_unit::CartesianRetainedUnits.RetainedUnitRecord
    right_unit::CartesianRetainedUnits.RetainedUnitRecord
    left_index::Int
    right_index::Int
    left_unit_key::Symbol
    right_unit_key::Symbol
    left_unit_kind::Symbol
    right_unit_kind::Symbol
    route_core_pair_sidecar::Union{CartesianRouteCore.UnitPair,Nothing}
    materialized::Bool
    metadata::NamedTuple
end

"""
    UnitPairIndexTable

Lightweight upper-triangular view over retained-unit pairs. It stores retained
units once and computes pair records on iteration for compatibility with older
pair consumers that still expect `UnitPairRecord`s.
"""
struct UnitPairIndexTable
    units::Vector{CartesianRetainedUnits.RetainedUnitRecord}
    metadata::NamedTuple
end

function unit_pair_index_table(
    retained_unit_plan::CartesianRetainedUnits.RetainedUnitPlan;
    metadata = (;),
)
    return unit_pair_index_table(
        CartesianRetainedUnits.units(retained_unit_plan);
        metadata,
    )
end

function unit_pair_index_table(retained_units; metadata = (;))
    unit_vector =
        retained_units isa Vector{CartesianRetainedUnits.RetainedUnitRecord} ?
        retained_units :
        CartesianRetainedUnits.RetainedUnitRecord[unit for unit in retained_units]
    return UnitPairIndexTable(unit_vector, NamedTuple(metadata))
end

Base.length(table::UnitPairIndexTable) =
    div(length(table.units) * (length(table.units) + 1), 2)
Base.isempty(table::UnitPairIndexTable) = isempty(table.units)
Base.IteratorSize(::Type{UnitPairIndexTable}) = Base.HasLength()
Base.eltype(::Type{UnitPairIndexTable}) = UnitPairRecord

function Base.iterate(table::UnitPairIndexTable)
    isempty(table) && return nothing
    return _unit_pair_index_table_iterate(table, 1, 1, 1)
end

function Base.iterate(table::UnitPairIndexTable, state)
    left_index, right_index, pair_index = state
    left_index > length(table.units) && return nothing
    return _unit_pair_index_table_iterate(
        table,
        left_index,
        right_index,
        pair_index,
    )
end

function _unit_pair_index_table_iterate(
    table::UnitPairIndexTable,
    left_index::Int,
    right_index::Int,
    pair_index::Int,
)
    pair = _unit_pair_record_from_index_table(
        table,
        left_index,
        right_index,
        pair_index,
    )
    next_right = right_index + 1
    next_left = left_index
    if next_right > length(table.units)
        next_left += 1
        next_right = next_left
    end
    return pair, (next_left, next_right, pair_index + 1)
end

function _unit_pair_record_from_index_table(
    table::UnitPairIndexTable,
    left_index::Int,
    right_index::Int,
    pair_index::Int,
)
    left = table.units[left_index]
    right = table.units[right_index]
    metadata = merge(
        table.metadata,
        (;
            retained_pair_source = :upper_triangular_unit_index_table,
            rich_unit_pair_record_stored = false,
        ),
    )
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
        nothing,
        false,
        metadata,
    )
end

"""
    UnitPairPlan

Metadata-only pair inventory for one retained-unit plan.
"""
struct UnitPairPlan
    policy::UnitPairPolicy
    retained_unit_plan::CartesianRetainedUnits.RetainedUnitPlan
    pairs::Union{UnitPairIndexTable,Vector{UnitPairRecord}}
    route_core_pair_inventory::Union{CartesianRouteCore.UnitPairInventory,Nothing}
    summary::NamedTuple
    metadata::NamedTuple

    function UnitPairPlan(
        policy::UnitPairPolicy,
        retained_unit_plan::CartesianRetainedUnits.RetainedUnitPlan,
        pairs,
        route_core_pair_inventory,
        summary,
        metadata,
    )
        pair_storage =
            pairs isa UnitPairIndexTable ? pairs :
            pairs isa Vector{UnitPairRecord} ? pairs :
            UnitPairRecord[pair for pair in pairs]
        return new(
            policy,
            retained_unit_plan,
            pair_storage,
            route_core_pair_inventory,
            NamedTuple(summary),
            NamedTuple(metadata),
        )
    end
end

unit_pairs(plan::UnitPairPlan) = plan.pairs
unit_pairs(table::UnitPairIndexTable) = table
summary(plan::UnitPairPlan) = plan.summary
route_core_pair_inventory(plan::UnitPairPlan) = plan.route_core_pair_inventory

function _pair_family(
    left::CartesianRetainedUnits.RetainedUnitRecord,
    right::CartesianRetainedUnits.RetainedUnitRecord,
)
    return Symbol(String(left.unit_kind), "__", String(right.unit_kind))
end

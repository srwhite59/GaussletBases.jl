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
    left_unit::CRU.RetainedUnitRecord
    right_unit::CRU.RetainedUnitRecord
    left_index::Int
    right_index::Int
    left_unit_key::Symbol
    right_unit_key::Symbol
    left_unit_kind::Symbol
    right_unit_kind::Symbol
    route_core_pair_sidecar::Union{CRC.UnitPair,Nothing}
    materialized::Bool
    metadata::NamedTuple
end

"""
    UnitPairPlan

Metadata-only pair inventory for one retained-unit plan.
"""
struct UnitPairPlan
    policy::UnitPairPolicy
    retained_unit_plan::CRU.RetainedUnitPlan
    pairs::Tuple{Vararg{UnitPairRecord}}
    route_core_pair_inventory::Union{CRC.UnitPairInventory,Nothing}
    summary::NamedTuple
    metadata::NamedTuple
end

unit_pairs(plan::UnitPairPlan) = plan.pairs
summary(plan::UnitPairPlan) = plan.summary
route_core_pair_inventory(plan::UnitPairPlan) = plan.route_core_pair_inventory

function _pair_family(left::CRU.RetainedUnitRecord, right::CRU.RetainedUnitRecord)
    return Symbol(String(left.unit_kind), "__", String(right.unit_kind))
end

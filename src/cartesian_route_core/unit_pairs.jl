# Unit-pair inventory. Pair planning starts from final retained units, not from
# shellification regions or CPBs directly.

"""
    UnitPair

Upper-triangular pair entry connecting two `FinalRetainedUnit` objects for
later operator-block planning.
"""
struct UnitPair
    pair_key::Tuple{Symbol,Symbol}
    left::FinalRetainedUnit
    right::FinalRetainedUnit
    left_index::Int
    right_index::Int
    symmetry::Symbol
    metadata::NamedTuple
end

"""
    UnitPairInventory

Bookkeeping inventory of final retained units and their upper-triangular pair
entries.
"""
struct UnitPairInventory
    units::Tuple{Vararg{FinalRetainedUnit}}
    pairs::Tuple{Vararg{UnitPair}}
    symmetry::Symbol
    metadata::NamedTuple
end

"""
    unit_pair_inventory(units; symmetry = :symmetric_upper_triangle, metadata = (;))

Build upper-triangular unit pairs from final retained units.

This is bookkeeping for later operator-block planning. The input must be
`FinalRetainedUnit` objects, not CPBs or shellification regions.
"""
function unit_pair_inventory(units; symmetry::Symbol = :symmetric_upper_triangle, metadata = (;))
    unit_tuple = Tuple(units)
    isempty(unit_tuple) && throw(ArgumentError("unit_pair_inventory requires at least one FinalRetainedUnit"))
    all(unit -> unit isa FinalRetainedUnit, unit_tuple) ||
        throw(ArgumentError("unit_pair_inventory accepts FinalRetainedUnit objects only"))

    keys = Tuple(unit.unit_key for unit in unit_tuple)
    length(unique(keys)) == length(keys) ||
        throw(ArgumentError("FinalRetainedUnit unit_key values must be unique"))

    pairs = UnitPair[]
    for left_index in eachindex(unit_tuple)
        for right_index in left_index:length(unit_tuple)
            left = unit_tuple[left_index]
            right = unit_tuple[right_index]
            push!(
                pairs,
                UnitPair(
                    (left.unit_key, right.unit_key),
                    left,
                    right,
                    left_index,
                    right_index,
                    symmetry,
                    NamedTuple(metadata),
                ),
            )
        end
    end

    return UnitPairInventory(unit_tuple, Tuple(pairs), symmetry, NamedTuple(metadata))
end

"""
    final_units(inventory)

Return the final retained units used to form a unit-pair inventory.
"""
final_units(inventory::UnitPairInventory) = inventory.units

"""
    pair_entries(inventory)

Return the upper-triangular `UnitPair` entries in a unit-pair inventory.
"""
pair_entries(inventory::UnitPairInventory) = inventory.pairs

"""
    unit_keys(inventory)

Return final retained unit keys in inventory order.
"""
unit_keys(inventory::UnitPairInventory) = Tuple(unit.unit_key for unit in inventory.units)

"""
    pair_keys(inventory)

Return pair keys in inventory order.
"""
pair_keys(inventory::UnitPairInventory) = Tuple(pair.pair_key for pair in inventory.pairs)

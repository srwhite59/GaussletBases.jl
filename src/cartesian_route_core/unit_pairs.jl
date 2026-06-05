# Unit-pair inventory. Pair planning starts from final retained units, not from
# shellification regions or CPBs directly.

struct UnitPair
    pair_key::Tuple{Symbol,Symbol}
    left::FinalRetainedUnit
    right::FinalRetainedUnit
    left_index::Int
    right_index::Int
    symmetry::Symbol
    metadata::NamedTuple
end

struct UnitPairInventory
    units::Tuple{Vararg{FinalRetainedUnit}}
    pairs::Tuple{Vararg{UnitPair}}
    symmetry::Symbol
    metadata::NamedTuple
end

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

final_units(inventory::UnitPairInventory) = inventory.units
pair_entries(inventory::UnitPairInventory) = inventory.pairs
unit_keys(inventory::UnitPairInventory) = Tuple(unit.unit_key for unit in inventory.units)
pair_keys(inventory::UnitPairInventory) = Tuple(pair.pair_key for pair in inventory.pairs)

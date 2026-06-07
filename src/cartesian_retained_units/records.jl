# Retained-unit planning records.

"""
    RetainedUnitPolicy

Policy marker for converting selected terminal-lowering contracts into planned
retained-unit metadata.
"""
abstract type RetainedUnitPolicy end

"""
    MetadataOnlyRetainedUnits()

Default retained-unit policy. It creates metadata-only retained-unit records
and RouteCore sidecars where possible, without materializing matrices or
column ranges.
"""
struct MetadataOnlyRetainedUnits <: RetainedUnitPolicy
end

policy_kind(::MetadataOnlyRetainedUnits) = :metadata_only_retained_units

"""
    RetainedUnitRecord

Metadata-only retained unit implied by one selected terminal-lowering contract.
This is the first layer where selected lowering contracts become column-owning
unit records. Dimensions and column ranges remain unavailable until later
construction materializes them.
"""
struct RetainedUnitRecord
    unit_key::Symbol
    unit_index::Int
    unit_kind::Symbol
    source_contract_key::Symbol
    terminal_region_key::Symbol
    terminal_region_role::Symbol
    terminal_region_kind::Symbol
    lowering_kind::Symbol
    retained_rule::Symbol
    realization_rule::Union{Symbol,Nothing}
    owned_support::Any
    source_cpbs::Tuple
    source_cpb_index::Union{Int,Nothing}
    dimension_status::Symbol
    dimension::Union{Int,Nothing}
    column_range_status::Symbol
    column_range::Union{UnitRange{Int},Nothing}
    route_core_final_unit::Union{CRC.FinalRetainedUnit,Nothing}
    materialized::Bool
    metadata::NamedTuple
end

"""
    RetainedUnitPlan

Metadata-only retained-unit plan for one terminal-lowering plan.
"""
struct RetainedUnitPlan
    policy::RetainedUnitPolicy
    lowering_plan::CTL.TerminalLoweringPlan
    units::Tuple{Vararg{RetainedUnitRecord}}
    summary::NamedTuple
    metadata::NamedTuple
end

units(plan::RetainedUnitPlan) = plan.units
summary(plan::RetainedUnitPlan) = plan.summary

function route_core_final_units(plan::RetainedUnitPlan)
    return Tuple(
        unit.route_core_final_unit
        for unit in plan.units
        if !isnothing(unit.route_core_final_unit)
    )
end

function _metadata_value(metadata::NamedTuple, key::Symbol, default = nothing)
    return haskey(metadata, key) ? getfield(metadata, key) : default
end

function _merge_metadata(parts...)
    merged = (;)
    for part in parts
        merged = merge(merged, NamedTuple(part))
    end
    return merged
end

function _unit_key(contract::CTL.TerminalLoweringContract, suffix::Symbol)
    return Symbol(String(contract.contract_key), "_", String(suffix))
end

function _unit_key(contract::CTL.TerminalLoweringContract, suffix::Symbol, index::Int)
    return Symbol(String(contract.contract_key), "_", String(suffix), "_", index)
end

function _stratum_kind(cpb::CPB.CoordinateProductBox)
    return _metadata_value(cpb.metadata, :stratum_kind, :unknown_cpb)
end

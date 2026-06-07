# Retained-unit transform-contract records.

"""
    RetainedUnitTransformContractPolicy

Policy marker for converting retained-unit records into planned
transform-contract metadata.
"""
abstract type RetainedUnitTransformContractPolicy end

"""
    MetadataOnlyRetainedUnitTransformContracts()

Default transform-contract policy. It creates one metadata-only transform
contract per retained unit and does not materialize transforms or coefficient
maps.
"""
struct MetadataOnlyRetainedUnitTransformContracts <: RetainedUnitTransformContractPolicy end

policy_kind(::MetadataOnlyRetainedUnitTransformContracts) =
    :metadata_only_retained_unit_transform_contracts

"""
    RetainedUnitTransformContract

Metadata-only transform contract for one retained unit.
"""
struct RetainedUnitTransformContract
    unit_key::Symbol
    unit_index::Int
    unit_kind::Symbol
    lowering_kind::Symbol
    retained_rule::Symbol
    realization_rule::Union{Symbol,Nothing}
    source_cpbs::Tuple
    transform_path::Symbol
    realization_path::Symbol
    dimension_status::Symbol
    column_range_status::Symbol
    materialized::Bool
    blocker::Union{Symbol,String,Nothing}
    metadata::NamedTuple
end

"""
    RetainedUnitTransformContractPlan

Metadata-only transform-contract plan for one retained-unit plan.
"""
struct RetainedUnitTransformContractPlan
    policy::RetainedUnitTransformContractPolicy
    retained_unit_plan::CRU.RetainedUnitPlan
    contracts::Tuple{Vararg{RetainedUnitTransformContract}}
    summary::NamedTuple
    metadata::NamedTuple
end

transform_contracts(plan::RetainedUnitTransformContractPlan) = plan.contracts
summary(plan::RetainedUnitTransformContractPlan) = plan.summary

function _merge_metadata(parts...)
    merged = (;)
    for part in parts
        merged = merge(merged, NamedTuple(part))
    end
    return merged
end

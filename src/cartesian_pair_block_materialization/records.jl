# Metadata-only pair-block materialization records.

"""
    PairBlockMaterializationPolicy

Policy marker for converting pair-operator plans into pair-block
materialization readiness records.
"""
abstract type PairBlockMaterializationPolicy end

"""
    MetadataOnlyPairBlockMaterialization()

Default pair-block materialization policy. It records readiness for future
numerical block construction without building matrices or operator blocks.
"""
struct MetadataOnlyPairBlockMaterialization <: PairBlockMaterializationPolicy end

policy_kind(::MetadataOnlyPairBlockMaterialization) =
    :metadata_only_pair_block_materialization

"""
    PairBlockMaterializationRecord

Metadata-only readiness record for one pair-operator plan record.
"""
struct PairBlockMaterializationRecord
    pair_key::Tuple{Symbol,Symbol}
    pair_index::Int
    pair_family::Symbol
    source_operator_path::Symbol
    transform_path::NamedTuple
    realization_path::NamedTuple
    final_block_path::Symbol
    supported_terms::Tuple{Vararg{Symbol}}
    materialization_path::Symbol
    readiness_status::Symbol
    blocker::Union{Symbol,String,Nothing}
    materialized::Bool
    metadata::NamedTuple
end

"""
    PairBlockMaterializationPlan

Metadata-only pair-block materialization readiness plan for one pair-operator
plan.
"""
struct PairBlockMaterializationPlan
    policy::PairBlockMaterializationPolicy
    pair_operator_plan::CPOP.PairOperatorPlan
    records::Tuple{Vararg{PairBlockMaterializationRecord}}
    summary::NamedTuple
    metadata::NamedTuple
end

pair_block_materialization_records(plan::PairBlockMaterializationPlan) = plan.records
summary(plan::PairBlockMaterializationPlan) = plan.summary

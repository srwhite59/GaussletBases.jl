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
    PairBlockMaterializationResult

Numerical result for one materialized pair block. This is still local pair-block
data, not a global operator or Hamiltonian assembly.
"""
struct PairBlockMaterializationResult
    term::Symbol
    pair_key::Tuple{Symbol,Symbol}
    block::Matrix{Float64}
    materialized::Bool
    source_operator_blocks_materialized::Bool
    final_pair_blocks_materialized::Bool
    operator_blocks_materialized::Bool
    hamiltonian_data_materialized::Bool
    artifacts_materialized::Bool
    metadata::NamedTuple
end

"""
    PairBlockMaterializationBatchResult

Compact result for a narrow plan-level materialization pass. It carries local
pair-block results and skipped-record summaries only.
"""
struct PairBlockMaterializationBatchResult
    term::Symbol
    materialized_results::Tuple{Vararg{PairBlockMaterializationResult}}
    skipped_records::Tuple{Vararg{NamedTuple}}
    materialized_count::Int
    skipped_count::Int
    materialized::Bool
    source_operator_blocks_materialized::Bool
    final_pair_blocks_materialized::Bool
    operator_blocks_materialized::Bool
    hamiltonian_data_materialized::Bool
    artifacts_materialized::Bool
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

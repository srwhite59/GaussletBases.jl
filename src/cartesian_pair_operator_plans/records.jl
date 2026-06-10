# Metadata-only retained-unit pair-operator planning records.

"""
    PairOperatorPlanPolicy

Policy marker for converting retained-unit pair plans into pair-operator plan
metadata.
"""
abstract type PairOperatorPlanPolicy end

"""
    MetadataOnlyPairOperatorPlans()

Default pair-operator planning policy. It classifies source-operator,
realization, and final-block paths without materializing numerical blocks.
"""
struct MetadataOnlyPairOperatorPlans <: PairOperatorPlanPolicy
end

policy_kind(::MetadataOnlyPairOperatorPlans) = :metadata_only_pair_operator_plans

"""
    PairOperatorPlanRecord

Metadata-only operator plan for one retained-unit pair.
"""
struct PairOperatorPlanRecord
    pair_key::Tuple{Symbol,Symbol}
    pair_index::Int
    pair_family::Symbol
    left_unit_key::Symbol
    right_unit_key::Symbol
    left_unit_kind::Symbol
    right_unit_kind::Symbol
    source_operator_path::Symbol
    transform_path::NamedTuple
    realization_path::NamedTuple
    final_block_path::Symbol
    supported_terms::Tuple{Vararg{Symbol}}
    status::Symbol
    blocker::Union{Symbol,String,Nothing}
    route_core_pair_operator_sidecar::Union{CartesianRouteCore.PairOperatorPlan,Nothing}
    materialized::Bool
    metadata::NamedTuple
end

"""
    PairOperatorPlan

Metadata-only pair-operator plan inventory for one `CartesianUnitPairs`
unit-pair plan.
"""
struct PairOperatorPlan
    policy::PairOperatorPlanPolicy
    unit_pair_plan::CartesianUnitPairs.UnitPairPlan
    transform_contract_plan::CartesianRetainedUnitTransformContracts.RetainedUnitTransformContractPlan
    records::Tuple{Vararg{PairOperatorPlanRecord}}
    route_core_pair_operator_plan_inventory::Union{
        CartesianRouteCore.PairOperatorPlanInventory,
        Nothing,
    }
    summary::NamedTuple
    metadata::NamedTuple
end

pair_operator_records(plan::PairOperatorPlan) = plan.records
summary(plan::PairOperatorPlan) = plan.summary
route_core_pair_operator_plan_inventory(plan::PairOperatorPlan) =
    plan.route_core_pair_operator_plan_inventory

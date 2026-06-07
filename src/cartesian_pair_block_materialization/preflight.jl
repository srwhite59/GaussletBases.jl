# Pair-operator plan to metadata-only pair-block materialization readiness.

"""
    pair_block_materialization_plan(pair_operator_plan; policy)

Build metadata-only pair-block materialization readiness records from
pair-operator plan metadata. This does not build numerical blocks.
"""
function pair_block_materialization_plan(
    pair_operator_plan::CPOP.PairOperatorPlan;
    policy::PairBlockMaterializationPolicy = MetadataOnlyPairBlockMaterialization(),
    metadata = (;),
)
    unit_pairs = CUP.unit_pairs(pair_operator_plan.unit_pair_plan)
    pair_operator_records = CPOP.pair_operator_records(pair_operator_plan)
    records = Tuple(
        _pair_block_materialization_record(record, unit_pairs[index], policy)
        for (index, record) in pairs(pair_operator_records)
    )
    plan_summary =
        _pair_block_materialization_plan_summary(policy, pair_operator_plan, records)
    return PairBlockMaterializationPlan(
        policy,
        pair_operator_plan,
        records,
        plan_summary,
        NamedTuple(metadata),
    )
end

function _pair_block_materialization_record(
    record::CPOP.PairOperatorPlanRecord,
    unit_pair::CUP.UnitPairRecord,
    ::MetadataOnlyPairBlockMaterialization,
)
    materialization_path, readiness_status, blocker =
        _pair_block_materialization_preflight(record)
    return PairBlockMaterializationRecord(
        record.pair_key,
        record.pair_index,
        record.pair_family,
        record.source_operator_path,
        record.transform_path,
        record.realization_path,
        record.final_block_path,
        record.supported_terms,
        materialization_path,
        readiness_status,
        blocker,
        false,
        _pair_block_materialization_record_metadata(record, unit_pair),
    )
end

function _pair_block_materialization_preflight(record::CPOP.PairOperatorPlanRecord)
    if !isnothing(record.blocker)
        return (
            :blocked_pair_block_materialization_path,
            :blocked_pair_operator_plan_not_ready,
            record.blocker,
        )
    end

    _is_direct_direct_pair_block_pilot(record) && return (
        :direct_direct_pair_block_materialization_pilot,
        :ready_metadata_only_not_materialized,
        nothing,
    )

    return (
        :deferred_pair_block_materialization_path,
        :blocked_pair_block_materialization_not_implemented,
        :non_direct_direct_pair_block_materialization_not_implemented,
    )
end

function _is_direct_direct_pair_block_pilot(record::CPOP.PairOperatorPlanRecord)
    return record.source_operator_path === :direct_identity_cpb_path &&
           record.transform_path.left === :direct_identity_transform_contract &&
           record.transform_path.right === :direct_identity_transform_contract &&
           record.realization_path.left === :identity_or_trivial_embedding &&
           record.realization_path.right === :identity_or_trivial_embedding &&
           record.final_block_path === :source_block_direct_to_final_block
end

function _pair_block_materialization_record_metadata(
    record::CPOP.PairOperatorPlanRecord,
    unit_pair::CUP.UnitPairRecord,
)
    base_metadata = (;
        pair_operator_status = record.status,
        pair_operator_blocker = record.blocker,
    )
    _is_direct_direct_pair_block_pilot(record) || return base_metadata

    return merge(
        base_metadata,
        (;
            direct_source_metadata_status = :available_from_retained_unit_source_cpbs,
            left_source_cpbs = unit_pair.left_unit.source_cpbs,
            right_source_cpbs = unit_pair.right_unit.source_cpbs,
        ),
    )
end

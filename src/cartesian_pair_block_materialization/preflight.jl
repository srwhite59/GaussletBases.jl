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
    transform_lookup =
        _pair_block_materialization_transform_contract_lookup(
            pair_operator_plan.transform_contract_plan,
        )
    records = Tuple(
        _pair_block_materialization_record(
            record,
            unit_pairs[index],
            transform_lookup,
            policy,
        )
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
    transform_lookup,
    ::MetadataOnlyPairBlockMaterialization,
)
    materialization_path, readiness_status, blocker =
        _pair_block_materialization_preflight(record, unit_pair, transform_lookup)
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
        _pair_block_materialization_record_metadata(record, unit_pair, transform_lookup),
    )
end

function _pair_block_materialization_transform_contract_lookup(
    transform_contract_plan::CRTC.RetainedUnitTransformContractPlan,
)
    contracts_by_key = Dict{Symbol,CRTC.RetainedUnitTransformContract}()
    for contract in CRTC.transform_contracts(transform_contract_plan)
        contracts_by_key[contract.unit_key] = contract
    end
    return contracts_by_key
end

function _pair_block_transform_contract(transform_lookup, unit_key::Symbol)
    return get(transform_lookup, unit_key, nothing)
end

function _pair_block_materialization_preflight(
    record::CPOP.PairOperatorPlanRecord,
    unit_pair::CUP.UnitPairRecord,
    transform_lookup,
)
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

    if _is_pqs_pqs_source_pair_preflight(record)
        readiness_status, blocker =
            _pqs_source_pair_preflight_status(record, unit_pair, transform_lookup)
        return (
            :pqs_source_pair_preflight,
            readiness_status,
            blocker,
        )
    end

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

function _is_pqs_pqs_source_pair_preflight(record::CPOP.PairOperatorPlanRecord)
    return record.left_unit_kind === :pqs_shell_retained_unit &&
           record.right_unit_kind === :pqs_shell_retained_unit &&
           record.source_operator_path === :pqs_source_cpb_1d_factor_path
end

function _pqs_source_pair_preflight_status(
    _record::CPOP.PairOperatorPlanRecord,
    unit_pair::CUP.UnitPairRecord,
    transform_lookup,
)
    left_contract =
        _pair_block_transform_contract(transform_lookup, unit_pair.left_unit_key)
    right_contract =
        _pair_block_transform_contract(transform_lookup, unit_pair.right_unit_key)

    _raw_product_source_plan_available(left_contract) || return (
        :blocked_missing_raw_product_source_plan,
        :missing_left_raw_product_source_plan,
    )
    _raw_product_source_plan_available(right_contract) || return (
        :blocked_missing_raw_product_source_plan,
        :missing_right_raw_product_source_plan,
    )

    left_summary = _raw_product_source_summary(left_contract)
    right_summary = _raw_product_source_summary(right_contract)
    _raw_product_source_facts_complete(left_summary) || return (
        :blocked_missing_raw_product_source_plan,
        :missing_left_raw_product_source_plan,
    )
    _raw_product_source_facts_complete(right_summary) || return (
        :blocked_missing_raw_product_source_plan,
        :missing_right_raw_product_source_plan,
    )

    left_summary.source_mode_ordering == right_summary.source_mode_ordering || return (
        :blocked_incompatible_raw_product_source_ordering,
        :incompatible_raw_product_source_ordering,
    )

    return :ready_metadata_only_not_materialized, nothing
end

function _raw_product_source_summary(contract)
    isnothing(contract) && return nothing
    haskey(contract.metadata, :raw_product_source_summary) ||
        return nothing
    return contract.metadata.raw_product_source_summary
end

function _raw_product_source_plan_status(contract)
    isnothing(contract) && return nothing
    haskey(contract.metadata, :raw_product_source_plan_status) ||
        return nothing
    return contract.metadata.raw_product_source_plan_status
end

function _raw_product_source_plan_available(contract)
    return _raw_product_source_plan_status(contract) === :available_raw_product_box_plan &&
           haskey(contract.metadata, :raw_product_source_plan) &&
           !isnothing(contract.metadata.raw_product_source_plan) &&
           !isnothing(_raw_product_source_summary(contract))
end

function _raw_product_source_facts_complete(summary)
    isnothing(summary) && return false
    for field in (:source_mode_dims, :source_mode_count, :source_mode_ordering)
        haskey(summary, field) || return false
        isnothing(getfield(summary, field)) && return false
    end
    return true
end

function _transform_contract_key(contract)
    isnothing(contract) && return nothing
    return contract.unit_key
end

function _source_contract_key(contract)
    isnothing(contract) && return nothing
    haskey(contract.metadata, :source_contract_key) ||
        return nothing
    return contract.metadata.source_contract_key
end

function _pair_block_materialization_record_metadata(
    record::CPOP.PairOperatorPlanRecord,
    unit_pair::CUP.UnitPairRecord,
    transform_lookup,
)
    base_metadata = (;
        pair_operator_status = record.status,
        pair_operator_blocker = record.blocker,
    )
    if _is_direct_direct_pair_block_pilot(record)
        return merge(
            base_metadata,
            (;
                direct_source_metadata_status = :available_from_retained_unit_source_cpbs,
                left_source_cpbs = unit_pair.left_unit.source_cpbs,
                right_source_cpbs = unit_pair.right_unit.source_cpbs,
            ),
        )
    end

    _is_pqs_pqs_source_pair_preflight(record) ||
        return base_metadata

    left_contract =
        _pair_block_transform_contract(transform_lookup, unit_pair.left_unit_key)
    right_contract =
        _pair_block_transform_contract(transform_lookup, unit_pair.right_unit_key)
    left_summary = _raw_product_source_summary(left_contract)
    right_summary = _raw_product_source_summary(right_contract)
    left_ordering =
        _raw_product_source_facts_complete(left_summary) ?
        left_summary.source_mode_ordering :
        nothing
    right_ordering =
        _raw_product_source_facts_complete(right_summary) ?
        right_summary.source_mode_ordering :
        nothing
    source_mode_ordering =
        left_ordering == right_ordering ? left_ordering : nothing
    return merge(
        base_metadata,
        (;
            transform_contract_keys = (;
                left = _transform_contract_key(left_contract),
                right = _transform_contract_key(right_contract),
            ),
            source_contract_keys = (;
                left = _source_contract_key(left_contract),
                right = _source_contract_key(right_contract),
            ),
            left_raw_product_source_plan_status =
                _raw_product_source_plan_status(left_contract),
            right_raw_product_source_plan_status =
                _raw_product_source_plan_status(right_contract),
            left_raw_product_source_summary = left_summary,
            right_raw_product_source_summary = right_summary,
            left_source_mode_dims =
                _raw_product_source_facts_complete(left_summary) ?
                left_summary.source_mode_dims :
                nothing,
            right_source_mode_dims =
                _raw_product_source_facts_complete(right_summary) ?
                right_summary.source_mode_dims :
                nothing,
            left_source_mode_count =
                _raw_product_source_facts_complete(left_summary) ?
                left_summary.source_mode_count :
                nothing,
            right_source_mode_count =
                _raw_product_source_facts_complete(right_summary) ?
                right_summary.source_mode_count :
                nothing,
            left_source_mode_ordering = left_ordering,
            right_source_mode_ordering = right_ordering,
            source_mode_ordering,
            source_operator_blocks_materialized = false,
            final_pair_blocks_materialized = false,
        ),
    )
end

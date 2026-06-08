# Metadata-only mixed one-body dispatch classifier.
#
# This skeleton classifies which existing family selector would be used later
# and whether caller-supplied inputs are present. It deliberately does not call
# numerical family selectors.

function _one_body_pair_block_dispatch_summary(
    record::PairBlockMaterializationRecord,
    term::Symbol;
    inputs = (;),
    provider = nothing,
    unit_pair = nothing,
)
    term_descriptor = _one_body_term_descriptor(term)
    selector_family = _one_body_record_selector_family(record)
    parent_axis_counts_required =
        _one_body_dispatch_parent_axis_counts_required(selector_family)
    factor_input_summary = _one_body_factor_input_summary(
        term;
        inputs,
        provider,
        selector_family =
            selector_family === :unsupported ? nothing : selector_family,
        parent_axis_counts_required,
    )
    term_supported_by_record = term in record.supported_terms
    record_ready = record.readiness_status === :ready_metadata_only_not_materialized &&
                   isnothing(record.blocker)
    unit_pair_requirement =
        _one_body_unit_pair_requirement(record, selector_family, unit_pair)
    blockers = _one_body_dispatch_blockers(
        selector_family,
        record_ready,
        term_supported_by_record,
        factor_input_summary,
        unit_pair_requirement,
    )
    dispatch_ready = isempty(blockers)

    return (;
        object_kind = :cartesian_pair_block_mixed_one_body_dispatch_summary,
        status =
            dispatch_ready ?
            :ready_metadata_only_not_materialized :
            :blocked_mixed_one_body_dispatch,
        blocker = dispatch_ready ? nothing : first(blockers),
        blockers,
        pair_key = record.pair_key,
        pair_index = record.pair_index,
        pair_family = record.pair_family,
        selector_family,
        selector_family_status =
            selector_family === :unsupported ? :unsupported : :available,
        requested_term = term,
        term_descriptor,
        term_supported_by_record,
        record_materialization_path = record.materialization_path,
        materialization_path = factor_input_summary.materialization_path,
        materialization_path_status =
            factor_input_summary.materialization_path_status,
        record_readiness_status = record.readiness_status,
        record_blocker = record.blocker,
        record_ready,
        factor_input_status = factor_input_summary.status,
        factor_input_blocker = factor_input_summary.blocker,
        factor_input_blockers = factor_input_summary.blockers,
        parent_axis_counts_required =
            factor_input_summary.parent_axis_counts_required,
        parent_axis_counts_status =
            factor_input_summary.parent_axis_counts_status,
        required_factor_names = factor_input_summary.required_factor_names,
        present_factor_names = factor_input_summary.present_factor_names,
        missing_factor_names = factor_input_summary.missing_factor_names,
        required_factors_available =
            factor_input_summary.required_factors_available,
        unit_pair_required = unit_pair_requirement.required,
        unit_pair_available = unit_pair_requirement.available,
        unit_pair_key = unit_pair_requirement.pair_key,
        unit_pair_requirement_status = unit_pair_requirement.status,
        unit_pair_blocker = unit_pair_requirement.blocker,
        factors_constructed = false,
        numerical_blocks_materialized = false,
        source_operator_blocks_materialized = false,
        final_pair_blocks_materialized = false,
        operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
        mixed_dispatcher_materialized = false,
        route_driver_wiring = false,
        pqs_lowdin_materialized = false,
        full_white_lindsey_route_assembled = false,
    )
end

function _one_body_pair_block_dispatch_summary(record, term; kwargs...)
    throw(
        ArgumentError(
            "one-body dispatch summary requires a PairBlockMaterializationRecord",
        ),
    )
end

function _one_body_pair_block_plan_dispatch_summary(
    plan::PairBlockMaterializationPlan,
    term::Symbol;
    inputs = (;),
    provider = nothing,
)
    unit_pair_lookup = _one_body_unit_pair_lookup(plan)
    record_summaries = Tuple(
        _one_body_pair_block_dispatch_summary(
            record,
            term;
            inputs,
            provider,
            unit_pair = _one_body_dispatch_unit_pair(record, unit_pair_lookup),
        ) for record in pair_block_materialization_records(plan)
    )
    ready_count = count(
        summary -> summary.status === :ready_metadata_only_not_materialized,
        record_summaries,
    )
    record_count = length(record_summaries)
    blocked_count = record_count - ready_count

    return (;
        object_kind =
            :cartesian_pair_block_mixed_one_body_plan_dispatch_summary,
        status =
            blocked_count == 0 ?
            :ready_metadata_only_not_materialized :
            ready_count == 0 ?
            :blocked_mixed_one_body_plan_dispatch :
            :partially_ready_mixed_one_body_plan_dispatch,
        blocker =
            blocked_count == 0 ?
            nothing :
            :blocked_mixed_one_body_dispatch_records,
        requested_term = term,
        term_descriptor = _one_body_term_descriptor(term),
        record_count,
        ready_count,
        blocked_count,
        selector_family_counts =
            _one_body_count_by(record_summaries, :selector_family),
        dispatch_status_counts = _one_body_count_by(record_summaries, :status),
        blocker_counts = _one_body_blocker_counts(record_summaries),
        record_materialization_path_counts =
            _one_body_count_by(record_summaries, :record_materialization_path),
        materialization_path_counts =
            _one_body_count_by(record_summaries, :materialization_path),
        white_lindsey_unit_pair_required_count = count(
            summary -> summary.unit_pair_required,
            record_summaries,
        ),
        white_lindsey_unit_pair_available_count = count(
            summary -> summary.unit_pair_required && summary.unit_pair_available,
            record_summaries,
        ),
        record_summaries,
        factors_constructed = false,
        numerical_blocks_materialized = false,
        source_operator_blocks_materialized = false,
        final_pair_blocks_materialized = false,
        operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
        mixed_dispatcher_materialized = false,
        route_driver_wiring = false,
        pqs_lowdin_materialized = false,
        full_white_lindsey_route_assembled = false,
    )
end

function _one_body_pair_block_plan_dispatch_summary(plan, term; kwargs...)
    throw(
        ArgumentError(
            "one-body plan dispatch summary requires a PairBlockMaterializationPlan",
        ),
    )
end

function _one_body_pair_block_set_preflight_summary(
    plan::PairBlockMaterializationPlan;
    terms = _ONE_BODY_TERMS,
    inputs = (;),
    provider = nothing,
    parent_axis_counts_required =
        _one_body_block_set_parent_axis_counts_required(plan),
)
    term_set_descriptor = _one_body_term_set_descriptor(terms)
    input_summary = _one_body_term_set_factor_input_summary(
        term_set_descriptor.terms;
        inputs,
        provider,
        parent_axis_counts_required,
    )
    term_summaries = Tuple(
        _one_body_pair_block_set_preflight_term_summary(
            _one_body_pair_block_plan_dispatch_summary(
                plan,
                term;
                inputs,
                provider,
            ),
        ) for term in term_set_descriptor.terms
    )
    total_ready_record_count = sum(summary.ready_count for summary in term_summaries)
    total_blocked_record_count =
        sum(summary.blocked_count for summary in term_summaries)

    return (;
        object_kind =
            :cartesian_pair_block_mixed_one_body_block_set_preflight_summary,
        status = _one_body_block_set_preflight_status(
            input_summary,
            total_ready_record_count,
            total_blocked_record_count,
        ),
        blocker = _one_body_block_set_preflight_blocker(
            input_summary,
            total_blocked_record_count,
        ),
        term_set_descriptor,
        input_summary,
        requested_terms = term_set_descriptor.terms,
        terms = term_set_descriptor.terms,
        term_count = term_set_descriptor.term_count,
        plan_record_count = length(pair_block_materialization_records(plan)),
        term_set_input_status = input_summary.status,
        term_set_input_blocker = input_summary.blocker,
        term_set_input_blockers = input_summary.blockers,
        input_source = input_summary.input_source,
        factor_provider_scope = input_summary.factor_provider_scope,
        parent_axis_counts_required = input_summary.parent_axis_counts_required,
        parent_axis_counts_status = input_summary.parent_axis_counts_status,
        required_factor_names = input_summary.required_factor_names,
        present_factor_names = input_summary.present_factor_names,
        missing_factor_names = input_summary.missing_factor_names,
        term_statuses = Tuple(
            (;
                term = term_summary.term,
                status = term_summary.status,
                ready_count = term_summary.ready_count,
                blocked_count = term_summary.blocked_count,
            ) for term_summary in term_summaries
        ),
        term_status_counts = _one_body_count_by(term_summaries, :status),
        total_ready_record_count,
        total_blocked_record_count,
        selector_family_counts = _one_body_block_set_aggregate_counts(
            term_summaries,
            :selector_family_counts,
            :selector_family,
        ),
        dispatch_status_counts = _one_body_block_set_aggregate_counts(
            term_summaries,
            :dispatch_status_counts,
            :status,
        ),
        record_materialization_path_counts =
            _one_body_block_set_aggregate_counts(
                term_summaries,
                :record_materialization_path_counts,
                :record_materialization_path,
            ),
        materialization_path_counts = _one_body_block_set_aggregate_counts(
            term_summaries,
            :materialization_path_counts,
            :materialization_path,
        ),
        blocker_counts = _one_body_block_set_aggregate_counts(
            term_summaries,
            :blocker_counts,
            :blocker,
        ),
        term_summaries,
        result_terms_remain_separated = true,
        block_set_results_summed = false,
        factor_values_stored = false,
        factors_constructed = false,
        numerical_blocks_materialized = false,
        materialized = false,
        source_operator_blocks_materialized = false,
        final_pair_blocks_materialized = false,
        operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
        global_operator_blocks_materialized = false,
        global_hamiltonian_data_materialized = false,
        global_artifacts_materialized = false,
        mixed_dispatcher_materialized = false,
        route_driver_wiring = false,
        coulomb_materialized = false,
        density_density_materialized = false,
        ida_mwg_data_materialized = false,
        pqs_lowdin_materialized = false,
        full_white_lindsey_route_assembled = false,
    )
end

function _one_body_pair_block_set_preflight_summary(plan; kwargs...)
    throw(
        ArgumentError(
            "one-body block-set preflight summary requires a PairBlockMaterializationPlan",
        ),
    )
end

function _one_body_pair_block_set_preflight_term_summary(plan_dispatch_summary)
    return (;
        object_kind =
            :cartesian_pair_block_mixed_one_body_block_set_preflight_term_summary,
        term = plan_dispatch_summary.requested_term,
        status = plan_dispatch_summary.status,
        blocker = plan_dispatch_summary.blocker,
        record_count = plan_dispatch_summary.record_count,
        ready_count = plan_dispatch_summary.ready_count,
        blocked_count = plan_dispatch_summary.blocked_count,
        selector_family_counts = plan_dispatch_summary.selector_family_counts,
        dispatch_status_counts = plan_dispatch_summary.dispatch_status_counts,
        blocker_counts = plan_dispatch_summary.blocker_counts,
        record_materialization_path_counts =
            plan_dispatch_summary.record_materialization_path_counts,
        materialization_path_counts =
            plan_dispatch_summary.materialization_path_counts,
        white_lindsey_unit_pair_required_count =
            plan_dispatch_summary.white_lindsey_unit_pair_required_count,
        white_lindsey_unit_pair_available_count =
            plan_dispatch_summary.white_lindsey_unit_pair_available_count,
        factors_constructed = false,
        numerical_blocks_materialized = false,
        materialized = false,
        source_operator_blocks_materialized = false,
        final_pair_blocks_materialized = false,
        operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
    )
end

function _one_body_block_set_parent_axis_counts_required(
    plan::PairBlockMaterializationPlan,
)
    return any(
        record -> _one_body_dispatch_parent_axis_counts_required(
            _one_body_record_selector_family(record),
        ),
        pair_block_materialization_records(plan),
    )
end

function _one_body_block_set_preflight_status(
    input_summary,
    total_ready_record_count::Int,
    total_blocked_record_count::Int,
)
    input_summary.status === :available_one_body_term_set_factor_inputs ||
        return :blocked_mixed_one_body_block_set_preflight
    total_blocked_record_count == 0 &&
        return :ready_metadata_only_not_materialized
    total_ready_record_count > 0 &&
        return :partially_ready_mixed_one_body_block_set_preflight
    return :blocked_mixed_one_body_block_set_preflight
end

function _one_body_block_set_preflight_blocker(
    input_summary,
    total_blocked_record_count::Int,
)
    input_summary.status === :available_one_body_term_set_factor_inputs ||
        return input_summary.blocker
    total_blocked_record_count == 0 && return nothing
    return :blocked_mixed_one_body_dispatch_records
end

function _one_body_pair_block(
    record::PairBlockMaterializationRecord,
    term::Symbol;
    inputs = (;),
    provider = nothing,
    unit_pair = nothing,
    materialize_selector_families = (
        :direct_direct,
        :pqs_source_pair,
        :white_lindsey_boundary_stratum,
    ),
)
    dispatch_summary = _one_body_pair_block_dispatch_summary(
        record,
        term;
        inputs,
        provider,
        unit_pair,
    )
    dispatch_summary.status === :ready_metadata_only_not_materialized ||
        return _one_body_skipped_pair_block_summary(dispatch_summary)

    dispatch_summary.selector_family in materialize_selector_families ||
        return _one_body_skipped_pair_block_summary(
            dispatch_summary,
            :mixed_one_body_selector_family_not_materialized_yet,
        )

    result_or_skip =
        dispatch_summary.selector_family === :direct_direct ?
        _one_body_direct_direct_pair_block(record, term; inputs, provider) :
        dispatch_summary.selector_family === :pqs_source_pair ?
        _one_body_pqs_source_pair_block(record, term; inputs, provider, dispatch_summary) :
        dispatch_summary.selector_family === :white_lindsey_boundary_stratum ?
        _one_body_white_lindsey_pair_block(
            record,
            term;
            inputs,
            provider,
            unit_pair,
            dispatch_summary,
        ) :
        return _one_body_skipped_pair_block_summary(
            dispatch_summary,
            :mixed_one_body_selector_family_not_materialized_yet,
        )
    result_or_skip isa PairBlockMaterializationResult || return result_or_skip
    return _one_body_result_with_mixed_metadata(result_or_skip, dispatch_summary)
end

function _one_body_pair_block(record, term; kwargs...)
    throw(
        ArgumentError(
            "mixed one-body pair block requires a PairBlockMaterializationRecord",
        ),
    )
end

function _one_body_pair_blocks(
    plan::PairBlockMaterializationPlan,
    term::Symbol;
    inputs = (;),
    provider = nothing,
    materialize_selector_families = (
        :direct_direct,
        :pqs_source_pair,
        :white_lindsey_boundary_stratum,
    ),
)
    unit_pair_lookup = _one_body_unit_pair_lookup(plan)
    materialized_results = PairBlockMaterializationResult[]
    skipped_records = NamedTuple[]
    for record in pair_block_materialization_records(plan)
        result_or_skip = _one_body_pair_block(
            record,
            term;
            inputs,
            provider,
            unit_pair = _one_body_dispatch_unit_pair(record, unit_pair_lookup),
            materialize_selector_families,
        )
        if result_or_skip isa PairBlockMaterializationResult
            push!(materialized_results, result_or_skip)
        else
            push!(skipped_records, NamedTuple(result_or_skip))
        end
    end

    result_tuple = Tuple(materialized_results)
    skipped_tuple = Tuple(skipped_records)
    any_materialized = !isempty(result_tuple)
    any_source_materialized =
        any(result -> result.source_operator_blocks_materialized, result_tuple)
    any_final_materialized =
        any(result -> result.final_pair_blocks_materialized, result_tuple)
    return PairBlockMaterializationBatchResult(
        term,
        result_tuple,
        skipped_tuple,
        length(result_tuple),
        length(skipped_tuple),
        any_materialized,
        any_source_materialized,
        any_final_materialized,
        false,
        false,
        false,
        (;
            materialization_path = :mixed_one_body_pair_block_batch_selector,
            mixed_one_body_dispatcher =
                _one_body_plan_dispatcher_name(materialize_selector_families),
            pair_block_record_count =
                length(pair_block_materialization_records(plan)),
            materialized_selector_family_counts =
                _one_body_materialized_selector_family_counts(result_tuple),
            skipped_selector_family_counts =
                _one_body_count_by(skipped_tuple, :selector_family),
            skipped_blocker_counts =
                _one_body_count_by(skipped_tuple, :blocker),
            skipped_status_counts =
                _one_body_count_by(skipped_tuple, :status),
            factors_constructed = false,
            numerical_dispatch_scope =
                _one_body_numerical_dispatch_scope(materialize_selector_families),
            pqs_source_pair_materialized =
                _one_body_result_selector_family_materialized(
                    result_tuple,
                    :pqs_source_pair,
                ),
            white_lindsey_materialized =
                _one_body_result_selector_family_materialized(
                    result_tuple,
                    :white_lindsey_boundary_stratum,
                ),
            route_driver_wiring = false,
            hamiltonian_data_materialized = false,
            artifacts_materialized = false,
        ),
    )
end

function _one_body_plan_dispatcher_name(materialize_selector_families)
    selector_families = Tuple(materialize_selector_families)
    selector_families == (:direct_direct,) &&
        return :direct_direct_only_plan_dispatcher
    selector_families == (:direct_direct, :pqs_source_pair) &&
        return :direct_pqs_source_plan_dispatcher
    selector_families == (
        :direct_direct,
        :pqs_source_pair,
        :white_lindsey_boundary_stratum,
    ) && return :direct_pqs_source_lw_plan_dispatcher
    return :mixed_one_body_plan_dispatcher
end

function _one_body_numerical_dispatch_scope(materialize_selector_families)
    selector_families = Tuple(materialize_selector_families)
    selector_families == (:direct_direct,) && return :direct_direct_only
    selector_families == (:direct_direct, :pqs_source_pair) &&
        return :direct_direct_and_pqs_source_pair
    selector_families == (
        :direct_direct,
        :pqs_source_pair,
        :white_lindsey_boundary_stratum,
    ) && return :direct_direct_pqs_source_pair_and_white_lindsey_boundary_stratum
    return :custom_selector_family_scope
end

function _one_body_direct_direct_pair_block(
    record::PairBlockMaterializationRecord,
    term::Symbol;
    inputs,
    provider,
)
    factor_input_summary = _one_body_factor_input_summary(
        term;
        inputs,
        provider,
        selector_family = :direct_direct,
        parent_axis_counts_required = true,
    )
    return direct_direct_one_body_block(
        record,
        term;
        parent_axis_counts = factor_input_summary.parent_axis_counts,
        overlap_1d = getproperty(factor_input_summary.factor_values, :overlap_1d),
        position_1d =
            _one_body_optional_factor(factor_input_summary.factor_values, :position_1d),
        x2_1d =
            _one_body_optional_factor(factor_input_summary.factor_values, :x2_1d),
        kinetic_1d =
            _one_body_optional_factor(factor_input_summary.factor_values, :kinetic_1d),
    )
end

function _one_body_pqs_source_pair_block(
    record::PairBlockMaterializationRecord,
    term::Symbol;
    inputs,
    provider,
    dispatch_summary,
)
    _one_body_pqs_source_metadata_ready(record) ||
        return _one_body_skipped_pair_block_summary(
            dispatch_summary,
            :missing_pqs_source_mode_metadata,
        )
    factor_input_summary = _one_body_factor_input_summary(
        term;
        inputs,
        provider,
        selector_family = :pqs_source_pair,
        parent_axis_counts_required = false,
    )
    return pqs_source_pair_one_body_block(
        record,
        term;
        overlap_1d = getproperty(factor_input_summary.factor_values, :overlap_1d),
        position_1d =
            _one_body_optional_factor(factor_input_summary.factor_values, :position_1d),
        x2_1d =
            _one_body_optional_factor(factor_input_summary.factor_values, :x2_1d),
        kinetic_1d =
            _one_body_optional_factor(factor_input_summary.factor_values, :kinetic_1d),
    )
end

function _one_body_white_lindsey_pair_block(
    record::PairBlockMaterializationRecord,
    term::Symbol;
    inputs,
    provider,
    unit_pair,
    dispatch_summary,
)
    unit_pair isa CUP.UnitPairRecord ||
        return _one_body_skipped_pair_block_summary(
            dispatch_summary,
            :missing_white_lindsey_unit_pair,
        )
    factor_input_summary = _one_body_factor_input_summary(
        term;
        inputs,
        provider,
        selector_family = :white_lindsey_boundary_stratum,
        parent_axis_counts_required = true,
    )
    return white_lindsey_boundary_stratum_one_body_block(
        unit_pair,
        term;
        parent_axis_counts = factor_input_summary.parent_axis_counts,
        overlap_1d = getproperty(factor_input_summary.factor_values, :overlap_1d),
        position_1d =
            _one_body_optional_factor(factor_input_summary.factor_values, :position_1d),
        x2_1d =
            _one_body_optional_factor(factor_input_summary.factor_values, :x2_1d),
        kinetic_1d =
            _one_body_optional_factor(factor_input_summary.factor_values, :kinetic_1d),
    )
end

function _one_body_pair_blocks(plan, term; kwargs...)
    throw(
        ArgumentError(
            "mixed one-body pair blocks require a PairBlockMaterializationPlan",
        ),
    )
end

function _one_body_pair_block_batch_summary(
    batch_result::PairBlockMaterializationBatchResult,
)
    result_tuple = batch_result.materialized_results
    skipped_tuple = batch_result.skipped_records
    return (;
        object_kind = :cartesian_pair_block_mixed_one_body_batch_summary,
        status = _one_body_batch_summary_status(batch_result),
        term = batch_result.term,
        materialized_count = batch_result.materialized_count,
        skipped_count = batch_result.skipped_count,
        materialized_selector_family_counts =
            _one_body_materialized_selector_family_counts(result_tuple),
        skipped_selector_family_counts =
            _one_body_count_optional_by(skipped_tuple, :selector_family),
        skipped_blocker_counts =
            _one_body_count_optional_by(skipped_tuple, :blocker),
        source_space_only_result_count =
            _one_body_source_space_only_result_count(result_tuple),
        final_local_block_result_count =
            _one_body_final_local_block_result_count(result_tuple),
        direct_direct_materialized =
            _one_body_result_selector_family_materialized(
                result_tuple,
                :direct_direct,
            ),
        pqs_source_pair_materialized =
            _one_body_result_selector_family_materialized(
                result_tuple,
                :pqs_source_pair,
            ),
        white_lindsey_boundary_stratum_materialized =
            _one_body_result_selector_family_materialized(
                result_tuple,
                :white_lindsey_boundary_stratum,
            ),
        white_lindsey_materialized =
            _one_body_result_selector_family_materialized(
                result_tuple,
                :white_lindsey_boundary_stratum,
            ),
        materialized = batch_result.materialized,
        source_operator_blocks_materialized =
            batch_result.source_operator_blocks_materialized,
        final_pair_blocks_materialized = batch_result.final_pair_blocks_materialized,
        operator_blocks_materialized = batch_result.operator_blocks_materialized,
        hamiltonian_data_materialized =
            batch_result.hamiltonian_data_materialized,
        artifacts_materialized = batch_result.artifacts_materialized,
        global_operator_blocks_materialized =
            batch_result.operator_blocks_materialized,
        global_hamiltonian_data_materialized =
            batch_result.hamiltonian_data_materialized,
        global_artifacts_materialized = batch_result.artifacts_materialized,
        materialization_path =
            _one_body_get_metadata(batch_result, :materialization_path, nothing),
        mixed_one_body_dispatcher =
            _one_body_get_metadata(batch_result, :mixed_one_body_dispatcher, nothing),
        numerical_dispatch_scope =
            _one_body_get_metadata(batch_result, :numerical_dispatch_scope, nothing),
        pair_block_record_count =
            _one_body_get_metadata(batch_result, :pair_block_record_count, nothing),
        factors_constructed =
            _one_body_get_metadata(batch_result, :factors_constructed, false),
        route_driver_wiring =
            _one_body_get_metadata(batch_result, :route_driver_wiring, false),
    )
end

function _one_body_pair_block_batch_summary(batch_result)
    throw(
        ArgumentError(
            "one-body pair block batch summary requires a PairBlockMaterializationBatchResult",
        ),
    )
end

function _one_body_pair_block_set_summary(
    plan::PairBlockMaterializationPlan;
    terms = _ONE_BODY_TERMS,
    term_batch_results = (;),
)
    term_set_descriptor = _one_body_term_set_descriptor(terms)
    batch_result_lookup =
        _one_body_block_set_batch_result_lookup(term_batch_results)
    term_summaries = Tuple(
        _one_body_pair_block_set_term_summary(
            term,
            _one_body_block_set_batch_result_for_term(
                term,
                batch_result_lookup,
            ),
        ) for term in term_set_descriptor.terms
    )
    supplied_term_count =
        count(summary -> summary.batch_result_supplied, term_summaries)
    deferred_term_count = length(term_summaries) - supplied_term_count
    total_materialized_count =
        sum(summary.materialized_count for summary in term_summaries)
    total_skipped_count = sum(summary.skipped_count for summary in term_summaries)
    any_source_materialized =
        any(summary -> summary.source_operator_blocks_materialized, term_summaries)
    any_final_materialized =
        any(summary -> summary.final_pair_blocks_materialized, term_summaries)
    any_operator_materialized =
        any(summary -> summary.operator_blocks_materialized, term_summaries)
    any_hamiltonian_materialized =
        any(summary -> summary.hamiltonian_data_materialized, term_summaries)
    any_artifacts_materialized =
        any(summary -> summary.artifacts_materialized, term_summaries)

    return (;
        object_kind = :cartesian_pair_block_mixed_one_body_block_set_summary,
        status = _one_body_block_set_summary_status(
            supplied_term_count,
            deferred_term_count,
            total_materialized_count,
            total_skipped_count,
        ),
        plan_record_count = length(pair_block_materialization_records(plan)),
        plan_summary_status = _one_body_namedtuple_value(
            summary(plan),
            :status,
            nothing,
        ),
        term_set_descriptor,
        requested_terms = term_set_descriptor.terms,
        terms = term_set_descriptor.terms,
        term_count = term_set_descriptor.term_count,
        supplied_term_count,
        deferred_term_count,
        term_statuses = Tuple(
            (;
                term = term_summary.term,
                status = term_summary.status,
                batch_result_supplied = term_summary.batch_result_supplied,
            ) for term_summary in term_summaries
        ),
        term_status_counts = _one_body_count_by(term_summaries, :status),
        total_materialized_count,
        total_skipped_count,
        materialized_selector_family_counts =
            _one_body_block_set_aggregate_counts(
                term_summaries,
                :materialized_selector_family_counts,
                :selector_family,
            ),
        skipped_selector_family_counts = _one_body_block_set_aggregate_counts(
            term_summaries,
            :skipped_selector_family_counts,
            :selector_family,
        ),
        skipped_blocker_counts = _one_body_block_set_aggregate_counts(
            term_summaries,
            :skipped_blocker_counts,
            :blocker,
        ),
        term_summaries,
        result_terms_remain_separated = true,
        block_set_results_summed = false,
        term_batch_results_stored_in_summary = false,
        factor_provider_scope = term_set_descriptor.factor_provider_scope,
        factors_constructed = false,
        materialized = total_materialized_count > 0,
        source_operator_blocks_materialized = any_source_materialized,
        final_pair_blocks_materialized = any_final_materialized,
        operator_blocks_materialized = any_operator_materialized,
        hamiltonian_data_materialized = any_hamiltonian_materialized,
        artifacts_materialized = any_artifacts_materialized,
        global_operator_blocks_materialized = any_operator_materialized,
        global_hamiltonian_data_materialized = any_hamiltonian_materialized,
        global_artifacts_materialized = any_artifacts_materialized,
        mixed_dispatcher_materialized = false,
        route_driver_wiring = false,
        coulomb_materialized = false,
        density_density_materialized = false,
        ida_mwg_data_materialized = false,
        pqs_lowdin_materialized = false,
        full_white_lindsey_route_assembled = false,
    )
end

function _one_body_pair_block_set_summary(plan; kwargs...)
    throw(
        ArgumentError(
            "one-body pair block-set summary requires a PairBlockMaterializationPlan",
        ),
    )
end

function _one_body_pair_block_set_term_summary(
    term::Symbol,
    batch_result::Nothing,
)
    return (;
        object_kind =
            :cartesian_pair_block_mixed_one_body_block_set_term_summary,
        status = :deferred_metadata_only_mixed_one_body_pair_block_batch,
        blocker = nothing,
        term,
        batch_result_supplied = false,
        materialized_count = 0,
        skipped_count = 0,
        materialized_selector_family_counts = (),
        skipped_selector_family_counts = (),
        skipped_blocker_counts = (),
        materialized = false,
        source_operator_blocks_materialized = false,
        final_pair_blocks_materialized = false,
        operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
    )
end

function _one_body_pair_block_set_term_summary(
    term::Symbol,
    batch_result::PairBlockMaterializationBatchResult,
)
    batch_result.term === term || throw(
        ArgumentError(
            "one-body block-set batch result term $(batch_result.term) does not match requested term $(term)",
        ),
    )
    batch_summary = _one_body_pair_block_batch_summary(batch_result)
    return (;
        object_kind =
            :cartesian_pair_block_mixed_one_body_block_set_term_summary,
        status = batch_summary.status,
        blocker = nothing,
        term,
        batch_result_supplied = true,
        materialized_count = batch_summary.materialized_count,
        skipped_count = batch_summary.skipped_count,
        materialized_selector_family_counts =
            batch_summary.materialized_selector_family_counts,
        skipped_selector_family_counts =
            batch_summary.skipped_selector_family_counts,
        skipped_blocker_counts = batch_summary.skipped_blocker_counts,
        materialized = batch_summary.materialized,
        source_operator_blocks_materialized =
            batch_summary.source_operator_blocks_materialized,
        final_pair_blocks_materialized =
            batch_summary.final_pair_blocks_materialized,
        operator_blocks_materialized = batch_summary.operator_blocks_materialized,
        hamiltonian_data_materialized =
            batch_summary.hamiltonian_data_materialized,
        artifacts_materialized = batch_summary.artifacts_materialized,
    )
end

function _one_body_block_set_batch_result_lookup(term_batch_results::NamedTuple)
    return term_batch_results
end

function _one_body_block_set_batch_result_lookup(
    batch_result::PairBlockMaterializationBatchResult,
)
    return NamedTuple{(batch_result.term,)}((batch_result,))
end

function _one_body_block_set_batch_result_lookup(batch_results::Tuple)
    all(result -> result isa PairBlockMaterializationBatchResult, batch_results) ||
        throw(
            ArgumentError(
                "one-body block-set term_batch_results tuple must contain PairBlockMaterializationBatchResult values",
            ),
        )
    terms = Tuple(result.term for result in batch_results)
    length(unique(terms)) == length(terms) || throw(
        ArgumentError("one-body block-set term_batch_results contain duplicate terms"),
    )
    return NamedTuple{terms}(batch_results)
end

function _one_body_block_set_batch_result_lookup(term_batch_results)
    throw(
        ArgumentError(
            "one-body block-set term_batch_results must be a NamedTuple, PairBlockMaterializationBatchResult, or tuple of PairBlockMaterializationBatchResult values",
        ),
    )
end

function _one_body_block_set_batch_result_for_term(
    term::Symbol,
    batch_result_lookup::NamedTuple,
)
    hasproperty(batch_result_lookup, term) || return nothing
    batch_result = getproperty(batch_result_lookup, term)
    batch_result isa PairBlockMaterializationBatchResult || throw(
        ArgumentError(
            "one-body block-set result for term $(term) must be a PairBlockMaterializationBatchResult",
        ),
    )
    return batch_result
end

function _one_body_block_set_summary_status(
    supplied_term_count::Int,
    deferred_term_count::Int,
    total_materialized_count::Int,
    total_skipped_count::Int,
)
    supplied_term_count == 0 &&
        return :deferred_metadata_only_mixed_one_body_block_set
    deferred_term_count > 0 &&
        return :partially_deferred_mixed_one_body_block_set
    total_materialized_count == 0 &&
        return :skipped_mixed_one_body_block_set
    total_skipped_count == 0 &&
        return :materialized_mixed_one_body_block_set
    return :partially_materialized_mixed_one_body_block_set
end

function _one_body_block_set_aggregate_counts(
    term_summaries::Tuple,
    count_field::Symbol,
    value_field::Symbol,
)
    counts = Dict{Any,Int}()
    ordered_values = Any[]
    for term_summary in term_summaries
        for entry in getproperty(term_summary, count_field)
            value = getproperty(entry, value_field)
            if !haskey(counts, value)
                counts[value] = 0
                push!(ordered_values, value)
            end
            counts[value] += entry.count
        end
    end
    return Tuple(
        NamedTuple{(value_field, :count)}((value, counts[value]))
        for value in ordered_values
    )
end

function _one_body_namedtuple_value(nt::NamedTuple, key::Symbol, default)
    haskey(nt, key) || return default
    return getproperty(nt, key)
end

function _one_body_pair_block_consumption(
    plan::PairBlockMaterializationPlan,
    term::Symbol;
    inputs = (;),
    provider = nothing,
    materialize_selector_families = (
        :direct_direct,
        :pqs_source_pair,
        :white_lindsey_boundary_stratum,
    ),
)
    batch_result = _one_body_pair_blocks(
        plan,
        term;
        inputs,
        provider,
        materialize_selector_families,
    )
    summary = _one_body_pair_block_batch_summary(batch_result)
    return (;
        object_kind = :cartesian_pair_block_mixed_one_body_consumption,
        status = summary.status,
        term = summary.term,
        batch_result,
        summary,
        materialized_count = summary.materialized_count,
        skipped_count = summary.skipped_count,
        direct_direct_materialized = summary.direct_direct_materialized,
        pqs_source_pair_materialized = summary.pqs_source_pair_materialized,
        white_lindsey_materialized = summary.white_lindsey_materialized,
        source_space_only_result_count = summary.source_space_only_result_count,
        final_local_block_result_count = summary.final_local_block_result_count,
        materialized = summary.materialized,
        source_operator_blocks_materialized =
            summary.source_operator_blocks_materialized,
        final_pair_blocks_materialized = summary.final_pair_blocks_materialized,
        operator_blocks_materialized = summary.global_operator_blocks_materialized,
        hamiltonian_data_materialized =
            summary.global_hamiltonian_data_materialized,
        artifacts_materialized = summary.global_artifacts_materialized,
        route_driver_wiring = summary.route_driver_wiring,
        factors_constructed = summary.factors_constructed,
        materialization_path = summary.materialization_path,
        mixed_one_body_dispatcher = summary.mixed_one_body_dispatcher,
        numerical_dispatch_scope = summary.numerical_dispatch_scope,
    )
end

function _one_body_pair_block_consumption(plan, term; kwargs...)
    throw(
        ArgumentError(
            "one-body pair block consumption requires a PairBlockMaterializationPlan",
        ),
    )
end

function _one_body_batch_summary_status(
    batch_result::PairBlockMaterializationBatchResult,
)
    batch_result.materialized_count == 0 &&
        return :skipped_mixed_one_body_pair_block_batch
    batch_result.skipped_count == 0 &&
        return :materialized_mixed_one_body_pair_block_batch
    return :partially_materialized_mixed_one_body_pair_block_batch
end

function _one_body_source_space_only_result_count(
    results::Tuple{Vararg{PairBlockMaterializationResult}},
)
    return count(
        result -> result.source_operator_blocks_materialized &&
                  !result.final_pair_blocks_materialized,
        results,
    )
end

function _one_body_final_local_block_result_count(
    results::Tuple{Vararg{PairBlockMaterializationResult}},
)
    return count(result -> result.final_pair_blocks_materialized, results)
end

function _one_body_record_selector_family(record::PairBlockMaterializationRecord)
    record.materialization_path === :direct_direct_pair_block_materialization_pilot &&
        return :direct_direct
    record.materialization_path === :pqs_source_pair_preflight &&
        return :pqs_source_pair
    record.materialization_path ===
        :white_lindsey_boundary_stratum_adapter_preflight &&
        return :white_lindsey_boundary_stratum
    return :unsupported
end

function _one_body_unit_pair_lookup(plan::PairBlockMaterializationPlan)
    lookup = Dict{Tuple{Symbol,Symbol},Any}()
    for unit_pair in CUP.unit_pairs(plan.pair_operator_plan.unit_pair_plan)
        lookup[unit_pair.pair_key] = unit_pair
    end
    return lookup
end

function _one_body_dispatch_unit_pair(
    record::PairBlockMaterializationRecord,
    unit_pair_lookup,
)
    record.materialization_path ===
        :white_lindsey_boundary_stratum_adapter_preflight || return nothing
    return get(unit_pair_lookup, record.pair_key, nothing)
end

function _one_body_optional_factor(factor_values::NamedTuple, name::Symbol)
    hasproperty(factor_values, name) || return nothing
    return getproperty(factor_values, name)
end

function _one_body_pqs_source_metadata_ready(record::PairBlockMaterializationRecord)
    required_keys = (
        :left_source_mode_dims,
        :right_source_mode_dims,
        :left_source_mode_count,
        :right_source_mode_count,
        :source_mode_ordering,
    )
    return all(
        key -> haskey(record.metadata, key) && !isnothing(getfield(record.metadata, key)),
        required_keys,
    )
end

function _one_body_result_with_mixed_metadata(
    result::PairBlockMaterializationResult,
    dispatch_summary,
)
    return PairBlockMaterializationResult(
        result.term,
        result.pair_key,
        result.block,
        result.materialized,
        result.source_operator_blocks_materialized,
        result.final_pair_blocks_materialized,
        result.operator_blocks_materialized,
        result.hamiltonian_data_materialized,
        result.artifacts_materialized,
        merge(
            result.metadata,
            (;
                mixed_one_body_dispatcher =
                    _one_body_record_dispatcher_name(dispatch_summary.selector_family),
                mixed_one_body_dispatch_status = dispatch_summary.status,
                selector_family = dispatch_summary.selector_family,
                factor_input_status = dispatch_summary.factor_input_status,
                numerical_family_selector =
                    _one_body_numerical_family_selector(dispatch_summary.selector_family),
            ),
        ),
    )
end

function _one_body_record_dispatcher_name(selector_family)
    selector_family === :direct_direct &&
        return :direct_direct_only_record_dispatcher
    selector_family === :pqs_source_pair &&
        return :pqs_source_pair_record_dispatcher
    selector_family === :white_lindsey_boundary_stratum &&
        return :white_lindsey_boundary_stratum_record_dispatcher
    return :mixed_one_body_record_dispatcher
end

function _one_body_numerical_family_selector(selector_family)
    selector_family === :direct_direct && return :direct_direct_one_body_block
    selector_family === :pqs_source_pair && return :pqs_source_pair_one_body_block
    selector_family === :white_lindsey_boundary_stratum &&
        return :white_lindsey_boundary_stratum_one_body_block
    return nothing
end

function _one_body_skipped_pair_block_summary(
    dispatch_summary,
    blocker = dispatch_summary.blocker,
)
    return (;
        object_kind = :cartesian_pair_block_mixed_one_body_skipped_summary,
        status = :skipped_mixed_one_body_pair_block,
        blocker,
        dispatch_status = dispatch_summary.status,
        dispatch_blocker = dispatch_summary.blocker,
        dispatch_blockers = dispatch_summary.blockers,
        pair_key = dispatch_summary.pair_key,
        pair_index = dispatch_summary.pair_index,
        pair_family = dispatch_summary.pair_family,
        selector_family = dispatch_summary.selector_family,
        requested_term = dispatch_summary.requested_term,
        materialization_path = dispatch_summary.materialization_path,
        record_materialization_path = dispatch_summary.record_materialization_path,
        record_readiness_status = dispatch_summary.record_readiness_status,
        record_blocker = dispatch_summary.record_blocker,
        factor_input_status = dispatch_summary.factor_input_status,
        factor_input_blocker = dispatch_summary.factor_input_blocker,
        factor_input_blockers = dispatch_summary.factor_input_blockers,
        unit_pair_requirement_status =
            dispatch_summary.unit_pair_requirement_status,
        materialized = false,
        factors_constructed = false,
        numerical_blocks_materialized = false,
        source_operator_blocks_materialized = false,
        final_pair_blocks_materialized = false,
        operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
        mixed_dispatcher_materialized = false,
        route_driver_wiring = false,
        pqs_lowdin_materialized = false,
        full_white_lindsey_route_assembled = false,
    )
end

function _one_body_count_by(record_summaries::Tuple, field::Symbol)
    values = Any[]
    for summary in record_summaries
        value = getproperty(summary, field)
        isnothing(value) && continue
        push!(values, value)
    end
    return _one_body_count_values(values, field)
end

function _one_body_materialized_selector_family_counts(
    results::Tuple{Vararg{PairBlockMaterializationResult}},
)
    selector_families = Any[]
    for result in results
        haskey(result.metadata, :selector_family) ||
            continue
        push!(selector_families, result.metadata.selector_family)
    end
    return _one_body_count_values(selector_families, :selector_family)
end

function _one_body_result_selector_family_materialized(
    results::Tuple{Vararg{PairBlockMaterializationResult}},
    selector_family::Symbol,
)
    return any(
        result ->
            haskey(result.metadata, :selector_family) &&
            result.metadata.selector_family === selector_family,
        results,
    )
end

function _one_body_blocker_counts(record_summaries::Tuple)
    blockers = Any[]
    for summary in record_summaries
        append!(blockers, summary.blockers)
    end
    return _one_body_count_values(blockers, :blocker)
end

function _one_body_count_optional_by(record_summaries::Tuple, field::Symbol)
    values = Any[]
    for summary in record_summaries
        haskey(summary, field) || continue
        value = getproperty(summary, field)
        isnothing(value) && continue
        push!(values, value)
    end
    return _one_body_count_values(values, field)
end

function _one_body_count_values(values, field::Symbol)
    counts = Dict{Any,Int}()
    ordered_values = Any[]
    for value in values
        if !haskey(counts, value)
            counts[value] = 0
            push!(ordered_values, value)
        end
        counts[value] += 1
    end
    return Tuple(
        NamedTuple{(field, :count)}((value, counts[value]))
        for value in ordered_values
    )
end

function _one_body_get_metadata(
    batch_result::PairBlockMaterializationBatchResult,
    key::Symbol,
    default,
)
    haskey(batch_result.metadata, key) || return default
    return getproperty(batch_result.metadata, key)
end

function _one_body_dispatch_parent_axis_counts_required(selector_family)
    return selector_family === :direct_direct ||
           selector_family === :white_lindsey_boundary_stratum
end

function _one_body_unit_pair_requirement(
    record::PairBlockMaterializationRecord,
    selector_family,
    unit_pair,
)
    required = selector_family === :white_lindsey_boundary_stratum
    if !required
        return (;
            required = false,
            available = !isnothing(unit_pair),
            pair_key = _one_body_unit_pair_key(unit_pair),
            status = :not_required,
            blocker = nothing,
        )
    end

    if isnothing(unit_pair)
        return (;
            required = true,
            available = false,
            pair_key = nothing,
            status = :missing_white_lindsey_unit_pair,
            blocker = :missing_white_lindsey_unit_pair,
        )
    end

    pair_key = _one_body_unit_pair_key(unit_pair)
    aligned = pair_key == record.pair_key
    return (;
        required = true,
        available = true,
        pair_key,
        status =
            aligned ?
            :available_aligned_white_lindsey_unit_pair :
            :mismatched_white_lindsey_unit_pair,
        blocker = aligned ? nothing : :mismatched_white_lindsey_unit_pair,
    )
end

function _one_body_unit_pair_key(unit_pair)
    isnothing(unit_pair) && return nothing
    hasproperty(unit_pair, :pair_key) && return getproperty(unit_pair, :pair_key)
    return nothing
end

function _one_body_dispatch_blockers(
    selector_family,
    record_ready::Bool,
    term_supported_by_record::Bool,
    factor_input_summary,
    unit_pair_requirement,
)
    blockers = Symbol[]
    selector_family === :unsupported &&
        push!(blockers, :unsupported_pair_block_materialization_path)
    record_ready || push!(blockers, :pair_block_record_not_ready)
    term_supported_by_record ||
        push!(blockers, :unsupported_one_body_term_for_record)
    factor_input_summary.status === :available_one_body_factor_inputs ||
        push!(blockers, :missing_one_body_factor_inputs)
    isnothing(unit_pair_requirement.blocker) ||
        push!(blockers, unit_pair_requirement.blocker)
    return Tuple(blockers)
end

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

function _one_body_pair_block(
    record::PairBlockMaterializationRecord,
    term::Symbol;
    inputs = (;),
    provider = nothing,
    unit_pair = nothing,
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

    dispatch_summary.selector_family === :direct_direct ||
        return _one_body_skipped_pair_block_summary(
            dispatch_summary,
            :mixed_one_body_selector_family_not_materialized_yet,
        )

    factor_input_summary = _one_body_factor_input_summary(
        term;
        inputs,
        provider,
        selector_family = :direct_direct,
        parent_axis_counts_required = true,
    )
    result = direct_direct_one_body_block(
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
    return _one_body_result_with_mixed_metadata(result, dispatch_summary)
end

function _one_body_pair_block(record, term; kwargs...)
    throw(
        ArgumentError(
            "mixed one-body pair block requires a PairBlockMaterializationRecord",
        ),
    )
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
                    :direct_direct_only_record_dispatcher,
                mixed_one_body_dispatch_status = dispatch_summary.status,
                selector_family = dispatch_summary.selector_family,
                factor_input_status = dispatch_summary.factor_input_status,
                numerical_family_selector = :direct_direct_one_body_block,
            ),
        ),
    )
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

function _one_body_blocker_counts(record_summaries::Tuple)
    blockers = Any[]
    for summary in record_summaries
        append!(blockers, summary.blockers)
    end
    return _one_body_count_values(blockers, :blocker)
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

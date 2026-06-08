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

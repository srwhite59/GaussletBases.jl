# Caller-supplied one-body factor input convention for future mixed dispatch.
#
# This layer only normalizes input metadata. It never constructs numerical
# factors or calls family materializers.

function _one_body_term_set_factor_input_summary(
    terms = _ONE_BODY_TERMS;
    inputs = (;),
    provider = nothing,
    parent_axis_counts_required::Bool = false,
)
    input_source = _one_body_factor_input_source(inputs, provider)
    term_set_descriptor = _one_body_term_set_descriptor(terms)
    _, present_factor_names, missing_factor_names =
        _one_body_required_factor_values(
            term_set_descriptor.required_factor_names,
            inputs,
            provider,
        )
    parent_axis_counts_available, parent_axis_counts =
        _one_body_input_value(:parent_axis_counts, inputs, provider)
    parent_axis_counts_status = _one_body_parent_axis_counts_status(
        parent_axis_counts_required,
        parent_axis_counts_available,
    )
    required_factors_available = isempty(missing_factor_names)
    parent_axis_counts_ready =
        !parent_axis_counts_required || parent_axis_counts_available
    blockers = _one_body_factor_input_blockers(
        required_factors_available,
        parent_axis_counts_ready,
    )
    ready_for_block_set = isempty(blockers)

    return (;
        object_kind =
            :cartesian_pair_block_one_body_term_set_factor_input_summary,
        status =
            ready_for_block_set ?
            :available_one_body_term_set_factor_inputs :
            :blocked_missing_one_body_term_set_inputs,
        blocker = ready_for_block_set ? nothing : first(blockers),
        blockers,
        term_set_descriptor,
        requested_terms = term_set_descriptor.terms,
        terms = term_set_descriptor.terms,
        term_count = term_set_descriptor.term_count,
        required_factor_names = term_set_descriptor.required_factor_names,
        present_factor_names,
        missing_factor_names,
        required_factors_available,
        input_source,
        parent_axis_counts_required,
        parent_axis_counts_status,
        parent_axis_counts,
        factor_provider_scope = term_set_descriptor.factor_provider_scope,
        factor_values_stored = false,
        factors_constructed = false,
        numerical_blocks_materialized = false,
        mixed_dispatcher_materialized = false,
        route_driver_wiring = false,
        global_operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
        coulomb_materialized = false,
        density_density_materialized = false,
        ida_mwg_data_materialized = false,
        pqs_lowdin_materialized = false,
        full_white_lindsey_route_assembled = false,
    )
end

function _one_body_factor_input_summary(
    term::Symbol;
    inputs = (;),
    provider = nothing,
    selector_family = nothing,
    materialization_path = nothing,
    parent_axis_counts_required::Bool = false,
)
    input_source = _one_body_factor_input_source(inputs, provider)
    descriptor = _one_body_term_descriptor(term)
    factor_values, present_factor_names, missing_factor_names =
        _one_body_required_factor_values(
            descriptor.required_factor_names,
            inputs,
            provider,
        )
    parent_axis_counts_available, parent_axis_counts =
        _one_body_input_value(:parent_axis_counts, inputs, provider)
    parent_axis_counts_status = _one_body_parent_axis_counts_status(
        parent_axis_counts_required,
        parent_axis_counts_available,
    )
    required_factors_available = isempty(missing_factor_names)
    parent_axis_counts_ready =
        !parent_axis_counts_required || parent_axis_counts_available
    blockers = _one_body_factor_input_blockers(
        required_factors_available,
        parent_axis_counts_ready,
    )
    ready_for_family_dispatch = isempty(blockers)
    resolved_materialization_path = _one_body_materialization_path(
        selector_family,
        materialization_path,
    )

    return (;
        object_kind = :cartesian_pair_block_one_body_factor_input_summary,
        status =
            ready_for_family_dispatch ?
            :available_one_body_factor_inputs :
            :blocked_missing_one_body_inputs,
        blocker = ready_for_family_dispatch ? nothing : first(blockers),
        blockers,
        requested_term = term,
        term_descriptor = descriptor,
        selector_family,
        selector_family_status =
            isnothing(selector_family) ? :not_specified : :specified,
        materialization_path = resolved_materialization_path,
        materialization_path_status =
            isnothing(resolved_materialization_path) ? :unknown : :available,
        input_source,
        parent_axis_counts_required,
        parent_axis_counts_status,
        parent_axis_counts,
        required_factor_names = descriptor.required_factor_names,
        present_factor_names,
        missing_factor_names,
        required_factors_available,
        factor_values,
        factors_constructed = false,
        numerical_blocks_materialized = false,
        mixed_dispatcher_materialized = false,
        route_driver_wiring = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
    )
end

function _one_body_factor_input_summary(term; kwargs...)
    throw(
        ArgumentError(
            "one-body factor input term must be a Symbol, got $(typeof(term))",
        ),
    )
end

function _one_body_factor_input_source(inputs, provider)
    inputs isa NamedTuple ||
        throw(ArgumentError("one-body factor inputs must be a NamedTuple"))
    if isnothing(provider)
        return :named_tuple
    end
    isempty(keys(inputs)) || throw(
        ArgumentError(
            "provide one-body factor inputs as either a NamedTuple or a provider callback",
        ),
    )
    provider isa Function ||
        throw(ArgumentError("one-body factor provider must be callable"))
    return :provider_callback
end

function _one_body_required_factor_values(
    required_factor_names::Tuple,
    inputs::NamedTuple,
    provider,
)
    present_names = Symbol[]
    missing_names = Symbol[]
    values = Any[]
    for name in required_factor_names
        available, value = _one_body_input_value(name, inputs, provider)
        if available
            push!(present_names, name)
            push!(values, value)
        else
            push!(missing_names, name)
        end
    end
    present_tuple = Tuple(present_names)
    return (
        NamedTuple{present_tuple}(Tuple(values)),
        present_tuple,
        Tuple(missing_names),
    )
end

function _one_body_input_value(name::Symbol, inputs::NamedTuple, provider)
    if isnothing(provider)
        hasproperty(inputs, name) || return false, nothing
        value = getproperty(inputs, name)
        return isnothing(value) ? (false, nothing) : (true, value)
    end
    value = provider(name)
    return isnothing(value) ? (false, nothing) : (true, value)
end

function _one_body_parent_axis_counts_status(required::Bool, available::Bool)
    if required
        return available ? :available_parent_axis_counts : :missing_parent_axis_counts
    end
    return available ? :available_parent_axis_counts_not_required : :not_required
end

function _one_body_factor_input_blockers(
    required_factors_available::Bool,
    parent_axis_counts_ready::Bool,
)
    blockers = Symbol[]
    required_factors_available ||
        push!(blockers, :missing_required_one_body_factors)
    parent_axis_counts_ready || push!(blockers, :missing_parent_axis_counts)
    return Tuple(blockers)
end

function _one_body_materialization_path(selector_family, materialization_path)
    !isnothing(materialization_path) && return materialization_path
    selector_family === :direct_direct &&
        return :direct_direct_one_body_selector
    selector_family === :pqs_source_pair &&
        return :pqs_source_pair_one_body_selector
    return nothing
end

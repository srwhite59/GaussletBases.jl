# Private route-shaped global one-body matrix adapter.
#
# This bridge composes the existing route-local one-body collection, local
# placement plan, and dense global safe one-body matrix pilots. It currently
# supports overlap, kinetic, and position_x/y/z only and does not assemble
# Hamiltonians, Coulomb data, IDA/MWG data, exports, artifacts, or PQS
# shell/Lowdin realizations.

function route_global_one_body_matrix(
    source;
    term::Symbol = :overlap,
    global_dimension = nothing,
    inputs = (;),
    provider = nothing,
    factors = nothing,
    factor_provider = nothing,
    metadata = (;),
)
    term in _route_global_one_body_supported_terms() ||
        return _route_global_one_body_blocked_result(
            term,
            :unsupported_route_global_one_body_term;
            metadata,
        )

    if isnothing(global_dimension)
        return _route_global_one_body_blocked_result(
            term,
            :missing_global_dimension;
            metadata,
        )
    end

    resolved_inputs = _route_global_one_body_inputs(inputs, factors)
    resolved_provider = _route_global_one_body_provider(provider, factor_provider)
    local_adapter = route_local_one_body_block_collection(
        source;
        terms = (term,),
        inputs = resolved_inputs,
        provider = resolved_provider,
        materialize_terms = (term,),
    )
    placement_plan = _route_global_one_body_placement_plan(
        term,
        local_adapter.local_block_collection;
        global_dimension,
    )
    global_matrix_result =
        _route_global_one_body_global_matrix(term, placement_plan)

    return _route_global_one_body_result(
        term,
        local_adapter,
        placement_plan,
        global_matrix_result;
        metadata,
    )
end

function route_global_overlap_matrix(source; kwargs...)
    return route_global_one_body_matrix(source; term = :overlap, kwargs...)
end

function route_global_kinetic_matrix(source; kwargs...)
    return route_global_one_body_matrix(source; term = :kinetic, kwargs...)
end

function route_global_position_matrix(source; axis, kwargs...)
    return route_global_one_body_matrix(
        source;
        term = _route_global_position_term(axis),
        kwargs...,
    )
end

function route_global_position_x_matrix(source; kwargs...)
    return route_global_one_body_matrix(source; term = :position_x, kwargs...)
end

function route_global_position_y_matrix(source; kwargs...)
    return route_global_one_body_matrix(source; term = :position_y, kwargs...)
end

function route_global_position_z_matrix(source; kwargs...)
    return route_global_one_body_matrix(source; term = :position_z, kwargs...)
end

function _route_global_one_body_supported_terms()
    return (:overlap, :kinetic, :position_x, :position_y, :position_z)
end

function _route_global_one_body_placement_plan(
    term::Symbol,
    collection::NamedTuple;
    global_dimension,
)
    term === :overlap && return one_body_overlap_placement_plan(
        collection;
        global_dimension,
    )
    term === :kinetic && return one_body_kinetic_placement_plan(
        collection;
        global_dimension,
    )
    term === :position_x && return one_body_position_x_placement_plan(
        collection;
        global_dimension,
    )
    term === :position_y && return one_body_position_y_placement_plan(
        collection;
        global_dimension,
    )
    term === :position_z && return one_body_position_z_placement_plan(
        collection;
        global_dimension,
    )
    throw(ArgumentError("unsupported route global one-body term: $(term)"))
end

function _route_global_one_body_global_matrix(
    term::Symbol,
    placement_plan::NamedTuple,
)
    term === :overlap && return one_body_global_overlap_matrix(placement_plan)
    term === :kinetic && return one_body_global_kinetic_matrix(placement_plan)
    term === :position_x &&
        return one_body_global_position_x_matrix(placement_plan)
    term === :position_y &&
        return one_body_global_position_y_matrix(placement_plan)
    term === :position_z &&
        return one_body_global_position_z_matrix(placement_plan)
    throw(ArgumentError("unsupported route global one-body term: $(term)"))
end

function _route_global_position_term(axis)
    axis === :x && return :position_x
    axis === :y && return :position_y
    axis === :z && return :position_z
    throw(ArgumentError("route global position axis must be :x, :y, or :z"))
end

function _route_global_one_body_inputs(inputs::NamedTuple, factors)
    isnothing(factors) && return inputs
    factors isa NamedTuple || throw(
        ArgumentError("route global one-body factors must be a NamedTuple"),
    )
    isempty(keys(inputs)) || throw(
        ArgumentError("provide route global one-body inputs or factors, not both"),
    )
    return factors
end

function _route_global_one_body_inputs(inputs, factors)
    throw(
        ArgumentError("route global one-body inputs must be a NamedTuple"),
    )
end

function _route_global_one_body_provider(provider, factor_provider)
    isnothing(factor_provider) && return provider
    isnothing(provider) || throw(
        ArgumentError(
            "provide route global one-body provider or factor_provider, not both",
        ),
    )
    factor_provider isa Function || throw(
        ArgumentError("route global one-body factor_provider must be callable"),
    )
    return factor_provider
end

function _route_global_one_body_result(
    term::Symbol,
    local_adapter::NamedTuple,
    placement_plan::NamedTuple,
    global_matrix_result::NamedTuple;
    metadata,
)
    materialized =
        _route_global_one_body_value(
            global_matrix_result,
            :operator_matrix_materialized,
            false,
        ) === true
    return (;
        object_kind = :cartesian_pair_block_route_global_one_body_matrix_adapter,
        term,
        status =
            materialized ?
            _route_global_one_body_materialized_status(term) :
            _route_global_one_body_blocked_status(term),
        blocker =
            materialized ?
            nothing :
            _route_global_one_body_value(global_matrix_result, :blocker, nothing),
        pair_block_materialization_plan =
            local_adapter.pair_block_materialization_plan,
        block_set_consumption = local_adapter.block_set_consumption,
        block_set_consumption_summary =
            local_adapter.block_set_consumption_summary,
        local_collection = local_adapter.local_block_collection,
        local_block_collection = local_adapter.local_block_collection,
        local_block_collection_summary =
            local_adapter.local_block_collection_summary,
        placement_plan,
        global_matrix_result,
        metadata = NamedTuple(metadata),
        local_adapter_status = local_adapter.status,
        local_adapter_blocker = local_adapter.blocker,
        placement_plan_status = placement_plan.status,
        placement_plan_blocker = placement_plan.blocker,
        global_matrix_status = global_matrix_result.status,
        global_matrix_blocker = global_matrix_result.blocker,
        materialized_local_block_count = local_adapter.materialized_entry_count,
        skipped_local_block_count = local_adapter.skipped_entry_count,
        placeable_record_count = placement_plan.placeable_count,
        blocked_placement_count = placement_plan.blocked_count,
        placed_block_count =
            _route_global_one_body_value(
                global_matrix_result,
                :placed_block_count,
                0,
            ),
        skipped_block_count =
            _route_global_one_body_value(
                global_matrix_result,
                :skipped_block_count,
                0,
            ),
        global_dimension =
            _route_global_one_body_value(
                placement_plan,
                :global_dimension,
                nothing,
            ),
        global_one_body_term_matrix_materialized = materialized,
        global_overlap_matrix_materialized =
            _route_global_one_body_value(
                global_matrix_result,
                :global_overlap_matrix_materialized,
                false,
            ),
        global_kinetic_matrix_materialized =
            _route_global_one_body_value(
                global_matrix_result,
                :global_kinetic_matrix_materialized,
                false,
            ),
        global_position_x_matrix_materialized =
            _route_global_one_body_value(
                global_matrix_result,
                :global_position_x_matrix_materialized,
                false,
            ),
        global_position_y_matrix_materialized =
            _route_global_one_body_value(
                global_matrix_result,
                :global_position_y_matrix_materialized,
                false,
            ),
        global_position_z_matrix_materialized =
            _route_global_one_body_value(
                global_matrix_result,
                :global_position_z_matrix_materialized,
                false,
            ),
        operator_matrix_materialized =
            _route_global_one_body_value(
                global_matrix_result,
                :operator_matrix_materialized,
                false,
            ),
        global_operator_assembled =
            _route_global_one_body_value(
                global_matrix_result,
                :global_operator_assembled,
                false,
            ),
        _route_global_one_body_nonclaim_flags()...,
    )
end

function _route_global_one_body_blocked_result(
    term::Symbol,
    blocker::Symbol;
    metadata,
)
    return (;
        object_kind = :cartesian_pair_block_route_global_one_body_matrix_adapter,
        term,
        status = _route_global_one_body_blocked_status(term),
        blocker,
        pair_block_materialization_plan = nothing,
        block_set_consumption = nothing,
        block_set_consumption_summary = nothing,
        local_collection = nothing,
        local_block_collection = nothing,
        local_block_collection_summary = nothing,
        placement_plan = nothing,
        global_matrix_result = nothing,
        metadata = NamedTuple(metadata),
        local_adapter_status = nothing,
        local_adapter_blocker = nothing,
        placement_plan_status = nothing,
        placement_plan_blocker = nothing,
        global_matrix_status = nothing,
        global_matrix_blocker = blocker,
        materialized_local_block_count = 0,
        skipped_local_block_count = 0,
        placeable_record_count = 0,
        blocked_placement_count = 0,
        placed_block_count = 0,
        skipped_block_count = 0,
        global_dimension = nothing,
        global_one_body_term_matrix_materialized = false,
        global_overlap_matrix_materialized = false,
        global_kinetic_matrix_materialized = false,
        global_position_x_matrix_materialized = false,
        global_position_y_matrix_materialized = false,
        global_position_z_matrix_materialized = false,
        operator_matrix_materialized = false,
        global_operator_assembled = false,
        _route_global_one_body_nonclaim_flags()...,
    )
end

function _route_global_one_body_materialized_status(term::Symbol)
    term in _route_global_one_body_supported_terms() &&
        return Symbol("materialized_route_global_", String(term), "_matrix")
    return :materialized_route_global_one_body_matrix
end

function _route_global_one_body_blocked_status(term::Symbol)
    term in _route_global_one_body_supported_terms() &&
        return Symbol("blocked_route_global_", String(term), "_matrix")
    return :blocked_route_global_one_body_matrix
end

function _route_global_one_body_nonclaim_flags()
    return (;
        route_driver_wiring = false,
        operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
        exports_materialized = false,
        global_operator_blocks_materialized = false,
        global_hamiltonian_data_materialized = false,
        global_artifacts_materialized = false,
        coulomb_materialized = false,
        density_density_materialized = false,
        ida_mwg_data_materialized = false,
        pqs_lowdin_materialized = false,
        pqs_shell_projection_materialized = false,
        full_white_lindsey_route_assembled = false,
    )
end

function _route_global_one_body_value(
    source::NamedTuple,
    key::Symbol,
    default = nothing,
)
    return haskey(source, key) ? getfield(source, key) : default
end

function _route_global_one_body_value(source, _key::Symbol, default = nothing)
    return default
end

# Private route-shaped global one-body matrix adapter.
#
# This bridge composes the existing route-local one-body collection, local
# placement plan, and dense global safe one-body matrix pilots. It currently
# supports the safe one-body terms only and does not assemble Hamiltonians,
# Coulomb data, IDA/MWG data, exports, artifacts, or PQS shell/Lowdin
# realizations.

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

function route_global_one_body_matrix_set(
    source;
    terms = route_global_safe_one_body_terms(),
    global_dimension = nothing,
    inputs = (;),
    provider = nothing,
    factors = nothing,
    factor_provider = nothing,
    metadata = (;),
)
    term_tuple = _route_global_one_body_term_tuple(terms)
    term_results = Tuple(
        route_global_one_body_matrix(
            source;
            term,
            global_dimension,
            inputs,
            provider,
            factors,
            factor_provider,
            metadata,
        ) for term in term_tuple
    )
    summary = _route_global_one_body_matrix_set_summary(
        term_tuple,
        term_results,
    )
    return (;
        object_kind =
            :cartesian_pair_block_route_global_safe_one_body_matrix_set,
        status = summary.status,
        blocker = summary.blocker,
        terms = term_tuple,
        term_results,
        summary,
        metadata = NamedTuple(metadata),
        all_terms_materialized = summary.all_terms_materialized,
        any_terms_materialized = summary.any_terms_materialized,
        term_count = summary.term_count,
        materialized_term_count = summary.materialized_term_count,
        blocked_term_count = summary.blocked_term_count,
        global_one_body_term_matrices_materialized =
            summary.all_terms_materialized,
        any_global_one_body_term_matrix_materialized =
            summary.any_terms_materialized,
        _route_global_one_body_nonclaim_flags()...,
    )
end

function route_global_safe_one_body_matrices(source; kwargs...)
    return route_global_one_body_matrix_set(source; kwargs...)
end

function route_global_one_body_matrix_set_result(
    matrix_set::NamedTuple,
    term::Symbol,
)
    for result in matrix_set.term_results
        result.term === term && return result
    end
    return nothing
end

function route_global_overlap_matrix(source; kwargs...)
    return route_global_one_body_matrix(source; term = :overlap, kwargs...)
end

function route_state_global_overlap_matrix(
    source;
    global_dimension = nothing,
    inputs = (;),
    provider = nothing,
    factors = nothing,
    factor_provider = nothing,
    metadata = (;),
)
    return _route_state_global_one_body_matrix(
        source,
        :overlap;
        global_dimension,
        inputs,
        provider,
        factors,
        factor_provider,
        metadata,
    )
end

function driver_global_overlap_result(
    source;
    global_dimension = nothing,
    inputs = (;),
    provider = nothing,
    factors = nothing,
    factor_provider = nothing,
    metadata = (;),
)
    return route_state_global_overlap_matrix(
        source;
        global_dimension,
        inputs,
        provider,
        factors,
        factor_provider,
        metadata,
    )
end

function route_global_kinetic_matrix(source; kwargs...)
    return route_global_one_body_matrix(source; term = :kinetic, kwargs...)
end

function route_state_global_kinetic_matrix(
    source;
    global_dimension = nothing,
    inputs = (;),
    provider = nothing,
    factors = nothing,
    factor_provider = nothing,
    metadata = (;),
)
    return _route_state_global_one_body_matrix(
        source,
        :kinetic;
        global_dimension,
        inputs,
        provider,
        factors,
        factor_provider,
        metadata,
    )
end

function route_global_position_matrix(source; axis, kwargs...)
    return route_global_one_body_matrix(
        source;
        term = _route_global_position_term(axis),
        kwargs...,
    )
end

function route_state_global_position_matrix(source; axis, kwargs...)
    return _route_state_global_one_body_matrix(
        source,
        _route_global_position_term(axis);
        kwargs...,
    )
end

function route_global_position_x_matrix(source; kwargs...)
    return route_global_one_body_matrix(source; term = :position_x, kwargs...)
end

function route_state_global_position_x_matrix(source; kwargs...)
    return _route_state_global_one_body_matrix(
        source,
        :position_x;
        kwargs...,
    )
end

function route_global_position_y_matrix(source; kwargs...)
    return route_global_one_body_matrix(source; term = :position_y, kwargs...)
end

function route_state_global_position_y_matrix(source; kwargs...)
    return _route_state_global_one_body_matrix(
        source,
        :position_y;
        kwargs...,
    )
end

function route_global_position_z_matrix(source; kwargs...)
    return route_global_one_body_matrix(source; term = :position_z, kwargs...)
end

function route_state_global_position_z_matrix(source; kwargs...)
    return _route_state_global_one_body_matrix(
        source,
        :position_z;
        kwargs...,
    )
end

function route_global_x2_matrix(source; axis, kwargs...)
    return route_global_one_body_matrix(
        source;
        term = _route_global_x2_term(axis),
        kwargs...,
    )
end

function route_state_global_x2_matrix(source; axis, kwargs...)
    return _route_state_global_one_body_matrix(
        source,
        _route_global_x2_term(axis);
        kwargs...,
    )
end

function route_global_x2_x_matrix(source; kwargs...)
    return route_global_one_body_matrix(source; term = :x2_x, kwargs...)
end

function route_state_global_x2_x_matrix(source; kwargs...)
    return _route_state_global_one_body_matrix(
        source,
        :x2_x;
        kwargs...,
    )
end

function route_global_x2_y_matrix(source; kwargs...)
    return route_global_one_body_matrix(source; term = :x2_y, kwargs...)
end

function route_state_global_x2_y_matrix(source; kwargs...)
    return _route_state_global_one_body_matrix(
        source,
        :x2_y;
        kwargs...,
    )
end

function route_global_x2_z_matrix(source; kwargs...)
    return route_global_one_body_matrix(source; term = :x2_z, kwargs...)
end

function route_state_global_x2_z_matrix(source; kwargs...)
    return _route_state_global_one_body_matrix(
        source,
        :x2_z;
        kwargs...,
    )
end

function route_global_electron_nuclear_by_center_matrix(
    source;
    parent_axis_counts = nothing,
    parent_axis_bundle_object = nothing,
    coulomb_expansion = nothing,
    center_record = nothing,
    global_dimension = nothing,
    inputs = (;),
    provider = nothing,
    factors = nothing,
    factor_provider = nothing,
    metadata = (;),
)
    if _route_global_one_body_pair_block_source(source)
        resolved_inputs = _route_global_electron_nuclear_by_center_inputs(
            _route_global_one_body_inputs(inputs, factors),
            parent_axis_bundle_object,
            coulomb_expansion,
            center_record,
        )
        return route_global_one_body_matrix(
            source;
            term = :electron_nuclear_by_center,
            global_dimension,
            inputs = resolved_inputs,
            provider,
            factor_provider,
            metadata,
        )
    end

    return _route_global_electron_nuclear_by_center_matrix_from_inventory(
        source;
        parent_axis_counts,
        parent_axis_bundle_object,
        coulomb_expansion,
        center_record,
        metadata,
    )
end

function route_global_electron_nuclear_by_center_matrices(
    source;
    parent_axis_counts = nothing,
    parent_axis_bundle_object = nothing,
    coulomb_expansion = nothing,
    center_records = (),
    metadata = (;),
)
    return @timeg "decomposed_wl.electron_nuclear_by_center.total" begin
        inventory = _route_global_electron_nuclear_by_center_inventory(
            source;
            parent_axis_counts,
            parent_axis_bundle_object,
        )
        if inventory.status !== :available_white_lindsey_decomposed_unit_pair_inventory
            return _route_global_electron_nuclear_by_center_matrices_result(
                inventory,
                (),
                :blocked_route_global_electron_nuclear_by_center_matrix_set,
                inventory.blocker;
                metadata,
            )
        end

        center_tuple = _route_global_electron_nuclear_by_center_records(center_records)
        isempty(center_tuple) &&
            return _route_global_electron_nuclear_by_center_matrices_result(
                inventory,
                (),
                :blocked_route_global_electron_nuclear_by_center_matrix_set,
                :missing_electron_nuclear_center_records;
                metadata,
            )

        matrix_results = Tuple(
            _route_global_electron_nuclear_by_center_matrix_from_inventory(
                inventory;
                parent_axis_counts,
                parent_axis_bundle_object,
                coulomb_expansion,
                center_record,
                metadata,
            ) for center_record in center_tuple
        )
        all_materialized = all(
            result -> result.global_electron_nuclear_by_center_matrix_materialized,
            matrix_results,
        )
        return _route_global_electron_nuclear_by_center_matrices_result(
            inventory,
            matrix_results,
            all_materialized ?
            :materialized_route_global_electron_nuclear_by_center_matrix_set :
            :blocked_route_global_electron_nuclear_by_center_matrix_set,
            all_materialized ?
            nothing :
            _route_global_electron_nuclear_by_center_first_blocker(matrix_results);
            metadata,
        )
    end
end

function route_state_global_electron_nuclear_by_center_matrix(source; kwargs...)
    return route_global_electron_nuclear_by_center_matrix(source; kwargs...)
end

function route_global_decomposed_wl_one_body_matrix(
    source;
    term::Symbol,
    parent_axis_counts = nothing,
    parent_axis_bundle_object = nothing,
    overlap_1d = nothing,
    position_1d = nothing,
    x2_1d = nothing,
    kinetic_1d = nothing,
    metadata = (;),
)
    timing_label =
        term === :overlap ? "decomposed_wl.overlap.total" :
        term === :kinetic ? "decomposed_wl.kinetic.total" :
        term in (:position_x, :position_y, :position_z) ?
        "decomposed_wl.position.total" :
        term in (:x2_x, :x2_y, :x2_z) ? "decomposed_wl.x2.total" :
        "decomposed_wl.one_body.total"
    return @timeg timing_label begin
    term in (
        :overlap,
        :kinetic,
        :position_x,
        :position_y,
        :position_z,
        :x2_x,
        :x2_y,
        :x2_z,
    ) ||
        return _route_global_decomposed_wl_one_body_blocked_result(
            term,
            :unsupported_decomposed_wl_one_body_term;
            metadata,
        )

    inventory = _route_global_electron_nuclear_by_center_inventory(
        source;
        parent_axis_counts,
        parent_axis_bundle_object,
    )
    if inventory.status !== :available_white_lindsey_decomposed_unit_pair_inventory
        return _route_global_decomposed_wl_one_body_blocked_result(
            term,
            inventory.blocker;
            inventory,
            metadata,
        )
    end

    if inventory.unit_pairs isa CUP.UnitPairIndexTable
        factorized = _route_global_decomposed_wl_factorized_one_body_matrix(
            inventory,
            term;
            parent_axis_counts,
            overlap_1d,
            position_1d,
            x2_1d,
            kinetic_1d,
            metadata,
        )
        if factorized.status === :materialized_decomposed_wl_factorized_one_body_matrix
            return _route_global_decomposed_wl_one_body_result(
                inventory,
                factorized.local_batch,
                factorized.local_collection,
                factorized.placement_plan,
                factorized.global_matrix_result;
                term,
                metadata,
            )
        end
        streaming = _route_global_decomposed_wl_streaming_one_body_matrix(
            inventory,
            term;
            parent_axis_counts,
            overlap_1d,
            position_1d,
            x2_1d,
            kinetic_1d,
            metadata,
        )
        return _route_global_decomposed_wl_one_body_result(
            inventory,
            streaming.local_batch,
            streaming.local_collection,
            streaming.placement_plan,
            streaming.global_matrix_result;
            term,
            metadata,
        )
    end

    local_batch = @timeg "decomposed_wl.one_body.local_batch" begin
        white_lindsey_boundary_stratum_one_body_blocks(
            inventory.unit_pairs,
            term;
            parent_axis_counts,
            overlap_1d,
            position_1d,
            x2_1d,
            kinetic_1d,
        )
    end
    local_collection = _route_global_decomposed_wl_one_body_local_collection(
        local_batch,
        inventory.unit_pairs,
        term,
    )
    placement_plan = nothing
    global_matrix_result = nothing
    @timeg "decomposed_wl.one_body.global_placement" begin
        placement_plan = _route_global_one_body_placement_plan(
            term,
            local_collection;
            global_dimension = inventory.retained_dimension,
        )
        global_matrix_result =
            _route_global_one_body_global_matrix(term, placement_plan)
    end
    return _route_global_decomposed_wl_one_body_result(
        inventory,
        local_batch,
        local_collection,
        placement_plan,
        global_matrix_result;
        term,
        metadata,
    )
    end
end

function route_global_decomposed_wl_overlap_matrix(source; kwargs...)
    return route_global_decomposed_wl_one_body_matrix(
        source;
        term = :overlap,
        kwargs...,
    )
end

function route_global_decomposed_wl_kinetic_matrix(source; kwargs...)
    return route_global_decomposed_wl_one_body_matrix(
        source;
        term = :kinetic,
        kwargs...,
    )
end

function route_global_decomposed_wl_position_x_matrix(source; kwargs...)
    return route_global_decomposed_wl_one_body_matrix(
        source;
        term = :position_x,
        kwargs...,
    )
end

function route_global_decomposed_wl_position_y_matrix(source; kwargs...)
    return route_global_decomposed_wl_one_body_matrix(
        source;
        term = :position_y,
        kwargs...,
    )
end

function route_global_decomposed_wl_position_z_matrix(source; kwargs...)
    return route_global_decomposed_wl_one_body_matrix(
        source;
        term = :position_z,
        kwargs...,
    )
end

function route_global_decomposed_wl_x2_x_matrix(source; kwargs...)
    return route_global_decomposed_wl_one_body_matrix(
        source;
        term = :x2_x,
        kwargs...,
    )
end

function route_global_decomposed_wl_x2_y_matrix(source; kwargs...)
    return route_global_decomposed_wl_one_body_matrix(
        source;
        term = :x2_y,
        kwargs...,
    )
end

function route_global_decomposed_wl_x2_z_matrix(source; kwargs...)
    return route_global_decomposed_wl_one_body_matrix(
        source;
        term = :x2_z,
        kwargs...,
    )
end

function white_lindsey_decomposed_one_electron_hamiltonian(
    kinetic_result::NamedTuple,
    electron_nuclear_by_center_results::NamedTuple;
    metadata = (;),
)
    kinetic_matrix = _route_global_one_body_value(kinetic_result, :matrix, nothing)
    kinetic_materialized =
        _route_global_one_body_value(
            kinetic_result,
            :global_kinetic_matrix_materialized,
            false,
        ) === true &&
        kinetic_matrix isa AbstractMatrix
    kinetic_materialized ||
        return _white_lindsey_decomposed_one_electron_hamiltonian_blocked_result(
            :missing_route_global_kinetic_matrix;
            metadata,
        )

    nuclear_materialized =
        _route_global_one_body_value(
            electron_nuclear_by_center_results,
            :global_electron_nuclear_by_center_matrices_materialized,
            false,
        ) === true
    nuclear_materialized ||
        return _white_lindsey_decomposed_one_electron_hamiltonian_blocked_result(
            :missing_route_global_electron_nuclear_by_center_matrices;
            metadata,
        )

    if _route_global_one_body_value(
        electron_nuclear_by_center_results,
        :nuclear_charge_applied,
        false,
    ) === true
        return _white_lindsey_decomposed_one_electron_hamiltonian_blocked_result(
            :nuclear_charge_already_applied_to_by_center_matrices;
            metadata,
        )
    end
    if _route_global_one_body_value(
        electron_nuclear_by_center_results,
        :centers_summed,
        false,
    ) === true
        return _white_lindsey_decomposed_one_electron_hamiltonian_blocked_result(
            :centers_already_summed_in_by_center_matrices;
            metadata,
        )
    end

    matrix_results = _route_global_one_body_value(
        electron_nuclear_by_center_results,
        :matrix_results,
        (),
    )
    dimension = size(kinetic_matrix)
    hamiltonian = Matrix{Float64}(kinetic_matrix)
    applied_center_count = 0
    center_indices = Any[]
    nuclear_charges = Float64[]
    for result in matrix_results
        by_center_matrix = _route_global_one_body_value(result, :matrix, nothing)
        by_center_matrix isa AbstractMatrix ||
            return _white_lindsey_decomposed_one_electron_hamiltonian_blocked_result(
                :missing_route_global_electron_nuclear_by_center_matrix;
                metadata,
            )
        size(by_center_matrix) == dimension ||
            return _white_lindsey_decomposed_one_electron_hamiltonian_blocked_result(
                :route_global_electron_nuclear_matrix_dimension_mismatch;
                metadata,
            )
        charge = _route_global_one_body_value(result, :nuclear_charge, nothing)
        charge isa Real ||
            return _white_lindsey_decomposed_one_electron_hamiltonian_blocked_result(
                :missing_nuclear_charge_for_hamiltonian_assembly;
                metadata,
            )

        # The WL by-center result is already a unit-charge nuclear-attraction
        # matrix; Hamiltonian assembly applies only the recorded nuclear charge.
        hamiltonian .+= Float64(charge) .* by_center_matrix
        applied_center_count += 1
        push!(
            center_indices,
            _route_global_one_body_value(result, :center_index, nothing),
        )
        push!(nuclear_charges, Float64(charge))
    end

    return (;
        object_kind =
            :white_lindsey_decomposed_one_electron_hamiltonian_assembly,
        status = :materialized_white_lindsey_decomposed_one_electron_hamiltonian,
        blocker = nothing,
        term = :one_electron_hamiltonian,
        matrix = hamiltonian,
        matrix_shape = size(hamiltonian),
        retained_dimension = size(hamiltonian, 1),
        kinetic_status = kinetic_result.status,
        kinetic_matrix_materialized = true,
        electron_nuclear_by_center_status =
            electron_nuclear_by_center_results.status,
        by_center_matrix_count =
            electron_nuclear_by_center_results.by_center_matrix_count,
        applied_center_count,
        center_indices = Tuple(center_indices),
        nuclear_charges = Tuple(nuclear_charges),
        unit_charge_nuclear_attraction_matrices_consumed = true,
        nuclear_charge_application_convention =
            :multiply_unit_charge_nuclear_attraction_by_recorded_charge,
        charge_application_stage = :white_lindsey_hamiltonian_assembly,
        nuclear_charge_applied_at_hamiltonian_assembly = true,
        center_summation_stage = :white_lindsey_hamiltonian_assembly,
        centers_summed_at_hamiltonian_assembly = true,
        by_center_matrices_remain_separated =
            !electron_nuclear_by_center_results.centers_summed,
        overlap_matrix_consumed = false,
        overlap_used_as_solve_metric_only = true,
        full_parent_window_cpb_used = false,
        direct_cartesian_product_assembly_used = false,
        ordinary_cartesian_ida_operators_used = false,
        metadata = NamedTuple(metadata),
        _route_global_one_body_nonclaim_flags(
            hamiltonian_data_materialized = true,
            global_hamiltonian_data_materialized = true,
        )...,
    )
end

function _route_global_decomposed_wl_one_body_local_collection(
    batch_result::PairBlockMaterializationBatchResult,
    unit_pairs,
    term::Symbol,
)
    pair_lookup = Dict(pair.pair_key => pair for pair in unit_pairs)
    materialized_entries = Tuple(
        merge(
            _one_body_local_block_collection_entry(result; block_set_term = term),
            (;
                left_column_range =
                    pair_lookup[result.pair_key].left_unit.column_range,
                right_column_range =
                    pair_lookup[result.pair_key].right_unit.column_range,
            ),
        ) for result in batch_result.materialized_results
    )
    skipped_entries = Tuple(
        _one_body_local_block_collection_skipped_entry(skip; block_set_term = term)
        for skip in batch_result.skipped_records
    )
    entries = (materialized_entries..., skipped_entries...)
    return (;
        object_kind = :cartesian_pair_block_local_one_body_block_collection,
        status =
            isempty(materialized_entries) ?
            :blocked_local_one_body_block_collection :
            :materialized_local_one_body_block_collection,
        blocker =
            isempty(materialized_entries) ?
            Symbol("no_materialized_", String(term), "_blocks") :
            nothing,
        block_set_status =
            isempty(materialized_entries) ?
            :blocked_mixed_one_body_block_set_consumption :
            :materialized_mixed_one_body_block_set_consumption,
        block_set_blocker =
            isempty(materialized_entries) ?
            Symbol("no_materialized_", String(term), "_blocks") :
            nothing,
        terms = (term,),
        requested_terms = (term,),
        requested_materialize_terms = (term,),
        materialized_terms = isempty(materialized_entries) ? () : (term,),
        deferred_terms = (),
        term_statuses = (),
        entries,
        materialized_entries,
        skipped_entries,
        entry_count = length(entries),
        materialized_entry_count = length(materialized_entries),
        skipped_entry_count = length(skipped_entries),
        deferred_term_count = 0,
        source_space_entry_count = 0,
        final_local_entry_count = length(materialized_entries),
        term_separated_entries = true,
        pair_separated_entries = true,
        block_set_results_summed = false,
        block_matrices_copied_into_collection = false,
        _one_body_local_block_collection_summary_materialization_flags(
            source_operator_blocks_materialized =
                batch_result.source_operator_blocks_materialized,
            final_pair_blocks_materialized =
                batch_result.final_pair_blocks_materialized,
        )...,
    )
end

function _route_global_decomposed_wl_streaming_one_body_matrix(
    inventory,
    term::Symbol;
    parent_axis_counts,
    overlap_1d = nothing,
    position_1d = nothing,
    x2_1d = nothing,
    kinetic_1d = nothing,
    metadata = (;),
)
    dimension = inventory.retained_dimension
    matrix = zeros(Float64, dimension, dimension)
    coefficient_cache =
        _white_lindsey_decomposed_operator_unit_coefficient_cache(inventory)
    prepared_cache =
        _white_lindsey_decomposed_operator_prepared_unit_cache(inventory)
    stream_state = @timeg "decomposed_wl.one_body.local_batch" begin
        _route_global_decomposed_wl_streaming_fill!(
            matrix,
            inventory.unit_pairs,
            term;
            parent_axis_counts,
            overlap_1d,
            position_1d,
            x2_1d,
            kinetic_1d,
            unit_coefficient_cache = coefficient_cache,
            prepared_unit_cache = prepared_cache,
        )
    end
    global_matrix_result =
        _route_global_decomposed_wl_streaming_global_matrix_result(
            term,
            inventory,
            matrix,
            stream_state,
        )
    local_batch =
        _route_global_decomposed_wl_streaming_local_batch(term, stream_state)
    local_collection =
        _route_global_decomposed_wl_streaming_local_collection(term, stream_state)
    placement_plan =
        _route_global_decomposed_wl_streaming_placement_plan(
            term,
            inventory,
            stream_state,
        )
    return (;
        local_batch,
        local_collection,
        placement_plan,
        global_matrix_result,
        metadata = NamedTuple(metadata),
    )
end

function _route_global_decomposed_wl_streaming_fill!(
    matrix::AbstractMatrix{Float64},
    unit_pairs::CUP.UnitPairIndexTable,
    term::Symbol;
    parent_axis_counts,
    overlap_1d = nothing,
    position_1d = nothing,
    x2_1d = nothing,
    kinetic_1d = nothing,
    parent_axis_bundle_object = nothing,
    coulomb_expansion = nothing,
    center_record = nothing,
    electron_nuclear_axis_context = nothing,
    unit_coefficient_cache = nothing,
    prepared_unit_cache = nothing,
    electron_nuclear_scratch = nothing,
)
    cache =
        isnothing(unit_coefficient_cache) ?
        Dict{Symbol,Any}() :
        unit_coefficient_cache
    prepared_cache =
        isnothing(prepared_unit_cache) ?
        Dict{Any,Any}() :
        prepared_unit_cache
    axis_counts =
        (term in (
             :overlap,
             :kinetic,
             :position_x,
             :position_y,
             :position_z,
             :x2_x,
             :x2_y,
             :x2_z,
         ) ||
         term === :electron_nuclear_by_center) ?
        _axis_counts_tuple(parent_axis_counts) :
        nothing
    overlap_axes =
        term in (
            :overlap,
            :kinetic,
            :position_x,
            :position_y,
            :position_z,
            :x2_x,
            :x2_y,
            :x2_z,
        ) ? _overlap_1d_tuple(overlap_1d) : nothing
    position_axes = term in (:position_x, :position_y, :position_z) ?
                    _operator_1d_tuple(position_1d, "position_1d") :
                    nothing
    x2_axes = term in (:x2_x, :x2_y, :x2_z) ?
              _operator_1d_tuple(x2_1d, "x2_1d") :
              nothing
    kinetic_axes = term === :kinetic ? _operator_1d_tuple(
        kinetic_1d,
        "kinetic_1d",
    ) : nothing
    if term in (
        :overlap,
        :kinetic,
        :position_x,
        :position_y,
        :position_z,
        :x2_x,
        :x2_y,
        :x2_z,
    )
        _assert_overlap_axis_sizes(overlap_axes, axis_counts)
    end
    if term in (:position_x, :position_y, :position_z)
        _assert_operator_axis_sizes(position_axes, axis_counts, "position_1d")
    end
    if term in (:x2_x, :x2_y, :x2_z)
        _assert_operator_axis_sizes(x2_axes, axis_counts, "x2_1d")
    end
    if term === :kinetic
        _assert_operator_axis_sizes(kinetic_axes, axis_counts, "kinetic_1d")
    end
    materialized_count = 0
    skipped_count = 0
    placed_block_count = 0
    blocker = nothing
    center_index = nothing
    center_key = nothing
    center_location = nothing
    nuclear_charge = nothing
    nuclear_charge_applied = false
    for unit_pair in unit_pairs
        if term === :electron_nuclear_by_center &&
           !isnothing(electron_nuclear_axis_context)
            axis_context = electron_nuclear_axis_context
            left_coefficients =
                _white_lindsey_unit_coefficients_from_local_cache(
                    cache,
                    unit_pair.left_unit,
                )
            right_coefficients =
                _white_lindsey_unit_coefficients_from_local_cache(
                    cache,
                    unit_pair.right_unit,
                )
            left_ready =
                _white_lindsey_unit_coefficients_materialized(left_coefficients)
            right_ready =
                _white_lindsey_unit_coefficients_materialized(right_coefficients)
            if !(left_ready && right_ready)
                skipped_count += 1
                if isnothing(blocker)
                    _status, pair_blocker =
                        _white_lindsey_pair_unit_coefficients_status(
                            left_ready,
                            right_ready,
                        )
                    blocker = pair_blocker
                end
                continue
            end
            left_unit =
                _white_lindsey_prepared_unit_for_local_operator_from_cache(
                    prepared_cache,
                    cache,
                    unit_pair.left_unit,
                    axis_context.axis_counts,
                )
            right_unit =
                _white_lindsey_prepared_unit_for_local_operator_from_cache(
                    prepared_cache,
                    cache,
                    unit_pair.right_unit,
                    axis_context.axis_counts,
                )
            block =
                isnothing(electron_nuclear_scratch) ?
                _white_lindsey_electron_nuclear_pair_block_from_prepared_units(
                    left_unit,
                    right_unit,
                    axis_context,
                ) :
                _white_lindsey_fill_electron_nuclear_pair_block!(
                    electron_nuclear_scratch,
                    left_unit,
                    right_unit,
                    axis_context,
                )
            rows = unit_pair.left_unit.column_range
            columns = unit_pair.right_unit.column_range
            placement_blocker =
                _route_global_decomposed_wl_streaming_placement_blocker(
                    block,
                    rows,
                    columns,
                    size(matrix, 1),
                )
            if !isnothing(placement_blocker)
                skipped_count += 1
                isnothing(blocker) && (blocker = placement_blocker)
                continue
            end
            placed_block_count +=
                _white_lindsey_accumulate_streaming_block!(
                    matrix,
                    rows,
                    columns,
                    block,
                )
            materialized_count += 1
            if isnothing(center_index)
                center_summary = axis_context.center_summary
                center_index = center_summary.center_index
                center_key = center_summary.center_key
                center_location = center_summary.location
                nuclear_charge = center_summary.charge
                nuclear_charge_applied = false
            end
            continue
        end
        left_coefficients =
            _white_lindsey_unit_coefficients_from_local_cache(
                cache,
                unit_pair.left_unit,
            )
        right_coefficients =
            _white_lindsey_unit_coefficients_from_local_cache(
                cache,
                unit_pair.right_unit,
            )
        left_ready =
            _white_lindsey_unit_coefficients_materialized(left_coefficients)
        right_ready =
            _white_lindsey_unit_coefficients_materialized(right_coefficients)
        if !(left_ready && right_ready)
            skipped_count += 1
            if isnothing(blocker)
                _status, pair_blocker =
                    _white_lindsey_pair_unit_coefficients_status(
                        left_ready,
                        right_ready,
                    )
                blocker = pair_blocker
            end
            continue
        end
        left_unit =
            _white_lindsey_prepared_unit_for_local_operator_from_cache(
                prepared_cache,
                cache,
                unit_pair.left_unit,
                axis_counts,
            )
        right_unit =
            _white_lindsey_prepared_unit_for_local_operator_from_cache(
                prepared_cache,
                cache,
                unit_pair.right_unit,
                axis_counts,
            )
        block =
            term === :overlap ?
            _white_lindsey_pair_product_block_from_prepared_units(
                left_unit,
                right_unit,
                overlap_axes,
            ) :
            (
                term === :kinetic ?
                _white_lindsey_pair_kinetic_block_from_prepared_units(
                    left_unit,
                    right_unit,
                    overlap_axes,
                    kinetic_axes,
                ) :
                (
                    term in (:position_x, :position_y, :position_z) ?
                    _white_lindsey_pair_product_block_from_prepared_units(
                        left_unit,
                        right_unit,
                        _route_global_decomposed_wl_axis_product_axes(
                            term,
                            overlap_axes,
                            position_axes,
                        ),
                    ) :
                    (
                        term in (:x2_x, :x2_y, :x2_z) ?
                        _white_lindsey_pair_product_block_from_prepared_units(
                            left_unit,
                            right_unit,
                            _route_global_decomposed_wl_axis_product_axes(
                                term,
                                overlap_axes,
                                x2_axes,
                            ),
                        ) :
                        nothing
                    )
                )
            )
        isnothing(block) && throw(
            ArgumentError("unsupported streaming decomposed WL one-body term $(term)"),
        )
        rows = unit_pair.left_unit.column_range
        columns = unit_pair.right_unit.column_range
        placement_blocker =
            _route_global_decomposed_wl_streaming_placement_blocker(
                block,
                rows,
                columns,
                size(matrix, 1),
            )
        if !isnothing(placement_blocker)
            skipped_count += 1
            isnothing(blocker) && (blocker = placement_blocker)
            continue
        end
        placed_block_count += _white_lindsey_accumulate_streaming_block!(
            matrix,
            rows,
            columns,
            block,
        )
        materialized_count += 1
    end
    return (;
        materialized_count,
        skipped_count,
        placed_block_count,
        blocker,
        unit_coefficient_cache_entry_count = length(cache),
        electron_nuclear_prepared_unit_cache_entry_count =
            length(prepared_cache),
        pair_count = length(unit_pairs),
        center_index,
        center_key,
        center_location,
        nuclear_charge,
        nuclear_charge_recorded = term === :electron_nuclear_by_center &&
                                  !isnothing(center_index),
        nuclear_charge_applied,
        materialization_path =
            :white_lindsey_decomposed_wl_streaming_one_body_global_assembly,
        local_blocks_collected = false,
        pair_lookup_materialized = false,
        placement_records_materialized = false,
    )
end

function _white_lindsey_accumulate_streaming_block!(
    matrix::AbstractMatrix{Float64},
    rows::UnitRange{Int},
    columns::UnitRange{Int},
    block,
)
    row_first = first(rows)
    column_first = first(columns)
    @inbounds for local_column in axes(block, 2)
        global_column = column_first + local_column - 1
        for local_row in axes(block, 1)
            global_row = row_first + local_row - 1
            matrix[global_row, global_column] += block[local_row, local_column]
        end
    end
    rows == columns && return 1
    @inbounds for local_column in axes(block, 2)
        global_row = column_first + local_column - 1
        for local_row in axes(block, 1)
            global_column = row_first + local_row - 1
            matrix[global_row, global_column] += block[local_row, local_column]
        end
    end
    return 2
end

function _white_lindsey_decomposed_operator_cache(inventory)
    cache = _route_global_one_body_value(inventory, :operator_cache, nothing)
    isnothing(cache) && return _white_lindsey_decomposed_operator_cache()
    return cache
end

function _white_lindsey_decomposed_operator_unit_coefficient_cache(inventory)
    cache = _white_lindsey_decomposed_operator_cache(inventory)
    return cache.unit_coefficients
end

function _white_lindsey_decomposed_operator_prepared_unit_cache(inventory)
    cache = _white_lindsey_decomposed_operator_cache(inventory)
    return cache.prepared_units
end

function _white_lindsey_decomposed_operator_factorized_basis_cache(inventory)
    cache = _white_lindsey_decomposed_operator_cache(inventory)
    hasproperty(cache, :factorized_retained_basis) ||
        return Dict{Any,Any}()
    return cache.factorized_retained_basis
end

function _white_lindsey_decomposed_factorized_retained_basis(
    inventory,
    axis_counts::NTuple{3,Int},
)
    cache = _white_lindsey_decomposed_operator_factorized_basis_cache(inventory)
    cache_key = (:factorized_retained_basis, axis_counts)
    haskey(cache, cache_key) && return cache[cache_key]
    result = try
        coefficient_cache =
            _white_lindsey_decomposed_operator_unit_coefficient_cache(inventory)
        coefficient_matrix =
            _white_lindsey_decomposed_retained_coefficient_matrix(
                inventory,
                axis_counts,
                coefficient_cache,
            )
        factorized_basis = @timeg "decomposed_wl.factorized_retained_basis.extract" begin
            ParentGaussletBases._nested_extract_factorized_basis(
                coefficient_matrix,
                axis_counts;
                verify_reconstruction = false,
                timing_prefix =
                    "decomposed_wl.factorized_retained_basis.bridge",
            )
        end
        (;
            status = :materialized_decomposed_wl_factorized_retained_basis,
            blocker = nothing,
            factorized_basis,
            coefficient_matrix,
            retained_dimension = inventory.retained_dimension,
            coefficient_cache_entry_count = length(coefficient_cache),
            materialization_path =
                :shellification_retained_units_to_factorized_basis_bridge,
            old_fixed_block_matrix_authority_used = false,
            full_parent_window_cpb_used = false,
            direct_cartesian_product_assembly_used = false,
            ordinary_cartesian_ida_operators_used = false,
        )
    catch err
        (;
            status = :blocked_decomposed_wl_factorized_retained_basis,
            blocker = :decomposed_wl_factorized_basis_extraction_failed,
            error = sprint(showerror, err),
            factorized_basis = nothing,
            coefficient_matrix = nothing,
            retained_dimension =
                _route_global_one_body_value(inventory, :retained_dimension, nothing),
            coefficient_cache_entry_count = 0,
            materialization_path =
                :shellification_retained_units_to_factorized_basis_bridge,
            old_fixed_block_matrix_authority_used = false,
            full_parent_window_cpb_used = false,
            direct_cartesian_product_assembly_used = false,
            ordinary_cartesian_ida_operators_used = false,
        )
    end
    cache[cache_key] = result
    return result
end

function _white_lindsey_decomposed_retained_coefficient_matrix(
    inventory,
    axis_counts::NTuple{3,Int},
    coefficient_cache,
)
    parent_dimension = prod(axis_counts)
    retained_dimension = Int(inventory.retained_dimension)
    matrix = zeros(Float64, parent_dimension, retained_dimension)
    for unit in inventory.retained_units
        unit_coefficients =
            _white_lindsey_unit_coefficients_from_local_cache(
                coefficient_cache,
                unit,
            )
        _white_lindsey_unit_coefficients_materialized(unit_coefficients) ||
            throw(
                ArgumentError(
                    "factorized decomposed WL backend requires materialized retained-unit coefficients",
                ),
            )
        _white_lindsey_insert_unit_coefficients!(
            matrix,
            unit,
            unit_coefficients,
            parent_dimension,
        )
    end
    return matrix
end

function _white_lindsey_insert_unit_coefficients!(
    destination::Matrix{Float64},
    unit,
    unit_coefficients,
    parent_dimension::Int,
)
    columns = unit.column_range
    coefficient_matrix = unit_coefficients.coefficient_matrix
    size(coefficient_matrix, 2) == length(columns) || throw(
        ArgumentError("retained-unit coefficient column count does not match retained range"),
    )
    if unit_coefficients.coefficient_space === :parent_cartesian_sparse_adapter &&
       size(coefficient_matrix, 1) == parent_dimension
        destination[:, columns] .= coefficient_matrix
        return destination
    end
    support_indices = unit_coefficients.support_indices
    support_indices isa AbstractVector || throw(
        ArgumentError("retained-unit coefficient insertion requires parent support indices"),
    )
    size(coefficient_matrix, 1) == length(support_indices) || throw(
        ArgumentError("support-local retained-unit coefficient row count mismatch"),
    )
    @inbounds for local_column in axes(coefficient_matrix, 2)
        global_column = first(columns) + local_column - 1
        for local_row in axes(coefficient_matrix, 1)
            global_row = Int(support_indices[local_row])
            destination[global_row, global_column] =
                Float64(coefficient_matrix[local_row, local_column])
        end
    end
    return destination
end

function _route_global_decomposed_wl_factorized_one_body_matrix(
    inventory,
    term::Symbol;
    parent_axis_counts,
    overlap_1d = nothing,
    position_1d = nothing,
    x2_1d = nothing,
    kinetic_1d = nothing,
    metadata = (;),
)
    axis_counts = _axis_counts_tuple(parent_axis_counts)
    overlap_axes = _overlap_1d_tuple(overlap_1d)
    _assert_overlap_axis_sizes(overlap_axes, axis_counts)
    position_axes = term in (:position_x, :position_y, :position_z) ?
                    _operator_1d_tuple(position_1d, "position_1d") :
                    nothing
    term in (:position_x, :position_y, :position_z) &&
        _assert_operator_axis_sizes(position_axes, axis_counts, "position_1d")
    x2_axes = term in (:x2_x, :x2_y, :x2_z) ?
              _operator_1d_tuple(x2_1d, "x2_1d") :
              nothing
    term in (:x2_x, :x2_y, :x2_z) &&
        _assert_operator_axis_sizes(x2_axes, axis_counts, "x2_1d")
    kinetic_axes = term === :kinetic ? _operator_1d_tuple(
        kinetic_1d,
        "kinetic_1d",
    ) : nothing
    term === :kinetic &&
        _assert_operator_axis_sizes(kinetic_axes, axis_counts, "kinetic_1d")
    factorized =
        _white_lindsey_decomposed_factorized_retained_basis(
            inventory,
            axis_counts,
        )
    factorized.status === :materialized_decomposed_wl_factorized_retained_basis ||
        return _route_global_decomposed_wl_factorized_one_body_blocked(
            term,
            factorized.blocker,
            inventory;
            metadata,
        )
    matrix = @timeg "decomposed_wl.factorized_retained_basis.one_body_fill" begin
        axis_overlap_x = _white_lindsey_factorized_axis_matrix_table(
            factorized.factorized_basis.x_functions,
            overlap_axes[1],
        )
        axis_overlap_y = _white_lindsey_factorized_axis_matrix_table(
            factorized.factorized_basis.y_functions,
            overlap_axes[2],
        )
        axis_overlap_z = _white_lindsey_factorized_axis_matrix_table(
            factorized.factorized_basis.z_functions,
            overlap_axes[3],
        )
        if term === :overlap
            ParentGaussletBases._nested_factorized_product_matrix(
                factorized.factorized_basis,
                axis_overlap_x,
                axis_overlap_y,
                axis_overlap_z,
            )
        elseif term === :kinetic
            axis_kinetic_x = _white_lindsey_factorized_axis_matrix_table(
                factorized.factorized_basis.x_functions,
                kinetic_axes[1],
            )
            axis_kinetic_y = _white_lindsey_factorized_axis_matrix_table(
                factorized.factorized_basis.y_functions,
                kinetic_axes[2],
            )
            axis_kinetic_z = _white_lindsey_factorized_axis_matrix_table(
                factorized.factorized_basis.z_functions,
                kinetic_axes[3],
            )
            ParentGaussletBases._nested_factorized_sum_of_products(
                factorized.factorized_basis,
                (
                    (axis_kinetic_x, axis_overlap_y, axis_overlap_z),
                    (axis_overlap_x, axis_kinetic_y, axis_overlap_z),
                    (axis_overlap_x, axis_overlap_y, axis_kinetic_z),
                ),
            )
        elseif term in (:position_x, :position_y, :position_z)
            ParentGaussletBases._nested_factorized_product_matrix(
                factorized.factorized_basis,
                _route_global_decomposed_wl_axis_factor_table(
                    factorized.factorized_basis.x_functions,
                    _route_global_decomposed_wl_axis_product_axes(
                        term,
                        overlap_axes,
                        position_axes,
                    )[1],
                ),
                _route_global_decomposed_wl_axis_factor_table(
                    factorized.factorized_basis.y_functions,
                    _route_global_decomposed_wl_axis_product_axes(
                        term,
                        overlap_axes,
                        position_axes,
                    )[2],
                ),
                _route_global_decomposed_wl_axis_factor_table(
                    factorized.factorized_basis.z_functions,
                    _route_global_decomposed_wl_axis_product_axes(
                        term,
                        overlap_axes,
                        position_axes,
                    )[3],
                ),
            )
        elseif term in (:x2_x, :x2_y, :x2_z)
            ParentGaussletBases._nested_factorized_product_matrix(
                factorized.factorized_basis,
                _route_global_decomposed_wl_axis_factor_table(
                    factorized.factorized_basis.x_functions,
                    _route_global_decomposed_wl_axis_product_axes(
                        term,
                        overlap_axes,
                        x2_axes,
                    )[1],
                ),
                _route_global_decomposed_wl_axis_factor_table(
                    factorized.factorized_basis.y_functions,
                    _route_global_decomposed_wl_axis_product_axes(
                        term,
                        overlap_axes,
                        x2_axes,
                    )[2],
                ),
                _route_global_decomposed_wl_axis_factor_table(
                    factorized.factorized_basis.z_functions,
                    _route_global_decomposed_wl_axis_product_axes(
                        term,
                        overlap_axes,
                        x2_axes,
                    )[3],
                ),
            )
        else
            throw(ArgumentError("unsupported factorized decomposed WL one-body term $(term)"))
        end
    end
    stream_state =
        _route_global_decomposed_wl_factorized_stream_state(
            inventory,
            term,
            factorized,
            _white_lindsey_factorized_symmetric_placed_block_count(inventory),
        )
    global_matrix_result =
        _route_global_decomposed_wl_streaming_global_matrix_result(
            term,
            inventory,
            matrix,
            stream_state,
        )
    return (;
        status = :materialized_decomposed_wl_factorized_one_body_matrix,
        blocker = nothing,
        local_batch =
            _route_global_decomposed_wl_streaming_local_batch(term, stream_state),
        local_collection =
            _route_global_decomposed_wl_streaming_local_collection(term, stream_state),
        placement_plan =
            _route_global_decomposed_wl_streaming_placement_plan(
                term,
                inventory,
                stream_state,
            ),
        global_matrix_result,
        metadata = NamedTuple(metadata),
    )
end

function _route_global_decomposed_wl_axis_product_axes(
    term::Symbol,
    overlap_axes,
    selected_axes,
)
    term in (:position_x, :x2_x) &&
        return (selected_axes[1], overlap_axes[2], overlap_axes[3])
    term in (:position_y, :x2_y) &&
        return (overlap_axes[1], selected_axes[2], overlap_axes[3])
    term in (:position_z, :x2_z) &&
        return (overlap_axes[1], overlap_axes[2], selected_axes[3])
    throw(ArgumentError("unsupported axis-product decomposed WL term $(term)"))
end

function _route_global_decomposed_wl_axis_factor_table(axis_functions, operator)
    return _white_lindsey_factorized_axis_matrix_table(axis_functions, operator)
end

function _white_lindsey_factorized_axis_matrix_table(
    axis_functions::AbstractMatrix{<:Real},
    operator::AbstractMatrix{<:Real},
)
    scratch = Matrix{Float64}(
        undef,
        size(axis_functions, 2),
        size(axis_functions, 1),
    )
    return ParentGaussletBases._nested_factorized_axis_matrix_table(
        operator,
        axis_functions,
        scratch,
    )
end

function _route_global_decomposed_wl_factorized_one_body_blocked(
    term::Symbol,
    blocker,
    inventory;
    metadata,
)
    return (;
        status = :blocked_decomposed_wl_factorized_one_body_matrix,
        blocker,
        local_batch = nothing,
        local_collection = nothing,
        placement_plan = nothing,
        global_matrix_result = nothing,
        metadata = NamedTuple(metadata),
        inventory_source_kind =
            _route_global_one_body_value(inventory, :source_kind, nothing),
    )
end

function _route_global_decomposed_wl_factorized_stream_state(
    inventory,
    term::Symbol,
    factorized,
    placed_block_count::Int;
    center_summary = nothing,
)
    return (;
        materialized_count = Int(inventory.pair_count),
        skipped_count = 0,
        placed_block_count,
        blocker = nothing,
        unit_coefficient_cache_entry_count =
            factorized.coefficient_cache_entry_count,
        electron_nuclear_prepared_unit_cache_entry_count = 0,
        pair_count = Int(inventory.pair_count),
        center_index =
            isnothing(center_summary) ? nothing : center_summary.center_index,
        center_key =
            isnothing(center_summary) ? nothing : center_summary.center_key,
        center_location =
            isnothing(center_summary) ? nothing : center_summary.location,
        nuclear_charge =
            isnothing(center_summary) ? nothing : center_summary.charge,
        nuclear_charge_applied = false,
        materialization_path =
            term === :electron_nuclear_by_center ?
            :white_lindsey_decomposed_wl_factorized_electron_nuclear_global_assembly :
            :white_lindsey_decomposed_wl_factorized_one_body_global_assembly,
        local_blocks_collected = false,
        pair_lookup_materialized = false,
        placement_records_materialized = false,
        factorized_retained_basis_materialized = true,
        old_fixed_block_matrix_authority_used = false,
        full_parent_window_cpb_used = false,
        direct_cartesian_product_assembly_used = false,
        ordinary_cartesian_ida_operators_used = false,
    )
end

function _white_lindsey_factorized_symmetric_placed_block_count(inventory)
    unit_count = Int(inventory.unit_count)
    pair_count = Int(inventory.pair_count)
    return 2 * pair_count - unit_count
end

function _white_lindsey_pair_product_block_from_prepared_units(
    left_unit,
    right_unit,
    operator_axes,
)
    support_block =
        Matrix{Float64}(
            undef,
            length(left_unit.support_states),
            length(right_unit.support_states),
        )
    _fill_source_mode_product_block!(
        support_block,
        left_unit.support_states,
        right_unit.support_states,
        operator_axes[1],
        operator_axes[2],
        operator_axes[3],
    )
    return Matrix{Float64}(
        transpose(left_unit.support_coefficients) *
        support_block *
        right_unit.support_coefficients,
    )
end

function _white_lindsey_pair_kinetic_block_from_prepared_units(
    left_unit,
    right_unit,
    overlap_axes,
    kinetic_axes,
)
    block =
        _white_lindsey_pair_product_block_from_prepared_units(
            left_unit,
            right_unit,
            (kinetic_axes[1], overlap_axes[2], overlap_axes[3]),
        )
    block .+= _white_lindsey_pair_product_block_from_prepared_units(
        left_unit,
        right_unit,
        (overlap_axes[1], kinetic_axes[2], overlap_axes[3]),
    )
    block .+= _white_lindsey_pair_product_block_from_prepared_units(
        left_unit,
        right_unit,
        (overlap_axes[1], overlap_axes[2], kinetic_axes[3]),
    )
    return block
end

function _route_global_decomposed_wl_streaming_placement_blocker(
    block,
    rows,
    columns,
    dimension::Int,
)
    block isa AbstractMatrix{<:Real} || return :missing_local_block_result
    (
        _one_body_placement_valid_column_range(rows) &&
        _one_body_placement_valid_column_range(columns)
    ) || return :missing_column_ranges
    size(block) == (length(rows), length(columns)) ||
        return :local_block_shape_mismatch
    _one_body_placement_ranges_inside_global_dimension(rows, columns, dimension) ||
        return :target_ranges_outside_global_dimension
    return nothing
end

function _route_global_decomposed_wl_streaming_local_batch(
    term::Symbol,
    stream_state,
)
    materialized = stream_state.materialized_count > 0
    return (;
        object_kind = :white_lindsey_boundary_stratum_one_body_streaming_batch,
        term,
        materialized_count = stream_state.materialized_count,
        skipped_count = stream_state.skipped_count,
        materialized,
        source_operator_blocks_materialized = materialized,
        final_pair_blocks_materialized = materialized,
        operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
        metadata = (;
            materialization_path = stream_state.materialization_path,
            pair_input_kind = :unit_pair_index_table,
            pair_unit_coefficient_record_count = stream_state.pair_count,
            unit_coefficient_cache_entry_count =
                stream_state.unit_coefficient_cache_entry_count,
            local_blocks_collected = false,
            pair_lookup_materialized = false,
            placement_records_materialized = false,
            streaming_global_matrix_insertion = true,
        ),
    )
end

function _route_global_decomposed_wl_streaming_local_collection(
    term::Symbol,
    stream_state,
)
    materialized = stream_state.materialized_count > 0
    return (;
        object_kind = :cartesian_pair_block_local_one_body_streaming_collection,
        status =
            materialized ?
            :materialized_local_one_body_block_collection :
            :blocked_local_one_body_block_collection,
        blocker =
            materialized ? nothing : Symbol("no_materialized_", String(term), "_blocks"),
        terms = (term,),
        requested_terms = (term,),
        requested_materialize_terms = (term,),
        materialized_terms = materialized ? (term,) : (),
        deferred_terms = (),
        term_statuses = (),
        entries = (),
        materialized_entries = (),
        skipped_entries = (),
        entry_count = stream_state.materialized_count + stream_state.skipped_count,
        materialized_entry_count = stream_state.materialized_count,
        skipped_entry_count = stream_state.skipped_count,
        final_local_entry_count = stream_state.materialized_count,
        local_blocks_collected = false,
        pair_lookup_materialized = false,
        placement_records_materialized = false,
        streaming_global_matrix_insertion = true,
        block_matrices_copied_into_collection = false,
        source_operator_blocks_materialized = materialized,
        final_pair_blocks_materialized = materialized,
        operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
        global_operator_blocks_materialized = false,
        route_global_matrix_materialized = true,
    )
end

function _route_global_decomposed_wl_streaming_placement_plan(
    term::Symbol,
    inventory,
    stream_state,
)
    materialized = stream_state.materialized_count > 0
    return (;
        object_kind = :cartesian_pair_block_local_one_body_streaming_placement_plan,
        term,
        status =
            materialized ?
            :placeable_local_one_body_placement_plan :
            :blocked_local_one_body_placement_plan,
        blocker = materialized ? nothing : stream_state.blocker,
        global_dimension = inventory.retained_dimension,
        global_dimension_status = inventory.retained_dimension_status,
        placeable_records = (),
        blocked_records = (),
        placeable_count = stream_state.materialized_count,
        blocked_count = stream_state.skipped_count,
        record_count = stream_state.materialized_count + stream_state.skipped_count,
        placement_records_materialized = false,
        streaming_global_matrix_insertion = true,
    )
end

function _route_global_decomposed_wl_streaming_global_matrix_result(
    term::Symbol,
    inventory,
    matrix::Matrix{Float64},
    stream_state,
)
    placement_plan =
        _route_global_decomposed_wl_streaming_placement_plan(
            term,
            inventory,
            stream_state,
        )
    materialized = stream_state.materialized_count > 0 &&
                   stream_state.skipped_count == 0
    if term === :overlap
        return _one_body_global_overlap_result(
            placement_plan,
            materialized ?
            :materialized_global_overlap_matrix :
            :blocked_global_overlap_matrix,
            materialized ? nothing : stream_state.blocker,
            materialized ? matrix : nothing;
            placed_block_count = stream_state.placed_block_count,
            skipped_block_count = stream_state.skipped_count,
            materialized,
        )
    elseif term === :kinetic
        return _one_body_global_kinetic_result(
            placement_plan,
            materialized ?
            :materialized_global_kinetic_matrix :
            :blocked_global_kinetic_matrix,
            materialized ? nothing : stream_state.blocker,
            materialized ? matrix : nothing;
            placed_block_count = stream_state.placed_block_count,
            skipped_block_count = stream_state.skipped_count,
            materialized,
        )
    elseif term in (:position_x, :position_y, :position_z)
        return _one_body_global_position_result(
            placement_plan,
            term,
            materialized ?
            Symbol("materialized_global_", String(term), "_matrix") :
            Symbol("blocked_global_", String(term), "_matrix"),
            materialized ? nothing : stream_state.blocker,
            materialized ? matrix : nothing;
            placed_block_count = stream_state.placed_block_count,
            skipped_block_count = stream_state.skipped_count,
            materialized,
        )
    elseif term in (:x2_x, :x2_y, :x2_z)
        return _one_body_global_x2_result(
            placement_plan,
            term,
            materialized ?
            Symbol("materialized_global_", String(term), "_matrix") :
            Symbol("blocked_global_", String(term), "_matrix"),
            materialized ? nothing : stream_state.blocker,
            materialized ? matrix : nothing;
            placed_block_count = stream_state.placed_block_count,
            skipped_block_count = stream_state.skipped_count,
            materialized,
        )
    elseif term === :electron_nuclear_by_center
        return _one_body_global_electron_nuclear_by_center_result(
            placement_plan,
            materialized ?
            :materialized_global_electron_nuclear_by_center_matrix :
            :blocked_global_electron_nuclear_by_center_matrix,
            materialized ? nothing : stream_state.blocker,
            materialized ? matrix : nothing;
            center_index = stream_state.center_index,
            center_key = stream_state.center_key,
            center_location = stream_state.center_location,
            nuclear_charge = stream_state.nuclear_charge,
            placed_block_count = stream_state.placed_block_count,
            skipped_block_count = stream_state.skipped_count,
            materialized,
        )
    end
    throw(ArgumentError("unsupported streaming decomposed WL one-body term $(term)"))
end

function _route_global_decomposed_wl_one_body_result(
    inventory,
    local_batch,
    local_collection::NamedTuple,
    placement_plan::NamedTuple,
    global_matrix_result::NamedTuple;
    term::Symbol,
    metadata,
)
    materialized =
        _route_global_one_body_value(
            global_matrix_result,
            :operator_matrix_materialized,
            false,
        ) === true
    return (;
        object_kind =
            :cartesian_pair_block_route_global_decomposed_wl_one_body_matrix_adapter,
        term,
        status =
            materialized ?
            _route_global_one_body_materialized_status(term) :
            _route_global_one_body_blocked_status(term),
        blocker = materialized ? nothing : global_matrix_result.blocker,
        inventory_source_kind = inventory.source_kind,
        decomposed_unit_pair_inventory_status = inventory.status,
        retained_dimension = inventory.retained_dimension,
        retained_dimension_status = inventory.retained_dimension_status,
        unit_count = inventory.unit_count,
        pair_count = inventory.pair_count,
        local_batch,
        local_collection,
        placement_plan,
        global_matrix_result,
        matrix = global_matrix_result.matrix,
        local_pair_block_count = local_batch.materialized_count,
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
        operator_matrix_materialized = materialized,
        global_operator_assembled = materialized,
        full_parent_window_cpb_used = false,
        direct_cartesian_product_assembly_used = false,
        ordinary_cartesian_ida_operators_used = false,
        metadata = NamedTuple(metadata),
        _route_global_one_body_nonclaim_flags()...,
    )
end

function _route_global_decomposed_wl_one_body_blocked_result(
    term::Symbol,
    blocker::Symbol;
    inventory = nothing,
    metadata,
)
    return (;
        object_kind =
            :cartesian_pair_block_route_global_decomposed_wl_one_body_matrix_adapter,
        term,
        status = _route_global_one_body_blocked_status(term),
        blocker,
        inventory_source_kind =
            _route_global_one_body_value(inventory, :source_kind, nothing),
        decomposed_unit_pair_inventory_status =
            _route_global_one_body_value(inventory, :status, nothing),
        retained_dimension = nothing,
        retained_dimension_status = :not_available,
        unit_count = 0,
        pair_count = 0,
        local_batch = nothing,
        local_collection = nothing,
        placement_plan = nothing,
        global_matrix_result = nothing,
        matrix = nothing,
        local_pair_block_count = 0,
        placeable_record_count = 0,
        blocked_placement_count = 0,
        placed_block_count = 0,
        skipped_block_count = 0,
        global_one_body_term_matrix_materialized = false,
        global_overlap_matrix_materialized = false,
        global_kinetic_matrix_materialized = false,
        operator_matrix_materialized = false,
        global_operator_assembled = false,
        full_parent_window_cpb_used = false,
        direct_cartesian_product_assembly_used = false,
        ordinary_cartesian_ida_operators_used = false,
        metadata = NamedTuple(metadata),
        _route_global_one_body_nonclaim_flags()...,
    )
end

function _white_lindsey_decomposed_one_electron_hamiltonian_blocked_result(
    blocker::Symbol;
    metadata,
)
    return (;
        object_kind =
            :white_lindsey_decomposed_one_electron_hamiltonian_assembly,
        status = :blocked_white_lindsey_decomposed_one_electron_hamiltonian,
        blocker,
        term = :one_electron_hamiltonian,
        matrix = nothing,
        matrix_shape = :not_materialized,
        retained_dimension = nothing,
        kinetic_status = nothing,
        kinetic_matrix_materialized = false,
        electron_nuclear_by_center_status = nothing,
        by_center_matrix_count = 0,
        applied_center_count = 0,
        center_indices = (),
        nuclear_charges = (),
        unit_charge_nuclear_attraction_matrices_consumed = false,
        nuclear_charge_application_convention =
            :multiply_unit_charge_nuclear_attraction_by_recorded_charge,
        charge_application_stage = :white_lindsey_hamiltonian_assembly,
        nuclear_charge_applied_at_hamiltonian_assembly = false,
        center_summation_stage = :white_lindsey_hamiltonian_assembly,
        centers_summed_at_hamiltonian_assembly = false,
        by_center_matrices_remain_separated = false,
        overlap_matrix_consumed = false,
        overlap_used_as_solve_metric_only = true,
        full_parent_window_cpb_used = false,
        direct_cartesian_product_assembly_used = false,
        ordinary_cartesian_ida_operators_used = false,
        metadata = NamedTuple(metadata),
        _route_global_one_body_nonclaim_flags()...,
    )
end

function route_global_safe_one_body_terms()
    return _route_global_one_body_safe_terms()
end

function _route_global_one_body_safe_terms()
    return (
        :overlap,
        :kinetic,
        :position_x,
        :position_y,
        :position_z,
        :x2_x,
        :x2_y,
        :x2_z,
    )
end

function _route_global_one_body_supported_terms()
    return (_route_global_one_body_safe_terms()..., :electron_nuclear_by_center)
end

function _route_global_one_body_term_tuple(terms::Symbol)
    return (terms,)
end

function _route_global_one_body_term_tuple(terms)
    term_tuple = Tuple(terms)
    all(term -> term isa Symbol, term_tuple) ||
        throw(ArgumentError("route global one-body terms must be Symbols"))
    return term_tuple
end

function _route_global_one_body_matrix_set_summary(
    terms::Tuple,
    term_results::Tuple,
)
    materialized_terms = Tuple(
        result.term for result in term_results
        if result.global_one_body_term_matrix_materialized
    )
    blocked_results = Tuple(
        result for result in term_results
        if !result.global_one_body_term_matrix_materialized
    )
    blocked_terms = Tuple(result.term for result in blocked_results)
    materialized_count = length(materialized_terms)
    blocked_count = length(blocked_terms)
    term_count = length(terms)
    all_materialized = term_count > 0 && blocked_count == 0
    any_materialized = materialized_count > 0

    return (;
        object_kind =
            :cartesian_pair_block_route_global_safe_one_body_matrix_set_summary,
        status =
            all_materialized ?
            :materialized_route_global_safe_one_body_matrix_set :
            any_materialized ?
            :partial_route_global_safe_one_body_matrix_set :
            :blocked_route_global_safe_one_body_matrix_set,
        blocker =
            all_materialized ?
            nothing :
            any_materialized ?
            :some_route_global_safe_one_body_terms_blocked :
            _route_global_one_body_matrix_set_blocker(blocked_results),
        terms,
        materialized_terms,
        blocked_terms,
        term_count,
        materialized_term_count = materialized_count,
        blocked_term_count = blocked_count,
        all_terms_materialized = all_materialized,
        any_terms_materialized = any_materialized,
        term_status_counts =
            _count_by_value((result.status for result in term_results), :status),
        blocked_term_blocker_counts =
            _count_by_value(
                (result.blocker for result in blocked_results),
                :blocker,
            ),
        global_one_body_term_matrices_materialized = all_materialized,
        any_global_one_body_term_matrix_materialized = any_materialized,
        _route_global_one_body_nonclaim_flags()...,
    )
end

function _route_global_one_body_matrix_set_blocker(blocked_results::Tuple)
    isempty(blocked_results) && return nothing
    for result in blocked_results
        !isnothing(result.blocker) && return result.blocker
    end
    return :no_route_global_safe_one_body_terms_materialized
end

function _route_state_global_one_body_matrix(
    source,
    term::Symbol;
    global_dimension = nothing,
    inputs = (;),
    provider = nothing,
    factors = nothing,
    factor_provider = nothing,
    metadata = (;),
)
    plan = _route_state_global_pair_block_materialization_plan(source)
    isnothing(plan) &&
        return _route_global_one_body_blocked_result(
            term,
            :missing_pair_block_materialization_plan;
            metadata,
        )

    return route_global_one_body_matrix(
        plan;
        term,
        global_dimension,
        inputs,
        provider,
        factors,
        factor_provider,
        metadata,
    )
end

function _route_state_global_pair_block_materialization_plan(
    plan::PairBlockMaterializationPlan,
)
    return plan
end

function _route_state_global_pair_block_materialization_plan(source)
    if hasproperty(source, :pair_block_materialization_plan)
        plan = getproperty(source, :pair_block_materialization_plan)
        plan isa PairBlockMaterializationPlan && return plan
        return nothing
    end

    if hasproperty(source, :terminal_route_state)
        route_state = getproperty(source, :terminal_route_state)
        return _route_state_global_pair_block_materialization_plan(route_state)
    end

    return nothing
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
    term === :x2_x && return one_body_x2_x_placement_plan(
        collection;
        global_dimension,
    )
    term === :x2_y && return one_body_x2_y_placement_plan(
        collection;
        global_dimension,
    )
    term === :x2_z && return one_body_x2_z_placement_plan(
        collection;
        global_dimension,
    )
    term === :electron_nuclear_by_center &&
        return one_body_electron_nuclear_by_center_placement_plan(
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
    term === :x2_x && return one_body_global_x2_x_matrix(placement_plan)
    term === :x2_y && return one_body_global_x2_y_matrix(placement_plan)
    term === :x2_z && return one_body_global_x2_z_matrix(placement_plan)
    term === :electron_nuclear_by_center &&
        return one_body_global_electron_nuclear_by_center_matrix(placement_plan)
    throw(ArgumentError("unsupported route global one-body term: $(term)"))
end

function _route_global_one_body_pair_block_source(source)
    source isa PairBlockMaterializationPlan && return true
    hasproperty(source, :pair_block_materialization_plan) && return true
    hasproperty(source, :terminal_route_state) && return true
    return false
end

function _route_global_electron_nuclear_by_center_records(center_record::NamedTuple)
    return (center_record,)
end

function _route_global_electron_nuclear_by_center_records(center_records)
    isnothing(center_records) && return ()
    return Tuple(center_records)
end

function _route_global_electron_nuclear_by_center_inventory(
    source;
    parent_axis_counts = nothing,
    parent_axis_bundle_object = nothing,
)
    return @timeg "decomposed_wl.inventory" begin
        if _route_global_one_body_value(source, :object_kind, nothing) ===
           :white_lindsey_decomposed_unit_pair_inventory
            return source
        end
        return white_lindsey_decomposed_unit_pair_inventory(
            source;
            metadata = (;
                parent_axis_counts,
                parent_axis_bundle_object,
            ),
        )
    end
end

function _route_global_electron_nuclear_by_center_matrix_from_inventory(
    source;
    parent_axis_counts,
    parent_axis_bundle_object,
    coulomb_expansion,
    center_record,
    metadata,
)
    inventory = _route_global_electron_nuclear_by_center_inventory(
        source;
        parent_axis_counts,
        parent_axis_bundle_object,
    )
    if inventory.status !== :available_white_lindsey_decomposed_unit_pair_inventory
        return _route_global_electron_nuclear_by_center_blocked_inventory_result(
            inventory,
            inventory.blocker;
            metadata,
        )
    end

    if inventory.unit_pairs isa CUP.UnitPairIndexTable
        factorized =
            _route_global_decomposed_wl_factorized_electron_nuclear_by_center_matrix(
                inventory;
                parent_axis_counts,
                parent_axis_bundle_object,
                coulomb_expansion,
                center_record,
                metadata,
            )
        if factorized.status ===
           :materialized_decomposed_wl_factorized_electron_nuclear_by_center_matrix
            return _route_global_electron_nuclear_by_center_result(
                inventory,
                factorized.local_batch,
                factorized.local_collection,
                factorized.placement_plan,
                factorized.global_matrix_result;
                metadata,
            )
        end
        streaming =
            _route_global_decomposed_wl_streaming_electron_nuclear_by_center_matrix(
                inventory;
                parent_axis_counts,
                parent_axis_bundle_object,
                coulomb_expansion,
                center_record,
                metadata,
            )
        return _route_global_electron_nuclear_by_center_result(
            inventory,
            streaming.local_batch,
            streaming.local_collection,
            streaming.placement_plan,
            streaming.global_matrix_result;
            metadata,
        )
    end

    local_batch = @timeg "decomposed_wl.electron_nuclear_by_center.local_batch" begin
        white_lindsey_boundary_stratum_one_body_blocks(
            inventory.unit_pairs,
            :electron_nuclear_by_center;
            parent_axis_counts,
            parent_axis_bundle_object,
            coulomb_expansion,
            center_record,
        )
    end
    local_collection =
        _route_global_electron_nuclear_by_center_local_collection(
            local_batch,
            inventory.unit_pairs,
        )
    placement_plan = nothing
    global_matrix_result = nothing
    @timeg "decomposed_wl.electron_nuclear_by_center.global_placement" begin
        placement_plan = one_body_electron_nuclear_by_center_placement_plan(
            local_collection;
            global_dimension = inventory.retained_dimension,
        )
        global_matrix_result =
            one_body_global_electron_nuclear_by_center_matrix(placement_plan)
    end
    return _route_global_electron_nuclear_by_center_result(
        inventory,
        local_batch,
        local_collection,
        placement_plan,
        global_matrix_result;
        metadata,
    )
end

function _route_global_decomposed_wl_streaming_electron_nuclear_by_center_matrix(
    inventory;
    parent_axis_counts,
    parent_axis_bundle_object,
    coulomb_expansion,
    center_record,
    metadata = (;),
)
    dimension = inventory.retained_dimension
    axis_context = _white_lindsey_electron_nuclear_axis_context(
        parent_axis_counts,
        parent_axis_bundle_object,
        coulomb_expansion,
        center_record,
    )
    matrix = zeros(Float64, dimension, dimension)
    coefficient_cache =
        _white_lindsey_decomposed_operator_unit_coefficient_cache(inventory)
    prepared_cache =
        _white_lindsey_decomposed_operator_prepared_unit_cache(inventory)
    scratch = _white_lindsey_electron_nuclear_streaming_scratch(
        inventory,
        coefficient_cache,
        prepared_cache,
        axis_context.axis_counts,
    )
    stream_state =
        @timeg "decomposed_wl.electron_nuclear_by_center.local_batch" begin
            _route_global_decomposed_wl_streaming_fill!(
                matrix,
                inventory.unit_pairs,
                :electron_nuclear_by_center;
                parent_axis_counts,
                parent_axis_bundle_object,
                coulomb_expansion,
                center_record,
                electron_nuclear_axis_context = axis_context,
                unit_coefficient_cache = coefficient_cache,
                prepared_unit_cache = prepared_cache,
                electron_nuclear_scratch = scratch,
            )
        end
    global_matrix_result =
        _route_global_decomposed_wl_streaming_global_matrix_result(
            :electron_nuclear_by_center,
            inventory,
            matrix,
            stream_state,
        )
    local_batch =
        _route_global_decomposed_wl_streaming_local_batch(
            :electron_nuclear_by_center,
            stream_state,
        )
    local_collection =
        _route_global_decomposed_wl_streaming_local_collection(
            :electron_nuclear_by_center,
            stream_state,
        )
    placement_plan =
        _route_global_decomposed_wl_streaming_placement_plan(
            :electron_nuclear_by_center,
            inventory,
            stream_state,
        )
    return (;
        local_batch,
        local_collection,
        placement_plan,
        global_matrix_result,
        metadata = NamedTuple(metadata),
    )
end

function _route_global_decomposed_wl_factorized_electron_nuclear_by_center_matrix(
    inventory;
    parent_axis_counts,
    parent_axis_bundle_object,
    coulomb_expansion,
    center_record,
    metadata = (;),
)
    axis_counts = _axis_counts_tuple(parent_axis_counts)
    factorized =
        _white_lindsey_decomposed_factorized_retained_basis(
            inventory,
            axis_counts,
        )
    factorized.status === :materialized_decomposed_wl_factorized_retained_basis ||
        return (;
            status =
                :blocked_decomposed_wl_factorized_electron_nuclear_by_center_matrix,
            blocker = factorized.blocker,
            local_batch = nothing,
            local_collection = nothing,
            placement_plan = nothing,
            global_matrix_result = nothing,
            metadata = NamedTuple(metadata),
        )
    axis_context = _white_lindsey_electron_nuclear_axis_context(
        parent_axis_counts,
        parent_axis_bundle_object,
        coulomb_expansion,
        center_record,
    )
    matrix = @timeg "decomposed_wl.factorized_retained_basis.electron_nuclear_fill" begin
        axis_term_tables_x =
            ParentGaussletBases._nested_factorized_axis_term_tables(
                axis_context.axis_terms.x,
                factorized.factorized_basis.x_functions,
            )
        axis_term_tables_y =
            ParentGaussletBases._nested_factorized_axis_term_tables(
                axis_context.axis_terms.y,
                factorized.factorized_basis.y_functions,
            )
        axis_term_tables_z =
            ParentGaussletBases._nested_factorized_axis_term_tables(
                axis_context.axis_terms.z,
                factorized.factorized_basis.z_functions,
            )
        matrix_local = Matrix{Float64}(
            undef,
            inventory.retained_dimension,
            inventory.retained_dimension,
        )
        ParentGaussletBases._nested_fill_factorized_weighted_term_sum!(
            matrix_local,
            factorized.factorized_basis,
            axis_context.coefficients,
            axis_term_tables_x,
            axis_term_tables_y,
            axis_term_tables_z;
            include_basis_amplitudes = true,
        )
        matrix_local
    end
    stream_state =
        _route_global_decomposed_wl_factorized_stream_state(
            inventory,
            :electron_nuclear_by_center,
            factorized,
            _white_lindsey_factorized_symmetric_placed_block_count(inventory);
            center_summary = axis_context.center_summary,
        )
    global_matrix_result =
        _route_global_decomposed_wl_streaming_global_matrix_result(
            :electron_nuclear_by_center,
            inventory,
            matrix,
            stream_state,
        )
    return (;
        status =
            :materialized_decomposed_wl_factorized_electron_nuclear_by_center_matrix,
        blocker = nothing,
        local_batch =
            _route_global_decomposed_wl_streaming_local_batch(
                :electron_nuclear_by_center,
                stream_state,
            ),
        local_collection =
            _route_global_decomposed_wl_streaming_local_collection(
                :electron_nuclear_by_center,
                stream_state,
            ),
        placement_plan =
            _route_global_decomposed_wl_streaming_placement_plan(
                :electron_nuclear_by_center,
                inventory,
                stream_state,
            ),
        global_matrix_result,
        metadata = NamedTuple(metadata),
    )
end

function _route_global_electron_nuclear_by_center_local_collection(
    batch_result::PairBlockMaterializationBatchResult,
    unit_pairs,
)
    pair_lookup = Dict(pair.pair_key => pair for pair in unit_pairs)
    materialized_entries = Tuple(
        _route_global_electron_nuclear_by_center_collection_entry(
            result,
            pair_lookup,
        ) for result in batch_result.materialized_results
    )
    skipped_entries = Tuple(
        _one_body_local_block_collection_skipped_entry(
            skip;
            block_set_term = :electron_nuclear_by_center,
        ) for skip in batch_result.skipped_records
    )
    entries = (materialized_entries..., skipped_entries...)
    return (;
        object_kind = :cartesian_pair_block_local_one_body_block_collection,
        status =
            isempty(materialized_entries) ?
            :blocked_local_one_body_block_collection :
            :materialized_local_one_body_block_collection,
        blocker =
            isempty(materialized_entries) ?
            :no_materialized_electron_nuclear_by_center_blocks :
            nothing,
        block_set_status =
            isempty(materialized_entries) ?
            :blocked_mixed_one_body_block_set_consumption :
            :materialized_mixed_one_body_block_set_consumption,
        block_set_blocker =
            isempty(materialized_entries) ?
            :no_materialized_electron_nuclear_by_center_blocks :
            nothing,
        terms = (:electron_nuclear_by_center,),
        requested_terms = (:electron_nuclear_by_center,),
        requested_materialize_terms = (:electron_nuclear_by_center,),
        materialized_terms =
            isempty(materialized_entries) ? () : (:electron_nuclear_by_center,),
        deferred_terms = (),
        term_statuses = (),
        entries,
        materialized_entries,
        skipped_entries,
        entry_count = length(entries),
        materialized_entry_count = length(materialized_entries),
        skipped_entry_count = length(skipped_entries),
        deferred_term_count = 0,
        source_space_entry_count = 0,
        final_local_entry_count = length(materialized_entries),
        term_separated_entries = true,
        pair_separated_entries = true,
        block_set_results_summed = false,
        block_matrices_copied_into_collection = false,
        _one_body_local_block_collection_summary_materialization_flags(
            source_operator_blocks_materialized =
                batch_result.source_operator_blocks_materialized,
            final_pair_blocks_materialized =
                batch_result.final_pair_blocks_materialized,
        )...,
    )
end

function _route_global_electron_nuclear_by_center_collection_entry(
    result::PairBlockMaterializationResult,
    pair_lookup::Dict,
)
    pair = pair_lookup[result.pair_key]
    return merge(
        _one_body_local_block_collection_entry(
            result;
            block_set_term = :electron_nuclear_by_center,
        ),
        (;
            left_column_range = pair.left_unit.column_range,
            right_column_range = pair.right_unit.column_range,
            center_index = result.metadata.center_index,
            center_key = result.metadata.center_key,
            center_location = result.metadata.center_location,
            nuclear_charge_recorded = result.metadata.nuclear_charge_recorded,
            nuclear_charge = result.metadata.nuclear_charge,
            nuclear_charge_applied = result.metadata.nuclear_charge_applied,
        ),
    )
end

function _route_global_electron_nuclear_by_center_result(
    inventory,
    local_batch,
    local_collection::NamedTuple,
    placement_plan::NamedTuple,
    global_matrix_result::NamedTuple;
    metadata,
)
    materialized =
        _route_global_one_body_value(
            global_matrix_result,
            :global_electron_nuclear_by_center_matrix_materialized,
            false,
        ) === true
    return (;
        object_kind =
            :cartesian_pair_block_route_global_electron_nuclear_by_center_matrix_adapter,
        term = :electron_nuclear_by_center,
        status =
            materialized ?
            :materialized_route_global_electron_nuclear_by_center_matrix :
            :blocked_route_global_electron_nuclear_by_center_matrix,
        blocker = materialized ? nothing : global_matrix_result.blocker,
        inventory_source_kind = inventory.source_kind,
        decomposed_unit_pair_inventory_status = inventory.status,
        retained_dimension = inventory.retained_dimension,
        retained_dimension_status = inventory.retained_dimension_status,
        unit_count = inventory.unit_count,
        pair_count = inventory.pair_count,
        center_index = global_matrix_result.center_index,
        center_key = global_matrix_result.center_key,
        center_location = global_matrix_result.center_location,
        nuclear_charge = global_matrix_result.nuclear_charge,
        by_center = true,
        centers_summed = false,
        nuclear_charge_recorded = global_matrix_result.nuclear_charge_recorded,
        nuclear_charge_applied = false,
        charge_application_stage = :acceptance_or_hamiltonian_assembly,
        local_batch,
        local_collection,
        placement_plan,
        global_matrix_result,
        matrix = global_matrix_result.matrix,
        local_pair_block_count = local_batch.materialized_count,
        placeable_record_count = placement_plan.placeable_count,
        blocked_placement_count = placement_plan.blocked_count,
        placed_block_count = global_matrix_result.placed_block_count,
        skipped_block_count = global_matrix_result.skipped_block_count,
        global_electron_nuclear_by_center_matrix_materialized = materialized,
        global_one_body_term_matrix_materialized = materialized,
        operator_matrix_materialized = materialized,
        global_operator_assembled = materialized,
        full_parent_window_cpb_used = false,
        direct_cartesian_product_assembly_used = false,
        ordinary_cartesian_ida_operators_used = false,
        metadata = NamedTuple(metadata),
        _route_global_one_body_nonclaim_flags()...,
    )
end

function _route_global_electron_nuclear_by_center_blocked_inventory_result(
    inventory,
    blocker;
    metadata,
)
    return (;
        object_kind =
            :cartesian_pair_block_route_global_electron_nuclear_by_center_matrix_adapter,
        term = :electron_nuclear_by_center,
        status = :blocked_route_global_electron_nuclear_by_center_matrix,
        blocker,
        inventory_source_kind =
            _route_global_one_body_value(inventory, :source_kind, nothing),
        decomposed_unit_pair_inventory_status =
            _route_global_one_body_value(inventory, :status, nothing),
        retained_dimension = nothing,
        retained_dimension_status = :not_available,
        unit_count = 0,
        pair_count = 0,
        center_index = nothing,
        center_key = nothing,
        center_location = nothing,
        nuclear_charge = nothing,
        by_center = true,
        centers_summed = false,
        nuclear_charge_recorded = false,
        nuclear_charge_applied = false,
        charge_application_stage = :acceptance_or_hamiltonian_assembly,
        local_batch = nothing,
        local_collection = nothing,
        placement_plan = nothing,
        global_matrix_result = nothing,
        matrix = nothing,
        local_pair_block_count = 0,
        placeable_record_count = 0,
        blocked_placement_count = 0,
        placed_block_count = 0,
        skipped_block_count = 0,
        global_electron_nuclear_by_center_matrix_materialized = false,
        global_one_body_term_matrix_materialized = false,
        operator_matrix_materialized = false,
        global_operator_assembled = false,
        full_parent_window_cpb_used = false,
        direct_cartesian_product_assembly_used = false,
        ordinary_cartesian_ida_operators_used = false,
        metadata = NamedTuple(metadata),
        _route_global_one_body_nonclaim_flags()...,
    )
end

function _route_global_electron_nuclear_by_center_matrices_result(
    inventory,
    matrix_results::Tuple,
    status::Symbol,
    blocker;
    metadata,
)
    materialized = status ===
                   :materialized_route_global_electron_nuclear_by_center_matrix_set
    return (;
        object_kind =
            :cartesian_pair_block_route_global_electron_nuclear_by_center_matrix_set_adapter,
        term = :electron_nuclear_by_center,
        status,
        blocker,
        inventory_source_kind =
            _route_global_one_body_value(inventory, :source_kind, nothing),
        decomposed_unit_pair_inventory_status =
            _route_global_one_body_value(inventory, :status, nothing),
        retained_dimension =
            _route_global_one_body_value(inventory, :retained_dimension, nothing),
        retained_dimension_status =
            _route_global_one_body_value(
                inventory,
                :retained_dimension_status,
                :not_available,
            ),
        unit_count = _route_global_one_body_value(inventory, :unit_count, 0),
        pair_count = _route_global_one_body_value(inventory, :pair_count, 0),
        center_count = length(matrix_results),
        by_center_matrix_count = count(
            result ->
                result.global_electron_nuclear_by_center_matrix_materialized,
            matrix_results,
        ),
        center_indices = Tuple(result.center_index for result in matrix_results),
        center_keys = Tuple(result.center_key for result in matrix_results),
        nuclear_charges = Tuple(result.nuclear_charge for result in matrix_results),
        matrix_results,
        by_center = true,
        centers_summed = false,
        nuclear_charge_recorded = !isempty(matrix_results),
        nuclear_charge_applied = false,
        charge_application_stage = :acceptance_or_hamiltonian_assembly,
        route_global_by_center_matrices_materialized = materialized,
        global_electron_nuclear_by_center_matrices_materialized = materialized,
        global_one_body_term_matrix_materialized = materialized,
        hamiltonian_data_materialized = false,
        full_parent_window_cpb_used = false,
        direct_cartesian_product_assembly_used = false,
        ordinary_cartesian_ida_operators_used = false,
        metadata = NamedTuple(metadata),
        _route_global_one_body_nonclaim_flags()...,
    )
end

function _route_global_electron_nuclear_by_center_first_blocker(results::Tuple)
    for result in results
        !isnothing(result.blocker) && return result.blocker
    end
    return :route_global_electron_nuclear_by_center_matrix_blocked
end

function _route_global_electron_nuclear_by_center_inputs(
    inputs::NamedTuple,
    parent_axis_bundle_object,
    coulomb_expansion,
    center_record,
)
    extra = NamedTuple()
    isnothing(parent_axis_bundle_object) ||
        (extra = merge(extra, (; parent_axis_bundle_object)))
    isnothing(coulomb_expansion) ||
        (extra = merge(extra, (; coulomb_expansion)))
    isnothing(center_record) ||
        (extra = merge(extra, (; center_record)))
    return merge(inputs, extra)
end

function _route_global_position_term(axis)
    axis === :x && return :position_x
    axis === :y && return :position_y
    axis === :z && return :position_z
    throw(ArgumentError("route global position axis must be :x, :y, or :z"))
end

function _route_global_x2_term(axis)
    axis === :x && return :x2_x
    axis === :y && return :x2_y
    axis === :z && return :x2_z
    throw(ArgumentError("route global x2 axis must be :x, :y, or :z"))
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
        global_x2_x_matrix_materialized =
            _route_global_one_body_value(
                global_matrix_result,
                :global_x2_x_matrix_materialized,
                false,
            ),
        global_x2_y_matrix_materialized =
            _route_global_one_body_value(
                global_matrix_result,
                :global_x2_y_matrix_materialized,
                false,
            ),
        global_x2_z_matrix_materialized =
            _route_global_one_body_value(
                global_matrix_result,
                :global_x2_z_matrix_materialized,
                false,
            ),
        global_electron_nuclear_by_center_matrix_materialized =
            _route_global_one_body_value(
                global_matrix_result,
                :global_electron_nuclear_by_center_matrix_materialized,
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
        global_x2_x_matrix_materialized = false,
        global_x2_y_matrix_materialized = false,
        global_x2_z_matrix_materialized = false,
        global_electron_nuclear_by_center_matrix_materialized = false,
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

function _route_global_one_body_nonclaim_flags(;
    hamiltonian_data_materialized = false,
    global_hamiltonian_data_materialized = false,
)
    return (;
        route_driver_wiring = false,
        operator_blocks_materialized = false,
        hamiltonian_data_materialized,
        artifacts_materialized = false,
        exports_materialized = false,
        global_operator_blocks_materialized = false,
        global_hamiltonian_data_materialized,
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

# Private local one-body block collection entry vocabulary.
#
# These entries are a local view layer over existing mixed one-body results and
# skip summaries. They do not assemble terms, place global blocks, or copy
# matrix data; materialized entries keep only compact block-shape fields beside
# the original result reference.

function _one_body_local_block_collection_entry_materialization_flags(;
    source_operator_blocks_materialized = false,
    final_pair_blocks_materialized = false,
    operator_blocks_materialized = false,
    hamiltonian_data_materialized = false,
    artifacts_materialized = false,
)
    return (;
        source_operator_blocks_materialized,
        final_pair_blocks_materialized,
        operator_blocks_materialized,
        hamiltonian_data_materialized,
        artifacts_materialized,
        block_copied_into_entry = false,
        local_operator_assembled = false,
        global_operator_assembled = false,
        route_driver_wiring = false,
        coulomb_materialized = false,
        ida_mwg_data_materialized = false,
        pqs_lowdin_materialized = false,
        pqs_shell_projection_materialized = false,
        full_white_lindsey_route_assembled = false,
    )
end

function _one_body_local_block_collection_summary_materialization_flags(;
    source_operator_blocks_materialized = false,
    final_pair_blocks_materialized = false,
)
    return (;
        local_operator_assembled = false,
        global_operator_assembled = false,
        route_driver_wiring = false,
        source_operator_blocks_materialized,
        final_pair_blocks_materialized,
        operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
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

function _one_body_local_block_collection_entry(
    result::PairBlockMaterializationResult,
    ;
    block_set_term = nothing,
)
    block_space = _one_body_collection_result_block_space(result)
    return (;
        object_kind = :cartesian_pair_block_local_one_body_block_collection_entry,
        entry_kind = :materialized_result,
        term = result.term,
        block_set_term,
        result_term = result.term,
        source_space_term = block_space === :source_space ? result.term : nothing,
        pair_key = result.pair_key,
        pair_index = _one_body_collection_metadata_value(
            result.metadata,
            :pair_index,
            nothing,
        ),
        selector_family = _one_body_collection_metadata_value(
            result.metadata,
            :selector_family,
            nothing,
        ),
        materialization_path = _one_body_collection_metadata_value(
            result.metadata,
            :materialization_path,
            nothing,
        ),
        block_space,
        block_shape = size(result.block),
        block_size = length(result.block),
        status = :materialized_local_one_body_block_collection_entry,
        blocker = nothing,
        materialized = result.materialized,
        result_available = true,
        skipped_record_available = false,
        result,
        skipped_record = nothing,
        _one_body_local_block_collection_entry_materialization_flags(
            source_operator_blocks_materialized =
                result.source_operator_blocks_materialized,
            final_pair_blocks_materialized =
                result.final_pair_blocks_materialized,
            operator_blocks_materialized = result.operator_blocks_materialized,
            hamiltonian_data_materialized = result.hamiltonian_data_materialized,
            artifacts_materialized = result.artifacts_materialized,
        )...,
    )
end

function _one_body_local_block_collection_skipped_entry(
    skip::NamedTuple;
    block_set_term = _one_body_collection_value(skip, :requested_term, nothing),
)
    return (;
        object_kind = :cartesian_pair_block_local_one_body_block_collection_entry,
        entry_kind = :skipped_record,
        term = _one_body_collection_value(skip, :requested_term, nothing),
        block_set_term,
        result_term = nothing,
        source_space_term = nothing,
        pair_key = _one_body_collection_value(skip, :pair_key, nothing),
        pair_index = _one_body_collection_value(skip, :pair_index, nothing),
        selector_family = _one_body_collection_value(skip, :selector_family, nothing),
        materialization_path =
            _one_body_collection_value(skip, :materialization_path, nothing),
        block_space = :not_materialized,
        block_shape = nothing,
        block_size = 0,
        status = :skipped_local_one_body_block_collection_entry,
        blocker = _one_body_collection_value(skip, :blocker, nothing),
        materialized = false,
        result_available = false,
        skipped_record_available = true,
        result = nothing,
        skipped_record = skip,
        _one_body_local_block_collection_entry_materialization_flags()...,
    )
end

function _one_body_local_block_collection(consumption::NamedTuple)
    summary = _one_body_pair_block_set_consumption_summary(consumption)
    materialized_entries = NamedTuple[]
    skipped_entries = NamedTuple[]

    for batch_result in _one_body_block_set_consumption_batch_results(consumption)
        block_set_term = batch_result.term
        append!(
            materialized_entries,
            _one_body_local_block_collection_entry(
                result;
                block_set_term,
            ) for result in batch_result.materialized_results
        )
        append!(
            skipped_entries,
            _one_body_local_block_collection_skipped_entry(
                skip;
                block_set_term,
            ) for skip in batch_result.skipped_records
        )
    end

    materialized_tuple = Tuple(materialized_entries)
    skipped_tuple = Tuple(skipped_entries)
    entries = (materialized_tuple..., skipped_tuple...)
    return (;
        object_kind = :cartesian_pair_block_local_one_body_block_collection,
        status = _one_body_local_block_collection_status(
            summary,
            materialized_tuple,
            skipped_tuple,
        ),
        blocker = _one_body_local_block_collection_blocker(summary),
        block_set_status = summary.status,
        block_set_blocker = summary.blocker,
        terms = summary.requested_terms,
        requested_terms = summary.requested_terms,
        requested_materialize_terms = summary.requested_materialize_terms,
        materialized_terms = summary.materialized_terms,
        deferred_terms = summary.deferred_terms,
        term_statuses = summary.term_statuses,
        entries,
        materialized_entries = materialized_tuple,
        skipped_entries = skipped_tuple,
        entry_count = length(entries),
        materialized_entry_count = length(materialized_tuple),
        skipped_entry_count = length(skipped_tuple),
        deferred_term_count = length(summary.deferred_terms),
        source_space_entry_count = count(
            entry -> entry.block_space === :source_space,
            materialized_tuple,
        ),
        final_local_entry_count = count(
            entry -> entry.block_space === :final_local_space,
            materialized_tuple,
        ),
        term_separated_entries = true,
        pair_separated_entries = true,
        block_set_results_summed = false,
        _one_body_local_block_collection_summary_materialization_flags(
            source_operator_blocks_materialized = any(
                entry -> entry.source_operator_blocks_materialized,
                materialized_tuple,
            ),
            final_pair_blocks_materialized = any(
                entry -> entry.final_pair_blocks_materialized,
                materialized_tuple,
            ),
        )...,
    )
end

function _one_body_local_block_collection_summary(collection::NamedTuple)
    _one_body_assert_local_block_collection(collection)
    return (;
        object_kind = :cartesian_pair_block_local_one_body_block_collection_summary,
        collection_object_kind = collection.object_kind,
        status = collection.status,
        blocker = collection.blocker,
        terms = collection.terms,
        requested_terms = collection.requested_terms,
        requested_materialize_terms = collection.requested_materialize_terms,
        materialized_terms = collection.materialized_terms,
        deferred_terms = collection.deferred_terms,
        term_count = length(collection.terms),
        materialized_term_count = length(collection.materialized_terms),
        deferred_term_count = length(collection.deferred_terms),
        entry_count = collection.entry_count,
        materialized_entry_count = collection.materialized_entry_count,
        skipped_entry_count = collection.skipped_entry_count,
        source_space_entry_count = collection.source_space_entry_count,
        final_local_entry_count = collection.final_local_entry_count,
        not_materialized_entry_count = count(
            entry -> entry.block_space === :not_materialized,
            collection.skipped_entries,
        ),
        selector_family_counts =
            _one_body_count_optional_by(collection.entries, :selector_family),
        materialized_selector_family_counts =
            _one_body_count_optional_by(
                collection.materialized_entries,
                :selector_family,
            ),
        skipped_selector_family_counts =
            _one_body_count_optional_by(
                collection.skipped_entries,
                :selector_family,
            ),
        skipped_blocker_counts =
            _one_body_count_optional_by(collection.skipped_entries, :blocker),
        block_space_counts =
            _one_body_count_optional_by(collection.entries, :block_space),
        term_separated_entries = collection.term_separated_entries,
        pair_separated_entries = collection.pair_separated_entries,
        block_set_results_summed = collection.block_set_results_summed,
        _one_body_local_block_collection_summary_materialization_flags(
            source_operator_blocks_materialized =
                collection.source_operator_blocks_materialized,
            final_pair_blocks_materialized =
                collection.final_pair_blocks_materialized,
        )...,
    )
end

function _one_body_local_block_collection_entries_for_term(
    collection::NamedTuple,
    term::Symbol,
)
    _one_body_assert_local_block_collection(collection)
    return Tuple(entry for entry in collection.entries if entry.block_set_term === term)
end

function _one_body_local_block_collection_materialized_entries_for_term(
    collection::NamedTuple,
    term::Symbol,
)
    _one_body_assert_local_block_collection(collection)
    return Tuple(
        entry for entry in collection.materialized_entries
        if entry.block_set_term === term
    )
end

function _one_body_local_block_collection_skipped_entries_for_term(
    collection::NamedTuple,
    term::Symbol,
)
    _one_body_assert_local_block_collection(collection)
    return Tuple(
        entry for entry in collection.skipped_entries
        if entry.block_set_term === term
    )
end

function _one_body_local_block_collection_term_status(
    collection::NamedTuple,
    term::Symbol,
)
    _one_body_assert_local_block_collection(collection)
    materialized_entries =
        _one_body_local_block_collection_materialized_entries_for_term(
            collection,
            term,
        )
    skipped_entries =
        _one_body_local_block_collection_skipped_entries_for_term(
            collection,
            term,
        )
    entries = (materialized_entries..., skipped_entries...)
    status = _one_body_local_block_collection_term_status_symbol(
        collection,
        term,
        materialized_entries,
        skipped_entries,
    )
    return (;
        object_kind = :cartesian_pair_block_local_one_body_collection_term_status,
        term,
        status,
        blocker = _one_body_local_block_collection_term_blocker(
            status,
            skipped_entries,
        ),
        collection_status = collection.status,
        collection_blocker = collection.blocker,
        requested = term in collection.requested_terms,
        requested_materialization = term in collection.requested_materialize_terms,
        materialized_term = term in collection.materialized_terms,
        deferred_term = term in collection.deferred_terms,
        materialized_entry_count = length(materialized_entries),
        skipped_entry_count = length(skipped_entries),
        entry_count = length(entries),
        selector_family_counts =
            _one_body_count_optional_by(entries, :selector_family),
        materialized_selector_family_counts =
            _one_body_count_optional_by(materialized_entries, :selector_family),
        skipped_selector_family_counts =
            _one_body_count_optional_by(skipped_entries, :selector_family),
        skipped_blocker_counts =
            _one_body_count_optional_by(skipped_entries, :blocker),
        block_space_counts = _one_body_count_optional_by(entries, :block_space),
        term_separated_entries = true,
        pair_separated_entries = true,
        block_set_results_summed = false,
        _one_body_local_block_collection_summary_materialization_flags(
            source_operator_blocks_materialized = any(
                entry -> entry.source_operator_blocks_materialized,
                materialized_entries,
            ),
            final_pair_blocks_materialized = any(
                entry -> entry.final_pair_blocks_materialized,
                materialized_entries,
            ),
        )...,
    )
end

function _one_body_local_block_collection_entries_for_pair(
    collection::NamedTuple,
    pair_key::Tuple{Symbol,Symbol},
)
    _one_body_assert_local_block_collection(collection)
    return Tuple(entry for entry in collection.entries if entry.pair_key == pair_key)
end

function _one_body_local_block_collection_materialized_entries_for_pair(
    collection::NamedTuple,
    pair_key::Tuple{Symbol,Symbol},
)
    _one_body_assert_local_block_collection(collection)
    return Tuple(
        entry for entry in collection.materialized_entries
        if entry.pair_key == pair_key
    )
end

function _one_body_local_block_collection_skipped_entries_for_pair(
    collection::NamedTuple,
    pair_key::Tuple{Symbol,Symbol},
)
    _one_body_assert_local_block_collection(collection)
    return Tuple(
        entry for entry in collection.skipped_entries
        if entry.pair_key == pair_key
    )
end

function _one_body_local_block_collection_pair_status(
    collection::NamedTuple,
    pair_key::Tuple{Symbol,Symbol},
)
    _one_body_assert_local_block_collection(collection)
    materialized_entries =
        _one_body_local_block_collection_materialized_entries_for_pair(
            collection,
            pair_key,
        )
    skipped_entries =
        _one_body_local_block_collection_skipped_entries_for_pair(
            collection,
            pair_key,
        )
    entries = (materialized_entries..., skipped_entries...)
    status = _one_body_local_block_collection_pair_status_symbol(
        materialized_entries,
        skipped_entries,
    )
    return (;
        object_kind = :cartesian_pair_block_local_one_body_collection_pair_status,
        pair_key,
        status,
        blocker = _one_body_local_block_collection_pair_blocker(
            status,
            skipped_entries,
        ),
        collection_status = collection.status,
        collection_blocker = collection.blocker,
        block_set_terms = _one_body_collection_unique_values(
            entries,
            :block_set_term,
        ),
        result_terms = _one_body_collection_unique_values(entries, :result_term),
        source_space_terms =
            _one_body_collection_unique_values(entries, :source_space_term),
        materialized_entry_count = length(materialized_entries),
        skipped_entry_count = length(skipped_entries),
        entry_count = length(entries),
        selector_family_counts =
            _one_body_count_optional_by(entries, :selector_family),
        materialized_selector_family_counts =
            _one_body_count_optional_by(materialized_entries, :selector_family),
        skipped_selector_family_counts =
            _one_body_count_optional_by(skipped_entries, :selector_family),
        skipped_blocker_counts =
            _one_body_count_optional_by(skipped_entries, :blocker),
        block_space_counts = _one_body_count_optional_by(entries, :block_space),
        term_separated_entries = true,
        pair_separated_entries = true,
        block_set_results_summed = false,
        _one_body_local_block_collection_summary_materialization_flags(
            source_operator_blocks_materialized = any(
                entry -> entry.source_operator_blocks_materialized,
                materialized_entries,
            ),
            final_pair_blocks_materialized = any(
                entry -> entry.final_pair_blocks_materialized,
                materialized_entries,
            ),
        )...,
    )
end

function _one_body_local_block_collection_lookup(
    collection::NamedTuple,
    term::Symbol,
    pair_key::Tuple{Symbol,Symbol},
)
    _one_body_assert_local_block_collection(collection)
    materialized_entries = Tuple(
        entry for entry in
        _one_body_local_block_collection_materialized_entries_for_term(
            collection,
            term,
        ) if entry.pair_key == pair_key
    )
    skipped_entries = Tuple(
        entry for entry in _one_body_local_block_collection_skipped_entries_for_term(
            collection,
            term,
        ) if entry.pair_key == pair_key
    )
    entries = (materialized_entries..., skipped_entries...)
    status = _one_body_local_block_collection_lookup_status(
        collection,
        term,
        materialized_entries,
        skipped_entries,
    )
    base = (;
        object_kind = :cartesian_pair_block_local_one_body_collection_lookup,
        block_set_term = term,
        term,
        pair_key,
        status,
        blocker = _one_body_local_block_collection_lookup_blocker(
            status,
            skipped_entries,
        ),
        collection_status = collection.status,
        collection_blocker = collection.blocker,
        entry_available = !isempty(entries),
        materialized_entry_available = !isempty(materialized_entries),
        skipped_entry_available = !isempty(skipped_entries),
        materialized_entry_count = length(materialized_entries),
        skipped_entry_count = length(skipped_entries),
        entry_count = length(entries),
        selector_family =
            _one_body_local_block_collection_single_entry_value(
                entries,
                :selector_family,
            ),
        block_space =
            _one_body_local_block_collection_single_entry_value(
                entries,
                :block_space,
            ),
        result_term =
            _one_body_local_block_collection_single_entry_value(
                entries,
                :result_term,
            ),
        source_space_term =
            _one_body_local_block_collection_single_entry_value(
                entries,
                :source_space_term,
            ),
        block_set_terms = _one_body_collection_unique_values(
            entries,
            :block_set_term,
        ),
        result_terms = _one_body_collection_unique_values(entries, :result_term),
        source_space_terms =
            _one_body_collection_unique_values(entries, :source_space_term),
        term_requested = term in collection.requested_terms,
        requested_materialization = term in collection.requested_materialize_terms,
        term_materialized = term in collection.materialized_terms,
        term_deferred = term in collection.deferred_terms,
        term_separated_entries = true,
        pair_separated_entries = true,
        lookup_chose_between_multiple_entries = false,
        block_set_results_summed = false,
        _one_body_local_block_collection_summary_materialization_flags(
            source_operator_blocks_materialized = any(
                entry -> entry.source_operator_blocks_materialized,
                materialized_entries,
            ),
            final_pair_blocks_materialized = any(
                entry -> entry.final_pair_blocks_materialized,
                materialized_entries,
            ),
        )...,
    )
    length(entries) == 1 && return merge(base, (; entry = only(entries)))
    return base
end

function _one_body_collection_result_block_space(
    result::PairBlockMaterializationResult,
)
    result.final_pair_blocks_materialized && return :final_local_space
    result.source_operator_blocks_materialized && return :source_space
    return :not_materialized
end

function _one_body_local_block_collection_status(
    summary,
    materialized_entries::Tuple,
    skipped_entries::Tuple,
)
    isnothing(_one_body_local_block_collection_blocker(summary)) ||
        return :blocked_local_one_body_block_collection
    !isempty(materialized_entries) && !isempty(skipped_entries) &&
        return :partially_materialized_local_one_body_block_collection
    !isempty(materialized_entries) &&
        return :materialized_local_one_body_block_collection
    !isempty(skipped_entries) && return :skipped_local_one_body_block_collection
    !isempty(summary.deferred_terms) &&
        return :deferred_metadata_only_local_one_body_block_collection
    return :empty_local_one_body_block_collection
end

function _one_body_local_block_collection_blocker(summary)
    return _one_body_collection_value(summary, :blocker, nothing)
end

function _one_body_local_block_collection_term_status_symbol(
    collection::NamedTuple,
    term::Symbol,
    materialized_entries::Tuple,
    skipped_entries::Tuple,
)
    term in collection.requested_terms || return :term_not_requested
    term in collection.deferred_terms &&
        return :deferred_metadata_only_local_one_body_collection_term
    !isempty(materialized_entries) && !isempty(skipped_entries) &&
        return :partially_materialized_local_one_body_collection_term
    !isempty(materialized_entries) &&
        return :materialized_local_one_body_collection_term
    !isempty(skipped_entries) && return :skipped_local_one_body_collection_term
    return :empty_requested_local_one_body_collection_term
end

function _one_body_local_block_collection_term_blocker(
    status::Symbol,
    skipped_entries::Tuple,
)
    status === :term_not_requested && return :term_not_requested
    isempty(skipped_entries) && return nothing
    return _one_body_collection_value(first(skipped_entries), :blocker, nothing)
end

function _one_body_local_block_collection_pair_status_symbol(
    materialized_entries::Tuple,
    skipped_entries::Tuple,
)
    !isempty(materialized_entries) && !isempty(skipped_entries) &&
        return :partially_materialized_local_one_body_collection_pair
    !isempty(materialized_entries) &&
        return :materialized_local_one_body_collection_pair
    !isempty(skipped_entries) && return :skipped_local_one_body_collection_pair
    return :pair_key_not_found
end

function _one_body_local_block_collection_pair_blocker(
    status::Symbol,
    skipped_entries::Tuple,
)
    status === :pair_key_not_found && return :pair_key_not_found
    isempty(skipped_entries) && return nothing
    return _one_body_collection_value(first(skipped_entries), :blocker, nothing)
end

function _one_body_local_block_collection_lookup_status(
    collection::NamedTuple,
    term::Symbol,
    materialized_entries::Tuple,
    skipped_entries::Tuple,
)
    term in collection.requested_terms || return :term_not_requested
    term in collection.deferred_terms &&
        return :deferred_metadata_only_local_one_body_collection_lookup
    !isempty(materialized_entries) && !isempty(skipped_entries) &&
        return :partially_materialized_local_one_body_collection_lookup
    !isempty(materialized_entries) &&
        return :materialized_local_one_body_collection_lookup
    !isempty(skipped_entries) && return :skipped_local_one_body_collection_lookup
    return :pair_key_not_found
end

function _one_body_local_block_collection_lookup_blocker(
    status::Symbol,
    skipped_entries::Tuple,
)
    status === :term_not_requested && return :term_not_requested
    status === :pair_key_not_found && return :pair_key_not_found
    isempty(skipped_entries) && return nothing
    return _one_body_collection_value(first(skipped_entries), :blocker, nothing)
end

function _one_body_local_block_collection_single_entry_value(
    entries::Tuple,
    field::Symbol,
)
    length(entries) == 1 || return nothing
    return _one_body_collection_value(first(entries), field, nothing)
end

function _one_body_collection_unique_values(entries::Tuple, field::Symbol)
    values = Any[]
    for entry in entries
        haskey(entry, field) || continue
        value = getfield(entry, field)
        isnothing(value) && continue
        value in values && continue
        push!(values, value)
    end
    return Tuple(values)
end

function _one_body_assert_local_block_collection(collection::NamedTuple)
    _one_body_collection_value(collection, :object_kind, nothing) ===
    :cartesian_pair_block_local_one_body_block_collection || throw(
        ArgumentError(
            "local one-body block collection accessor requires a local block collection object",
        ),
    )
    return nothing
end

function _one_body_collection_metadata_value(metadata::NamedTuple, key::Symbol, default)
    return haskey(metadata, key) ? getfield(metadata, key) : default
end

function _one_body_collection_value(record::NamedTuple, key::Symbol, default)
    return haskey(record, key) ? getfield(record, key) : default
end

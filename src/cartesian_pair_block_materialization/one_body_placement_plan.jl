# Metadata-only placement plan for local one-body block collections.
#
# This layer converts local one-body collection entries into placement records
# for a future global retained-operator assembly step. It does not allocate a
# global matrix, copy block data, build Hamiltonians, materialize Coulomb or
# IDA/MWG data, export artifacts, or perform PQS shell projection/Lowdin
# realization.

function one_body_overlap_placement_plan(collection; global_dimension = nothing)
    return one_body_placement_plan(collection; term = :overlap, global_dimension)
end

function one_body_kinetic_placement_plan(collection; global_dimension = nothing)
    return one_body_placement_plan(collection; term = :kinetic, global_dimension)
end

function one_body_electron_nuclear_by_center_placement_plan(
    collection;
    global_dimension = nothing,
)
    return one_body_placement_plan(
        collection;
        term = :electron_nuclear_by_center,
        global_dimension,
    )
end

function one_body_position_x_placement_plan(collection; global_dimension = nothing)
    return one_body_placement_plan(collection; term = :position_x, global_dimension)
end

function one_body_position_y_placement_plan(collection; global_dimension = nothing)
    return one_body_placement_plan(collection; term = :position_y, global_dimension)
end

function one_body_position_z_placement_plan(collection; global_dimension = nothing)
    return one_body_placement_plan(collection; term = :position_z, global_dimension)
end

function one_body_x2_x_placement_plan(collection; global_dimension = nothing)
    return one_body_placement_plan(collection; term = :x2_x, global_dimension)
end

function one_body_x2_y_placement_plan(collection; global_dimension = nothing)
    return one_body_placement_plan(collection; term = :x2_y, global_dimension)
end

function one_body_x2_z_placement_plan(collection; global_dimension = nothing)
    return one_body_placement_plan(collection; term = :x2_z, global_dimension)
end

function one_body_placement_plan(
    collection::NamedTuple;
    term::Symbol = :overlap,
    global_dimension = nothing,
)
    term in _one_body_supported_global_placement_terms() || throw(
        ArgumentError(
            "local one-body placement planning currently supports :overlap, :kinetic, :position_x/:position_y/:position_z, :x2_x/:x2_y/:x2_z, and :electron_nuclear_by_center only",
        ),
    )
    _one_body_assert_local_block_collection(collection)
    global_dimension_status =
        _one_body_placement_global_dimension_status(global_dimension)
    placement_records = Tuple(
        _one_body_placement_record(entry; term, global_dimension)
        for entry in collection.entries
        if _one_body_placement_entry_matches_term(entry, term)
    )
    placeable_records = Tuple(
        record for record in placement_records
        if record.placeable_in_global_retained_operator
    )
    blocked_records = Tuple(
        record for record in placement_records
        if !record.placeable_in_global_retained_operator
    )
    blocker_counts = _one_body_placement_count_by(blocked_records, :blocker)

    return (;
        object_kind = :cartesian_pair_block_local_one_body_placement_plan,
        status = _one_body_placement_plan_status(
            placement_records,
            placeable_records,
            blocked_records,
        ),
        blocker = _one_body_placement_plan_blocker(
            placeable_records,
            blocked_records,
        ),
        term,
        placement_records,
        placeable_records,
        blocked_records,
        record_count = length(placement_records),
        placeable_count = length(placeable_records),
        blocked_count = length(blocked_records),
        blocker_counts,
        global_dimension,
        global_dimension_status,
        target_space = :global_retained_operator,
        global_matrix_assembly_status = :not_materialized,
        records_store_existing_collection_entry_references = true,
        records_copy_block_matrices = false,
        _one_body_placement_nonclaim_flags()...,
    )
end

function _one_body_supported_global_placement_terms()
    return (
        :overlap,
        :kinetic,
        :position_x,
        :position_y,
        :position_z,
        :x2_x,
        :x2_y,
        :x2_z,
        :electron_nuclear_by_center,
    )
end

function one_body_placement_plan(collection; kwargs...)
    throw(
        ArgumentError(
            "local one-body placement plan requires a local block collection NamedTuple",
        ),
    )
end

function _one_body_placement_record(
    entry::NamedTuple;
    term::Symbol,
    global_dimension,
)
    left_column_range = _one_body_placement_column_range(entry, :left)
    right_column_range = _one_body_placement_column_range(entry, :right)
    blocker = _one_body_placement_record_blocker(
        entry,
        left_column_range,
        right_column_range,
        global_dimension,
    )
    placeable = isnothing(blocker)
    target_ranges =
        placeable ?
        (; rows = left_column_range, columns = right_column_range) :
        nothing
    return (;
        object_kind = :cartesian_pair_block_local_one_body_placement_record,
        term,
        block_set_term = _one_body_placement_value(entry, :block_set_term, nothing),
        result_term = _one_body_placement_value(entry, :result_term, nothing),
        source_space_term =
            _one_body_placement_value(entry, :source_space_term, nothing),
        center_index = _one_body_placement_entry_metadata_value(
            entry,
            :center_index,
            nothing,
        ),
        center_key = _one_body_placement_entry_metadata_value(
            entry,
            :center_key,
            nothing,
        ),
        center_location = _one_body_placement_entry_metadata_value(
            entry,
            :center_location,
            nothing,
        ),
        nuclear_charge_recorded = _one_body_placement_entry_metadata_value(
            entry,
            :nuclear_charge_recorded,
            false,
        ),
        nuclear_charge = _one_body_placement_entry_metadata_value(
            entry,
            :nuclear_charge,
            nothing,
        ),
        nuclear_charge_applied = _one_body_placement_entry_metadata_value(
            entry,
            :nuclear_charge_applied,
            false,
        ),
        pair_key = _one_body_placement_value(entry, :pair_key, nothing),
        pair_index = _one_body_placement_value(entry, :pair_index, nothing),
        selector_family =
            _one_body_placement_value(entry, :selector_family, nothing),
        materialization_path =
            _one_body_placement_value(entry, :materialization_path, nothing),
        entry_kind = _one_body_placement_value(entry, :entry_kind, nothing),
        block_space = _one_body_placement_value(entry, :block_space, nothing),
        block_shape = _one_body_placement_value(entry, :block_shape, nothing),
        source_entry_status = _one_body_placement_value(entry, :status, nothing),
        status =
            placeable ?
            :placeable_local_one_body_placement_record :
            :blocked_local_one_body_placement_record,
        blocker,
        not_placeable_reason = blocker,
        left_column_range,
        right_column_range,
        target_ranges,
        target_slice = target_ranges,
        global_dimension,
        placeable_in_global_retained_operator = placeable,
        source_collection_entry = entry,
        result = _one_body_placement_value(entry, :result, nothing),
        skipped_record = _one_body_placement_value(entry, :skipped_record, nothing),
        record_copies_block_matrix = false,
        _one_body_placement_nonclaim_flags()...,
    )
end

function _one_body_placement_entry_matches_term(entry::NamedTuple, term::Symbol)
    return _one_body_placement_value(entry, :block_set_term, nothing) === term
end

function _one_body_placement_record_blocker(
    entry::NamedTuple,
    left_column_range,
    right_column_range,
    global_dimension,
)
    block_space = _one_body_placement_value(entry, :block_space, nothing)
    block_space === :source_space &&
        return :source_space_block_requires_shell_realization
    block_space === :not_materialized &&
        return _one_body_placement_value(entry, :blocker, :not_materialized)
    block_space === :final_local_space ||
        return :non_final_local_block_not_placeable

    (
        _one_body_placement_valid_column_range(left_column_range) &&
        _one_body_placement_valid_column_range(right_column_range)
    ) || return :missing_column_ranges

    block_shape = _one_body_placement_value(entry, :block_shape, nothing)
    _one_body_placement_column_ranges_match_shape(
        left_column_range,
        right_column_range,
        block_shape,
    ) || return :column_ranges_do_not_match_block_shape

    _one_body_placement_ranges_inside_global_dimension(
        left_column_range,
        right_column_range,
        global_dimension,
    ) || return :column_ranges_outside_global_dimension

    return nothing
end

function _one_body_placement_column_range(entry::NamedTuple, side::Symbol)
    side in (:left, :right) ||
        throw(ArgumentError("placement column range side must be :left or :right"))
    keys =
        side === :left ?
        (:left_column_range, :left_final_column_range) :
        (:right_column_range, :right_final_column_range)
    for key in keys
        value = _one_body_placement_value(entry, key, nothing)
        isnothing(value) || return value
    end

    metadata = _one_body_placement_entry_metadata(entry)
    for key in keys
        value = _one_body_placement_value(metadata, key, nothing)
        isnothing(value) || return value
    end
    return nothing
end

function _one_body_placement_entry_metadata(entry::NamedTuple)
    result = _one_body_placement_value(entry, :result, nothing)
    result isa PairBlockMaterializationResult && return result.metadata

    skipped_record = _one_body_placement_value(entry, :skipped_record, nothing)
    skipped_record isa NamedTuple && return skipped_record

    return (;)
end

function _one_body_placement_entry_metadata_value(
    entry::NamedTuple,
    key::Symbol,
    default = nothing,
)
    value = _one_body_placement_value(entry, key, nothing)
    isnothing(value) || return value
    metadata = _one_body_placement_entry_metadata(entry)
    return _one_body_placement_value(metadata, key, default)
end

function _one_body_placement_valid_column_range(value)
    isnothing(value) && return false
    value isa AbstractUnitRange{<:Integer} && return true
    return false
end

function _one_body_placement_column_ranges_match_shape(
    left_column_range,
    right_column_range,
    block_shape,
)
    block_shape isa Tuple || return false
    length(block_shape) == 2 || return false
    return length(left_column_range) == block_shape[1] &&
           length(right_column_range) == block_shape[2]
end

function _one_body_placement_ranges_inside_global_dimension(
    left_column_range,
    right_column_range,
    global_dimension,
)
    isnothing(global_dimension) && return true
    dimension = _one_body_placement_global_dimension(global_dimension)
    return first(left_column_range) >= 1 &&
           first(right_column_range) >= 1 &&
           last(left_column_range) <= dimension &&
           last(right_column_range) <= dimension
end

function _one_body_placement_global_dimension_status(global_dimension)
    isnothing(global_dimension) && return :not_supplied
    _one_body_placement_global_dimension(global_dimension)
    return :available
end

function _one_body_placement_global_dimension(global_dimension)
    global_dimension isa Integer ||
        throw(ArgumentError("global_dimension must be an integer or nothing"))
    dimension = Int(global_dimension)
    dimension > 0 ||
        throw(ArgumentError("global_dimension must be positive"))
    return dimension
end

function _one_body_placement_plan_status(
    placement_records::Tuple,
    placeable_records::Tuple,
    blocked_records::Tuple,
)
    isempty(placement_records) && return :empty_local_one_body_placement_plan
    !isempty(placeable_records) && isempty(blocked_records) &&
        return :placeable_local_one_body_placement_plan
    !isempty(placeable_records) && !isempty(blocked_records) &&
        return :partially_placeable_local_one_body_placement_plan
    return :blocked_local_one_body_placement_plan
end

function _one_body_placement_plan_blocker(
    placeable_records::Tuple,
    blocked_records::Tuple,
)
    isempty(blocked_records) && return nothing
    !isempty(placeable_records) && return :blocked_placement_records_present
    blockers = Tuple(record.blocker for record in blocked_records)
    length(unique(blockers)) == 1 && return only(blockers)
    return :multiple_blockers
end

function _one_body_placement_count_by(records::Tuple, field::Symbol)
    values = Symbol[]
    for record in records
        value = _one_body_placement_value(record, field, nothing)
        value isa Symbol || continue
        push!(values, value)
    end
    ordered_values = Symbol[]
    for value in values
        value in ordered_values && continue
        push!(ordered_values, value)
    end
    return Tuple(
        (; field => value, count = count(==(value), values))
        for value in ordered_values
    )
end

function _one_body_placement_value(source::NamedTuple, key::Symbol, default = nothing)
    return haskey(source, key) ? getfield(source, key) : default
end

function _one_body_placement_value(source, _key::Symbol, default = nothing)
    return default
end

function _one_body_placement_nonclaim_flags()
    return (;
        operator_matrix_materialized = false,
        operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
        exports_materialized = false,
        global_operator_assembled = false,
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

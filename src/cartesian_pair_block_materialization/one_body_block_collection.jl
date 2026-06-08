# Private local one-body block collection entry vocabulary.
#
# These entries are a local view layer over existing mixed one-body results and
# skip summaries. They do not assemble terms, place global blocks, or copy
# matrix data; materialized entries keep only compact block-shape fields beside
# the original result reference.

function _one_body_local_block_collection_entry(
    result::PairBlockMaterializationResult,
)
    return (;
        object_kind = :cartesian_pair_block_local_one_body_block_collection_entry,
        entry_kind = :materialized_result,
        term = result.term,
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
        block_space = _one_body_collection_result_block_space(result),
        block_shape = size(result.block),
        block_size = length(result.block),
        status = :materialized_local_one_body_block_collection_entry,
        blocker = nothing,
        materialized = result.materialized,
        result_available = true,
        skipped_record_available = false,
        result,
        skipped_record = nothing,
        source_operator_blocks_materialized =
            result.source_operator_blocks_materialized,
        final_pair_blocks_materialized = result.final_pair_blocks_materialized,
        operator_blocks_materialized = result.operator_blocks_materialized,
        hamiltonian_data_materialized = result.hamiltonian_data_materialized,
        artifacts_materialized = result.artifacts_materialized,
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

function _one_body_local_block_collection_entry(result)
    throw(
        ArgumentError(
            "local one-body block collection result entry requires a PairBlockMaterializationResult",
        ),
    )
end

function _one_body_local_block_collection_skipped_entry(skip::NamedTuple)
    return (;
        object_kind = :cartesian_pair_block_local_one_body_block_collection_entry,
        entry_kind = :skipped_record,
        term = _one_body_collection_value(skip, :requested_term, nothing),
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
        source_operator_blocks_materialized = false,
        final_pair_blocks_materialized = false,
        operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
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

function _one_body_local_block_collection_skipped_entry(skip)
    throw(
        ArgumentError(
            "local one-body block collection skipped entry requires a NamedTuple",
        ),
    )
end

function _one_body_collection_result_block_space(
    result::PairBlockMaterializationResult,
)
    result.final_pair_blocks_materialized && return :final_local_space
    result.source_operator_blocks_materialized && return :source_space
    return :not_materialized
end

function _one_body_collection_metadata_value(metadata::NamedTuple, key::Symbol, default)
    return haskey(metadata, key) ? getfield(metadata, key) : default
end

function _one_body_collection_value(record::NamedTuple, key::Symbol, default)
    return haskey(record, key) ? getfield(record, key) : default
end

# PQS source-space block to shell-realization bridge metadata.

"""
    pqs_source_pair_shell_realization_bridge_summary(result)

Return a compact metadata-only bridge summary for a materialized PQS/PQS raw
source-space block. The summary records the source block facts and the planned
shell-realization path without building shell projection, Lowdin objects, final
retained pair blocks, Hamiltonian data, or artifacts.
"""
function pqs_source_pair_shell_realization_bridge_summary(
    result::PairBlockMaterializationResult,
)
    status, blocker = _pqs_source_shell_realization_bridge_status(result)
    metadata = result.metadata
    return (;
        object_kind = :pqs_source_pair_shell_realization_bridge_summary,
        status,
        blocker,
        pair_key = result.pair_key,
        source_block_term = result.term,
        source_block_status = _pqs_source_block_status(result),
        source_block_materialized = result.materialized,
        block_space = _pair_block_metadata_value(metadata, :block_space),
        left_source_mode_dims =
            _pair_block_metadata_value(metadata, :left_source_mode_dims),
        right_source_mode_dims =
            _pair_block_metadata_value(metadata, :right_source_mode_dims),
        left_source_mode_count =
            _pair_block_metadata_value(metadata, :left_source_mode_count),
        right_source_mode_count =
            _pair_block_metadata_value(metadata, :right_source_mode_count),
        source_mode_ordering =
            _pair_block_metadata_value(metadata, :source_mode_ordering),
        transform_contract_keys =
            _pair_block_metadata_value(metadata, :transform_contract_keys),
        source_contract_keys =
            _pair_block_metadata_value(metadata, :source_contract_keys),
        transform_paths = _pair_block_metadata_value(metadata, :transform_paths),
        realization_paths = _pair_block_metadata_value(metadata, :realization_paths),
        source_operator_blocks_materialized =
            result.source_operator_blocks_materialized,
        final_pair_blocks_materialized = false,
        shell_realization_materialized = false,
        operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
    )
end

function _pqs_source_shell_realization_bridge_status(
    result::PairBlockMaterializationResult,
)
    metadata = result.metadata
    _pair_block_metadata_value(metadata, :block_space) === :raw_product_source_modes ||
        return (
            :blocked_not_pqs_source_space_block,
            :not_raw_product_source_modes,
        )
    result.materialized && result.source_operator_blocks_materialized ||
        return (
            :blocked_source_block_not_materialized,
            :source_operator_block_not_materialized,
        )

    for field in (
        :left_source_mode_dims,
        :right_source_mode_dims,
        :left_source_mode_count,
        :right_source_mode_count,
        :source_mode_ordering,
    )
        isnothing(_pair_block_metadata_value(metadata, field)) &&
            return (:blocked_missing_source_block_metadata, Symbol(:missing_, field))
    end

    transform_contract_keys =
        _pair_block_metadata_value(metadata, :transform_contract_keys)
    _complete_left_right_metadata(transform_contract_keys) ||
        return (
            :blocked_missing_shell_realization_facts,
            :missing_transform_contract_keys,
        )

    realization_paths = _pair_block_metadata_value(metadata, :realization_paths)
    _complete_left_right_metadata(realization_paths) ||
        return (
            :blocked_missing_shell_realization_facts,
            :missing_shell_realization_path,
        )

    return :available_metadata_only_shell_realization_bridge, nothing
end

function _pqs_source_block_status(result::PairBlockMaterializationResult)
    return result.materialized && result.source_operator_blocks_materialized ?
           :source_operator_block_materialized :
           :source_operator_block_not_materialized
end

function _pair_block_metadata_value(metadata::NamedTuple, key::Symbol, default = nothing)
    return haskey(metadata, key) ? getfield(metadata, key) : default
end

function _complete_left_right_metadata(value)
    value isa NamedTuple || return false
    haskey(value, :left) && haskey(value, :right) || return false
    isnothing(value.left) && return false
    isnothing(value.right) && return false
    return true
end

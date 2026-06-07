# PQS source bridge to final retained pair-block readiness metadata.

"""
    pqs_source_pair_final_block_readiness_summary(bridge_summary)

Return a compact metadata-only readiness summary for attempting a future final
PQS retained pair block from a PQS source shell-realization bridge summary.
This does not build shell projection, Lowdin objects, final pair blocks,
Hamiltonian data, exports, artifacts, IDA/MWG data, or Coulomb blocks.
"""
function pqs_source_pair_final_block_readiness_summary(bridge_summary::NamedTuple)
    object_kind = _bridge_summary_value(bridge_summary, :object_kind)
    if object_kind === :pqs_source_pair_shell_realization_bridge_summary
        return _pqs_source_pair_final_block_readiness_single(bridge_summary)
    elseif object_kind === :pqs_source_pair_shell_realization_bridge_batch_summary
        return _pqs_source_pair_final_block_readiness_batch(bridge_summary)
    end
    return (;
        object_kind = :pqs_source_pair_final_block_readiness_summary,
        status = :blocked_final_pqs_pair_block_not_ready,
        blocker = :not_pqs_source_shell_realization_bridge_summary,
        bridge_object_kind = object_kind,
        source_operator_blocks_materialized = false,
        shell_realization_materialized = false,
        final_pair_blocks_materialized = false,
        operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
    )
end

function _pqs_source_pair_final_block_readiness_single(bridge_summary::NamedTuple)
    status, blocker = _pqs_source_pair_final_block_readiness_status(bridge_summary)
    return (;
        object_kind = :pqs_source_pair_final_block_readiness_summary,
        status,
        blocker,
        bridge_object_kind = bridge_summary.object_kind,
        bridge_status = _bridge_summary_value(bridge_summary, :status),
        bridge_blocker = _bridge_summary_value(bridge_summary, :blocker),
        pair_key = _bridge_summary_value(bridge_summary, :pair_key),
        source_block_term = _bridge_summary_value(bridge_summary, :source_block_term),
        source_block_status =
            _bridge_summary_value(bridge_summary, :source_block_status),
        block_space = _bridge_summary_value(bridge_summary, :block_space),
        left_source_mode_dims =
            _bridge_summary_value(bridge_summary, :left_source_mode_dims),
        right_source_mode_dims =
            _bridge_summary_value(bridge_summary, :right_source_mode_dims),
        left_source_mode_count =
            _bridge_summary_value(bridge_summary, :left_source_mode_count),
        right_source_mode_count =
            _bridge_summary_value(bridge_summary, :right_source_mode_count),
        source_mode_ordering =
            _bridge_summary_value(bridge_summary, :source_mode_ordering),
        transform_contract_keys =
            _bridge_summary_value(bridge_summary, :transform_contract_keys),
        source_contract_keys =
            _bridge_summary_value(bridge_summary, :source_contract_keys),
        realization_paths = _bridge_summary_value(bridge_summary, :realization_paths),
        source_operator_blocks_materialized =
            _bridge_summary_value(
                bridge_summary,
                :source_operator_blocks_materialized,
                false,
            ),
        shell_realization_materialized =
            _bridge_summary_value(
                bridge_summary,
                :shell_realization_materialized,
                false,
            ),
        final_pair_blocks_materialized = false,
        operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
    )
end

function _pqs_source_pair_final_block_readiness_batch(bridge_summary::NamedTuple)
    status, blocker = _pqs_source_pair_final_block_readiness_status(bridge_summary)
    return (;
        object_kind = :pqs_source_pair_final_block_readiness_batch_summary,
        status,
        blocker,
        bridge_object_kind = bridge_summary.object_kind,
        bridge_status = _bridge_summary_value(bridge_summary, :status),
        bridge_blocker = _bridge_summary_value(bridge_summary, :blocker),
        term = _bridge_summary_value(bridge_summary, :term),
        term_counts = _bridge_summary_value(bridge_summary, :term_counts, ()),
        result_count = _bridge_summary_value(bridge_summary, :result_count, 0),
        available_count = _bridge_summary_value(bridge_summary, :available_count, 0),
        blocked_count = _bridge_summary_value(bridge_summary, :blocked_count, 0),
        bridge_status_counts =
            _bridge_summary_value(bridge_summary, :bridge_status_counts, ()),
        blocker_counts = _bridge_summary_value(bridge_summary, :blocker_counts, ()),
        skipped_record_count =
            _bridge_summary_value(bridge_summary, :skipped_record_count, 0),
        skipped_blocker_counts =
            _bridge_summary_value(bridge_summary, :skipped_blocker_counts, ()),
        source_mode_ordering_status =
            _bridge_summary_value(bridge_summary, :source_mode_ordering_status),
        source_mode_ordering =
            _bridge_summary_value(bridge_summary, :source_mode_ordering),
        source_mode_ordering_counts =
            _bridge_summary_value(
                bridge_summary,
                :source_mode_ordering_counts,
                (),
            ),
        source_operator_blocks_materialized =
            _bridge_summary_value(
                bridge_summary,
                :source_operator_blocks_materialized,
                false,
            ),
        shell_realization_materialized =
            _bridge_summary_value(
                bridge_summary,
                :shell_realization_materialized,
                false,
            ),
        final_pair_blocks_materialized = false,
        operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
    )
end

function _pqs_source_pair_final_block_readiness_status(bridge_summary::NamedTuple)
    bridge_status = _bridge_summary_value(bridge_summary, :status)
    if bridge_status === :available_metadata_only_shell_realization_bridge ||
       bridge_status === :available_metadata_only_shell_realization_bridge_batch
        return _pqs_source_pair_final_block_ready_or_blocked(bridge_summary)
    end
    bridge_blocker = _bridge_summary_value(
        bridge_summary,
        :blocker,
        :source_shell_realization_bridge_not_available,
    )
    return :blocked_final_pqs_pair_block_not_ready, bridge_blocker
end

function _pqs_source_pair_final_block_ready_or_blocked(bridge_summary::NamedTuple)
    _bridge_summary_value(bridge_summary, :source_operator_blocks_materialized, false) ||
        return (
            :blocked_final_pqs_pair_block_not_ready,
            :source_operator_block_not_materialized,
        )

    if haskey(bridge_summary, :result_count)
        _bridge_summary_value(bridge_summary, :result_count, 0) > 0 ||
            return (
                :blocked_final_pqs_pair_block_not_ready,
                :no_pqs_source_block_results,
            )
    else
        for field in (
            :left_source_mode_dims,
            :right_source_mode_dims,
            :left_source_mode_count,
            :right_source_mode_count,
            :source_mode_ordering,
        )
            isnothing(_bridge_summary_value(bridge_summary, field)) &&
                return (
                    :blocked_final_pqs_pair_block_not_ready,
                    Symbol(:missing_, field),
                )
        end

        _complete_left_right_metadata(
            _bridge_summary_value(bridge_summary, :transform_contract_keys),
        ) || return (
            :blocked_final_pqs_pair_block_not_ready,
            :missing_transform_contract_keys,
        )
        _complete_left_right_metadata(
            _bridge_summary_value(bridge_summary, :realization_paths),
        ) || return (
            :blocked_final_pqs_pair_block_not_ready,
            :missing_shell_realization_path,
        )
    end

    _bridge_summary_value(bridge_summary, :shell_realization_materialized, false) ||
        return (
            :blocked_final_pqs_pair_block_not_ready,
            :shell_realization_not_materialized,
        )

    return :ready_metadata_only_final_pqs_pair_block, nothing
end

function _bridge_summary_value(
    bridge_summary::NamedTuple,
    key::Symbol,
    default = nothing,
)
    return haskey(bridge_summary, key) ? getfield(bridge_summary, key) : default
end

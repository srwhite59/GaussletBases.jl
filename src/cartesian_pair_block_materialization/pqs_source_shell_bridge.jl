# PQS source-space block to shell-realization bridge metadata.

"""
    pqs_source_pair_shell_realization_bridge_summary(result)
    pqs_source_pair_shell_realization_bridge_summary(batch_result)
    pqs_source_pair_shell_realization_bridge_summary(results)

Return compact metadata-only bridge summaries for materialized PQS/PQS raw
source-space blocks. The summaries record source block facts and the planned
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

function pqs_source_pair_shell_realization_bridge_summary(
    batch_result::PairBlockMaterializationBatchResult,
)
    return _pqs_source_pair_shell_realization_bridge_batch_summary(
        batch_result.materialized_results;
        term = batch_result.term,
        skipped_records = batch_result.skipped_records,
    )
end

function pqs_source_pair_shell_realization_bridge_summary(
    results::Tuple{Vararg{PairBlockMaterializationResult}},
)
    return _pqs_source_pair_shell_realization_bridge_batch_summary(results)
end

function pqs_source_pair_shell_realization_bridge_summary(
    results::AbstractVector{<:PairBlockMaterializationResult},
)
    return _pqs_source_pair_shell_realization_bridge_batch_summary(Tuple(results))
end

function _pqs_source_pair_shell_realization_bridge_batch_summary(
    results::Tuple{Vararg{PairBlockMaterializationResult}};
    term = _bridge_summary_term(results),
    skipped_records = (),
)
    bridge_summaries = Tuple(
        pqs_source_pair_shell_realization_bridge_summary(result)
        for result in results
    )
    available_count =
        count(
            summary ->
                summary.status === :available_metadata_only_shell_realization_bridge,
            bridge_summaries,
        )
    blocked_count = length(bridge_summaries) - available_count
    status, blocker =
        _pqs_source_pair_shell_realization_bridge_batch_status(
            length(bridge_summaries),
            blocked_count,
            bridge_summaries,
        )
    ordering_status, ordering = _bridge_source_mode_ordering_status(bridge_summaries)

    return (;
        object_kind = :pqs_source_pair_shell_realization_bridge_batch_summary,
        status,
        blocker,
        term,
        term_counts = _bridge_counts(
            Tuple(summary.source_block_term for summary in bridge_summaries),
            :term,
        ),
        result_count = length(bridge_summaries),
        available_count,
        blocked_count,
        bridge_status_counts = _bridge_counts(
            Tuple(summary.status for summary in bridge_summaries),
            :status,
        ),
        blocker_counts = _bridge_blocker_counts(bridge_summaries),
        skipped_record_count = length(skipped_records),
        skipped_blocker_counts = _skipped_record_blocker_counts(skipped_records),
        source_mode_ordering_status = ordering_status,
        source_mode_ordering = ordering,
        source_mode_ordering_counts =
            _source_mode_ordering_counts(bridge_summaries),
        source_operator_blocks_materialized =
            any(
                summary -> summary.source_operator_blocks_materialized,
                bridge_summaries,
            ),
        final_pair_blocks_materialized = false,
        shell_realization_materialized = false,
        operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
        bridge_summaries,
    )
end

function _pqs_source_pair_shell_realization_bridge_batch_status(
    result_count::Int,
    blocked_count::Int,
    bridge_summaries,
)
    result_count == 0 && return (
        :blocked_no_pqs_source_block_results,
        :no_pqs_source_block_results,
    )
    blocked_count == 0 &&
        return :available_metadata_only_shell_realization_bridge_batch, nothing
    return (
        :blocked_pqs_source_shell_realization_bridge_batch,
        first(
            summary.blocker for summary in bridge_summaries
            if !isnothing(summary.blocker)
        ),
    )
end

function _pqs_source_shell_realization_bridge_status(
    result::PairBlockMaterializationResult,
)
    metadata = result.metadata
    block_space = _pair_block_metadata_value(metadata, :block_space)
    isnothing(block_space) &&
        return (:blocked_missing_source_block_metadata, :missing_block_space)
    block_space === :raw_product_source_modes ||
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

function _bridge_summary_term(results::Tuple{Vararg{PairBlockMaterializationResult}})
    isempty(results) && return nothing
    terms = Tuple(result.term for result in results)
    all(term -> term === first(terms), terms) && return first(terms)
    return :mixed_source_terms
end

function _bridge_source_mode_ordering_status(bridge_summaries)
    orderings = Tuple(
        summary.source_mode_ordering for summary in bridge_summaries
        if !isnothing(summary.source_mode_ordering)
    )
    isempty(orderings) && return :unavailable_source_mode_ordering, nothing

    unique_orderings = unique(orderings)
    if length(unique_orderings) == 1 && length(orderings) == length(bridge_summaries)
        return :uniform_source_mode_ordering, only(unique_orderings)
    end
    length(unique_orderings) == 1 &&
        return :partially_unavailable_source_mode_ordering, only(unique_orderings)
    return :mixed_source_mode_ordering, nothing
end

function _source_mode_ordering_counts(bridge_summaries)
    orderings = Tuple(
        summary.source_mode_ordering for summary in bridge_summaries
        if !isnothing(summary.source_mode_ordering)
    )
    return _bridge_counts(orderings, :source_mode_ordering)
end

function _bridge_blocker_counts(bridge_summaries)
    blockers = Tuple(
        summary.blocker for summary in bridge_summaries
        if !isnothing(summary.blocker)
    )
    return _bridge_counts(blockers, :blocker)
end

function _skipped_record_blocker_counts(skipped_records)
    blockers = Tuple(
        skipped.blocker for skipped in skipped_records
        if haskey(skipped, :blocker) && !isnothing(skipped.blocker)
    )
    return _bridge_counts(blockers, :blocker)
end

function _pqs_source_block_status(result::PairBlockMaterializationResult)
    return result.materialized && result.source_operator_blocks_materialized ?
           :source_operator_block_materialized :
           :source_operator_block_not_materialized
end

function _pair_block_metadata_value(metadata::NamedTuple, key::Symbol, default = nothing)
    return haskey(metadata, key) ? getfield(metadata, key) : default
end

function _bridge_counts(values, field::Symbol)
    counts = Dict{Any,Int}()
    order = Any[]
    for value in values
        if !haskey(counts, value)
            push!(order, value)
            counts[value] = 0
        end
        counts[value] += 1
    end
    return Tuple((; field => value, count = counts[value]) for value in order)
end

function _complete_left_right_metadata(value)
    value isa NamedTuple || return false
    haskey(value, :left) && haskey(value, :right) || return false
    isnothing(value.left) && return false
    isnothing(value.right) && return false
    return true
end

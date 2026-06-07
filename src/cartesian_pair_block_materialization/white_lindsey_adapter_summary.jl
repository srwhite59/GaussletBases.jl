# Metadata-only old-kernel reuse summary for White--Lindsey boundary strata.

"""
    white_lindsey_boundary_stratum_adapter_summary(record)

Return compact metadata describing which existing White--Lindsey numerical
kernels a future boundary-stratum adapter should reuse. This does not call the
old kernels and does not build LW numerical blocks, coefficient maps, doside
transforms, Hamiltonian data, exports, artifacts, IDA/MWG data, or Coulomb.
"""
function white_lindsey_boundary_stratum_adapter_summary(
    record::PairBlockMaterializationRecord,
)
    record.materialization_path ===
        :white_lindsey_boundary_stratum_adapter_preflight ||
        return _white_lindsey_boundary_stratum_blocked_summary(
            record,
            :blocked_white_lindsey_boundary_stratum_adapter_summary,
            :not_white_lindsey_boundary_stratum_adapter_preflight,
        )

    left_stratum_kind = _white_lindsey_record_metadata_value(
        record,
        :left_stratum_kind,
    )
    right_stratum_kind = _white_lindsey_record_metadata_value(
        record,
        :right_stratum_kind,
    )
    left_plan = _white_lindsey_stratum_kernel_plan(left_stratum_kind)
    right_plan = _white_lindsey_stratum_kernel_plan(right_stratum_kind)

    if left_plan.status !== :available_white_lindsey_stratum_kernel_plan
        return _white_lindsey_boundary_stratum_blocked_summary(
            record,
            :blocked_white_lindsey_boundary_stratum_adapter_summary,
            left_plan.blocker,
            left_stratum_kind,
            right_stratum_kind,
            left_plan,
            right_plan,
        )
    elseif right_plan.status !== :available_white_lindsey_stratum_kernel_plan
        return _white_lindsey_boundary_stratum_blocked_summary(
            record,
            :blocked_white_lindsey_boundary_stratum_adapter_summary,
            right_plan.blocker,
            left_stratum_kind,
            right_stratum_kind,
            left_plan,
            right_plan,
        )
    end

    return _white_lindsey_boundary_stratum_adapter_summary(
        record,
        record.readiness_status,
        record.blocker,
        left_stratum_kind,
        right_stratum_kind,
        left_plan,
        right_plan,
    )
end

function white_lindsey_boundary_stratum_adapter_summary(
    records::Tuple{Vararg{PairBlockMaterializationRecord}},
)
    return _white_lindsey_boundary_stratum_adapter_batch_summary(records)
end

function white_lindsey_boundary_stratum_adapter_summary(
    records::AbstractVector{<:PairBlockMaterializationRecord},
)
    return _white_lindsey_boundary_stratum_adapter_batch_summary(Tuple(records))
end

function white_lindsey_boundary_stratum_adapter_summary(
    plan::PairBlockMaterializationPlan,
)
    return _white_lindsey_boundary_stratum_adapter_batch_summary(
        pair_block_materialization_records(plan),
    )
end

function _white_lindsey_boundary_stratum_adapter_summary(
    record::PairBlockMaterializationRecord,
    status::Symbol,
    blocker,
    left_stratum_kind,
    right_stratum_kind,
    left_plan,
    right_plan,
)
    return (;
        object_kind = :white_lindsey_boundary_stratum_adapter_summary,
        status,
        blocker,
        pair_key = record.pair_key,
        pair_index = record.pair_index,
        materialization_path = record.materialization_path,
        readiness_status = record.readiness_status,
        left_stratum_kind,
        right_stratum_kind,
        left_planned_kernel = left_plan.kernel,
        right_planned_kernel = right_plan.kernel,
        left_side_1d_helper = left_plan.side_1d_helper,
        right_side_1d_helper = right_plan.side_1d_helper,
        shared_1d_helper =
            _white_lindsey_shared_1d_helper(left_plan, right_plan),
        legacy_kernel_reuse_map = (;
            facet = :_nested_face_product,
            face = :_nested_face_product,
            edge = :_nested_edge_product,
            corner = :_nested_corner_piece,
            side_1d = :_nested_doside_1d,
        ),
        source_operator_blocks_materialized = false,
        final_pair_blocks_materialized = false,
        operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
    )
end

function _white_lindsey_boundary_stratum_adapter_batch_summary(
    records::Tuple{Vararg{PairBlockMaterializationRecord}},
)
    adapter_records = Tuple(
        record for record in records
        if record.materialization_path ===
           :white_lindsey_boundary_stratum_adapter_preflight
    )
    skipped_count = length(records) - length(adapter_records)
    adapter_summaries = Tuple(
        white_lindsey_boundary_stratum_adapter_summary(record)
        for record in adapter_records
    )
    available_count = count(
        _white_lindsey_adapter_reuse_metadata_available,
        adapter_summaries,
    )
    blocked_count = length(adapter_summaries) - available_count
    status, blocker =
        _white_lindsey_boundary_stratum_adapter_batch_status(
            length(adapter_summaries),
            blocked_count,
            adapter_summaries,
        )

    left_kernels = Tuple(
        summary.left_planned_kernel for summary in adapter_summaries
        if !isnothing(summary.left_planned_kernel)
    )
    right_kernels = Tuple(
        summary.right_planned_kernel for summary in adapter_summaries
        if !isnothing(summary.right_planned_kernel)
    )
    strata = Tuple(
        stratum for summary in adapter_summaries
        for stratum in (summary.left_stratum_kind, summary.right_stratum_kind)
        if !isnothing(stratum)
    )
    shared_helpers = Tuple(
        summary.shared_1d_helper for summary in adapter_summaries
        if !isnothing(summary.shared_1d_helper)
    )

    return (;
        object_kind = :white_lindsey_boundary_stratum_adapter_batch_summary,
        status,
        blocker,
        input_record_count = length(records),
        record_count = length(adapter_summaries),
        available_count,
        blocked_count,
        skipped_record_count = skipped_count,
        skipped_blocker_counts =
            skipped_count == 0 ?
            () :
            ((;
                blocker = :not_white_lindsey_boundary_stratum_adapter_preflight,
                count = skipped_count,
            ),),
        status_counts = _white_lindsey_count_by_value(
            Tuple(summary.status for summary in adapter_summaries),
            :status,
        ),
        blocker_counts = _white_lindsey_count_by_value(
            Tuple(
                summary.blocker for summary in adapter_summaries
                if !isnothing(summary.blocker)
            ),
            :blocker,
        ),
        left_planned_kernel_counts =
            _white_lindsey_count_by_value(left_kernels, :planned_kernel),
        right_planned_kernel_counts =
            _white_lindsey_count_by_value(right_kernels, :planned_kernel),
        all_planned_kernel_counts =
            _white_lindsey_count_by_value(
                (left_kernels..., right_kernels...),
                :planned_kernel,
            ),
        stratum_kind_counts =
            _white_lindsey_count_by_value(strata, :stratum_kind),
        shared_1d_helper_counts =
            _white_lindsey_count_by_value(shared_helpers, :helper),
        source_operator_blocks_materialized = false,
        final_pair_blocks_materialized = false,
        operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
    )
end

function _white_lindsey_boundary_stratum_adapter_batch_status(
    record_count::Int,
    blocked_count::Int,
    adapter_summaries,
)
    record_count == 0 && return (
        :blocked_white_lindsey_boundary_stratum_adapter_batch_summary,
        :no_white_lindsey_boundary_stratum_adapter_records,
    )
    blocked_count == 0 &&
        return :available_metadata_only_white_lindsey_adapter_reuse_batch, nothing
    return (
        :blocked_white_lindsey_boundary_stratum_adapter_batch_summary,
        first(
            summary.blocker for summary in adapter_summaries
            if !_white_lindsey_adapter_reuse_metadata_available(summary) &&
               !isnothing(summary.blocker)
        ),
    )
end

function _white_lindsey_adapter_reuse_metadata_available(summary)
    return !isnothing(summary.left_planned_kernel) &&
           !isnothing(summary.right_planned_kernel)
end

function _white_lindsey_boundary_stratum_blocked_summary(
    record::PairBlockMaterializationRecord,
    status::Symbol,
    blocker,
    left_stratum_kind = _white_lindsey_record_metadata_value(
        record,
        :left_stratum_kind,
    ),
    right_stratum_kind = _white_lindsey_record_metadata_value(
        record,
        :right_stratum_kind,
    ),
    left_plan = _white_lindsey_stratum_kernel_plan(left_stratum_kind),
    right_plan = _white_lindsey_stratum_kernel_plan(right_stratum_kind),
)
    return _white_lindsey_boundary_stratum_adapter_summary(
        record,
        status,
        blocker,
        left_stratum_kind,
        right_stratum_kind,
        left_plan,
        right_plan,
    )
end

function _white_lindsey_record_metadata_value(
    record::PairBlockMaterializationRecord,
    key::Symbol,
    default = nothing,
)
    return haskey(record.metadata, key) ? getfield(record.metadata, key) : default
end

function _white_lindsey_stratum_kernel_plan(stratum_kind)
    isnothing(stratum_kind) &&
        return _white_lindsey_unavailable_stratum_kernel_plan(
            :missing_white_lindsey_stratum_kind,
        )
    stratum_kind in (:facet_cpb, :face_cpb) && return (;
        status = :available_white_lindsey_stratum_kernel_plan,
        blocker = nothing,
        kernel = :_nested_face_product,
        side_1d_helper = :_nested_doside_1d,
    )
    stratum_kind === :edge_cpb && return (;
        status = :available_white_lindsey_stratum_kernel_plan,
        blocker = nothing,
        kernel = :_nested_edge_product,
        side_1d_helper = :_nested_doside_1d,
    )
    stratum_kind === :corner_cpb && return (;
        status = :available_white_lindsey_stratum_kernel_plan,
        blocker = nothing,
        kernel = :_nested_corner_piece,
        side_1d_helper = nothing,
    )
    return _white_lindsey_unavailable_stratum_kernel_plan(
        :unknown_white_lindsey_stratum_kind,
    )
end

function _white_lindsey_unavailable_stratum_kernel_plan(blocker::Symbol)
    return (;
        status = :blocked_white_lindsey_stratum_kernel_plan,
        blocker,
        kernel = nothing,
        side_1d_helper = nothing,
    )
end

function _white_lindsey_shared_1d_helper(left_plan, right_plan)
    left_plan.side_1d_helper === :_nested_doside_1d && return :_nested_doside_1d
    right_plan.side_1d_helper === :_nested_doside_1d && return :_nested_doside_1d
    return nothing
end

function _white_lindsey_count_by_value(values, field::Symbol)
    counts = Dict{Any,Int}()
    order = Any[]
    for value in values
        if !haskey(counts, value)
            counts[value] = 0
            push!(order, value)
        end
        counts[value] += 1
    end
    return Tuple((; field => value, count = counts[value]) for value in order)
end

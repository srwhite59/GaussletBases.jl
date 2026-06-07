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

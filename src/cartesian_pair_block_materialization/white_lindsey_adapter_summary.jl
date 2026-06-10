# Metadata-only old-kernel reuse summary for White--Lindsey boundary strata.

const _WHITE_LINDSEY_UNIT_DESCRIPTOR_AXES = (:x, :y, :z)

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

"""
    white_lindsey_boundary_stratum_unit_adapter_descriptor(unit)

Return compact metadata describing the source-CPB and old-kernel inputs for one
White--Lindsey boundary-stratum retained unit. This does not call old kernels
and does not build coefficient maps, LW numerical blocks, doside transforms,
Hamiltonian data, exports, artifacts, IDA/MWG data, or Coulomb.
"""
function white_lindsey_boundary_stratum_unit_adapter_descriptor(unit)
    stratum_kind = _white_lindsey_unit_metadata_value(unit, :stratum_kind)
    plan = _white_lindsey_stratum_kernel_plan(stratum_kind)
    source_cpb_count = length(unit.source_cpbs)
    source_cpb = source_cpb_count == 1 ? only(unit.source_cpbs) : nothing

    status, blocker = _white_lindsey_unit_descriptor_status(
        unit,
        stratum_kind,
        plan,
        source_cpb,
        source_cpb_count,
    )

    return (;
        object_kind =
            :white_lindsey_boundary_stratum_unit_adapter_descriptor,
        status,
        blocker,
        unit_key = unit.unit_key,
        unit_index = unit.unit_index,
        unit_kind = unit.unit_kind,
        lowering_kind = unit.lowering_kind,
        retained_rule = unit.retained_rule,
        source_contract_key = unit.source_contract_key,
        terminal_region_key = unit.terminal_region_key,
        stratum_kind,
        source_cpb_index = unit.source_cpb_index,
        source_cpb_count,
        source_cpb_roles = Tuple(CPB.role(cpb) for cpb in unit.source_cpbs),
        source_cpb_role =
            isnothing(source_cpb) ? nothing : CPB.role(source_cpb),
        source_cpb_codimension =
            isnothing(source_cpb) ? nothing : CPB.codimension(source_cpb),
        source_cpb_shape =
            isnothing(source_cpb) ? nothing : CPB.shape(source_cpb),
        source_cpb_intervals =
            isnothing(source_cpb) ? nothing : CPB.intervals(source_cpb),
        source_axis_intervals =
            _white_lindsey_source_cpb_axis_intervals(source_cpb),
        active_product_axes =
            _white_lindsey_source_cpb_active_product_axes(source_cpb),
        active_product_axis_intervals =
            _white_lindsey_source_cpb_active_product_axis_intervals(source_cpb),
        free_axis = _white_lindsey_source_cpb_free_axis(source_cpb),
        free_axis_interval =
            _white_lindsey_source_cpb_free_axis_interval(source_cpb),
        fixed_axes = _white_lindsey_source_cpb_fixed_axes(source_cpb),
        fixed_axis_coordinates =
            _white_lindsey_source_cpb_fixed_axis_coordinates(source_cpb),
        fixed_side_metadata =
            _white_lindsey_source_cpb_fixed_side_metadata(source_cpb),
        retained_count_status = unit.dimension_status,
        retained_count = unit.dimension,
        retained_counts = _white_lindsey_unit_metadata_value(
            unit,
            :retained_counts,
        ),
        parent_dims = _white_lindsey_unit_metadata_value(unit, :parent_dims),
        doside_source_1d =
            _white_lindsey_unit_metadata_value(unit, :doside_source_1d),
        kernel_plan_status = plan.status,
        planned_old_kernel = plan.kernel,
        planned_1d_helper = plan.side_1d_helper,
        coefficient_maps_materialized = false,
        source_operator_blocks_materialized = false,
        final_pair_blocks_materialized = false,
        operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
    )
end

"""
    white_lindsey_boundary_stratum_pair_adapter_descriptor(record)
    white_lindsey_boundary_stratum_pair_adapter_descriptor(record, unit_pair)

Return compact metadata describing one White--Lindsey boundary-stratum adapter
pair. The `unit_pair` overload uses unit-level descriptors when available; the
record-only overload reports record-derived facts only. This does not call old
kernels or build coefficient maps, LW numerical blocks, Hamiltonian data,
exports, artifacts, IDA/MWG data, or Coulomb.
"""
function white_lindsey_boundary_stratum_pair_adapter_descriptor(
    record::PairBlockMaterializationRecord,
)
    return _white_lindsey_boundary_stratum_pair_adapter_descriptor(
        record,
        nothing,
        nothing,
        :pair_block_materialization_record_metadata,
    )
end

function white_lindsey_boundary_stratum_pair_adapter_descriptor(
    record::PairBlockMaterializationRecord,
    unit_pair::CUP.UnitPairRecord,
)
    if record.pair_key != unit_pair.pair_key
        return _white_lindsey_boundary_stratum_pair_adapter_descriptor(
            record,
            nothing,
            nothing,
            :unit_pair_context_mismatch,
            :blocked_white_lindsey_boundary_stratum_pair_adapter_descriptor,
            :white_lindsey_pair_record_unit_pair_mismatch,
        )
    end

    left_descriptor =
        white_lindsey_boundary_stratum_unit_adapter_descriptor(unit_pair.left_unit)
    right_descriptor =
        white_lindsey_boundary_stratum_unit_adapter_descriptor(unit_pair.right_unit)

    return _white_lindsey_boundary_stratum_pair_adapter_descriptor(
        record,
        left_descriptor,
        right_descriptor,
        :unit_pair_retained_unit_descriptors,
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
        reuse_metadata_available_count = available_count,
        reuse_metadata_blocked_count = blocked_count,
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

function _white_lindsey_boundary_stratum_pair_adapter_descriptor(
    record::PairBlockMaterializationRecord,
    left_descriptor,
    right_descriptor,
    descriptor_source::Symbol,
    status_override = nothing,
    blocker_override = nothing,
)
    left_stratum_kind = _white_lindsey_pair_descriptor_stratum_kind(
        record,
        left_descriptor,
        :left_stratum_kind,
    )
    right_stratum_kind = _white_lindsey_pair_descriptor_stratum_kind(
        record,
        right_descriptor,
        :right_stratum_kind,
    )
    left_plan = _white_lindsey_stratum_kernel_plan(left_stratum_kind)
    right_plan = _white_lindsey_stratum_kernel_plan(right_stratum_kind)
    pair_family_classification =
        _white_lindsey_pair_family_classification(
            left_stratum_kind,
            right_stratum_kind,
        )
    status, blocker =
        isnothing(status_override) ?
        _white_lindsey_pair_descriptor_status(
            record,
            left_descriptor,
            right_descriptor,
            left_plan,
            right_plan,
            pair_family_classification,
        ) :
        (status_override, blocker_override)

    return (;
        object_kind =
            :white_lindsey_boundary_stratum_pair_adapter_descriptor,
        status,
        blocker,
        descriptor_source,
        pair_key = record.pair_key,
        pair_index = record.pair_index,
        pair_family = record.pair_family,
        pair_family_classification,
        left_unit_key = record.pair_key[1],
        right_unit_key = record.pair_key[2],
        left_unit_kind = _white_lindsey_pair_descriptor_unit_kind(
            record,
            left_descriptor,
            :left_unit_kind,
        ),
        right_unit_kind = _white_lindsey_pair_descriptor_unit_kind(
            record,
            right_descriptor,
            :right_unit_kind,
        ),
        left_unit_descriptor_status =
            isnothing(left_descriptor) ? nothing : left_descriptor.status,
        right_unit_descriptor_status =
            isnothing(right_descriptor) ? nothing : right_descriptor.status,
        left_stratum_kind,
        right_stratum_kind,
        left_planned_old_kernel = left_plan.kernel,
        right_planned_old_kernel = right_plan.kernel,
        left_planned_1d_helper = left_plan.side_1d_helper,
        right_planned_1d_helper = right_plan.side_1d_helper,
        shared_1d_helper =
            _white_lindsey_shared_1d_helper(left_plan, right_plan),
        materialization_path = record.materialization_path,
        pair_block_readiness_status = record.readiness_status,
        pair_block_blocker = record.blocker,
        coefficient_maps_materialized = false,
        source_operator_blocks_materialized = false,
        final_pair_blocks_materialized = false,
        operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
    )
end

function _white_lindsey_pair_descriptor_status(
    record::PairBlockMaterializationRecord,
    left_descriptor,
    right_descriptor,
    left_plan,
    right_plan,
    pair_family_classification,
)
    record.materialization_path ===
        :white_lindsey_boundary_stratum_adapter_preflight || return (
        :blocked_white_lindsey_boundary_stratum_pair_adapter_descriptor,
        :not_white_lindsey_boundary_stratum_adapter_preflight,
    )
    if !isnothing(left_descriptor) &&
       left_descriptor.status !==
       :available_metadata_only_white_lindsey_unit_adapter_descriptor
        return (
            :blocked_white_lindsey_boundary_stratum_pair_adapter_descriptor,
            left_descriptor.blocker,
        )
    elseif !isnothing(right_descriptor) &&
           right_descriptor.status !==
           :available_metadata_only_white_lindsey_unit_adapter_descriptor
        return (
            :blocked_white_lindsey_boundary_stratum_pair_adapter_descriptor,
            right_descriptor.blocker,
        )
    elseif left_plan.status !== :available_white_lindsey_stratum_kernel_plan
        return (
            :blocked_white_lindsey_boundary_stratum_pair_adapter_descriptor,
            left_plan.blocker,
        )
    elseif right_plan.status !== :available_white_lindsey_stratum_kernel_plan
        return (
            :blocked_white_lindsey_boundary_stratum_pair_adapter_descriptor,
            right_plan.blocker,
        )
    elseif isnothing(pair_family_classification)
        return (
            :blocked_white_lindsey_boundary_stratum_pair_adapter_descriptor,
            :unknown_white_lindsey_pair_family_classification,
        )
    end
    return :available_metadata_only_white_lindsey_pair_adapter_descriptor, nothing
end

function _white_lindsey_pair_descriptor_stratum_kind(
    record::PairBlockMaterializationRecord,
    unit_descriptor,
    metadata_key::Symbol,
)
    isnothing(unit_descriptor) ||
        return unit_descriptor.stratum_kind
    return _white_lindsey_record_metadata_value(record, metadata_key)
end

function _white_lindsey_pair_descriptor_unit_kind(
    record::PairBlockMaterializationRecord,
    unit_descriptor,
    metadata_key::Symbol,
)
    isnothing(unit_descriptor) ||
        return unit_descriptor.unit_kind
    return _white_lindsey_record_metadata_value(record, metadata_key)
end

function _white_lindsey_unit_metadata_value(unit, key::Symbol, default = nothing)
    return haskey(unit.metadata, key) ? getfield(unit.metadata, key) : default
end

function _white_lindsey_unit_descriptor_status(
    unit,
    stratum_kind,
    plan,
    source_cpb,
    source_cpb_count::Int,
)
    unit.unit_kind === :white_lindsey_boundary_stratum_retained_unit || return (
        :blocked_white_lindsey_boundary_stratum_unit_adapter_descriptor,
        :not_white_lindsey_boundary_stratum_retained_unit,
    )
    plan.status === :available_white_lindsey_stratum_kernel_plan || return (
        :blocked_white_lindsey_boundary_stratum_unit_adapter_descriptor,
        plan.blocker,
    )
    source_cpb_count == 1 || return (
        :blocked_white_lindsey_boundary_stratum_unit_adapter_descriptor,
        :white_lindsey_unit_source_cpb_count_not_one,
    )
    _white_lindsey_stratum_codimension_matches(
        stratum_kind,
        CPB.codimension(source_cpb),
    ) || return (
        :blocked_white_lindsey_boundary_stratum_unit_adapter_descriptor,
        :white_lindsey_stratum_codimension_mismatch,
    )
    return :available_metadata_only_white_lindsey_unit_adapter_descriptor, nothing
end

function _white_lindsey_stratum_codimension_matches(stratum_kind, codimension)
    stratum_kind === :direct_core && return codimension == 0
    stratum_kind in (:facet_cpb, :face_cpb) && return codimension == 1
    stratum_kind === :edge_cpb && return codimension == 2
    stratum_kind === :corner_cpb && return codimension == 3
    return false
end

function _white_lindsey_pair_family_classification(
    left_stratum_kind,
    right_stratum_kind,
)
    left_family = _white_lindsey_stratum_family(left_stratum_kind)
    right_family = _white_lindsey_stratum_family(right_stratum_kind)
    isnothing(left_family) && return nothing
    isnothing(right_family) && return nothing
    left_rank = _white_lindsey_stratum_family_rank(left_family)
    right_rank = _white_lindsey_stratum_family_rank(right_family)
    first_family, second_family =
        left_rank <= right_rank ?
        (left_family, right_family) :
        (right_family, left_family)
    return Symbol(String(first_family), "_", String(second_family))
end

function _white_lindsey_stratum_family(stratum_kind)
    stratum_kind === :direct_core && return :direct_core
    stratum_kind in (:facet_cpb, :face_cpb) && return :facet
    stratum_kind === :edge_cpb && return :edge
    stratum_kind === :corner_cpb && return :corner
    return nothing
end

function _white_lindsey_stratum_family_rank(stratum_family::Symbol)
    stratum_family === :direct_core && return 0
    stratum_family === :facet && return 1
    stratum_family === :edge && return 2
    stratum_family === :corner && return 3
    return typemax(Int)
end

function _white_lindsey_source_cpb_active_product_axes(source_cpb)
    isnothing(source_cpb) && return nothing
    return Tuple(
        axis for (axis, interval) in zip(
            _WHITE_LINDSEY_UNIT_DESCRIPTOR_AXES,
            CPB.intervals(source_cpb),
        )
        if length(interval) > 1
    )
end

function _white_lindsey_source_cpb_axis_intervals(source_cpb)
    isnothing(source_cpb) && return nothing
    intervals = CPB.intervals(source_cpb)
    return (;
        x = intervals[1],
        y = intervals[2],
        z = intervals[3],
    )
end

function _white_lindsey_source_cpb_active_product_axis_intervals(source_cpb)
    isnothing(source_cpb) && return nothing
    return Tuple(
        (; axis, interval)
        for (axis, interval) in zip(
            _WHITE_LINDSEY_UNIT_DESCRIPTOR_AXES,
            CPB.intervals(source_cpb),
        )
        if length(interval) > 1
    )
end

function _white_lindsey_source_cpb_free_axis(source_cpb)
    active_axes = _white_lindsey_source_cpb_active_product_axes(source_cpb)
    isnothing(active_axes) && return nothing
    length(active_axes) == 1 || return nothing
    return only(active_axes)
end

function _white_lindsey_source_cpb_free_axis_interval(source_cpb)
    active_intervals =
        _white_lindsey_source_cpb_active_product_axis_intervals(source_cpb)
    isnothing(active_intervals) && return nothing
    length(active_intervals) == 1 || return nothing
    return only(active_intervals).interval
end

function _white_lindsey_source_cpb_fixed_axes(source_cpb)
    isnothing(source_cpb) && return nothing
    return Tuple(
        axis for (axis, interval) in zip(
            _WHITE_LINDSEY_UNIT_DESCRIPTOR_AXES,
            CPB.intervals(source_cpb),
        )
        if length(interval) == 1
    )
end

function _white_lindsey_source_cpb_fixed_side_metadata(source_cpb)
    isnothing(source_cpb) && return nothing
    metadata = source_cpb.metadata
    if haskey(metadata, :axis) && haskey(metadata, :side)
        return ((; axis = metadata.axis, side = metadata.side),)
    elseif haskey(metadata, :fixed_axes) && haskey(metadata, :sides)
        return Tuple(
            (; axis, side)
            for (axis, side) in zip(metadata.fixed_axes, metadata.sides)
        )
    elseif haskey(metadata, :sides)
        return Tuple(
            (; axis, side)
            for (axis, side) in zip(
                _WHITE_LINDSEY_UNIT_DESCRIPTOR_AXES,
                metadata.sides,
            )
        )
    end
    return nothing
end

function _white_lindsey_source_cpb_fixed_axis_coordinates(source_cpb)
    isnothing(source_cpb) && return nothing
    return Tuple(
        (; axis, coordinate = first(interval))
        for (axis, interval) in zip(
            _WHITE_LINDSEY_UNIT_DESCRIPTOR_AXES,
            CPB.intervals(source_cpb),
        )
        if length(interval) == 1
    )
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
    stratum_kind === :direct_core && return (;
        status = :available_white_lindsey_stratum_kernel_plan,
        blocker = nothing,
        kernel = :direct_core_identity_support,
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

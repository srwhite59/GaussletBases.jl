# First tiny direct/direct overlap pair-block pilot.

"""
    direct_direct_overlap_block(record; parent_axis_counts, overlap_1d)

Materialize one direct/direct overlap pair block from explicit 1D overlap
matrices. This is a local pair-block pilot: it does not assemble global
operators, Hamiltonians, exports, or artifacts.
"""
function direct_direct_overlap_block(
    record::PairBlockMaterializationRecord;
    parent_axis_counts,
    overlap_1d,
)
    _assert_direct_direct_overlap_record(record)
    axis_counts = _axis_counts_tuple(parent_axis_counts)
    overlap_x, overlap_y, overlap_z = _overlap_1d_tuple(overlap_1d)
    _assert_overlap_axis_sizes((overlap_x, overlap_y, overlap_z), axis_counts)

    left_cpb = _direct_source_cpb(record, :left)
    right_cpb = _direct_source_cpb(record, :right)
    _assert_cpb_inside_parent(left_cpb, axis_counts, :left)
    _assert_cpb_inside_parent(right_cpb, axis_counts, :right)

    left_states = _cpb_support_states(left_cpb)
    right_states = _cpb_support_states(right_cpb)
    block = Matrix{Float64}(undef, length(left_states), length(right_states))
    _fill_direct_direct_overlap_block!(
        block,
        left_states,
        right_states,
        overlap_x,
        overlap_y,
        overlap_z,
    )

    return PairBlockMaterializationResult(
        :overlap,
        record.pair_key,
        block,
        true,
        true,
        true,
        false,
        false,
        false,
        (;
            materialization_path = record.materialization_path,
            readiness_status_before_materialization = record.readiness_status,
            parent_axis_counts = axis_counts,
            left_source_shape = CPB.shape(left_cpb),
            right_source_shape = CPB.shape(right_cpb),
        ),
    )
end

"""
    direct_direct_overlap_blocks(plan; parent_axis_counts, overlap_1d)

Materialize overlap only for ready direct/direct records in a pair-block
materialization plan. Unsupported or blocked records are returned as compact
skipped summaries.
"""
function direct_direct_overlap_blocks(
    plan::PairBlockMaterializationPlan;
    parent_axis_counts,
    overlap_1d,
)
    results = PairBlockMaterializationResult[]
    skipped = NamedTuple[]

    for record in pair_block_materialization_records(plan)
        if _is_ready_direct_direct_overlap_record(record)
            push!(
                results,
                direct_direct_overlap_block(
                    record;
                    parent_axis_counts,
                    overlap_1d,
                ),
            )
        else
            push!(skipped, _skipped_overlap_record_summary(record))
        end
    end

    result_tuple = Tuple(results)
    skipped_tuple = Tuple(skipped)
    any_materialized = !isempty(result_tuple)
    return PairBlockMaterializationBatchResult(
        :overlap,
        result_tuple,
        skipped_tuple,
        length(result_tuple),
        length(skipped_tuple),
        any_materialized,
        any_materialized,
        any_materialized,
        false,
        false,
        false,
        (;
            materialization_path = :ready_direct_direct_overlap_blocks_only,
            pair_block_record_count = length(pair_block_materialization_records(plan)),
        ),
    )
end

function _assert_direct_direct_overlap_record(record::PairBlockMaterializationRecord)
    record.materialization_path === :direct_direct_pair_block_materialization_pilot ||
        throw(ArgumentError("direct/direct overlap requires a direct/direct materialization record"))
    record.readiness_status === :ready_metadata_only_not_materialized ||
        throw(ArgumentError("direct/direct overlap requires a ready metadata-only record"))
    isnothing(record.blocker) ||
        throw(ArgumentError("direct/direct overlap record is blocked"))
    return nothing
end

function _is_ready_direct_direct_overlap_record(record::PairBlockMaterializationRecord)
    return record.materialization_path === :direct_direct_pair_block_materialization_pilot &&
           record.readiness_status === :ready_metadata_only_not_materialized
end

function _skipped_overlap_record_summary(record::PairBlockMaterializationRecord)
    return (;
        pair_key = record.pair_key,
        pair_index = record.pair_index,
        pair_family = record.pair_family,
        materialization_path = record.materialization_path,
        readiness_status = record.readiness_status,
        blocker =
            isnothing(record.blocker) ?
            :unsupported_direct_direct_overlap_materialization_record :
            record.blocker,
    )
end

function _axis_counts_tuple(parent_axis_counts)
    if parent_axis_counts isa Integer
        count = Int(parent_axis_counts)
        count > 0 ||
            throw(ArgumentError("parent_axis_counts must be positive"))
        return (count, count, count)
    end

    length(parent_axis_counts) == 3 ||
        throw(ArgumentError("parent_axis_counts must have three axes"))
    counts = Tuple(Int(parent_axis_counts[index]) for index in 1:3)
    all(count -> count > 0, counts) ||
        throw(ArgumentError("parent_axis_counts must be positive"))
    return counts
end

function _overlap_1d_tuple(overlap_1d)
    return _operator_1d_tuple(overlap_1d, "overlap_1d")
end

function _operator_1d_tuple(operator_1d, name::AbstractString)
    if operator_1d isa NamedTuple &&
       all(key -> haskey(operator_1d, key), (:x, :y, :z))
        return (operator_1d.x, operator_1d.y, operator_1d.z)
    end

    length(operator_1d) == 3 ||
        throw(ArgumentError("$(name) must provide x, y, and z matrices"))
    return (operator_1d[1], operator_1d[2], operator_1d[3])
end

function _assert_overlap_axis_sizes(overlap_axes, axis_counts::NTuple{3,Int})
    return _assert_operator_axis_sizes(overlap_axes, axis_counts, "overlap_1d")
end

function _assert_operator_axis_sizes(
    operator_axes,
    axis_counts::NTuple{3,Int},
    name::AbstractString,
)
    for (axis_index, overlap) in pairs(operator_axes)
        overlap isa AbstractMatrix{<:Real} ||
            throw(ArgumentError("$(name) entries must be real matrices"))
        size(overlap, 1) >= axis_counts[axis_index] &&
            size(overlap, 2) >= axis_counts[axis_index] ||
            throw(ArgumentError("$(name) matrix is smaller than parent_axis_counts"))
    end
    return nothing
end

function _direct_source_cpb(record::PairBlockMaterializationRecord, side::Symbol)
    source_key =
        side === :left ? :left_source_cpbs :
        side === :right ? :right_source_cpbs :
        throw(ArgumentError("direct source side must be :left or :right"))
    haskey(record.metadata, source_key) ||
        throw(ArgumentError("direct/direct overlap record is missing source CPB metadata"))
    source_cpbs = getfield(record.metadata, source_key)
    length(source_cpbs) == 1 ||
        throw(ArgumentError("direct/direct overlap requires exactly one source CPB per side"))
    source_cpb = only(source_cpbs)
    source_cpb isa CPB.CoordinateProductBox ||
        throw(ArgumentError("direct/direct overlap source metadata must contain CPBs"))
    return source_cpb
end

function _assert_cpb_inside_parent(
    source_cpb::CPB.CoordinateProductBox,
    axis_counts::NTuple{3,Int},
    side::Symbol,
)
    for axis_index in 1:3
        interval = CPB.intervals(source_cpb)[axis_index]
        first(interval) >= 1 && last(interval) <= axis_counts[axis_index] ||
            throw(ArgumentError("$(side) direct source CPB lies outside parent_axis_counts"))
    end
    return nothing
end

function _cpb_support_states(source_cpb::CPB.CoordinateProductBox)
    ix, iy, iz = CPB.intervals(source_cpb)
    return Tuple((x, y, z) for x in ix for y in iy for z in iz)
end

function _fill_direct_direct_overlap_block!(
    block::Matrix{Float64},
    left_states,
    right_states,
    overlap_x,
    overlap_y,
    overlap_z,
)
    @inbounds for (left_index, left_state) in pairs(left_states)
        lx, ly, lz = left_state
        for (right_index, right_state) in pairs(right_states)
            rx, ry, rz = right_state
            block[left_index, right_index] =
                Float64(overlap_x[lx, rx] * overlap_y[ly, ry] * overlap_z[lz, rz])
        end
    end
    return block
end

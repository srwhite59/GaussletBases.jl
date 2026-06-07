# Direct/direct position pair-block pilot.

"""
    direct_direct_position_block(record; axis, parent_axis_counts, overlap_1d, position_1d)

Materialize one direct/direct position pair block for `axis in (:x, :y, :z)`.
This is local pair-block data only, not global operator or Hamiltonian assembly.
"""
function direct_direct_position_block(
    record::PairBlockMaterializationRecord;
    axis,
    parent_axis_counts,
    overlap_1d,
    position_1d,
)
    term = _position_term(axis)
    _assert_direct_direct_overlap_record(record)
    axis_counts = _axis_counts_tuple(parent_axis_counts)
    overlap_x, overlap_y, overlap_z = _overlap_1d_tuple(overlap_1d)
    position_x, position_y, position_z =
        _operator_1d_tuple(position_1d, "position_1d")
    _assert_overlap_axis_sizes((overlap_x, overlap_y, overlap_z), axis_counts)
    _assert_operator_axis_sizes(
        (position_x, position_y, position_z),
        axis_counts,
        "position_1d",
    )

    operator_x, operator_y, operator_z =
        axis === :x ? (position_x, overlap_y, overlap_z) :
        axis === :y ? (overlap_x, position_y, overlap_z) :
        (overlap_x, overlap_y, position_z)

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
        operator_x,
        operator_y,
        operator_z,
    )

    return PairBlockMaterializationResult(
        term,
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
            position_axis = axis,
            left_source_shape = CPB.shape(left_cpb),
            right_source_shape = CPB.shape(right_cpb),
        ),
    )
end

"""
    direct_direct_position_blocks(plan; axis, parent_axis_counts, overlap_1d, position_1d)

Materialize position blocks only for ready direct/direct records in a plan.
Unsupported or blocked records are returned as compact skipped summaries.
"""
function direct_direct_position_blocks(
    plan::PairBlockMaterializationPlan;
    axis,
    parent_axis_counts,
    overlap_1d,
    position_1d,
)
    term = _position_term(axis)
    results = PairBlockMaterializationResult[]
    skipped = NamedTuple[]

    for record in pair_block_materialization_records(plan)
        if _is_ready_direct_direct_overlap_record(record)
            push!(
                results,
                direct_direct_position_block(
                    record;
                    axis,
                    parent_axis_counts,
                    overlap_1d,
                    position_1d,
                ),
            )
        else
            push!(skipped, _skipped_position_record_summary(record))
        end
    end

    result_tuple = Tuple(results)
    skipped_tuple = Tuple(skipped)
    any_materialized = !isempty(result_tuple)
    return PairBlockMaterializationBatchResult(
        term,
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
            materialization_path = :ready_direct_direct_position_blocks_only,
            position_axis = axis,
            pair_block_record_count = length(pair_block_materialization_records(plan)),
        ),
    )
end

function _position_term(axis)
    axis === :x && return :position_x
    axis === :y && return :position_y
    axis === :z && return :position_z
    throw(ArgumentError("position axis must be :x, :y, or :z"))
end

function _skipped_position_record_summary(record::PairBlockMaterializationRecord)
    return (;
        pair_key = record.pair_key,
        pair_index = record.pair_index,
        pair_family = record.pair_family,
        materialization_path = record.materialization_path,
        readiness_status = record.readiness_status,
        blocker =
            isnothing(record.blocker) ?
            :unsupported_direct_direct_position_materialization_record :
            record.blocker,
    )
end

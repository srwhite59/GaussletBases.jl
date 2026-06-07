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
    return _direct_direct_product_result(
        record,
        term,
        axis_counts,
        (operator_x, operator_y, operator_z),
        (; position_axis = axis),
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
    return _direct_direct_batch_results(
        record -> direct_direct_position_block(
            record;
            axis,
            parent_axis_counts,
            overlap_1d,
            position_1d,
        ),
        plan,
        term,
        :ready_direct_direct_position_blocks_only,
        _skipped_position_record_summary,
        (; position_axis = axis),
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

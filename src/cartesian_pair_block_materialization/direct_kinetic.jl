# Direct/direct kinetic pair-block pilot.

"""
    direct_direct_kinetic_block(record; parent_axis_counts, overlap_1d, kinetic_1d)

Materialize one direct/direct kinetic pair block from explicit 1D kinetic
matrices. The sign and prefactor convention belong to the supplied `kinetic_1d`
matrices.
"""
function direct_direct_kinetic_block(
    record::PairBlockMaterializationRecord;
    parent_axis_counts,
    overlap_1d,
    kinetic_1d,
)
    _assert_direct_direct_overlap_record(record)
    axis_counts = _axis_counts_tuple(parent_axis_counts)
    overlap_x, overlap_y, overlap_z = _overlap_1d_tuple(overlap_1d)
    kinetic_x, kinetic_y, kinetic_z =
        _operator_1d_tuple(kinetic_1d, "kinetic_1d")
    _assert_overlap_axis_sizes((overlap_x, overlap_y, overlap_z), axis_counts)
    _assert_operator_axis_sizes(
        (kinetic_x, kinetic_y, kinetic_z),
        axis_counts,
        "kinetic_1d",
    )

    kinetic_x_result = _direct_direct_product_result(
        record,
        :kinetic,
        axis_counts,
        (kinetic_x, overlap_y, overlap_z),
        (;),
    )
    kinetic_y_result = _direct_direct_product_result(
        record,
        :kinetic,
        axis_counts,
        (overlap_x, kinetic_y, overlap_z),
        (;),
    )
    kinetic_z_result = _direct_direct_product_result(
        record,
        :kinetic,
        axis_counts,
        (overlap_x, overlap_y, kinetic_z),
        (;),
    )

    return PairBlockMaterializationResult(
        :kinetic,
        record.pair_key,
        kinetic_x_result.block + kinetic_y_result.block + kinetic_z_result.block,
        true,
        true,
        true,
        false,
        false,
        false,
        kinetic_x_result.metadata,
    )
end

"""
    direct_direct_kinetic_blocks(plan; parent_axis_counts, overlap_1d, kinetic_1d)

Materialize kinetic blocks only for ready direct/direct records in a plan.
Unsupported or blocked records are returned as compact skipped summaries.
"""
function direct_direct_kinetic_blocks(
    plan::PairBlockMaterializationPlan;
    parent_axis_counts,
    overlap_1d,
    kinetic_1d,
)
    return _direct_direct_batch_results(
        record -> direct_direct_kinetic_block(
            record;
            parent_axis_counts,
            overlap_1d,
            kinetic_1d,
        ),
        plan,
        :kinetic,
        :ready_direct_direct_kinetic_blocks_only,
        _skipped_kinetic_record_summary,
        (;),
    )
end

function _skipped_kinetic_record_summary(record::PairBlockMaterializationRecord)
    return (;
        pair_key = record.pair_key,
        pair_index = record.pair_index,
        pair_family = record.pair_family,
        materialization_path = record.materialization_path,
        readiness_status = record.readiness_status,
        blocker =
            isnothing(record.blocker) ?
            :unsupported_direct_direct_kinetic_materialization_record :
            record.blocker,
    )
end

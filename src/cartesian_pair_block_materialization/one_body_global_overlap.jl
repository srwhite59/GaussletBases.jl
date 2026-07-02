# Dense global overlap assembly pilot from local placement records.
#
# This is deliberately overlap-only and consumes only placement records already
# marked placeable in final retained space. It ignores blocked records except
# for counts. It does not assemble Hamiltonians, kinetic/position/x2,
# Coulomb/IDA/MWG data, exports, artifacts, or PQS shell/Lowdin realization.

function one_body_global_overlap_matrix(placement_plan::NamedTuple)
    _one_body_assert_overlap_placement_plan_object(placement_plan)

    placement_plan.term === :overlap ||
        return _one_body_global_overlap_blocked_result(
            placement_plan,
            :non_overlap_placement_plan,
        )

    if _one_body_placement_value(placement_plan, :global_dimension_status, nothing) !==
       :available ||
       isnothing(_one_body_placement_value(placement_plan, :global_dimension, nothing))
        return _one_body_global_overlap_blocked_result(
            placement_plan,
            :missing_global_dimension,
        )
    end

    dimension = _one_body_placement_global_dimension(placement_plan.global_dimension)
    isempty(placement_plan.placeable_records) &&
        return _one_body_global_overlap_blocked_result(
            placement_plan,
            :no_placeable_overlap_blocks,
        )

    validation_blocker =
        _one_body_global_overlap_validation_blocker(placement_plan, dimension)
    isnothing(validation_blocker) ||
        return _one_body_global_overlap_blocked_result(
            placement_plan,
            validation_blocker,
        )

    placed = _one_body_global_symmetric_matrix_and_count(placement_plan, dimension)

    return _one_body_global_overlap_result(
        placement_plan,
        :materialized_global_overlap_matrix,
        nothing,
        placed.matrix;
        placed_block_count = placed.placed_block_count,
        skipped_block_count = placement_plan.blocked_count,
        materialized = true,
    )
end

function one_body_global_overlap_matrix(placement_plan)
    throw(
        ArgumentError(
            "global overlap matrix assembly requires a local one-body placement plan NamedTuple",
        ),
    )
end

const global_overlap_matrix = one_body_global_overlap_matrix

function _one_body_assert_overlap_placement_plan_object(placement_plan::NamedTuple)
    return _one_body_assert_global_placement_plan_object(placement_plan, :overlap)
end

function _one_body_global_overlap_validation_blocker(
    placement_plan::NamedTuple,
    dimension::Int,
)
    return _one_body_global_symmetric_validation_blocker(
        placement_plan,
        dimension,
        :non_symmetric_overlap_placement_record,
    )
end

function _one_body_global_overlap_blocked_result(
    placement_plan::NamedTuple,
    blocker::Symbol,
)
    return _one_body_global_overlap_result(
        placement_plan,
        :blocked_global_overlap_matrix,
        blocker,
        nothing;
        placed_block_count = 0,
        skipped_block_count =
            _one_body_placement_value(placement_plan, :record_count, 0),
        materialized = false,
    )
end

function _one_body_global_overlap_result(
    placement_plan::NamedTuple,
    status::Symbol,
    blocker,
    matrix;
    placed_block_count::Int,
    skipped_block_count::Int,
    materialized::Bool,
)
    return (;
        object_kind = :cartesian_pair_block_global_overlap_matrix,
        status,
        blocker,
        term = :overlap,
        global_dimension =
            _one_body_placement_value(placement_plan, :global_dimension, nothing),
        matrix,
        placeable_record_count =
            _one_body_placement_value(placement_plan, :placeable_count, 0),
        blocked_record_count =
            _one_body_placement_value(placement_plan, :blocked_count, 0),
        placed_block_count,
        skipped_block_count,
        symmetry = :symmetric,
        placement_plan_status =
            _one_body_placement_value(placement_plan, :status, nothing),
        placement_plan_blocker =
            _one_body_placement_value(placement_plan, :blocker, nothing),
        operator_matrix_materialized = materialized,
        global_overlap_matrix_materialized = materialized,
        global_operator_assembled = materialized,
    )
end

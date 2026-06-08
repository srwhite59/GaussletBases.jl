# Shared private helpers for dense global one-body matrix pilots.
#
# These helpers cover only behavior-neutral placement mechanics: local
# placement-plan object checks, placeable-record validation, and symmetric
# insertion of already-placeable final-local blocks. Term-specific entry points,
# result object kinds, statuses, blockers, materialization flags, and supported
# term sets stay owned by the individual overlap/kinetic/position/x2 files.

function _one_body_assert_global_placement_plan_object(
    placement_plan::NamedTuple,
    term::Symbol,
)
    _one_body_placement_value(placement_plan, :object_kind, nothing) ===
    :cartesian_pair_block_local_one_body_placement_plan || throw(
        ArgumentError(
            "global $(term) matrix assembly requires a local one-body placement plan",
        ),
    )
    return nothing
end

function _one_body_global_symmetric_validation_blocker(
    placement_plan::NamedTuple,
    dimension::Int,
    non_symmetric_blocker::Symbol,
)
    for record in placement_plan.placeable_records
        symmetry = _one_body_placement_value(record, :symmetry, :symmetric)
        symmetry === :symmetric || return non_symmetric_blocker

        result = _one_body_placement_value(record, :result, nothing)
        result isa PairBlockMaterializationResult ||
            return :missing_local_block_result

        block = result.block
        block isa AbstractMatrix{<:Real} ||
            return :missing_local_block_result

        target_ranges = _one_body_placement_value(record, :target_ranges, nothing)
        target_ranges isa NamedTuple ||
            return :missing_column_ranges
        haskey(target_ranges, :rows) && haskey(target_ranges, :columns) ||
            return :missing_column_ranges
        rows = target_ranges.rows
        columns = target_ranges.columns
        (
            _one_body_placement_valid_column_range(rows) &&
            _one_body_placement_valid_column_range(columns)
        ) || return :missing_column_ranges

        size(block) == (length(rows), length(columns)) ||
            return :local_block_shape_mismatch

        _one_body_placement_ranges_inside_global_dimension(
            rows,
            columns,
            dimension,
        ) || return :target_ranges_outside_global_dimension
    end
    return nothing
end

function _one_body_global_symmetric_matrix_and_count(
    placement_plan::NamedTuple,
    dimension::Int,
)
    matrix = zeros(Float64, dimension, dimension)
    placed_block_count = 0
    for record in placement_plan.placeable_records
        block = record.result.block
        rows = record.target_ranges.rows
        columns = record.target_ranges.columns
        matrix[rows, columns] .+= block
        placed_block_count += 1
        rows == columns && continue
        matrix[columns, rows] .+= transpose(block)
        placed_block_count += 1
    end
    return (; matrix, placed_block_count)
end

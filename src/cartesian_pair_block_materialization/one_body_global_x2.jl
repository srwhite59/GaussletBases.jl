# Dense global x2-axis assembly pilots from local placement records.
#
# This is deliberately limited to x2_x. It consumes only placement records
# already marked placeable in final retained space and keeps x2 axis terms
# separated. It ignores blocked records except for counts. It does not assemble
# Hamiltonians, overlap/kinetic/position, x2_y/z, Coulomb/IDA/MWG data, exports,
# artifacts, or PQS shell/Lowdin realization.

function one_body_global_x2_x_matrix(placement_plan::NamedTuple)
    _one_body_assert_x2_placement_plan_object(placement_plan, :x2_x)

    placement_plan.term === :x2_x ||
        return _one_body_global_x2_blocked_result(
            placement_plan,
            :x2_x,
            :non_x2_x_placement_plan,
        )

    if _one_body_placement_value(placement_plan, :global_dimension_status, nothing) !==
       :available ||
       isnothing(_one_body_placement_value(placement_plan, :global_dimension, nothing))
        return _one_body_global_x2_blocked_result(
            placement_plan,
            :x2_x,
            :missing_global_dimension,
        )
    end

    dimension = _one_body_placement_global_dimension(placement_plan.global_dimension)
    isempty(placement_plan.placeable_records) &&
        return _one_body_global_x2_blocked_result(
            placement_plan,
            :x2_x,
            :no_placeable_x2_x_blocks,
        )

    validation_blocker =
        _one_body_global_x2_validation_blocker(placement_plan, :x2_x, dimension)
    isnothing(validation_blocker) ||
        return _one_body_global_x2_blocked_result(
            placement_plan,
            :x2_x,
            validation_blocker,
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

    return _one_body_global_x2_result(
        placement_plan,
        :x2_x,
        :materialized_global_x2_x_matrix,
        nothing,
        matrix;
        placed_block_count,
        skipped_block_count = placement_plan.blocked_count,
        materialized = true,
    )
end

function one_body_global_x2_x_matrix(placement_plan)
    throw(
        ArgumentError(
            "global x2_x matrix assembly requires a local one-body placement plan NamedTuple",
        ),
    )
end

const global_x2_x_matrix = one_body_global_x2_x_matrix

function _one_body_assert_x2_placement_plan_object(
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

function _one_body_global_x2_validation_blocker(
    placement_plan::NamedTuple,
    term::Symbol,
    dimension::Int,
)
    for record in placement_plan.placeable_records
        symmetry = _one_body_placement_value(record, :symmetry, :symmetric)
        symmetry === :symmetric ||
            return Symbol("non_symmetric_", String(term), "_placement_record")

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

function _one_body_global_x2_blocked_result(
    placement_plan::NamedTuple,
    term::Symbol,
    blocker::Symbol,
)
    return _one_body_global_x2_result(
        placement_plan,
        term,
        Symbol("blocked_global_", String(term), "_matrix"),
        blocker,
        nothing;
        placed_block_count = 0,
        skipped_block_count =
            _one_body_placement_value(placement_plan, :record_count, 0),
        materialized = false,
    )
end

function _one_body_global_x2_result(
    placement_plan::NamedTuple,
    term::Symbol,
    status::Symbol,
    blocker,
    matrix;
    placed_block_count::Int,
    skipped_block_count::Int,
    materialized::Bool,
)
    term_string = String(term)
    term_flag = NamedTuple{
        (Symbol("global_", term_string, "_matrix_materialized"),),
    }((materialized,))
    base = (;
        object_kind = Symbol("cartesian_pair_block_global_", term_string, "_matrix"),
        status,
        blocker,
        term,
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
    )
    nonclaims = (;
        global_operator_assembled = materialized,
        operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
        exports_materialized = false,
        global_operator_blocks_materialized = false,
        global_hamiltonian_data_materialized = false,
        global_artifacts_materialized = false,
        coulomb_materialized = false,
        density_density_materialized = false,
        ida_mwg_data_materialized = false,
        pqs_lowdin_materialized = false,
        pqs_shell_projection_materialized = false,
        full_white_lindsey_route_assembled = false,
    )
    return merge(base, term_flag, nonclaims)
end

# Dense global position-axis assembly pilots from local placement records.
#
# This is deliberately limited to position_x, position_y, and position_z. It
# consumes only placement records already marked placeable in final retained
# space and keeps axis terms separated. It ignores blocked records except for
# counts. It does not assemble Hamiltonians, overlap/kinetic/x2,
# Coulomb/IDA/MWG data, exports, artifacts, or PQS shell/Lowdin realization.

function one_body_global_position_x_matrix(placement_plan::NamedTuple)
    return _one_body_global_position_matrix(placement_plan, :position_x)
end

function one_body_global_position_y_matrix(placement_plan::NamedTuple)
    return _one_body_global_position_matrix(placement_plan, :position_y)
end

function one_body_global_position_z_matrix(placement_plan::NamedTuple)
    return _one_body_global_position_matrix(placement_plan, :position_z)
end

function one_body_global_position_x_matrix(placement_plan)
    _one_body_global_position_throw_requires_plan(:position_x)
end

function one_body_global_position_y_matrix(placement_plan)
    _one_body_global_position_throw_requires_plan(:position_y)
end

function one_body_global_position_z_matrix(placement_plan)
    _one_body_global_position_throw_requires_plan(:position_z)
end

const global_position_x_matrix = one_body_global_position_x_matrix
const global_position_y_matrix = one_body_global_position_y_matrix
const global_position_z_matrix = one_body_global_position_z_matrix

function _one_body_global_position_matrix(
    placement_plan::NamedTuple,
    term::Symbol,
)
    _one_body_assert_position_term(term)
    _one_body_assert_position_placement_plan_object(placement_plan, term)

    placement_plan.term === term ||
        return _one_body_global_position_blocked_result(
            placement_plan,
            term,
            Symbol("non_", String(term), "_placement_plan"),
        )

    if _one_body_placement_value(placement_plan, :global_dimension_status, nothing) !==
       :available ||
       isnothing(_one_body_placement_value(placement_plan, :global_dimension, nothing))
        return _one_body_global_position_blocked_result(
            placement_plan,
            term,
            :missing_global_dimension,
        )
    end

    dimension = _one_body_placement_global_dimension(placement_plan.global_dimension)
    isempty(placement_plan.placeable_records) &&
        return _one_body_global_position_blocked_result(
            placement_plan,
            term,
            Symbol("no_placeable_", String(term), "_blocks"),
        )

    validation_blocker =
        _one_body_global_position_validation_blocker(
            placement_plan,
            term,
            dimension,
        )
    isnothing(validation_blocker) ||
        return _one_body_global_position_blocked_result(
            placement_plan,
            term,
            validation_blocker,
        )

    placed = _one_body_global_symmetric_matrix_and_count(placement_plan, dimension)

    return _one_body_global_position_result(
        placement_plan,
        term,
        Symbol("materialized_global_", String(term), "_matrix"),
        nothing,
        placed.matrix;
        placed_block_count = placed.placed_block_count,
        skipped_block_count = placement_plan.blocked_count,
        materialized = true,
    )
end

function _one_body_global_position_throw_requires_plan(term::Symbol)
    throw(
        ArgumentError(
            "global $(term) matrix assembly requires a local one-body placement plan NamedTuple",
        ),
    )
end

function _one_body_assert_position_term(term::Symbol)
    term in (:position_x, :position_y, :position_z) ||
        throw(ArgumentError("global position matrix assembly supports position_x/y/z only"))
    return nothing
end

function _one_body_assert_position_placement_plan_object(
    placement_plan::NamedTuple,
    term::Symbol,
)
    return _one_body_assert_global_placement_plan_object(placement_plan, term)
end

function _one_body_global_position_validation_blocker(
    placement_plan::NamedTuple,
    term::Symbol,
    dimension::Int,
)
    return _one_body_global_symmetric_validation_blocker(
        placement_plan,
        dimension,
        Symbol("non_symmetric_", String(term), "_placement_record"),
    )
end

function _one_body_global_position_blocked_result(
    placement_plan::NamedTuple,
    term::Symbol,
    blocker::Symbol,
)
    return _one_body_global_position_result(
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

function _one_body_global_position_result(
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
        global_operator_assembled = materialized,
    )
    return merge(base, term_flag)
end

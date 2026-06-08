# Dense global x2-axis assembly pilots from local placement records.
#
# This is deliberately limited to x2_x, x2_y, and x2_z. It consumes only
# placement records already marked placeable in final retained space and keeps
# x2 axis terms separated. It ignores blocked records except for counts. It does
# not assemble Hamiltonians, overlap/kinetic/position, Coulomb/IDA/MWG data,
# exports, artifacts, or PQS shell/Lowdin realization.

function one_body_global_x2_x_matrix(placement_plan::NamedTuple)
    return _one_body_global_x2_matrix(placement_plan, :x2_x)
end

function one_body_global_x2_y_matrix(placement_plan::NamedTuple)
    return _one_body_global_x2_matrix(placement_plan, :x2_y)
end

function one_body_global_x2_z_matrix(placement_plan::NamedTuple)
    return _one_body_global_x2_matrix(placement_plan, :x2_z)
end

function one_body_global_x2_x_matrix(placement_plan)
    _one_body_global_x2_throw_requires_plan(:x2_x)
end

function one_body_global_x2_y_matrix(placement_plan)
    _one_body_global_x2_throw_requires_plan(:x2_y)
end

function one_body_global_x2_z_matrix(placement_plan)
    _one_body_global_x2_throw_requires_plan(:x2_z)
end

const global_x2_x_matrix = one_body_global_x2_x_matrix
const global_x2_y_matrix = one_body_global_x2_y_matrix
const global_x2_z_matrix = one_body_global_x2_z_matrix

function _one_body_global_x2_matrix(placement_plan::NamedTuple, term::Symbol)
    _one_body_assert_x2_term(term)
    _one_body_assert_x2_placement_plan_object(placement_plan, term)

    placement_plan.term === term ||
        return _one_body_global_x2_blocked_result(
            placement_plan,
            term,
            Symbol("non_", String(term), "_placement_plan"),
        )

    if _one_body_placement_value(placement_plan, :global_dimension_status, nothing) !==
       :available ||
       isnothing(_one_body_placement_value(placement_plan, :global_dimension, nothing))
        return _one_body_global_x2_blocked_result(
            placement_plan,
            term,
            :missing_global_dimension,
        )
    end

    dimension = _one_body_placement_global_dimension(placement_plan.global_dimension)
    isempty(placement_plan.placeable_records) &&
        return _one_body_global_x2_blocked_result(
            placement_plan,
            term,
            Symbol("no_placeable_", String(term), "_blocks"),
        )

    validation_blocker =
        _one_body_global_x2_validation_blocker(placement_plan, term, dimension)
    isnothing(validation_blocker) ||
        return _one_body_global_x2_blocked_result(
            placement_plan,
            term,
            validation_blocker,
        )

    placed = _one_body_global_symmetric_matrix_and_count(placement_plan, dimension)

    return _one_body_global_x2_result(
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

function _one_body_global_x2_throw_requires_plan(term::Symbol)
    throw(
        ArgumentError(
            "global $(term) matrix assembly requires a local one-body placement plan NamedTuple",
        ),
    )
end

function _one_body_assert_x2_term(term::Symbol)
    term in (:x2_x, :x2_y, :x2_z) ||
        throw(ArgumentError("global x2 matrix assembly supports x2_x/y/z only"))
    return nothing
end

function _one_body_assert_x2_placement_plan_object(
    placement_plan::NamedTuple,
    term::Symbol,
)
    return _one_body_assert_global_placement_plan_object(placement_plan, term)
end

function _one_body_global_x2_validation_blocker(
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

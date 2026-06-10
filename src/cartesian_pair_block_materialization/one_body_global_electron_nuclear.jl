# Dense global by-center electron-nuclear assembly pilot from local placement
# records.
#
# This consumes only placement records already marked placeable in final
# retained space for one center. It keeps center records separated and does not
# apply nuclear charges, sum centers, assemble Hamiltonians, export artifacts,
# or use IDA/MWG/PQS semantics.

function one_body_global_electron_nuclear_by_center_matrix(
    placement_plan::NamedTuple,
)
    _one_body_assert_electron_nuclear_by_center_placement_plan_object(
        placement_plan,
    )

    placement_plan.term === :electron_nuclear_by_center ||
        return _one_body_global_electron_nuclear_by_center_blocked_result(
            placement_plan,
            :non_electron_nuclear_by_center_placement_plan,
        )

    if _one_body_placement_value(placement_plan, :global_dimension_status, nothing) !==
       :available ||
       isnothing(_one_body_placement_value(placement_plan, :global_dimension, nothing))
        return _one_body_global_electron_nuclear_by_center_blocked_result(
            placement_plan,
            :missing_global_dimension,
        )
    end

    dimension = _one_body_placement_global_dimension(placement_plan.global_dimension)
    isempty(placement_plan.placeable_records) &&
        return _one_body_global_electron_nuclear_by_center_blocked_result(
            placement_plan,
            :no_placeable_electron_nuclear_by_center_blocks,
        )

    validation_blocker =
        _one_body_global_electron_nuclear_by_center_validation_blocker(
            placement_plan,
            dimension,
        )
    isnothing(validation_blocker) ||
        return _one_body_global_electron_nuclear_by_center_blocked_result(
            placement_plan,
            validation_blocker,
        )

    placed = _one_body_global_symmetric_matrix_and_count(placement_plan, dimension)
    center_index, center_key, center_location, nuclear_charge =
        _one_body_global_electron_nuclear_by_center_identity(placement_plan)

    return _one_body_global_electron_nuclear_by_center_result(
        placement_plan,
        :materialized_global_electron_nuclear_by_center_matrix,
        nothing,
        placed.matrix;
        center_index,
        center_key,
        center_location,
        nuclear_charge,
        placed_block_count = placed.placed_block_count,
        skipped_block_count = placement_plan.blocked_count,
        materialized = true,
    )
end

function one_body_global_electron_nuclear_by_center_matrix(placement_plan)
    throw(
        ArgumentError(
            "global electron-nuclear by-center matrix assembly requires a local one-body placement plan NamedTuple",
        ),
    )
end

function _one_body_assert_electron_nuclear_by_center_placement_plan_object(
    placement_plan::NamedTuple,
)
    return _one_body_assert_global_placement_plan_object(
        placement_plan,
        :electron_nuclear_by_center,
    )
end

function _one_body_global_electron_nuclear_by_center_validation_blocker(
    placement_plan::NamedTuple,
    dimension::Int,
)
    center_count =
        length(_one_body_global_electron_nuclear_by_center_indices(placement_plan))
    center_count == 1 ||
        return :mixed_electron_nuclear_centers_in_by_center_matrix
    any(record -> record.nuclear_charge_applied, placement_plan.placeable_records) &&
        return :nuclear_charge_already_applied_to_by_center_blocks
    return _one_body_global_symmetric_validation_blocker(
        placement_plan,
        dimension,
        :non_symmetric_electron_nuclear_by_center_placement_record,
    )
end

function _one_body_global_electron_nuclear_by_center_indices(
    placement_plan::NamedTuple,
)
    values = Any[]
    for record in placement_plan.placeable_records
        value = _one_body_placement_value(record, :center_index, nothing)
        isnothing(value) && continue
        value in values && continue
        push!(values, value)
    end
    return Tuple(values)
end

function _one_body_global_electron_nuclear_by_center_identity(
    placement_plan::NamedTuple,
)
    record = first(placement_plan.placeable_records)
    return (
        _one_body_placement_value(record, :center_index, nothing),
        _one_body_placement_value(record, :center_key, nothing),
        _one_body_placement_value(record, :center_location, nothing),
        _one_body_placement_value(record, :nuclear_charge, nothing),
    )
end

function _one_body_global_electron_nuclear_by_center_blocked_result(
    placement_plan::NamedTuple,
    blocker::Symbol,
)
    return _one_body_global_electron_nuclear_by_center_result(
        placement_plan,
        :blocked_global_electron_nuclear_by_center_matrix,
        blocker,
        nothing;
        center_index = nothing,
        center_key = nothing,
        center_location = nothing,
        nuclear_charge = nothing,
        placed_block_count = 0,
        skipped_block_count =
            _one_body_placement_value(placement_plan, :record_count, 0),
        materialized = false,
    )
end

function _one_body_global_electron_nuclear_by_center_result(
    placement_plan::NamedTuple,
    status::Symbol,
    blocker,
    matrix;
    center_index,
    center_key,
    center_location,
    nuclear_charge,
    placed_block_count::Int,
    skipped_block_count::Int,
    materialized::Bool,
)
    return (;
        object_kind = :cartesian_pair_block_global_electron_nuclear_by_center_matrix,
        status,
        blocker,
        term = :electron_nuclear_by_center,
        center_index,
        center_key,
        center_location,
        nuclear_charge,
        by_center = true,
        centers_summed = false,
        nuclear_charge_recorded = !isnothing(center_index),
        nuclear_charge_applied = false,
        charge_application_stage = :acceptance_or_hamiltonian_assembly,
        matrix,
        global_dimension =
            _one_body_placement_value(placement_plan, :global_dimension, nothing),
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
        global_electron_nuclear_by_center_matrix_materialized = materialized,
        global_operator_assembled = materialized,
        operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
        exports_materialized = false,
        global_operator_blocks_materialized = false,
        global_hamiltonian_data_materialized = false,
        global_artifacts_materialized = false,
        density_density_materialized = false,
        ida_mwg_data_materialized = false,
        pqs_lowdin_materialized = false,
        pqs_shell_projection_materialized = false,
        full_white_lindsey_route_assembled = false,
    )
end

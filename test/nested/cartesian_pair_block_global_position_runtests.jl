# Runtime role: tiny smoke / global position_x assembly pilot.
#
# This validates position_x-only dense global assembly from already-placeable
# final-local placement records. It does not cover Hamiltonians, overlap,
# kinetic/x2, Coulomb, IDA/MWG, PQS realization, exports, or White-Lindsey
# oracle paths.

using Test
using GaussletBases

const CPBMGlobalPosition = GaussletBases.CartesianPairBlockMaterialization

function _global_position_result(
    term::Symbol,
    pair_key::Tuple{Symbol,Symbol},
    block;
    source_operator_blocks_materialized::Bool,
    final_pair_blocks_materialized::Bool,
    metadata = (;),
)
    return CPBMGlobalPosition.PairBlockMaterializationResult(
        term,
        pair_key,
        Matrix{Float64}(block),
        true,
        source_operator_blocks_materialized,
        final_pair_blocks_materialized,
        false,
        false,
        false,
        NamedTuple(metadata),
    )
end

function _global_position_entry(result; block_set_term = :position_x)
    return CPBMGlobalPosition._one_body_local_block_collection_entry(
        result;
        block_set_term,
    )
end

function _global_position_collection(; include_missing_ranges::Bool = true)
    diagonal = _global_position_entry(
        _global_position_result(
            :position_x,
            (:diag_left, :diag_left),
            [1.5 0.4; 0.4 2.5];
            source_operator_blocks_materialized = true,
            final_pair_blocks_materialized = true,
            metadata = (;
                pair_index = 1,
                selector_family = :direct_direct,
                materialization_path =
                    :direct_direct_pair_block_materialization_pilot,
                left_column_range = 1:2,
                right_column_range = 1:2,
            ),
        ),
    )
    off_diagonal = _global_position_entry(
        _global_position_result(
            :position_x,
            (:left, :right),
            [3.0 4.0; 5.0 6.0];
            source_operator_blocks_materialized = true,
            final_pair_blocks_materialized = true,
            metadata = (;
                pair_index = 2,
                selector_family = :direct_direct,
                materialization_path =
                    :direct_direct_pair_block_materialization_pilot,
                left_column_range = 1:2,
                right_column_range = 3:4,
            ),
        ),
    )
    pqs_source = _global_position_entry(
        _global_position_result(
            :source_position_x,
            (:pqs_left, :pqs_right),
            [10.0 11.0; 12.0 13.0];
            source_operator_blocks_materialized = true,
            final_pair_blocks_materialized = false,
            metadata = (;
                pair_index = 3,
                selector_family = :pqs_source_pair,
                materialization_path = :pqs_source_pair_preflight,
            ),
        ),
    )
    materialized_entries = (diagonal, off_diagonal, pqs_source)
    if include_missing_ranges
        missing_ranges = _global_position_entry(
            _global_position_result(
                :position_x,
                (:missing_left, :missing_right),
                [2.0 0.1; 0.1 2.2];
                source_operator_blocks_materialized = true,
                final_pair_blocks_materialized = true,
                metadata = (;
                    pair_index = 4,
                    selector_family = :direct_direct,
                    materialization_path =
                        :direct_direct_pair_block_materialization_pilot,
                ),
            ),
        )
        materialized_entries = (materialized_entries..., missing_ranges)
    end
    return (;
        object_kind = :cartesian_pair_block_local_one_body_block_collection,
        terms = (:position_x,),
        requested_terms = (:position_x,),
        materialized_terms = (:position_x,),
        deferred_terms = (),
        entries = materialized_entries,
        materialized_entries,
        skipped_entries = (),
    )
end

function _global_position_plan(;
    global_dimension = 4,
    include_missing_ranges::Bool = true,
)
    return CPBMGlobalPosition.one_body_position_x_placement_plan(
        _global_position_collection(; include_missing_ranges);
        global_dimension,
    )
end

function _global_position_count(counts, field::Symbol, value)
    matches = Tuple(entry for entry in counts if getproperty(entry, field) == value)
    isempty(matches) && return 0
    return sum(entry -> hasproperty(entry, :count) ? entry.count : 1, matches)
end

@testset "CartesianPairBlockMaterialization global position_x matrix" begin
    placement_plan = _global_position_plan()
    result = CPBMGlobalPosition.one_body_global_position_x_matrix(placement_plan)
    alias_result = CPBMGlobalPosition.global_position_x_matrix(placement_plan)

    expected = [
        1.5 0.4 3.0 4.0
        0.4 2.5 5.0 6.0
        3.0 5.0 0.0 0.0
        4.0 6.0 0.0 0.0
    ]
    @test result.object_kind === :cartesian_pair_block_global_position_x_matrix
    @test result.status === :materialized_global_position_x_matrix
    @test isnothing(result.blocker)
    @test result.term === :position_x
    @test result.global_dimension == 4
    @test result.matrix == expected
    @test alias_result.matrix == expected
    @test result.placeable_record_count == 2
    @test result.blocked_record_count == 2
    @test result.placed_block_count == 3
    @test result.skipped_block_count == 2
    @test result.symmetry === :symmetric

    @test _global_position_count(
        placement_plan.blocker_counts,
        :blocker,
        :source_space_block_requires_shell_realization,
    ) == 1
    @test _global_position_count(
        placement_plan.blocker_counts,
        :blocker,
        :missing_column_ranges,
    ) == 1
    @test !any(value -> value >= 10.0, result.matrix)

    missing_dimension_plan = CPBMGlobalPosition.one_body_position_x_placement_plan(
        _global_position_collection(),
    )
    missing_dimension =
        CPBMGlobalPosition.one_body_global_position_x_matrix(
            missing_dimension_plan,
        )
    @test missing_dimension.status === :blocked_global_position_x_matrix
    @test missing_dimension.blocker === :missing_global_dimension
    @test isnothing(missing_dimension.matrix)

    non_position = CPBMGlobalPosition.one_body_global_position_x_matrix(
        merge(placement_plan, (; term = :overlap)),
    )
    @test non_position.status === :blocked_global_position_x_matrix
    @test non_position.blocker === :non_position_x_placement_plan

    @test_throws ArgumentError CPBMGlobalPosition.one_body_placement_plan(
        _global_position_collection();
        term = :position_y,
    )

    placeable = only(
        record for record in placement_plan.placeable_records
        if record.pair_key == (:left, :right)
    )
    shape_mismatch_record = merge(
        placeable,
        (; target_ranges = (; rows = 1:3, columns = 3:4)),
    )
    shape_mismatch_plan = merge(
        placement_plan,
        (;
            placeable_records = (shape_mismatch_record,),
            placeable_count = 1,
            blocked_count = 0,
            blocked_records = (),
            record_count = 1,
        ),
    )
    shape_mismatch =
        CPBMGlobalPosition.one_body_global_position_x_matrix(
            shape_mismatch_plan,
        )
    @test shape_mismatch.status === :blocked_global_position_x_matrix
    @test shape_mismatch.blocker === :local_block_shape_mismatch

    outside_record = merge(
        placeable,
        (; target_ranges = (; rows = 1:2, columns = 4:5)),
    )
    outside_plan = merge(
        placement_plan,
        (;
            placeable_records = (outside_record,),
            placeable_count = 1,
            blocked_count = 0,
            blocked_records = (),
            record_count = 1,
        ),
    )
    outside = CPBMGlobalPosition.one_body_global_position_x_matrix(outside_plan)
    @test outside.status === :blocked_global_position_x_matrix
    @test outside.blocker === :target_ranges_outside_global_dimension

    empty_plan = merge(
        placement_plan,
        (;
            placeable_records = (),
            placeable_count = 0,
            blocked_count = placement_plan.record_count,
        ),
    )
    empty = CPBMGlobalPosition.one_body_global_position_x_matrix(empty_plan)
    @test empty.status === :blocked_global_position_x_matrix
    @test empty.blocker === :no_placeable_position_x_blocks

    @test result.operator_matrix_materialized
    @test result.global_position_x_matrix_materialized
    @test result.global_operator_assembled
    @test !result.operator_blocks_materialized
    @test !result.hamiltonian_data_materialized
    @test !result.artifacts_materialized
    @test !result.exports_materialized
    @test !result.global_operator_blocks_materialized
    @test !result.global_hamiltonian_data_materialized
    @test !result.global_artifacts_materialized
    @test !result.coulomb_materialized
    @test !result.density_density_materialized
    @test !result.ida_mwg_data_materialized
    @test !result.pqs_lowdin_materialized
    @test !result.pqs_shell_projection_materialized
    @test !result.full_white_lindsey_route_assembled
end

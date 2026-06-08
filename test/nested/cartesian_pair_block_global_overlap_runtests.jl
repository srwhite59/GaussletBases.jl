# Runtime role: tiny smoke / global overlap assembly pilot.
#
# This validates overlap-only dense global assembly from already-placeable
# final-local placement records. It does not cover Hamiltonians, Coulomb,
# IDA/MWG, PQS realization, exports, or White-Lindsey oracle paths.

using Test
using GaussletBases

const CPBMGlobalOverlap = GaussletBases.CartesianPairBlockMaterialization

function _global_overlap_result(
    term::Symbol,
    pair_key::Tuple{Symbol,Symbol},
    block;
    source_operator_blocks_materialized::Bool,
    final_pair_blocks_materialized::Bool,
    metadata = (;),
)
    return CPBMGlobalOverlap.PairBlockMaterializationResult(
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

function _global_overlap_entry(result; block_set_term = :overlap)
    return CPBMGlobalOverlap._one_body_local_block_collection_entry(
        result;
        block_set_term,
    )
end

function _global_overlap_collection()
    diagonal = _global_overlap_entry(
        _global_overlap_result(
            :overlap,
            (:diag_left, :diag_left),
            [1.0 0.2; 0.2 1.1];
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
    off_diagonal = _global_overlap_entry(
        _global_overlap_result(
            :overlap,
            (:left, :right),
            [2.0 3.0; 4.0 5.0];
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
    pqs_source = _global_overlap_entry(
        _global_overlap_result(
            :source_overlap,
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
    return (;
        object_kind = :cartesian_pair_block_local_one_body_block_collection,
        terms = (:overlap,),
        requested_terms = (:overlap,),
        materialized_terms = (:overlap,),
        deferred_terms = (),
        entries = materialized_entries,
        materialized_entries,
        skipped_entries = (),
    )
end

function _global_overlap_plan(; global_dimension = 4)
    return CPBMGlobalOverlap.one_body_overlap_placement_plan(
        _global_overlap_collection();
        global_dimension,
    )
end

@testset "CartesianPairBlockMaterialization global overlap matrix" begin
    placement_plan = _global_overlap_plan()
    result = CPBMGlobalOverlap.one_body_global_overlap_matrix(placement_plan)
    alias_result = CPBMGlobalOverlap.global_overlap_matrix(placement_plan)

    expected = [
        1.0 0.2 2.0 3.0
        0.2 1.1 4.0 5.0
        2.0 4.0 0.0 0.0
        3.0 5.0 0.0 0.0
    ]
    @test result.object_kind === :cartesian_pair_block_global_overlap_matrix
    @test result.status === :materialized_global_overlap_matrix
    @test isnothing(result.blocker)
    @test result.term === :overlap
    @test result.global_dimension == 4
    @test result.matrix == expected
    @test alias_result.matrix == expected
    @test result.placeable_record_count == 2
    @test result.blocked_record_count == 1
    @test result.placed_block_count == 3
    @test result.skipped_block_count == 1
    @test result.symmetry === :symmetric

    pqs_record = only(
        record for record in placement_plan.blocked_records
        if record.selector_family === :pqs_source_pair
    )
    @test pqs_record.blocker === :source_space_block_requires_shell_realization
    @test !any(value -> value >= 10.0, result.matrix)

    missing_dimension_plan = CPBMGlobalOverlap.one_body_overlap_placement_plan(
        _global_overlap_collection(),
    )
    missing_dimension =
        CPBMGlobalOverlap.one_body_global_overlap_matrix(missing_dimension_plan)
    @test missing_dimension.status === :blocked_global_overlap_matrix
    @test missing_dimension.blocker === :missing_global_dimension
    @test isnothing(missing_dimension.matrix)

    non_overlap = CPBMGlobalOverlap.one_body_global_overlap_matrix(
        merge(placement_plan, (; term = :kinetic)),
    )
    @test non_overlap.status === :blocked_global_overlap_matrix
    @test non_overlap.blocker === :non_overlap_placement_plan

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
        CPBMGlobalOverlap.one_body_global_overlap_matrix(shape_mismatch_plan)
    @test shape_mismatch.status === :blocked_global_overlap_matrix
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
    outside = CPBMGlobalOverlap.one_body_global_overlap_matrix(outside_plan)
    @test outside.status === :blocked_global_overlap_matrix
    @test outside.blocker === :target_ranges_outside_global_dimension

    nonsymmetric_record = merge(placeable, (; symmetry = :nonsymmetric))
    nonsymmetric_plan = merge(
        placement_plan,
        (;
            placeable_records = (nonsymmetric_record,),
            placeable_count = 1,
            blocked_count = 0,
            blocked_records = (),
            record_count = 1,
        ),
    )
    nonsymmetric =
        CPBMGlobalOverlap.one_body_global_overlap_matrix(nonsymmetric_plan)
    @test nonsymmetric.status === :blocked_global_overlap_matrix
    @test nonsymmetric.blocker === :non_symmetric_overlap_placement_record

    empty_plan = merge(
        placement_plan,
        (;
            placeable_records = (),
            placeable_count = 0,
            blocked_count = placement_plan.record_count,
        ),
    )
    empty = CPBMGlobalOverlap.one_body_global_overlap_matrix(empty_plan)
    @test empty.status === :blocked_global_overlap_matrix
    @test empty.blocker === :no_placeable_overlap_blocks

    @test result.operator_matrix_materialized
    @test result.global_overlap_matrix_materialized
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

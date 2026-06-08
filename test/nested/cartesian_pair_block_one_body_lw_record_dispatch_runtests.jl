using Test

include("cartesian_white_lindsey_adapter_fixture_helpers.jl")

const CPBMMixedLWDispatch = GaussletBases.CartesianPairBlockMaterialization

function _mixed_lw_dispatch_record(unit_pair)
    return CPBMMixedLWDispatch.PairBlockMaterializationRecord(
        unit_pair.pair_key,
        unit_pair.pair_index,
        unit_pair.pair_family,
        :synthetic_white_lindsey_source_operator_path,
        (; left = :white_lindsey_left_transform, right = :white_lindsey_right_transform),
        (; left = :white_lindsey_left_realization, right = :white_lindsey_right_realization),
        :synthetic_white_lindsey_final_block_path,
        (
            :overlap,
            :position_x,
            :position_y,
            :position_z,
            :x2_x,
            :x2_y,
            :x2_z,
            :kinetic,
        ),
        :white_lindsey_boundary_stratum_adapter_preflight,
        :ready_metadata_only_not_materialized,
        nothing,
        false,
        (; fixture = :mixed_lw_record_dispatch),
    )
end

function _check_mixed_lw_result(result, expected, term::Symbol)
    @test result isa CPBMMixedLWDispatch.PairBlockMaterializationResult
    @test result.term == term
    @test result.pair_key == expected.pair_key
    @test result.block == expected.block
    @test result.materialized
    @test result.source_operator_blocks_materialized
    @test result.final_pair_blocks_materialized
    @test !result.operator_blocks_materialized
    @test !result.hamiltonian_data_materialized
    @test !result.artifacts_materialized
    @test result.metadata.selector_family == :white_lindsey_boundary_stratum
    @test result.metadata.mixed_one_body_dispatcher ==
          :white_lindsey_boundary_stratum_record_dispatcher
    @test result.metadata.numerical_family_selector ==
          :white_lindsey_boundary_stratum_one_body_block
    @test !result.metadata.operator_blocks_materialized
    @test !result.metadata.hamiltonian_data_materialized
    @test !result.metadata.artifacts_materialized
    @test result.metadata.coefficient_source ==
          :white_lindsey_boundary_stratum_pair_unit_coefficients
end

@testset "CartesianPairBlockMaterialization mixed White-Lindsey record dispatch" begin
    fixture = _lw_adapter_prepared_facet_edge_fixture(;
        prefix = "mixed_lw_record_dispatch",
    )
    unit_pair = _lw_adapter_unit_pair(
        fixture.real_units.real_facet_unit,
        fixture.real_units.real_edge_unit,
        1,
    )
    record = _mixed_lw_dispatch_record(unit_pair)
    parent_axis_counts = (7, 7, 7)
    overlap_1d = fixture.factors.overlap_1d
    position_1d = fixture.factors.position_1d

    overlap = CPBMMixedLWDispatch._one_body_pair_block(
        record,
        :overlap;
        inputs = (; parent_axis_counts, overlap_1d),
        unit_pair,
    )
    expected_overlap =
        CPBMMixedLWDispatch.white_lindsey_boundary_stratum_one_body_block(
            unit_pair,
            :overlap;
            parent_axis_counts,
            overlap_1d,
        )
    _check_mixed_lw_result(overlap, expected_overlap, :overlap)
    @test overlap.metadata.materialization_path ==
          :white_lindsey_boundary_stratum_overlap_adapter

    position_y = CPBMMixedLWDispatch._one_body_pair_block(
        record,
        :position_y;
        inputs = (; parent_axis_counts, overlap_1d, position_1d),
        unit_pair,
    )
    expected_position_y =
        CPBMMixedLWDispatch.white_lindsey_boundary_stratum_one_body_block(
            unit_pair,
            :position_y;
            parent_axis_counts,
            overlap_1d,
            position_1d,
        )
    _check_mixed_lw_result(position_y, expected_position_y, :position_y)
    @test position_y.metadata.position_axis == :y

    missing_unit_pair = CPBMMixedLWDispatch._one_body_pair_block(
        record,
        :overlap;
        inputs = (; parent_axis_counts, overlap_1d),
    )
    @test missing_unit_pair.status == :skipped_mixed_one_body_pair_block
    @test missing_unit_pair.blocker == :missing_white_lindsey_unit_pair
    @test missing_unit_pair.selector_family == :white_lindsey_boundary_stratum
    @test !missing_unit_pair.materialized

    mismatched_unit_pair = CPBMMixedLWDispatch._one_body_pair_block(
        record,
        :overlap;
        inputs = (; parent_axis_counts, overlap_1d),
        unit_pair = (; pair_key = (:wrong_left, :wrong_right)),
    )
    @test mismatched_unit_pair.status == :skipped_mixed_one_body_pair_block
    @test mismatched_unit_pair.blocker == :mismatched_white_lindsey_unit_pair
    @test mismatched_unit_pair.selector_family == :white_lindsey_boundary_stratum
    @test !mismatched_unit_pair.materialized
end

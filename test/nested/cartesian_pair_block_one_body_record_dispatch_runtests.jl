using Test
using GaussletBases

const CPBMRecordDispatch = GaussletBases.CartesianPairBlockMaterialization
const CPBRecordDispatch = GaussletBases.CartesianCPB

function _record_dispatch_direct_record()
    left_source = CPBRecordDispatch.cpb(
        1:1,
        1:2,
        1:1;
        role = :mixed_one_body_direct_left_source,
    )
    right_source = CPBRecordDispatch.cpb(
        2:3,
        1:1,
        1:1;
        role = :mixed_one_body_direct_right_source,
    )
    return _record_dispatch_record(
        :direct_direct_pair_block_materialization_pilot;
        pair_key = (:direct_left, :direct_right),
        pair_family = :direct_direct,
        metadata = (;
            left_source_cpbs = (left_source,),
            right_source_cpbs = (right_source,),
        ),
    )
end

function _record_dispatch_pqs_record()
    return _record_dispatch_record(
        :pqs_source_pair_preflight;
        pair_key = (:pqs_left, :pqs_right),
        pair_family = :pqs_pqs,
        metadata = (;
            left_source_mode_dims = (2, 2, 1),
            right_source_mode_dims = (3, 1, 2),
            left_source_mode_count = 4,
            right_source_mode_count = 6,
            source_mode_ordering = :x_major_y_major_z_fast,
            left_source_mode_ordering = :x_major_y_major_z_fast,
            right_source_mode_ordering = :x_major_y_major_z_fast,
        ),
    )
end

function _record_dispatch_record(
    materialization_path::Symbol;
    pair_key = (:left_unit, :right_unit),
    pair_family = :synthetic_pair_family,
    metadata = (;),
)
    return CPBMRecordDispatch.PairBlockMaterializationRecord(
        pair_key,
        1,
        pair_family,
        :synthetic_source_operator_path,
        (; left = :synthetic_left_transform, right = :synthetic_right_transform),
        (; left = :synthetic_left_realization, right = :synthetic_right_realization),
        :synthetic_final_block_path,
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
        materialization_path,
        :ready_metadata_only_not_materialized,
        nothing,
        false,
        NamedTuple(metadata),
    )
end

function _record_dispatch_overlap_1d()
    return (;
        x = [
            1.0 0.1 0.2 0.3
            0.1 1.1 0.4 0.5
            0.2 0.4 1.2 0.6
            0.3 0.5 0.6 1.3
        ],
        y = [
            1.4 0.7 0.8 0.9
            0.7 1.5 1.0 1.1
            0.8 1.0 1.6 1.2
            0.9 1.1 1.2 1.7
        ],
        z = [
            1.8 1.3 1.4 1.5
            1.3 1.9 1.6 1.7
            1.4 1.6 2.0 1.8
            1.5 1.7 1.8 2.1
        ],
    )
end

function _record_dispatch_position_1d()
    return (;
        x = [
            2.0 2.1 2.2 2.3
            2.1 2.4 2.5 2.6
            2.2 2.5 2.7 2.8
            2.3 2.6 2.8 2.9
        ],
        y = [
            3.0 3.1 3.2 3.3
            3.1 3.4 3.5 3.6
            3.2 3.5 3.7 3.8
            3.3 3.6 3.8 3.9
        ],
        z = [
            4.0 4.1 4.2 4.3
            4.1 4.4 4.5 4.6
            4.2 4.5 4.7 4.8
            4.3 4.6 4.8 4.9
        ],
    )
end

function _record_dispatch_pqs_overlap_1d()
    return (;
        x = [
            1.0 0.1 0.2
            0.3 1.1 0.4
        ],
        y = reshape([1.2, 1.3], 2, 1),
        z = reshape([1.4, 1.5], 1, 2),
    )
end

function _record_dispatch_pqs_position_1d()
    return (;
        x = [
            2.0 2.1 2.2
            2.3 2.4 2.5
        ],
        y = reshape([3.0, 3.1], 2, 1),
        z = reshape([4.0, 4.1], 1, 2),
    )
end

function _check_record_dispatch_skipped_flags(summary)
    @test summary.status == :skipped_mixed_one_body_pair_block
    @test !summary.materialized
    @test !summary.factors_constructed
    @test !summary.numerical_blocks_materialized
    @test !summary.source_operator_blocks_materialized
    @test !summary.final_pair_blocks_materialized
    @test !summary.operator_blocks_materialized
    @test !summary.hamiltonian_data_materialized
    @test !summary.artifacts_materialized
    @test !summary.mixed_dispatcher_materialized
    @test !summary.route_driver_wiring
    @test !summary.pqs_lowdin_materialized
    @test !summary.full_white_lindsey_route_assembled
end

@testset "CartesianPairBlockMaterialization mixed direct record dispatch" begin
    record = _record_dispatch_direct_record()
    parent_axis_counts = (4, 4, 4)
    overlap_1d = _record_dispatch_overlap_1d()
    position_1d = _record_dispatch_position_1d()

    overlap = CPBMRecordDispatch._one_body_pair_block(
        record,
        :overlap;
        inputs = (; parent_axis_counts, overlap_1d),
    )
    expected_overlap = CPBMRecordDispatch.direct_direct_one_body_block(
        record,
        :overlap;
        parent_axis_counts,
        overlap_1d,
    )
    @test overlap isa CPBMRecordDispatch.PairBlockMaterializationResult
    @test overlap.term == :overlap
    @test overlap.pair_key == (:direct_left, :direct_right)
    @test overlap.materialized
    @test overlap.source_operator_blocks_materialized
    @test overlap.final_pair_blocks_materialized
    @test !overlap.operator_blocks_materialized
    @test overlap.block == expected_overlap.block
    @test overlap.metadata.mixed_one_body_dispatcher ==
          :direct_direct_only_record_dispatcher
    @test overlap.metadata.selector_family == :direct_direct
    @test overlap.metadata.numerical_family_selector ==
          :direct_direct_one_body_block

    position_y = CPBMRecordDispatch._one_body_pair_block(
        record,
        :position_y;
        inputs = (; parent_axis_counts, overlap_1d, position_1d),
    )
    expected_position_y = CPBMRecordDispatch.direct_direct_one_body_block(
        record,
        :position_y;
        parent_axis_counts,
        overlap_1d,
        position_1d,
    )
    @test position_y isa CPBMRecordDispatch.PairBlockMaterializationResult
    @test position_y.term == :position_y
    @test position_y.block == expected_position_y.block
    @test position_y.metadata.position_axis == :y
    @test position_y.metadata.selector_family == :direct_direct
end

@testset "CartesianPairBlockMaterialization mixed PQS source record dispatch" begin
    record = _record_dispatch_pqs_record()
    overlap_1d = _record_dispatch_pqs_overlap_1d()
    position_1d = _record_dispatch_pqs_position_1d()

    overlap = CPBMRecordDispatch._one_body_pair_block(
        record,
        :overlap;
        inputs = (; overlap_1d),
    )
    expected_overlap = CPBMRecordDispatch.pqs_source_pair_one_body_block(
        record,
        :overlap;
        overlap_1d,
    )
    @test overlap isa CPBMRecordDispatch.PairBlockMaterializationResult
    @test overlap.term == :source_overlap
    @test overlap.pair_key == (:pqs_left, :pqs_right)
    @test overlap.materialized
    @test overlap.source_operator_blocks_materialized
    @test !overlap.final_pair_blocks_materialized
    @test !overlap.operator_blocks_materialized
    @test !overlap.hamiltonian_data_materialized
    @test !overlap.artifacts_materialized
    @test overlap.block == expected_overlap.block
    @test overlap.metadata.block_space == :raw_product_source_modes
    @test overlap.metadata.source_mode_ordering == :x_major_y_major_z_fast
    @test overlap.metadata.left_source_mode_dims == (2, 2, 1)
    @test overlap.metadata.right_source_mode_dims == (3, 1, 2)
    @test overlap.metadata.left_source_mode_count == 4
    @test overlap.metadata.right_source_mode_count == 6
    @test !overlap.metadata.final_pair_blocks_materialized
    @test !overlap.metadata.shell_realization_materialized
    @test !overlap.metadata.operator_blocks_materialized
    @test !overlap.metadata.hamiltonian_data_materialized
    @test !overlap.metadata.artifacts_materialized
    @test overlap.metadata.mixed_one_body_dispatcher ==
          :pqs_source_pair_record_dispatcher
    @test overlap.metadata.selector_family == :pqs_source_pair
    @test overlap.metadata.numerical_family_selector ==
          :pqs_source_pair_one_body_block

    position_y = CPBMRecordDispatch._one_body_pair_block(
        record,
        :position_y;
        inputs = (; overlap_1d, position_1d),
    )
    expected_position_y = CPBMRecordDispatch.pqs_source_pair_one_body_block(
        record,
        :position_y;
        overlap_1d,
        position_1d,
    )
    @test position_y isa CPBMRecordDispatch.PairBlockMaterializationResult
    @test position_y.term == :source_position_y
    @test position_y.block == expected_position_y.block
    @test position_y.metadata.position_axis == :y
    @test !position_y.final_pair_blocks_materialized
    @test !position_y.metadata.shell_realization_materialized
end

@testset "CartesianPairBlockMaterialization mixed record dispatch skipped cases" begin
    direct_record = _record_dispatch_direct_record()
    overlap_1d = _record_dispatch_overlap_1d()
    missing_inputs = CPBMRecordDispatch._one_body_pair_block(
        direct_record,
        :kinetic;
        inputs = (; overlap_1d),
    )
    @test missing_inputs.blocker == :missing_one_body_factor_inputs
    @test missing_inputs.selector_family == :direct_direct
    @test missing_inputs.factor_input_blockers ==
          (:missing_required_one_body_factors, :missing_parent_axis_counts)
    _check_record_dispatch_skipped_flags(missing_inputs)

    pqs_record = _record_dispatch_record(
        :pqs_source_pair_preflight;
        pair_key = (:pqs_left, :pqs_right),
        pair_family = :pqs_pqs,
    )
    pqs_skipped = CPBMRecordDispatch._one_body_pair_block(
        pqs_record,
        :overlap;
        inputs = (; overlap_1d),
    )
    @test pqs_skipped.blocker == :missing_pqs_source_mode_metadata
    @test pqs_skipped.selector_family == :pqs_source_pair
    @test pqs_skipped.dispatch_status == :ready_metadata_only_not_materialized
    _check_record_dispatch_skipped_flags(pqs_skipped)

    lw_record = _record_dispatch_record(
        :white_lindsey_boundary_stratum_adapter_preflight;
        pair_key = (:lw_left, :lw_right),
        pair_family = :white_lindsey_boundary_stratum,
    )
    lw_skipped = CPBMRecordDispatch._one_body_pair_block(
        lw_record,
        :overlap;
        inputs = (; parent_axis_counts = (4, 4, 4), overlap_1d),
        unit_pair = (; pair_key = (:lw_left, :lw_right)),
    )
    @test lw_skipped.blocker ==
          :mixed_one_body_selector_family_not_materialized_yet
    @test lw_skipped.selector_family == :white_lindsey_boundary_stratum
    @test lw_skipped.unit_pair_requirement_status ==
          :available_aligned_white_lindsey_unit_pair
    _check_record_dispatch_skipped_flags(lw_skipped)
end

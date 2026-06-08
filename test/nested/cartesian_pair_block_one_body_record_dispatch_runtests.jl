# Runtime role: contract.
#
# Record-level numerical mixed dispatch contract across all safe one-body terms
# for direct/direct and PQS source-space families. Use for selector semantics;
# the tiny consumer smoke is preferred for routine mixed-consumer edits.

using Test
using GaussletBases

const CPBMRecordDispatch = GaussletBases.CartesianPairBlockMaterialization
const CPBRecordDispatch = GaussletBases.CartesianCPB

const _RECORD_DISPATCH_SAFE_TERMS = (
    :overlap,
    :position_x,
    :position_y,
    :position_z,
    :x2_x,
    :x2_y,
    :x2_z,
    :kinetic,
)

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

function _record_dispatch_x2_1d()
    return (;
        x = [
            5.0 5.1 5.2 5.3
            5.1 5.4 5.5 5.6
            5.2 5.5 5.7 5.8
            5.3 5.6 5.8 5.9
        ],
        y = [
            6.0 6.1 6.2 6.3
            6.1 6.4 6.5 6.6
            6.2 6.5 6.7 6.8
            6.3 6.6 6.8 6.9
        ],
        z = [
            7.0 7.1 7.2 7.3
            7.1 7.4 7.5 7.6
            7.2 7.5 7.7 7.8
            7.3 7.6 7.8 7.9
        ],
    )
end

function _record_dispatch_kinetic_1d()
    return (;
        x = [
            8.0 8.1 8.2 8.3
            8.1 8.4 8.5 8.6
            8.2 8.5 8.7 8.8
            8.3 8.6 8.8 8.9
        ],
        y = [
            9.0 9.1 9.2 9.3
            9.1 9.4 9.5 9.6
            9.2 9.5 9.7 9.8
            9.3 9.6 9.8 9.9
        ],
        z = [
            10.0 10.1 10.2 10.3
            10.1 10.4 10.5 10.6
            10.2 10.5 10.7 10.8
            10.3 10.6 10.8 10.9
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

function _record_dispatch_pqs_x2_1d()
    return (;
        x = [
            5.0 5.1 5.2
            5.3 5.4 5.5
        ],
        y = reshape([6.0, 6.1], 2, 1),
        z = reshape([7.0, 7.1], 1, 2),
    )
end

function _record_dispatch_pqs_kinetic_1d()
    return (;
        x = [
            8.0 8.1 8.2
            8.3 8.4 8.5
        ],
        y = reshape([9.0, 9.1], 2, 1),
        z = reshape([10.0, 10.1], 1, 2),
    )
end

function _check_record_dispatch_materialized_flags(result)
    @test result isa CPBMRecordDispatch.PairBlockMaterializationResult
    @test result.materialized
    @test !result.operator_blocks_materialized
    @test !result.hamiltonian_data_materialized
    @test !result.artifacts_materialized
end

function _check_record_dispatch_direct_result(result, expected)
    _check_record_dispatch_materialized_flags(result)
    @test result.term == expected.term
    @test result.pair_key == (:direct_left, :direct_right)
    @test result.block == expected.block
    @test result.source_operator_blocks_materialized
    @test result.final_pair_blocks_materialized
    @test result.metadata.mixed_one_body_dispatcher ==
          :direct_direct_only_record_dispatcher
    @test result.metadata.selector_family == :direct_direct
    @test result.metadata.numerical_family_selector ==
          :direct_direct_one_body_block
end

function _check_record_dispatch_pqs_result(result, expected)
    _check_record_dispatch_materialized_flags(result)
    @test result.term == expected.term
    @test result.pair_key == (:pqs_left, :pqs_right)
    @test result.block == expected.block
    @test result.source_operator_blocks_materialized
    @test !result.final_pair_blocks_materialized
    @test result.metadata.block_space == :raw_product_source_modes
    @test result.metadata.source_mode_ordering == :x_major_y_major_z_fast
    @test result.metadata.left_source_mode_dims == (2, 2, 1)
    @test result.metadata.right_source_mode_dims == (3, 1, 2)
    @test result.metadata.left_source_mode_count == 4
    @test result.metadata.right_source_mode_count == 6
    @test !result.metadata.final_pair_blocks_materialized
    @test !result.metadata.shell_realization_materialized
    @test !result.metadata.operator_blocks_materialized
    @test !result.metadata.hamiltonian_data_materialized
    @test !result.metadata.artifacts_materialized
    @test result.metadata.mixed_one_body_dispatcher ==
          :pqs_source_pair_record_dispatcher
    @test result.metadata.selector_family == :pqs_source_pair
    @test result.metadata.numerical_family_selector ==
          :pqs_source_pair_one_body_block
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
    x2_1d = _record_dispatch_x2_1d()
    kinetic_1d = _record_dispatch_kinetic_1d()
    inputs = (; parent_axis_counts, overlap_1d, position_1d, x2_1d, kinetic_1d)

    for term in _RECORD_DISPATCH_SAFE_TERMS
        result = CPBMRecordDispatch._one_body_pair_block(record, term; inputs)
        expected = CPBMRecordDispatch.direct_direct_one_body_block(
            record,
            term;
            parent_axis_counts,
            overlap_1d,
            position_1d,
            x2_1d,
            kinetic_1d,
        )
        _check_record_dispatch_direct_result(result, expected)
    end
end

@testset "CartesianPairBlockMaterialization mixed PQS source record dispatch" begin
    record = _record_dispatch_pqs_record()
    overlap_1d = _record_dispatch_pqs_overlap_1d()
    position_1d = _record_dispatch_pqs_position_1d()
    x2_1d = _record_dispatch_pqs_x2_1d()
    kinetic_1d = _record_dispatch_pqs_kinetic_1d()
    inputs = (; overlap_1d, position_1d, x2_1d, kinetic_1d)

    for term in _RECORD_DISPATCH_SAFE_TERMS
        result = CPBMRecordDispatch._one_body_pair_block(record, term; inputs)
        expected = CPBMRecordDispatch.pqs_source_pair_one_body_block(
            record,
            term;
            overlap_1d,
            position_1d,
            x2_1d,
            kinetic_1d,
        )
        _check_record_dispatch_pqs_result(result, expected)
    end
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
        materialize_selector_families = (:direct_direct, :pqs_source_pair),
    )
    @test lw_skipped.blocker ==
          :mixed_one_body_selector_family_not_materialized_yet
    @test lw_skipped.selector_family == :white_lindsey_boundary_stratum
    @test lw_skipped.unit_pair_requirement_status ==
          :available_aligned_white_lindsey_unit_pair
    _check_record_dispatch_skipped_flags(lw_skipped)
end

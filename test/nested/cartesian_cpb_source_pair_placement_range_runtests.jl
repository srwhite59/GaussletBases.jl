# Runtime role: metadata-only CPB source-pair placement range contract.
#
# This test validates retained column range and dimension metadata only. It does
# not apply transforms, place matrices, or assemble route/global overlap.

using Test
using GaussletBases

const CPBRangeCarry = GaussletBases.CartesianCPBBlockProviders

function _range_transform_carry(;
    side = :left,
    target_retained_column_range = 1:2,
    transform_object = [1.0 0.0; 0.0 1.0],
)
    source_shape = (x = 1, y = 2, z = 1)
    return CPBRangeCarry.cpb_retained_transform_carry(
        side,
        (:tiny_left, :tiny_right),
        (; object_kind = :test_cpb_source_summary, shape = source_shape),
        source_shape,
        :parent_compatible_x_slowest_z_fastest,
        target_retained_column_range;
        transform_object,
        transform_convention = :test_local_to_retained_columns,
        transform_provenance = :test_fixture,
    )
end

function _placement_range(; kwargs...)
    return CPBRangeCarry.cpb_source_pair_placement_range(
        (:tiny_left, :tiny_right);
        left_column_range = 1:2,
        right_column_range = 3:4,
        global_dimension = 4,
        global_dimension_source = :test_retained_layout,
        range_source = :test_source_pair_ranges,
        range_provenance = :test_fixture,
        kwargs...,
    )
end

@testset "CPB source pair placement range metadata" begin
    missing_left = _placement_range(; left_column_range = nothing)
    missing_left_summary = CPBRangeCarry.summary(missing_left)
    @test missing_left_summary.object_kind ===
          :cartesian_cpb_source_pair_placement_range_summary
    @test missing_left_summary.status ===
          :blocked_cpb_source_pair_placement_range
    @test missing_left_summary.blocker === :missing_left_column_range
    @test missing_left_summary.block_key === (:tiny_left, :tiny_right)
    @test missing_left_summary.left_column_count === :unavailable
    @test missing_left_summary.right_column_count == 2
    @test missing_left_summary.global_dimension == 4
    @test missing_left_summary.global_dimension_source === :test_retained_layout
    @test missing_left_summary.range_source === :test_source_pair_ranges
    @test missing_left_summary.range_provenance === :test_fixture
    @test missing_left_summary.placement_status === :unassigned
    @test missing_left_summary.global_matrix_materialized === false
    @test missing_left_summary.route_driver_wiring === false

    missing_right = _placement_range(; right_column_range = nothing)
    @test CPBRangeCarry.summary(missing_right).status ===
          :blocked_cpb_source_pair_placement_range
    @test CPBRangeCarry.summary(missing_right).blocker ===
          :missing_right_column_range

    missing_dimension = _placement_range(; global_dimension = nothing)
    @test CPBRangeCarry.summary(missing_dimension).status ===
          :blocked_cpb_source_pair_placement_range
    @test CPBRangeCarry.summary(missing_dimension).blocker ===
          :missing_global_dimension

    available = _placement_range()
    available_summary = CPBRangeCarry.summary(available)
    @test available_summary.status === :available_cpb_source_pair_placement_range
    @test available_summary.blocker === nothing
    @test available_summary.left_column_range == 1:2
    @test available_summary.right_column_range == 3:4
    @test available_summary.left_column_count == 2
    @test available_summary.right_column_count == 2
    @test available_summary.left_transform_status === :unavailable
    @test available_summary.right_transform_status === :unavailable
    @test available_summary.left_transform_target_count === :unavailable
    @test available_summary.right_transform_target_count === :unavailable
    @test available_summary.global_matrix_materialized === false
    @test available_summary.route_driver_wiring === false

    left_outside = _placement_range(; left_column_range = 0:1)
    @test CPBRangeCarry.summary(left_outside).status ===
          :blocked_cpb_source_pair_placement_range
    @test CPBRangeCarry.summary(left_outside).blocker ===
          :left_column_range_dimension_mismatch

    right_outside = _placement_range(; right_column_range = 4:5)
    @test CPBRangeCarry.summary(right_outside).status ===
          :blocked_cpb_source_pair_placement_range
    @test CPBRangeCarry.summary(right_outside).blocker ===
          :right_column_range_dimension_mismatch

    left_transform_mismatch = _placement_range(;
        left_column_range = 1:3,
        global_dimension = 5,
        left_transform_carry = _range_transform_carry(;
            side = :left,
            target_retained_column_range = 1:2,
        ),
    )
    @test CPBRangeCarry.summary(left_transform_mismatch).status ===
          :blocked_cpb_source_pair_placement_range
    @test CPBRangeCarry.summary(left_transform_mismatch).blocker ===
          :left_column_range_dimension_mismatch
    @test CPBRangeCarry.summary(left_transform_mismatch).left_transform_status ===
          :available_cpb_retained_transform_carry
    @test CPBRangeCarry.summary(left_transform_mismatch).left_transform_target_count ==
          2

    right_transform_mismatch = _placement_range(;
        right_column_range = 3:5,
        global_dimension = 5,
        right_transform_carry = _range_transform_carry(;
            side = :right,
            target_retained_column_range = 1:2,
        ),
    )
    @test CPBRangeCarry.summary(right_transform_mismatch).status ===
          :blocked_cpb_source_pair_placement_range
    @test CPBRangeCarry.summary(right_transform_mismatch).blocker ===
          :right_column_range_dimension_mismatch
    @test CPBRangeCarry.summary(right_transform_mismatch).right_transform_status ===
          :available_cpb_retained_transform_carry
    @test CPBRangeCarry.summary(right_transform_mismatch).right_transform_target_count ==
          2

    with_transform_carries = _placement_range(;
        left_transform_carry = _range_transform_carry(; side = :left),
        right_transform_carry = _range_transform_carry(; side = :right),
    )
    with_transform_summary = CPBRangeCarry.summary(with_transform_carries)
    @test with_transform_summary.status ===
          :available_cpb_source_pair_placement_range
    @test with_transform_summary.left_transform_status ===
          :available_cpb_retained_transform_carry
    @test with_transform_summary.right_transform_status ===
          :available_cpb_retained_transform_carry
    @test with_transform_summary.left_transform_target_count == 2
    @test with_transform_summary.right_transform_target_count == 2
    @test with_transform_summary.placement_status === :unassigned
    @test with_transform_summary.global_matrix_materialized === false
    @test with_transform_summary.route_driver_wiring === false
end

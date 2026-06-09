# Runtime role: metadata-only CPB retained transform carry contract.
#
# This test validates transform carry shape/count/order metadata only. It does
# not apply retained transforms or place local CPB overlap blocks globally.

using Test
using GaussletBases

const CPBTransformCarry = GaussletBases.CartesianCPBBlockProviders

function _transform_carry_source_summary(; shape = (x = 1, y = 2, z = 1))
    return (;
        object_kind = :test_cpb_source_summary,
        shape,
        support_count = shape.x * shape.y * shape.z,
    )
end

function _transform_carry(;
    source_shape = (x = 1, y = 2, z = 1),
    source_ordering = :parent_compatible_x_slowest_z_fastest,
    target_retained_column_range = 1:2,
    transform_object = [1.0 0.0; 0.0 1.0],
)
    return CPBTransformCarry.cpb_retained_transform_carry(
        :left,
        (:tiny_left, :tiny_right),
        _transform_carry_source_summary(; shape = source_shape),
        source_shape,
        source_ordering,
        target_retained_column_range;
        transform_object,
        transform_convention = :test_local_to_retained_columns,
        transform_provenance = :test_fixture,
    )
end

@testset "CPB retained transform carry metadata" begin
    missing_transform = _transform_carry(; transform_object = nothing)
    missing_summary = CPBTransformCarry.summary(missing_transform)
    @test missing_summary.object_kind ===
          :cartesian_cpb_retained_transform_carry_summary
    @test missing_summary.status === :blocked_cpb_retained_transform_carry
    @test missing_summary.blocker === :missing_retained_transform
    @test missing_summary.side === :left
    @test missing_summary.block_key === (:tiny_left, :tiny_right)
    @test missing_summary.source_shape == (x = 1, y = 2, z = 1)
    @test missing_summary.source_support_count == 2
    @test missing_summary.target_retained_column_range == 1:2
    @test missing_summary.target_retained_column_count == 2
    @test missing_summary.transform_available === false
    @test missing_summary.transform_shape === :unavailable
    @test missing_summary.transform_reference_kind === :unavailable
    @test missing_summary.placement_status === :unassigned
    @test missing_summary.transform_applied === false
    @test missing_summary.global_matrix_materialized === false
    @test missing_summary.route_driver_wiring === false

    available = _transform_carry()
    available_summary = CPBTransformCarry.summary(available)
    @test available_summary.status === :available_cpb_retained_transform_carry
    @test available_summary.blocker === nothing
    @test available_summary.transform_available
    @test available_summary.transform_shape == (2, 2)
    @test available_summary.transform_reference_kind === :matrix
    @test available_summary.transform_convention ===
          :test_local_to_retained_columns
    @test available_summary.transform_provenance === :test_fixture
    @test available_summary.transform_applied === false
    @test available_summary.global_matrix_materialized === false
    @test available_summary.route_driver_wiring === false

    source_shape_mismatch = _transform_carry(;
        source_shape = (x = 1, y = 3, z = 1),
        transform_object = [1.0 0.0; 0.0 1.0],
    )
    @test CPBTransformCarry.summary(source_shape_mismatch).status ===
          :blocked_cpb_retained_transform_carry
    @test CPBTransformCarry.summary(source_shape_mismatch).blocker ===
          :retained_transform_source_shape_mismatch

    target_count_mismatch = _transform_carry(;
        target_retained_column_range = 1:3,
        transform_object = [1.0 0.0; 0.0 1.0],
    )
    @test CPBTransformCarry.summary(target_count_mismatch).status ===
          :blocked_cpb_retained_transform_carry
    @test CPBTransformCarry.summary(target_count_mismatch).blocker ===
          :retained_transform_target_count_mismatch

    ordering_mismatch = _transform_carry(;
        source_ordering = :z_slowest_x_fastest,
        transform_object = [1.0 0.0; 0.0 1.0],
    )
    @test CPBTransformCarry.summary(ordering_mismatch).status ===
          :blocked_cpb_retained_transform_carry
    @test CPBTransformCarry.summary(ordering_mismatch).blocker ===
          :retained_transform_ordering_mismatch
    @test CPBTransformCarry.summary(ordering_mismatch).transform_applied === false
    @test CPBTransformCarry.summary(ordering_mismatch).global_matrix_materialized ===
          false
    @test CPBTransformCarry.summary(ordering_mismatch).route_driver_wiring === false
end

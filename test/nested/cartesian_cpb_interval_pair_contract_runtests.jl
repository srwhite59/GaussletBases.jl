# Runtime role: tiny CPB interval-pair contract test.
#
# This validates parent-box CPB pair metadata for future matrix kernels. It
# does not consume parent factor packets, slice operator matrices, wire route
# drivers, or exercise Hamiltonian/IDA/MWG/PQS export paths.

using Test
using GaussletBases

const CPBPairTest = GaussletBases.CartesianCPB
const CPGBPairTest = GaussletBases.CartesianParentGaussletBases
const CKPairTest = GaussletBases.CartesianCPBBlockProviders

function _interval_pair_test_parent()
    axis_x = build_basis(MappedUniformBasisSpec(
        :G10;
        count = 4,
        mapping = IdentityMapping(),
        reference_spacing = 1.0,
    ))
    axis_y = build_basis(MappedUniformBasisSpec(
        :G10;
        count = 5,
        mapping = IdentityMapping(),
        reference_spacing = 1.0,
    ))
    axis_z = build_basis(MappedUniformBasisSpec(
        :G10;
        count = 6,
        mapping = IdentityMapping(),
        reference_spacing = 1.0,
    ))
    return CPGBPairTest.CartesianParentGaussletBasis3D(axis_x, axis_y, axis_z)
end

@testset "Cartesian CPB interval pair contract" begin
    parent = _interval_pair_test_parent()
    left_filled = CPBPairTest.filled_cpb(
        1:2,
        2:4,
        3:6;
        role = :left_filled,
    )
    right_filled = CPBPairTest.filled_cpb(
        2:4,
        1:3,
        1:2;
        role = :right_filled,
    )
    filled_pair =
        CKPairTest.cpb_interval_pair(parent, left_filled, right_filled)
    filled_summary = CKPairTest.summary(filled_pair)

    @test filled_pair.parent === parent
    @test filled_pair.left_cpb === left_filled
    @test filled_pair.right_cpb === right_filled
    @test filled_summary.status === :available_cpb_interval_pair
    @test filled_summary.blocker === nothing
    @test filled_summary.parent_axis_counts == (4, 5, 6)
    @test filled_summary.parent_intervals == (x = 1:4, y = 1:5, z = 1:6)
    @test filled_summary.left_intervals == (x = 1:2, y = 2:4, z = 3:6)
    @test filled_summary.right_intervals == (x = 2:4, y = 1:3, z = 1:2)
    @test filled_summary.left_shape == (x = 2, y = 3, z = 4)
    @test filled_summary.right_shape == (x = 3, y = 3, z = 2)
    @test filled_summary.left_support_count ==
          CPBPairTest.support_count(left_filled)
    @test filled_summary.right_support_count ==
          CPBPairTest.support_count(right_filled)
    @test filled_summary.left_support_count ==
          filled_summary.left_shape.x *
          filled_summary.left_shape.y *
          filled_summary.left_shape.z
    @test filled_summary.axis_order === (:x, :y, :z)
    @test filled_summary.local_ordering ===
          :parent_compatible_x_slowest_z_fastest
    @test isempty(filled_summary.outside_parent_intervals)
    @test filled_summary.operator_matrices_sliced === false
    @test filled_summary.parent_factor_packet_consumed === false
    @test filled_summary.route_driver_wiring === false

    left_slab = CPBPairTest.slab_cpb(
        1:4,
        3:3,
        2:6;
        role = :left_y_slab,
    )
    right_slab = CPBPairTest.slab_cpb(
        2:2,
        1:5,
        1:6;
        role = :right_x_slab,
    )
    slab_pair = CKPairTest.cpb_interval_pair(parent, left_slab, right_slab)
    slab_summary = CKPairTest.summary(slab_pair)

    @test slab_summary.status === :available_cpb_interval_pair
    @test slab_summary.blocker === nothing
    @test slab_summary.left_intervals == (x = 1:4, y = 3:3, z = 2:6)
    @test slab_summary.right_intervals == (x = 2:2, y = 1:5, z = 1:6)
    @test slab_summary.left_shape == (x = 4, y = 1, z = 5)
    @test slab_summary.right_shape == (x = 1, y = 5, z = 6)
    @test slab_summary.left_codimension == 1
    @test slab_summary.right_codimension == 1
    @test slab_summary.left_support_count ==
          CPBPairTest.support_count(left_slab)
    @test slab_summary.right_support_count ==
          CPBPairTest.support_count(right_slab)
    @test slab_summary.axis_order === (:x, :y, :z)
    @test slab_summary.local_ordering ===
          :parent_compatible_x_slowest_z_fastest

    outside_left = CPBPairTest.cpb(
        0:1,
        2:4,
        3:6;
        role = :outside_left_x,
    )
    outside_right = CPBPairTest.cpb(
        2:4,
        1:3,
        1:7;
        role = :outside_right_z,
    )
    blocked_pair =
        CKPairTest.cpb_interval_pair(parent, outside_left, outside_right)
    blocked_summary = CKPairTest.summary(blocked_pair)

    @test blocked_summary.status === :blocked_cpb_interval_pair
    @test blocked_summary.blocker === :left_x_interval_outside_parent
    @test blocked_summary.parent_axis_counts == (4, 5, 6)
    @test blocked_summary.left_intervals.x == 0:1
    @test blocked_summary.right_intervals.z == 1:7
    @test length(blocked_summary.outside_parent_intervals) == 2
    @test blocked_summary.outside_parent_intervals[1] == (;
        side = :left,
        axis = :x,
        interval = 0:1,
        parent_interval = 1:4,
    )
    @test blocked_summary.outside_parent_intervals[2] == (;
        side = :right,
        axis = :z,
        interval = 1:7,
        parent_interval = 1:6,
    )
    @test blocked_summary.operator_matrices_sliced === false
    @test blocked_summary.parent_factor_packet_consumed === false
    @test blocked_summary.route_driver_wiring === false
end

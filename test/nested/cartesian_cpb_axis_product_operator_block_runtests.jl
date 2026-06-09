# Runtime role: CPB-local axis-product operator block contract.
#
# This test validates the generic CPB-local axis-product dense block primitive
# and the overlap wrapper around it. It does not add WL/PQS realization,
# retained transforms, route/global placement, driver wiring, kinetic,
# position, x2, Coulomb, Hamiltonian, IDA/MWG, PQS Lowdin/projection, exports,
# or artifacts.

using Test
using GaussletBases

const CPBOperatorBlock = GaussletBases.CartesianCPB
const CPGBOperatorBlock = GaussletBases.CartesianParentGaussletBases
const CBPOperatorBlock = GaussletBases.CartesianCPBBlockProviders

function _axis_product_expected_dense(axis_ops)
    nx_left, nx_right = size(axis_ops.x)
    ny_left, ny_right = size(axis_ops.y)
    nz_left, nz_right = size(axis_ops.z)
    dense = Matrix{Float64}(
        undef,
        nx_left * ny_left * nz_left,
        nx_right * ny_right * nz_right,
    )
    for ix_left in 1:nx_left, iy_left in 1:ny_left, iz_left in 1:nz_left
        left_index =
            (ix_left - 1) * ny_left * nz_left +
            (iy_left - 1) * nz_left +
            iz_left
        for ix_right in 1:nx_right, iy_right in 1:ny_right, iz_right in 1:nz_right
            right_index =
                (ix_right - 1) * ny_right * nz_right +
                (iy_right - 1) * nz_right +
                iz_right
            dense[left_index, right_index] =
                axis_ops.x[ix_left, ix_right] *
                axis_ops.y[iy_left, iy_right] *
                axis_ops.z[iz_left, iz_right]
        end
    end
    return dense
end

function _axis_product_case(axis_ops)
    block = CBPOperatorBlock.cpb_axis_product_operator_block(
        axis_ops;
        term = :test_axis_product,
        factor_space = :test_factor_space,
        factor_convention = :test_factor_convention,
        index_domain = :parent_axis_indices,
    )
    block_summary = CBPOperatorBlock.summary(block)
    @test block_summary.status === :materialized_cpb_axis_product_operator_block
    @test block_summary.blocker === nothing
    @test block_summary.term === :test_axis_product
    @test block_summary.representation === :dense_local_cpb_product_space
    @test block_summary.local_ordering ===
          :parent_compatible_x_slowest_z_fastest
    @test block_summary.left_shape == (
        x = size(axis_ops.x, 1),
        y = size(axis_ops.y, 1),
        z = size(axis_ops.z, 1),
    )
    @test block_summary.right_shape == (
        x = size(axis_ops.x, 2),
        y = size(axis_ops.y, 2),
        z = size(axis_ops.z, 2),
    )
    @test block_summary.dense_block_shape ==
          (block_summary.left_support_count, block_summary.right_support_count)
    @test block_summary.dense_block_eltype === Float64
    @test block_summary.provider_level_local_matrix_materialized === true
    @test block_summary.realization_status === :unrealized
    @test block_summary.route_global_status === :unassigned
    @test block_summary.route_driver_wiring === false
    @test block_summary.route_global_matrix_materialized === false
    @test block_summary.global_matrix_materialized === false
    @test !hasproperty(block_summary, :dense_block)
    @test block.dense_block == _axis_product_expected_dense(axis_ops)
    return block
end

function _operator_block_parent(; count = 3)
    axis = build_basis(MappedUniformBasisSpec(
        :G10;
        count,
        mapping = IdentityMapping(),
        reference_spacing = 1.0,
    ))
    return CPGBOperatorBlock.CartesianParentGaussletBasis3D(axis)
end

function _operator_block_overlap_1d()
    return (;
        x = [
            1.0 0.1 0.2
            0.3 1.1 0.4
            0.5 0.6 1.2
        ],
        y = [
            2.0 0.7 0.8
            0.9 2.1 1.0
            1.1 1.2 2.2
        ],
        z = [
            3.0 1.3 1.4
            1.5 3.1 1.6
            1.7 1.8 3.2
        ],
    )
end

function _operator_block_axis_bundle(overlap_1d)
    return (;
        x = (; pgdg_intermediate = (; overlap = overlap_1d.x)),
        y = (; pgdg_intermediate = (; overlap = overlap_1d.y)),
        z = (; pgdg_intermediate = (; overlap = overlap_1d.z)),
    )
end

@testset "CPB axis-product operator block" begin
    point_ops = (;
        x = fill(2.0, 1, 1),
        y = fill(3.0, 1, 1),
        z = fill(5.0, 1, 1),
    )
    point_block = _axis_product_case(point_ops)
    @test point_block.dense_block == fill(30.0, 1, 1)

    edge_ops = (;
        x = [1.0 2.0; 3.0 4.0],
        y = fill(5.0, 1, 1),
        z = fill(7.0, 1, 1),
    )
    edge_block = _axis_product_case(edge_ops)
    @test edge_block.dense_block == 35.0 .* edge_ops.x

    face_ops = (;
        x = reshape([1.0, 2.0], 2, 1),
        y = [3.0 4.0; 5.0 6.0],
        z = fill(7.0, 1, 1),
    )
    _axis_product_case(face_ops)

    cube_ops = (;
        x = [1.0 2.0; 3.0 4.0],
        y = reshape([5.0, 6.0], 2, 1),
        z = [7.0 8.0; 9.0 10.0],
    )
    cube_block = _axis_product_case(cube_ops)
    @test CBPOperatorBlock.summary(cube_block).left_shape == (x = 2, y = 2, z = 2)
    @test CBPOperatorBlock.summary(cube_block).right_shape == (x = 2, y = 1, z = 2)

    parent = _operator_block_parent()
    overlap_1d = _operator_block_overlap_1d()
    packet = CPGBOperatorBlock.parent_overlap_axis_factor_packet(
        parent,
        _operator_block_axis_bundle(overlap_1d),
    )
    left = CPBOperatorBlock.cpb(1:2, 2:3, 1:3; role = :operator_left)
    right = CPBOperatorBlock.cpb(2:3, 1:2, 2:2; role = :operator_right)
    interval_pair = CBPOperatorBlock.cpb_interval_pair(parent, left, right)
    overlap_operator =
        CBPOperatorBlock.cpb_overlap_operator_block(packet, interval_pair)
    overlap_summary = CBPOperatorBlock.summary(overlap_operator)
    axis_block_set = CBPOperatorBlock.cpb_overlap_axis_blocks(packet, interval_pair)
    existing_dense = CBPOperatorBlock.cpb_overlap_dense_block(axis_block_set)

    @test overlap_summary.status === :materialized_cpb_overlap_operator_block
    @test overlap_summary.blocker === nothing
    @test overlap_summary.term === :overlap
    @test overlap_summary.representation === :dense_local_cpb_product_space
    @test overlap_summary.realization_status === :unrealized
    @test overlap_summary.route_global_status === :unassigned
    @test overlap_summary.route_driver_wiring === false
    @test overlap_summary.route_global_matrix_materialized === false
    @test overlap_summary.provider_level_local_matrix_materialized === true
    @test overlap_summary.factor_convention === :axis_bundle_one_body_overlap
    @test !isnothing(overlap_operator.axis_product_block)
    @test overlap_operator.axis_product_block.dense_block == existing_dense.dense_block
    @test CBPOperatorBlock.summary(overlap_operator.axis_product_block).term === :overlap
    @test !hasproperty(overlap_summary, :dense_block)
    @test !hasproperty(overlap_summary, :global_overlap_matrix)
end

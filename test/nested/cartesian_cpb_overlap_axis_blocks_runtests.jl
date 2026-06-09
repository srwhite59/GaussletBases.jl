# Runtime role: tiny CPB overlap axis-block provider contract test.
#
# This validates axis-level overlap slices from a parent factor packet and a
# CPB interval pair, plus local dense CPB overlap materialization. It does not
# wire route/global placement or exercise kinetic/position/x2/Coulomb,
# Hamiltonian, IDA/MWG, PQS Lowdin/projection, exports, or artifacts.

using Test
using GaussletBases

const CPBOverlapBlocks = GaussletBases.CartesianCPB
const CPGBOverlapBlocks = GaussletBases.CartesianParentGaussletBases
const CBPOverlapBlocks = GaussletBases.CartesianCPBBlockProviders

function _overlap_blocks_parent(; count = 3)
    axis = build_basis(MappedUniformBasisSpec(
        :G10;
        count,
        mapping = IdentityMapping(),
        reference_spacing = 1.0,
    ))
    return CPGBOverlapBlocks.CartesianParentGaussletBasis3D(axis)
end

function _overlap_blocks_overlap_1d()
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

function _overlap_blocks_driver_fixture_overlap_1d()
    return (;
        x = [1.0 0.2; 0.2 1.1],
        y = [1.2 0.3; 0.3 1.3],
        z = [1.4 0.4; 0.4 1.5],
    )
end

function _overlap_blocks_driver_fixture_expected_matrix()
    overlap_1d = _overlap_blocks_driver_fixture_overlap_1d()
    scale = overlap_1d.x[1, 1] * overlap_1d.z[1, 1]
    return scale .* overlap_1d.y
end

function _overlap_blocks_axis_bundle(overlap_1d = _overlap_blocks_overlap_1d())
    return (;
        x = (; pgdg_intermediate = (; overlap = overlap_1d.x)),
        y = (; pgdg_intermediate = (; overlap = overlap_1d.y)),
        z = (; pgdg_intermediate = (; overlap = overlap_1d.z)),
    )
end

function _overlap_blocks_interval_pair(parent)
    left = CPBOverlapBlocks.cpb(1:2, 2:3, 1:3; role = :left_overlap_slice)
    right = CPBOverlapBlocks.cpb(2:3, 1:2, 2:2; role = :right_overlap_slice)
    return CBPOverlapBlocks.cpb_interval_pair(parent, left, right)
end

function _overlap_blocks_packet_with_metadata(packet, metadata)
    return CPGBOverlapBlocks.CartesianParentAxisFactorPacket3D(
        packet.parent,
        packet.overlap_1d,
        metadata,
    )
end

function _overlap_blocks_expected_dense(axis_blocks)
    nx_left, nx_right = size(axis_blocks.x)
    ny_left, ny_right = size(axis_blocks.y)
    nz_left, nz_right = size(axis_blocks.z)
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
                axis_blocks.x[ix_left, ix_right] *
                axis_blocks.y[iy_left, iy_right] *
                axis_blocks.z[iz_left, iz_right]
        end
    end
    return dense
end

@testset "Cartesian CPB overlap axis-block slices" begin
    parent = _overlap_blocks_parent()
    overlap_1d = _overlap_blocks_overlap_1d()
    packet = CPGBOverlapBlocks.parent_overlap_axis_factor_packet(
        parent,
        _overlap_blocks_axis_bundle(overlap_1d),
    )
    interval_pair = _overlap_blocks_interval_pair(parent)
    block_set = CBPOverlapBlocks.cpb_overlap_axis_blocks(packet, interval_pair)
    block_summary = CBPOverlapBlocks.summary(block_set)

    @test block_set.overlap_packet === packet
    @test block_set.interval_pair === interval_pair
    @test block_summary.status === :available_cpb_overlap_axis_blocks
    @test block_summary.blocker === nothing
    @test Matrix(block_set.axis_blocks.x) == overlap_1d.x[1:2, 2:3]
    @test Matrix(block_set.axis_blocks.y) == overlap_1d.y[2:3, 1:2]
    @test Matrix(block_set.axis_blocks.z) == overlap_1d.z[1:3, 2:2]
    @test block_summary.axis_block_shapes == (x = (2, 2), y = (2, 2), z = (3, 1))
    @test block_summary.factor_space === :parent_axis_bundle_pgdg_intermediate
    @test block_summary.factor_convention === :axis_bundle_one_body_overlap
    @test block_summary.normalization_convention ===
          :not_separate_from_axis_bundle_one_body_overlap
    @test block_summary.index_domain === :parent_axis_indices
    @test block_summary.index_domain_source === :axis_bundle_contract
    @test block_summary.index_domain_status ===
          :assumed_parent_axis_indexed_by_current_axis_bundle_contract
    @test block_summary.axis_order === (:x, :y, :z)
    @test block_summary.bra_ket_order === (:bra, :ket)
    @test block_summary.local_ordering ===
          :parent_compatible_x_slowest_z_fastest
    @test block_summary.blocks_are_views
    @test block_summary.blocks_are_copies === false
    @test block_summary.dense_local_block_materialized === false
    @test block_summary.parent_factor_packet_consumed
    @test block_summary.route_driver_wiring === false
    @test block_summary.hamiltonian_data_materialized === false
    @test block_summary.coulomb_data_materialized === false
    @test block_summary.ida_mwg_semantics === false
    @test block_summary.exports_or_artifacts === false
    @test block_summary.overlap_packet_summary.status ===
          :available_parent_overlap_axis_factors
    @test block_summary.interval_pair_summary.status ===
          :available_cpb_interval_pair

    dense_block = CBPOverlapBlocks.cpb_overlap_dense_block(block_set)
    dense_summary = CBPOverlapBlocks.summary(dense_block)

    @test dense_block.axis_block_set === block_set
    @test dense_summary.status === :materialized_cpb_overlap_dense_block
    @test dense_summary.blocker === nothing
    @test dense_summary.dense_block_available
    @test dense_block.dense_block ==
          _overlap_blocks_expected_dense(block_set.axis_blocks)
    @test dense_summary.dense_block_shape == (
        block_summary.interval_pair_summary.left_support_count,
        block_summary.interval_pair_summary.right_support_count,
    )
    @test dense_summary.dense_block_eltype === Float64
    @test dense_summary.left_shape == (x = 2, y = 2, z = 3)
    @test dense_summary.right_shape == (x = 2, y = 2, z = 1)
    @test dense_summary.local_ordering ===
          :parent_compatible_x_slowest_z_fastest
    @test dense_summary.source_axis_block_summary === block_summary
    @test dense_summary.factor_space ===
          :parent_axis_bundle_pgdg_intermediate
    @test dense_summary.factor_convention === :axis_bundle_one_body_overlap
    @test dense_summary.index_domain === :parent_axis_indices
    @test dense_summary.dense_local_block_materialized
    @test dense_summary.route_driver_wiring === false
    @test dense_summary.global_matrix_materialized === false
    @test dense_summary.hamiltonian_data_materialized === false
    @test dense_summary.coulomb_data_materialized === false
    @test dense_summary.ida_mwg_semantics === false
    @test dense_summary.exports_or_artifacts === false

    blocked_packet = CPGBOverlapBlocks.parent_overlap_axis_factor_packet(
        parent,
        (;
            x = (; pgdg_intermediate = (; overlap = overlap_1d.x)),
            y = (; pgdg_intermediate = (; overlap = overlap_1d.y)),
            z = (; pgdg_intermediate = (;)),
        ),
    )
    blocked_packet_blocks =
        CBPOverlapBlocks.cpb_overlap_axis_blocks(blocked_packet, interval_pair)
    blocked_packet_summary = CBPOverlapBlocks.summary(blocked_packet_blocks)

    @test blocked_packet_blocks.axis_blocks === nothing
    @test blocked_packet_summary.status === :blocked_cpb_overlap_axis_blocks
    @test blocked_packet_summary.blocker ===
          :missing_parent_axis_bundle_overlap_factors
    @test blocked_packet_summary.axis_blocks_available === false
    @test blocked_packet_summary.axis_block_shapes === :unavailable
    @test blocked_packet_summary.factor_space === :unavailable
    @test blocked_packet_summary.index_domain === :unavailable
    @test blocked_packet_summary.blocks_are_views === false
    @test blocked_packet_summary.blocks_are_copies === false
    @test blocked_packet_summary.dense_local_block_materialized === false

    blocked_dense =
        CBPOverlapBlocks.cpb_overlap_dense_block(blocked_packet_blocks)
    blocked_dense_summary = CBPOverlapBlocks.summary(blocked_dense)

    @test blocked_dense.axis_block_set === blocked_packet_blocks
    @test blocked_dense.dense_block === nothing
    @test blocked_dense_summary.status === :blocked_cpb_overlap_dense_block
    @test blocked_dense_summary.blocker ===
          :missing_parent_axis_bundle_overlap_factors
    @test blocked_dense_summary.dense_block_available === false
    @test blocked_dense_summary.dense_block_shape === :unavailable
    @test blocked_dense_summary.dense_block_eltype === :unavailable
    @test blocked_dense_summary.source_axis_block_summary ===
          blocked_packet_summary
    @test blocked_dense_summary.factor_space === :unavailable
    @test blocked_dense_summary.index_domain === :unavailable
    @test blocked_dense_summary.dense_local_block_materialized === false
    @test blocked_dense_summary.route_driver_wiring === false
    @test blocked_dense_summary.global_matrix_materialized === false

    outside_pair = CBPOverlapBlocks.cpb_interval_pair(
        parent,
        CPBOverlapBlocks.cpb(0:1, 1:2, 1:2; role = :outside_left),
        CPBOverlapBlocks.cpb(1:2, 1:2, 1:2; role = :inside_right),
    )
    blocked_interval_blocks =
        CBPOverlapBlocks.cpb_overlap_axis_blocks(packet, outside_pair)
    blocked_interval_summary = CBPOverlapBlocks.summary(blocked_interval_blocks)

    @test blocked_interval_blocks.axis_blocks === nothing
    @test blocked_interval_summary.status === :blocked_cpb_overlap_axis_blocks
    @test blocked_interval_summary.blocker === :left_x_interval_outside_parent
    @test blocked_interval_summary.factor_space ===
          :parent_axis_bundle_pgdg_intermediate
    @test blocked_interval_summary.index_domain === :parent_axis_indices
    @test blocked_interval_summary.blocks_are_views === false
    @test blocked_interval_summary.dense_local_block_materialized === false

    nonsliceable_packet = _overlap_blocks_packet_with_metadata(
        packet,
        merge(CPGBOverlapBlocks.summary(packet), (; sliceable_by_cpb = false)),
    )
    nonsliceable_blocks =
        CBPOverlapBlocks.cpb_overlap_axis_blocks(nonsliceable_packet, interval_pair)
    nonsliceable_summary = CBPOverlapBlocks.summary(nonsliceable_blocks)

    @test nonsliceable_blocks.axis_blocks === nothing
    @test nonsliceable_summary.status === :blocked_cpb_overlap_axis_blocks
    @test nonsliceable_summary.blocker === :overlap_packet_not_cpb_sliceable
    @test nonsliceable_summary.axis_blocks_available === false
    @test nonsliceable_summary.axis_block_shapes === :unavailable
    @test nonsliceable_summary.blocks_are_views === false
    @test nonsliceable_summary.dense_local_block_materialized === false

    other_parent = _overlap_blocks_parent()
    mismatch_pair = _overlap_blocks_interval_pair(other_parent)
    mismatch_blocks = CBPOverlapBlocks.cpb_overlap_axis_blocks(packet, mismatch_pair)
    mismatch_summary = CBPOverlapBlocks.summary(mismatch_blocks)

    @test packet.parent !== mismatch_pair.parent
    @test mismatch_blocks.axis_blocks === nothing
    @test mismatch_summary.status === :blocked_cpb_overlap_axis_blocks
    @test mismatch_summary.blocker === :parent_mismatch
    @test mismatch_summary.factor_space ===
          :parent_axis_bundle_pgdg_intermediate
    @test mismatch_summary.index_domain === :parent_axis_indices
    @test mismatch_summary.blocks_are_views === false
    @test mismatch_summary.dense_local_block_materialized === false
end

@testset "Cartesian CPB local overlap bridge comparison" begin
    parent = _overlap_blocks_parent(count = 2)
    overlap_1d = _overlap_blocks_driver_fixture_overlap_1d()
    packet = CPGBOverlapBlocks.parent_overlap_axis_factor_packet(
        parent,
        _overlap_blocks_axis_bundle(overlap_1d),
    )
    left = CPBOverlapBlocks.cpb(1:1, 1:2, 1:1; role = :driver_fixture_left)
    right = CPBOverlapBlocks.cpb(1:1, 1:2, 1:1; role = :driver_fixture_right)
    interval_pair = CBPOverlapBlocks.cpb_interval_pair(parent, left, right)
    axis_block_set = CBPOverlapBlocks.cpb_overlap_axis_blocks(packet, interval_pair)
    dense_block = CBPOverlapBlocks.cpb_overlap_dense_block(axis_block_set)
    dense_summary = CBPOverlapBlocks.summary(dense_block)
    packet_summary = CPGBOverlapBlocks.summary(packet)

    @test Matrix(dense_block.dense_block) ==
          _overlap_blocks_driver_fixture_expected_matrix()
    @test dense_summary.dense_block_shape == (2, 2)
    @test dense_summary.dense_local_block_materialized
    @test dense_summary.global_matrix_materialized === false
    @test dense_summary.route_driver_wiring === false
    @test dense_summary.factor_space === packet_summary.factor_space
    @test dense_summary.factor_convention === packet_summary.factor_convention
    @test dense_summary.normalization_convention ===
          packet_summary.normalization_convention
    @test dense_summary.index_domain === packet_summary.index_domain
    @test dense_summary.index_domain_source === packet_summary.index_domain_source
    @test dense_summary.index_domain_status === packet_summary.index_domain_status
    @test dense_summary.local_ordering ===
          :parent_compatible_x_slowest_z_fastest
    @test dense_summary.source_axis_block_summary.status ===
          :available_cpb_overlap_axis_blocks
end

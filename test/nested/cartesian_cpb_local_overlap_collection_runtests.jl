# Runtime role: tiny CPB local overlap block collection contract test.
#
# This validates local overlap provider block records and collections. It does
# not assign placement, wire route/global overlap, or exercise kinetic,
# position, x2, Coulomb, Hamiltonian, IDA/MWG, PQS Lowdin/projection, exports,
# or artifacts.

using Test
using GaussletBases

const CPBOverlapCollection = GaussletBases.CartesianCPB
const CPGBOverlapCollection = GaussletBases.CartesianParentGaussletBases
const CBPOverlapCollection = GaussletBases.CartesianCPBBlockProviders

function _collection_parent(; count = 3)
    axis = build_basis(MappedUniformBasisSpec(
        :G10;
        count,
        mapping = IdentityMapping(),
        reference_spacing = 1.0,
    ))
    return CPGBOverlapCollection.CartesianParentGaussletBasis3D(axis)
end

function _collection_overlap_1d()
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

function _collection_axis_bundle(overlap_1d = _collection_overlap_1d())
    return (;
        x = (; pgdg_intermediate = (; overlap = overlap_1d.x)),
        y = (; pgdg_intermediate = (; overlap = overlap_1d.y)),
        z = (; pgdg_intermediate = (; overlap = overlap_1d.z)),
    )
end

function _collection_packet(parent)
    return CPGBOverlapCollection.parent_overlap_axis_factor_packet(
        parent,
        _collection_axis_bundle(),
    )
end

function _collection_dense_block(parent, packet, left, right)
    interval_pair = CBPOverlapCollection.cpb_interval_pair(parent, left, right)
    axis_blocks =
        CBPOverlapCollection.cpb_overlap_axis_blocks(packet, interval_pair)
    return CBPOverlapCollection.cpb_overlap_dense_block(axis_blocks)
end

@testset "Cartesian CPB local overlap block collection" begin
    parent = _collection_parent()
    packet = _collection_packet(parent)
    left_a = CPBOverlapCollection.cpb(1:2, 1:2, 1:1; role = :left_a)
    right_a = CPBOverlapCollection.cpb(2:3, 1:2, 1:1; role = :right_a)
    dense_a = _collection_dense_block(parent, packet, left_a, right_a)
    record_a = CBPOverlapCollection.cpb_local_overlap_block_record(
        dense_a;
        block_key = (:product_a, :product_b),
    )
    record_summary = CBPOverlapCollection.summary(record_a)

    @test record_a.block_key === (:product_a, :product_b)
    @test record_a.source_block === dense_a
    @test record_summary.object_kind ===
          :cartesian_cpb_local_overlap_block_record_summary
    @test record_summary.term === :overlap
    @test record_summary.status === :available_cpb_local_overlap_block_record
    @test record_summary.blocker === nothing
    @test record_summary.source_kind === :cpb_overlap_dense_block
    @test record_summary.interval_pair_summary.status ===
          :available_cpb_interval_pair
    @test record_summary.axis_block_summary.status ===
          :available_cpb_overlap_axis_blocks
    @test record_summary.dense_block_summary.status ===
          :materialized_cpb_overlap_dense_block
    @test record_summary.dense_block_available
    @test record_summary.dense_block_shape == size(dense_a.dense_block)
    @test record_summary.local_ordering ===
          :parent_compatible_x_slowest_z_fastest
    @test record_summary.factor_space ===
          :parent_axis_bundle_pgdg_intermediate
    @test record_summary.factor_convention === :axis_bundle_one_body_overlap
    @test record_summary.index_domain === :parent_axis_indices
    @test record_summary.placement_status === :unassigned
    @test record_summary.global_matrix_materialized === false
    @test record_summary.route_driver_wiring === false
    @test !hasproperty(record_summary, :dense_block)
    @test !hasproperty(record_summary.dense_block_summary, :dense_block)

    collection_a =
        CBPOverlapCollection.cpb_local_overlap_block_collection((record_a,))
    collection_a_summary = CBPOverlapCollection.summary(collection_a)

    @test collection_a.records === (record_a,)
    @test collection_a_summary.status ===
          :available_cpb_local_overlap_block_collection
    @test collection_a_summary.blocker === nothing
    @test collection_a_summary.record_count == 1
    @test collection_a_summary.terms === (:overlap,)
    @test collection_a_summary.block_keys === ((:product_a, :product_b),)
    @test collection_a_summary.dense_block_count == 1
    @test collection_a_summary.placement_status === :unassigned
    @test collection_a_summary.global_matrix_materialized === false
    @test collection_a_summary.route_driver_wiring === false
    @test !hasproperty(collection_a_summary, :dense_block)

    left_b = CPBOverlapCollection.cpb(1:1, 2:3, 1:2; role = :left_b)
    right_b = CPBOverlapCollection.cpb(1:2, 1:1, 2:3; role = :right_b)
    dense_b = _collection_dense_block(parent, packet, left_b, right_b)
    collection_b =
        CBPOverlapCollection.cpb_local_overlap_block_collection((dense_a, dense_b))
    collection_b_summary = CBPOverlapCollection.summary(collection_b)

    @test length(collection_b.records) == 2
    @test collection_b.records[1].source_block === dense_a
    @test collection_b.records[2].source_block === dense_b
    @test collection_b_summary.status ===
          :available_cpb_local_overlap_block_collection
    @test collection_b_summary.record_count == 2
    @test collection_b_summary.block_keys ===
          (:local_overlap_block_1, :local_overlap_block_2)
    @test collection_b_summary.dense_block_count == 2
    @test collection_b_summary.global_matrix_materialized === false
    @test collection_b_summary.route_driver_wiring === false

    outside_pair = CBPOverlapCollection.cpb_interval_pair(
        parent,
        CPBOverlapCollection.cpb(0:1, 1:2, 1:1; role = :outside_left),
        right_a,
    )
    blocked_axis_blocks =
        CBPOverlapCollection.cpb_overlap_axis_blocks(packet, outside_pair)
    blocked_record = CBPOverlapCollection.cpb_local_overlap_block_record(
        blocked_axis_blocks;
        block_key = :blocked_outside_left,
    )
    blocked_record_summary = CBPOverlapCollection.summary(blocked_record)

    @test blocked_record.source_block === blocked_axis_blocks
    @test blocked_record_summary.status ===
          :blocked_cpb_local_overlap_block_record
    @test blocked_record_summary.blocker === :left_x_interval_outside_parent
    @test blocked_record_summary.source_kind === :cpb_overlap_axis_block_set
    @test blocked_record_summary.dense_block_summary === nothing
    @test blocked_record_summary.dense_block_available === false
    @test blocked_record_summary.dense_block_shape === :not_materialized
    @test blocked_record_summary.placement_status === :unassigned
    @test blocked_record_summary.global_matrix_materialized === false
    @test blocked_record_summary.route_driver_wiring === false
    @test !hasproperty(blocked_record_summary, :dense_block)

    blocked_collection =
        CBPOverlapCollection.cpb_local_overlap_block_collection((
            record_a,
            blocked_record,
        ))
    blocked_collection_summary =
        CBPOverlapCollection.summary(blocked_collection)

    @test blocked_collection_summary.status ===
          :blocked_cpb_local_overlap_block_collection
    @test blocked_collection_summary.blocker ===
          :blocked_cpb_local_overlap_block_records
    @test blocked_collection_summary.record_count == 2
    @test blocked_collection_summary.dense_block_count == 1
    @test blocked_collection_summary.blocked_record_count == 1
    @test blocked_collection_summary.blocked_record_blockers ===
          (:left_x_interval_outside_parent,)
    @test blocked_collection_summary.placement_status === :unassigned
    @test blocked_collection_summary.global_matrix_materialized === false
    @test blocked_collection_summary.route_driver_wiring === false
end

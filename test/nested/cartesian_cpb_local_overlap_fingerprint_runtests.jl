# Runtime role: tiny private-overlap audit fingerprint for CPB local overlap.
#
# This test verifies that the private overlap side can recognize the new
# CPB-provider local overlap source as local structured data. It does not make
# the route-global overlap input-facts path consume that local source.

using Test
using GaussletBases

const CPBLocalOverlapFingerprint = GaussletBases.CartesianCPB
const CPGBLocalOverlapFingerprint = GaussletBases.CartesianParentGaussletBases
const CBPLocalOverlapFingerprint = GaussletBases.CartesianCPBBlockProviders

function _local_overlap_fingerprint_parent()
    axis = build_basis(MappedUniformBasisSpec(
        :G10;
        count = 2,
        mapping = IdentityMapping(),
        reference_spacing = 1.0,
    ))
    return CPGBLocalOverlapFingerprint.CartesianParentGaussletBasis3D(axis)
end

function _local_overlap_fingerprint_overlap_1d()
    return (;
        x = [1.0 0.2; 0.2 1.1],
        y = [1.2 0.3; 0.3 1.3],
        z = [1.4 0.4; 0.4 1.5],
    )
end

function _local_overlap_fingerprint_axis_bundle()
    overlap_1d = _local_overlap_fingerprint_overlap_1d()
    return (;
        x = (; pgdg_intermediate = (; overlap = overlap_1d.x)),
        y = (; pgdg_intermediate = (; overlap = overlap_1d.y)),
        z = (; pgdg_intermediate = (; overlap = overlap_1d.z)),
    )
end

function _local_overlap_fingerprint_source()
    parent = _local_overlap_fingerprint_parent()
    packet = CPGBLocalOverlapFingerprint.parent_overlap_axis_factor_packet(
        parent,
        _local_overlap_fingerprint_axis_bundle(),
    )
    left = CPBLocalOverlapFingerprint.cpb(
        1:1,
        1:2,
        1:1;
        role = :private_overlap_local_left,
    )
    right = CPBLocalOverlapFingerprint.cpb(
        1:1,
        1:2,
        1:1;
        role = :private_overlap_local_right,
    )
    interval_pair =
        CBPLocalOverlapFingerprint.cpb_interval_pair(parent, left, right)
    axis_blocks =
        CBPLocalOverlapFingerprint.cpb_overlap_axis_blocks(packet, interval_pair)
    dense_block =
        CBPLocalOverlapFingerprint.cpb_overlap_dense_block(axis_blocks)
    return (; parent, packet, interval_pair, axis_blocks, dense_block)
end

@testset "Private overlap CPB local source fingerprint" begin
    source = _local_overlap_fingerprint_source()
    axis_block_summary =
        CBPLocalOverlapFingerprint.summary(source.axis_blocks)
    dense_summary =
        CBPLocalOverlapFingerprint.summary(source.dense_block)
    fingerprint =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_local_source_fingerprint(
            source.dense_block,
        )

    @test fingerprint.status ===
          :available_private_global_overlap_local_source_fingerprint
    @test fingerprint.blocker === nothing
    @test fingerprint.source_kind === :cpb_provider_local_overlap
    @test fingerprint.axis_blocks_available
    @test fingerprint.dense_local_block_materialized
    @test fingerprint.global_matrix_materialized === false
    @test fingerprint.global_overlap_matrix_materialized === false
    @test fingerprint.route_driver_wiring === false
    @test fingerprint.private_global_overlap_input_facts_available === false
    @test fingerprint.route_global_overlap_stage_source === false
    @test fingerprint.factor_space === :parent_axis_bundle_pgdg_intermediate
    @test fingerprint.factor_convention === :axis_bundle_one_body_overlap
    @test fingerprint.normalization_convention ===
          :not_separate_from_axis_bundle_one_body_overlap
    @test fingerprint.index_domain === :parent_axis_indices
    @test fingerprint.index_domain_source === :axis_bundle_contract
    @test fingerprint.index_domain_status ===
          :assumed_parent_axis_indexed_by_current_axis_bundle_contract
    @test fingerprint.local_ordering ===
          :parent_compatible_x_slowest_z_fastest
    @test fingerprint.left_shape == (x = 1, y = 2, z = 1)
    @test fingerprint.right_shape == (x = 1, y = 2, z = 1)
    @test fingerprint.dense_block_shape == (2, 2)
    @test fingerprint.dense_block_shape == dense_summary.dense_block_shape

    axis_fingerprint =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_local_source_fingerprint(
            source.axis_blocks,
        )
    @test axis_fingerprint.status ===
          :available_private_global_overlap_local_source_fingerprint
    @test axis_fingerprint.axis_blocks_available
    @test axis_fingerprint.dense_local_block_materialized === false
    @test axis_fingerprint.dense_block_shape === :not_materialized
    @test axis_fingerprint.left_shape ==
          axis_block_summary.interval_pair_summary.left_shape

    facts_report = (;
        retained_dimension = 2,
        parent_axis_counts = (2, 2, 2),
        cpb_provider_local_overlap = source.dense_block,
    )
    facts =
        GaussletBases._pqs_source_box_route_driver_private_global_overlap_input_facts(
            facts_report,
        )
    @test facts.status === :blocked_private_global_overlap_input_facts
    @test facts.blocker === :missing_parent_axis_bundle_overlap_factors
    @test facts.input_source === :unavailable
    @test facts.overlap_1d === nothing
    @test facts.global_dimension == 2
    @test facts.factor_space === :unavailable
    @test facts.factor_convention === :unavailable
end

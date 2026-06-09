# Runtime role: tiny parent-owned overlap axis factor packet contract test.
#
# This validates the first overlap-only parent factor packet. It does not
# exercise CPB slicing, route-driver wiring, Hamiltonian assembly, Coulomb,
# IDA/MWG, PQS Lowdin/projection, exports, or artifacts.

using Test
using GaussletBases

const CPAFPacketTest = GaussletBases.CartesianParentAxisFactors
const CPGBPacketTest = GaussletBases.CartesianParentGaussletBases

function _packet_test_parent()
    axis = build_basis(MappedUniformBasisSpec(
        :G10;
        count = 2,
        mapping = IdentityMapping(),
        reference_spacing = 1.0,
    ))
    return CPGBPacketTest.CartesianParentGaussletBasis3D(axis)
end

function _packet_test_overlap_1d()
    return (;
        x = [1.0 0.2; 0.2 1.1],
        y = [1.2 0.3; 0.3 1.3],
        z = [1.4 0.4; 0.4 1.5],
    )
end

function _packet_test_axis_bundle(overlap_1d = _packet_test_overlap_1d())
    return (;
        x = (; pgdg_intermediate = (; overlap = overlap_1d.x)),
        y = (; pgdg_intermediate = (; overlap = overlap_1d.y)),
        z = (; pgdg_intermediate = (; overlap = overlap_1d.z)),
    )
end

@testset "Cartesian parent overlap axis factor packet" begin
    parent = _packet_test_parent()
    overlap_1d = _packet_test_overlap_1d()
    packet = CPAFPacketTest.parent_overlap_axis_factor_packet(
        parent,
        _packet_test_axis_bundle(overlap_1d),
    )
    packet_summary = CPAFPacketTest.summary(packet)

    @test packet.parent === parent
    @test packet.overlap_1d.x === overlap_1d.x
    @test packet.overlap_1d.y === overlap_1d.y
    @test packet.overlap_1d.z === overlap_1d.z
    @test packet_summary.status === :available_parent_overlap_axis_factors
    @test packet_summary.blocker === nothing
    @test packet_summary.parent_axis_counts == (2, 2, 2)
    @test packet_summary.parent_axis_counts_source ===
          :parent_object_parent_axis_counts
    @test packet_summary.factor_space ===
          :parent_axis_bundle_pgdg_intermediate
    @test packet_summary.factor_convention ===
          :axis_bundle_one_body_overlap
    @test packet_summary.normalization_convention ===
          :not_separate_from_axis_bundle_one_body_overlap
    @test packet_summary.index_domain === :parent_axis_indices
    @test packet_summary.axis_order === (:x, :y, :z)
    @test packet_summary.bra_ket_order === (:bra, :ket)
    @test packet_summary.sliceable_by_cpb
    @test packet_summary.category_availability.overlap ===
          :available_parent_overlap_axis_factors
    @test packet_summary.category_availability.kinetic ===
          :not_requested_parent_kinetic_axis_factors
    @test packet_summary.category_availability.position ===
          :not_requested_parent_position_axis_factors
    @test packet_summary.category_availability.x2 ===
          :not_requested_parent_x2_axis_factors
    @test packet_summary.category_availability.coulomb ===
          :not_requested_parent_coulomb_axis_factors
    @test packet_summary.full_3d_parent_matrices === false
    @test packet_summary.cpb_slicing_implemented === false
    @test packet_summary.route_driver_wiring === false
    @test packet_summary.hamiltonian_data_materialized === false
    @test packet_summary.coulomb_data_materialized === false
    @test packet_summary.ida_mwg_semantics === false
    @test packet_summary.exports_or_artifacts === false

    incomplete_bundle = (;
        x = (; pgdg_intermediate = (; overlap = overlap_1d.x)),
        y = (; pgdg_intermediate = (; overlap = overlap_1d.y)),
        z = (; pgdg_intermediate = (;)),
    )
    blocked_packet = CPAFPacketTest.parent_overlap_axis_factor_packet(
        parent,
        incomplete_bundle,
    )
    blocked_summary = CPAFPacketTest.summary(blocked_packet)

    @test blocked_packet.parent === parent
    @test blocked_packet.overlap_1d === nothing
    @test blocked_summary.status === :blocked_parent_overlap_axis_factors
    @test blocked_summary.blocker ===
          :missing_parent_axis_bundle_overlap_factors
    @test blocked_summary.parent_axis_counts == (2, 2, 2)
    @test blocked_summary.factor_space === :unavailable
    @test blocked_summary.factor_convention === :unavailable
    @test blocked_summary.normalization_convention === :unavailable
    @test blocked_summary.index_domain === :unavailable
    @test blocked_summary.bra_ket_order === :unavailable
    @test blocked_summary.sliceable_by_cpb === false
    @test blocked_summary.category_availability.overlap ===
          :missing_parent_axis_bundle_overlap_factors
    @test blocked_summary.category_availability.kinetic ===
          :not_requested_parent_kinetic_axis_factors
    @test blocked_summary.category_availability.position ===
          :not_requested_parent_position_axis_factors
    @test blocked_summary.category_availability.x2 ===
          :not_requested_parent_x2_axis_factors
    @test blocked_summary.category_availability.coulomb ===
          :not_requested_parent_coulomb_axis_factors
end

# Runtime role: tiny parent-owned overlap axis factor packet contract test.
#
# This validates parent overlap and optional kinetic axis factor packets. It
# does not exercise CPB slicing, route-driver wiring, Hamiltonian assembly,
# Coulomb, IDA/MWG, PQS Lowdin/projection, exports, or artifacts.

using Test
using GaussletBases

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

function _packet_test_kinetic_1d()
    return (;
        x = [2.0 -0.2; -0.2 2.1],
        y = [2.2 -0.3; -0.3 2.3],
        z = [2.4 -0.4; -0.4 2.5],
    )
end

function _packet_test_axis_bundle(
    overlap_1d = _packet_test_overlap_1d();
    kinetic_1d = nothing,
    kinetic_source = :missing,
)
    axis_bundle(axis, overlap_matrix) =
        if kinetic_source === :pgdg_intermediate
            (; pgdg_intermediate = (;
                overlap = overlap_matrix,
                kinetic = getproperty(kinetic_1d, axis),
            ))
        elseif kinetic_source === :axis_property
            (;
                pgdg_intermediate = (; overlap = overlap_matrix),
                kinetic = getproperty(kinetic_1d, axis),
            )
        else
            (; pgdg_intermediate = (; overlap = overlap_matrix))
        end
    return (;
        x = axis_bundle(:x, overlap_1d.x),
        y = axis_bundle(:y, overlap_1d.y),
        z = axis_bundle(:z, overlap_1d.z),
    )
end

@testset "Cartesian parent overlap axis factor packet" begin
    parent = _packet_test_parent()
    overlap_1d = _packet_test_overlap_1d()
    packet = CPGBPacketTest.parent_overlap_axis_factor_packet(
        parent,
        _packet_test_axis_bundle(overlap_1d),
    )
    packet_summary = CPGBPacketTest.summary(packet)

    @test packet.parent === parent
    @test packet.overlap_1d.x === overlap_1d.x
    @test packet.overlap_1d.y === overlap_1d.y
    @test packet.overlap_1d.z === overlap_1d.z
    @test isnothing(packet.kinetic_1d)
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
    @test packet_summary.index_domain_source === :axis_bundle_contract
    @test packet_summary.index_domain_status ===
          :assumed_parent_axis_indexed_by_current_axis_bundle_contract
    @test packet_summary.axis_order === (:x, :y, :z)
    @test packet_summary.bra_ket_order === (:bra, :ket)
    @test packet_summary.sliceable_by_cpb
    @test packet_summary.sliceability_source === :index_domain_contract
    @test packet_summary.sliceability_status ===
          :sliceable_by_cpb_parent_axis_index_contract
    @test packet_summary.category_availability.overlap ===
          :available_parent_overlap_axis_factors
    @test packet_summary.category_availability.kinetic ===
          :missing_parent_axis_bundle_kinetic_factors
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
    blocked_packet = CPGBPacketTest.parent_overlap_axis_factor_packet(
        parent,
        incomplete_bundle,
    )
    blocked_summary = CPGBPacketTest.summary(blocked_packet)

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
    @test blocked_summary.index_domain_source === :unavailable
    @test blocked_summary.index_domain_status === :unavailable
    @test blocked_summary.bra_ket_order === :unavailable
    @test blocked_summary.sliceable_by_cpb === false
    @test blocked_summary.sliceability_source === :unavailable
    @test blocked_summary.sliceability_status === :unavailable
    @test blocked_summary.category_availability.overlap ===
          :missing_parent_axis_bundle_overlap_factors
    @test blocked_summary.category_availability.kinetic ===
          :missing_parent_axis_bundle_kinetic_factors
    @test blocked_summary.category_availability.position ===
          :not_requested_parent_position_axis_factors
    @test blocked_summary.category_availability.x2 ===
          :not_requested_parent_x2_axis_factors
    @test blocked_summary.category_availability.coulomb ===
          :not_requested_parent_coulomb_axis_factors

    nonmatrix_bundle = _packet_test_axis_bundle((
        x = :not_a_matrix,
        y = overlap_1d.y,
        z = overlap_1d.z,
    ))
    nonmatrix_packet = CPGBPacketTest.parent_overlap_axis_factor_packet(
        parent,
        nonmatrix_bundle,
    )
    nonmatrix_summary = CPGBPacketTest.summary(nonmatrix_packet)

    @test nonmatrix_packet.parent === parent
    @test nonmatrix_packet.overlap_1d === nothing
    @test nonmatrix_summary.status === :blocked_parent_overlap_axis_factors
    @test nonmatrix_summary.blocker === :x_overlap_axis_factor_not_matrix
    @test nonmatrix_summary.index_domain === :unavailable
    @test nonmatrix_summary.index_domain_source === :unavailable
    @test nonmatrix_summary.index_domain_status === :unavailable
    @test nonmatrix_summary.sliceable_by_cpb === false
    @test nonmatrix_summary.sliceability_source === :unavailable
    @test nonmatrix_summary.sliceability_status === :unavailable

    size_mismatch_bundle = _packet_test_axis_bundle((
        x = [1.0 0.2 0.0; 0.2 1.1 0.0],
        y = overlap_1d.y,
        z = overlap_1d.z,
    ))
    size_mismatch_packet = CPGBPacketTest.parent_overlap_axis_factor_packet(
        parent,
        size_mismatch_bundle,
    )
    size_mismatch_summary = CPGBPacketTest.summary(size_mismatch_packet)

    @test size_mismatch_packet.parent === parent
    @test size_mismatch_packet.overlap_1d === nothing
    @test size_mismatch_summary.status === :blocked_parent_overlap_axis_factors
    @test size_mismatch_summary.blocker ===
          :x_overlap_axis_factor_size_mismatch
    @test size_mismatch_summary.index_domain === :unavailable
    @test size_mismatch_summary.index_domain_source === :unavailable
    @test size_mismatch_summary.index_domain_status === :unavailable
    @test size_mismatch_summary.sliceable_by_cpb === false
    @test size_mismatch_summary.sliceability_source === :unavailable
    @test size_mismatch_summary.sliceability_status === :unavailable

    kinetic_1d = _packet_test_kinetic_1d()
    kinetic_packet = CPGBPacketTest.parent_overlap_axis_factor_packet(
        parent,
        _packet_test_axis_bundle(
            overlap_1d;
            kinetic_1d,
            kinetic_source = :pgdg_intermediate,
        ),
    )
    kinetic_summary = CPGBPacketTest.summary(kinetic_packet)

    @test kinetic_packet.parent === parent
    @test kinetic_packet.overlap_1d.x === overlap_1d.x
    @test kinetic_packet.kinetic_1d.x === kinetic_1d.x
    @test kinetic_packet.kinetic_1d.y === kinetic_1d.y
    @test kinetic_packet.kinetic_1d.z === kinetic_1d.z
    @test kinetic_summary.status === :available_parent_overlap_axis_factors
    @test kinetic_summary.blocker === nothing
    @test kinetic_summary.packet_kind === :overlap_kinetic_parent_axis_factors
    @test kinetic_summary.overlap_1d_available === true
    @test kinetic_summary.kinetic_status ===
          :available_parent_kinetic_axis_factors
    @test kinetic_summary.kinetic_1d_available === true
    @test kinetic_summary.kinetic_factor_space ===
          :parent_axis_bundle_pgdg_intermediate
    @test kinetic_summary.kinetic_factor_convention ===
          :axis_bundle_one_body_kinetic
    @test kinetic_summary.kinetic_index_domain === :parent_axis_indices
    @test kinetic_summary.kinetic_index_domain_source === :axis_bundle_contract
    @test kinetic_summary.kinetic_index_domain_status ===
          :assumed_parent_axis_indexed_by_current_axis_bundle_contract
    @test kinetic_summary.kinetic_sliceable_by_cpb === true
    @test kinetic_summary.kinetic_sliceability_source === :index_domain_contract
    @test kinetic_summary.kinetic_sliceability_status ===
          :sliceable_by_cpb_parent_axis_index_contract
    @test kinetic_summary.category_availability.kinetic ===
          :available_parent_kinetic_axis_factors
    @test kinetic_summary.full_3d_parent_matrices === false
    @test kinetic_summary.route_driver_wiring === false
    @test kinetic_summary.hamiltonian_data_materialized === false
    @test kinetic_summary.coulomb_data_materialized === false

    top_level_kinetic_packet = CPGBPacketTest.parent_overlap_axis_factor_packet(
        parent,
        _packet_test_axis_bundle(
            overlap_1d;
            kinetic_1d,
            kinetic_source = :axis_property,
        ),
    )
    top_level_kinetic_summary = CPGBPacketTest.summary(top_level_kinetic_packet)

    @test top_level_kinetic_packet.kinetic_1d.x === kinetic_1d.x
    @test top_level_kinetic_summary.kinetic_status ===
          :available_parent_kinetic_axis_factors

    nonmatrix_kinetic_packet = CPGBPacketTest.parent_overlap_axis_factor_packet(
        parent,
        _packet_test_axis_bundle(
            overlap_1d;
            kinetic_1d = (x = :not_a_matrix, y = kinetic_1d.y, z = kinetic_1d.z),
            kinetic_source = :pgdg_intermediate,
        ),
    )
    nonmatrix_kinetic_summary =
        CPGBPacketTest.summary(nonmatrix_kinetic_packet)

    @test nonmatrix_kinetic_packet.overlap_1d.x === overlap_1d.x
    @test isnothing(nonmatrix_kinetic_packet.kinetic_1d)
    @test nonmatrix_kinetic_summary.status ===
          :available_parent_overlap_axis_factors
    @test nonmatrix_kinetic_summary.kinetic_status ===
          :x_kinetic_axis_factor_not_matrix
    @test nonmatrix_kinetic_summary.category_availability.kinetic ===
          :x_kinetic_axis_factor_not_matrix

    size_mismatch_kinetic_packet =
        CPGBPacketTest.parent_overlap_axis_factor_packet(
            parent,
            _packet_test_axis_bundle(
                overlap_1d;
                kinetic_1d = (
                    x = [2.0 -0.2 0.0; -0.2 2.1 0.0],
                    y = kinetic_1d.y,
                    z = kinetic_1d.z,
                ),
                kinetic_source = :pgdg_intermediate,
            ),
        )
    size_mismatch_kinetic_summary =
        CPGBPacketTest.summary(size_mismatch_kinetic_packet)

    @test isnothing(size_mismatch_kinetic_packet.kinetic_1d)
    @test size_mismatch_kinetic_summary.kinetic_status ===
          :x_kinetic_axis_factor_size_mismatch
    @test size_mismatch_kinetic_summary.kinetic_sliceable_by_cpb === false
end

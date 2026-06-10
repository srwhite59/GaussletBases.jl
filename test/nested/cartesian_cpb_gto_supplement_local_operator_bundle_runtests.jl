# Runtime role: CPB provider-level GTO supplement local operator bundle.
#
# This test validates packaging of existing mixed CPB/GTO, GTO/GTO, and
# by-center nuclear provider blocks into one compact supplement operator bundle.
# It does not add route/global placement, WL/PQS realization, Hamiltonian
# assembly, driver options, exports, artifacts, or IDA/MWG/PQS semantics.

using Test
using GaussletBases

const CPBGTOSupplementBundle = GaussletBases.CartesianCPB
const CPBGTOSupplementBundleParent = GaussletBases.CartesianParentGaussletBases
const CPBGTOSupplementBundleProvider = GaussletBases.CartesianCPBBlockProviders

function _gto_bundle_parent(; count = 3)
    axis = build_basis(MappedUniformBasisSpec(
        :G10;
        count,
        mapping = IdentityMapping(),
        reference_spacing = 1.0,
    ))
    return CPBGTOSupplementBundleParent.CartesianParentGaussletBasis3D(axis)
end

function _gto_bundle_fixture()
    parent = _gto_bundle_parent()
    cpb = CPBGTOSupplementBundle.cpb(
        1:2,
        2:3,
        2:2;
        role = :gto_supplement_bundle_cpb,
    )
    supplement = basis_representation(
        legacy_atomic_gaussian_supplement("He", "cc-pVTZ"; lmax = 0),
    )
    expansion = CoulombGaussianExpansion(
        [1.0, 0.35],
        [0.3, 0.9];
        del = 1.0,
        s = 0.16,
        c = 0.01,
        maxu = 135.0,
    )
    center_records = (
        (;
            center_key = :left_center,
            center_index = 1,
            nuclear_charge = 2.0,
            location = (0.0, 0.0, -0.3),
        ),
        (;
            center_key = :right_center,
            center_index = 2,
            nuclear_charge = 2.0,
            location = (0.0, 0.0, 0.4),
        ),
    )
    return (; parent, cpb, supplement, expansion, center_records)
end

function _check_gto_bundle_summary_compact(summary)
    @test !hasproperty(summary, :dense_block)
    @test !hasproperty(summary, :mixed_blocks)
    @test !hasproperty(summary, :gto_blocks)
    @test !hasproperty(summary, :mixed_nuclear_by_center_blocks)
    @test !hasproperty(summary, :gto_nuclear_by_center_blocks)
    @test !hasproperty(summary, :axis_tables)
    @test !hasproperty(summary, :primitive_tables)
    @test !hasproperty(summary, :oracle_matrix)
    @test !hasproperty(summary, :route_global_matrix)
    @test !hasproperty(summary, :hamiltonian_matrix)
    return nothing
end

@testset "CPB GTO supplement local operator bundle" begin
    fixture = _gto_bundle_fixture()
    orbital_count = length(fixture.supplement.orbitals)
    support_count = CPBGTOSupplementBundle.support_count(fixture.cpb)

    bundle =
        CPBGTOSupplementBundleProvider.cpb_gto_supplement_local_operator_bundle(
            fixture.parent,
            fixture.cpb,
            fixture.supplement;
            expansion = fixture.expansion,
            center_records = fixture.center_records,
        )
    bundle_summary = CPBGTOSupplementBundleProvider.summary(bundle)

    @test bundle_summary.status ==
          :materialized_cpb_gto_supplement_local_operator_bundle
    @test isnothing(bundle_summary.blocker)
    @test bundle_summary.supplement_kind == fixture.supplement.supplement_kind
    @test bundle_summary.cpb_support_count == support_count
    @test bundle_summary.gto_orbital_count == orbital_count
    @test bundle_summary.mixed_one_body_term_count == 8
    @test bundle_summary.gto_one_body_term_count == 8
    @test bundle_summary.by_center_nuclear_count == 2
    @test bundle_summary.mixed_by_center_nuclear_count == 2
    @test bundle_summary.gto_by_center_nuclear_count == 2
    @test isempty(bundle_summary.missing_term_list)
    for term in (
        :mixed_gto_overlap,
        :mixed_gto_kinetic,
        :mixed_gto_position_x,
        :mixed_gto_position_y,
        :mixed_gto_position_z,
        :mixed_gto_x2_x,
        :mixed_gto_x2_y,
        :mixed_gto_x2_z,
        :gto_overlap,
        :gto_kinetic,
        :gto_position_x,
        :gto_position_y,
        :gto_position_z,
        :gto_x2_x,
        :gto_x2_y,
        :gto_x2_z,
        :mixed_gto_electron_nuclear_by_center,
        :gto_electron_nuclear_by_center,
    )
        @test term in bundle_summary.included_term_list
    end
    @test bundle_summary.nuclear_charge_recorded == true
    @test bundle_summary.nuclear_charge_applied == false
    @test bundle_summary.center_summation == false
    @test bundle_summary.centers_summed == false
    @test bundle_summary.center_indices == (1, 2)
    @test bundle_summary.nuclear_charges == (2.0, 2.0)
    @test bundle_summary.provider_level_local_blocks_materialized == true
    @test bundle_summary.provider_level_pilot == true
    @test bundle_summary.route_driver_wiring == false
    @test bundle_summary.route_global_matrix_materialized == false
    @test bundle_summary.global_matrix_materialized == false
    @test bundle_summary.hamiltonian_assembly == false
    @test bundle_summary.hamiltonian_data_materialized == false
    @test bundle_summary.ida_mwg_semantics == false
    @test bundle_summary.exports_or_artifacts == false
    _check_gto_bundle_summary_compact(bundle_summary)

    direct_mixed_overlap =
        CPBGTOSupplementBundleProvider.cpb_mixed_gto_overlap_block(
            fixture.parent,
            fixture.cpb,
            fixture.supplement,
        )
    direct_mixed_kinetic =
        CPBGTOSupplementBundleProvider.cpb_mixed_gto_kinetic_operator_block(
            fixture.parent,
            fixture.cpb,
            fixture.supplement,
        )
    direct_gto_overlap =
        CPBGTOSupplementBundleProvider.cpb_gto_overlap_operator_block(
            fixture.supplement,
        )
    direct_gto_kinetic =
        CPBGTOSupplementBundleProvider.cpb_gto_kinetic_operator_block(
            fixture.supplement,
        )
    @test bundle.mixed_blocks.overlap.dense_block ≈
          direct_mixed_overlap.dense_block atol = 0.0 rtol = 0.0
    @test bundle.mixed_blocks.kinetic.dense_block ≈
          direct_mixed_kinetic.dense_block atol = 0.0 rtol = 0.0
    @test bundle.gto_blocks.overlap.dense_block ≈
          direct_gto_overlap.dense_block atol = 0.0 rtol = 0.0
    @test bundle.gto_blocks.kinetic.dense_block ≈
          direct_gto_kinetic.dense_block atol = 0.0 rtol = 0.0

    @test size(bundle.mixed_blocks.overlap.dense_block) ==
          (support_count, orbital_count)
    @test size(bundle.gto_blocks.overlap.dense_block) ==
          (orbital_count, orbital_count)
    @test length(bundle.mixed_nuclear_by_center_blocks) == 2
    @test length(bundle.gto_nuclear_by_center_blocks) == 2
    @test !isapprox(
        bundle.mixed_nuclear_by_center_blocks[1].dense_block,
        bundle.mixed_nuclear_by_center_blocks[2].dense_block;
        atol = 1.0e-12,
        rtol = 1.0e-12,
    )
    @test !isapprox(
        bundle.gto_nuclear_by_center_blocks[1].dense_block,
        bundle.gto_nuclear_by_center_blocks[2].dense_block;
        atol = 1.0e-12,
        rtol = 1.0e-12,
    )
    direct_mixed_nuclear =
        CPBGTOSupplementBundleProvider.cpb_mixed_gto_nuclear_by_center_block(
            fixture.parent,
            fixture.cpb,
            fixture.supplement,
            fixture.expansion,
            fixture.center_records[1],
        )
    direct_gto_nuclear =
        CPBGTOSupplementBundleProvider.cpb_gto_nuclear_by_center_block(
            fixture.supplement,
            fixture.expansion,
            fixture.center_records[1],
        )
    @test bundle.mixed_nuclear_by_center_blocks[1].dense_block ≈
          direct_mixed_nuclear.dense_block atol = 0.0 rtol = 0.0
    @test bundle.gto_nuclear_by_center_blocks[1].dense_block ≈
          direct_gto_nuclear.dense_block atol = 0.0 rtol = 0.0
    for block in (
        bundle.mixed_nuclear_by_center_blocks...,
        bundle.gto_nuclear_by_center_blocks...,
    )
        block_summary = CPBGTOSupplementBundleProvider.summary(block)
        @test block_summary.by_center == true
        @test block_summary.center_summation == false
        @test block_summary.centers_summed == false
        @test block_summary.nuclear_charge_recorded == true
        @test block_summary.nuclear_charge_applied == false
        @test block_summary.oracle_convention == :qw_nuclear_by_center_uncharged
    end

    blocked_bundle =
        CPBGTOSupplementBundleProvider.cpb_gto_supplement_local_operator_bundle(
            fixture.parent,
            fixture.cpb,
            fixture.supplement;
            center_records = (fixture.center_records[1],),
        )
    blocked_summary = CPBGTOSupplementBundleProvider.summary(blocked_bundle)
    @test blocked_summary.status ==
          :blocked_cpb_gto_supplement_local_operator_bundle
    @test blocked_summary.blocker == :missing_coulomb_gaussian_expansion
    @test :mixed_gto_electron_nuclear_by_center in blocked_summary.missing_term_list
    @test :gto_electron_nuclear_by_center in blocked_summary.missing_term_list
    @test :mixed_gto_overlap in blocked_summary.included_term_list
    @test blocked_summary.by_center_nuclear_count == 1
    @test blocked_summary.nuclear_charge_recorded == true
    @test blocked_summary.nuclear_charge_applied == false
    @test blocked_summary.center_summation == false
    @test blocked_summary.provider_level_local_blocks_materialized == false
    @test blocked_summary.route_driver_wiring == false
    @test blocked_summary.hamiltonian_data_materialized == false
    @test blocked_summary.exports_or_artifacts == false
    _check_gto_bundle_summary_compact(blocked_summary)
end

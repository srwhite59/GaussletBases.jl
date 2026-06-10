# Runtime role: CPB-local GTO nuclear-attraction by-center pilots.
#
# This test validates provider-level mixed CPB/GTO and GTO/GTO by-center
# nuclear-attraction blocks against the existing QW/GTO by-center oracle. It
# keeps centers separate and does not assemble a Hamiltonian, place route/global
# matrices, add IDA/MWG/PQS semantics, or export artifacts.

using Test
using LinearAlgebra
using GaussletBases

const CPBGTONuclear = GaussletBases.CartesianCPB
const CPBGTONuclearProvider = GaussletBases.CartesianCPBBlockProviders
const CPBGTONuclearParent = GaussletBases.CartesianParentGaussletBases

function _gto_nuclear_fixture()
    basis = bond_aligned_homonuclear_qw_basis(
        bond_length = 1.2;
        core_spacing = 0.5,
        xmax_parallel = 1.5,
        xmax_transverse = 1.0,
    )
    expansion = CoulombGaussianExpansion(
        [1.0, 0.35],
        [0.3, 0.9];
        del = 1.0,
        s = 0.16,
        c = 0.01,
        maxu = 135.0,
    )
    bundle_x = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        basis.basis_x;
        exponents = expansion.exponents,
        center = 0.0,
        backend = :numerical_reference,
    )
    bundle_y = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        basis.basis_y;
        exponents = expansion.exponents,
        center = 0.0,
        backend = :numerical_reference,
    )
    bundle_z = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        basis.basis_z;
        exponents = expansion.exponents,
        center = 0.0,
        backend = :numerical_reference,
    )
    bundles = GaussletBases._CartesianNestedAxisBundles3D(
        bundle_x,
        bundle_y,
        bundle_z,
    )
    legacy_supplement = legacy_bond_aligned_diatomic_gaussian_supplement(
        "H",
        "cc-pVTZ",
        basis.nuclei;
        lmax = 0,
        max_width = 1.0,
    )
    supplement = basis_representation(legacy_supplement)
    internal_supplement =
        GaussletBases._bond_aligned_diatomic_cartesian_shell_supplement_3d(
            legacy_supplement,
        )
    oracle = GaussletBases._qwrg_diatomic_cartesian_shell_blocks_3d(
        bundles,
        internal_supplement,
        basis,
        expansion,
        basis.nuclear_charges,
    )
    parent = CPBGTONuclearParent.CartesianParentGaussletBasis3D(basis)
    cpb = CPBGTONuclear.cpb(
        1:2,
        2:3,
        3:4;
        role = :mixed_gto_nuclear_by_center_cpb,
    )
    center_records = Tuple(
        (;
            center_key = Symbol("nucleus_", index),
            center_index = index,
            nuclear_charge = basis.nuclear_charges[index],
            location = basis.nuclei[index],
        ) for index in eachindex(basis.nuclei)
    )
    return (;
        basis,
        expansion,
        supplement,
        oracle,
        parent,
        cpb,
        center_records,
    )
end

function _gto_nuclear_rows(parent, cpb)
    intervals = CPBGTONuclear.intervals(cpb)
    rows = Int[]
    for ix in intervals[1], iy in intervals[2], iz in intervals[3]
        push!(rows, CPBGTONuclearParent.parent_flat_index(parent, ix, iy, iz))
    end
    return rows
end

function _check_gto_nuclear_common_summary(summary; expected_shape, expected_orbital_count)
    @test summary.status in (
        :materialized_cpb_mixed_gto_nuclear_by_center_block,
        :materialized_cpb_gto_nuclear_by_center_block,
    )
    @test isnothing(summary.blocker)
    @test summary.term == :electron_nuclear_by_center
    @test summary.by_center == true
    @test summary.centers_summed == false
    @test summary.nuclear_charge_applied == false
    @test summary.charge_application_stage == :hamiltonian_or_center_summation
    @test summary.nuclear_attraction_sign_applied == true
    @test summary.galerkin_operator == true
    @test summary.cpb_integral_weights_applied == false
    @test summary.provider_level_pilot == true
    @test summary.gaussian_expansion_loop == :inner_local_contraction
    @test summary.gaussian_term_count == 2
    @test summary.orbital_count == expected_orbital_count
    @test summary.dense_block_available == true
    @test summary.dense_block_shape == expected_shape
    @test summary.provider_level_local_matrix_materialized == true
    @test summary.route_driver_wiring == false
    @test summary.route_global_matrix_materialized == false
    @test summary.global_matrix_materialized == false
    @test summary.hamiltonian_assembly == false
    @test summary.hamiltonian_data_materialized == false
    @test summary.coulomb_data_materialized == false
    @test summary.ida_mwg_semantics == false
    @test summary.pqs_lowdin_materialized == false
    @test summary.exports_or_artifacts == false
    @test !hasproperty(summary, :dense_block)
    @test !hasproperty(summary, :axis_tables)
    @test !hasproperty(summary, :primitive_tables)
    @test !hasproperty(summary, :oracle_matrix)
    @test !hasproperty(summary, :nuclear_ga_by_center)
    @test !hasproperty(summary, :nuclear_aa_by_center)
    @test !hasproperty(summary, :route_global_matrix)
    @test !hasproperty(summary, :retained_matrix)
    @test !hasproperty(summary, :payload)
    return nothing
end

@testset "CPB provider GTO nuclear by-center blocks" begin
    fixture = _gto_nuclear_fixture()
    rows = _gto_nuclear_rows(fixture.parent, fixture.cpb)
    orbital_count = length(fixture.supplement.orbitals)

    mixed_blocks = [
        CPBGTONuclearProvider.cpb_mixed_gto_nuclear_by_center_block(
            fixture.parent,
            fixture.cpb,
            fixture.supplement,
            fixture.expansion,
            center_record,
        ) for center_record in fixture.center_records
    ]
    self_blocks = [
        CPBGTONuclearProvider.cpb_gto_nuclear_by_center_block(
            fixture.supplement,
            fixture.expansion,
            center_record,
        ) for center_record in fixture.center_records
    ]

    @test length(mixed_blocks) == 2
    @test length(self_blocks) == 2
    @test !isapprox(
        mixed_blocks[1].dense_block,
        mixed_blocks[2].dense_block;
        atol = 1.0e-12,
        rtol = 1.0e-12,
    )
    @test !isapprox(
        self_blocks[1].dense_block,
        self_blocks[2].dense_block;
        atol = 1.0e-12,
        rtol = 1.0e-12,
    )

    for center_index in eachindex(fixture.center_records)
        mixed = mixed_blocks[center_index]
        mixed_summary = CPBGTONuclearProvider.summary(mixed)
        _check_gto_nuclear_common_summary(
            mixed_summary;
            expected_shape = (CPBGTONuclear.support_count(fixture.cpb), orbital_count),
            expected_orbital_count = orbital_count,
        )
        @test mixed_summary.source_kind ==
              :mixed_gausslet_gto_supplement_nuclear_by_center
        @test mixed_summary.source_oracle_helper == :nuclear_ga_by_center
        @test mixed_summary.factor_source_path ==
              :_qwrg_atomic_axis_factor_cross_data
        @test mixed_summary.mixed_gto_pilot == true
        @test mixed_summary.gto_supplement_self_block == false
        @test mixed_summary.center_index == center_index
        @test mixed_summary.center_coordinates ==
              fixture.center_records[center_index].location
        @test mixed_summary.nuclear_charge ==
              fixture.center_records[center_index].nuclear_charge
        @test isapprox(
            mixed.dense_block,
            fixture.oracle.nuclear_ga_by_center[center_index][rows, :];
            atol = 1.0e-11,
            rtol = 1.0e-11,
        )

        self = self_blocks[center_index]
        self_summary = CPBGTONuclearProvider.summary(self)
        _check_gto_nuclear_common_summary(
            self_summary;
            expected_shape = (orbital_count, orbital_count),
            expected_orbital_count = orbital_count,
        )
        @test self_summary.source_kind == :gto_supplement_nuclear_by_center
        @test self_summary.source_oracle_helper == :nuclear_aa_by_center
        @test self_summary.factor_source_path ==
              :_qwrg_atomic_axis_factor_aa_data
        @test self_summary.mixed_gto_pilot == false
        @test self_summary.gto_supplement_self_block == true
        @test self_summary.center_index == center_index
        @test self_summary.center_coordinates ==
              fixture.center_records[center_index].location
        @test self_summary.nuclear_charge ==
              fixture.center_records[center_index].nuclear_charge
        @test isapprox(
            self.dense_block,
            fixture.oracle.nuclear_aa_by_center[center_index];
            atol = 1.0e-11,
            rtol = 1.0e-11,
        )
    end

    missing_expansion = CPBGTONuclearProvider.cpb_gto_nuclear_by_center_block(
        fixture.supplement,
        nothing,
        fixture.center_records[1],
    )
    missing_expansion_summary = CPBGTONuclearProvider.summary(missing_expansion)
    @test missing_expansion_summary.status ==
          :blocked_cpb_gto_nuclear_by_center_block
    @test missing_expansion_summary.blocker == :missing_coulomb_gaussian_expansion
    @test isnothing(missing_expansion.dense_block)
    @test missing_expansion_summary.dense_block_available == false
    @test missing_expansion_summary.route_driver_wiring == false
    @test missing_expansion_summary.hamiltonian_data_materialized == false

    missing_center = CPBGTONuclearProvider.cpb_mixed_gto_nuclear_by_center_block(
        fixture.parent,
        fixture.cpb,
        fixture.supplement,
        fixture.expansion,
        nothing,
    )
    missing_center_summary = CPBGTONuclearProvider.summary(missing_center)
    @test missing_center_summary.status ==
          :blocked_cpb_mixed_gto_nuclear_by_center_block
    @test missing_center_summary.blocker ==
          :missing_electron_nuclear_center_record
    @test isnothing(missing_center.dense_block)
    @test missing_center_summary.dense_block_available == false
    @test missing_center_summary.route_driver_wiring == false
    @test missing_center_summary.hamiltonian_data_materialized == false
end

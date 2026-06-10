# Runtime role: CPB provider-level GTO/GTO one-body supplement blocks.
#
# This test validates whole-supplement GTO/GTO overlap, kinetic, position, and
# x2 Galerkin blocks against existing supplement self oracles. It does not add
# mixed CPB-row blocks, nuclear/Coulomb GTO terms, route/global placement,
# WL/PQS realization, Hamiltonian assembly, exports, artifacts, or
# IDA/MWG/PQS semantics.

using Test
using GaussletBases

const CBPGTOSelf = GaussletBases.CartesianCPBBlockProviders

function _gto_self_internal_orbital(orbital)
    return GaussletBases._AtomicCartesianShellOrbital3D(
        orbital.label,
        orbital.angular_powers[1],
        orbital.angular_powers[2],
        orbital.angular_powers[3],
        orbital.exponents,
        orbital.coefficients,
        orbital.center,
    )
end

function _gto_self_qw_oracle(supplement)
    expansion = CoulombGaussianExpansion(
        [1.0],
        [1.0];
        del = 1.0,
        s = 1.0,
        c = 1.0,
        maxu = 1.0,
    )
    internal_supplement = (;
        orbitals = [
            _gto_self_internal_orbital(orbital) for orbital in supplement.orbitals
        ],
    )
    return GaussletBases._qwrg_cartesian_shell_self_moment_blocks_3d(
        internal_supplement,
        expansion;
        include_factor_terms = false,
    )
end

function _check_gto_self_summary_common(
    block_summary,
    supplement;
    expected_status,
    expected_term,
    expected_source_kind,
)
    orbital_count = length(supplement.orbitals)
    @test block_summary.status == expected_status
    @test isnothing(block_summary.blocker)
    @test block_summary.term == expected_term
    @test block_summary.source_kind == expected_source_kind
    @test block_summary.supplement_representation_kind ==
          :cartesian_gaussian_shell_supplement_representation
    @test block_summary.supplement_kind == supplement.supplement_kind
    @test block_summary.source_metadata_summary.source_kind ==
          supplement.metadata.source_kind
    @test block_summary.source_metadata_summary.lmax == supplement.metadata.lmax
    @test block_summary.source_metadata_summary.uncontracted ==
          supplement.metadata.uncontracted
    @test block_summary.orbital_count == orbital_count
    @test block_summary.orbital_labels ==
          Tuple(orbital.label for orbital in supplement.orbitals)
    @test block_summary.orbital_angular_powers ==
          Tuple(orbital.angular_powers for orbital in supplement.orbitals)
    @test block_summary.orbital_centers ==
          Tuple(orbital.center for orbital in supplement.orbitals)
    @test block_summary.primitive_normalization_convention ==
          :axiswise_normalized_cartesian_gaussian
    @test block_summary.contraction_convention ==
          :orbital_coefficients_contract_primitive_axis_tables
    @test block_summary.center_shift_convention ==
          :explicit_axis_center_coordinates
    @test block_summary.shell_power_order == (:x, :y, :z)
    @test block_summary.axis_order == (:x, :y, :z)
    @test block_summary.left_basis_kind ==
          :cartesian_gaussian_shell_supplement_representation
    @test block_summary.right_basis_kind ==
          :cartesian_gaussian_shell_supplement_representation
    @test block_summary.left_orbital_count == orbital_count
    @test block_summary.right_orbital_count == orbital_count
    @test block_summary.left_shape == (gto = orbital_count,)
    @test block_summary.right_shape == (gto = orbital_count,)
    @test block_summary.dense_block_available == true
    @test block_summary.dense_block_shape == (orbital_count, orbital_count)
    @test block_summary.representation ==
          :dense_provider_level_gto_supplement_galerkin_block
    @test block_summary.galerkin_operator == true
    @test block_summary.provider_level_local_matrix_materialized == true
    @test block_summary.provider_level_pilot == true
    @test block_summary.gto_supplement_self_block == true
    @test block_summary.mixed_gto_pilot == false
    @test block_summary.parent_one_body_factor_packet_consumed == false
    @test block_summary.gto_axis_integral_source ==
          :qw_polynomial_gaussian_primitive_tables
    @test block_summary.route_driver_wiring == false
    @test block_summary.route_global_matrix_materialized == false
    @test block_summary.global_matrix_materialized == false
    @test block_summary.hamiltonian_data_materialized == false
    @test block_summary.coulomb_data_materialized == false
    @test block_summary.ida_mwg_semantics == false
    @test block_summary.exports_or_artifacts == false
    @test !hasproperty(block_summary, :dense_block)
    @test !hasproperty(block_summary, :axis_tables)
    @test !hasproperty(block_summary, :primitive_tables)
    @test !hasproperty(block_summary, :oracle_matrix)
    @test !hasproperty(block_summary, :route_global_matrix)
    @test !hasproperty(block_summary, :retained_matrix)
    @test !hasproperty(block_summary, :payload)
    return nothing
end

@testset "CPB provider GTO/GTO supplement one-body blocks" begin
    supplement = basis_representation(
        legacy_atomic_gaussian_supplement("He", "cc-pVTZ"; lmax = 0),
    )
    @test length(supplement.orbitals) == 3

    overlap = CBPGTOSelf.cpb_gto_overlap_operator_block(supplement)
    overlap_summary = CBPGTOSelf.summary(overlap)
    _check_gto_self_summary_common(
        overlap_summary,
        supplement;
        expected_status = :materialized_cpb_gto_supplement_overlap_block,
        expected_term = :gto_overlap,
        expected_source_kind = :gto_supplement_self_overlap,
    )
    @test overlap_summary.oracle_source == :_cartesian_supplement_cross_overlap
    overlap_oracle = GaussletBases._cartesian_supplement_cross_overlap(
        supplement,
        supplement,
    )
    @test isapprox(overlap.dense_block, overlap_oracle; atol = 1.0e-12, rtol = 1.0e-12)

    moment_oracle = _gto_self_qw_oracle(supplement)
    for axis in (:x, :y, :z)
        position = CBPGTOSelf.cpb_gto_position_operator_block(supplement; axis)
        position_summary = CBPGTOSelf.summary(position)
        _check_gto_self_summary_common(
            position_summary,
            supplement;
            expected_status =
                :materialized_cpb_gto_supplement_coordinate_moment_block,
            expected_term = Symbol("gto_position_", axis),
            expected_source_kind = :gto_supplement_self_coordinate_moment,
        )
        @test position_summary.coordinate_moment == :position
        @test position_summary.active_axis == axis
        @test position_summary.xpower == 1
        @test position_summary.oracle_source ==
              :_qwrg_cartesian_shell_self_moment_blocks_3d
        @test isapprox(
            position.dense_block,
            getproperty(moment_oracle, Symbol("position_", axis, "_aa"));
            atol = 1.0e-12,
            rtol = 1.0e-12,
        )

        x2 = CBPGTOSelf.cpb_gto_x2_operator_block(supplement; axis)
        x2_summary = CBPGTOSelf.summary(x2)
        _check_gto_self_summary_common(
            x2_summary,
            supplement;
            expected_status =
                :materialized_cpb_gto_supplement_coordinate_moment_block,
            expected_term = Symbol("gto_x2_", axis),
            expected_source_kind = :gto_supplement_self_coordinate_moment,
        )
        @test x2_summary.coordinate_moment == :x2
        @test x2_summary.active_axis == axis
        @test x2_summary.xpower == 2
        @test x2_summary.oracle_source ==
              :_qwrg_cartesian_shell_self_moment_blocks_3d
        @test isapprox(
            x2.dense_block,
            getproperty(moment_oracle, Symbol("x2_", axis, "_aa"));
            atol = 1.0e-12,
            rtol = 1.0e-12,
        )
    end

    kinetic = CBPGTOSelf.cpb_gto_kinetic_operator_block(supplement)
    kinetic_summary = CBPGTOSelf.summary(kinetic)
    _check_gto_self_summary_common(
        kinetic_summary,
        supplement;
        expected_status = :materialized_cpb_gto_supplement_kinetic_block,
        expected_term = :gto_kinetic,
        expected_source_kind = :gto_supplement_self_kinetic,
    )
    @test kinetic_summary.oracle_source ==
          :_qwrg_cartesian_shell_self_moment_blocks_3d
    @test isapprox(
        kinetic.dense_block,
        moment_oracle.kinetic_aa;
        atol = 1.0e-12,
        rtol = 1.0e-12,
    )

    uncontracted_supplement = basis_representation(
        legacy_atomic_gaussian_supplement(
            "He",
            "cc-pVTZ";
            lmax = 0,
            uncontracted = true,
        ),
    )
    @test length(uncontracted_supplement.orbitals) == 6
    uncontracted_overlap =
        CBPGTOSelf.cpb_gto_overlap_operator_block(uncontracted_supplement)
    uncontracted_summary = CBPGTOSelf.summary(uncontracted_overlap)
    _check_gto_self_summary_common(
        uncontracted_summary,
        uncontracted_supplement;
        expected_status = :materialized_cpb_gto_supplement_overlap_block,
        expected_term = :gto_overlap,
        expected_source_kind = :gto_supplement_self_overlap,
    )
    @test uncontracted_summary.source_metadata_summary.uncontracted == true
    uncontracted_oracle = GaussletBases._cartesian_supplement_cross_overlap(
        uncontracted_supplement,
        uncontracted_supplement,
    )
    @test isapprox(
        uncontracted_overlap.dense_block,
        uncontracted_oracle;
        atol = 1.0e-12,
        rtol = 1.0e-12,
    )

    unsupported_axis =
        CBPGTOSelf.cpb_gto_position_operator_block(supplement; axis = :q)
    unsupported_summary = CBPGTOSelf.summary(unsupported_axis)
    @test unsupported_summary.status ==
          :blocked_cpb_gto_supplement_coordinate_moment_block
    @test unsupported_summary.blocker == :unsupported_gto_position_axis
    @test unsupported_summary.dense_block_available == false
    @test isnothing(unsupported_axis.dense_block)
    @test unsupported_summary.route_driver_wiring == false
    @test unsupported_summary.hamiltonian_data_materialized == false
    @test unsupported_summary.coulomb_data_materialized == false
end

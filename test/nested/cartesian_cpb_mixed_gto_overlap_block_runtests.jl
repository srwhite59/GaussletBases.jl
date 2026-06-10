# Runtime role: CPB-local mixed gausslet/GTO overlap pilot.
#
# This test validates tiny provider-level mixed gausslet/GTO overlap,
# coordinate-moment, and kinetic blocks against existing Cartesian/GTO and
# QW/GTO oracles. It does not add nuclear/Coulomb terms, WL/PQS realization,
# retained transforms, route/global placement, driver wiring, Hamiltonian
# assembly, IDA/MWG, PQS Lowdin/projection, exports, or artifacts.

using Test
using GaussletBases

const CPBMixedGTO = GaussletBases.CartesianCPB
const CPGBMixedGTO = GaussletBases.CartesianParentGaussletBases
const CBPMixedGTO = GaussletBases.CartesianCPBBlockProviders

function _mixed_gto_parent(; count = 3)
    axis = build_basis(MappedUniformBasisSpec(
        :G10;
        count,
        mapping = IdentityMapping(),
        reference_spacing = 1.0,
    ))
    return axis, CPGBMixedGTO.CartesianParentGaussletBasis3D(axis)
end

function _mixed_gto_orbital(;
    label = "test_px_dz",
    angular_powers = (1, 0, 2),
    center = (0.1, -0.2, 0.3),
    exponents = [0.7, 1.3],
    coefficients = [0.8, -0.4],
    primitive_normalization = :axiswise_normalized_cartesian_gaussian,
)
    return CartesianGaussianShellOrbitalRepresentation3D(
        label,
        angular_powers,
        center,
        Float64[Float64(value) for value in exponents],
        Float64[Float64(value) for value in coefficients],
        primitive_normalization,
    )
end

function _mixed_gto_oracle_rows(parent, cpb)
    intervals = CPBMixedGTO.intervals(cpb)
    rows = Int[]
    for ix in intervals[1], iy in intervals[2], iz in intervals[3]
        push!(rows, CPGBMixedGTO.parent_flat_index(parent, ix, iy, iz))
    end
    return rows
end

function _mixed_gto_internal_orbital(orbital)
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

function _mixed_gto_qw_cross_moment_oracle(axis_basis, orbital)
    expansion = CoulombGaussianExpansion(
        [1.0],
        [1.0];
        del = 1.0,
        s = 1.0,
        c = 1.0,
        maxu = 1.0,
    )
    bundle = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        axis_basis;
        exponents = expansion.exponents,
        center = 0.0,
        backend = :numerical_reference,
    )
    proxy = GaussletBases._qwrg_mapped_supplement_proxy_layer(axis_basis, bundle)
    supplement = (; orbitals = [ _mixed_gto_internal_orbital(orbital) ])
    return GaussletBases._qwrg_cartesian_shell_cross_moment_blocks_3d(
        (x = proxy, y = proxy, z = proxy),
        supplement,
        expansion,
        length(axis_basis)^3;
        include_factor_terms = false,
    )
end

function _check_blocked_mixed_gto(
    block,
    expected_blocker;
    expected_status = :blocked_cpb_mixed_gto_local_overlap_block,
)
    block_summary = CBPMixedGTO.summary(block)
    @test block_summary.status == expected_status
    @test block_summary.blocker == expected_blocker
    @test isnothing(block.dense_block)
    @test block_summary.dense_block_available == false
    @test block_summary.dense_block_shape == :unavailable
    @test block_summary.provider_level_local_matrix_materialized == false
    @test block_summary.route_driver_wiring == false
    @test block_summary.route_global_matrix_materialized == false
    @test block_summary.global_matrix_materialized == false
    @test block_summary.hamiltonian_data_materialized == false
    @test block_summary.coulomb_data_materialized == false
    @test block_summary.ida_mwg_semantics == false
    @test block_summary.exports_or_artifacts == false
    @test !hasproperty(block_summary, :dense_block)
    @test !hasproperty(block_summary, :axis_ops)
    @test !hasproperty(block_summary, :axis_tables)
    @test !hasproperty(block_summary, :primitive_tables)
    @test !hasproperty(block_summary, :oracle_matrix)
    @test !hasproperty(block_summary, :global_matrix)
    @test !hasproperty(block_summary, :global_overlap_matrix)
    @test !hasproperty(block_summary, :retained_blocks)
    return nothing
end

function _check_mixed_gto_summary_common(
    block_summary;
    expected_shape = (x = 2, y = 2, z = 1),
    expected_dense_shape = (4, 1),
)
    @test isnothing(block_summary.blocker)
    @test block_summary.supplement_representation_kind ==
          :cartesian_gaussian_shell_orbital
    @test block_summary.primitive_normalization ==
          :axiswise_normalized_cartesian_gaussian
    @test block_summary.formula_source ==
          :GaussianAnalyticIntegrals_polynomial_gaussian
    @test block_summary.local_shape == expected_shape
    @test block_summary.dense_block_shape == expected_dense_shape
    @test block_summary.local_ordering ==
          :parent_compatible_x_slowest_z_fastest
    @test block_summary.galerkin_operator == true
    @test block_summary.provider_level_local_matrix_materialized == true
    @test block_summary.route_driver_wiring == false
    @test block_summary.route_global_matrix_materialized == false
    @test block_summary.global_matrix_materialized == false
    @test block_summary.hamiltonian_data_materialized == false
    @test block_summary.coulomb_data_materialized == false
    @test block_summary.ida_mwg_semantics == false
    @test block_summary.exports_or_artifacts == false
    @test !hasproperty(block_summary, :dense_block)
    @test !hasproperty(block_summary, :axis_ops)
    @test !hasproperty(block_summary, :axis_tables)
    @test !hasproperty(block_summary, :primitive_tables)
    @test !hasproperty(block_summary, :oracle_matrix)
    @test !hasproperty(block_summary, :global_matrix)
    @test !hasproperty(block_summary, :global_overlap_matrix)
    @test !hasproperty(block_summary, :retained_blocks)
    return nothing
end

function _check_mixed_gto_one_body_wrapper_blockers(
    parent,
    cpb,
    orbital,
    expected_blocker,
)
    _check_blocked_mixed_gto(
        CBPMixedGTO.cpb_mixed_gto_position_operator_block(parent, cpb, orbital),
        expected_blocker;
        expected_status = :blocked_cpb_mixed_gto_coordinate_moment_local_block,
    )
    _check_blocked_mixed_gto(
        CBPMixedGTO.cpb_mixed_gto_x2_operator_block(parent, cpb, orbital),
        expected_blocker;
        expected_status = :blocked_cpb_mixed_gto_coordinate_moment_local_block,
    )
    _check_blocked_mixed_gto(
        CBPMixedGTO.cpb_mixed_gto_kinetic_operator_block(parent, cpb, orbital),
        expected_blocker;
        expected_status = :blocked_cpb_mixed_gto_kinetic_local_block,
    )
    return nothing
end

function _mixed_gto_moment_oracle_column(moment_oracle, moment::Symbol, axis::Symbol, rows)
    field =
        moment == :position ? Symbol("position_$(axis)_ga") :
        moment == :x2 ? Symbol("x2_$(axis)_ga") :
        error("unsupported test moment")
    return getproperty(moment_oracle, field)[rows, 1]
end

@testset "CPB mixed GTO overlap block" begin
    axis, parent = _mixed_gto_parent()
    orbital = _mixed_gto_orbital()
    cpb = CPBMixedGTO.cpb(1:2, 2:3, 3:3; role = :mixed_gto_overlap_cpb)

    block = CBPMixedGTO.cpb_mixed_gto_overlap_block(parent, cpb, orbital)
    block_summary = CBPMixedGTO.summary(block)

    supplement = CartesianGaussianShellSupplementRepresentation3D(
        :test_cartesian_shell,
        CartesianGaussianShellOrbitalRepresentation3D[orbital],
        (; source_kind = :test_mixed_gto_overlap_oracle),
    )
    parent_representation = GaussletBases._cartesian_direct_product_representation(axis)
    oracle = GaussletBases._cartesian_basis_supplement_cross(
        parent_representation,
        supplement,
    )
    oracle_rows = _mixed_gto_oracle_rows(parent, cpb)

    @test block_summary.status == :materialized_cpb_mixed_gto_local_overlap_block
    @test isnothing(block_summary.blocker)
    @test block_summary.term == :mixed_gto_overlap
    @test block_summary.source_kind == :mixed_gausslet_gto_supplement_overlap
    @test block_summary.supplement_representation_kind ==
          :cartesian_gaussian_shell_orbital
    @test block_summary.orbital_label == orbital.label
    @test block_summary.angular_powers == orbital.angular_powers
    @test block_summary.center == orbital.center
    @test block_summary.primitive_count == length(orbital.exponents)
    @test block_summary.primitive_normalization ==
          :axiswise_normalized_cartesian_gaussian
    @test block_summary.formula_source ==
          :GaussianAnalyticIntegrals_polynomial_gaussian
    @test block_summary.axis_kernel_source ==
          :existing_cartesian_basis_supplement_axis_cross
    _check_mixed_gto_summary_common(block_summary)
    @test block.dense_block[:, 1] ≈ oracle[oracle_rows, 1] atol = 1.0e-12 rtol = 1.0e-12

    moment_oracle = _mixed_gto_qw_cross_moment_oracle(axis, orbital)

    for active_axis in (:x, :y, :z)
        position_block = CBPMixedGTO.cpb_mixed_gto_position_operator_block(
            parent,
            cpb,
            orbital;
            axis = active_axis,
        )
        position_summary = CBPMixedGTO.summary(position_block)
        @test position_summary.status ==
              :materialized_cpb_mixed_gto_coordinate_moment_local_block
        @test position_summary.term == Symbol("mixed_gto_position_$(active_axis)")
        @test position_summary.coordinate_moment == :position
        @test position_summary.active_axis == active_axis
        @test position_summary.xpower == 1
        @test position_summary.source_kind ==
              :mixed_gausslet_gto_supplement_coordinate_moment
        @test position_summary.orbital_label == orbital.label
        @test position_summary.axis_kernel_source ==
              :existing_qw_polynomial_gaussian_axis_integrals
        _check_mixed_gto_summary_common(position_summary)
        @test position_block.dense_block[:, 1] ≈
              _mixed_gto_moment_oracle_column(
                  moment_oracle,
                  :position,
                  active_axis,
                  oracle_rows,
              ) atol = 1.0e-12 rtol = 1.0e-12

        x2_block = CBPMixedGTO.cpb_mixed_gto_x2_operator_block(
            parent,
            cpb,
            orbital;
            axis = active_axis,
        )
        x2_summary = CBPMixedGTO.summary(x2_block)
        @test x2_summary.status ==
              :materialized_cpb_mixed_gto_coordinate_moment_local_block
        @test x2_summary.term == Symbol("mixed_gto_x2_$(active_axis)")
        @test x2_summary.coordinate_moment == :x2
        @test x2_summary.active_axis == active_axis
        @test x2_summary.xpower == 2
        @test x2_summary.source_kind ==
              :mixed_gausslet_gto_supplement_coordinate_moment
        @test x2_summary.axis_kernel_source ==
              :existing_qw_polynomial_gaussian_axis_integrals
        _check_mixed_gto_summary_common(x2_summary)
        @test x2_block.dense_block[:, 1] ≈
              _mixed_gto_moment_oracle_column(
                  moment_oracle,
                  :x2,
                  active_axis,
                  oracle_rows,
              ) atol = 1.0e-12 rtol = 1.0e-12
    end

    kinetic_block = CBPMixedGTO.cpb_mixed_gto_kinetic_operator_block(
        parent,
        cpb,
        orbital,
    )
    kinetic_summary = CBPMixedGTO.summary(kinetic_block)
    @test kinetic_summary.status == :materialized_cpb_mixed_gto_kinetic_local_block
    @test kinetic_summary.term == :mixed_gto_kinetic
    @test kinetic_summary.source_kind == :mixed_gausslet_gto_supplement_kinetic
    @test kinetic_summary.kinetic_factor_form == :sum_of_axis_products
    @test kinetic_summary.kinetic_component_terms ==
          (:kinetic_x_component, :kinetic_y_component, :kinetic_z_component)
    @test kinetic_summary.axis_kernel_source ==
          :existing_qw_polynomial_gaussian_axis_integrals
    _check_mixed_gto_summary_common(kinetic_summary)
    @test kinetic_block.dense_block[:, 1] ≈
          moment_oracle.kinetic_ga[oracle_rows, 1] atol = 1.0e-12 rtol = 1.0e-12

    unsupported_axis_block =
        CBPMixedGTO.cpb_mixed_gto_position_operator_block(
            parent,
            cpb,
            orbital;
            axis = :q,
        )
    _check_blocked_mixed_gto(
        unsupported_axis_block,
        :unsupported_mixed_gto_position_axis;
        expected_status = :blocked_cpb_mixed_gto_coordinate_moment_local_block,
    )

    unsupported_x2_axis_block =
        CBPMixedGTO.cpb_mixed_gto_x2_operator_block(
            parent,
            cpb,
            orbital;
            axis = :q,
        )
    _check_blocked_mixed_gto(
        unsupported_x2_axis_block,
        :unsupported_mixed_gto_x2_axis;
        expected_status = :blocked_cpb_mixed_gto_coordinate_moment_local_block,
    )

    unsupported_normalization = _mixed_gto_orbital(
        primitive_normalization = :unsupported_normalization,
    )
    _check_blocked_mixed_gto(
        CBPMixedGTO.cpb_mixed_gto_overlap_block(
            parent,
            cpb,
            unsupported_normalization,
        ),
        :unsupported_supplement_primitive_normalization,
    )
    _check_mixed_gto_one_body_wrapper_blockers(
        parent,
        cpb,
        unsupported_normalization,
        :unsupported_supplement_primitive_normalization,
    )

    invalid_powers = _mixed_gto_orbital(angular_powers = (-1, 0, 0))
    _check_blocked_mixed_gto(
        CBPMixedGTO.cpb_mixed_gto_overlap_block(parent, cpb, invalid_powers),
        :invalid_supplement_angular_powers,
    )

    missing_primitives = _mixed_gto_orbital(exponents = Float64[], coefficients = Float64[])
    _check_blocked_mixed_gto(
        CBPMixedGTO.cpb_mixed_gto_overlap_block(parent, cpb, missing_primitives),
        :missing_supplement_primitive_data,
    )
    _check_mixed_gto_one_body_wrapper_blockers(
        parent,
        cpb,
        missing_primitives,
        :missing_supplement_primitive_data,
    )

    mismatched_primitives = _mixed_gto_orbital(
        exponents = [0.7, 1.3],
        coefficients = [1.0],
    )
    _check_blocked_mixed_gto(
        CBPMixedGTO.cpb_mixed_gto_overlap_block(parent, cpb, mismatched_primitives),
        :supplement_primitive_count_mismatch,
    )
    _check_mixed_gto_one_body_wrapper_blockers(
        parent,
        cpb,
        mismatched_primitives,
        :supplement_primitive_count_mismatch,
    )

    invalid_exponents = _mixed_gto_orbital(exponents = [0.7, -1.3])
    _check_blocked_mixed_gto(
        CBPMixedGTO.cpb_mixed_gto_overlap_block(parent, cpb, invalid_exponents),
        :invalid_supplement_exponents,
    )

    invalid_coefficients = _mixed_gto_orbital(coefficients = [0.8, Inf])
    _check_blocked_mixed_gto(
        CBPMixedGTO.cpb_mixed_gto_overlap_block(parent, cpb, invalid_coefficients),
        :invalid_supplement_coefficients,
    )

    outside_cpb = CPBMixedGTO.cpb(1:2, 2:3, 3:4; role = :outside_mixed_gto_cpb)
    _check_blocked_mixed_gto(
        CBPMixedGTO.cpb_mixed_gto_overlap_block(parent, outside_cpb, orbital),
        :cpb_z_interval_outside_parent,
    )
    _check_mixed_gto_one_body_wrapper_blockers(
        parent,
        outside_cpb,
        orbital,
        :cpb_z_interval_outside_parent,
    )
end

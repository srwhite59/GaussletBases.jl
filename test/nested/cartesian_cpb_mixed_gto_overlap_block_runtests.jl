# Runtime role: CPB-local mixed gausslet/GTO overlap pilot.
#
# This test validates the tiny provider-level mixed gausslet/GTO overlap block
# against the existing Cartesian/GTO overlap oracle. It does not add kinetic,
# position, x2, nuclear/Coulomb terms, WL/PQS realization, retained transforms,
# route/global placement, driver wiring, Hamiltonian assembly, IDA/MWG,
# PQS Lowdin/projection, exports, or artifacts.

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

function _check_blocked_mixed_gto(block, expected_blocker)
    block_summary = CBPMixedGTO.summary(block)
    @test block_summary.status == :blocked_cpb_mixed_gto_local_overlap_block
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
    @test !hasproperty(block_summary, :global_matrix)
    @test !hasproperty(block_summary, :global_overlap_matrix)
    @test !hasproperty(block_summary, :retained_blocks)
    return nothing
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
    @test block_summary.local_shape == (x = 2, y = 2, z = 1)
    @test block_summary.dense_block_shape == (4, 1)
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
    @test !hasproperty(block_summary, :global_matrix)
    @test !hasproperty(block_summary, :global_overlap_matrix)
    @test !hasproperty(block_summary, :retained_blocks)
    @test block.dense_block[:, 1] ≈ oracle[oracle_rows, 1] atol = 1.0e-12 rtol = 1.0e-12

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

    mismatched_primitives = _mixed_gto_orbital(
        exponents = [0.7, 1.3],
        coefficients = [1.0],
    )
    _check_blocked_mixed_gto(
        CBPMixedGTO.cpb_mixed_gto_overlap_block(parent, cpb, mismatched_primitives),
        :supplement_primitive_count_mismatch,
    )

    outside_cpb = CPBMixedGTO.cpb(1:2, 2:3, 3:4; role = :outside_mixed_gto_cpb)
    _check_blocked_mixed_gto(
        CBPMixedGTO.cpb_mixed_gto_overlap_block(parent, outside_cpb, orbital),
        :cpb_z_interval_outside_parent,
    )
end

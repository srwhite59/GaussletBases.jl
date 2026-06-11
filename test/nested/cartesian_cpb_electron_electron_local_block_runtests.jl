# Runtime role: tiny CPB-local electron-electron pair-factor pilot test.
#
# This validates one provider-level CPB-local Coulomb block against the
# existing Cartesian/White-Lindsey direct-product interaction oracle. It does
# not wire route/global placement, Hamiltonian assembly, IDA/MWG/PQS semantics,
# PQS Lowdin/projection, electron-nuclear kernels, exports, or artifacts.

using Test
using GaussletBases

const CPBEEPilot = GaussletBases.CartesianCPB
const CPEEPilot = GaussletBases.CartesianCPBBlockProviders
const CPGEEPilot = GaussletBases.CartesianParentGaussletBases

function _ee_pilot_fixture()
    expansion = CoulombGaussianExpansion(
        [1.0, 0.5],
        [0.2, 0.7];
        del = 1.0,
        s = 0.16,
        c = 0.01,
        maxu = 135.0,
    )
    axis = build_basis(MappedUniformBasisSpec(
        :G10;
        count = 2,
        mapping = IdentityMapping(),
        reference_spacing = 1.0,
    ))
    bundle = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        axis;
        exponents = expansion.exponents,
    )
    parent = CPGEEPilot.CartesianParentGaussletBasis3D(axis)
    axis_bundle = (x = bundle, y = bundle, z = bundle)
    left_cpb = CPBEEPilot.cpb(
        1:2,
        1:1,
        1:2;
        role = :left_electron_electron_window,
    )
    right_cpb = CPBEEPilot.cpb(
        2:2,
        1:2,
        1:2;
        role = :right_electron_electron_window,
    )
    interval_pair = CPEEPilot.cpb_interval_pair(parent, left_cpb, right_cpb)
    return (;
        expansion,
        bundle,
        axis_bundle,
        interval_pair,
    )
end

function _ee_pilot_product_indices(intervals, dims)
    indices = Int[]
    for ix in intervals.x, iy in intervals.y, iz in intervals.z
        push!(indices, (ix - 1) * dims.y * dims.z + (iy - 1) * dims.z + iz)
    end
    return indices
end

@testset "CPB electron-electron local pair-factor interaction block pilot" begin
    fixture = _ee_pilot_fixture()
    interval_summary = CPEEPilot.summary(fixture.interval_pair)

    local_block = CPEEPilot.cpb_electron_electron_local_block(
        fixture.axis_bundle,
        fixture.expansion,
        fixture.interval_pair,
    )
    local_summary = CPEEPilot.summary(local_block)

    oracle = GaussletBases._qwrg_diatomic_interaction_matrix(
        fixture.bundle,
        fixture.bundle,
        fixture.bundle,
        fixture.expansion,
    )
    parent_dims = (x = 2, y = 2, z = 2)
    left_indices = _ee_pilot_product_indices(
        interval_summary.left_intervals,
        parent_dims,
    )
    right_indices = _ee_pilot_product_indices(
        interval_summary.right_intervals,
        parent_dims,
    )
    expected = oracle[left_indices, right_indices]

    @test local_summary.status ===
          :materialized_cpb_electron_electron_local_block
    @test isnothing(local_summary.blocker)
    @test local_summary.term === :electron_electron_pair_factor_interaction
    @test local_summary.source_kind === :parent_axis_bundle_pair_factor_terms
    @test local_summary.factor_source_path ===
          :axis_pgdg_intermediate_pair_factor_terms
    @test local_summary.representation ===
          :dense_local_cpb_electron_electron_pair_factor_interaction
    @test local_summary.pair_factor_source ===
          :axis_pgdg_intermediate_pair_factor_terms
    @test local_summary.pair_factor_weighting ===
          :pgdg_weight_divided_pair_factor_storage
    @test local_summary.pair_factor_normalization ===
          :source_weight_divided_storage
    @test !local_summary.density_normalized_pair_factors
    @test !local_summary.raw_weighted_pair_factors
    @test local_summary.source_weight_division_owner ===
          :pgdg_auxiliary_source_weights
    @test local_summary.source_weight_division_shape ===
          :axis_pair_weight_outer
    @test !local_summary.source_weight_division_applied_by_provider
    @test local_summary.source_weight_division_stage ===
          :parent_pgdg_intermediate_construction
    @test local_summary.axis_integral_weights_applied === false
    @test local_summary.axis_integral_weights_deferred === false
    @test local_summary.weight_application_stage ===
          :not_final_ida_density_interaction
    @test local_summary.retained_pqs_weights === false
    @test local_summary.ida_mwg_semantics === false
    @test local_summary.gaussian_expansion_loop === :inner_local_contraction
    @test local_summary.gaussian_term_count == length(fixture.expansion.coefficients)
    @test local_summary.left_shape == (x = 2, y = 1, z = 2)
    @test local_summary.right_shape == (x = 1, y = 2, z = 2)
    @test local_summary.dense_block_shape == (4, 4)
    @test local_summary.dense_block_eltype === Float64
    @test local_summary.factor_space === :parent_axis_bundle_pgdg_intermediate
    @test local_summary.factor_convention ===
          :axis_bundle_electron_electron_pair_factor_terms
    @test local_summary.local_ordering ===
          :parent_compatible_x_slowest_z_fastest
    @test local_summary.provider_level_local_matrix_materialized
    @test local_summary.provider_level_coulomb_block_materialized
    @test local_summary.provider_level_electron_electron_block_materialized
    @test local_summary.cpb_local_coulomb_kernel_implemented
    @test local_summary.cpb_local_electron_electron_kernel_implemented
    # The local CPB pilot consumes WL pair_factor_terms directly to match the
    # existing WL oracle. It is not the final retained IDA density-interaction
    # boundary, so retained density weights must not be applied here.
    @test local_block.dense_block ≈ expected
    @test !hasproperty(local_summary, :dense_block)
    @test !hasproperty(local_summary, :axis_ops)
    @test !hasproperty(local_summary, :global_overlap_matrix)
    @test local_summary.route_global_matrix_materialized === false
    @test local_summary.global_matrix_materialized === false
    @test local_summary.hamiltonian_assembly === false
    @test local_summary.hamiltonian_data_materialized === false
    @test local_summary.route_driver_wiring === false
    @test local_summary.exports_or_artifacts === false

    missing_expansion = CPEEPilot.cpb_electron_electron_local_block(
        fixture.axis_bundle,
        nothing,
        fixture.interval_pair,
    )
    missing_summary = CPEEPilot.summary(missing_expansion)

    @test missing_summary.status === :blocked_cpb_electron_electron_local_block
    @test missing_summary.blocker === :missing_coulomb_gaussian_expansion
    @test isnothing(missing_expansion.dense_block)
    @test missing_summary.provider_level_local_matrix_materialized === false
    @test missing_summary.route_global_matrix_materialized === false
    @test missing_summary.hamiltonian_assembly === false

    missing_pair_axis_bundle = merge(
        fixture.axis_bundle,
        (x = (; pgdg_intermediate = (;)),),
    )
    missing_pair = CPEEPilot.cpb_electron_electron_local_block(
        missing_pair_axis_bundle,
        fixture.expansion,
        fixture.interval_pair,
    )
    missing_pair_summary = CPEEPilot.summary(missing_pair)

    @test missing_pair_summary.status ===
          :blocked_cpb_electron_electron_local_block
    @test missing_pair_summary.blocker === :missing_x_pair_factor_terms
    @test isnothing(missing_pair.dense_block)
    @test missing_pair_summary.axis_integral_weights_applied === false
    @test missing_pair_summary.axis_integral_weights_deferred === false
    @test missing_pair_summary.pair_factor_normalization ===
          :source_weight_divided_storage
    @test missing_pair_summary.retained_pqs_weights === false
    @test missing_pair_summary.ida_mwg_semantics === false
    @test missing_pair_summary.provider_level_local_matrix_materialized === false
    @test missing_pair_summary.route_global_matrix_materialized === false
    @test missing_pair_summary.hamiltonian_assembly === false
end

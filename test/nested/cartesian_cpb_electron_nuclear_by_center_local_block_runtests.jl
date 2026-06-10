# Runtime role: tiny CPB-local electron-nuclear by-center Galerkin pilot test.
#
# This validates one provider-level CPB-local by-center nuclear attraction
# block against the existing Cartesian/White-Lindsey by-center oracle. It does
# not sum centers, build a Hamiltonian, apply CPB integral weights, route/global
# place matrices, add IDA/MWG/PQS semantics, or export artifacts.

using Test
using GaussletBases

const CPBENPilot = GaussletBases.CartesianCPB
const CPBProviderENPilot = GaussletBases.CartesianCPBBlockProviders
const CPGBENPilot = GaussletBases.CartesianParentGaussletBases

function _en_pilot_fixture()
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
        center = 0.0,
    )
    parent = CPGBENPilot.CartesianParentGaussletBasis3D(axis)
    axis_bundle = (x = bundle, y = bundle, z = bundle)
    left_cpb = CPBENPilot.cpb(
        1:2,
        1:1,
        1:2;
        role = :left_electron_nuclear_window,
    )
    right_cpb = CPBENPilot.cpb(
        2:2,
        1:2,
        1:2;
        role = :right_electron_nuclear_window,
    )
    interval_pair = CPBProviderENPilot.cpb_interval_pair(
        parent,
        left_cpb,
        right_cpb,
    )
    center_record = (;
        center_key = :nucleus_1,
        center_index = 1,
        nuclear_charge = 2.0,
        location = (0.0, 0.0, 0.0),
    )
    oracle_basis = BondAlignedDiatomicQWBasis3D(
        :x,
        axis,
        axis,
        axis,
        [(0.0, 0.0, 0.0)],
        [center_record.nuclear_charge],
        1.0,
    )
    return (;
        expansion,
        axis,
        bundle,
        axis_bundle,
        interval_pair,
        center_record,
        oracle_basis,
    )
end

function _en_pilot_product_indices(intervals, dims)
    indices = Int[]
    for ix in intervals.x, iy in intervals.y, iz in intervals.z
        push!(indices, (ix - 1) * dims.y * dims.z + (iy - 1) * dims.z + iz)
    end
    return indices
end

function _test_en_nonclaims(summary)
    @test summary.cpb_integral_weights_applied === false
    @test summary.ida_mwg_semantics === false
    @test summary.route_driver_wiring === false
    @test summary.route_global_matrix_materialized === false
    @test summary.global_matrix_materialized === false
    @test summary.hamiltonian_assembly === false
    @test summary.hamiltonian_data_materialized === false
    @test summary.exports_or_artifacts === false
end

@testset "CPB electron-nuclear by-center local Galerkin block pilot" begin
    fixture = _en_pilot_fixture()
    interval_summary = CPBProviderENPilot.summary(fixture.interval_pair)

    local_block =
        CPBProviderENPilot.cpb_electron_nuclear_by_center_local_block(
            fixture.axis_bundle,
            fixture.expansion,
            fixture.center_record,
            fixture.interval_pair,
        )
    local_summary = CPBProviderENPilot.summary(local_block)

    oracle_by_center = GaussletBases._qwrg_diatomic_nuclear_one_body_by_center(
        fixture.oracle_basis,
        fixture.bundle,
        fixture.bundle,
        fixture.bundle,
        fixture.expansion,
    )
    parent_dims = (x = 2, y = 2, z = 2)
    left_indices = _en_pilot_product_indices(
        interval_summary.left_intervals,
        parent_dims,
    )
    right_indices = _en_pilot_product_indices(
        interval_summary.right_intervals,
        parent_dims,
    )
    expected =
        fixture.center_record.nuclear_charge .*
        oracle_by_center[1][left_indices, right_indices]

    @test local_summary.status ===
          :materialized_cpb_electron_nuclear_by_center_local_block
    @test isnothing(local_summary.blocker)
    @test local_summary.term === :electron_nuclear_by_center_galerkin
    @test local_summary.source_kind ===
          :parent_axis_bundle_per_center_gaussian_factor_terms
    @test local_summary.center_key === :nucleus_1
    @test local_summary.center_index == 1
    @test local_summary.charge == 2.0
    @test local_summary.center_location == (0.0, 0.0, 0.0)
    @test local_summary.by_center
    @test local_summary.centers_summed === false
    @test local_summary.galerkin_operator
    @test local_summary.factor_source_path ===
          :axis_pgdg_intermediate_gaussian_factor_terms
    @test local_summary.gaussian_expansion_loop === :inner_local_contraction
    @test local_summary.gaussian_term_count == length(fixture.expansion.coefficients)
    @test local_summary.left_shape == (x = 2, y = 1, z = 2)
    @test local_summary.right_shape == (x = 1, y = 2, z = 2)
    @test local_summary.dense_block_shape == (4, 4)
    @test local_summary.dense_block_eltype === Float64
    @test local_summary.factor_space === :parent_axis_bundle_pgdg_intermediate
    @test local_summary.factor_convention ===
          :axis_bundle_electron_nuclear_per_center_gaussian_factor_terms
    @test local_summary.local_ordering ===
          :parent_compatible_x_slowest_z_fastest
    @test local_summary.provider_level_local_matrix_materialized
    @test local_summary.provider_level_coulomb_block_materialized
    @test local_summary.provider_level_electron_nuclear_block_materialized
    @test local_summary.cpb_local_coulomb_kernel_implemented
    @test local_summary.cpb_local_electron_nuclear_kernel_implemented
    @test local_block.dense_block ≈ expected
    @test !hasproperty(local_summary, :dense_block)
    @test !hasproperty(local_summary, :axis_ops)
    @test !hasproperty(local_summary, :global_overlap_matrix)
    _test_en_nonclaims(local_summary)

    missing_expansion =
        CPBProviderENPilot.cpb_electron_nuclear_by_center_local_block(
            fixture.axis_bundle,
            nothing,
            fixture.center_record,
            fixture.interval_pair,
        )
    missing_expansion_summary = CPBProviderENPilot.summary(missing_expansion)

    @test missing_expansion_summary.status ===
          :blocked_cpb_electron_nuclear_by_center_local_block
    @test missing_expansion_summary.blocker === :missing_coulomb_gaussian_expansion
    @test isnothing(missing_expansion.dense_block)
    @test missing_expansion_summary.provider_level_local_matrix_materialized === false
    _test_en_nonclaims(missing_expansion_summary)

    missing_center =
        CPBProviderENPilot.cpb_electron_nuclear_by_center_local_block(
            fixture.axis_bundle,
            fixture.expansion,
            (; center_key = :bad_center, nuclear_charge = 2.0),
            fixture.interval_pair,
        )
    missing_center_summary = CPBProviderENPilot.summary(missing_center)

    @test missing_center_summary.status ===
          :blocked_cpb_electron_nuclear_by_center_local_block
    @test missing_center_summary.blocker ===
          :missing_electron_nuclear_center_location
    @test isnothing(missing_center.dense_block)
    _test_en_nonclaims(missing_center_summary)

    missing_axis_bundle = merge(
        fixture.axis_bundle,
        (x = (; pgdg_intermediate = (;)),),
    )
    missing_terms =
        CPBProviderENPilot.cpb_electron_nuclear_by_center_local_block(
            missing_axis_bundle,
            fixture.expansion,
            fixture.center_record,
            fixture.interval_pair,
        )
    missing_terms_summary = CPBProviderENPilot.summary(missing_terms)

    @test missing_terms_summary.status ===
          :blocked_cpb_electron_nuclear_by_center_local_block
    @test missing_terms_summary.blocker ===
          :missing_x_electron_nuclear_axis_factor_terms
    @test isnothing(missing_terms.dense_block)
    _test_en_nonclaims(missing_terms_summary)
end

# Runtime role: metadata-only parent Coulomb source fingerprint test.
#
# This validates a compact source summary for future CPB-local Coulomb kernels.
# It does not build numerical Coulomb blocks, CPB kernels, WL/PQS realization,
# route/global placement, Hamiltonians, IDA/MWG/PQS semantics, exports, or
# artifacts.

using Test
using GaussletBases

const CPGBSourceTest = GaussletBases.CartesianParentGaussletBases

function _coulomb_source_test_parent()
    axis = build_basis(MappedUniformBasisSpec(
        :G10;
        count = 2,
        mapping = IdentityMapping(),
        reference_spacing = 1.0,
    ))
    return CPGBSourceTest.CartesianParentGaussletBasis3D(axis)
end

function _coulomb_source_test_axis_terms()
    return Array{Float64,3}(reshape(collect(1.0:8.0), 2, 2, 2))
end

function _coulomb_source_test_axis_bundle(; include_pair_terms::Bool = true)
    terms = _coulomb_source_test_axis_terms()
    function axis_bundle()
        pgdg = (;
            gaussian_factor_terms = terms,
            gaussian_factors = [[1.0 0.1; 0.1 1.2], [0.8 0.0; 0.0 0.9]],
            pair_factor_terms_raw = terms,
            exponents = [0.2, 0.7],
            center = 0.0,
            weights = [1.0, 1.0],
            centers = [-0.5, 0.5],
        )
        if include_pair_terms
            pgdg = merge(pgdg, (;
                pair_factor_terms = terms,
                pair_factors = [[1.1 0.2; 0.2 1.3], [0.7 0.1; 0.1 0.6]],
            ))
        end
        return (; pgdg_intermediate = pgdg)
    end
    return (x = axis_bundle(), y = axis_bundle(), z = axis_bundle())
end

function _coulomb_source_test_expansion()
    return CoulombGaussianExpansion(
        [1.0, 0.5],
        [0.2, 0.7];
        del = 1.0,
        s = 0.16,
        c = 0.01,
        maxu = 135.0,
    )
end

@testset "Cartesian parent Coulomb source summary" begin
    parent = _coulomb_source_test_parent()
    expansion = _coulomb_source_test_expansion()
    center_table = (;
        center_index = 1,
        center_key = :center_1,
        atom_symbol = "He",
        nuclear_charge = 2.0,
        location = (0.0, 0.0, 0.0),
    )

    source = CPGBSourceTest.parent_coulomb_axis_source_summary(
        parent,
        _coulomb_source_test_axis_bundle(),
        expansion;
        center_table = (center_table,),
    )
    source_summary = CPGBSourceTest.summary(source)

    @test source_summary.object_kind ===
          :cartesian_parent_coulomb_axis_source_summary
    @test source_summary.status ===
          :available_parent_coulomb_axis_source_summary
    @test isnothing(source_summary.blocker)
    @test source_summary.parent_axis_counts == (2, 2, 2)
    @test source_summary.axis_bundle_available
    @test source_summary.expansion_available
    @test source_summary.expansion_coefficients_available
    @test source_summary.expansion_exponents_available
    @test source_summary.gaussian_expansion_source ===
          :coulomb_gaussian_expansion_object
    @test source_summary.electron_nuclear_center_metadata_available
    @test source_summary.electron_nuclear_center_metadata_source ===
          :center_table
    @test source_summary.electron_nuclear_axis_terms_available === false
    @test source_summary.electron_nuclear_axis_terms_status ===
          :missing_per_center_electron_nuclear_axis_term_tables
    @test source_summary.electron_electron_pair_axis_terms_available
    @test source_summary.electron_electron_pair_axis_terms_status ===
          :available_parent_electron_electron_pair_axis_terms
    @test source_summary.electron_electron_cpb_pair_pair_source_record_available ===
          false
    @test source_summary.gaussian_factor_terms_available
    @test source_summary.gaussian_factors_available
    @test source_summary.pair_factor_terms_available
    @test source_summary.pair_factors_available
    @test source_summary.pair_factor_terms_raw_available
    @test source_summary.axis_exponents_available
    @test source_summary.source_paths.electron_nuclear_gaussian_factor_terms ===
          :axis_pgdg_intermediate_gaussian_factor_terms
    @test source_summary.source_paths.electron_electron_pair_factor_terms ===
          :axis_pgdg_intermediate_pair_factor_terms
    @test source_summary.source_paths.expansion_coefficients ===
          :expansion_coefficients
    @test source_summary.source_paths.expansion_exponents ===
          :expansion_exponents
    @test source_summary.source_paths.center_table === :center_table
    @test source_summary.numerical_coulomb_blocks_materialized === false
    @test source_summary.cpb_local_coulomb_kernel_implemented === false
    @test source_summary.wl_pqs_realization === false
    @test source_summary.route_global_placement === false
    @test source_summary.route_driver_wiring === false
    @test source_summary.hamiltonian_assembly === false
    @test source_summary.ida_mwg_pqs_semantics === false
    @test source_summary.exports_or_artifacts === false

    missing_expansion = CPGBSourceTest.parent_coulomb_axis_source_summary(
        parent,
        _coulomb_source_test_axis_bundle(),
        nothing;
        nuclear_charges = (2.0,),
        atom_locations = ((0.0, 0.0, 0.0),),
    )
    missing_expansion_summary = CPGBSourceTest.summary(missing_expansion)

    @test missing_expansion_summary.status ===
          :blocked_parent_coulomb_axis_source_summary
    @test missing_expansion_summary.blocker ===
          :missing_coulomb_gaussian_expansion
    @test :missing_coulomb_expansion_coefficients in
          missing_expansion_summary.missing_sources
    @test :missing_coulomb_expansion_exponents in
          missing_expansion_summary.missing_sources
    @test missing_expansion_summary.electron_nuclear_center_metadata_source ===
          :nuclear_charges_and_atom_locations
    @test missing_expansion_summary.electron_electron_pair_axis_terms_available
    @test missing_expansion_summary.numerical_coulomb_blocks_materialized === false

    missing_pair_terms = CPGBSourceTest.parent_coulomb_axis_source_summary(
        parent,
        _coulomb_source_test_axis_bundle(; include_pair_terms = false),
        expansion;
        parent_qw_basis_object = (;
            nuclei = ((0.0, 0.0, 0.0),),
            nuclear_charges = (2.0,),
        ),
    )
    missing_pair_terms_summary = CPGBSourceTest.summary(missing_pair_terms)

    @test missing_pair_terms_summary.status ===
          :blocked_parent_coulomb_axis_source_summary
    @test missing_pair_terms_summary.blocker ===
          :missing_electron_electron_pair_axis_terms
    @test missing_pair_terms_summary.electron_nuclear_center_metadata_source ===
          :parent_qw_basis_object_nuclei_and_charges
    @test missing_pair_terms_summary.gaussian_factor_terms_available
    @test missing_pair_terms_summary.pair_factor_terms_available === false
    @test missing_pair_terms_summary.pair_factors_available === false
    @test missing_pair_terms_summary.pair_factor_terms_raw_available
    @test missing_pair_terms_summary.electron_nuclear_axis_terms_available === false
    @test missing_pair_terms_summary.cpb_local_coulomb_kernel_implemented === false
end

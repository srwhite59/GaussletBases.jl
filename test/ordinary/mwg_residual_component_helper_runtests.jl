using Logging
using LinearAlgebra
using Test

using GaussletBases

@testset "Final-residual MWG component block helper" begin
    Z = 2.0
    expansion = GaussletBases.coulomb_gaussian_expansion(doacc = false)
    basis = GaussletBases.build_basis(GaussletBases.MappedUniformBasisSpec(
        :G10;
        count = 9,
        mapping = GaussletBases.white_lindsey_atomic_mapping(
            Z = Z,
            d = 0.2,
            tail_spacing = 10.0,
        ),
        reference_spacing = 1.0,
    ))
    supplement = GaussletBases.legacy_atomic_gaussian_supplement(
        "He",
        "cc-pVTZ";
        lmax = 0,
    )
    fixed_block = @test_logs min_level = Logging.Warn GaussletBases.one_center_atomic_full_parent_fixed_block(
        basis;
        expansion,
        nside = 3,
    )
    operators = @test_logs min_level = Logging.Warn GaussletBases.ordinary_cartesian_qiu_white_operators(
        fixed_block,
        supplement;
        expansion,
        Z = Z,
        interaction_treatment = :mwg,
    )
    fixed_interaction =
        GaussletBases._qwrg_fixed_block_interaction_matrix(fixed_block, expansion)
    bundle = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        fixed_block.parent_basis;
        exponents = expansion.exponents,
        center = 0.0,
        backend = fixed_block.gausslet_backend,
    )

    component_blocks = GaussletBases._qwrg_final_residual_mwg_component_blocks(
        fixed_interaction,
        fixed_block.coefficient_matrix,
        bundle,
        bundle,
        bundle,
        expansion,
        operators.residual_centers,
        operators.residual_widths,
    )
    final_interaction = GaussletBases._qwrg_fixed_block_interaction_matrix_mwg(
        fixed_interaction,
        fixed_block.coefficient_matrix,
        bundle,
        bundle,
        bundle,
        expansion,
        operators.residual_centers,
        operators.residual_widths,
    )

    nfixed = size(fixed_interaction, 1)
    nresidual = operators.residual_count
    @test nresidual > 0
    @test component_blocks.final_interaction == final_interaction
    @test component_blocks.final_interaction == operators.interaction_matrix
    @test component_blocks.fixed_fixed == fixed_interaction
    @test size(component_blocks.fixed_fixed) == (nfixed, nfixed)
    @test size(component_blocks.fixed_residual) == (nfixed, nresidual)
    @test size(component_blocks.residual_residual) == (nresidual, nresidual)
    @test size(component_blocks.final_interaction) ==
          (nfixed + nresidual, nfixed + nresidual)

    diagnostics = component_blocks.diagnostics
    @test diagnostics.fixed_fixed_block_source == :existing_fixed_gausslet_ida_path
    @test diagnostics.fixed_residual_component_source == :_qwrg_mwg_interaction_components
    @test diagnostics.residual_residual_component_source == :_qwrg_mwg_interaction_components
    @test diagnostics.residual_centers_finite
    @test diagnostics.residual_widths_finite
    @test diagnostics.residual_widths_positive
    @test diagnostics.raw_gto_rows_role == :residual_construction_inputs_only
    @test diagnostics.raw_gto_gto_mwg_interaction_blocks_used == false
    @test diagnostics.fixed_raw_gto_mwg_interaction_blocks_used == false
    @test diagnostics.raw_gto_interaction_blocks_used == false
    @test diagnostics.owner_semantics_source == :caller_residual_metadata
    @test diagnostics.owner_semantics_generalized_by_helper == false
    @test diagnostics.production_route_adoption == false
    @test diagnostics.packet_fixed_block_qw_hamiltonian_adoption == false
    @test diagnostics.public_default_route == false
    @test diagnostics.ecp == false
    @test diagnostics.be2_cr2_science_claim == false
    @test diagnostics.mwg_ida_semantic_change == false
end

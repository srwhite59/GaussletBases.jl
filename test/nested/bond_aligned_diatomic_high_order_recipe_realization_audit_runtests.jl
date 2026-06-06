@testset "Bond-aligned diatomic high-order recipe realization audit" begin
    plan29 = GaussletBases._nested_bond_aligned_diatomic_atom_growth_construction_plan(
        (1:29, 1:29, 1:29);
        bond_axis = :z,
        atom_axis_indices = (10, 20),
        protected_atom_side_count = 5,
    )
    policy = GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        plan29;
        q_min = 4,
    )
    audit =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_realization_audit(
            policy,
        )
    diagnostics =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_realization_diagnostics(
            audit,
        )

    @test audit isa GaussletBases._BondAlignedDiatomicHighOrderRecipeRealizationAudit3D
    @test diagnostics.recipe_label == :mixed_atom_cubic_shared_endcap_panel
    @test diagnostics.region_count == 9
    @test diagnostics.mapped_region_count == 9
    @test diagnostics.missing_region_count == 0
    @test diagnostics.buildable_without_mapped_primitive_count == 0
    @test diagnostics.active_builder_consumed_region_count == 0
    @test diagnostics.descriptor_region_count == 2
    @test diagnostics.exact_descriptor_region_count == 2
    @test diagnostics.ready_for_opt_in_builder == true
    @test diagnostics.active_builder_uses_policy == false
    @test diagnostics.support_coverage.coverage_ok

    realizations = diagnostics.region_realizations
    @test realizations[1].role == :outer_mismatch_shared_molecular_shell
    @test realizations[1].mapped_primitive ==
        :_nested_diatomic_high_order_outer_mismatch_descriptor
    @test realizations[1].mapped_primitive_status == :construction_piece_descriptor
    @test isnothing(realizations[1].missing_implementation)
    @test realizations[1].realization_descriptor.parent_support_count ==
        length(plan29.regions[1].support_indices)
    @test realizations[1].realization_descriptor.owned_unit_count == 4
    @test realizations[1].realization_descriptor.primitive_family ==
        :outer_mismatch_boundary_slab_set
    @test realizations[1].realization_descriptor.support_coverage.coverage_ok
    @test realizations[1].realization_descriptor.exact_full_coverage
    @test all(
        piece.primitive_family == :outer_mismatch_boundary_slab
        for piece in realizations[1].realization_descriptor.pieces
    )
    @test all(
        realization.mapped_primitive == :_nested_endcap_panel_shell_layer &&
        realization.mapped_primitive_status == :existing_internal_primitive &&
        isnothing(realization.missing_implementation) &&
        realization.existing_opt_in_route == :shared_shell_layer_policy_endcap_panel_owned
        for realization in realizations
        if realization.region_category == :shared_exterior
    )
    shared_realizations = [
        realization for realization in realizations
        if realization.region_category == :shared_exterior
    ]
    @test all(
        !isnothing(realization.metadata.realization_notes.projected_q_shell_candidate)
        for realization in shared_realizations
    )
    @test all(
        realization.metadata.realization_notes.projected_q_shell_candidate.mapped_primitive ==
        :_nested_projected_q_shell_layer
        for realization in shared_realizations
    )
    @test all(
        realization.metadata.realization_notes.projected_q_shell_candidate.primitive_family ==
        :projected_q_shell
        for realization in shared_realizations
    )
    @test all(
        realization.metadata.realization_notes.projected_q_shell_candidate.realization_primitive ==
        :projected_q_shell_boundary_comx_product_modes
        for realization in shared_realizations
    )
    @test all(
        realization.metadata.realization_notes.projected_q_shell_candidate.support_contract ==
        :projected_q_shell_raw_boundary
        for realization in shared_realizations
    )
    @test all(
        realization.metadata.realization_notes.projected_q_shell_candidate.coefficient_contract ==
        :full_block_boundary_comx_product_mode_projection
        for realization in shared_realizations
    )
    @test all(
        realization.metadata.realization_notes.projected_q_shell_candidate.cleanup_contract ==
        :full_rank_symmetric_lowdin
        for realization in shared_realizations
    )
    @test all(
        realization.metadata.realization_notes.projected_q_shell_candidate.mode_selection_rule ==
        :any_axis_mode_index_first_or_last
        for realization in shared_realizations
    )
    @test all(
        realization.metadata.realization_notes.projected_q_shell_candidate.sidecar_status ==
        :not_yet_optimized_product_staged_for_pqs
        for realization in shared_realizations
    )
    @test all(
        !realization.metadata.realization_notes.projected_q_shell_candidate.active_builder_consumes &&
        !realization.metadata.realization_notes.projected_q_shell_candidate.source_builder_consumes &&
        !realization.metadata.realization_notes.projected_q_shell_candidate.fixed_block_consumes &&
        !realization.metadata.realization_notes.projected_q_shell_candidate.qw_consumes &&
        !realization.metadata.realization_notes.projected_q_shell_candidate.hamiltonian_consumes
        for realization in shared_realizations
    )
    @test all(
        realization.metadata.realization_notes.projected_q_shell_candidate.current_transitional_implementation ==
        :_nested_endcap_panel_shell_layer
        for realization in shared_realizations
    )
    @test all(
        realization.mapped_primitive == :_nested_bond_aligned_diatomic_sequence_for_box &&
        realization.mapped_primitive_status == :existing_internal_primitive &&
        isnothing(realization.missing_implementation)
        for realization in realizations
        if realization.region_category == :atom_local
    )
    @test realizations[end].role == :contact_cap
    @test realizations[end].mapped_primitive ==
        :_nested_diatomic_high_order_contact_cap_descriptor
    @test realizations[end].mapped_primitive_status == :construction_piece_descriptor
    @test isnothing(realizations[end].missing_implementation)
    @test realizations[end].realization_descriptor.parent_support_count ==
        length(plan29.regions[end].support_indices)
    @test realizations[end].realization_descriptor.owned_unit_count == 1
    @test realizations[end].realization_descriptor.primitive_family == :contact_cap_owned_slab
    @test realizations[end].realization_descriptor.support_coverage.coverage_ok
    @test realizations[end].realization_descriptor.exact_full_coverage
    @test all(
        !(
            realization.buildability_status == :buildable_now &&
            isnothing(realization.mapped_primitive)
        )
        for realization in realizations
    )
    @test all(!realization.active_builder_consumes for realization in realizations)

    annulus_policy = GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        plan29;
        q_min = 4,
        atom_q = 4,
        shared_q = 5,
        contact_q = 4,
        outer_mismatch_q = 4,
        shared_exterior_family = :transverse_annulus_exterior,
    )
    annulus_audit =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_realization_audit(
            annulus_policy,
        )
    annulus_diagnostics =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_realization_diagnostics(
            annulus_audit,
        )
    @test annulus_diagnostics.recipe_label == :mixed_atom_cubic_shared_transverse_annulus
    @test annulus_diagnostics.mapped_region_count == 4
    @test annulus_diagnostics.missing_region_count == 5
    @test annulus_diagnostics.buildable_without_mapped_primitive_count == 0
    @test annulus_diagnostics.descriptor_region_count == 2
    @test annulus_diagnostics.exact_descriptor_region_count == 2
    @test annulus_diagnostics.ready_for_opt_in_builder == false
    @test all(
        realization.recipe_family == :transverse_annulus_exterior &&
        realization.q == 5 &&
        realization.mapped_primitive_status == :missing_experimental_primitive &&
        realization.missing_implementation == :transverse_annulus_owned_unit_producer &&
        !(:projected_q_shell_candidate in propertynames(realization.metadata.realization_notes))
        for realization in annulus_diagnostics.region_realizations
        if realization.region_category == :shared_exterior
    )
end

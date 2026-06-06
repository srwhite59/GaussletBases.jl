@testset "Bond-aligned diatomic high-order recipe policy metadata" begin
    plan29 = GaussletBases._nested_bond_aligned_diatomic_atom_growth_construction_plan(
        (1:29, 1:29, 1:29);
        bond_axis = :z,
        atom_axis_indices = (10, 20),
        protected_atom_side_count = 5,
    )
    policy = GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        plan29;
        q_min = 4,
        atom_q = 4,
        shared_q = 4,
    )
    diagnostics =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy_diagnostics(
            policy,
        )

    @test policy isa GaussletBases._BondAlignedDiatomicHighOrderRecipePolicy3D
    @test policy.recipe_label == :mixed_atom_cubic_shared_endcap_panel
    @test policy.q_min == 4
    @test policy.q_region_counts == Dict(4 => 9)
    @test policy.buildable_region_count == 7
    @test policy.planned_region_count == 2
    @test policy.experimental_region_count == 0
    @test policy.metadata.active_builder_uses_policy == false
    @test policy.metadata.source_builder_changed == false
    @test diagnostics.region_count == length(plan29.regions)
    @test diagnostics.active_builder_uses_policy == false
    @test diagnostics.support_coverage.coverage_ok
    @test diagnostics.outer_mismatch_is_outermost
    @test diagnostics.middle_contact_clean

    roles = [choice.region_role for choice in policy.region_choices]
    families = [choice.recipe_family for choice in policy.region_choices]
    categories = [choice.region_category for choice in policy.region_choices]
    statuses = [choice.buildability_status for choice in policy.region_choices]
    support_counts = [choice.support_count for choice in policy.region_choices]
    @test roles == plan29.region_order
    @test families[1] == :outermost_mismatch_shared_molecular_shell
    @test all(==(:shared_endcap_panel_exterior), families[2:6])
    @test all(==(:protected_atom_cubic_shell), families[7:8])
    @test families[9] == :shared_contact_cap
    @test categories[1] == :outer_mismatch
    @test all(==(:shared_exterior), categories[2:6])
    @test all(==(:atom_local), categories[7:8])
    @test categories[9] == :contact_cap
    @test statuses[1] == :planned_only
    @test all(==(:buildable_now), statuses[2:8])
    @test statuses[9] == :planned_only
    @test support_counts == [length(region.support_indices) for region in plan29.regions]
    @test all(choice.q == 4 && choice.order == 4 for choice in policy.region_choices)

    annulus_policy = GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        plan29;
        q_min = 4,
        atom_q = 4,
        shared_q = 5,
        contact_q = 4,
        outer_mismatch_q = 4,
        shared_exterior_family = :transverse_annulus_exterior,
    )
    annulus_diagnostics =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy_diagnostics(
            annulus_policy,
        )
    @test annulus_policy.recipe_label == :mixed_atom_cubic_shared_transverse_annulus
    @test annulus_policy.q_region_counts == Dict(4 => 4, 5 => 5)
    @test annulus_policy.buildable_region_count == 2
    @test annulus_policy.planned_region_count == 2
    @test annulus_policy.experimental_region_count == 5
    @test all(
        choice.recipe_family == :transverse_annulus_exterior &&
        choice.q == 5 &&
        choice.implementation_status == :promising_experimental &&
        choice.buildability_status == :planned_experimental
        for choice in annulus_policy.region_choices
        if choice.region_category == :shared_exterior
    )
    @test annulus_diagnostics.region_choices[2].recipe_family ==
        :transverse_annulus_exterior
    @test annulus_diagnostics.region_choices[2].q == 5

    @test_throws ArgumentError GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        plan29;
        q_min = 5,
        atom_q = 4,
    )
end

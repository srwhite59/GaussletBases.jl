@testset "Bond-aligned diatomic atom-growth construction plan" begin
    parent29 = (1:29, 1:29, 1:29)
    anatomy29 = GaussletBases._nested_bond_aligned_diatomic_atom_growth_anatomy(
        parent29;
        bond_axis = :z,
        atom_axis_indices = (10, 20),
        protected_atom_side_count = 5,
    )
    plan29 = GaussletBases._nested_bond_aligned_diatomic_atom_growth_construction_plan(
        anatomy29,
    )

    @test plan29 isa GaussletBases._BondAlignedDiatomicAtomGrowthConstructionPlan3D
    @test plan29.region_order == [
        :outer_mismatch_shared_molecular_shell,
        :regular_shared_molecular_shell,
        :regular_shared_molecular_shell,
        :regular_shared_molecular_shell,
        :regular_shared_molecular_shell,
        :regular_shared_molecular_shell,
        :left_atom_box,
        :right_atom_box,
        :contact_cap,
    ]
    @test [region.order_index for region in plan29.regions] == collect(1:9)
    @test plan29.outer_mismatch_is_outermost
    @test plan29.middle_contact_clean
    @test plan29.support_coverage.expected_support_count == 29^3
    @test plan29.support_coverage.region_support_count == 29^3
    @test plan29.support_coverage.covered_support_count == 29^3
    @test plan29.support_coverage.duplicate_count == 0
    @test plan29.support_coverage.missing_count == 0
    @test plan29.support_coverage.outside_count == 0
    @test plan29.support_coverage.coverage_ok
    @test length(plan29.regions[1].support_indices) == 29^3 - 19 * 19 * 29
    @test [
        length(region.support_indices)
        for region in plan29.regions
        if region.role == :regular_shared_molecular_shell
    ] == [2666, 2178, 1738, 1346, 1002]
    @test length(plan29.regions[end - 2].support_indices) == 9^3
    @test length(plan29.regions[end - 1].support_indices) == 9^3
    @test length(plan29.regions[end].support_indices) == 9 * 9
    @test plan29.regions[1].metadata.mismatch_absorption_policy ==
        :outermost_shared_molecular_shell
    @test plan29.regions[2].metadata.shell_offset == 5
    @test plan29.regions[6].metadata.shell_offset == 1
    @test plan29.regions[end].metadata.contact_policy == :single_shared_contact_cap
    CCP = GaussletBases.CartesianContractedParents
    atom_region = CCP.cartesian_shell_region(
        plan29.regions[end - 2];
        parent_dimension = 29^3,
    )
    shared_region = CCP.cartesian_shell_region(
        plan29.regions[2];
        parent_dimension = 29^3,
    )
    contact_region = CCP.cartesian_shell_region(
        plan29.regions[end];
        parent_dimension = 29^3,
    )
    mismatch_region = CCP.cartesian_shell_region(
        plan29.regions[1];
        parent_dimension = 29^3,
    )
    @test atom_region isa CCP.CartesianShellRegion3D
    @test atom_region.region_family == :atom_core_cube
    @test atom_region.status == :clean
    @test atom_region.role == :left_atom_box
    @test atom_region.support_summary.entry_count == 9^3
    @test atom_region.support_summary.outside_count == 0
    @test atom_region.ownership_coverage_contract == :disjoint_partition_piece
    @test atom_region.retention.retention_rule == :protected_atom_cubic_shell
    @test atom_region.retention.preferred_contraction_rule == :complete_shell_sequence
    @test atom_region.retention.metric_capability == :support_local_product
    @test :coefficient_matrix in atom_region.retention.missing_payload_fields
    @test !atom_region.current_route_consumes
    @test !atom_region.descriptor_drives_builder
    @test atom_region.descriptor_only
    @test shared_region.region_family == :rectangular_molecular_shell
    @test shared_region.retention.retention_rule == :policy_selected_shared_exterior
    @test shared_region.retention.metric_capability == :metadata_only_policy_dependent
    @test contact_region.region_family == :shared_midpoint_slab_cap
    @test contact_region.status == :clean
    @test contact_region.retention.preferred_contraction_rule == :contact_cap_owned_slab
    @test !contact_region.current_route_consumes
    @test mismatch_region.region_family == :outer_boundary_shell
    @test mismatch_region.retention.retention_rule ==
          :outermost_mismatch_shared_molecular_shell
    @test mismatch_region.retention.preferred_contraction_rule ==
          :outer_mismatch_boundary_slab_set
    @test mismatch_region.geometry.mismatch_low_counts == (5, 5, 0)
    @test !mismatch_region.current_route_consumes

    even_plan = GaussletBases._nested_bond_aligned_diatomic_atom_growth_construction_plan(
        (1:20, 1:20, 1:20);
        bond_axis = :z,
        atom_axis_indices = (7, 14),
        protected_atom_side_count = 4,
    )
    @test even_plan.support_coverage.coverage_ok
    @test even_plan.middle_contact_clean
    @test !(:contact_cap in even_plan.region_order)
    @test even_plan.region_order[end - 1:end] == [:left_atom_box, :right_atom_box]

    expansion = coulomb_gaussian_expansion(doacc = false)
    real_basis = bond_aligned_homonuclear_qw_basis(
        family = :G10,
        bond_length = 5.0,
        core_spacing = 0.7,
        xmax_parallel = 8.0,
        xmax_transverse = 4.0,
        bond_axis = :z,
    )
    real_bundles = GaussletBases._qwrg_bond_aligned_axis_bundles(real_basis, expansion)
    real_anatomy = GaussletBases._nested_bond_aligned_diatomic_atom_growth_anatomy(
        real_basis,
        real_bundles;
        protected_atom_side_count = 5,
    )
    real_plan = GaussletBases._nested_bond_aligned_diatomic_atom_growth_construction_plan(
        real_anatomy,
    )
    diagnostic_anatomy = bond_aligned_diatomic_nested_geometry_diagnostics(
        real_basis;
        nside = 5,
    ).atom_growth_anatomy.anatomy
    diagnostic_plan = GaussletBases._nested_bond_aligned_diatomic_atom_growth_construction_plan(
        diagnostic_anatomy,
    )
    @test real_plan.region_order == [
        :outer_mismatch_shared_molecular_shell,
        :regular_shared_molecular_shell,
        :left_atom_box,
        :right_atom_box,
        :contact_cap,
    ]
    @test [length(region.support_indices) for region in real_plan.regions] ==
        [98, 362, 125, 125, 25]
    @test real_plan.support_coverage.expected_support_count == 7 * 7 * 15
    @test real_plan.support_coverage.region_support_count == 7 * 7 * 15
    @test real_plan.support_coverage.covered_support_count == 7 * 7 * 15
    @test real_plan.support_coverage.duplicate_count == 0
    @test real_plan.support_coverage.missing_count == 0
    @test real_plan.support_coverage.outside_count == 0
    @test real_plan.support_coverage.coverage_ok
    @test real_plan.outer_mismatch_is_outermost
    @test real_plan.middle_contact_clean
    @test real_plan.regions[1].metadata.low_counts == (0, 0, 1)
    @test real_plan.regions[1].metadata.high_counts == (0, 0, 1)
    @test diagnostic_plan.region_order == real_plan.region_order
    @test diagnostic_plan.support_coverage.coverage_ok
end

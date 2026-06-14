# Integration/slow test. Do not include in default nested runner.

@testset "Bond-aligned diatomic high-order recipe opt-in source construction" begin
    basis = bond_aligned_homonuclear_qw_basis(
        family = :G10,
        bond_length = 5.0,
        core_spacing = 0.7,
        xmax_parallel = 8.0,
        xmax_transverse = 4.0,
        bond_axis = :z,
    )
    expansion = coulomb_gaussian_expansion(doacc = false)
    bundles = GaussletBases._qwrg_bond_aligned_axis_bundles(basis, expansion)
    policy = GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        basis,
        bundles;
        protected_atom_side_count = 5,
        q_min = 4,
    )
    realization_diagnostics =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_realization_diagnostics(
            GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_realization_audit(
                policy,
            ),
        )
    @test realization_diagnostics.ready_for_opt_in_builder
    @test !realization_diagnostics.active_builder_uses_policy
    @test realization_diagnostics.active_builder_consumed_region_count == 0

    construction =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction(
            basis,
            bundles,
            policy;
            nside = 5,
            term_coefficients = Float64.(expansion.coefficients),
            packet_kernel = :factorized_direct,
            build_sequence_packet = false,
        )
    diagnostics =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction_diagnostics(
            construction,
        )

    @test construction isa
          GaussletBases._BondAlignedDiatomicHighOrderRecipeSourceConstruction3D
    @test diagnostics.recipe_label == :mixed_atom_cubic_shared_endcap_panel
    @test diagnostics.active_builder_consumes
    @test diagnostics.active_builder_uses_policy
    @test diagnostics.metadata.default_source_builder_changed == false
    @test diagnostics.consumed_region_count == diagnostics.region_count == 5
    @test diagnostics.unsupported_region_count == 0
    @test diagnostics.parent_dimension == 7 * 7 * 15
    @test diagnostics.fixed_dimension == 469
    @test diagnostics.support_coverage.expected_support_count == 7 * 7 * 15
    @test diagnostics.support_coverage.coverage_ok
    @test isnothing(construction.sequence.packet)

    region_builds = diagnostics.region_builds
    @test [build.role for build in region_builds] == [
        :outer_mismatch_shared_molecular_shell,
        :regular_shared_molecular_shell,
        :left_atom_box,
        :right_atom_box,
        :contact_cap,
    ]
    @test [build.primitive_family for build in region_builds] == [
        :outer_mismatch_boundary_slab_set,
        :shared_endcap_panel_shell_layer,
        :atom_local_complete_shell_sequence,
        :atom_local_complete_shell_sequence,
        :contact_cap_owned_slab,
    ]
    @test [build.built_support_count for build in region_builds] == [98, 362, 125, 125, 25]
    @test [build.retained_count for build in region_builds] == [98, 96, 125, 125, 25]
    @test [build.column_range for build in region_builds] ==
          [1:98, 374:469, 99:223, 224:348, 349:373]
    @test all(build.built && build.active_builder_consumes for build in region_builds)
    @test all(build.support_coverage.coverage_ok for build in region_builds)
    @test region_builds[2].metadata.support_contract == :thin_endcap_box_perimeter
    @test region_builds[2].metadata.coefficient_contract == :product_doside
    @test region_builds[5].metadata.descriptor_scope == :middle_contact_cap
    CCP = GaussletBases.CartesianContractedParents
    inventory = CCP.cartesian_shell_region_inventory(
        construction;
        parent_dimension = diagnostics.parent_dimension,
    )
    @test inventory isa CCP.CartesianShellRegionInventory3D
    @test inventory.region_count == diagnostics.region_count == length(region_builds)
    @test inventory.region_order == [build.role for build in region_builds]
    @test [region.role for region in inventory.regions] == inventory.region_order
    @test [region.status for region in inventory.regions] ==
          [:clean, :transitional, :clean, :clean, :clean]
    @test [region.retention.preferred_contraction_rule for region in inventory.regions] == [
        :outer_mismatch_boundary_slab_set,
        :old_endcap_panel_product_split,
        :complete_shell_sequence,
        :complete_shell_sequence,
        :contact_cap_owned_slab,
    ]
    @test inventory.status_counts.clean == 4
    @test inventory.status_counts.transitional == 1
    @test inventory.current_route_consumes_count == diagnostics.region_count
    @test inventory.descriptor_only_count == 0
    @test inventory.support_summary.region_support_entry_count ==
          diagnostics.support_coverage.covered_support_count
    @test inventory.support_summary.support_complete_by_region_counts
    @test inventory.support_summary.count_only_summaries_for_all_regions
    @test all(isnothing(region.support_indices) for region in inventory.regions)
    @test all(region.current_route_consumes for region in inventory.regions)
    @test all(!region.descriptor_drives_builder for region in inventory.regions)
    @test all(!region.descriptor_only for region in inventory.regions)
    @test inventory.regions[2].ownership_coverage_contract == :boundary_only
    @test inventory.regions[2].retention.metric_capability ==
          :product_staged_metric_contraction
    @test isempty(inventory.regions[2].retention.missing_payload_fields)
    @test inventory.regions[5].ownership_coverage_contract == :disjoint_partition_piece
    @test isempty(inventory.regions[5].retention.missing_payload_fields)
    @test diagnostics.metadata.q_policy == :atom_growth_endcap_panel_shared_q_variable
    @test diagnostics.metadata.q4_acceptance_fixture
    @test diagnostics.metadata.non_shared_q_policy == :fixed_q4_order4
    @test diagnostics.metadata.shared_q_values == (4,)
    @test diagnostics.metadata.shared_order_values == (4,)
    @test diagnostics.metadata.shared_shell_realization == :endcap_panel_owned
    @test !diagnostics.metadata.projected_q_shell_opt_in

    be2_basis = bond_aligned_homonuclear_qw_basis(
        family = :G10,
        bond_length = 5.0,
        core_spacing = 0.15,
        xmax_parallel = 10.5,
        xmax_transverse = 8.0,
        bond_axis = :z,
        nuclear_charge = 4.0,
    )
    be2_bundles = GaussletBases._qwrg_bond_aligned_axis_bundles(be2_basis, expansion)
    @test GaussletBases._nested_axis_lengths(be2_bundles) == (15, 15, 27)
    be2_policy0 = GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        be2_basis,
        be2_bundles;
        protected_atom_side_count = 5,
        q_min = 4,
    )
    be2_policy = GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        be2_policy0.construction_plan;
        q_min = 4,
        atom_q = 4,
        atom_order = 4,
        shared_q = 5,
        shared_order = 5,
        contact_q = 4,
        contact_order = 4,
        outer_mismatch_q = 4,
        outer_mismatch_order = 4,
    )
    be2_retention = GaussletBases._nested_resolve_complete_shell_retention(5)
    be2_shared_dimensions = [
        GaussletBases._nested_diatomic_projected_q_shell_adaptive_source_dimensions(
            be2_basis,
            be2_bundles,
            region,
            be2_retention;
            bond_axis = :z,
            nside = 5,
            selected_q = choice.q,
            shared_shell_angular_resolution_scale = 1.4,
        )
        for (choice, region) in zip(
            be2_policy.region_choices,
            be2_policy.construction_plan.regions,
        )
        if choice.recipe_family == :shared_endcap_panel_exterior
    ]
    be2_source_box_dimension_plans = [
        GaussletBases._nested_diatomic_source_box_dimension_plan(
            be2_basis,
            be2_bundles,
            region.box,
            something(region.inner_exclusion_box),
            be2_retention;
            bond_axis = :z,
            nside = 5,
            selected_q = choice.q,
            shared_shell_angular_resolution_scale = 1.4,
            support_count = length(region.support_indices),
        )
        for (choice, region) in zip(
            be2_policy.region_choices,
            be2_policy.construction_plan.regions,
        )
        if choice.recipe_family == :shared_endcap_panel_exterior
    ]
    be2_shared_regions = [
        region
        for (choice, region) in zip(
            be2_policy.region_choices,
            be2_policy.construction_plan.regions,
        )
        if choice.recipe_family == :shared_endcap_panel_exterior
    ]
    pqs_plan = GaussletBases._nested_projected_q_shell_source_mode_plan(
        (5, 5, 6);
        bond_axis = :z,
        selected_q = 5,
        physical_box_lengths = (11, 11, 21),
        support_count = 1002,
    )
    @test pqs_plan.raw_source_dims == (5, 5, 6)
    @test pqs_plan.source_mode_dims == (5, 5, 6)
    @test pqs_plan.axis_selector_retained_counts == (3, 3, 4)
    @test pqs_plan.raw_q == 5
    @test pqs_plan.raw_L == 6
    @test pqs_plan.raw_q_matches_selected_q
    @test pqs_plan.physical_box_lengths == (11, 11, 21)
    @test pqs_plan.support_count == 1002
    @test pqs_plan.pqs_retained_count == 114
    @test pqs_plan.decomposition_status == :adaptive_broad_support_q_local_modes
    @test !pqs_plan.broad_parent_boundary_reference
    @test !pqs_plan.excluded_from_mvp_gate
    pqs_mismatch_plan = GaussletBases._nested_projected_q_shell_source_mode_plan(
        (4, 4, 5);
        bond_axis = :z,
        selected_q = 5,
        physical_box_lengths = (9, 9, 9),
        support_count = 488,
    )
    @test !pqs_mismatch_plan.raw_q_matches_selected_q
    @test pqs_mismatch_plan.decomposition_status == :adaptive_raw_q_mismatch
    @test pqs_mismatch_plan.excluded_from_mvp_gate
    @test !(:adaptive_retain in propertynames(pqs_plan))
    @test [length.(region.box) for region in be2_shared_regions] ==
          [(15, 15, 25), (13, 13, 23), (11, 11, 21)]
    @test [dims.raw_source_dims for dims in be2_shared_dimensions] ==
          [(5, 5, 5), (5, 5, 5), (5, 5, 6)]
    @test [plan.dimension_policy for plan in be2_source_box_dimension_plans] ==
          fill(:diatomic_adaptive_angular_source_box, 3)
    @test [plan.source_mode_dims for plan in be2_source_box_dimension_plans] ==
          [(5, 5, 5), (5, 5, 5), (5, 5, 6)]
    @test [plan.side_dimensions for plan in be2_source_box_dimension_plans] ==
          [(5, 5, 5), (5, 5, 5), (5, 5, 6)]
    @test [plan.raw_L for plan in be2_source_box_dimension_plans] == [5, 5, 6]
    @test all(plan -> plan.total_source_dimensions_primary, be2_source_box_dimension_plans)
    @test all(
        plan -> plan.diagnostics.source_mode_dims_are_total_lengths,
        be2_source_box_dimension_plans,
    )
    @test all(
        plan -> plan.diagnostics.axis_selector_retained_counts_are_diagnostic,
        be2_source_box_dimension_plans,
    )
    @test [dims.pqs_retained_count for dims in be2_shared_dimensions] ==
          [98, 98, 114]
    @test [
        dims.source_box_dimension_plan.source_mode_dims
        for dims in be2_shared_dimensions
    ] == [(5, 5, 5), (5, 5, 5), (5, 5, 6)]
    @test all(
        dims -> !(:adaptive_retain in propertynames(dims)),
        be2_shared_dimensions,
    )
    @test all(dims -> dims.raw_q == 5, be2_shared_dimensions)
    @test all(dims -> dims.raw_q_matches_selected_q, be2_shared_dimensions)
    @test all(
        dims -> dims.decomposition_status == :adaptive_broad_support_q_local_modes,
        be2_shared_dimensions,
    )
    @test all(dims -> !dims.broad_parent_boundary_reference, be2_shared_dimensions)
    @test all(dims -> !dims.excluded_from_mvp_gate, be2_shared_dimensions)
    @test [
        GaussletBases._nested_diatomic_projected_q_shell_retained_count(length.(region.box))
        for region in be2_shared_regions
    ] == [1738, 1346, 1002]
    @test [
        dims.pqs_retained_count !=
        GaussletBases._nested_diatomic_projected_q_shell_retained_count(length.(region.box))
        for (dims, region) in zip(be2_shared_dimensions, be2_shared_regions)
    ] == [true, true, true]

    pqs_construction =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction(
            basis,
            bundles,
            policy;
            nside = 5,
            term_coefficients = Float64.(expansion.coefficients),
            packet_kernel = :support_reference,
            shared_shell_realization = :projected_q_shell,
        )
    pqs_diagnostics =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction_diagnostics(
            pqs_construction,
        )
    @test pqs_diagnostics.metadata.shared_shell_realization == :projected_q_shell
    @test pqs_diagnostics.metadata.projected_q_shell_opt_in
    @test pqs_diagnostics.metadata.default_source_builder_changed == false
    @test pqs_diagnostics.metadata.packet_kernel == :support_reference
    @test pqs_diagnostics.metadata.build_sequence_packet
    @test pqs_diagnostics.support_coverage.coverage_ok
    @test pqs_diagnostics.support_coverage.expected_support_count == 7 * 7 * 15
    @test pqs_diagnostics.support_coverage.covered_support_count == 7 * 7 * 15
    @test pqs_diagnostics.support_coverage.duplicate_count == 0
    @test pqs_diagnostics.support_coverage.missing_count == 0
    @test pqs_diagnostics.support_coverage.outside_count == 0
    @test [build.role for build in pqs_diagnostics.region_builds] ==
          [build.role for build in region_builds]
    @test [build.primitive_family for build in pqs_diagnostics.region_builds] == [
        :outer_mismatch_boundary_slab_set,
        :projected_q_shell,
        :atom_local_complete_shell_sequence,
        :atom_local_complete_shell_sequence,
        :contact_cap_owned_slab,
    ]
    @test [
        build.built_support_count for build in pqs_diagnostics.region_builds
    ] == [build.built_support_count for build in region_builds]
    @test pqs_diagnostics.region_builds[2].mapped_primitive ==
          :_nested_projected_q_shell_layer
    @test pqs_diagnostics.region_builds[2].metadata.support_contract ==
          :projected_q_shell_raw_boundary
    @test pqs_diagnostics.region_builds[2].metadata.coefficient_contract ==
          :full_block_boundary_comx_product_mode_projection
    @test pqs_diagnostics.region_builds[2].metadata.seed_contract ==
          :raw_boundary_projection_of_boundary_comx_product_modes_from_full_local_block_transform
    @test pqs_diagnostics.region_builds[2].metadata.cleanup_contract ==
          :full_rank_symmetric_lowdin
    @test pqs_diagnostics.region_builds[2].metadata.cleanup_method ==
          :projected_boundary_symmetric_lowdin
    @test !pqs_diagnostics.region_builds[2].metadata.pqs_product_staged_sidecar_available
    @test !pqs_diagnostics.region_builds[2].metadata.factorized_direct_allowed
    @test !pqs_diagnostics.region_builds[2].metadata.active_default_builder_changed
    @test pqs_diagnostics.region_builds[2].metadata.selected_q == 4
    @test pqs_diagnostics.region_builds[2].metadata.raw_source_dims == (5, 5, 6)
    @test pqs_diagnostics.region_builds[2].metadata.source_mode_dims == (5, 5, 6)
    @test !(
        :adaptive_retain in
        propertynames(pqs_diagnostics.region_builds[2].metadata)
    )
    @test pqs_diagnostics.region_builds[2].metadata.raw_q == 5
    @test pqs_diagnostics.region_builds[2].metadata.raw_L == 6
    @test !pqs_diagnostics.region_builds[2].metadata.raw_q_matches_selected_q
    @test pqs_diagnostics.region_builds[2].metadata.physical_box_lengths == (7, 7, 13)
    @test pqs_diagnostics.region_builds[2].metadata.support_count == 362
    @test pqs_diagnostics.region_builds[2].metadata.pqs_retained_count == 114
    @test pqs_diagnostics.region_builds[2].metadata.decomposition_status ==
          :adaptive_raw_q_mismatch
    @test !pqs_diagnostics.region_builds[2].metadata.broad_parent_boundary_reference
    @test pqs_diagnostics.region_builds[2].metadata.excluded_from_mvp_gate
    @test pqs_diagnostics.region_builds[2].metadata.policy_q == 4
    @test pqs_diagnostics.region_builds[2].metadata.policy_order == 4
    pqs_source_descriptor =
        pqs_diagnostics.region_builds[2].metadata.pqs_staged_unit_descriptor
    @test pqs_source_descriptor.kind == :projected_q_shell
    @test length.(pqs_source_descriptor.current_box) == (7, 7, 13)
    @test all(
        axis ->
            first(pqs_source_descriptor.inner_box[axis]) ==
            first(pqs_source_descriptor.current_box[axis]) + 1 &&
            last(pqs_source_descriptor.inner_box[axis]) ==
            last(pqs_source_descriptor.current_box[axis]) - 1,
        1:3,
    )
    @test pqs_source_descriptor.bond_axis == :z
    @test pqs_source_descriptor.q == 5
    @test pqs_source_descriptor.L == 6
    @test pqs_source_descriptor.support_count == 362
    @test pqs_source_descriptor.mode_count == 114
    @test pqs_source_descriptor.retained_count == 114
    @test pqs_source_descriptor.cleanup_method == :projected_boundary_symmetric_lowdin
    @test pqs_source_descriptor.cleanup_matrix_size == (114, 114)
    @test pqs_source_descriptor.cleanup_rank_count == 114
    @test pqs_source_descriptor.cleanup_rank_drop_count == 0
    @test pqs_source_descriptor.selection_rule == :any_axis_mode_index_first_or_last
    @test all(
        mode -> any(axis -> mode[axis] == 1 || mode[axis] == (5, 5, 6)[axis], 1:3),
        pqs_source_descriptor.boundary_mode_indices,
    )
    pqs_route_fact_diagnostic =
        CCPM._pqs_pqs_product_route_descriptor_diagnostic(
            pqs_construction,
            _pqs_axis_metrics(bundles),
        )
    @test pqs_route_fact_diagnostic.status == :descriptor_unavailable
    @test pqs_route_fact_diagnostic.descriptor === nothing
    @test :second_pqs_raw_plan in pqs_route_fact_diagnostic.missing
    @test :middle_product_doside_unit in pqs_route_fact_diagnostic.missing
    @test pqs_route_fact_diagnostic.diagnostics.pqs_descriptor_count == 1
    @test pqs_route_fact_diagnostic.diagnostics.pqs_raw_plan_convertible_count == 1
    @test pqs_route_fact_diagnostic.diagnostics.product_doside_unit_count == 0
    @test pqs_route_fact_diagnostic.diagnostics.direct_or_support_body_piece_count == 4
    @test !pqs_route_fact_diagnostic.diagnostics.descriptor_emitted
    @test pqs_route_fact_diagnostic.diagnostics.raw_plan_convertibility_checked
    @test pqs_route_fact_diagnostic.diagnostics.raw_plan_conversion_failure_count == 0
    @test pqs_route_fact_diagnostic.diagnostics.current_route_contains_pqs_descriptors
    @test !pqs_route_fact_diagnostic.diagnostics.current_route_contains_explicit_product_doside_body_unit
    @test !pqs_route_fact_diagnostic.diagnostics.packet_adoption
    @test !pqs_route_fact_diagnostic.diagnostics.fixed_block_construction_changed
    @test !pqs_route_fact_diagnostic.diagnostics.qwhamiltonian_changed
    @test !pqs_route_fact_diagnostic.diagnostics.shell_projection_used
    @test !pqs_route_fact_diagnostic.diagnostics.lowdin_cleanup_used
    @test !pqs_route_fact_diagnostic.diagnostics.support_local_pqs_oracle_used
    @test pqs_route_fact_diagnostic.diagnostics.retained_weight_semantics ==
          :not_positive_quadrature_weights
    @test !pqs_route_fact_diagnostic.diagnostics.ida_weight_division_allowed
    @test !pqs_route_fact_diagnostic.diagnostics.direct_support_reinterpreted_as_product_doside
    @test any(
        mismatch ->
            hasproperty(mismatch, :reason) &&
            mismatch.reason == :direct_or_support_piece_not_product_doside &&
            hasproperty(mismatch, :primitive_family) &&
            mismatch.primitive_family == :contact_cap_owned_slab,
        pqs_route_fact_diagnostic.mismatches,
    )
    @test any(
        mismatch ->
            hasproperty(mismatch, :reason) &&
            mismatch.reason == :shared_pqs_descriptors_are_not_route_left_right_group,
        pqs_route_fact_diagnostic.mismatches,
    )
    pqs_retained_unit_audit =
        CCPM._pqs_route_retained_unit_fact_audit(pqs_construction)
    @test pqs_retained_unit_audit.object_kind ==
          :pqs_route_retained_unit_fact_audit
    @test pqs_retained_unit_audit.status == :audit_only
    @test length(pqs_retained_unit_audit.unit_facts) == 5
    @test pqs_retained_unit_audit.summary.classification_counts.product_box_constructible == 2
    @test pqs_retained_unit_audit.summary.classification_counts.needs_direct_support_retained_unit_kind == 2
    @test pqs_retained_unit_audit.summary.classification_counts.out_of_scope == 1
    @test pqs_retained_unit_audit.summary.product_box_constructible_slab_rule_count == 3
    pqs_retained_facts_by_role =
        Dict(fact.role => fact for fact in pqs_retained_unit_audit.unit_facts)
    contact_fact = pqs_retained_facts_by_role[:contact_cap]
    @test contact_fact.classification == :product_box_constructible
    @test contact_fact.primitive_family == :contact_cap_owned_slab
    @test !contact_fact.raw_product_box_operator_contract
    @test contact_fact.product_box_construction_rule_available
    @test !contact_fact.product_doside_unit
    @test contact_fact.safe_term_capability ==
          :product_box_rule_available_not_instantiated
    @test contact_fact.coefficient_scope == :support_local_direct_rows
    @test contact_fact.construction_rule.rule_kind ==
          :identity_selector_product_doside_slab
    @test contact_fact.construction_rule.fixed_axis == :z
    @test contact_fact.construction_rule.fixed_index == 8
    @test contact_fact.construction_rule.active_axes == (:x, :y)
    @test contact_fact.construction_rule.active_intervals == (2:6, 2:6)
    @test contact_fact.construction_rule.retained_count == 25
    outer_fact =
        pqs_retained_facts_by_role[:outer_mismatch_shared_molecular_shell]
    @test outer_fact.classification == :product_box_constructible
    @test outer_fact.primitive_family == :outer_mismatch_boundary_slab_set
    @test !outer_fact.raw_product_box_operator_contract
    @test outer_fact.product_box_construction_rule_available
    @test !outer_fact.product_doside_unit
    @test outer_fact.construction_rule.rule_kind ==
          :identity_selector_product_doside_boundary_slab_set
    @test outer_fact.construction_rule.boundary_slab_set
    @test outer_fact.construction_rule.slab_piece_count == 2
    @test all(
        piece_rule -> piece_rule.fixed_axis == :z,
        outer_fact.construction_rule.slab_piece_rules,
    )
    @test sum(
        piece_rule -> piece_rule.support_count,
        outer_fact.construction_rule.slab_piece_rules,
    ) == outer_fact.support_count
    outer_mismatch_product_units =
        CCPM._pqs_outer_mismatch_product_doside_units(pqs_construction)
    @test outer_mismatch_product_units.object_kind ==
          :pqs_outer_mismatch_product_doside_units_fixture
    @test outer_mismatch_product_units.status == :private_diagnostic_only
    @test outer_mismatch_product_units.fact.role ==
          :outer_mismatch_shared_molecular_shell
    @test outer_mismatch_product_units.fact.primitive_family ==
          :outer_mismatch_boundary_slab_set
    @test !outer_mismatch_product_units.fact.raw_product_box_operator_contract
    @test outer_mismatch_product_units.fact.product_box_construction_rule_available
    @test length(outer_mismatch_product_units.units) == 2
    @test map(unit -> unit.kind, outer_mismatch_product_units.units) ==
          (:product_doside, :product_doside)
    @test map(unit -> unit.role, outer_mismatch_product_units.units) ==
          (:outer_mismatch_z_low_slab, :outer_mismatch_z_high_slab)
    @test map(unit -> unit.column_range, outer_mismatch_product_units.units) ==
          (1:49, 50:98)
    @test map(unit -> length(unit.support_indices), outer_mismatch_product_units.units) ==
          (49, 49)
    outer_mismatch_piece_support_indices = vcat(
        outer_mismatch_product_units.units[1].support_indices,
        outer_mismatch_product_units.units[2].support_indices,
    )
    @test sort(outer_mismatch_piece_support_indices) ==
          sort(outer_mismatch_product_units.fact.support_indices)
    @test map(
        unit -> map(axis -> axis.kind, unit.axes),
        outer_mismatch_product_units.units,
    ) == ((:active, :active, :fixed), (:active, :active, :fixed))
    @test outer_mismatch_product_units.units[1].axes[1].interval == 1:7
    @test outer_mismatch_product_units.units[1].axes[2].interval == 1:7
    @test outer_mismatch_product_units.units[1].axes[3].fixed_index == 1
    @test outer_mismatch_product_units.units[2].axes[1].interval == 1:7
    @test outer_mismatch_product_units.units[2].axes[2].interval == 1:7
    @test outer_mismatch_product_units.units[2].axes[3].fixed_index == 15
    @test all(
        unit -> Matrix(unit.coefficient_matrix) == Matrix{Float64}(I, 49, 49),
        outer_mismatch_product_units.units,
    )
    @test all(
        unit -> unit.axis_function_indices ==
                GaussletBases._nested_product_axis_function_indices(3, 1, 7, 2, 7),
        outer_mismatch_product_units.units,
    )
    @test all(
        equivalence -> equivalence.support_indices_match,
        outer_mismatch_product_units.piece_equivalences,
    )
    @test all(
        equivalence -> equivalence.retained_count_match,
        outer_mismatch_product_units.piece_equivalences,
    )
    @test map(
        equivalence -> equivalence.column_range,
        outer_mismatch_product_units.piece_equivalences,
    ) == (1:49, 50:98)
    @test all(
        equivalence -> equivalence.coefficient_matrix_matches_direct_selector,
        outer_mismatch_product_units.piece_equivalences,
    )
    @test all(
        equivalence -> equivalence.max_parent_coefficient_error == 0.0,
        outer_mismatch_product_units.piece_equivalences,
    )
    @test outer_mismatch_product_units.aggregate_equivalence.support_indices_match
    @test outer_mismatch_product_units.aggregate_equivalence.audited_support_set_match
    @test outer_mismatch_product_units.aggregate_equivalence.retained_count_match
    @test outer_mismatch_product_units.aggregate_equivalence.column_range_partition
    @test outer_mismatch_product_units.aggregate_equivalence.coefficient_matrix_matches_direct_selector
    @test outer_mismatch_product_units.aggregate_equivalence.max_parent_coefficient_error == 0.0
    @test outer_mismatch_product_units.diagnostics.outer_mismatch_only
    @test outer_mismatch_product_units.diagnostics.boundary_slab_set
    @test outer_mismatch_product_units.diagnostics.product_doside_units_created
    @test outer_mismatch_product_units.diagnostics.unit_count == 2
    @test outer_mismatch_product_units.diagnostics.slab_piece_count == 2
    @test !outer_mismatch_product_units.diagnostics.route_descriptor_emitted
    @test !outer_mismatch_product_units.diagnostics.construction_mutated
    @test !outer_mismatch_product_units.diagnostics.sidecar_installation
    @test !outer_mismatch_product_units.diagnostics.packet_adoption
    @test !outer_mismatch_product_units.diagnostics.fixed_block_construction_changed
    @test !outer_mismatch_product_units.diagnostics.qwhamiltonian_changed
    @test !outer_mismatch_product_units.diagnostics.ida_weight_division_allowed
    @test outer_mismatch_product_units.diagnostics.retained_weight_semantics ==
          :not_positive_quadrature_weights
    @test !outer_mismatch_product_units.diagnostics.input_fact_raw_product_box_operator_contract
    @test outer_mismatch_product_units.diagnostics.created_units_raw_product_box_operator_contract
    @test outer_mismatch_product_units.diagnostics.descriptor_piece_order_defines_columns
    @test outer_mismatch_product_units.diagnostics.audited_support_checked_as_set
    @test outer_mismatch_product_units.diagnostics.product_box_construction_rule_available
    @test !outer_mismatch_product_units.diagnostics.local_ecp_gaussian_mwg_interaction_changed
    outer_mismatch_safe_term_metrics = _pqs_axis_metrics(bundles)
    outer_mismatch_safe_terms = (
        :overlap,
        :position_x,
        :position_y,
        :position_z,
        :x2_x,
        :x2_y,
        :x2_z,
        :kinetic,
    )
    outer_mismatch_safe_term_comparison =
        CCPM._pqs_outer_mismatch_safe_term_operator_comparison(
            pqs_construction,
            outer_mismatch_safe_term_metrics,
        )
    @test outer_mismatch_safe_term_comparison.object_kind ==
          :pqs_outer_mismatch_safe_term_operator_comparison
    @test outer_mismatch_safe_term_comparison.status == :private_diagnostic_only
    @test outer_mismatch_safe_term_comparison.terms == outer_mismatch_safe_terms
    @test length(outer_mismatch_safe_term_comparison.fixture.units) == 2
    @test outer_mismatch_safe_term_comparison.max_block_error <= 1.0e-12
    @test outer_mismatch_safe_term_comparison.diagnostics.source ==
          :pqs_outer_mismatch_safe_term_operator_comparison
    @test outer_mismatch_safe_term_comparison.diagnostics.outer_mismatch_only
    @test outer_mismatch_safe_term_comparison.diagnostics.boundary_slab_set
    @test outer_mismatch_safe_term_comparison.diagnostics.private_diagnostic_only
    @test outer_mismatch_safe_term_comparison.diagnostics.terms_checked ==
          outer_mismatch_safe_terms
    @test outer_mismatch_safe_term_comparison.diagnostics.supported_terms ==
          outer_mismatch_safe_terms
    @test :weights in outer_mismatch_safe_term_comparison.diagnostics.unsupported_terms
    @test outer_mismatch_safe_term_comparison.diagnostics.product_path ==
          :_product_doside_source_box_reference_block
    @test outer_mismatch_safe_term_comparison.diagnostics.direct_oracle_path ==
          :support_local_direct_selector_contract_pair_block
    @test outer_mismatch_safe_term_comparison.diagnostics.product_doside_units_created
    @test outer_mismatch_safe_term_comparison.diagnostics.complete_slab_set_block_assembled
    @test outer_mismatch_safe_term_comparison.diagnostics.cross_slab_blocks_included
    @test outer_mismatch_safe_term_comparison.diagnostics.direct_support_oracle_compared
    @test !outer_mismatch_safe_term_comparison.diagnostics.route_descriptor_emitted
    @test !outer_mismatch_safe_term_comparison.diagnostics.construction_mutated
    @test !outer_mismatch_safe_term_comparison.diagnostics.sidecar_installation
    @test !outer_mismatch_safe_term_comparison.diagnostics.packet_adoption
    @test !outer_mismatch_safe_term_comparison.diagnostics.fixed_block_construction_changed
    @test !outer_mismatch_safe_term_comparison.diagnostics.qwhamiltonian_changed
    @test !outer_mismatch_safe_term_comparison.diagnostics.ida_weight_division_allowed
    @test outer_mismatch_safe_term_comparison.diagnostics.retained_weight_semantics ==
          :not_positive_quadrature_weights
    @test !outer_mismatch_safe_term_comparison.diagnostics.local_ecp_gaussian_mwg_interaction_changed
    @test outer_mismatch_safe_term_comparison.diagnostics.operator_factor_source ==
          :explicit_metric_operator_data
    @test outer_mismatch_safe_term_comparison.diagnostics.operator_metric_sources ==
          (:nested_pgdg_axis, :nested_pgdg_axis, :nested_pgdg_axis)
    @test !outer_mismatch_safe_term_comparison.diagnostics.input_metric_operator_data_pgdg_checked
    @test !outer_mismatch_safe_term_comparison.diagnostics.pgdg_analytic_operator_provenance_claimed
    @test !outer_mismatch_safe_term_comparison.diagnostics.numerical_reference_fallback
    @test outer_mismatch_safe_term_comparison.diagnostics.product_source_box_reference_compared
    @test outer_mismatch_safe_term_comparison.diagnostics.direct_support_oracle_entries_built
    @test outer_mismatch_safe_term_comparison.diagnostics.retained_count == 98
    @test outer_mismatch_safe_term_comparison.diagnostics.support_count == 98
    @test outer_mismatch_safe_term_comparison.diagnostics.unit_count == 2
    @test outer_mismatch_safe_term_comparison.diagnostics.max_block_error <= 1.0e-12
    @test outer_mismatch_safe_term_comparison.diagnostics.output_finite
    for term in outer_mismatch_safe_terms
        product_block = outer_mismatch_safe_term_comparison.product_blocks[term]
        oracle_block = outer_mismatch_safe_term_comparison.direct_oracle_blocks[term]
        @test size(product_block) == (98, 98)
        @test size(oracle_block) == (98, 98)
        @test all(isfinite, product_block)
        @test all(isfinite, oracle_block)
        @test product_block ≈ oracle_block atol = 1.0e-12 rtol = 0.0
        @test outer_mismatch_safe_term_comparison.term_errors[term] <= 1.0e-12
        @test outer_mismatch_safe_term_comparison.diagnostics.pair_block_counts[term] == 4
        @test outer_mismatch_safe_term_comparison.diagnostics.cross_slab_pair_counts[term] == 2
        @test length(outer_mismatch_safe_term_comparison.product_references[term].pair_references) == 4
    end
    @test_throws ArgumentError CCPM._pqs_outer_mismatch_safe_term_operator_comparison(
        pqs_construction,
        outer_mismatch_safe_term_metrics;
        terms = (:weights,),
    )
    left_atom_fact = pqs_retained_facts_by_role[:left_atom_box]
    right_atom_fact = pqs_retained_facts_by_role[:right_atom_box]
    @test left_atom_fact.classification ==
          :needs_direct_support_retained_unit_kind
    @test right_atom_fact.classification ==
          :needs_direct_support_retained_unit_kind
    @test !left_atom_fact.product_doside_unit
    @test !right_atom_fact.product_doside_unit
    @test !left_atom_fact.raw_product_box_operator_contract
    @test !right_atom_fact.raw_product_box_operator_contract
    @test !left_atom_fact.product_box_construction_rule_available
    @test !right_atom_fact.product_box_construction_rule_available
    @test left_atom_fact.safe_term_capability == :support_local_reference_only
    @test right_atom_fact.safe_term_capability == :support_local_reference_only
    atom_box_support_dense_units =
        CCPM._pqs_atom_box_support_dense_units(pqs_construction)
    @test atom_box_support_dense_units.object_kind ==
          :pqs_atom_box_support_dense_units_fixture
    @test atom_box_support_dense_units.status == :private_diagnostic_only
    @test length(atom_box_support_dense_units.units) == 2
    @test map(unit -> unit.role, atom_box_support_dense_units.units) ==
          (:left_atom_box, :right_atom_box)
    @test map(unit -> unit.kind, atom_box_support_dense_units.units) ==
          (:support_dense, :support_dense)
    @test map(unit -> unit.column_range, atom_box_support_dense_units.units) ==
          (99:223, 224:348)
    @test map(unit -> length(unit.support_indices), atom_box_support_dense_units.units) ==
          (125, 125)
    @test map(unit -> size(unit.coefficient_matrix), atom_box_support_dense_units.units) ==
          ((125, 125), (125, 125))
    @test all(
        unit -> map(axis -> axis.kind, unit.axes) == (:fixed, :fixed, :fixed),
        atom_box_support_dense_units.units,
    )
    @test all(
        unit -> all(index -> index == (1, 1, 1), unit.axis_function_indices),
        atom_box_support_dense_units.units,
    )
    @test all(
        equivalence -> equivalence.support_indices_match,
        atom_box_support_dense_units.equivalences,
    )
    @test all(
        equivalence -> equivalence.audited_support_set_match,
        atom_box_support_dense_units.equivalences,
    )
    @test all(
        equivalence -> equivalence.retained_count_match,
        atom_box_support_dense_units.equivalences,
    )
    @test all(
        equivalence -> equivalence.coefficient_matrix_matches_direct_support,
        atom_box_support_dense_units.equivalences,
    )
    @test all(
        equivalence -> equivalence.max_parent_coefficient_error == 0.0,
        atom_box_support_dense_units.equivalences,
    )
    @test all(
        equivalence -> equivalence.local_identity_error <= 1.0e-12,
        atom_box_support_dense_units.equivalences,
    )
    @test atom_box_support_dense_units.diagnostics.atom_box_only
    @test atom_box_support_dense_units.diagnostics.support_dense_direct_support_units_created
    @test !atom_box_support_dense_units.diagnostics.product_doside_units_created
    @test !atom_box_support_dense_units.diagnostics.raw_product_box_operator_contract
    @test atom_box_support_dense_units.diagnostics.support_local_reference_only
    @test !atom_box_support_dense_units.diagnostics.product_box_construction_rule_available
    @test !atom_box_support_dense_units.diagnostics.route_descriptor_emitted
    @test !atom_box_support_dense_units.diagnostics.construction_mutated
    @test !atom_box_support_dense_units.diagnostics.sidecar_installation
    @test !atom_box_support_dense_units.diagnostics.packet_adoption
    @test !atom_box_support_dense_units.diagnostics.fixed_block_construction_changed
    @test !atom_box_support_dense_units.diagnostics.qwhamiltonian_changed
    @test !atom_box_support_dense_units.diagnostics.ida_weight_division_allowed
    @test atom_box_support_dense_units.diagnostics.retained_weight_semantics ==
          :not_positive_quadrature_weights
    @test !atom_box_support_dense_units.diagnostics.local_ecp_gaussian_mwg_interaction_changed
    @test atom_box_support_dense_units.diagnostics.max_parent_coefficient_error == 0.0
    @test !atom_box_support_dense_units.diagnostics.local_identity_is_product_box_claim
    @test !atom_box_support_dense_units.diagnostics.safe_term_operator_comparison_added
    atom_box_safe_term_metrics = _pqs_axis_metrics(bundles)
    atom_box_safe_terms = (
        :overlap,
        :position_x,
        :position_y,
        :position_z,
        :x2_x,
        :x2_y,
        :x2_z,
        :kinetic,
    )
    atom_box_safe_term_comparison =
        CCPM._pqs_atom_box_safe_term_operator_comparison(
            pqs_construction,
            atom_box_safe_term_metrics,
        )
    @test atom_box_safe_term_comparison.object_kind ==
          :pqs_atom_box_safe_term_operator_comparison
    @test atom_box_safe_term_comparison.status == :private_diagnostic_only
    @test atom_box_safe_term_comparison.terms == atom_box_safe_terms
    @test length(atom_box_safe_term_comparison.fixture.units) == 2
    @test atom_box_safe_term_comparison.max_block_error <= 1.0e-12
    @test atom_box_safe_term_comparison.diagnostics.source ==
          :pqs_atom_box_safe_term_operator_comparison
    @test atom_box_safe_term_comparison.diagnostics.atom_box_only
    @test atom_box_safe_term_comparison.diagnostics.support_dense_direct_support_units_created
    @test atom_box_safe_term_comparison.diagnostics.support_local_fallback_operator_comparison
    @test !atom_box_safe_term_comparison.diagnostics.product_doside_units_created
    @test !atom_box_safe_term_comparison.diagnostics.raw_product_box_operator_contract
    @test !atom_box_safe_term_comparison.diagnostics.product_box_construction_rule_available
    @test atom_box_safe_term_comparison.diagnostics.complete_atom_box_block_assembled
    @test atom_box_safe_term_comparison.diagnostics.cross_atom_blocks_included
    @test atom_box_safe_term_comparison.diagnostics.direct_support_oracle_compared
    @test !atom_box_safe_term_comparison.diagnostics.route_descriptor_emitted
    @test !atom_box_safe_term_comparison.diagnostics.construction_mutated
    @test !atom_box_safe_term_comparison.diagnostics.sidecar_installation
    @test !atom_box_safe_term_comparison.diagnostics.packet_adoption
    @test !atom_box_safe_term_comparison.diagnostics.fixed_block_construction_changed
    @test !atom_box_safe_term_comparison.diagnostics.qwhamiltonian_changed
    @test !atom_box_safe_term_comparison.diagnostics.ida_weight_division_allowed
    @test atom_box_safe_term_comparison.diagnostics.retained_weight_semantics ==
          :not_positive_quadrature_weights
    @test !atom_box_safe_term_comparison.diagnostics.local_ecp_gaussian_mwg_interaction_changed
    @test atom_box_safe_term_comparison.diagnostics.operator_factor_source ==
          :explicit_metric_operator_data
    @test atom_box_safe_term_comparison.diagnostics.operator_metric_sources ==
          (:nested_pgdg_axis, :nested_pgdg_axis, :nested_pgdg_axis)
    @test !atom_box_safe_term_comparison.diagnostics.input_metric_operator_data_pgdg_checked
    @test !atom_box_safe_term_comparison.diagnostics.pgdg_analytic_operator_provenance_claimed
    @test !atom_box_safe_term_comparison.diagnostics.numerical_reference_fallback
    @test atom_box_safe_term_comparison.diagnostics.retained_count == 250
    @test atom_box_safe_term_comparison.diagnostics.support_count == 250
    @test atom_box_safe_term_comparison.diagnostics.unit_count == 2
    @test atom_box_safe_term_comparison.diagnostics.max_block_error <= 1.0e-12
    @test atom_box_safe_term_comparison.diagnostics.output_finite
    for term in atom_box_safe_terms
        support_dense_block = atom_box_safe_term_comparison.support_dense_blocks[term]
        oracle_block = atom_box_safe_term_comparison.direct_oracle_blocks[term]
        @test size(support_dense_block) == (250, 250)
        @test size(oracle_block) == (250, 250)
        @test all(isfinite, support_dense_block)
        @test all(isfinite, oracle_block)
        @test support_dense_block ≈ oracle_block atol = 1.0e-12 rtol = 0.0
        @test atom_box_safe_term_comparison.term_errors[term] <= 1.0e-12
        @test atom_box_safe_term_comparison.diagnostics.pair_block_counts[term] == 4
        @test atom_box_safe_term_comparison.diagnostics.cross_atom_pair_counts[term] == 2
    end
    @test_throws ArgumentError CCPM._pqs_atom_box_safe_term_operator_comparison(
        pqs_construction,
        atom_box_safe_term_metrics;
        terms = (:weights,),
    )
    shared_pqs_fact =
        pqs_retained_facts_by_role[:regular_shared_molecular_shell]
    @test shared_pqs_fact.classification == :out_of_scope
    @test shared_pqs_fact.primitive_family == :projected_q_shell
    @test !shared_pqs_fact.product_doside_unit
    @test !shared_pqs_fact.raw_product_box_operator_contract
    @test !shared_pqs_fact.product_box_construction_rule_available
    @test shared_pqs_fact.safe_term_capability == :not_body_retained_unit
    @test shared_pqs_fact.notes.current_single_pqs_descriptor
    @test pqs_retained_unit_audit.diagnostics.private_diagnostic_only
    @test !pqs_retained_unit_audit.diagnostics.descriptor_emitted
    @test !pqs_retained_unit_audit.diagnostics.packet_adoption
    @test !pqs_retained_unit_audit.diagnostics.fixed_block_construction_changed
    @test !pqs_retained_unit_audit.diagnostics.qwhamiltonian_changed
    @test !pqs_retained_unit_audit.diagnostics.sidecar_mutation
    @test !pqs_retained_unit_audit.diagnostics.sidecar_installation
    @test !pqs_retained_unit_audit.diagnostics.direct_support_reinterpreted_as_product_doside
    @test pqs_retained_unit_audit.diagnostics.retained_weight_semantics ==
          :not_positive_quadrature_weights
    @test !pqs_retained_unit_audit.diagnostics.ida_weight_division_allowed
    @test !pqs_retained_unit_audit.diagnostics.local_ecp_gaussian_mwg_interaction_changed
    contact_cap_product_unit =
        CCPM._pqs_contact_cap_product_doside_unit(pqs_construction)
    @test contact_cap_product_unit.object_kind ==
          :pqs_contact_cap_product_doside_unit_fixture
    @test contact_cap_product_unit.status == :private_diagnostic_only
    @test contact_cap_product_unit.fact.role == :contact_cap
    @test contact_cap_product_unit.fact.classification ==
          :product_box_constructible
    @test !contact_cap_product_unit.fact.raw_product_box_operator_contract
    @test contact_cap_product_unit.fact.product_box_construction_rule_available
    @test contact_cap_product_unit.unit.role == :contact_cap_slab
    @test contact_cap_product_unit.unit.kind == :product_doside
    @test contact_cap_product_unit.unit.column_range == 349:373
    @test contact_cap_product_unit.unit.support_indices ==
          contact_cap_product_unit.fact.support_indices
    @test contact_cap_product_unit.unit.support_states == [
        GaussletBases._cartesian_unflat_index(index, (7, 7, 15))
        for index in contact_cap_product_unit.unit.support_indices
    ]
    @test Matrix(contact_cap_product_unit.unit.coefficient_matrix) ==
          Matrix{Float64}(I, 25, 25)
    @test map(axis -> axis.kind, contact_cap_product_unit.unit.axes) ==
          (:active, :active, :fixed)
    @test contact_cap_product_unit.unit.axes[1].interval == 2:6
    @test contact_cap_product_unit.unit.axes[2].interval == 2:6
    @test contact_cap_product_unit.unit.axes[3].fixed_index == 8
    @test Matrix(contact_cap_product_unit.unit.axes[1].coefficient_matrix) ==
          Matrix{Float64}(I, 5, 5)
    @test Matrix(contact_cap_product_unit.unit.axes[2].coefficient_matrix) ==
          Matrix{Float64}(I, 5, 5)
    @test contact_cap_product_unit.unit.axis_function_indices ==
          GaussletBases._nested_product_axis_function_indices(3, 1, 5, 2, 5)
    @test contact_cap_product_unit.equivalence.support_indices_match
    @test contact_cap_product_unit.equivalence.support_states_match
    @test contact_cap_product_unit.equivalence.retained_count_match
    @test contact_cap_product_unit.equivalence.column_range_match
    @test contact_cap_product_unit.equivalence.coefficient_matrix_matches_direct_selector
    @test contact_cap_product_unit.equivalence.max_parent_coefficient_error == 0.0
    @test contact_cap_product_unit.diagnostics.contact_cap_only
    @test contact_cap_product_unit.diagnostics.product_doside_unit_created
    @test !contact_cap_product_unit.diagnostics.route_descriptor_emitted
    @test !contact_cap_product_unit.diagnostics.construction_mutated
    @test !contact_cap_product_unit.diagnostics.sidecar_installation
    @test !contact_cap_product_unit.diagnostics.packet_adoption
    @test !contact_cap_product_unit.diagnostics.fixed_block_construction_changed
    @test !contact_cap_product_unit.diagnostics.qwhamiltonian_changed
    @test !contact_cap_product_unit.diagnostics.ida_weight_division_allowed
    @test contact_cap_product_unit.diagnostics.retained_weight_semantics ==
          :not_positive_quadrature_weights
    @test contact_cap_product_unit.diagnostics.product_box_construction_rule_available
    @test !contact_cap_product_unit.diagnostics.input_fact_raw_product_box_operator_contract
    @test contact_cap_product_unit.diagnostics.created_unit_raw_product_box_operator_contract
    contact_safe_term_metrics = _pqs_axis_metrics(bundles)
    contact_safe_term_comparison =
        CCPM._pqs_contact_cap_safe_term_operator_comparison(
            pqs_construction,
            contact_safe_term_metrics,
        )
    contact_safe_terms = (
        :overlap,
        :position_x,
        :position_y,
        :position_z,
        :x2_x,
        :x2_y,
        :x2_z,
        :kinetic,
    )
    @test contact_safe_term_comparison.object_kind ==
          :pqs_contact_cap_safe_term_operator_comparison
    @test contact_safe_term_comparison.status == :private_diagnostic_only
    @test contact_safe_term_comparison.terms == contact_safe_terms
    @test contact_safe_term_comparison.fixture.unit.kind == :product_doside
    @test contact_safe_term_comparison.fixture.unit.column_range == 349:373
    @test contact_safe_term_comparison.max_block_error <= 1.0e-12
    @test contact_safe_term_comparison.diagnostics.source ==
          :pqs_contact_cap_safe_term_operator_comparison
    @test contact_safe_term_comparison.diagnostics.contact_cap_only
    @test contact_safe_term_comparison.diagnostics.private_diagnostic_only
    @test contact_safe_term_comparison.diagnostics.terms_checked ==
          contact_safe_terms
    @test contact_safe_term_comparison.diagnostics.supported_terms ==
          contact_safe_terms
    @test :weights in contact_safe_term_comparison.diagnostics.unsupported_terms
    @test contact_safe_term_comparison.diagnostics.product_path ==
          :_product_doside_source_box_reference_block
    @test contact_safe_term_comparison.diagnostics.direct_oracle_path ==
          :support_local_direct_selector_contract_pair_block
    @test contact_safe_term_comparison.diagnostics.current_direct_support_selector_compared
    @test contact_safe_term_comparison.diagnostics.product_doside_unit_created
    @test !contact_safe_term_comparison.diagnostics.input_fact_raw_product_box_operator_contract
    @test contact_safe_term_comparison.diagnostics.created_unit_raw_product_box_operator_contract
    @test contact_safe_term_comparison.diagnostics.product_box_construction_rule_available
    @test !contact_safe_term_comparison.diagnostics.route_descriptor_emitted
    @test !contact_safe_term_comparison.diagnostics.construction_mutated
    @test !contact_safe_term_comparison.diagnostics.sidecar_installation
    @test !contact_safe_term_comparison.diagnostics.packet_adoption
    @test !contact_safe_term_comparison.diagnostics.fixed_block_construction_changed
    @test !contact_safe_term_comparison.diagnostics.qwhamiltonian_changed
    @test !contact_safe_term_comparison.diagnostics.ida_weight_division_allowed
    @test contact_safe_term_comparison.diagnostics.retained_weight_semantics ==
          :not_positive_quadrature_weights
    @test !contact_safe_term_comparison.diagnostics.local_ecp_gaussian_mwg_interaction_changed
    @test contact_safe_term_comparison.diagnostics.operator_factor_source ==
          :explicit_metric_operator_data
    @test contact_safe_term_comparison.diagnostics.operator_metric_sources ==
          (:nested_pgdg_axis, :nested_pgdg_axis, :nested_pgdg_axis)
    @test !contact_safe_term_comparison.diagnostics.input_metric_operator_data_pgdg_checked
    @test !contact_safe_term_comparison.diagnostics.pgdg_analytic_operator_provenance_claimed
    @test !contact_safe_term_comparison.diagnostics.numerical_reference_fallback
    @test contact_safe_term_comparison.diagnostics.product_source_box_reference_compared
    @test contact_safe_term_comparison.diagnostics.direct_support_oracle_entries_built
    @test contact_safe_term_comparison.diagnostics.retained_count == 25
    @test contact_safe_term_comparison.diagnostics.support_count == 25
    @test contact_safe_term_comparison.diagnostics.column_range == 349:373
    @test contact_safe_term_comparison.diagnostics.max_block_error <= 1.0e-12
    @test contact_safe_term_comparison.diagnostics.output_finite
    for term in contact_safe_terms
        product_block = contact_safe_term_comparison.product_blocks[term]
        oracle_block = contact_safe_term_comparison.direct_oracle_blocks[term]
        @test size(product_block) == (25, 25)
        @test size(oracle_block) == (25, 25)
        @test all(isfinite, product_block)
        @test all(isfinite, oracle_block)
        @test product_block ≈ oracle_block atol = 1.0e-12 rtol = 0.0
        @test contact_safe_term_comparison.term_errors[term] <= 1.0e-12
        @test contact_safe_term_comparison.product_references[term].diagnostics.authoritative_block_compared
    end
    @test_throws ArgumentError CCPM._pqs_contact_cap_safe_term_operator_comparison(
        pqs_construction,
        contact_safe_term_metrics;
        terms = (:weights,),
    )
    current_route_inventory =
        CCPM._pqs_current_route_retained_unit_inventory(pqs_construction)
    @test current_route_inventory.object_kind ==
          :pqs_current_route_retained_unit_inventory_fixture
    @test current_route_inventory.status == :private_diagnostic_only
    @test length(current_route_inventory.units) == 6
    @test current_route_inventory.coverage.ordered_roles == (
        :outer_mismatch_z_low_slab,
        :outer_mismatch_z_high_slab,
        :left_atom_box,
        :right_atom_box,
        :contact_cap_slab,
        :regular_shared_molecular_shell,
    )
    @test current_route_inventory.coverage.ordered_column_ranges == (
        1:49,
        50:98,
        99:223,
        224:348,
        349:373,
        374:487,
    )
    @test current_route_inventory.coverage.first_column == 1
    @test current_route_inventory.coverage.last_column == 487
    @test current_route_inventory.coverage.represented_count == 487
    @test current_route_inventory.coverage.contiguous
    @test current_route_inventory.coverage.non_overlapping
    @test current_route_inventory.coverage.covers_every_column_once
    @test map(unit -> unit.category, current_route_inventory.units) == (
        :product_doside,
        :product_doside,
        :support_dense,
        :support_dense,
        :product_doside,
        :shell_realized_pqs_fixture,
    )
    @test map(unit -> unit.retained_count, current_route_inventory.units) ==
          (49, 49, 125, 125, 25, 114)
    @test current_route_inventory.by_role[:regular_shared_molecular_shell].kind ==
          :projected_q_shell
    @test current_route_inventory.by_role[:regular_shared_molecular_shell].active_representation_stage ==
          :shell_realized_pqs_fixture
    @test !current_route_inventory.by_role[:regular_shared_molecular_shell].raw_product_box_operator_contract
    @test current_route_inventory.by_role[:regular_shared_molecular_shell].safe_term_capability ==
          :support_local_oracle_for_shell_realization
    @test current_route_inventory.by_role[:regular_shared_molecular_shell].support_count == 362
    @test current_route_inventory.by_role[:regular_shared_molecular_shell].retained_count == 114
    @test current_route_inventory.by_role[:regular_shared_molecular_shell].raw_box_auxiliary_metadata.available
    @test current_route_inventory.by_role[:regular_shared_molecular_shell].raw_box_auxiliary_metadata.reference_only
    @test !current_route_inventory.by_role[:regular_shared_molecular_shell].raw_box_auxiliary_metadata.active_current_route_contract
    shared_pqs_unit = current_route_inventory.by_role[:regular_shared_molecular_shell]
    shared_transform_fact = shared_pqs_unit.shell_realization_transform_fact
    @test shared_transform_fact.object_kind ==
          :pqs_current_route_shell_realization_transform_fact
    @test shared_transform_fact.status == :metadata_precursor
    @test shared_transform_fact.representation_stage == :shell_realized_pqs_fixture
    @test shared_transform_fact.source_box.source_mode_dims ==
          shared_pqs_unit.raw_box_auxiliary_metadata.source_mode_dims
    @test shared_transform_fact.source_box.source_mode_count ==
          prod(shared_transform_fact.source_box.source_mode_dims)
    @test shared_transform_fact.boundary_selection.mode_count == shared_pqs_unit.retained_count
    @test shared_transform_fact.shell_projection.support_count == shared_pqs_unit.support_count
    @test shared_transform_fact.shell_projection.matrix_shape ==
          (shared_pqs_unit.support_count, shared_pqs_unit.retained_count)
    @test shared_transform_fact.lowdin_cleanup.transform_shape ==
          (shared_pqs_unit.retained_count, shared_pqs_unit.retained_count)
    @test shared_transform_fact.retained_columns.support_local_coefficient_shape ==
          size(shared_pqs_unit.support_local_coefficient_matrix)
    @test shared_transform_fact.retained_columns.coefficient_matches_descriptor_realization
    @test shared_transform_fact.retained_columns.max_support_local_coefficient_error <=
          1.0e-12
    @test shared_transform_fact.shell_realization.shell_projection_used
    @test shared_transform_fact.shell_realization.lowdin_cleanup_used
    @test !shared_transform_fact.compact_source_space_transform.available
    @test !shared_transform_fact.source_box_operator_application_ready
    @test shared_transform_fact.diagnostics.support_local_oracle_used
    @test shared_transform_fact.diagnostics.shell_row_oracle_only
    @test shared_transform_fact.diagnostics.metadata_precursor
    shared_transform_fact_checked =
        CCPM._pqs_current_route_shell_realization_transform_fact(
            shared_pqs_unit;
            metrics = contact_safe_term_metrics,
        )
    @test shared_transform_fact_checked.shell_realization.isometry_checked
    @test shared_transform_fact_checked.shell_realization.isometry_error <= 1.0e-8
    @test shared_transform_fact_checked.shell_realization.isometric
    @test current_route_inventory.by_role[:regular_shared_molecular_shell].equivalence.coefficient_matrix_matches_active_shell
    @test current_route_inventory.by_role[:regular_shared_molecular_shell].equivalence.max_parent_coefficient_error == 0.0
    @test current_route_inventory.by_role[:left_atom_box].category == :support_dense
    @test current_route_inventory.by_role[:left_atom_box].safe_term_capability ==
          :support_local_fallback_safe_terms
    @test current_route_inventory.by_role[:contact_cap_slab].category == :product_doside
    @test current_route_inventory.by_role[:outer_mismatch_z_low_slab].category ==
          :product_doside
    @test map(policy -> policy.pair_type, current_route_inventory.pair_policies) == (
        :product_product,
        :support_support,
        :support_product,
        :shell_realized_pqs_product,
        :shell_realized_pqs_support,
        :shell_realized_pqs_pqs,
        :raw_box_pqs_helpers,
    )
    @test current_route_inventory.pair_policies[1].policy ==
          :product_doside_source_box_path
    @test current_route_inventory.pair_policies[4].policy ==
          :support_local_oracle_for_shell_realization
    @test !current_route_inventory.pair_policies[4].active_current_route
    @test !current_route_inventory.pair_policies[4].active_algorithmic_policy
    @test !current_route_inventory.pair_policies[4].source_box_algorithm_available
    @test current_route_inventory.pair_policies[4].support_local_oracle_used
    @test current_route_inventory.pair_policies[4].shell_row_oracle_only
    @test !current_route_inventory.pair_policies[end].active_current_route
    @test current_route_inventory.diagnostics.private_diagnostic_only
    @test current_route_inventory.diagnostics.current_route_inventory
    @test !current_route_inventory.diagnostics.route_descriptor_emitted
    @test !current_route_inventory.diagnostics.construction_mutated
    @test !current_route_inventory.diagnostics.sidecar_installation
    @test !current_route_inventory.diagnostics.packet_adoption
    @test !current_route_inventory.diagnostics.fixed_block_construction_changed
    @test !current_route_inventory.diagnostics.qwhamiltonian_changed
    @test !current_route_inventory.diagnostics.ida_weight_division_allowed
    @test current_route_inventory.diagnostics.retained_weight_semantics ==
          :not_positive_quadrature_weights
    @test !current_route_inventory.diagnostics.local_ecp_gaussian_mwg_interaction_changed
    @test current_route_inventory.diagnostics.shared_pqs_active_representation ==
          :shell_realized_pqs_fixture
    @test !current_route_inventory.diagnostics.shared_pqs_raw_box_operator_contract
    @test current_route_inventory.diagnostics.shared_pqs_shell_realization_transform_fact_count == 1
    @test current_route_inventory.diagnostics.shared_pqs_source_box_operator_application_ready_count == 0
    @test current_route_inventory.diagnostics.raw_box_pqs_auxiliary_reference_available
    @test !current_route_inventory.diagnostics.whole_route_safe_term_matrix_consumer
    @test current_route_inventory.diagnostics.fixed_dimension == 487
    @test current_route_inventory.diagnostics.coverage_complete
    current_route_pair_inventory =
        CCPM._pqs_current_route_retained_pair_inventory(current_route_inventory)
    @test current_route_pair_inventory.object_kind ==
          :pqs_current_route_retained_pair_inventory_fixture
    @test current_route_pair_inventory.status == :private_diagnostic_only
    @test current_route_pair_inventory.unit_inventory === current_route_inventory
    @test length(current_route_pair_inventory.pairs) == 21
    expected_pair_roles = Tuple(
        (current_route_inventory.units[left].role, current_route_inventory.units[right].role)
        for left in 1:length(current_route_inventory.units)
        for right in left:length(current_route_inventory.units)
    )
    expected_pair_shapes = Tuple(
        (
            current_route_inventory.units[left].retained_count,
            current_route_inventory.units[right].retained_count,
        ) for left in 1:length(current_route_inventory.units)
        for right in left:length(current_route_inventory.units)
    )
    @test map(
        pair -> (pair.left_role, pair.right_role),
        current_route_pair_inventory.pairs,
    ) == expected_pair_roles
    @test map(pair -> pair.pair_shape, current_route_pair_inventory.pairs) ==
          expected_pair_shapes
    @test current_route_pair_inventory.counts.pair_count == 21
    @test current_route_pair_inventory.counts.product_product == 6
    @test current_route_pair_inventory.counts.support_support == 3
    @test current_route_pair_inventory.counts.support_product == 6
    @test current_route_pair_inventory.counts.shell_realized_pqs_product == 3
    @test current_route_pair_inventory.counts.shell_realized_pqs_support == 2
    @test current_route_pair_inventory.counts.shell_realized_pqs_pqs == 1
    @test current_route_pair_inventory.counts.raw_box_pqs_active == 0
    @test current_route_pair_inventory.counts.active_algorithmic_policy == 15
    @test current_route_pair_inventory.counts.source_box_algorithm_available == 6
    @test current_route_pair_inventory.counts.product_doside_source_box_path == 6
    @test current_route_pair_inventory.counts.support_local_fallback == 9
    @test current_route_pair_inventory.counts.support_local_oracle_for_shell_realization == 6
    @test current_route_pair_inventory.counts.support_local_oracle_pair_count == 6
    @test current_route_pair_inventory.counts.shell_row_oracle_pair_count == 6
    @test current_route_pair_inventory.pairs[1].pair_group == :product_product
    @test current_route_pair_inventory.pairs[1].policy ==
          :product_doside_source_box_path
    @test current_route_pair_inventory.pairs[3].pair_group == :support_product
    @test current_route_pair_inventory.pairs[3].policy == :support_local_fallback
    @test current_route_pair_inventory.pairs[6].pair_group ==
          :shell_realized_pqs_product
    @test current_route_pair_inventory.pairs[6].policy ==
          :support_local_oracle_for_shell_realization
    @test current_route_pair_inventory.pairs[6].shell_row_oracle_only
    @test current_route_pair_inventory.pairs[6].support_local_oracle_used
    @test !current_route_pair_inventory.pairs[6].active_algorithmic_policy
    @test !current_route_pair_inventory.pairs[6].source_box_algorithm_available
    @test current_route_pair_inventory.pairs[15].pair_group ==
          :shell_realized_pqs_support
    @test current_route_pair_inventory.pairs[15].policy ==
          :support_local_oracle_for_shell_realization
    @test current_route_pair_inventory.pairs[end].pair_group ==
          :shell_realized_pqs_pqs
    @test current_route_pair_inventory.pairs[end].policy ==
          :support_local_oracle_for_shell_realization
    @test current_route_pair_inventory.pairs[end].shell_row_oracle_only
    @test current_route_pair_inventory.pairs[end].support_local_oracle_used
    @test count(pair -> pair.active_current_route, current_route_pair_inventory.pairs) == 15
    @test count(pair -> pair.active_algorithmic_policy, current_route_pair_inventory.pairs) == 15
    @test count(pair -> pair.shell_row_oracle_only, current_route_pair_inventory.pairs) == 6
    @test all(
        pair -> !pair.raw_box_pqs_active_pair_policy,
        current_route_pair_inventory.pairs,
    )
    @test current_route_pair_inventory.diagnostics.private_diagnostic_only
    @test current_route_pair_inventory.diagnostics.current_route_pair_inventory
    @test current_route_pair_inventory.diagnostics.unit_inventory_complete
    @test current_route_pair_inventory.diagnostics.upper_triangular_pairs
    @test current_route_pair_inventory.diagnostics.pair_count == 21
    @test current_route_pair_inventory.diagnostics.raw_box_pqs_active_pair_policy_count == 0
    @test current_route_pair_inventory.diagnostics.active_algorithmic_policy_pair_count == 15
    @test current_route_pair_inventory.diagnostics.source_box_algorithm_available_pair_count == 6
    @test current_route_pair_inventory.diagnostics.support_local_oracle_for_shell_realization_pair_count == 6
    @test current_route_pair_inventory.diagnostics.shell_row_oracle_pair_count == 6
    @test current_route_pair_inventory.diagnostics.shell_realized_pqs_pairs_are_oracle_only
    @test !current_route_pair_inventory.diagnostics.route_descriptor_emitted
    @test !current_route_pair_inventory.diagnostics.construction_mutated
    @test !current_route_pair_inventory.diagnostics.sidecar_installation
    @test !current_route_pair_inventory.diagnostics.packet_adoption
    @test !current_route_pair_inventory.diagnostics.fixed_block_construction_changed
    @test !current_route_pair_inventory.diagnostics.qwhamiltonian_changed
    @test !current_route_pair_inventory.diagnostics.ida_weight_division_allowed
    @test current_route_pair_inventory.diagnostics.retained_weight_semantics ==
          :not_positive_quadrature_weights
    @test !current_route_pair_inventory.diagnostics.local_ecp_gaussian_mwg_interaction_changed
    @test !current_route_pair_inventory.diagnostics.whole_route_safe_term_matrix_consumer
    be2_inventory_timed = @timed begin
        be2_pqs_construction =
            GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction(
                be2_basis,
                be2_bundles,
                be2_policy;
                nside = 5,
                term_coefficients = Float64.(expansion.coefficients),
                packet_kernel = :support_reference,
                shared_shell_realization = :projected_q_shell,
            )
        be2_current_route_inventory =
            CCPM._pqs_current_route_retained_unit_inventory(be2_pqs_construction)
        be2_current_route_pair_inventory =
            CCPM._pqs_current_route_retained_pair_inventory(
                be2_current_route_inventory,
            )
        (
            construction = be2_pqs_construction,
            inventory = be2_current_route_inventory,
            pair_inventory = be2_current_route_pair_inventory,
        )
    end
    be2_inventory_payload = be2_inventory_timed.value
    be2_current_route_inventory = be2_inventory_payload.inventory
    be2_current_route_pair_inventory = be2_inventory_payload.pair_inventory
    be2_shared_pqs_units = Tuple(
        unit for unit in be2_current_route_inventory.units
        if unit.category == :shell_realized_pqs_fixture
    )
    @test be2_inventory_timed.time >= 0.0
    @test be2_inventory_timed.bytes >= 0
    @test be2_current_route_inventory.object_kind ==
          :pqs_current_route_retained_unit_inventory_fixture
    @test be2_current_route_inventory.status == :private_diagnostic_only
    @test length(be2_current_route_inventory.units) == 8
    @test be2_current_route_inventory.diagnostics.unit_count == 8
    @test be2_current_route_inventory.diagnostics.shared_pqs_unit_count == 3
    @test be2_current_route_inventory.diagnostics.shared_pqs_roles == (
        :regular_shared_molecular_shell_1,
        :regular_shared_molecular_shell_2,
        :regular_shared_molecular_shell_3,
    )
    @test be2_current_route_inventory.diagnostics.shared_pqs_original_roles ==
          (:regular_shared_molecular_shell, :regular_shared_molecular_shell, :regular_shared_molecular_shell)
    @test map(unit -> unit.role, be2_shared_pqs_units) ==
          be2_current_route_inventory.diagnostics.shared_pqs_roles
    @test map(unit -> unit.original_role, be2_shared_pqs_units) ==
          be2_current_route_inventory.diagnostics.shared_pqs_original_roles
    @test map(unit -> unit.retained_count, be2_shared_pqs_units) == (98, 98, 114)
    @test map(unit -> unit.support_count, be2_shared_pqs_units) ==
          (1738, 1346, 1002)
    @test map(unit -> unit.column_range, be2_shared_pqs_units) ==
          (1174:1271, 1272:1369, 1370:1483)
    @test map(unit -> unit.original_column_range, be2_shared_pqs_units) ==
          (1174:1271, 1272:1369, 1370:1483)
    @test be2_current_route_inventory.diagnostics.shared_pqs_column_ranges ==
          (1174:1271, 1272:1369, 1370:1483)
    @test be2_current_route_inventory.diagnostics.shared_pqs_original_column_ranges ==
          (1174:1271, 1272:1369, 1370:1483)
    @test be2_current_route_inventory.diagnostics.shared_pqs_shell_realization_transform_fact_count == 3
    @test be2_current_route_inventory.diagnostics.shared_pqs_source_box_operator_application_ready_count == 0
    @test be2_current_route_inventory.coverage.first_column == 1
    @test be2_current_route_inventory.coverage.last_column == 1483
    @test be2_current_route_inventory.coverage.represented_count == 1483
    @test be2_current_route_inventory.coverage.covers_every_column_once
    @test be2_current_route_inventory.diagnostics.fixed_dimension == 1483
    @test be2_current_route_inventory.diagnostics.coverage_complete
    @test !be2_current_route_inventory.diagnostics.q4_single_shared_role_order_preserved
    @test all(
        unit -> unit.active_representation_stage == :shell_realized_pqs_fixture,
        be2_shared_pqs_units,
    )
    @test all(
        unit -> unit.safe_term_capability == :support_local_oracle_for_shell_realization,
        be2_shared_pqs_units,
    )
    @test all(unit -> !unit.raw_product_box_operator_contract, be2_shared_pqs_units)
    @test all(
        unit -> unit.raw_box_auxiliary_metadata.available,
        be2_shared_pqs_units,
    )
    @test all(
        unit -> unit.raw_box_auxiliary_metadata.reference_only,
        be2_shared_pqs_units,
    )
    @test all(
        unit -> !unit.raw_box_auxiliary_metadata.active_current_route_contract,
        be2_shared_pqs_units,
    )
    @test all(
        unit -> unit.diagnostics.representation_stage == :shell_realized_pqs_fixture,
        be2_shared_pqs_units,
    )
    @test all(
        unit -> unit.diagnostics.shell_projection_lowdin_realization,
        be2_shared_pqs_units,
    )
    @test all(
        unit -> !unit.diagnostics.raw_product_box_operator_contract,
        be2_shared_pqs_units,
    )
    @test all(
        unit -> !unit.diagnostics.ida_weight_division_allowed,
        be2_shared_pqs_units,
    )
    @test all(
        unit -> unit.diagnostics.retained_weight_semantics ==
                :not_positive_quadrature_weights,
        be2_shared_pqs_units,
    )
    @test length(be2_current_route_inventory.source_fixtures.shared_pqs) == 3
    @test be2_current_route_pair_inventory.object_kind ==
          :pqs_current_route_retained_pair_inventory_fixture
    @test be2_current_route_pair_inventory.unit_inventory ===
          be2_current_route_inventory
    @test length(be2_current_route_pair_inventory.pairs) == 36
    @test be2_current_route_pair_inventory.counts.pair_count == 36
    @test be2_current_route_pair_inventory.diagnostics.pair_count == 36
    @test be2_current_route_pair_inventory.diagnostics.expected_pair_count == 36
    @test be2_current_route_pair_inventory.diagnostics.unit_count == 8
    @test be2_current_route_pair_inventory.counts.raw_box_pqs_active == 0
    @test be2_current_route_pair_inventory.diagnostics.raw_box_pqs_active_pair_policy_count == 0
    @test be2_current_route_pair_inventory.counts.support_local_oracle_for_shell_realization == 21
    @test be2_current_route_pair_inventory.diagnostics.support_local_oracle_for_shell_realization_pair_count == 21
    @test be2_current_route_pair_inventory.diagnostics.shell_row_oracle_pair_count == 21
    @test be2_current_route_pair_inventory.diagnostics.shell_realized_pqs_pairs_are_oracle_only
    @test all(
        pair -> !pair.raw_box_pqs_active_pair_policy,
        be2_current_route_pair_inventory.pairs,
    )
    @test count(pair -> pair.shell_row_oracle_only, be2_current_route_pair_inventory.pairs) == 21
    @test count(pair -> pair.active_algorithmic_policy, be2_current_route_pair_inventory.pairs) == 15
    @test !be2_current_route_pair_inventory.diagnostics.whole_route_safe_term_matrix_consumer
    current_route_safe_terms = CCPM._pqs_current_route_safe_term_matrices(
        pqs_construction,
        contact_safe_term_metrics;
        inventory = current_route_inventory,
        pair_inventory = current_route_pair_inventory,
    )
    @test current_route_safe_terms.object_kind ==
          :pqs_current_route_safe_term_matrices_fixture
    @test current_route_safe_terms.status == :private_diagnostic_only
    @test current_route_safe_terms.terms == (
        :overlap,
        :position_x,
        :position_y,
        :position_z,
        :x2_x,
        :x2_y,
        :x2_z,
        :kinetic,
    )
    @test current_route_safe_terms.global_max_error <= 1.0e-12
    @test current_route_safe_terms.diagnostics.private_diagnostic_only
    @test current_route_safe_terms.diagnostics.finite_output
    @test :product_doside_unit in pqs_source_descriptor.non_contracts
    @test :dense_full_parent_fallback in pqs_source_descriptor.non_contracts
    @test pqs_source_descriptor.diagnostics.metadata_only
    @test !pqs_source_descriptor.active_consumption.fixed_block_sidecar_installed
    @test !pqs_source_descriptor.active_consumption.metric_packet_consumes
    @test !pqs_source_descriptor.active_consumption.by_center_consumes
    @test pqs_diagnostics.region_builds[2].retained_count == 114
    @test pqs_diagnostics.fixed_dimension == 487
    @test !isnothing(pqs_construction.sequence.packet)
    @test all(isfinite, pqs_construction.sequence.packet.overlap)
    @test all(isfinite, pqs_construction.sequence.packet.kinetic)
    @test all(isfinite, pqs_construction.sequence.packet.weights)
    @test all(isfinite, pqs_construction.sequence.packet.gaussian_sum)
    @test all(isfinite, pqs_construction.sequence.packet.pair_sum)
    @test norm(pqs_construction.sequence.packet.overlap - I, Inf) < 1.0e-8
    pqs_fixed_block =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_fixed_block(
            pqs_construction,
        )
    @test size(pqs_fixed_block.coefficient_matrix) == (7 * 7 * 15, 487)
    @test all(isfinite, pqs_fixed_block.overlap)
    @test all(isfinite, pqs_fixed_block.kinetic)
    @test all(isfinite, pqs_fixed_block.weights)
    @test all(isfinite, pqs_fixed_block.gaussian_sum)
    @test all(isfinite, pqs_fixed_block.pair_sum)
    @test norm(pqs_fixed_block.overlap - I, Inf) < 1.0e-8
    @test pqs_fixed_block.staged_by_center_sidecar[] === nothing
    current_route_authority_comparison =
        CCPM._pqs_current_route_safe_term_authority_comparison(
            pqs_construction,
            contact_safe_term_metrics;
            inventory = current_route_inventory,
            pair_inventory = current_route_pair_inventory,
            safe_terms = current_route_safe_terms,
            fixed_block = pqs_fixed_block,
        )
    @test current_route_authority_comparison.object_kind ==
          :pqs_current_route_safe_term_authority_comparison_fixture
    @test current_route_authority_comparison.status == :private_diagnostic_only
    @test current_route_authority_comparison.terms == current_route_safe_terms.terms
    @test isempty(current_route_authority_comparison.unavailable_terms)
    @test current_route_authority_comparison.max_authority_error <= 1.0e-8
    @test current_route_authority_comparison.diagnostics.private_diagnostic_only
    @test current_route_authority_comparison.diagnostics.finite_output
    QWCS = GaussletBases.CartesianQWOperatorCarriedSpaces
    pqs_receipt = @test_logs min_level = Logging.Warn QWCS.cartesian_qw_operator_construction_receipt(
        pqs_fixed_block;
        nuclear_charges = [1.0, 1.0],
        nuclear_term_storage = :total_only,
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :pgdg_localized_experimental,
        expansion = expansion,
    )
    pqs_operators = QWCS.qw_operator_construction_receipt_operators(pqs_receipt)
    @test pqs_operators.gausslet_backend == :pgdg_localized_experimental
    @test pqs_operators.interaction_treatment == :ggt_nearest
    @test pqs_operators.nuclear_term_storage == :total_only
    @test pqs_operators.gausslet_count == 487
    @test pqs_operators.residual_count == 0
    @test size(pqs_operators.overlap) == (487, 487)
    @test size(pqs_operators.one_body_hamiltonian) == (487, 487)
    @test size(pqs_operators.interaction_matrix) == (487, 487)
    @test all(isfinite, pqs_operators.overlap)
    @test all(isfinite, pqs_operators.one_body_hamiltonian)
    @test all(isfinite, pqs_operators.interaction_matrix)
    @test norm(pqs_operators.overlap - I, Inf) < 1.0e-8
    @test norm(
        pqs_operators.one_body_hamiltonian - transpose(pqs_operators.one_body_hamiltonian),
        Inf,
    ) < 1.0e-12
    @test norm(
        pqs_operators.interaction_matrix - transpose(pqs_operators.interaction_matrix),
        Inf,
    ) < 1.0e-12
    @test_throws ArgumentError GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction(
        basis,
        bundles,
        policy;
        nside = 5,
        term_coefficients = Float64.(expansion.coefficients),
        packet_kernel = :factorized_direct,
        shared_shell_realization = :projected_q_shell,
    )

    shared_q5_policy = GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        policy.construction_plan;
        q_min = 4,
        atom_q = 4,
        atom_order = 4,
        shared_q = 5,
        shared_order = 5,
        contact_q = 4,
        contact_order = 4,
        outer_mismatch_q = 4,
        outer_mismatch_order = 4,
    )
    shared_q5_construction =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction(
            basis,
            bundles,
            shared_q5_policy;
            nside = 5,
            term_coefficients = Float64.(expansion.coefficients),
            packet_kernel = :factorized_direct,
            build_sequence_packet = false,
        )
    shared_q5_diagnostics =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction_diagnostics(
            shared_q5_construction,
        )
    shared_q5_region_builds = shared_q5_diagnostics.region_builds
    @test shared_q5_diagnostics.recipe_label == :mixed_atom_cubic_shared_endcap_panel
    @test shared_q5_diagnostics.active_builder_consumes
    @test shared_q5_diagnostics.parent_dimension == 7 * 7 * 15
    @test shared_q5_diagnostics.fixed_dimension == 523
    @test shared_q5_diagnostics.support_coverage.coverage_ok
    @test [build.q for build in shared_q5_region_builds] == [4, 5, 4, 4, 4]
    @test [build.order for build in shared_q5_region_builds] == [4, 5, 4, 4, 4]
    @test [build.retained_count for build in shared_q5_region_builds] == [98, 150, 125, 125, 25]
    @test [build.column_range for build in shared_q5_region_builds] ==
          [1:98, 374:523, 99:223, 224:348, 349:373]
    @test shared_q5_region_builds[2].metadata.coefficient_contract == :product_doside
    @test shared_q5_region_builds[5].metadata.descriptor_scope == :middle_contact_cap
    @test shared_q5_diagnostics.metadata.q_policy == :atom_growth_endcap_panel_shared_q_variable
    @test !shared_q5_diagnostics.metadata.q4_acceptance_fixture
    @test shared_q5_diagnostics.metadata.non_shared_q_policy == :fixed_q4_order4
    @test shared_q5_diagnostics.metadata.shared_q_values == (5,)
    @test shared_q5_diagnostics.metadata.shared_order_values == (5,)
    @test isnothing(shared_q5_construction.sequence.packet)

    no_packet_readiness =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_readiness(
            construction,
        )
    no_packet_diagnostics =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_readiness_diagnostics(
            no_packet_readiness,
        )
    @test no_packet_diagnostics.parent_dimension == 7 * 7 * 15
    @test no_packet_diagnostics.fixed_dimension == 469
    @test no_packet_diagnostics.support_coverage.coverage_ok
    @test !no_packet_diagnostics.sequence_packet_available
    @test !no_packet_diagnostics.can_produce_fixed_block
    @test no_packet_diagnostics.fixed_block_missing_fields == [:sequence_packet]
    @test !no_packet_diagnostics.can_produce_nested_source
    @test :split_geometry in no_packet_diagnostics.nested_source_missing_fields
    @test no_packet_diagnostics.default_builders_unchanged

    packeted_construction =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction(
            basis,
            bundles,
            policy;
            nside = 5,
            term_coefficients = Float64.(expansion.coefficients),
            packet_kernel = :factorized_direct,
            build_sequence_packet = true,
        )
    readiness =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_readiness(
            packeted_construction;
            build_fixed_block = true,
        )
    readiness_diagnostics =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_readiness_diagnostics(
            readiness,
        )
    @test readiness isa
          GaussletBases._BondAlignedDiatomicHighOrderRecipeSourceReadiness3D
    @test readiness_diagnostics.sequence_packet_available
    @test readiness_diagnostics.overlap_available
    @test readiness_diagnostics.weights_available
    @test readiness_diagnostics.overlap_error < 1.0e-8
    @test readiness_diagnostics.can_produce_fixed_block
    @test isempty(readiness_diagnostics.fixed_block_missing_fields)
    @test readiness_diagnostics.fixed_block_built
    @test readiness_diagnostics.fixed_block_backend == :unknown
    @test readiness_diagnostics.fixed_block_dimension == 469
    @test readiness_diagnostics.fixed_block_support_count == 7 * 7 * 15
    @test readiness_diagnostics.staged_by_center_sidecar_available
    @test !readiness_diagnostics.can_produce_nested_source
    @test :child_sequences in readiness_diagnostics.nested_source_missing_fields
    fixed_block =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_fixed_block(
            packeted_construction,
        )
    @test fixed_block isa GaussletBases._NestedFixedBlock3D
    @test size(fixed_block.coefficient_matrix) == (7 * 7 * 15, 469)
    @test norm(fixed_block.overlap - I, Inf) < 1.0e-8
    @test all(isfinite, fixed_block.weights)
    @test fixed_block.staged_by_center_sidecar[] isa
          GaussletBases._CartesianNestedProductStagedByCenterSidecar3D

    receipt = @test_logs min_level = Logging.Warn QWCS.cartesian_qw_operator_construction_receipt(
        fixed_block;
        nuclear_charges = [1.0, 1.0],
        nuclear_term_storage = :total_only,
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :pgdg_localized_experimental,
        expansion = expansion,
    )
    operators = QWCS.qw_operator_construction_receipt_operators(receipt)
    @test operators.gausslet_backend == :pgdg_localized_experimental
    @test operators.interaction_treatment == :ggt_nearest
    @test operators.gausslet_count == 469
    @test operators.residual_count == 0
    @test size(operators.overlap) == (469, 469)
    @test size(operators.one_body_hamiltonian) == (469, 469)
    @test size(operators.interaction_matrix) == (469, 469)
    @test all(isfinite, operators.overlap)
    @test all(isfinite, operators.one_body_hamiltonian)
    @test all(isfinite, operators.interaction_matrix)
    @test norm(operators.overlap - I, Inf) < 1.0e-8
    @test norm(
        operators.one_body_hamiltonian - transpose(operators.one_body_hamiltonian),
        Inf,
    ) < 1.0e-12
    @test norm(
        operators.interaction_matrix - transpose(operators.interaction_matrix),
        Inf,
    ) < 1.0e-12

    shared_q5_packeted_construction =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction(
            basis,
            bundles,
            shared_q5_policy;
            nside = 5,
            term_coefficients = Float64.(expansion.coefficients),
            packet_kernel = :factorized_direct,
            build_sequence_packet = true,
        )
    shared_q5_readiness =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_readiness(
            shared_q5_packeted_construction;
            build_fixed_block = true,
        )
    shared_q5_readiness_diagnostics =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_readiness_diagnostics(
            shared_q5_readiness,
        )
    @test shared_q5_readiness_diagnostics.sequence_packet_available
    @test shared_q5_readiness_diagnostics.overlap_available
    @test shared_q5_readiness_diagnostics.weights_available
    @test shared_q5_readiness_diagnostics.overlap_error < 1.0e-8
    @test shared_q5_readiness_diagnostics.can_produce_fixed_block
    @test isempty(shared_q5_readiness_diagnostics.fixed_block_missing_fields)
    @test shared_q5_readiness_diagnostics.fixed_block_built
    @test shared_q5_readiness_diagnostics.fixed_block_dimension == 523
    @test shared_q5_readiness_diagnostics.fixed_block_support_count == 7 * 7 * 15
    @test shared_q5_readiness_diagnostics.staged_by_center_sidecar_available
    shared_q5_fixed_block =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_fixed_block(
            shared_q5_packeted_construction,
        )
    @test shared_q5_fixed_block isa GaussletBases._NestedFixedBlock3D
    @test size(shared_q5_fixed_block.coefficient_matrix) == (7 * 7 * 15, 523)
    @test norm(shared_q5_fixed_block.overlap - I, Inf) < 1.0e-8
    @test all(isfinite, shared_q5_fixed_block.weights)
    @test shared_q5_fixed_block.staged_by_center_sidecar[] isa
          GaussletBases._CartesianNestedProductStagedByCenterSidecar3D

    shared_q5_receipt = @test_logs min_level = Logging.Warn QWCS.cartesian_qw_operator_construction_receipt(
        shared_q5_fixed_block;
        nuclear_charges = [1.0, 1.0],
        nuclear_term_storage = :total_only,
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :pgdg_localized_experimental,
        expansion = expansion,
    )
    shared_q5_operators = QWCS.qw_operator_construction_receipt_operators(shared_q5_receipt)
    @test shared_q5_operators.gausslet_backend == :pgdg_localized_experimental
    @test shared_q5_operators.interaction_treatment == :ggt_nearest
    @test shared_q5_operators.gausslet_count == 523
    @test shared_q5_operators.residual_count == 0
    @test size(shared_q5_operators.overlap) == (523, 523)
    @test size(shared_q5_operators.one_body_hamiltonian) == (523, 523)
    @test size(shared_q5_operators.interaction_matrix) == (523, 523)
    @test all(isfinite, shared_q5_operators.overlap)
    @test all(isfinite, shared_q5_operators.one_body_hamiltonian)
    @test all(isfinite, shared_q5_operators.interaction_matrix)
    @test norm(shared_q5_operators.overlap - I, Inf) < 1.0e-8
    @test norm(
        shared_q5_operators.one_body_hamiltonian -
        transpose(shared_q5_operators.one_body_hamiltonian),
        Inf,
    ) < 1.0e-12
    @test norm(
        shared_q5_operators.interaction_matrix -
        transpose(shared_q5_operators.interaction_matrix),
        Inf,
    ) < 1.0e-12

    shared_q6_policy = GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        policy.construction_plan;
        q_min = 4,
        atom_q = 4,
        atom_order = 4,
        shared_q = 6,
        shared_order = 6,
        contact_q = 4,
        contact_order = 4,
        outer_mismatch_q = 4,
        outer_mismatch_order = 4,
    )
    shared_q6_packeted_construction =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction(
            basis,
            bundles,
            shared_q6_policy;
            nside = 5,
            term_coefficients = Float64.(expansion.coefficients),
            packet_kernel = :factorized_direct,
            build_sequence_packet = true,
        )
    shared_q6_diagnostics =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction_diagnostics(
            shared_q6_packeted_construction,
        )
    shared_q6_region_builds = shared_q6_diagnostics.region_builds
    @test shared_q6_diagnostics.recipe_label == :mixed_atom_cubic_shared_endcap_panel
    @test shared_q6_diagnostics.active_builder_consumes
    @test shared_q6_diagnostics.parent_dimension == 7 * 7 * 15
    @test shared_q6_diagnostics.fixed_dimension == 589
    @test shared_q6_diagnostics.support_coverage.coverage_ok
    @test [build.q for build in shared_q6_region_builds] == [4, 6, 4, 4, 4]
    @test [build.order for build in shared_q6_region_builds] == [4, 6, 4, 4, 4]
    @test [build.retained_count for build in shared_q6_region_builds] ==
          [98, 216, 125, 125, 25]
    @test [build.column_range for build in shared_q6_region_builds] ==
          [1:98, 374:589, 99:223, 224:348, 349:373]
    @test shared_q6_region_builds[2].metadata.coefficient_contract == :product_doside
    @test all(
        build.retained_count == build.built_support_count for
        build in shared_q6_region_builds[[1, 3, 4, 5]]
    )
    @test shared_q6_diagnostics.metadata.non_shared_q_policy == :fixed_q4_order4
    @test shared_q6_diagnostics.metadata.shared_q_values == (6,)
    @test shared_q6_diagnostics.metadata.shared_order_values == (6,)

    shared_q6_readiness =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_readiness(
            shared_q6_packeted_construction;
            build_fixed_block = true,
        )
    shared_q6_readiness_diagnostics =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_readiness_diagnostics(
            shared_q6_readiness,
        )
    @test shared_q6_readiness_diagnostics.sequence_packet_available
    @test shared_q6_readiness_diagnostics.overlap_available
    @test shared_q6_readiness_diagnostics.weights_available
    @test shared_q6_readiness_diagnostics.overlap_error < 1.0e-8
    @test shared_q6_readiness_diagnostics.can_produce_fixed_block
    @test isempty(shared_q6_readiness_diagnostics.fixed_block_missing_fields)
    @test shared_q6_readiness_diagnostics.fixed_block_built
    @test shared_q6_readiness_diagnostics.fixed_block_dimension == 589
    @test shared_q6_readiness_diagnostics.fixed_block_support_count == 7 * 7 * 15
    @test shared_q6_readiness_diagnostics.staged_by_center_sidecar_available
    shared_q6_fixed_block = shared_q6_readiness.fixed_block
    @test shared_q6_fixed_block isa GaussletBases._NestedFixedBlock3D
    @test size(shared_q6_fixed_block.coefficient_matrix) == (7 * 7 * 15, 589)
    @test norm(shared_q6_fixed_block.overlap - I, Inf) < 1.0e-8
    @test all(isfinite, shared_q6_fixed_block.weights)
    @test shared_q6_fixed_block.staged_by_center_sidecar[] isa
          GaussletBases._CartesianNestedProductStagedByCenterSidecar3D

    shared_q6_receipt = @test_logs min_level = Logging.Warn QWCS.cartesian_qw_operator_construction_receipt(
        shared_q6_fixed_block;
        nuclear_charges = [1.0, 1.0],
        nuclear_term_storage = :total_only,
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :pgdg_localized_experimental,
        expansion = expansion,
    )
    shared_q6_operators = QWCS.qw_operator_construction_receipt_operators(shared_q6_receipt)
    @test shared_q6_operators.gausslet_backend == :pgdg_localized_experimental
    @test shared_q6_operators.interaction_treatment == :ggt_nearest
    @test shared_q6_operators.gausslet_count == 589
    @test shared_q6_operators.residual_count == 0
    @test size(shared_q6_operators.overlap) == (589, 589)
    @test size(shared_q6_operators.one_body_hamiltonian) == (589, 589)
    @test size(shared_q6_operators.interaction_matrix) == (589, 589)
    @test all(isfinite, shared_q6_operators.overlap)
    @test all(isfinite, shared_q6_operators.one_body_hamiltonian)
    @test all(isfinite, shared_q6_operators.interaction_matrix)
    @test norm(shared_q6_operators.overlap - I, Inf) < 1.0e-8
    @test norm(
        shared_q6_operators.one_body_hamiltonian -
        transpose(shared_q6_operators.one_body_hamiltonian),
        Inf,
    ) < 1.0e-12
    @test norm(
        shared_q6_operators.interaction_matrix -
        transpose(shared_q6_operators.interaction_matrix),
        Inf,
    ) < 1.0e-12

    shared_q7_policy = GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        policy.construction_plan;
        q_min = 4,
        atom_q = 4,
        atom_order = 4,
        shared_q = 7,
        shared_order = 7,
        contact_q = 4,
        contact_order = 4,
        outer_mismatch_q = 4,
        outer_mismatch_order = 4,
    )
    shared_q7_error = try
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction(
            basis,
            bundles,
            shared_q7_policy;
            nside = 5,
            term_coefficients = Float64.(expansion.coefficients),
            packet_kernel = :factorized_direct,
            build_sequence_packet = false,
        )
        nothing
    catch err
        err
    end
    @test shared_q7_error isa ArgumentError
    @test occursin(
        "nested doside retained_count must not exceed the interval size",
        sprint(showerror, shared_q7_error),
    )

    q_row_expectations = (
        (q = 4, fixed_dimension = 469, retained = (98, 96, 125, 125, 25), shared = 96),
        (q = 5, fixed_dimension = 523, retained = (98, 150, 125, 125, 25), shared = 150),
        (q = 6, fixed_dimension = 589, retained = (98, 216, 125, 125, 25), shared = 216),
    )
    multi_shared_region_builds = (
        (role = :outer_mismatch_shared_molecular_shell, retained_count = 12),
        (role = :regular_shared_molecular_shell, retained_count = 24),
        (role = :regular_shared_molecular_shell, retained_count = 36),
        (role = :left_atom_box, retained_count = 64),
    )
    @test QWCS._nested_q_row_shared_retained_counts(multi_shared_region_builds) ==
          (24, 36)
    @test QWCS._nested_q_row_shared_retained_count(multi_shared_region_builds) ===
          nothing
    @test QWCS._nested_q_row_shared_retained_counts(()) == ()
    @test QWCS._nested_q_row_shared_retained_count(()) === nothing
    for expected in q_row_expectations
        q_row_receipt = @test_logs min_level = Logging.Warn QWCS._nested_bond_aligned_diatomic_high_order_q_row_route_receipt(
            basis;
            shared_q = expected.q,
            shared_order = expected.q,
            protected_atom_side_count = 5,
            q_min = 4,
            nside = 5,
            expansion = expansion,
            nuclear_charges = [1.0, 1.0],
            nuclear_term_storage = :total_only,
            interaction_treatment = :ggt_nearest,
            gausslet_backend = :pgdg_localized_experimental,
        )
        q_row_diagnostics =
            QWCS._nested_bond_aligned_diatomic_high_order_q_row_route_diagnostics(
                q_row_receipt,
            )
        @test q_row_diagnostics.route_label ==
              :bond_aligned_diatomic_high_order_q_row_route
        @test q_row_diagnostics.shared_q == expected.q
        @test q_row_diagnostics.shared_order == expected.q
        @test q_row_diagnostics.non_shared_q_policy == :fixed_q4_order4
        @test q_row_diagnostics.parent_dimension == 7 * 7 * 15
        @test q_row_diagnostics.fixed_dimension == expected.fixed_dimension
        @test q_row_diagnostics.retained_counts_by_region == expected.retained
        @test q_row_diagnostics.shared_retained_counts == (expected.shared,)
        @test q_row_diagnostics.shared_retained_count == expected.shared
        @test q_row_diagnostics.overlap_error < 1.0e-8
        @test q_row_diagnostics.staged_sidecar_available
        @test q_row_diagnostics.backend == :pgdg_localized_experimental
        @test q_row_diagnostics.residual_count == 0
        @test q_row_diagnostics.gausslet_count == expected.fixed_dimension
        @test q_row_diagnostics.dense_parent_matrix_used == false
    end

    for expected in q_row_expectations
        fixture_receipt = @test_logs min_level = Logging.Warn QWCS._nested_bond_aligned_homonuclear_high_order_q_row_fixture_receipt(
            bond_length = 5.0,
            core_spacing = 0.7,
            xmax_parallel = 8.0,
            xmax_transverse = 4.0,
            shared_q = expected.q,
            shared_order = expected.q,
            family = :G10,
            bond_axis = :z,
            nuclear_charge = 1.0,
            reference_spacing = 1.0,
            tail_spacing = 10.0,
            protected_atom_side_count = 5,
            q_min = 4,
            nside = 5,
            expansion = expansion,
            nuclear_term_storage = :total_only,
            interaction_treatment = :ggt_nearest,
            gausslet_backend = :pgdg_localized_experimental,
        )
        fixture_diagnostics =
            QWCS._nested_bond_aligned_homonuclear_high_order_q_row_fixture_diagnostics(
                fixture_receipt,
            )
        fixture_provenance =
            QWCS._nested_bond_aligned_homonuclear_high_order_q_row_fixture_provenance(
                fixture_receipt,
            )
        route_diagnostics = fixture_diagnostics.q_row_route_diagnostics
        @test fixture_diagnostics.route_label ==
              :bond_aligned_homonuclear_high_order_q_row_fixture
        @test fixture_diagnostics.receipt_contract ==
              :construct_homonuclear_basis_then_delegate_q_row_route
        @test fixture_diagnostics.basis_constructor ==
              :bond_aligned_homonuclear_qw_basis
        @test fixture_diagnostics.family == :G10
        @test fixture_diagnostics.bond_length == 5.0
        @test fixture_diagnostics.core_spacing == 0.7
        @test fixture_diagnostics.xmax_parallel == 8.0
        @test fixture_diagnostics.xmax_transverse == 4.0
        @test fixture_diagnostics.bond_axis == :z
        @test fixture_diagnostics.reference_spacing == 1.0
        @test fixture_diagnostics.tail_spacing == 10.0
        @test fixture_diagnostics.nuclear_charge == 1.0
        @test fixture_diagnostics.basis_nuclear_charges == (1.0, 1.0)
        @test fixture_receipt.basis.nuclear_charges == [1.0, 1.0]
        @test fixture_diagnostics.basis_nuclei == ((0.0, 0.0, -2.5), (0.0, 0.0, 2.5))
        @test fixture_diagnostics.parent_axis_counts == (7, 7, 15)
        @test fixture_diagnostics.parent_dimension == 7 * 7 * 15
        @test fixture_diagnostics.flat_index_convention.one_based
        @test fixture_diagnostics.flat_index_convention.fastest_axis == :z
        @test fixture_diagnostics.flat_index_convention.slowest_axis == :x
        @test fixture_diagnostics.shared_q == expected.q
        @test fixture_diagnostics.shared_order == expected.q
        @test fixture_diagnostics.q_min == 4
        @test fixture_diagnostics.protected_atom_side_count == 5
        @test fixture_diagnostics.nside == 5
        @test fixture_diagnostics.fixed_dimension == expected.fixed_dimension
        @test fixture_diagnostics.retained_counts_by_region == expected.retained
        @test fixture_diagnostics.shared_retained_counts == (expected.shared,)
        @test fixture_diagnostics.shared_retained_count == expected.shared
        @test route_diagnostics.route_label ==
              :bond_aligned_diatomic_high_order_q_row_route
        @test route_diagnostics.parent_dimension == fixture_diagnostics.parent_dimension
        @test route_diagnostics.fixed_dimension == expected.fixed_dimension
        @test route_diagnostics.backend == :pgdg_localized_experimental
        @test fixture_provenance.charge_policy == :basis_nuclear_charges_only
        @test fixture_provenance.homonuclear_only
        @test !fixture_provenance.heteronuclear_support
        @test !fixture_provenance.public_api
        @test !fixture_provenance.science_validation
    end

    @test_throws ArgumentError QWCS._nested_bond_aligned_diatomic_high_order_q_row_route_receipt(
        basis;
        shared_q = 4,
        expansion = expansion,
        gausslet_backend = :numerical_reference,
    )
    q_row_q7_error = try
        QWCS._nested_bond_aligned_diatomic_high_order_q_row_route_receipt(
            basis;
            shared_q = 7,
            shared_order = 7,
            expansion = expansion,
            gausslet_backend = :pgdg_localized_experimental,
        )
        nothing
    catch err
        err
    end
    @test q_row_q7_error isa ArgumentError
    @test occursin(
        "nested doside retained_count must not exceed the interval size",
        sprint(showerror, q_row_q7_error),
    )

    @test_throws ArgumentError QWCS._nested_bond_aligned_homonuclear_high_order_q_row_fixture_receipt(
        bond_length = 5.0,
        core_spacing = 0.7,
        xmax_parallel = 8.0,
        xmax_transverse = 4.0,
        shared_q = 4,
        expansion = expansion,
        gausslet_backend = :numerical_reference,
    )
    @test_throws UndefKeywordError QWCS._nested_bond_aligned_homonuclear_high_order_q_row_fixture_receipt(
        bond_length = 5.0,
        core_spacing = 0.7,
        xmax_transverse = 4.0,
        shared_q = 4,
        expansion = expansion,
    )
    fixture_q7_error = try
        QWCS._nested_bond_aligned_homonuclear_high_order_q_row_fixture_receipt(
            bond_length = 5.0,
            core_spacing = 0.7,
            xmax_parallel = 8.0,
            xmax_transverse = 4.0,
            shared_q = 7,
            shared_order = 7,
            expansion = expansion,
            gausslet_backend = :pgdg_localized_experimental,
        )
        nothing
    catch err
        err
    end
    @test fixture_q7_error isa ArgumentError
    @test occursin(
        "nested doside retained_count must not exceed the interval size",
        sprint(showerror, fixture_q7_error),
    )

    annulus_policy = GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        policy.construction_plan;
        q_min = 4,
        atom_q = 4,
        shared_q = 5,
        contact_q = 4,
        outer_mismatch_q = 4,
        shared_exterior_family = :transverse_annulus_exterior,
    )
    @test_throws ArgumentError GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction(
        basis,
        bundles,
        annulus_policy;
        nside = 5,
        term_coefficients = Float64.(expansion.coefficients),
        packet_kernel = :factorized_direct,
        build_sequence_packet = false,
    )

    atom_q5_policy = GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        policy.construction_plan;
        q_min = 4,
        atom_q = 5,
        atom_order = 5,
        shared_q = 5,
        shared_order = 5,
        contact_q = 4,
        outer_mismatch_q = 4,
    )
    @test_throws ArgumentError GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction(
        basis,
        bundles,
        atom_q5_policy;
        nside = 5,
        term_coefficients = Float64.(expansion.coefficients),
        packet_kernel = :factorized_direct,
        build_sequence_packet = false,
    )
end

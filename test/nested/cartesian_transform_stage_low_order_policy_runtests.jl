using Test
using GaussletBases

function _cartesian_transform_stage_low_order_policy_fixture(;
    probe_parent_axis_construction = :auto,
    atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
    parent_axis_counts = (x = 9, y = 7, z = 9),
)
    system_inputs = (;
        atom_symbols = ("Be", "Be"),
        nuclear_charges = (4, 4),
        atom_locations,
        radius = 15.0,
        parent_axis_counts,
        map_backend = :pgdg_localized_experimental,
    )
    spacing_inputs = (;
        q = 5,
        n_s = 5,
        reference_spacing = 1.0,
        tail_spacing = 10.0,
        q_to_core_spacing_rule = :standard_pqs_ns_equals_q,
        core_spacing = 0.15,
    )
    parent_inputs = (;
        probe_parent_axis_construction,
        parent_axis_probe_backend = :pgdg_localized_experimental,
        parent_axis_probe_family = :G10,
    )
    route_probe_inputs = (;
        probe_raw_product_box_plans = :auto,
        raw_product_box_probe_backend = :pgdg_localized_experimental,
    )
    route_inputs = (;
        route_family = :white_lindsey_low_order,
        route_kind = :be2_cartesian_transform_stage_low_order_policy,
        route_shape = (:pqs_left, :product, :pqs_right),
        product_body_rule = :centered_single_z_slab,
        pqs_retained_rule = :boundary_comx_product_mode_selection,
        product_retained_rule = :product_doside_retained_unit,
        terms = (
            :overlap,
            :position_x,
            :position_y,
            :position_z,
            :x2_x,
            :x2_y,
            :x2_z,
            :kinetic,
        ),
        pair_factor_normalization = :density_normalized,
        support_dense_direct_allowed = false,
        reference_only_authorities =
            (:support_row_oracle, :dense_parent_projection),
        white_lindsey_route_shape =
            (:standard_cartesian_units, :low_order_comx_coarsening),
        white_lindsey_mapping_rule = :standard_unit_backbone_mapping_family,
        white_lindsey_nesting_rule = :unit_box_low_order_comx_coarsening,
        white_lindsey_retained_rule = :low_order_unit_comx_retained_basis,
        white_lindsey_operator_rule = :low_order_unit_operator_blocks,
        white_lindsey_benchmark_role =
            :published_cartesian_baseline_for_pqs_comparison,
    )

    system = GaussletBases.cartesian_system(system_inputs)
    recipe = GaussletBases.cartesian_recipe(route_inputs)
    parent = GaussletBases.cartesian_parent(
        system,
        spacing_inputs,
        parent_inputs,
        recipe,
    )
    return (; parent, spacing_inputs, recipe, route_probe_inputs)
end

@testset "cartesian transform stage carries selected low-order policy" begin
    fixture = _cartesian_transform_stage_low_order_policy_fixture()

    default_shells = GaussletBases.cartesian_shells(
        fixture.parent,
        fixture.spacing_inputs,
        fixture.recipe,
    )
    default_units = GaussletBases.cartesian_units(
        fixture.parent,
        default_shells,
        fixture.route_probe_inputs,
        fixture.recipe,
    )
    default_transforms =
        GaussletBases.cartesian_transforms(default_units, fixture.recipe)
    default_summary = default_transforms.low_order_transforms
    @test default_transforms.object_kind == :cartesian_transforms
    @test default_summary.object_kind ==
          :cartesian_transform_stage_low_order_summary
    @test default_summary.low_order_shellization_policy_resolved ==
          :legacy_diatomic_source
    @test default_summary.shellization_source ==
          :route_configured_bond_aligned_diatomic_source
    @test default_summary.unit_route_kind ==
          :legacy_diatomic_source_low_order_units
    @test default_summary.transform_route_kind ==
          :legacy_diatomic_source_low_order_transforms
    @test default_summary.legacy_source_transforms_selected
    @test default_summary.active_source_authority
    @test !default_summary.atom_growth_transforms_selected
    @test !default_summary.terminal_shellification_transforms_selected
    @test !default_summary.atom_growth_transform_contracts_available
    @test !default_summary.terminal_shellification_transform_contracts_available
    @test !default_summary.transform_contract_inventory_available
    @test default_summary.transform_contract_source ==
          :legacy_diatomic_source_summary
    @test !default_summary.coefficient_transforms_materialized
    @test !default_summary.coefficient_maps_materialized
    @test default_summary.transform_materialization_status ==
          :deferred_legacy_diatomic_source_transform_materialization
    @test !default_summary.retained_unit_dimensions_known
    @test !default_summary.retained_unit_ranges_known
    @test !default_summary.retained_dimension_known
    @test default_summary.transform_fields_preserved
    @test default_transforms.retained_units === default_units.retained_units
    @test hasproperty(default_transforms, :retained_counts)
    @test hasproperty(default_transforms, :ranges)
    @test hasproperty(default_transforms, :retained_dimension)
    @test default_transforms.transform_stage == :unit_retained_transforms_described

    atom_growth_shells = GaussletBases.cartesian_shells(
        fixture.parent,
        fixture.spacing_inputs,
        fixture.recipe;
        low_order_shellization_policy = :atom_growth_complete_rectangular,
    )
    atom_growth_units = GaussletBases.cartesian_units(
        fixture.parent,
        atom_growth_shells,
        fixture.route_probe_inputs,
        fixture.recipe,
    )
    atom_growth_transforms =
        GaussletBases.cartesian_transforms(atom_growth_units, fixture.recipe)
    atom_growth_summary = atom_growth_transforms.low_order_transforms
    @test atom_growth_transforms.low_order_transform_route_kind ==
          :atom_growth_complete_rectangular_low_order_transforms
    @test atom_growth_transforms.atom_growth_transforms_selected
    @test atom_growth_transforms.atom_growth_transform_contracts_available
    @test atom_growth_transforms.transform_contract_inventory_available
    @test atom_growth_transforms.transform_contract_source ==
          :atom_growth_plan_unit_inventory
    @test atom_growth_transforms.transform_contract_status ==
          :available_atom_growth_transform_contract_inventory
    @test atom_growth_transforms.cpb_contract_stage ==
          :construction_transform_contract
    @test atom_growth_transforms.pqs_transform_prototype_available
    @test atom_growth_transforms.transform_contracts_derive_from_lowering
    @test atom_growth_transforms.final_retained_units_are_pair_planning_inputs
    @test !atom_growth_transforms.coefficient_transforms_materialized
    @test !atom_growth_transforms.coefficient_maps_materialized
    @test !atom_growth_transforms.active_source_authority
    @test atom_growth_summary.low_order_shellization_policy_resolved ==
          :atom_growth_complete_rectangular
    @test atom_growth_summary.shellization_source ==
          :bond_aligned_diatomic_atom_growth_construction_plan
    @test atom_growth_summary.shellization_kind ==
          :atom_growth_complete_rectangular
    @test atom_growth_summary.unit_route_kind ==
          :atom_growth_complete_rectangular_low_order_units
    @test atom_growth_summary.transform_route_kind ==
          :atom_growth_complete_rectangular_low_order_transforms
    @test atom_growth_summary.atom_growth_transforms_selected
    @test !atom_growth_summary.terminal_shellification_transforms_selected
    @test !atom_growth_summary.legacy_source_transforms_selected
    @test atom_growth_summary.atom_growth_transform_contracts_available
    @test !atom_growth_summary.terminal_shellification_transform_contracts_available
    @test atom_growth_summary.transform_contract_inventory_available
    @test atom_growth_summary.transform_contract_source ==
          :atom_growth_plan_unit_inventory
    @test atom_growth_summary.transform_contract_status ==
          :available_atom_growth_transform_contract_inventory
    @test atom_growth_summary.cpb_contract_stage ==
          :construction_transform_contract
    @test atom_growth_summary.pqs_transform_prototype_available
    @test atom_growth_summary.transform_contracts_derive_from_lowering
    @test atom_growth_summary.final_retained_units_are_pair_planning_inputs
    @test atom_growth_summary.plan_authority
    @test !atom_growth_summary.active_source_authority
    @test !atom_growth_summary.legacy_source_authority
    @test !atom_growth_summary.coefficient_transforms_materialized
    @test !atom_growth_summary.coefficient_maps_materialized
    @test atom_growth_summary.transform_materialization_status ==
          :deferred_atom_growth_complete_rectangular_transform_materialization
    @test !atom_growth_summary.retained_unit_dimensions_known
    @test !atom_growth_summary.retained_unit_ranges_known
    @test !atom_growth_summary.retained_dimension_known
    @test atom_growth_summary.retained_dimension === nothing
    @test !atom_growth_summary.summary_only
    unit_inventory = atom_growth_units.plan_unit_inventory
    contract_inventory = atom_growth_summary.transform_contract_inventory
    @test contract_inventory.object_kind ==
          :cartesian_atom_growth_transform_contract_inventory
    @test contract_inventory.transform_contract_source ==
          :atom_growth_plan_unit_inventory
    @test contract_inventory.unit_keys == unit_inventory.unit_keys
    @test contract_inventory.unit_roles == unit_inventory.unit_roles
    @test contract_inventory.contract_count == unit_inventory.unit_count
    expected_contract_names = Tuple(
        role == :outer_mismatch_shared_molecular_shell ?
        :outer_mismatch_boundary_slab_set :
        role == :regular_shared_molecular_shell ?
        :adaptive_complete_shell_layer :
        role in (:left_atom_box, :right_atom_box) ?
        :atom_local_child_shellification_sequence :
        role == :contact_cap ?
        :direct_identity_selector :
        :unexpected_atom_growth_transform_contract
        for role in contract_inventory.unit_roles
    )
    @test contract_inventory.contract_names == expected_contract_names
    @test !(:unexpected_atom_growth_transform_contract in expected_contract_names)
    @test contract_inventory.source_backed_contract_count == 0
    @test !contract_inventory.coefficient_transforms_materialized
    @test !contract_inventory.coefficient_maps_materialized
    @test contract_inventory.lw_complete_shell_cpb_enumeration_available
    @test contract_inventory.lw_complete_shell_region_count == 4
    @test contract_inventory.lw_complete_shell_cpb_count == 104
    @test contract_inventory.lw_complete_shell_cpb_family_counts ==
          (facet_cpb = 24, edge_cpb = 48, corner_cpb = 32)
    @test contract_inventory.lw_complete_shell_enumeration_policy ==
          :white_lindsey_complete_shell_boundary_strata
    @test !contract_inventory.lw_complete_shell_coefficient_maps_materialized
    @test !contract_inventory.lw_complete_shell_operator_blocks_materialized
    @test !contract_inventory.lw_complete_shell_pair_operator_blocks_materialized
    @test !contract_inventory.lw_complete_shell_hamiltonian_data_materialized
    @test atom_growth_summary.lw_complete_shell_cpb_enumeration_available
    @test atom_growth_summary.lw_complete_shell_region_count == 4
    @test atom_growth_summary.lw_complete_shell_cpb_count == 104
    @test atom_growth_summary.lw_complete_shell_cpb_family_counts ==
          contract_inventory.lw_complete_shell_cpb_family_counts
    @test atom_growth_transforms.lw_complete_shell_cpb_enumeration_available
    @test atom_growth_transforms.lw_complete_shell_region_count == 4
    @test atom_growth_transforms.lw_complete_shell_cpb_count == 104
    @test atom_growth_transforms.lw_complete_shell_cpb_family_counts ==
          contract_inventory.lw_complete_shell_cpb_family_counts
    @test !atom_growth_transforms.lw_complete_shell_coefficient_maps_materialized
    @test !atom_growth_transforms.lw_complete_shell_operator_blocks_materialized
    @test !atom_growth_transforms.lw_complete_shell_pair_operator_blocks_materialized
    @test !atom_growth_transforms.lw_complete_shell_hamiltonian_data_materialized
    @test contract_inventory.pqs_transform_prototype_available
    @test contract_inventory.source_lowering_prototype_unit_key ==
          unit_inventory.pqs_lowering_prototype_unit_key
    pqs_transform_prototype = contract_inventory.pqs_transform_prototype
    pqs_lowering_prototype = unit_inventory.pqs_lowering_prototype
    @test pqs_transform_prototype.object_kind ==
          :cartesian_pqs_transform_metadata_prototype
    @test pqs_transform_prototype.status == :metadata_only_planned
    @test pqs_transform_prototype.source_lowering_prototype ===
          pqs_lowering_prototype
    @test pqs_transform_prototype.source_lowering_prototype_unit_key ==
          pqs_lowering_prototype.unit_key
    @test pqs_transform_prototype.unit_key == pqs_lowering_prototype.unit_key
    @test pqs_transform_prototype.transform_plan == (
        :source_retained_modes,
        :shell_projection,
        :lowdin_cleanup,
        :final_retained_unit,
    )
    @test pqs_transform_prototype.source_cpb === pqs_lowering_prototype.source_cpb
    @test pqs_transform_prototype.source_cpb.support_count ==
          prod(pqs_transform_prototype.source_cpb.dimensions)
    @test pqs_transform_prototype.source_cpb_support_count ==
          pqs_lowering_prototype.source_cpb_support_count
    @test pqs_transform_prototype.owned_support ===
          pqs_lowering_prototype.owned_support
    @test pqs_transform_prototype.owned_support_count ==
          pqs_lowering_prototype.owned_support_count
    @test pqs_transform_prototype.source_cpb_support_count !=
          pqs_transform_prototype.owned_support_count
    @test pqs_transform_prototype.intermediate_retained_space ===
          pqs_lowering_prototype.intermediate_retained_space
    @test pqs_transform_prototype.shell_realization ===
          pqs_lowering_prototype.shell_realization
    @test pqs_transform_prototype.final_retained_unit.object_kind ==
          :cartesian_final_retained_unit_contract
    @test pqs_transform_prototype.final_retained_unit.unit_key ==
          pqs_lowering_prototype.unit_key
    @test !pqs_transform_prototype.coefficient_transform_materialized
    @test !pqs_transform_prototype.coefficient_maps_materialized
    @test !pqs_transform_prototype.numerical_transform_materialized
    @test !pqs_transform_prototype.source_operator_blocks_materialized
    @test !pqs_transform_prototype.operator_blocks_materialized
    @test !pqs_transform_prototype.pair_operator_blocks_materialized
    @test !pqs_transform_prototype.hamiltonian_data_materialized
    @test !pqs_transform_prototype.dense_parent_space_operators_are_algorithm
    @test atom_growth_summary.pqs_transform_prototype ===
          pqs_transform_prototype
    @test atom_growth_transforms.pqs_transform_prototype ===
          pqs_transform_prototype
    @test !contract_inventory.retained_unit_dimensions_known
    @test !contract_inventory.retained_unit_ranges_known
    @test !contract_inventory.retained_dimension_known
    @test all(
        contract -> contract.cpb_contract_stage ==
                    :construction_transform_contract,
        contract_inventory.transform_contracts,
    )
    @test all(
        contract -> contract.owned_support.object_kind ==
                    :cartesian_owned_support_region3d,
        contract_inventory.transform_contracts,
    )
    @test all(
        contract -> contract.lowering_recipe.object_kind ==
                    :cartesian_cpb_lowering_recipe,
        contract_inventory.transform_contracts,
    )
    @test all(
        contract -> contract.lowering_recipe.shellification_authority_scope ==
                    :owned_support_only,
        contract_inventory.transform_contracts,
    )
    @test all(
        contract -> !contract.lowering_recipe.shellification_region_is_lowering_source,
        contract_inventory.transform_contracts,
    )
    @test all(
        contract -> contract.lowering_recipe.lowering_source_authority ==
                    :lowering_recipe_cpbs,
        contract_inventory.transform_contracts,
    )
    @test all(
        contract -> contract.intermediate_retained_space.object_kind ==
                    :cartesian_intermediate_retained_space_contract,
        contract_inventory.transform_contracts,
    )
    @test all(
        contract -> contract.shell_realization.object_kind ==
                    :cartesian_shell_realization_contract,
        contract_inventory.transform_contracts,
    )
    @test all(
        contract -> contract.final_retained_unit.object_kind ==
                    :cartesian_final_retained_unit_contract,
        contract_inventory.transform_contracts,
    )
    @test all(
        contract -> contract.final_unit_downstream_of_lowering,
        contract_inventory.transform_contracts,
    )
    @test all(
        contract -> contract.final_retained_unit.pair_planning_input,
        contract_inventory.transform_contracts,
    )
    @test all(
        contract -> !contract.source_backed,
        contract_inventory.transform_contracts,
    )
    @test all(
        contract -> contract.coefficient_transform_materialized == false,
        contract_inventory.transform_contracts,
    )
    @test all(
        contract -> contract.coefficient_map_materialized == false,
        contract_inventory.transform_contracts,
    )
    @test all(
        contract -> contract.retained_count_known == false,
        contract_inventory.transform_contracts,
    )
    @test all(
        contract -> contract.retained_range_known == false,
        contract_inventory.transform_contracts,
    )
    @test atom_growth_summary.transform_fields_preserved
    @test atom_growth_summary.route_skeleton_transform_inventory_source ==
          :route_skeleton_compatibility_fields
    @test atom_growth_transforms.retained_units === atom_growth_units.retained_units
    @test hasproperty(atom_growth_transforms, :retained_counts)
    @test hasproperty(atom_growth_transforms, :ranges)
    @test hasproperty(atom_growth_transforms, :retained_dimension)
    @test atom_growth_transforms.transform_stage ==
          :unit_retained_transforms_described

    terminal_fixture =
        _cartesian_transform_stage_low_order_policy_fixture(
            probe_parent_axis_construction = false,
            atom_locations = ((-4.0, 0.0, 0.0), (4.0, 0.0, 0.0)),
            parent_axis_counts = (x = 13, y = 7, z = 7),
        )
    terminal_shells = GaussletBases.cartesian_shells(
        terminal_fixture.parent,
        terminal_fixture.spacing_inputs,
        terminal_fixture.recipe;
        low_order_shellization_policy =
            :terminal_cartesian_shellification_geometry,
    )
    terminal_units = GaussletBases.cartesian_units(
        terminal_fixture.parent,
        terminal_shells,
        terminal_fixture.route_probe_inputs,
        terminal_fixture.recipe,
    )
    terminal_transforms =
        GaussletBases.cartesian_transforms(terminal_units, terminal_fixture.recipe)
    terminal_summary = terminal_transforms.low_order_transforms
    @test terminal_transforms.low_order_transform_route_kind ==
          :terminal_shellification_low_order_transforms
    @test terminal_transforms.terminal_shellification_transforms_selected
    @test terminal_transforms.terminal_shellification_transform_summary_available
    @test terminal_transforms.terminal_shellification_scaffold_available
    @test terminal_transforms.terminal_shellification_scaffold ===
          terminal_units.terminal_shellification_scaffold
    @test terminal_transforms.terminal_shellification_region_count ==
          terminal_units.terminal_shellification_region_count
    @test terminal_transforms.terminal_shellification_unit_inventory_available
    @test terminal_transforms.terminal_shellification_unit_inventory ===
          terminal_units.terminal_shellification_unit_inventory
    @test terminal_transforms.terminal_shellification_unit_count ==
          terminal_units.terminal_shellification_unit_count
    @test terminal_transforms.terminal_shellification_unit_keys ==
          terminal_units.terminal_shellification_unit_keys
    @test terminal_transforms.terminal_shellification_unit_roles ==
          terminal_units.terminal_shellification_unit_roles
    @test terminal_transforms.terminal_shellification_unit_kinds ==
          terminal_units.terminal_shellification_unit_kinds
    @test terminal_transforms.terminal_shellification_unit_support_counts ==
          terminal_units.terminal_shellification_unit_support_counts
    @test terminal_transforms.terminal_shellification_lowering_contract_inventory_available
    @test terminal_transforms.terminal_shellification_lowering_contract_inventory_status ==
          terminal_units.terminal_shellification_lowering_contract_inventory_status
    @test terminal_transforms.terminal_shellification_lowering_contract_inventory ===
          terminal_units.terminal_shellification_lowering_contract_inventory
    @test terminal_transforms.terminal_shellification_lowering_contract_count ==
          terminal_units.terminal_shellification_lowering_contract_count
    @test terminal_transforms.terminal_shellification_lowering_contract_kinds ==
          terminal_units.terminal_shellification_lowering_contract_kinds
    @test terminal_transforms.terminal_shellification_lowering_contract_kind_counts ==
          terminal_units.terminal_shellification_lowering_contract_kind_counts
    @test terminal_transforms.terminal_shellification_contract_counts_by_unit ==
          terminal_units.terminal_shellification_contract_counts_by_unit
    @test terminal_transforms.terminal_shellification_selected_lowering_contract_inventory_available
    @test terminal_transforms.terminal_shellification_selected_lowering_contract_inventory_status ==
          terminal_units.terminal_shellification_selected_lowering_contract_inventory_status
    @test terminal_transforms.terminal_shellification_selected_lowering_contract_inventory ===
          terminal_units.terminal_shellification_selected_lowering_contract_inventory
    @test terminal_transforms.terminal_shellification_selected_lowering_family ==
          terminal_units.terminal_shellification_selected_lowering_family
    @test terminal_transforms.terminal_shellification_selected_contract_count ==
          terminal_units.terminal_shellification_selected_contract_count
    @test terminal_transforms.terminal_shellification_selected_contract_kinds ==
          terminal_units.terminal_shellification_selected_contract_kinds
    @test terminal_transforms.terminal_shellification_selected_contract_kind_counts ==
          terminal_units.terminal_shellification_selected_contract_kind_counts
    @test terminal_transforms.terminal_shellification_selected_contract_counts_by_unit ==
          terminal_units.terminal_shellification_selected_contract_counts_by_unit
    @test terminal_transforms.terminal_shellification_all_units_have_exactly_one_selected_contract ==
          terminal_units.terminal_shellification_all_units_have_exactly_one_selected_contract
    @test terminal_transforms.terminal_shellification_unselected_contract_count ==
          terminal_units.terminal_shellification_unselected_contract_count
    @test terminal_transforms.terminal_shellification_unselected_contract_kinds ==
          terminal_units.terminal_shellification_unselected_contract_kinds
    @test terminal_transforms.terminal_shellification_lw_complete_shell_cpb_count ==
          terminal_units.terminal_shellification_lw_complete_shell_cpb_count
    @test terminal_transforms.terminal_shellification_lw_complete_shell_cpb_family_counts ==
          terminal_units.terminal_shellification_lw_complete_shell_cpb_family_counts
    terminal_lowering_inventory =
        terminal_transforms.terminal_shellification_lowering_contract_inventory
    selected_terminal_lowering_inventory =
        terminal_transforms.terminal_shellification_selected_lowering_contract_inventory
    @test terminal_lowering_inventory.lowering_contract_count >=
          terminal_transforms.terminal_shellification_unit_count
    @test all(
        entry.lowering_contract_count >= 1
        for entry in terminal_transforms.terminal_shellification_contract_counts_by_unit
    )
    @test selected_terminal_lowering_inventory.selected_contract_count ==
          terminal_transforms.terminal_shellification_unit_count
    @test selected_terminal_lowering_inventory.route_lowering_family ==
          :white_lindsey_low_order
    @test selected_terminal_lowering_inventory.all_units_have_exactly_one_selected_contract
    @test !terminal_lowering_inventory.final_retained_unit_inventory_available
    @test !terminal_lowering_inventory.pair_inventory_available
    @test !terminal_lowering_inventory.coefficient_maps_materialized
    @test !terminal_lowering_inventory.transform_contracts_materialized
    @test !terminal_lowering_inventory.retained_spaces_materialized
    @test !terminal_lowering_inventory.operator_blocks_materialized
    @test !terminal_lowering_inventory.pair_operator_blocks_materialized
    @test !terminal_lowering_inventory.hamiltonian_data_materialized
    @test !terminal_lowering_inventory.artifacts_materialized
    @test !selected_terminal_lowering_inventory.final_retained_unit_inventory_available
    @test !selected_terminal_lowering_inventory.pair_inventory_available
    @test !selected_terminal_lowering_inventory.coefficient_maps_materialized
    @test !selected_terminal_lowering_inventory.transform_contracts_materialized
    @test !selected_terminal_lowering_inventory.retained_spaces_materialized
    @test !selected_terminal_lowering_inventory.operator_blocks_materialized
    @test !selected_terminal_lowering_inventory.pair_operator_blocks_materialized
    @test !selected_terminal_lowering_inventory.hamiltonian_data_materialized
    @test !selected_terminal_lowering_inventory.artifacts_materialized
    @test !terminal_transforms.terminal_shellification_final_retained_unit_inventory_available
    @test !terminal_transforms.terminal_shellification_pair_inventory_available
    @test !terminal_transforms.terminal_shellification_transform_contracts_available
    @test terminal_transforms.terminal_shellification_transform_materialization_status ==
          :deferred_terminal_shellification_transform_contracts
    @test terminal_transforms.transform_contract_source ==
          :terminal_shellification_scaffold
    @test terminal_transforms.transform_contract_status ==
          :deferred_terminal_shellification_transform_contracts
    @test !terminal_transforms.transform_contract_inventory_available
    @test terminal_transforms.transform_contract_inventory === nothing
    @test !terminal_transforms.coefficient_transforms_materialized
    @test !terminal_transforms.coefficient_maps_materialized
    @test !terminal_transforms.active_source_authority
    @test terminal_transforms.terminal_shellification_central_gap_region_count ==
          terminal_units.terminal_shellification_central_gap_region_count
    @test terminal_transforms.terminal_shellification_central_midpoint_slab_count ==
          terminal_units.terminal_shellification_central_midpoint_slab_count
    @test terminal_transforms.terminal_shellification_central_distorted_product_box_count ==
          terminal_units.terminal_shellification_central_distorted_product_box_count
    @test terminal_transforms.terminal_shellification_central_gap_region_count == 3
    @test terminal_transforms.terminal_shellification_central_midpoint_slab_count ==
          3
    @test terminal_transforms.terminal_shellification_central_distorted_product_box_count ==
          0
    @test terminal_summary.object_kind ==
          :cartesian_transform_stage_low_order_summary
    @test terminal_summary.status ==
          :deferred_terminal_shellification_transform_contracts
    @test terminal_summary.low_order_shellization_policy_resolved ==
          :terminal_cartesian_shellification_geometry
    @test terminal_summary.shellization_source ==
          :terminal_cartesian_shellification_geometry
    @test terminal_summary.shellization_kind ==
          :terminal_cartesian_shellification_geometry
    @test terminal_summary.unit_route_kind ==
          :terminal_shellification_low_order_units
    @test terminal_summary.transform_route_kind ==
          :terminal_shellification_low_order_transforms
    @test terminal_summary.terminal_shellification_transforms_selected
    @test !terminal_summary.atom_growth_transforms_selected
    @test !terminal_summary.legacy_source_transforms_selected
    @test terminal_summary.terminal_shellification_transform_summary_available
    @test terminal_summary.terminal_shellification_scaffold_available
    @test terminal_summary.terminal_shellification_scaffold ===
          terminal_units.terminal_shellification_scaffold
    @test terminal_summary.terminal_shellification_region_count ==
          terminal_units.terminal_shellification_region_count
    @test terminal_summary.terminal_shellification_unit_inventory_available
    @test terminal_summary.terminal_shellification_unit_inventory ===
          terminal_units.terminal_shellification_unit_inventory
    @test terminal_summary.terminal_shellification_unit_count ==
          terminal_units.terminal_shellification_unit_count
    @test terminal_summary.terminal_shellification_unit_keys ==
          terminal_units.terminal_shellification_unit_keys
    @test terminal_summary.terminal_shellification_unit_roles ==
          terminal_units.terminal_shellification_unit_roles
    @test terminal_summary.terminal_shellification_unit_kinds ==
          terminal_units.terminal_shellification_unit_kinds
    @test terminal_summary.terminal_shellification_unit_support_counts ==
          terminal_units.terminal_shellification_unit_support_counts
    @test terminal_summary.terminal_shellification_lowering_contract_inventory_available
    @test terminal_summary.terminal_shellification_lowering_contract_inventory_status ==
          terminal_units.terminal_shellification_lowering_contract_inventory_status
    @test terminal_summary.terminal_shellification_lowering_contract_inventory ===
          terminal_units.terminal_shellification_lowering_contract_inventory
    @test terminal_summary.terminal_shellification_lowering_contract_count ==
          terminal_units.terminal_shellification_lowering_contract_count
    @test terminal_summary.terminal_shellification_lowering_contract_kinds ==
          terminal_units.terminal_shellification_lowering_contract_kinds
    @test terminal_summary.terminal_shellification_lowering_contract_kind_counts ==
          terminal_units.terminal_shellification_lowering_contract_kind_counts
    @test terminal_summary.terminal_shellification_contract_counts_by_unit ==
          terminal_units.terminal_shellification_contract_counts_by_unit
    @test terminal_summary.terminal_shellification_selected_lowering_contract_inventory_available
    @test terminal_summary.terminal_shellification_selected_lowering_contract_inventory_status ==
          terminal_units.terminal_shellification_selected_lowering_contract_inventory_status
    @test terminal_summary.terminal_shellification_selected_lowering_contract_inventory ===
          selected_terminal_lowering_inventory
    @test terminal_summary.terminal_shellification_selected_lowering_family ==
          terminal_units.terminal_shellification_selected_lowering_family
    @test terminal_summary.terminal_shellification_selected_contract_count ==
          terminal_units.terminal_shellification_selected_contract_count
    @test terminal_summary.terminal_shellification_selected_contract_kinds ==
          terminal_units.terminal_shellification_selected_contract_kinds
    @test terminal_summary.terminal_shellification_selected_contract_kind_counts ==
          terminal_units.terminal_shellification_selected_contract_kind_counts
    @test terminal_summary.terminal_shellification_selected_contract_counts_by_unit ==
          terminal_units.terminal_shellification_selected_contract_counts_by_unit
    @test terminal_summary.terminal_shellification_all_units_have_exactly_one_selected_contract ==
          terminal_units.terminal_shellification_all_units_have_exactly_one_selected_contract
    @test terminal_summary.terminal_shellification_unselected_contract_count ==
          terminal_units.terminal_shellification_unselected_contract_count
    @test terminal_summary.terminal_shellification_unselected_contract_kinds ==
          terminal_units.terminal_shellification_unselected_contract_kinds
    @test terminal_summary.terminal_shellification_lw_complete_shell_cpb_count ==
          terminal_units.terminal_shellification_lw_complete_shell_cpb_count
    @test terminal_summary.terminal_shellification_lw_complete_shell_cpb_family_counts ==
          terminal_units.terminal_shellification_lw_complete_shell_cpb_family_counts
    @test !terminal_summary.terminal_shellification_final_retained_unit_inventory_available
    @test !terminal_summary.terminal_shellification_pair_inventory_available
    @test !terminal_summary.terminal_shellification_transform_contracts_available
    @test terminal_summary.terminal_shellification_transform_materialization_status ==
          :deferred_terminal_shellification_transform_contracts
    @test !terminal_summary.transform_contract_inventory_available
    @test terminal_summary.transform_contract_inventory === nothing
    @test terminal_summary.transform_contract_count == 0
    @test terminal_summary.transform_contract_source ==
          :terminal_shellification_scaffold
    @test terminal_summary.transform_contract_status ==
          :deferred_terminal_shellification_transform_contracts
    @test terminal_summary.plan_authority
    @test !terminal_summary.active_source_authority
    @test !terminal_summary.legacy_source_authority
    @test !terminal_summary.coefficient_transforms_materialized
    @test !terminal_summary.coefficient_maps_materialized
    @test terminal_summary.transform_materialization_status ==
          :deferred_terminal_shellification_transform_contracts
    @test !terminal_summary.retained_unit_dimensions_known
    @test !terminal_summary.retained_unit_ranges_known
    @test !terminal_summary.retained_dimension_known
    @test terminal_summary.retained_dimension === nothing
    @test terminal_summary.summary_only
    @test all(
        record -> !record.owned_support_is_cpb,
        terminal_summary.terminal_shellification_unit_inventory.terminal_region_units,
    )
    @test all(
        record -> !record.shellification_region_is_cpb,
        terminal_summary.terminal_shellification_unit_inventory.terminal_region_units,
    )
    direct_contracts =
        filter(
            contract ->
                contract.lowering_contract_kind in (
                    :direct_core_identity_cpb,
                    :direct_slab_identity_cpb,
                    :direct_boundary_slab_identity_cpb,
                ),
            terminal_lowering_inventory.lowering_contracts,
        )
    @test all(contract.identity_like_source_contract for contract in direct_contracts)
    lw_contracts =
        filter(
            contract ->
                contract.lowering_contract_kind == :white_lindsey_boundary_strata,
            terminal_lowering_inventory.lowering_contracts,
        )
    if !isempty(lw_contracts)
        @test all(
            contract.source_cpb_family_counts ==
            (facet_cpb = 6, edge_cpb = 12, corner_cpb = 8)
            for contract in lw_contracts
        )
        @test all(contract.source_cpb_count == 26 for contract in lw_contracts)
    end
    pqs_contracts =
        filter(
            contract -> contract.lowering_contract_kind == :pqs_filled_source_cpb,
            terminal_lowering_inventory.lowering_contracts,
        )
    if !isempty(pqs_contracts)
        @test all(contract.source_cpb_count == 1 for contract in pqs_contracts)
        @test all(
            !contract.face_edge_corner_decomposition_required
            for contract in pqs_contracts
        )
    end
    @test !terminal_summary.lw_complete_shell_cpb_enumeration_available
    @test terminal_summary.lw_complete_shell_cpb_count == 0
    @test !terminal_summary.pqs_transform_prototype_available
    @test terminal_summary.pqs_transform_prototype === nothing
    @test terminal_summary.transform_fields_preserved
    @test terminal_transforms.retained_units === terminal_units.retained_units
    @test terminal_transforms.transform_stage ==
          :unit_retained_transforms_described
end

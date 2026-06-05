using Test
using GaussletBases

function _cartesian_transform_stage_low_order_policy_fixture()
    system_inputs = (;
        atom_symbols = ("Be", "Be"),
        nuclear_charges = (4, 4),
        atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
        radius = 15.0,
        parent_axis_counts = (x = 9, y = 7, z = 9),
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
        probe_parent_axis_construction = :auto,
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
    @test !default_summary.atom_growth_transform_contracts_available
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
    @test !atom_growth_summary.legacy_source_transforms_selected
    @test atom_growth_summary.atom_growth_transform_contracts_available
    @test atom_growth_summary.transform_contract_inventory_available
    @test atom_growth_summary.transform_contract_source ==
          :atom_growth_plan_unit_inventory
    @test atom_growth_summary.transform_contract_status ==
          :available_atom_growth_transform_contract_inventory
    @test atom_growth_summary.cpb_contract_stage ==
          :construction_transform_contract
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
end

using Test
using GaussletBases

function _cartesian_pair_stage_low_order_policy_fixture()
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
        route_kind = :be2_cartesian_pair_stage_low_order_policy,
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

function _cartesian_pair_stage_low_order_policy_pairs(fixture; policy = nothing)
    shells = isnothing(policy) ?
        GaussletBases.cartesian_shells(
            fixture.parent,
            fixture.spacing_inputs,
            fixture.recipe,
        ) :
        GaussletBases.cartesian_shells(
            fixture.parent,
            fixture.spacing_inputs,
            fixture.recipe;
            low_order_shellization_policy = policy,
        )
    units = GaussletBases.cartesian_units(
        fixture.parent,
        shells,
        fixture.route_probe_inputs,
        fixture.recipe,
    )
    transforms = GaussletBases.cartesian_transforms(units, fixture.recipe)
    pairs = GaussletBases.cartesian_pair_terms(
        units,
        transforms,
        fixture.recipe,
    )
    return (; shells, units, transforms, pairs)
end

@testset "cartesian pair stage carries selected low-order policy" begin
    fixture = _cartesian_pair_stage_low_order_policy_fixture()

    default_stages = _cartesian_pair_stage_low_order_policy_pairs(fixture)
    default_pairs = default_stages.pairs
    default_summary = default_pairs.low_order_pairs
    @test default_pairs.object_kind == :cartesian_pair_terms
    @test default_summary.object_kind == :cartesian_pair_stage_low_order_summary
    @test default_summary.low_order_shellization_policy_resolved ==
          :legacy_diatomic_source
    @test default_summary.shellization_source ==
          :route_configured_bond_aligned_diatomic_source
    @test default_summary.transform_route_kind ==
          :legacy_diatomic_source_low_order_transforms
    @test default_summary.pair_route_kind ==
          :legacy_diatomic_source_low_order_pairs
    @test default_summary.legacy_source_pairs_selected
    @test default_summary.active_source_authority
    @test !default_summary.atom_growth_pairs_selected
    @test !default_summary.pair_operator_blocks_materialized
    @test default_summary.pair_inventory_known
    @test default_summary.pair_inventory_source ==
          :route_skeleton_pair_entries_only
    @test default_summary.route_skeleton_pair_entry_count ==
          length(default_pairs.pair_entries)
    @test default_summary.route_skeleton_pair_family_counts ==
          default_pairs.pair_family_counts
    @test !default_summary.independent_atom_growth_pair_inventory_available
    @test !default_summary.route_core_pair_inventory_available
    @test default_summary.route_core_pair_inventory_status ==
          :not_selected_legacy_source_pairs
    @test default_summary.route_core_pair_count == 0
    @test isempty(default_summary.route_core_pair_keys)
    @test default_summary.pair_stage_fields_preserved
    @test default_pairs.pair_entries === default_stages.units.route_skeleton.pair_entries
    @test default_pairs.pair_family_counts ===
          default_stages.units.route_skeleton.pair_family_counts
    @test default_pairs.helper_by_pair_family ===
          default_stages.units.route_skeleton.helper_by_pair_family
    @test default_pairs.pair_operator_helper_by_family ===
          default_stages.units.route_skeleton.helper_by_pair_family
    @test isnothing(default_pairs.pair_helper_status_by_family)
    @test !default_pairs.route_core_pair_inventory_available
    @test default_pairs.route_core_pair_inventory_status ==
          :not_selected_legacy_source_pairs
    @test default_pairs.pair_stage == :pair_operator_terms_described

    atom_growth_stages = _cartesian_pair_stage_low_order_policy_pairs(
        fixture;
        policy = :atom_growth_complete_rectangular,
    )
    atom_growth_pairs = atom_growth_stages.pairs
    atom_growth_summary = atom_growth_pairs.low_order_pairs
    @test atom_growth_pairs.low_order_pair_route_kind ==
          :atom_growth_complete_rectangular_low_order_pairs
    @test atom_growth_pairs.atom_growth_pairs_selected
    @test !atom_growth_pairs.pair_operator_blocks_materialized
    @test !atom_growth_pairs.operator_pairs_materialized
    @test atom_growth_pairs.pair_inventory_source ==
          :atom_growth_unit_inventory
    @test atom_growth_pairs.independent_atom_growth_pair_inventory_available
    @test !atom_growth_pairs.active_source_authority
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
    @test atom_growth_summary.pair_route_kind ==
          :atom_growth_complete_rectangular_low_order_pairs
    @test atom_growth_summary.atom_growth_pairs_selected
    @test !atom_growth_summary.legacy_source_pairs_selected
    @test atom_growth_summary.plan_authority
    @test !atom_growth_summary.active_source_authority
    @test !atom_growth_summary.legacy_source_authority
    @test !atom_growth_summary.pair_operator_blocks_materialized
    @test !atom_growth_summary.operator_pairs_materialized
    @test atom_growth_summary.pair_inventory_known
    @test atom_growth_summary.pair_inventory_source ==
          :atom_growth_unit_inventory
    @test atom_growth_summary.independent_atom_growth_pair_inventory_available
    @test atom_growth_summary.route_core_pair_inventory_available
    @test atom_growth_summary.route_core_pair_inventory_status ==
          :available_route_core_unit_pair_inventory
    @test atom_growth_summary.route_core_pair_inventory isa
          GaussletBases.CartesianRouteCore.UnitPairInventory
    @test atom_growth_summary.route_core_pair_count == 36
    @test atom_growth_summary.route_core_pair_order_comparison_source ==
          :atom_growth_pair_inventory
    @test atom_growth_summary.route_core_pair_family_count_source ==
          :crc_final_unit_lowering_recipes
    @test atom_growth_pairs.route_core_pair_inventory_available
    @test atom_growth_pairs.route_core_pair_inventory ===
          atom_growth_summary.route_core_pair_inventory
    pair_inventory = atom_growth_summary.pair_inventory
    unit_inventory = atom_growth_stages.units.plan_unit_inventory
    @test pair_inventory.object_kind == :cartesian_atom_growth_plan_pair_inventory
    @test pair_inventory.status == :available_atom_growth_pair_inventory
    @test pair_inventory.pair_inventory_source == :atom_growth_unit_inventory
    @test pair_inventory.unit_count == unit_inventory.unit_count
    @test pair_inventory.unit_count == 8
    @test pair_inventory.pair_count ==
          pair_inventory.unit_count * (pair_inventory.unit_count + 1) ÷ 2
    @test pair_inventory.pair_count == 36
    @test pair_inventory.pair_entries === atom_growth_pairs.pair_entries
    @test atom_growth_summary.pair_entries === atom_growth_pairs.pair_entries
    @test atom_growth_summary.pair_count == pair_inventory.pair_count
    @test atom_growth_summary.pair_family_counts ==
          pair_inventory.pair_family_counts
    @test pair_inventory.pair_family_counts.white_lindsey_low_order_atom_growth_unit_pair ==
          pair_inventory.pair_count
    @test all(
        pair -> pair.pair_family == :white_lindsey_low_order_atom_growth_unit_pair,
        pair_inventory.pair_entries,
    )
    @test all(
        pair -> pair.pair_contract == :planned_low_order_unit_pair_operator_block,
        pair_inventory.pair_entries,
    )
    @test all(
        pair -> pair.pair_inventory_source == :atom_growth_unit_inventory,
        pair_inventory.pair_entries,
    )
    @test all(
        pair -> pair.pair_planning_source == :final_retained_units,
        pair_inventory.pair_entries,
    )
    @test all(
        pair -> pair.left_final_retained_unit.object_kind ==
                :cartesian_final_retained_unit_contract,
        pair_inventory.pair_entries,
    )
    @test all(
        pair -> pair.right_final_retained_unit.object_kind ==
                :cartesian_final_retained_unit_contract,
        pair_inventory.pair_entries,
    )
    @test all(
        pair -> pair.left_final_retained_unit.downstream_of_lowering,
        pair_inventory.pair_entries,
    )
    @test all(
        pair -> pair.right_final_retained_unit.downstream_of_lowering,
        pair_inventory.pair_entries,
    )
    @test all(
        pair -> !pair.left_final_retained_unit.direct_shellification_region_alias,
        pair_inventory.pair_entries,
    )
    @test all(
        pair -> !pair.right_final_retained_unit.direct_shellification_region_alias,
        pair_inventory.pair_entries,
    )
    @test all(
        pair -> pair.left_owned_support.object_kind ==
                :cartesian_owned_support_region3d,
        pair_inventory.pair_entries,
    )
    @test all(
        pair -> pair.right_owned_support.object_kind ==
                :cartesian_owned_support_region3d,
        pair_inventory.pair_entries,
    )
    @test all(
        pair -> pair.left_unit_index <= pair.right_unit_index,
        pair_inventory.pair_entries,
    )
    expected_pair_keys = Tuple(
        (
            unit_inventory.plan_units[left_index].unit_key,
            unit_inventory.plan_units[right_index].unit_key,
        )
        for left_index in 1:unit_inventory.unit_count
        for right_index in left_index:unit_inventory.unit_count
    )
    @test Tuple(pair.pair_key for pair in pair_inventory.pair_entries) ==
          expected_pair_keys
    @test atom_growth_summary.route_core_pair_keys == expected_pair_keys
    @test atom_growth_summary.route_core_pair_order_matches_staged
    @test atom_growth_pairs.route_core_pair_keys == expected_pair_keys
    route_core_family_counts =
        atom_growth_summary.route_core_pair_family_counts
    @test !isempty(route_core_family_counts)
    @test all(
        record -> record.pair_family_source == :crc_final_unit_lowering_recipes,
        route_core_family_counts,
    )
    @test all(record -> length(record.pair_family) == 2, route_core_family_counts)
    @test sum(record.pair_count for record in route_core_family_counts) ==
          atom_growth_summary.route_core_pair_count
    @test all(
        pair -> !pair.operator_pair_block_materialized,
        pair_inventory.pair_entries,
    )
    @test all(pair -> !pair.operator_block_materialized, pair_inventory.pair_entries)
    @test !pair_inventory.operator_pairs_materialized
    @test !pair_inventory.pair_operator_blocks_materialized
    @test atom_growth_summary.route_skeleton_pair_entry_count ==
          length(atom_growth_stages.units.route_skeleton.pair_entries)
    @test atom_growth_summary.route_skeleton_pair_family_counts ==
          atom_growth_stages.units.route_skeleton.pair_family_counts
    @test atom_growth_summary.route_skeleton_pair_entries ===
          atom_growth_stages.units.route_skeleton.pair_entries
    @test atom_growth_summary.route_skeleton_pair_inventory_source ==
          :route_skeleton_compatibility_fields
    @test !atom_growth_summary.summary_only
    @test atom_growth_summary.pair_stage_fields_preserved
    @test atom_growth_pairs.pair_entries !==
          atom_growth_stages.units.route_skeleton.pair_entries
    @test atom_growth_pairs.route_skeleton_pair_entries ===
          atom_growth_stages.units.route_skeleton.pair_entries
    @test atom_growth_pairs.route_skeleton_pair_family_counts ===
          atom_growth_stages.units.route_skeleton.pair_family_counts
    expected_atom_growth_helper_status = (
        white_lindsey_low_order_atom_growth_unit_pair =
            :deferred_no_pair_operator_block_helper,
    )
    @test atom_growth_pairs.helper_by_pair_family ==
          expected_atom_growth_helper_status
    @test atom_growth_pairs.pair_operator_helper_by_family ==
          expected_atom_growth_helper_status
    @test atom_growth_pairs.pair_helper_status_by_family ==
          expected_atom_growth_helper_status
    @test !(:white_lindsey_low_order in keys(atom_growth_pairs.helper_by_pair_family))
    @test atom_growth_pairs.route_skeleton_helper_by_pair_family ===
          atom_growth_stages.units.route_skeleton.helper_by_pair_family
    @test atom_growth_pairs.pair_stage == :pair_operator_terms_described

    blocked_units = merge(
        atom_growth_stages.units,
        (; plan_unit_inventory = nothing),
    )
    blocked_summary =
        GaussletBases._pqs_source_box_route_driver_pair_stage_low_order_summary(
            blocked_units,
            atom_growth_stages.transforms,
            atom_growth_stages.units.route_skeleton,
        )
    @test blocked_summary.atom_growth_pairs_selected
    @test !blocked_summary.independent_atom_growth_pair_inventory_available
    @test !blocked_summary.pair_inventory_known
    @test blocked_summary.pair_inventory_source ==
          :blocked_missing_plan_unit_inventory
    @test blocked_summary.pair_entries == ()
    @test blocked_summary.pair_count == 0
    @test blocked_summary.pair_family_counts.white_lindsey_low_order_atom_growth_unit_pair ==
          0
    @test blocked_summary.route_core_pair_inventory_available
    @test blocked_summary.route_core_pair_inventory_status ==
          :available_route_core_unit_pair_inventory
    @test blocked_summary.route_core_pair_count == 36
    @test blocked_summary.route_core_pair_order_matches_staged
    @test blocked_summary.route_core_pair_order_comparison_source ==
          :route_core_sidecar_staged_pair_keys
    @test blocked_summary.helper_by_pair_family ==
          expected_atom_growth_helper_status
    @test blocked_summary.route_skeleton_pair_entries ===
          atom_growth_stages.units.route_skeleton.pair_entries
    @test blocked_summary.route_skeleton_pair_family_counts ===
          atom_growth_stages.units.route_skeleton.pair_family_counts
    @test blocked_summary.route_skeleton_helper_by_pair_family ===
          atom_growth_stages.units.route_skeleton.helper_by_pair_family
end

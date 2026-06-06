using Test
using GaussletBases

function _cartesian_pair_stage_low_order_policy_fixture(;
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

function _cartesian_pair_stage_count_by_field(entries, field, value)
    for entry in entries
        getproperty(entry, field) == value && return entry.pair_count
    end
    return 0
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
    @test !default_summary.terminal_shellification_pairs_selected
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
    @test !default_summary.route_core_pair_operator_ready
    @test default_summary.route_core_pair_operator_readiness_status ==
          :not_selected_legacy_source_pairs
    @test default_summary.route_core_pair_operator_blocker ==
          :not_selected_legacy_source_pairs
    @test default_summary.route_core_pair_operator_preflight_available
    @test default_summary.route_core_pair_operator_preflight_status ==
          :blocked_route_core_pair_operator_preflight
    @test default_summary.route_core_pair_operator_preflight_blocker ==
          :not_selected_legacy_source_pairs
    @test default_summary.route_core_pair_operator_preflight.operator_blocks_materialized ==
          false
    @test default_summary.route_core_pair_operator_plan_available
    @test default_summary.route_core_pair_operator_plan_status ==
          :blocked_route_core_pair_operator_plan
    @test default_summary.route_core_pair_operator_plan_blocker ==
          :not_selected_legacy_source_pairs
    @test !default_summary.route_core_pair_operator_plan.operator_blocks_materialized
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
    @test !default_pairs.route_core_pair_operator_ready
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
    @test !atom_growth_summary.terminal_shellification_pairs_selected
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
    @test atom_growth_summary.route_core_pair_operator_ready
    @test atom_growth_summary.route_core_pair_operator_readiness_status ==
          :ready_route_core_pair_operator_metadata
    @test isnothing(atom_growth_summary.route_core_pair_operator_blocker)
    @test atom_growth_summary.route_core_pair_operator_readiness_requirements ==
          (
              :complete_crc_final_unit_inventory,
              :available_crc_pair_inventory,
              :positive_crc_pair_count,
              :crc_pair_order_matches_staged,
              :crc_pair_family_metadata_available,
          )
    @test atom_growth_summary.route_core_pair_operator_preflight_available
    @test atom_growth_summary.route_core_pair_operator_preflight_status ==
          :ready_route_core_pair_operator_preflight
    @test isnothing(atom_growth_summary.route_core_pair_operator_preflight_blocker)
    atom_growth_preflight =
        atom_growth_summary.route_core_pair_operator_preflight
    @test atom_growth_preflight.route_core_final_unit_count == 8
    @test atom_growth_preflight.route_core_pair_count == 36
    @test atom_growth_preflight.route_core_pair_order_matches_staged
    @test !atom_growth_preflight.operator_blocks_materialized
    @test !atom_growth_preflight.pair_operator_blocks_materialized
    @test !atom_growth_preflight.source_operator_blocks_materialized
    @test atom_growth_summary.route_core_pair_operator_plan_available
    @test atom_growth_summary.route_core_pair_operator_plan_status ==
          :ready_route_core_pair_operator_plan
    @test isnothing(atom_growth_summary.route_core_pair_operator_plan_blocker)
    atom_growth_plan = atom_growth_summary.route_core_pair_operator_plan
    @test atom_growth_plan.route_core_final_unit_count == 8
    @test atom_growth_plan.route_core_pair_count == 36
    @test atom_growth_plan.route_core_pair_order_matches_staged
    @test !atom_growth_plan.source_operator_blocks_materialized
    @test !atom_growth_plan.pair_operator_blocks_materialized
    @test !atom_growth_plan.operator_blocks_materialized
    @test !atom_growth_plan.hamiltonian_matrices_materialized
    @test !isempty(atom_growth_plan.operator_block_family_plan)
    @test atom_growth_summary.route_core_typed_pair_operator_plan_inventory_available
    @test atom_growth_summary.route_core_typed_pair_operator_plan_inventory_status ==
          :blocked_route_core_pair_operator_plan_inventory
    @test atom_growth_summary.route_core_typed_pair_operator_plan_blocker ==
          :aggregate_subtree_operator_plan_required
    @test atom_growth_summary.route_core_typed_pair_operator_plan_count == 36
    @test atom_growth_summary.route_core_typed_pair_operator_plan_blocked_count == 15
    @test !atom_growth_summary.route_core_typed_pair_operator_plan_materialized
    @test _cartesian_pair_stage_count_by_field(
        atom_growth_summary.route_core_typed_pair_operator_source_path_counts,
        :source_operator_path,
        :aggregate_subtree_adapter_required,
    ) == 15
    @test _cartesian_pair_stage_count_by_field(
        atom_growth_summary.route_core_typed_pair_operator_materialization_status_counts,
        :materialization_status,
        :metadata_only_not_materialized,
    ) == 21
    @test _cartesian_pair_stage_count_by_field(
        atom_growth_summary.route_core_typed_pair_operator_materialization_status_counts,
        :materialization_status,
        :blocked_metadata_only_not_materialized,
    ) == 15
    @test _cartesian_pair_stage_count_by_field(
        atom_growth_summary.route_core_typed_pair_operator_blocker_counts,
        :blocker,
        nothing,
    ) == 21
    @test _cartesian_pair_stage_count_by_field(
        atom_growth_summary.route_core_typed_pair_operator_blocker_counts,
        :blocker,
        :aggregate_subtree_operator_plan_required,
    ) == 15
    @test sum(
        entry.pair_count for entry in
        atom_growth_summary.route_core_typed_pair_operator_plan_family_counts
    ) == 36
    @test all(
        entry -> !entry.materialized,
        atom_growth_summary.route_core_typed_pair_operator_plan_family_counts,
    )
    @test !atom_growth_summary.route_core_typed_pair_operator_materialization_ready
    @test atom_growth_summary.route_core_typed_pair_operator_materialization_readiness_status ==
          :blocked_pair_operator_materialization
    @test atom_growth_summary.route_core_typed_pair_operator_materialization_readiness_blocker ==
          :aggregate_subtree_operator_plan_required
    @test atom_growth_summary.route_core_typed_pair_operator_materialization_readiness_plan_count ==
          36
    @test atom_growth_summary.route_core_typed_pair_operator_materialization_readiness_blocked_count ==
          15
    @test atom_growth_summary.route_core_typed_pair_operator_materialization_readiness_materialized_count ==
          0
    @test atom_growth_pairs.route_core_pair_inventory_available
    @test atom_growth_pairs.route_core_pair_inventory ===
          atom_growth_summary.route_core_pair_inventory
    @test atom_growth_pairs.route_core_pair_operator_ready
    @test atom_growth_pairs.route_core_pair_operator_preflight ===
          atom_growth_preflight
    @test atom_growth_pairs.route_core_pair_operator_plan === atom_growth_plan
    @test atom_growth_pairs.route_core_typed_pair_operator_plan_count ==
          atom_growth_summary.route_core_typed_pair_operator_plan_count
    @test atom_growth_pairs.route_core_typed_pair_operator_plan_blocked_count ==
          atom_growth_summary.route_core_typed_pair_operator_plan_blocked_count
    @test atom_growth_pairs.route_core_typed_pair_operator_plan_materialized ==
          atom_growth_summary.route_core_typed_pair_operator_plan_materialized
    @test atom_growth_pairs.route_core_typed_pair_operator_materialization_ready ==
          atom_growth_summary.route_core_typed_pair_operator_materialization_ready
    @test atom_growth_pairs.route_core_typed_pair_operator_materialization_readiness_status ==
          atom_growth_summary.route_core_typed_pair_operator_materialization_readiness_status
    @test atom_growth_pairs.route_core_typed_pair_operator_materialization_readiness_blocker ==
          atom_growth_summary.route_core_typed_pair_operator_materialization_readiness_blocker
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

    terminal_fixture =
        _cartesian_pair_stage_low_order_policy_fixture(
            probe_parent_axis_construction = false,
            atom_locations = ((-4.0, 0.0, 0.0), (4.0, 0.0, 0.0)),
            parent_axis_counts = (x = 13, y = 7, z = 7),
        )
    terminal_stages = _cartesian_pair_stage_low_order_policy_pairs(
        terminal_fixture;
        policy = :terminal_cartesian_shellification_geometry,
    )
    terminal_pairs = terminal_stages.pairs
    terminal_summary = terminal_pairs.low_order_pairs
    @test terminal_pairs.low_order_pair_route_kind ==
          :terminal_shellification_low_order_pairs
    @test terminal_pairs.terminal_shellification_pairs_selected
    @test terminal_pairs.terminal_shellification_pair_summary_available
    @test terminal_pairs.terminal_shellification_scaffold_available
    @test terminal_pairs.terminal_shellification_scaffold ===
          terminal_stages.transforms.terminal_shellification_scaffold
    @test terminal_pairs.terminal_shellification_region_count ==
          terminal_stages.transforms.terminal_shellification_region_count
    @test terminal_pairs.terminal_shellification_unit_inventory_available
    @test terminal_pairs.terminal_shellification_unit_inventory ===
          terminal_stages.transforms.terminal_shellification_unit_inventory
    @test terminal_pairs.terminal_shellification_unit_count ==
          terminal_stages.transforms.terminal_shellification_unit_count
    @test terminal_pairs.terminal_shellification_unit_keys ==
          terminal_stages.transforms.terminal_shellification_unit_keys
    @test terminal_pairs.terminal_shellification_unit_roles ==
          terminal_stages.transforms.terminal_shellification_unit_roles
    @test terminal_pairs.terminal_shellification_unit_kinds ==
          terminal_stages.transforms.terminal_shellification_unit_kinds
    @test terminal_pairs.terminal_shellification_unit_support_counts ==
          terminal_stages.transforms.terminal_shellification_unit_support_counts
    @test !terminal_pairs.terminal_shellification_final_retained_unit_inventory_available
    @test !terminal_pairs.terminal_shellification_transform_contracts_available
    @test !terminal_pairs.terminal_shellification_pair_inventory_available
    @test terminal_pairs.terminal_shellification_pair_inventory_status ==
          :deferred_terminal_shellification_pair_inventory
    @test terminal_pairs.terminal_shellification_pair_materialization_status ==
          :deferred_terminal_shellification_pair_materialization
    @test !terminal_pairs.pair_operator_blocks_materialized
    @test !terminal_pairs.operator_pairs_materialized
    @test terminal_pairs.pair_inventory === nothing
    @test terminal_pairs.pair_inventory_source ==
          :terminal_shellification_scaffold
    @test terminal_pairs.pair_entries == ()
    @test terminal_pairs.pair_entries !==
          terminal_stages.units.route_skeleton.pair_entries
    @test terminal_pairs.pair_family_counts == ()
    @test terminal_pairs.pair_family_counts !==
          terminal_stages.units.route_skeleton.pair_family_counts
    @test terminal_pairs.helper_by_pair_family == ()
    @test terminal_pairs.pair_operator_helper_by_family == ()
    @test terminal_pairs.pair_helper_status_by_family == ()
    @test !terminal_pairs.route_core_pair_inventory_available
    @test terminal_pairs.route_core_pair_inventory_status ==
          :deferred_terminal_shellification_pair_inventory
    @test !terminal_pairs.route_core_pair_operator_ready
    @test terminal_pairs.route_core_pair_operator_readiness_status ==
          :deferred_terminal_shellification_pair_inventory
    @test terminal_pairs.route_core_pair_operator_blocker ==
          :deferred_terminal_shellification_pair_inventory
    @test terminal_pairs.route_core_typed_pair_operator_plan_inventory_status ==
          :deferred_terminal_shellification_typed_pair_operator_plan_inventory
    @test terminal_pairs.route_core_typed_pair_operator_plan_blocker ==
          :deferred_terminal_shellification_pair_inventory
    @test !terminal_pairs.route_core_typed_pair_operator_plan_inventory_available
    @test terminal_pairs.route_core_typed_pair_operator_plan_count == 0
    @test terminal_pairs.route_core_typed_pair_operator_plan_blocked_count == 0
    @test !terminal_pairs.route_core_typed_pair_operator_plan_materialized
    @test !terminal_pairs.route_core_typed_pair_operator_materialization_ready
    @test terminal_pairs.terminal_shellification_central_gap_region_count ==
          terminal_stages.transforms.terminal_shellification_central_gap_region_count
    @test terminal_pairs.terminal_shellification_central_midpoint_slab_count ==
          terminal_stages.transforms.terminal_shellification_central_midpoint_slab_count
    @test terminal_pairs.terminal_shellification_central_distorted_product_box_count ==
          terminal_stages.transforms.terminal_shellification_central_distorted_product_box_count
    @test terminal_pairs.terminal_shellification_central_gap_region_count == 3
    @test terminal_pairs.terminal_shellification_central_midpoint_slab_count == 3
    @test terminal_pairs.terminal_shellification_central_distorted_product_box_count ==
          0
    @test terminal_summary.object_kind == :cartesian_pair_stage_low_order_summary
    @test terminal_summary.status ==
          :deferred_terminal_shellification_pair_inventory
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
    @test terminal_summary.pair_route_kind ==
          :terminal_shellification_low_order_pairs
    @test terminal_summary.terminal_shellification_pairs_selected
    @test !terminal_summary.atom_growth_pairs_selected
    @test !terminal_summary.legacy_source_pairs_selected
    @test terminal_summary.terminal_shellification_pair_summary_available
    @test terminal_summary.terminal_shellification_scaffold_available
    @test terminal_summary.terminal_shellification_scaffold ===
          terminal_stages.transforms.terminal_shellification_scaffold
    @test terminal_summary.terminal_shellification_region_count ==
          terminal_stages.transforms.terminal_shellification_region_count
    @test terminal_summary.terminal_shellification_unit_inventory_available
    @test terminal_summary.terminal_shellification_unit_inventory ===
          terminal_stages.transforms.terminal_shellification_unit_inventory
    @test terminal_summary.terminal_shellification_unit_count ==
          terminal_stages.transforms.terminal_shellification_unit_count
    @test terminal_summary.terminal_shellification_unit_keys ==
          terminal_stages.transforms.terminal_shellification_unit_keys
    @test terminal_summary.terminal_shellification_unit_roles ==
          terminal_stages.transforms.terminal_shellification_unit_roles
    @test terminal_summary.terminal_shellification_unit_kinds ==
          terminal_stages.transforms.terminal_shellification_unit_kinds
    @test terminal_summary.terminal_shellification_unit_support_counts ==
          terminal_stages.transforms.terminal_shellification_unit_support_counts
    @test !terminal_summary.terminal_shellification_final_retained_unit_inventory_available
    @test !terminal_summary.terminal_shellification_transform_contracts_available
    @test !terminal_summary.terminal_shellification_pair_inventory_available
    @test terminal_summary.terminal_shellification_pair_inventory_status ==
          :deferred_terminal_shellification_pair_inventory
    @test terminal_summary.terminal_shellification_pair_materialization_status ==
          :deferred_terminal_shellification_pair_materialization
    @test !terminal_summary.pair_inventory_known
    @test terminal_summary.pair_inventory_source ==
          :terminal_shellification_scaffold
    @test terminal_summary.pair_inventory === nothing
    @test terminal_summary.pair_entries == ()
    @test terminal_summary.pair_entries !==
          terminal_stages.units.route_skeleton.pair_entries
    @test terminal_summary.pair_count == 0
    @test terminal_summary.pair_family_counts == ()
    @test terminal_summary.pair_family_counts !==
          terminal_stages.units.route_skeleton.pair_family_counts
    @test terminal_summary.helper_by_pair_family == ()
    @test terminal_summary.pair_operator_helper_by_family == ()
    @test terminal_summary.pair_helper_status_by_family == ()
    @test !terminal_summary.pair_operator_blocks_materialized
    @test !terminal_summary.operator_pairs_materialized
    @test terminal_summary.plan_authority
    @test !terminal_summary.active_source_authority
    @test !terminal_summary.legacy_source_authority
    @test terminal_summary.route_core_pair_inventory_status ==
          :deferred_terminal_shellification_pair_inventory
    @test terminal_summary.route_core_pair_inventory === nothing
    @test terminal_summary.route_core_pair_count == 0
    @test terminal_summary.route_core_pair_keys == ()
    @test terminal_summary.route_core_pair_family_counts == ()
    @test terminal_summary.route_core_summary_status ==
          :deferred_terminal_shellification_pair_inventory
    @test !terminal_summary.route_core_pair_operator_ready
    @test terminal_summary.route_core_pair_operator_preflight_available
    @test terminal_summary.route_core_pair_operator_preflight_status ==
          :blocked_route_core_pair_operator_preflight
    @test terminal_summary.route_core_pair_operator_preflight_blocker ==
          :deferred_terminal_shellification_pair_inventory
    @test !terminal_summary.route_core_pair_operator_preflight.operator_blocks_materialized
    @test terminal_summary.route_core_pair_operator_plan_available
    @test terminal_summary.route_core_pair_operator_plan_status ==
          :blocked_route_core_pair_operator_plan
    @test terminal_summary.route_core_pair_operator_plan_blocker ==
          :deferred_terminal_shellification_pair_inventory
    @test !terminal_summary.route_core_pair_operator_plan.operator_blocks_materialized
    @test !terminal_summary.route_core_typed_pair_operator_plan_inventory_available
    @test terminal_summary.route_core_typed_pair_operator_plan_inventory_status ==
          :deferred_terminal_shellification_typed_pair_operator_plan_inventory
    @test terminal_summary.route_core_typed_pair_operator_plan_blocker ==
          :deferred_terminal_shellification_pair_inventory
    @test terminal_summary.route_skeleton_pair_entries ===
          terminal_stages.units.route_skeleton.pair_entries
    @test terminal_summary.route_skeleton_pair_family_counts ===
          terminal_stages.units.route_skeleton.pair_family_counts
    @test terminal_summary.summary_only
    @test all(
        record -> !record.owned_support_is_cpb,
        terminal_summary.terminal_shellification_unit_inventory.terminal_region_units,
    )
    @test all(
        record -> !record.shellification_region_is_cpb,
        terminal_summary.terminal_shellification_unit_inventory.terminal_region_units,
    )
    @test terminal_summary.pair_stage_fields_preserved
    @test terminal_pairs.route_skeleton_pair_entries ===
          terminal_stages.units.route_skeleton.pair_entries
    @test terminal_pairs.route_skeleton_pair_family_counts ===
          terminal_stages.units.route_skeleton.pair_family_counts
    @test terminal_pairs.route_skeleton_helper_by_pair_family ===
          terminal_stages.units.route_skeleton.helper_by_pair_family
    @test terminal_pairs.pair_stage == :pair_operator_terms_described

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
    @test blocked_summary.route_core_pair_operator_ready
    @test blocked_summary.helper_by_pair_family ==
          expected_atom_growth_helper_status
    @test blocked_summary.route_skeleton_pair_entries ===
          atom_growth_stages.units.route_skeleton.pair_entries
    @test blocked_summary.route_skeleton_pair_family_counts ===
          atom_growth_stages.units.route_skeleton.pair_family_counts
    @test blocked_summary.route_skeleton_helper_by_pair_family ===
          atom_growth_stages.units.route_skeleton.helper_by_pair_family

    missing_crc_units = merge(
        atom_growth_stages.units,
        (; route_core_sidecar_inventory = nothing),
    )
    missing_crc_summary =
        GaussletBases._pqs_source_box_route_driver_pair_stage_low_order_summary(
            missing_crc_units,
            atom_growth_stages.transforms,
            atom_growth_stages.units.route_skeleton,
        )
    @test missing_crc_summary.independent_atom_growth_pair_inventory_available
    @test !missing_crc_summary.route_core_pair_inventory_available
    @test !missing_crc_summary.route_core_pair_operator_ready
    @test missing_crc_summary.route_core_pair_operator_readiness_status ==
          :blocked_missing_route_core_sidecar_inventory
    @test missing_crc_summary.route_core_pair_operator_blocker ==
          :missing_route_core_sidecar_inventory
    @test missing_crc_summary.route_core_pair_operator_preflight_available
    @test missing_crc_summary.route_core_pair_operator_preflight_status ==
          :blocked_route_core_pair_operator_preflight
    @test missing_crc_summary.route_core_pair_operator_preflight_blocker ==
          :missing_route_core_sidecar_inventory
    @test !missing_crc_summary.route_core_pair_operator_preflight.operator_blocks_materialized
    @test missing_crc_summary.route_core_pair_operator_plan_available
    @test missing_crc_summary.route_core_pair_operator_plan_status ==
          :blocked_route_core_pair_operator_plan
    @test missing_crc_summary.route_core_pair_operator_plan_blocker ==
          :missing_route_core_sidecar_inventory
    @test !missing_crc_summary.route_core_pair_operator_plan.operator_blocks_materialized
end

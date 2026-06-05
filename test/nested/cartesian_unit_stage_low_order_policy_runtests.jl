using Test
using GaussletBases

function _cartesian_unit_stage_low_order_policy_fixture()
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
        route_kind = :be2_cartesian_unit_stage_low_order_policy,
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

@testset "cartesian unit stage carries selected low-order policy" begin
    fixture = _cartesian_unit_stage_low_order_policy_fixture()

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
    default_summary = default_units.low_order_units
    @test default_units.object_kind == :cartesian_units
    @test default_summary.object_kind ==
          :cartesian_unit_stage_low_order_summary
    @test default_summary.low_order_shellization_policy_resolved ==
          :legacy_diatomic_source
    @test default_summary.shellization_source ==
          :route_configured_bond_aligned_diatomic_source
    @test default_summary.shellization_kind == :legacy_diatomic_source
    @test default_summary.unit_route_kind ==
          :legacy_diatomic_source_low_order_units
    @test default_summary.legacy_source_units_selected
    @test default_summary.active_source_authority
    @test !default_summary.atom_growth_units_selected
    @test !default_summary.atom_growth_unit_summary_available
    @test !default_summary.atom_growth_unit_inventory_available
    @test !default_summary.plan_unit_inventory_available
    @test default_summary.unit_inventory_source ==
          :legacy_diatomic_source_summary
    @test !default_summary.materialized_units_available
    @test default_summary.materialization_status ==
          :deferred_legacy_diatomic_source_unit_materialization
    @test !default_summary.retained_unit_dimensions_known
    @test !default_summary.retained_unit_ranges_known
    @test !default_summary.retained_dimension_known
    @test default_summary.route_skeleton_unit_fields_preserved
    @test default_units.source_boxes === default_units.route_skeleton.source_boxes
    @test default_units.source_dimensions ===
          default_units.route_skeleton.source_dimensions
    @test default_units.retained_units === default_units.route_skeleton.retained_units
    @test default_units.retained_unit_order ===
          default_units.route_skeleton.retained_unit_order

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
    atom_growth_summary = atom_growth_units.low_order_units
    @test atom_growth_units.low_order_unit_route_kind ==
          :atom_growth_complete_rectangular_low_order_units
    @test atom_growth_units.atom_growth_unit_summary_available
    @test atom_growth_units.atom_growth_units_selected
    @test atom_growth_units.atom_growth_unit_inventory_available
    @test atom_growth_units.plan_unit_inventory_available
    @test atom_growth_units.unit_inventory_source ==
          :atom_growth_shellification_plan
    @test atom_growth_units.unit_inventory_status ==
          :available_atom_growth_plan_unit_inventory
    @test !atom_growth_units.active_source_authority
    @test atom_growth_summary.low_order_shellization_policy_resolved ==
          :atom_growth_complete_rectangular
    @test atom_growth_summary.shellization_source ==
          :bond_aligned_diatomic_atom_growth_construction_plan
    @test atom_growth_summary.shellization_kind ==
          :atom_growth_complete_rectangular
    @test atom_growth_summary.unit_route_kind ==
          :atom_growth_complete_rectangular_low_order_units
    @test atom_growth_summary.atom_growth_units_selected
    @test !atom_growth_summary.legacy_source_units_selected
    @test atom_growth_summary.atom_growth_unit_summary_available
    @test atom_growth_summary.atom_growth_unit_inventory_available
    @test atom_growth_summary.plan_unit_inventory_available
    @test atom_growth_summary.unit_inventory_source ==
          :atom_growth_shellification_plan
    @test atom_growth_summary.unit_inventory_status ==
          :available_atom_growth_plan_unit_inventory
    @test atom_growth_summary.plan_authority
    @test !atom_growth_summary.active_source_authority
    @test !atom_growth_summary.legacy_source_authority
    @test !atom_growth_summary.materialized_units_available
    @test atom_growth_summary.materialization_status ==
          :deferred_atom_growth_complete_rectangular_unit_materialization
    @test !atom_growth_summary.retained_unit_dimensions_known
    @test !atom_growth_summary.retained_unit_ranges_known
    @test !atom_growth_summary.retained_dimension_known
    @test atom_growth_summary.retained_dimension === nothing
    @test !atom_growth_summary.summary_only
    plan_inventory = atom_growth_summary.plan_unit_inventory
    @test plan_inventory.object_kind ==
          :cartesian_atom_growth_plan_unit_inventory
    @test plan_inventory.unit_inventory_source ==
          :atom_growth_shellification_plan
    @test plan_inventory.unit_roles ==
          atom_growth_shells.low_order_shellization.atom_growth_scaffold.ordered_region_roles
    @test plan_inventory.unit_count == length(plan_inventory.plan_units)
    @test plan_inventory.unit_count ==
          atom_growth_shells.low_order_shellization.atom_growth_scaffold.region_count
    @test plan_inventory.source_backed_region_count == 0
    @test plan_inventory.source_backed_unit_count == 0
    @test plan_inventory.plan_lowerable_unit_count == plan_inventory.unit_count
    @test !plan_inventory.materialized_units_available
    @test !plan_inventory.retained_unit_dimensions_known
    @test !plan_inventory.retained_unit_ranges_known
    @test !plan_inventory.retained_dimension_known
    @test plan_inventory.retained_dimension === nothing
    @test plan_inventory.cpb_contract_stage ==
          :shellification_to_cpb_lowering_to_construction
    @test !plan_inventory.shellification_regions_are_cpbs
    @test plan_inventory.owned_support_available
    @test plan_inventory.lowering_source_cpbs_available
    @test plan_inventory.source_cpb_count == 3
    @test all(unit -> !unit.source_backed, plan_inventory.plan_units)
    @test all(unit -> unit.independently_lowerable, plan_inventory.plan_units)
    @test all(
        unit -> unit.retirement_target == :already_plan_lowered_region,
        plan_inventory.plan_units,
    )
    @test all(unit -> unit.support_count > 0, plan_inventory.plan_units)
    @test all(unit -> unit.retained_count === nothing, plan_inventory.plan_units)
    @test all(unit -> unit.retained_range === nothing, plan_inventory.plan_units)
    @test all(unit -> !unit.shellification_region_is_cpb, plan_inventory.plan_units)
    @test all(
        unit -> unit.owned_support.object_kind == :cartesian_owned_support_region3d,
        plan_inventory.plan_units,
    )
    @test all(
        unit -> !unit.owned_support.shellification_region_is_cpb,
        plan_inventory.plan_units,
    )
    @test all(
        unit -> unit.lowering_recipe.object_kind == :cartesian_cpb_lowering_recipe,
        plan_inventory.plan_units,
    )
    @test all(
        unit -> unit.lowering_recipe.lowering_stage ==
                :coordinate_product_box_lowering,
        plan_inventory.plan_units,
    )
    @test all(
        unit -> unit.intermediate_retained_space.object_kind ==
                :cartesian_intermediate_retained_space_contract,
        plan_inventory.plan_units,
    )
    @test all(
        unit -> unit.shell_realization.object_kind ==
                :cartesian_shell_realization_contract,
        plan_inventory.plan_units,
    )
    @test all(
        unit -> unit.final_retained_unit.object_kind ==
                :cartesian_final_retained_unit_contract,
        plan_inventory.plan_units,
    )
    @test all(
        unit -> unit.final_retained_unit.downstream_of_lowering,
        plan_inventory.plan_units,
    )
    @test all(
        unit -> !unit.final_retained_unit.direct_shellification_region_alias,
        plan_inventory.plan_units,
    )
    @test all(
        unit -> unit.cpb_contract.final_unit_downstream_of_lowering,
        plan_inventory.plan_units,
    )
    shell_difference_units = Tuple(
        unit for unit in plan_inventory.plan_units
        if unit.owned_support.owned_support_is_shell_difference
    )
    @test !isempty(shell_difference_units)
    @test all(
        unit -> unit.owned_support.support_kind == :shell_difference_owned_support,
        shell_difference_units,
    )
    @test all(unit -> !unit.owned_support.box_difference_is_cpb, shell_difference_units)
    @test all(
        unit -> !isnothing(unit.owned_support.inner_exclusion_cpb),
        shell_difference_units,
    )
    @test all(
        unit -> unit.owned_support.inner_exclusion_cpb.object_kind ==
                :cartesian_coordinate_product_box3d,
        shell_difference_units,
    )
    explicit_source_cpbs = reduce(
        vcat,
        (collect(unit.source_cpbs) for unit in plan_inventory.plan_units);
        init = NamedTuple[],
    )
    @test length(explicit_source_cpbs) == plan_inventory.source_cpb_count
    @test all(
        cpb -> cpb.object_kind == :cartesian_coordinate_product_box3d,
        explicit_source_cpbs,
    )
    @test all(cpb -> cpb.coordinate_product_box, explicit_source_cpbs)
    @test count(
        unit -> unit.lowering_family == :outer_mismatch_boundary_slab_set &&
                unit.source_cpb_count == 2 &&
                unit.lowering_recipe.source_cpb_families ==
                (:direct_boundary_slab_cpb,),
        plan_inventory.plan_units,
    ) == 1
    @test count(
        unit -> unit.materialization_dependency == :plan_lowerable_direct_slab &&
                unit.source_cpb_count == 1 &&
                unit.lowering_recipe.source_cpb_families == (:direct_slab_cpb,),
        plan_inventory.plan_units,
    ) == 1
    @test all(
        unit -> unit.lowering_recipe.source_cpb_families ==
                (:facet_cpb, :edge_cpb, :corner_cpb),
        (
            unit for unit in plan_inventory.plan_units
            if unit.lowering_family == :white_lindsey_adaptive_complete_shell
        ),
    )
    @test all(
        unit -> unit.lowering_recipe.source_cpb_enumeration_status ==
                :planned_cpb_families_not_enumerated,
        (
            unit for unit in plan_inventory.plan_units
            if unit.lowering_family == :white_lindsey_adaptive_complete_shell
        ),
    )
    @test all(
        unit -> unit.source_descriptor.final_column_ranges_available == false,
        plan_inventory.plan_units,
    )
    @test Tuple(keys(plan_inventory.support_counts)) == plan_inventory.unit_keys
    @test Tuple(values(plan_inventory.support_counts)) ==
          Tuple(unit.support_count for unit in plan_inventory.plan_units)
    @test atom_growth_summary.route_skeleton_unit_fields_preserved
    @test atom_growth_summary.route_skeleton_unit_inventory_source ==
          :route_skeleton_compatibility_fields
    @test atom_growth_units.source_boxes ===
          atom_growth_units.route_skeleton.source_boxes
    @test atom_growth_units.source_dimensions ===
          atom_growth_units.route_skeleton.source_dimensions
    @test atom_growth_units.retained_units ===
          atom_growth_units.route_skeleton.retained_units
    @test atom_growth_units.retained_unit_order ===
          atom_growth_units.route_skeleton.retained_unit_order
end

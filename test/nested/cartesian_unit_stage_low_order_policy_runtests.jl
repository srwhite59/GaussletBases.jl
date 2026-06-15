using Test
using GaussletBases

function _cartesian_unit_stage_low_order_policy_fixture(;
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

function _cartesian_unit_stage_box_points(box)
    return Set((x, y, z) for x in box[1] for y in box[2] for z in box[3])
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
    @test default_summary.legacy_source_units_selected
    @test default_summary.active_source_authority
    @test !default_summary.atom_growth_units_selected
    @test !default_summary.terminal_shellification_units_selected
    @test default_summary.route_skeleton_unit_fields_preserved
    @test length(default_units.source_boxes) ==
          length(default_units.route_skeleton.source_boxes)
    @test default_units.source_dimensions ==
          default_units.route_skeleton.source_dimensions
    @test length(default_units.retained_units) ==
          length(default_units.route_skeleton.retained_units)
    @test default_units.retained_unit_order ==
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
    @test !atom_growth_summary.terminal_shellification_units_selected
    @test !atom_growth_summary.legacy_source_units_selected
    @test atom_growth_summary.atom_growth_unit_summary_available
    @test !atom_growth_summary.terminal_shellification_unit_summary_available
    @test atom_growth_summary.atom_growth_unit_inventory_available
    @test !atom_growth_summary.terminal_shellification_unit_inventory_available
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
    @test plan_inventory.source_cpb_count == 107
    @test plan_inventory.pqs_lowering_prototype_available
    @test !isnothing(plan_inventory.pqs_lowering_prototype_unit_key)
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
        unit -> unit.owned_support.owned_support_authority ==
                :shellification_region,
        plan_inventory.plan_units,
    )
    @test all(
        unit -> unit.owned_support.shellification_authority_scope ==
                :owned_support_only,
        plan_inventory.plan_units,
    )
    @test all(
        unit -> !unit.owned_support.shellification_region_is_lowering_source,
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
        unit -> unit.lowering_recipe.owned_support_authority ==
                :shellification_region,
        plan_inventory.plan_units,
    )
    @test all(
        unit -> unit.lowering_recipe.shellification_authority_scope ==
                :owned_support_only,
        plan_inventory.plan_units,
    )
    @test all(
        unit -> !unit.lowering_recipe.shellification_region_is_lowering_source,
        plan_inventory.plan_units,
    )
    @test all(
        unit -> unit.lowering_recipe.lowering_source_authority ==
                :lowering_recipe_cpbs,
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
    @test all(unit -> !unit.owned_support.owned_support_is_cpb, shell_difference_units)
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
    complete_shell_units = Tuple(
        unit for unit in plan_inventory.plan_units
        if unit.lowering_family == :white_lindsey_adaptive_complete_shell
    )
    @test length(complete_shell_units) == 4
    @test all(unit -> unit.source_cpb_count == 26, complete_shell_units)
    @test plan_inventory.source_cpb_count ==
          3 + 26 * length(complete_shell_units)
    @test all(
        unit -> unit.lowering_recipe.source_cpb_families ==
                (:facet_cpb, :edge_cpb, :corner_cpb),
        complete_shell_units,
    )
    @test all(
        unit -> unit.lowering_recipe.source_cpb_enumeration_status ==
                :explicit_complete_shell_boundary_strata,
        complete_shell_units,
    )
    @test all(
        unit -> unit.lowering_recipe.source_cpb_enumeration_reason === nothing,
        complete_shell_units,
    )
    @test all(
        unit -> unit.lowering_recipe.enumeration_policy ==
                :white_lindsey_complete_shell_boundary_strata,
        complete_shell_units,
    )
    @test all(
        unit -> unit.lowering_recipe.complete_shell_condition_status ==
                :supported_complete_shell,
        complete_shell_units,
    )
    @test all(
        unit -> unit.lowering_recipe.complete_shell_cpb_family_counts ==
                (facet_cpb = 6, edge_cpb = 12, corner_cpb = 8),
        complete_shell_units,
    )
    @test all(
        unit -> !unit.lowering_recipe.coefficient_maps_materialized,
        complete_shell_units,
    )
    @test all(
        unit -> !unit.lowering_recipe.operator_blocks_materialized,
        complete_shell_units,
    )
    @test all(
        unit -> !unit.lowering_recipe.pair_operator_blocks_materialized,
        complete_shell_units,
    )
    @test all(
        unit -> !unit.lowering_recipe.hamiltonian_data_materialized,
        complete_shell_units,
    )
    complete_shell_cpbs = reduce(
        vcat,
        (collect(unit.source_cpbs) for unit in complete_shell_units);
        init = NamedTuple[],
    )
    complete_shell_roles = Set(cpb.role for cpb in complete_shell_cpbs)
    @test :x_low_facet in complete_shell_roles
    @test :y_high_z_low_edge in complete_shell_roles
    @test :x_low_y_high_z_high_corner in complete_shell_roles
    @test all(
        cpb -> cpb.metadata.enumeration_policy ==
               :white_lindsey_complete_shell_boundary_strata,
        complete_shell_cpbs,
    )
    @test all(
        cpb -> cpb.metadata.shellification_authority_scope ==
               :owned_support_only,
        complete_shell_cpbs,
    )
    @test all(cpb -> cpb.metadata.lowering_geometry_only, complete_shell_cpbs)
    @test all(
        cpb -> !cpb.metadata.coefficient_maps_materialized,
        complete_shell_cpbs,
    )
    @test all(
        cpb -> !cpb.metadata.operator_blocks_materialized,
        complete_shell_cpbs,
    )
    @test all(
        cpb -> !cpb.metadata.pair_operator_blocks_materialized,
        complete_shell_cpbs,
    )
    @test all(
        cpb -> !cpb.metadata.hamiltonian_data_materialized,
        complete_shell_cpbs,
    )
    for unit in complete_shell_units
        @test count(cpb -> cpb.cpb_family == :facet_cpb, unit.source_cpbs) == 6
        @test count(cpb -> cpb.cpb_family == :edge_cpb, unit.source_cpbs) == 12
        @test count(cpb -> cpb.cpb_family == :corner_cpb, unit.source_cpbs) == 8
        @test sum(cpb.support_count for cpb in unit.source_cpbs) ==
              unit.support_count
        stratum_points = Set{NTuple{3,Int}}()
        for cpb in unit.source_cpbs
            cpb_points = _cartesian_unit_stage_box_points(cpb.box)
            @test isempty(intersect(stratum_points, cpb_points))
            union!(stratum_points, cpb_points)
        end
        owned_shell_points = setdiff(
            _cartesian_unit_stage_box_points(unit.owned_support.outer_cpb.box),
            _cartesian_unit_stage_box_points(
                unit.owned_support.inner_exclusion_cpb.box,
            ),
        )
        @test stratum_points == owned_shell_points
        @test length(stratum_points) == unit.support_count
    end
    pqs_prototype = plan_inventory.pqs_lowering_prototype
    @test pqs_prototype.object_kind ==
          :cartesian_pqs_lowering_metadata_prototype
    @test pqs_prototype.status == :metadata_only_planned
    @test pqs_prototype.unit_key == plan_inventory.pqs_lowering_prototype_unit_key
    pqs_unit = only(
        unit for unit in plan_inventory.plan_units
        if unit.unit_key == pqs_prototype.unit_key
    )
    @test pqs_unit.unit_role == :regular_shared_molecular_shell
    @test pqs_prototype.owned_support_authority == :shellification_region
    @test pqs_prototype.shellification_authority_scope == :owned_support_only
    @test !pqs_prototype.shellification_region_is_lowering_source
    @test !pqs_prototype.owned_support_is_cpb
    @test pqs_prototype.owned_support_is_shell_difference
    @test pqs_prototype.source_cpb.object_kind ==
          :cartesian_coordinate_product_box3d
    @test pqs_prototype.source_cpb.cpb_family == :filled_source_cpb
    @test pqs_prototype.source_cpb.coordinate_product_box
    @test pqs_prototype.source_cpb.codimension == 0
    @test pqs_prototype.source_cpb.support_count ==
          prod(pqs_prototype.source_cpb.dimensions)
    @test pqs_prototype.source_cpb.support_count ==
          pqs_prototype.source_cpb_support_count
    @test pqs_prototype.source_cpb_support_count_source ==
          :filled_coordinate_product_box
    @test pqs_prototype.owned_support.support_count == pqs_unit.support_count
    @test pqs_prototype.owned_support_count == pqs_unit.support_count
    @test pqs_prototype.owned_support_count_source == :shellification_region
    @test pqs_prototype.source_cpb_support_count !=
          pqs_prototype.owned_support_count
    @test pqs_prototype.lowering_source_authority ==
          :pqs_lowering_recipe_filled_source_cpb
    @test pqs_prototype.lowering_recipe == :pqs_filled_source_cpb
    @test pqs_prototype.lowering_recipe_contract.retained_rule ==
          :boundary_comx_product_mode_selection
    @test pqs_prototype.intermediate_retained_space.object_kind ==
          :pqs_boundary_comx_product_intermediate_retained_space
    @test pqs_prototype.intermediate_retained_space.status == :planned_deferred
    @test pqs_prototype.intermediate_retained_space.selected_modes ==
          :boundary_comx_product_modes
    @test pqs_prototype.shell_realization.object_kind ==
          :pqs_shell_projection_lowdin_realization
    @test pqs_prototype.shell_realization.status == :planned_deferred
    @test pqs_prototype.shell_realization.shell_projection_planned
    @test pqs_prototype.shell_realization.lowdin_cleanup_planned
    @test !pqs_prototype.dense_parent_space_operators_are_algorithm
    @test !pqs_prototype.shell_row_operator_algorithm
    @test !pqs_prototype.coefficient_maps_materialized
    @test !pqs_prototype.operator_blocks_materialized
    @test !pqs_prototype.pair_operator_blocks_materialized
    @test !pqs_prototype.public_route_adoption
    @test pqs_unit.lowering_recipe.source_cpb_families ==
          (:facet_cpb, :edge_cpb, :corner_cpb)
    @test pqs_unit.source_cpb_count == 26
    @test pqs_unit.lowering_recipe.source_cpb_enumeration_status ==
          :explicit_complete_shell_boundary_strata
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
    @test atom_growth_units.pqs_lowering_prototype_available
    @test atom_growth_units.pqs_lowering_prototype.status ==
          plan_inventory.pqs_lowering_prototype.status
    @test atom_growth_units.pqs_lowering_prototype.unit_key ==
          plan_inventory.pqs_lowering_prototype.unit_key
    @test length(atom_growth_units.source_boxes) ==
          length(atom_growth_units.route_skeleton.source_boxes)
    @test atom_growth_units.source_dimensions ==
          atom_growth_units.route_skeleton.source_dimensions
    @test length(atom_growth_units.retained_units) ==
          length(atom_growth_units.route_skeleton.retained_units)
    @test atom_growth_units.retained_unit_order ==
          atom_growth_units.route_skeleton.retained_unit_order

    terminal_fixture =
        _cartesian_unit_stage_low_order_policy_fixture(
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
    terminal_summary = terminal_units.low_order_units
    @test terminal_units.low_order_unit_route_kind ==
          :terminal_shellification_low_order_units
    @test terminal_units.terminal_shellification_units_selected
    @test terminal_units.terminal_shellification_unit_summary_available
    @test terminal_units.terminal_shellification_scaffold_available
    @test terminal_units.terminal_shellification_region_count ==
          terminal_shells.terminal_shellification_region_count
    @test terminal_units.terminal_shellification_unit_inventory_available
    @test terminal_units.terminal_shellification_unit_inventory_status ==
          :available_terminal_region_unit_inventory
    terminal_inventory = terminal_units.terminal_shellification_unit_inventory
    @test terminal_inventory.object_kind ==
          :cartesian_terminal_region_unit_inventory
    @test terminal_units.terminal_shellification_unit_count ==
          terminal_inventory.unit_count
    @test terminal_units.terminal_shellification_unit_count ==
          terminal_units.terminal_shellification_region_count
    @test terminal_units.terminal_shellification_lowering_contract_inventory_available
    terminal_lowering_inventory =
        terminal_units.terminal_shellification_lowering_contract_inventory
    @test terminal_units.terminal_shellification_lowering_contract_count >=
          terminal_units.terminal_shellification_unit_count
    @test all(
        entry.lowering_contract_count >= 1
        for entry in terminal_units.terminal_shellification_contract_counts_by_unit
    )
    @test terminal_units.terminal_shellification_selected_lowering_contract_inventory_available
    selected_terminal_lowering_inventory =
        terminal_units.terminal_shellification_selected_lowering_contract_inventory
    @test selected_terminal_lowering_inventory.route_lowering_family ==
          :white_lindsey_low_order
    @test all(
        entry.selected_contract_count == 1
        for entry in terminal_units.terminal_shellification_selected_contract_counts_by_unit
    )
    @test !selected_terminal_lowering_inventory.final_retained_unit_inventory_available
    @test !selected_terminal_lowering_inventory.pair_inventory_available
    @test selected_terminal_lowering_inventory.pair_inventory_status ==
          :not_available_selected_lowering_metadata_only
    @test !selected_terminal_lowering_inventory.coefficient_maps_materialized
    @test !selected_terminal_lowering_inventory.transform_contracts_materialized
    @test !selected_terminal_lowering_inventory.retained_spaces_materialized
    @test !selected_terminal_lowering_inventory.operator_blocks_materialized
    @test !selected_terminal_lowering_inventory.pair_operator_blocks_materialized
    @test !selected_terminal_lowering_inventory.hamiltonian_data_materialized
    @test !selected_terminal_lowering_inventory.artifacts_materialized
    selected_crc_sidecars =
        terminal_summary.terminal_shellification_selected_crc_sidecar_summary
    @test selected_crc_sidecars.object_kind ==
          :cartesian_unit_stage_selected_terminal_lowering_crc_sidecar_summary
    @test selected_crc_sidecars.status in (
        :available_selected_terminal_lowering_crc_sidecar_inventory,
        :partial_selected_terminal_lowering_crc_sidecar_inventory,
    )
    @test selected_crc_sidecars.private_development_only
    @test terminal_units.terminal_shellification_lw_complete_shell_cpb_count ==
          terminal_lowering_inventory.lw_complete_shell_cpb_count
    @test terminal_units.terminal_shellification_lw_complete_shell_cpb_family_counts ==
          terminal_lowering_inventory.lw_complete_shell_cpb_family_counts
    @test !terminal_lowering_inventory.final_retained_unit_inventory_available
    @test !terminal_lowering_inventory.pair_inventory_available
    @test !terminal_lowering_inventory.coefficient_maps_materialized
    @test !terminal_lowering_inventory.transform_contracts_materialized
    @test !terminal_lowering_inventory.retained_spaces_materialized
    @test !terminal_lowering_inventory.operator_blocks_materialized
    @test !terminal_lowering_inventory.pair_operator_blocks_materialized
    @test !terminal_lowering_inventory.hamiltonian_data_materialized
    @test !terminal_lowering_inventory.artifacts_materialized
    @test !terminal_units.terminal_shellification_final_retained_unit_inventory_available
    @test !terminal_units.terminal_shellification_pair_inventory_available
    @test terminal_units.terminal_shellification_central_gap_region_count ==
          terminal_shells.terminal_shellification_central_gap_region_count
    @test terminal_units.terminal_shellification_central_midpoint_slab_count ==
          terminal_shells.terminal_shellification_central_midpoint_slab_count
    @test terminal_units.terminal_shellification_central_distorted_product_box_count ==
          terminal_shells.terminal_shellification_central_distorted_product_box_count
    @test terminal_units.terminal_shellification_central_gap_region_count == 3
    @test terminal_units.terminal_shellification_central_midpoint_slab_count == 3
    @test terminal_units.terminal_shellification_central_distorted_product_box_count ==
          0
    @test terminal_summary.object_kind ==
          :cartesian_unit_stage_low_order_summary
    @test terminal_summary.status ==
          :available_unit_stage_low_order_summary
    @test terminal_summary.low_order_shellization_policy_resolved ==
          :terminal_cartesian_shellification_geometry
    @test terminal_summary.shellization_source ==
          :terminal_cartesian_shellification_geometry
    @test terminal_summary.shellization_kind ==
          :terminal_cartesian_shellification_geometry
    @test terminal_summary.unit_route_kind ==
          :terminal_shellification_low_order_units
    @test terminal_summary.terminal_shellification_units_selected
    @test !terminal_summary.atom_growth_units_selected
    @test !terminal_summary.legacy_source_units_selected
    @test terminal_summary.terminal_shellification_unit_summary_available
    @test terminal_summary.terminal_shellification_scaffold_available
    @test terminal_summary.terminal_shellification_region_count ==
          terminal_shells.terminal_shellification_region_count
    @test terminal_summary.terminal_shellification_unit_inventory_available
    @test terminal_summary.terminal_shellification_unit_inventory_status ==
          :available_terminal_region_unit_inventory
    @test terminal_summary.terminal_shellification_unit_count ==
          terminal_inventory.unit_count
    @test terminal_summary.terminal_shellification_lowering_contract_inventory_available
    @test terminal_summary.terminal_shellification_lowering_contract_count ==
          terminal_lowering_inventory.lowering_contract_count
    @test terminal_summary.terminal_shellification_lw_complete_shell_cpb_count ==
          terminal_lowering_inventory.lw_complete_shell_cpb_count
    @test terminal_summary.terminal_shellification_lw_complete_shell_cpb_family_counts ==
          terminal_lowering_inventory.lw_complete_shell_cpb_family_counts
    @test !terminal_summary.terminal_shellification_final_retained_unit_inventory_available
    @test !terminal_summary.terminal_shellification_pair_inventory_available
    @test terminal_summary.plan_unit_inventory_available
    @test terminal_summary.plan_unit_count == terminal_inventory.unit_count
    @test terminal_summary.plan_unit_roles == terminal_inventory.unit_roles
    @test terminal_summary.plan_unit_keys == terminal_inventory.unit_keys
    @test terminal_summary.plan_unit_support_counts ==
          terminal_inventory.support_counts
    @test terminal_summary.unit_inventory_source ==
          :terminal_shellification_scaffold
    @test terminal_summary.unit_inventory_status ==
          :available_terminal_region_unit_inventory
    @test !terminal_summary.atom_growth_unit_inventory_available
    @test terminal_summary.plan_authority
    @test !terminal_summary.active_source_authority
    @test !terminal_summary.legacy_source_authority
    @test !terminal_summary.materialized_units_available
    @test terminal_summary.materialization_status ==
          :deferred_terminal_shellification_unit_materialization
    @test !terminal_summary.retained_unit_dimensions_known
    @test !terminal_summary.retained_unit_ranges_known
    @test !terminal_summary.retained_dimension_known
    @test terminal_summary.retained_dimension === nothing
    @test terminal_summary.summary_only
    @test !terminal_summary.shellification_regions_are_cpbs
    @test !terminal_summary.owned_support_available
    @test !terminal_summary.lowering_source_cpbs_available
    @test terminal_summary.source_cpb_count == 0
    @test !terminal_summary.pqs_lowering_prototype_available
    @test all(
        record.unit_key ∉ (:left_atom_box, :right_atom_box)
        for record in terminal_inventory.terminal_region_units
    )
    @test all(
        !record.owned_support_is_cpb
        for record in terminal_inventory.terminal_region_units
    )
    @test all(
        !record.shellification_region_is_cpb
        for record in terminal_inventory.terminal_region_units
    )
    direct_terminal_units =
        filter(
            record ->
                record.terminal_region_kind in (
                    :direct_core,
                    :direct_midpoint_slab,
                    :outer_mismatch_slab,
                ),
            terminal_inventory.terminal_region_units,
        )
    @test !isempty(direct_terminal_units)
    @test all(
        record.source_cpb_plan.source_cpb_plan_equals_owned_support
        for record in direct_terminal_units
    )
    direct_terminal_contracts =
        filter(
            contract ->
                contract.lowering_contract_kind in (
                    :direct_core_identity_cpb,
                    :direct_slab_identity_cpb,
                    :direct_boundary_slab_identity_cpb,
                ),
            terminal_lowering_inventory.lowering_contracts,
        )
    @test !isempty(direct_terminal_contracts)
    @test all(
        contract.identity_like_source_contract
        for contract in direct_terminal_contracts
    )
    @test all(
        contract.source_cpb_plan_equals_owned_support
        for contract in direct_terminal_contracts
    )
    selected_direct_terminal_contracts =
        filter(
            contract ->
                contract.terminal_region_kind in (
                    :direct_core,
                    :direct_midpoint_slab,
                    :outer_mismatch_slab,
                ),
            selected_terminal_lowering_inventory.selected_contracts,
        )
    @test !isempty(selected_direct_terminal_contracts)
    @test all(
        contract.identity_like_source_contract
        for contract in selected_direct_terminal_contracts
    )
    complete_shell_units =
        filter(
            record -> record.terminal_region_kind == :complete_shell,
            terminal_inventory.terminal_region_units,
        )
    @test all(
        record.lw_complete_shell_lowering.total_source_cpb_count == 26
        for record in complete_shell_units
    )
    @test all(
        record.pqs_complete_shell_lowering.source_cpb_plan_kind ==
        :filled_source_cpb
        for record in complete_shell_units
    )
    if !isempty(complete_shell_units)
        lw_contracts =
            filter(
                contract ->
                    contract.lowering_contract_kind ==
                    :white_lindsey_boundary_strata,
                terminal_lowering_inventory.lowering_contracts,
            )
        @test length(lw_contracts) == length(complete_shell_units)
        @test all(
            contract.source_cpb_family_counts ==
            (facet_cpb = 6, edge_cpb = 12, corner_cpb = 8)
            for contract in lw_contracts
        )
        @test all(contract.source_cpb_count == 26 for contract in lw_contracts)
        pqs_contracts =
            filter(
                contract ->
                    contract.lowering_contract_kind == :pqs_filled_source_cpb,
                terminal_lowering_inventory.lowering_contracts,
            )
        @test length(pqs_contracts) == length(complete_shell_units)
        @test all(contract.source_cpb_count == 1 for contract in pqs_contracts)
        @test all(
            contract.source_cpb_plan_kind == :filled_source_cpb
            for contract in pqs_contracts
        )
        @test all(
            !contract.face_edge_corner_decomposition_required
            for contract in pqs_contracts
        )
        selected_shell_contracts =
            filter(
                contract -> contract.terminal_region_kind == :complete_shell,
                selected_terminal_lowering_inventory.selected_contracts,
            )
        @test length(selected_shell_contracts) == length(complete_shell_units)
        @test all(
            contract.lowering_contract_kind ==
            :white_lindsey_boundary_strata
            for contract in selected_shell_contracts
        )
        @test !any(
            contract ->
                contract.lowering_contract_kind == :pqs_filled_source_cpb,
            selected_terminal_lowering_inventory.selected_contracts,
        )
        @test selected_terminal_lowering_inventory.selected_contract_kind_counts.pqs_filled_source_cpb_count ==
              0
    end
    @test terminal_summary.route_skeleton_unit_fields_preserved
    @test length(terminal_units.source_boxes) ==
          length(terminal_units.route_skeleton.source_boxes)
    @test terminal_units.source_dimensions ==
          terminal_units.route_skeleton.source_dimensions
    @test length(terminal_units.retained_units) ==
          length(terminal_units.route_skeleton.retained_units)
    @test terminal_units.retained_unit_order ==
          terminal_units.route_skeleton.retained_unit_order

    distorted_fixture =
        _cartesian_unit_stage_low_order_policy_fixture(
            probe_parent_axis_construction = false,
            atom_locations = ((-8.0, 0.0, 0.0), (8.0, 0.0, 0.0)),
            parent_axis_counts = (x = 21, y = 7, z = 7),
        )
    distorted_shells = GaussletBases.cartesian_shells(
        distorted_fixture.parent,
        distorted_fixture.spacing_inputs,
        distorted_fixture.recipe;
        low_order_shellization_policy =
            :terminal_cartesian_shellification_geometry,
    )
    distorted_units = GaussletBases.cartesian_units(
        distorted_fixture.parent,
        distorted_shells,
        distorted_fixture.route_probe_inputs,
        distorted_fixture.recipe,
    )
    distorted_inventory =
        distorted_units.terminal_shellification_unit_inventory
    @test distorted_units.terminal_shellification_lowering_contract_inventory_available
    distorted_lowering_inventory =
        distorted_units.terminal_shellification_lowering_contract_inventory
    @test distorted_units.terminal_shellification_selected_lowering_contract_inventory_available
    distorted_selected_lowering_inventory =
        distorted_units.terminal_shellification_selected_lowering_contract_inventory
    @test distorted_selected_lowering_inventory.route_lowering_family ==
          :white_lindsey_low_order
    @test distorted_selected_lowering_inventory.all_units_have_exactly_one_selected_contract
    distorted_crc_sidecars =
        distorted_units.low_order_units.terminal_shellification_selected_crc_sidecar_summary
    @test distorted_crc_sidecars.selected_contract_count ==
          distorted_selected_lowering_inventory.selected_contract_count
    @test :distorted_product_box_comx in
          distorted_crc_sidecars.missing_sidecar_kinds
    @test :distorted_product_box_comx_sidecar_not_yet_mapped in
          distorted_crc_sidecars.missing_sidecar_reasons
    @test !distorted_crc_sidecars.final_retained_unit_inventory_available
    @test !distorted_crc_sidecars.pair_inventory_available
    distorted_units_for_gap =
        filter(
            record ->
                record.terminal_region_kind == :central_distorted_product_box,
            distorted_inventory.terminal_region_units,
        )
    @test length(distorted_units_for_gap) == 1
    distorted_unit = only(distorted_units_for_gap)
    @test distorted_unit.unit_kind == :central_distorted_product_box_unit
    @test !distorted_unit.identity_lowering_planned
    @test !isnothing(distorted_unit.source_mode_shape)
    @test distorted_unit.source_cpb_plan.source_mode_shape ==
          distorted_unit.source_mode_shape
    distorted_contract = only(
        contract for contract in distorted_lowering_inventory.lowering_contracts
        if contract.lowering_contract_kind == :distorted_product_box_comx
    )
    @test !distorted_contract.identity_like_source_contract
    @test distorted_contract.source_mode_shape == distorted_unit.source_mode_shape
    @test distorted_contract.q == distorted_unit.q
    @test distorted_contract.L == distorted_unit.L
    @test distorted_contract.aspect_ratio == distorted_unit.aspect_ratio
    distorted_selected_contract = only(
        contract for contract in distorted_selected_lowering_inventory.selected_contracts
        if contract.lowering_contract_kind == :distorted_product_box_comx
    )
    @test !distorted_selected_contract.identity_like_source_contract
    @test distorted_selected_contract.contract_key == distorted_contract.contract_key
    @test !distorted_inventory.final_retained_unit_inventory_available
    @test !distorted_inventory.pair_inventory_available
    @test !distorted_lowering_inventory.final_retained_unit_inventory_available
    @test !distorted_lowering_inventory.pair_inventory_available
    @test !distorted_selected_lowering_inventory.final_retained_unit_inventory_available
    @test !distorted_selected_lowering_inventory.pair_inventory_available
end

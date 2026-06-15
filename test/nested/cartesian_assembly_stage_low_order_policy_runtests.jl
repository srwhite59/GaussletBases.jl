using Test
using GaussletBases

function _cartesian_assembly_stage_low_order_policy_fixture(;
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
        route_kind = :be2_cartesian_assembly_stage_low_order_policy,
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
    return (; system, parent, spacing_inputs, recipe, route_probe_inputs)
end

function _cartesian_assembly_stage_low_order_policy_assembly(
    fixture;
    policy = nothing,
)
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
    assembly = GaussletBases.cartesian_assembly(
        fixture.parent,
        shells,
        units,
        transforms,
        pairs,
        fixture.recipe,
    )
    return (; shells, units, transforms, pairs, assembly)
end

@testset "cartesian assembly stage carries selected low-order policy" begin
    fixture = _cartesian_assembly_stage_low_order_policy_fixture()

    default_stages =
        _cartesian_assembly_stage_low_order_policy_assembly(fixture)
    default_assembly = default_stages.assembly
    default_summary = default_assembly.low_order_assembly
    @test default_assembly.object_kind == :cartesian_assembly
    @test default_summary.object_kind ==
          :cartesian_assembly_stage_low_order_summary
    @test default_summary.legacy_source_assembly_selected
    @test default_summary.active_source_authority
    @test !default_summary.atom_growth_assembly_selected
    @test !default_summary.terminal_shellification_assembly_selected
    @test default_summary.pair_inventory_known
    @test !default_summary.independent_atom_growth_pair_inventory_available
    @test default_summary.pair_count == length(default_stages.pairs.pair_entries)
    @test default_summary.pair_family_counts ==
          default_stages.pairs.pair_family_counts
    @test default_summary.assembly_stage_fields_preserved
    @test !default_assembly.atom_growth_assembly_selected
    @test default_assembly.low_order_pair_inventory_known
    @test !default_assembly.low_order_independent_atom_growth_pair_inventory_available
    @test default_assembly.low_order_pair_count ==
          length(default_stages.pairs.pair_entries)
    @test default_assembly.assembly_requires_materialization
    @test default_assembly.active_source_authority
    @test default_assembly.route_skeleton === default_stages.shells.route_skeleton
    @test default_assembly.shells === default_stages.shells
    @test default_assembly.units === default_stages.units
    @test default_assembly.transforms === default_stages.transforms
    @test default_assembly.pairs === default_stages.pairs
    @test hasproperty(default_assembly, :route_facts)
    @test hasproperty(default_assembly, :contract)
    @test default_assembly.assembly_stage == :assembled_report_inputs

    atom_growth_stages =
        _cartesian_assembly_stage_low_order_policy_assembly(
            fixture;
            policy = :atom_growth_complete_rectangular,
        )
    atom_growth_assembly = atom_growth_stages.assembly
    atom_growth_summary = atom_growth_assembly.low_order_assembly
    @test atom_growth_assembly.low_order_assembly_route_kind ==
          :atom_growth_complete_rectangular_low_order_assembly
    @test atom_growth_assembly.low_order_assembly_source ==
          :atom_growth_complete_rectangular_low_order_pair_terms
    @test atom_growth_assembly.atom_growth_assembly_selected
    @test !atom_growth_assembly.hamiltonian_matrices_materialized
    @test !atom_growth_assembly.operator_matrices_materialized
    @test !atom_growth_assembly.pair_operator_blocks_materialized
    @test atom_growth_assembly.assembly_requires_materialization
    @test !atom_growth_assembly.active_source_authority
    @test atom_growth_summary.low_order_shellization_policy_resolved ==
          :atom_growth_complete_rectangular
    @test atom_growth_summary.shellization_source ==
          :bond_aligned_diatomic_atom_growth_construction_plan
    @test atom_growth_summary.shellization_kind ==
          :atom_growth_complete_rectangular
    @test atom_growth_summary.pair_route_kind ==
          :atom_growth_complete_rectangular_low_order_pairs
    @test atom_growth_summary.assembly_source ==
          :atom_growth_complete_rectangular_low_order_pair_terms
    @test atom_growth_summary.assembly_route_kind ==
          :atom_growth_complete_rectangular_low_order_assembly
    @test atom_growth_summary.assembly_kind ==
          :atom_growth_complete_rectangular_low_order
    @test atom_growth_summary.atom_growth_assembly_selected
    @test !atom_growth_summary.terminal_shellification_assembly_selected
    @test !atom_growth_summary.legacy_source_assembly_selected
    @test atom_growth_summary.plan_authority
    @test !atom_growth_summary.active_source_authority
    @test !atom_growth_summary.legacy_source_authority
    @test !atom_growth_summary.hamiltonian_matrices_materialized
    @test !atom_growth_summary.operator_matrices_materialized
    @test !atom_growth_summary.pair_operator_blocks_materialized
    @test !atom_growth_summary.pair_operator_blocks_available
    @test atom_growth_summary.pair_inventory_source ==
          :atom_growth_unit_inventory
    @test atom_growth_summary.pair_inventory_known
    @test atom_growth_summary.independent_atom_growth_pair_inventory_available
    @test atom_growth_summary.pair_count == 36
    @test atom_growth_summary.pair_family_counts.white_lindsey_low_order_atom_growth_unit_pair ==
          36
    @test !atom_growth_summary.assembly_can_proceed_from_current_staged_data
    @test atom_growth_summary.assembly_requires_materialization
    @test atom_growth_summary.assembly_materialization_status ==
          :deferred_atom_growth_complete_rectangular_pair_block_materialization
    @test atom_growth_summary.assembly_blocker ==
          :pair_operator_blocks_deferred
    @test atom_growth_summary.summary_only
    @test atom_growth_summary.assembly_stage_fields_preserved
    expected_atom_growth_helper_status = (
        white_lindsey_low_order_atom_growth_unit_pair =
            :deferred_no_pair_operator_block_helper,
    )
    @test atom_growth_summary.pair_helper_status_by_family ==
          expected_atom_growth_helper_status
    @test atom_growth_assembly.low_order_pair_inventory_source ==
          :atom_growth_unit_inventory
    @test atom_growth_assembly.low_order_pair_inventory_known
    @test atom_growth_assembly.low_order_independent_atom_growth_pair_inventory_available
    @test atom_growth_assembly.low_order_pair_count == 36
    @test atom_growth_assembly.low_order_pair_family_counts.white_lindsey_low_order_atom_growth_unit_pair ==
          36
    @test atom_growth_assembly.low_order_pair_helper_status_by_family ==
          expected_atom_growth_helper_status
    @test atom_growth_assembly.route_skeleton ===
          atom_growth_stages.shells.route_skeleton
    @test hasproperty(atom_growth_assembly, :route_facts)
    @test hasproperty(atom_growth_assembly, :contract)
    @test atom_growth_assembly.shells === atom_growth_stages.shells
    @test atom_growth_assembly.units === atom_growth_stages.units
    @test atom_growth_assembly.transforms === atom_growth_stages.transforms
    @test atom_growth_assembly.pairs === atom_growth_stages.pairs
    @test atom_growth_assembly.assembly_stage == :assembled_report_inputs

    terminal_fixture =
        _cartesian_assembly_stage_low_order_policy_fixture(
            probe_parent_axis_construction = false,
            atom_locations = ((-4.0, 0.0, 0.0), (4.0, 0.0, 0.0)),
            parent_axis_counts = (x = 13, y = 7, z = 7),
        )
    terminal_stages =
        _cartesian_assembly_stage_low_order_policy_assembly(
            terminal_fixture;
            policy = :terminal_cartesian_shellification_geometry,
        )
    terminal_assembly = terminal_stages.assembly
    terminal_summary = terminal_assembly.low_order_assembly
    @test terminal_assembly.low_order_assembly_route_kind ==
          :terminal_shellification_low_order_assembly
    @test terminal_assembly.low_order_assembly_source ==
          :terminal_shellification_pair_terms
    @test terminal_assembly.terminal_shellification_assembly_selected
    @test terminal_assembly.terminal_shellification_region_count ==
          terminal_stages.pairs.terminal_shellification_region_count
    @test terminal_assembly.terminal_shellification_unit_count ==
          terminal_stages.pairs.terminal_shellification_unit_count
    @test terminal_assembly.assembly_requires_materialization
    @test !terminal_assembly.active_source_authority
    @test terminal_assembly.terminal_shellification_central_gap_region_count == 3
    @test terminal_assembly.terminal_shellification_central_midpoint_slab_count ==
          3
    @test terminal_assembly.terminal_shellification_central_distorted_product_box_count ==
          0
    @test terminal_summary.object_kind ==
          :cartesian_assembly_stage_low_order_summary
    @test terminal_summary.status ==
          :deferred_terminal_shellification_assembly_materialization
    @test terminal_summary.low_order_shellization_policy_resolved ==
          :terminal_cartesian_shellification_geometry
    @test terminal_summary.shellization_source ==
          :terminal_cartesian_shellification_geometry
    @test terminal_summary.shellization_kind ==
          :terminal_cartesian_shellification_geometry
    @test terminal_summary.assembly_source ==
          :terminal_shellification_pair_terms
    @test terminal_summary.assembly_route_kind ==
          :terminal_shellification_low_order_assembly
    @test terminal_summary.assembly_kind ==
          :terminal_shellification_low_order
    @test terminal_summary.terminal_shellification_assembly_selected
    @test !terminal_summary.atom_growth_assembly_selected
    @test !terminal_summary.legacy_source_assembly_selected
    @test terminal_summary.terminal_shellification_region_count ==
          terminal_stages.pairs.terminal_shellification_region_count
    @test terminal_summary.terminal_shellification_unit_inventory_available
    @test !terminal_summary.pair_operator_blocks_materialized
    @test terminal_summary.pair_count == 0
    @test !terminal_summary.assembly_can_proceed_from_current_staged_data
    @test terminal_summary.assembly_requires_materialization
    @test terminal_summary.assembly_materialization_status ==
          :deferred_terminal_shellification_assembly_materialization
    @test terminal_summary.assembly_blocker ==
          :terminal_shellification_pair_blocks_deferred
    @test terminal_summary.plan_authority
    @test !terminal_summary.active_source_authority
    @test !terminal_summary.legacy_source_authority
    @test terminal_summary.summary_only
    terminal_units =
        terminal_summary.terminal_shellification_unit_inventory.terminal_region_units
    @test length(terminal_units) ==
          terminal_stages.pairs.terminal_shellification_unit_count
    @test all(
        record -> !record.owned_support_is_cpb,
        terminal_units,
    )
    @test all(
        record -> !record.shellification_region_is_cpb,
        terminal_units,
    )
    @test terminal_summary.assembly_stage_fields_preserved
    @test terminal_assembly.shells === terminal_stages.shells
    @test terminal_assembly.units === terminal_stages.units
    @test terminal_assembly.transforms === terminal_stages.transforms
    @test terminal_assembly.pairs === terminal_stages.pairs
    @test terminal_assembly.assembly_stage == :assembled_report_inputs
end

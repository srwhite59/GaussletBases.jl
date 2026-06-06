using Test
using GaussletBases

function _cartesian_report_stage_low_order_policy_fixture(;
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
        route_kind = :be2_cartesian_report_stage_low_order_policy,
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

function _cartesian_report_stage_low_order_policy_report(
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
    report = GaussletBases.cartesian_report(
        fixture.system,
        fixture.parent,
        assembly,
        fixture.recipe,
    )
    return (; shells, units, transforms, pairs, assembly, report)
end

function _cartesian_report_stage_count_by_field(entries, field, value)
    for entry in entries
        getproperty(entry, field) == value && return entry.pair_count
    end
    return 0
end

@testset "cartesian report stage carries selected low-order policy" begin
    fixture = _cartesian_report_stage_low_order_policy_fixture()

    default_stages = _cartesian_report_stage_low_order_policy_report(fixture)
    default_report = default_stages.report
    default_summary = default_report.low_order_route_summary
    @test default_report.object_kind ==
          :cartesian_nesting_route_driver_skeleton_report
    @test default_summary.object_kind ==
          :cartesian_report_stage_low_order_route_summary
    @test default_summary.low_order_shellization_policy_requested === nothing
    @test default_summary.low_order_shellization_policy_resolved ==
          :legacy_diatomic_source
    @test default_summary.low_order_shellization_policy_source ==
          :default_legacy_diatomic_source
    @test default_summary.shellization_source ==
          :route_configured_bond_aligned_diatomic_source
    @test default_summary.unit_route_kind ==
          :legacy_diatomic_source_low_order_units
    @test default_summary.transform_route_kind ==
          :legacy_diatomic_source_low_order_transforms
    @test default_summary.pair_route_kind ==
          :legacy_diatomic_source_low_order_pairs
    @test default_summary.assembly_source ==
          :legacy_diatomic_source_pair_terms
    @test default_summary.assembly_route_kind ==
          :legacy_diatomic_source_low_order_assembly
    @test default_summary.legacy_source_selected
    @test default_summary.active_source_authority
    @test !default_summary.atom_growth_selected
    @test !default_summary.terminal_shellification_selected
    @test !default_summary.plan_authority
    @test default_summary.materialization_required
    @test default_summary.materialization_status ==
          :deferred_legacy_diatomic_source_pair_block_materialization
    @test !default_summary.hamiltonian_matrices_materialized
    @test !default_summary.operator_matrices_materialized
    @test !default_summary.pair_operator_blocks_materialized
    @test default_summary.pair_inventory_source ==
          :route_skeleton_pair_entries_only
    @test default_summary.pair_inventory_known
    @test !default_summary.independent_atom_growth_pair_inventory_available
    @test default_summary.pair_count == length(default_stages.pairs.pair_entries)
    @test default_summary.pair_family_counts ==
          default_stages.pairs.pair_family_counts
    @test default_summary.route_core_summary_status == :not_selected_legacy_source_pairs
    @test default_summary.route_core_final_unit_count == 0
    @test !default_summary.route_core_pair_inventory_available
    @test default_summary.route_core_pair_count == 0
    @test !default_summary.route_core_pair_order_matches_staged
    @test isempty(default_summary.route_core_pair_family_counts)
    @test !default_summary.route_core_pair_operator_ready
    @test default_summary.route_core_pair_operator_readiness_status ==
          :not_selected_legacy_source_pairs
    @test default_summary.route_core_pair_operator_blocker ==
          :not_selected_legacy_source_pairs
    @test default_summary.report_stage_fields_preserved
    @test default_report.low_order_shellization_policy_resolved ==
          :legacy_diatomic_source
    @test default_report.low_order_shellization_policy_source ==
          :default_legacy_diatomic_source
    @test default_report.low_order_unit_route_kind ==
          :legacy_diatomic_source_low_order_units
    @test default_report.low_order_transform_route_kind ==
          :legacy_diatomic_source_low_order_transforms
    @test default_report.low_order_pair_route_kind ==
          :legacy_diatomic_source_low_order_pairs
    @test default_report.low_order_assembly_route_kind ==
          :legacy_diatomic_source_low_order_assembly
    @test default_report.legacy_source_low_order_route_selected
    @test default_report.low_order_active_source_authority
    @test !default_report.atom_growth_low_order_route_selected
    @test !default_report.low_order_terminal_shellification_selected
    @test !default_report.low_order_plan_authority
    @test default_report.low_order_materialization_required
    @test !default_report.low_order_hamiltonian_matrices_materialized
    @test !default_report.low_order_operator_matrices_materialized
    @test !default_report.low_order_pair_operator_blocks_materialized
    @test default_report.low_order_pair_inventory_source ==
          :route_skeleton_pair_entries_only
    @test default_report.low_order_pair_inventory_known
    @test !default_report.low_order_independent_atom_growth_pair_inventory_available
    @test default_report.low_order_pair_count ==
          length(default_stages.pairs.pair_entries)
    @test default_report.low_order_route_core_summary_status ==
          :not_selected_legacy_source_pairs
    @test default_report.low_order_route_core_final_unit_count == 0
    @test !default_report.low_order_route_core_pair_inventory_available
    @test default_report.low_order_route_core_pair_count == 0
    @test !default_report.low_order_route_core_pair_operator_ready
    @test default_report.low_order_route_core_pair_operator_readiness_status ==
          :not_selected_legacy_source_pairs
    @test default_report.low_order_route_core_pair_operator_preflight_available
    @test default_report.low_order_route_core_pair_operator_preflight_status ==
          :blocked_route_core_pair_operator_preflight
    @test default_report.low_order_route_core_pair_operator_preflight_blocker ==
          :not_selected_legacy_source_pairs
    @test !default_report.low_order_route_core_pair_operator_preflight.operator_blocks_materialized
    @test default_report.low_order_route_core_pair_operator_plan_available
    @test default_report.low_order_route_core_pair_operator_plan_status ==
          :blocked_route_core_pair_operator_plan
    @test default_report.low_order_route_core_pair_operator_plan_blocker ==
          :not_selected_legacy_source_pairs
    @test !default_report.low_order_route_core_pair_operator_plan.operator_blocks_materialized
    @test !default_summary.pqs_lowering_prototype_available
    @test !default_summary.pqs_transform_prototype_available
    @test default_summary.pqs_prototype_stage == :not_available
    @test !default_report.low_order_pqs_lowering_prototype_available
    @test !default_report.low_order_pqs_transform_prototype_available
    @test default_report.low_order_pqs_prototype_stage == :not_available
    @test hasproperty(default_report, :route_materializer_payload)
    @test hasproperty(default_report, :diagnostics)
    @test default_report.route_skeleton === default_stages.shells.route_skeleton
    @test default_report.pair_entries ===
          default_stages.shells.route_skeleton.pair_entries
    @test isnothing(
        GaussletBases._pqs_source_box_route_driver_terminal_shellification_print_line(
            default_report,
        ),
    )
    default_plan_line =
        GaussletBases._pqs_source_box_route_driver_crc_operator_plan_print_line(
            default_report,
        )
    @test occursin("blocked metadata plan", default_plan_line)
    @test occursin("not_selected_legacy_source_pairs", default_plan_line)
    @test occursin("operator blocks materialized no", default_plan_line)

    atom_growth_stages =
        _cartesian_report_stage_low_order_policy_report(
            fixture;
            policy = :atom_growth_complete_rectangular,
        )
    atom_growth_report = atom_growth_stages.report
    atom_growth_summary = atom_growth_report.low_order_route_summary
    @test atom_growth_summary.low_order_shellization_policy_requested ==
          :atom_growth_complete_rectangular
    @test atom_growth_summary.low_order_shellization_policy_resolved ==
          :atom_growth_complete_rectangular
    @test atom_growth_summary.low_order_shellization_policy_source ==
          :explicit_low_order_shellization_policy
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
    @test atom_growth_summary.assembly_source ==
          :atom_growth_complete_rectangular_low_order_pair_terms
    @test atom_growth_summary.assembly_route_kind ==
          :atom_growth_complete_rectangular_low_order_assembly
    @test atom_growth_summary.assembly_kind ==
          :atom_growth_complete_rectangular_low_order
    @test atom_growth_summary.atom_growth_selected
    @test !atom_growth_summary.terminal_shellification_selected
    @test !atom_growth_summary.legacy_source_selected
    @test atom_growth_summary.plan_authority
    @test !atom_growth_summary.active_source_authority
    @test !atom_growth_summary.legacy_source_authority
    @test atom_growth_summary.materialization_required
    @test atom_growth_summary.materialization_status ==
          :deferred_atom_growth_complete_rectangular_pair_block_materialization
    @test atom_growth_summary.materialization_blocker ==
          :pair_operator_blocks_deferred
    @test !atom_growth_summary.hamiltonian_matrices_materialized
    @test !atom_growth_summary.operator_matrices_materialized
    @test !atom_growth_summary.pair_operator_blocks_materialized
    @test atom_growth_summary.pair_inventory_source ==
          :atom_growth_unit_inventory
    @test atom_growth_summary.pair_inventory_known
    @test atom_growth_summary.independent_atom_growth_pair_inventory_available
    @test atom_growth_summary.pair_count == 36
    @test atom_growth_summary.pair_family_counts.white_lindsey_low_order_atom_growth_unit_pair ==
          36
    @test atom_growth_summary.route_core_summary_status ==
          :available_route_core_unit_pair_summary
    @test atom_growth_summary.route_core_final_unit_count == 8
    @test atom_growth_summary.route_core_pair_inventory_available
    @test atom_growth_summary.route_core_pair_inventory_status ==
          :available_route_core_unit_pair_inventory
    @test atom_growth_summary.route_core_pair_count == 36
    @test atom_growth_summary.route_core_pair_order_matches_staged
    @test atom_growth_summary.route_core_pair_order_comparison_source ==
          :atom_growth_pair_inventory
    @test sum(
        family.pair_count for family in atom_growth_summary.route_core_pair_family_counts
    ) == 36
    @test atom_growth_summary.route_core_pair_operator_ready
    @test atom_growth_summary.route_core_pair_operator_readiness_status ==
          :ready_route_core_pair_operator_metadata
    @test isnothing(atom_growth_summary.route_core_pair_operator_blocker)
    @test atom_growth_summary.summary_only
    @test atom_growth_report.low_order_shellization_policy_requested ==
          :atom_growth_complete_rectangular
    @test atom_growth_report.low_order_shellization_policy_source ==
          :explicit_low_order_shellization_policy
    @test atom_growth_report.low_order_shellization_source ==
          :bond_aligned_diatomic_atom_growth_construction_plan
    @test atom_growth_report.low_order_unit_route_kind ==
          :atom_growth_complete_rectangular_low_order_units
    @test atom_growth_report.low_order_transform_route_kind ==
          :atom_growth_complete_rectangular_low_order_transforms
    @test atom_growth_report.low_order_pair_route_kind ==
          :atom_growth_complete_rectangular_low_order_pairs
    @test atom_growth_report.low_order_assembly_source ==
          :atom_growth_complete_rectangular_low_order_pair_terms
    @test atom_growth_report.low_order_assembly_route_kind ==
          :atom_growth_complete_rectangular_low_order_assembly
    @test atom_growth_report.atom_growth_low_order_route_selected
    @test !atom_growth_report.low_order_terminal_shellification_selected
    @test !atom_growth_report.legacy_source_low_order_route_selected
    @test atom_growth_report.low_order_plan_authority
    @test !atom_growth_report.low_order_active_source_authority
    @test !atom_growth_report.low_order_legacy_source_authority
    @test atom_growth_report.low_order_materialization_required
    @test atom_growth_report.low_order_materialization_status ==
          :deferred_atom_growth_complete_rectangular_pair_block_materialization
    @test !atom_growth_report.low_order_hamiltonian_matrices_materialized
    @test !atom_growth_report.low_order_operator_matrices_materialized
    @test !atom_growth_report.low_order_pair_operator_blocks_materialized
    @test atom_growth_report.low_order_pair_inventory_source ==
          :atom_growth_unit_inventory
    @test atom_growth_report.low_order_pair_inventory_known
    @test atom_growth_report.low_order_independent_atom_growth_pair_inventory_available
    @test atom_growth_report.low_order_pair_count == 36
    @test atom_growth_report.low_order_pair_family_counts.white_lindsey_low_order_atom_growth_unit_pair ==
          36
    @test atom_growth_report.low_order_route_core_summary_status ==
          :available_route_core_unit_pair_summary
    @test atom_growth_report.low_order_route_core_final_unit_count == 8
    @test atom_growth_report.low_order_route_core_pair_inventory_available
    @test atom_growth_report.low_order_route_core_pair_inventory_status ==
          :available_route_core_unit_pair_inventory
    @test atom_growth_report.low_order_route_core_pair_count == 36
    @test atom_growth_report.low_order_route_core_pair_order_matches_staged
    @test atom_growth_report.low_order_route_core_pair_order_comparison_source ==
          :atom_growth_pair_inventory
    @test atom_growth_report.low_order_route_core_pair_family_counts ==
          atom_growth_summary.route_core_pair_family_counts
    @test atom_growth_report.low_order_route_core_pair_operator_ready
    @test atom_growth_report.low_order_route_core_pair_operator_readiness_status ==
          :ready_route_core_pair_operator_metadata
    @test isnothing(atom_growth_report.low_order_route_core_pair_operator_blocker)
    @test atom_growth_report.low_order_route_core_pair_operator_preflight_available
    @test atom_growth_report.low_order_route_core_pair_operator_preflight_status ==
          :ready_route_core_pair_operator_preflight
    @test isnothing(
        atom_growth_report.low_order_route_core_pair_operator_preflight_blocker,
    )
    @test atom_growth_report.low_order_route_core_pair_operator_preflight.route_core_final_unit_count ==
          8
    @test atom_growth_report.low_order_route_core_pair_operator_preflight.route_core_pair_count ==
          36
    @test !atom_growth_report.low_order_route_core_pair_operator_preflight.operator_blocks_materialized
    @test atom_growth_report.low_order_route_core_pair_operator_plan_available
    @test atom_growth_report.low_order_route_core_pair_operator_plan_status ==
          :ready_route_core_pair_operator_plan
    @test isnothing(atom_growth_report.low_order_route_core_pair_operator_plan_blocker)
    @test atom_growth_report.low_order_route_core_pair_operator_plan.route_core_final_unit_count ==
          8
    @test atom_growth_report.low_order_route_core_pair_operator_plan.route_core_pair_count ==
          36
    @test !atom_growth_report.low_order_route_core_pair_operator_plan.operator_blocks_materialized
    @test !atom_growth_report.low_order_route_core_pair_operator_plan.hamiltonian_matrices_materialized
    @test atom_growth_summary.route_core_typed_pair_operator_plan_inventory_available
    @test atom_growth_summary.route_core_typed_pair_operator_plan_inventory_status ==
          :blocked_route_core_pair_operator_plan_inventory
    @test atom_growth_summary.route_core_typed_pair_operator_plan_blocker ==
          :aggregate_subtree_operator_plan_required
    @test atom_growth_summary.route_core_typed_pair_operator_plan_count == 36
    @test atom_growth_summary.route_core_typed_pair_operator_plan_blocked_count == 15
    @test !atom_growth_summary.route_core_typed_pair_operator_plan_materialized
    @test _cartesian_report_stage_count_by_field(
        atom_growth_summary.route_core_typed_pair_operator_source_path_counts,
        :source_operator_path,
        :aggregate_subtree_adapter_required,
    ) == 15
    @test _cartesian_report_stage_count_by_field(
        atom_growth_summary.route_core_typed_pair_operator_materialization_status_counts,
        :materialization_status,
        :metadata_only_not_materialized,
    ) == 21
    @test _cartesian_report_stage_count_by_field(
        atom_growth_summary.route_core_typed_pair_operator_materialization_status_counts,
        :materialization_status,
        :blocked_metadata_only_not_materialized,
    ) == 15
    @test _cartesian_report_stage_count_by_field(
        atom_growth_summary.route_core_typed_pair_operator_blocker_counts,
        :blocker,
        nothing,
    ) == 21
    @test _cartesian_report_stage_count_by_field(
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
    @test atom_growth_report.low_order_route_core_typed_pair_operator_plan_count ==
          atom_growth_summary.route_core_typed_pair_operator_plan_count
    @test atom_growth_report.low_order_route_core_typed_pair_operator_plan_blocked_count ==
          atom_growth_summary.route_core_typed_pair_operator_plan_blocked_count
    @test !atom_growth_report.low_order_route_core_typed_pair_operator_plan_materialized
    @test atom_growth_report.low_order_route_core_typed_pair_operator_source_path_counts ==
          atom_growth_summary.route_core_typed_pair_operator_source_path_counts
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
    @test atom_growth_report.low_order_route_core_typed_pair_operator_materialization_ready ==
          atom_growth_summary.route_core_typed_pair_operator_materialization_ready
    @test atom_growth_report.low_order_route_core_typed_pair_operator_materialization_readiness_status ==
          atom_growth_summary.route_core_typed_pair_operator_materialization_readiness_status
    @test atom_growth_report.low_order_route_core_typed_pair_operator_materialization_readiness_blocker ==
          atom_growth_summary.route_core_typed_pair_operator_materialization_readiness_blocker
    @test atom_growth_report.low_order_route_core_typed_pair_operator_materialization_readiness_plan_count ==
          atom_growth_summary.route_core_typed_pair_operator_materialization_readiness_plan_count
    @test atom_growth_report.low_order_route_core_typed_pair_operator_materialization_readiness_blocked_count ==
          atom_growth_summary.route_core_typed_pair_operator_materialization_readiness_blocked_count
    @test atom_growth_report.low_order_route_core_typed_pair_operator_materialization_readiness_materialized_count ==
          atom_growth_summary.route_core_typed_pair_operator_materialization_readiness_materialized_count
    @test atom_growth_summary.lw_complete_shell_cpb_enumeration_available
    @test atom_growth_summary.lw_complete_shell_region_count == 4
    @test atom_growth_summary.lw_complete_shell_cpb_count == 104
    @test atom_growth_summary.lw_complete_shell_cpb_family_counts ==
          (facet_cpb = 24, edge_cpb = 48, corner_cpb = 32)
    @test atom_growth_summary.lw_complete_shell_enumeration_policy ==
          :white_lindsey_complete_shell_boundary_strata
    @test !atom_growth_summary.lw_complete_shell_coefficient_maps_materialized
    @test !atom_growth_summary.lw_complete_shell_operator_blocks_materialized
    @test !atom_growth_summary.lw_complete_shell_pair_operator_blocks_materialized
    @test !atom_growth_summary.lw_complete_shell_hamiltonian_data_materialized
    @test atom_growth_report.low_order_lw_complete_shell_cpb_enumeration_available
    @test atom_growth_report.low_order_lw_complete_shell_region_count == 4
    @test atom_growth_report.low_order_lw_complete_shell_cpb_count == 104
    @test atom_growth_report.low_order_lw_complete_shell_cpb_family_counts ==
          atom_growth_summary.lw_complete_shell_cpb_family_counts
    @test atom_growth_report.low_order_lw_complete_shell_enumeration_policy ==
          :white_lindsey_complete_shell_boundary_strata
    @test !atom_growth_report.low_order_lw_complete_shell_coefficient_maps_materialized
    @test !atom_growth_report.low_order_lw_complete_shell_operator_blocks_materialized
    @test !atom_growth_report.low_order_lw_complete_shell_pair_operator_blocks_materialized
    @test !atom_growth_report.low_order_lw_complete_shell_hamiltonian_data_materialized
    @test atom_growth_summary.pqs_lowering_prototype_available
    @test atom_growth_summary.pqs_transform_prototype_available
    @test atom_growth_summary.pqs_lowering_prototype ===
          atom_growth_stages.units.pqs_lowering_prototype
    @test atom_growth_summary.pqs_transform_prototype ===
          atom_growth_stages.transforms.pqs_transform_prototype
    @test atom_growth_summary.pqs_prototype_unit_key ==
          atom_growth_stages.units.pqs_lowering_prototype_unit_key
    @test atom_growth_summary.pqs_prototype_stage == :metadata_only
    @test atom_growth_summary.pqs_prototype_source ==
          :transform_stage_pqs_transform_prototype
    @test atom_growth_summary.pqs_prototype_source_cpb_kind ==
          :filled_source_cpb
    @test !atom_growth_summary.pqs_prototype_owned_support_is_cpb
    @test atom_growth_summary.pqs_prototype_intermediate_retained_space ==
          :boundary_comx_product_mode_selection
    @test atom_growth_summary.pqs_prototype_shell_realization ==
          :shell_projection_lowdin_deferred
    @test atom_growth_summary.pqs_prototype_source_cpb_support_count ==
          atom_growth_stages.units.pqs_lowering_prototype.source_cpb_support_count
    @test atom_growth_summary.pqs_prototype_owned_support_count ==
          atom_growth_stages.units.pqs_lowering_prototype.owned_support_count
    @test atom_growth_summary.pqs_prototype_source_count_distinct_from_owned_support_count
    @test !atom_growth_summary.pqs_prototype_coefficient_maps_materialized
    @test !atom_growth_summary.pqs_prototype_coefficient_transform_materialized
    @test !atom_growth_summary.pqs_prototype_numerical_transform_materialized
    @test !atom_growth_summary.pqs_prototype_source_operator_blocks_materialized
    @test !atom_growth_summary.pqs_prototype_operator_blocks_materialized
    @test !atom_growth_summary.pqs_prototype_pair_operator_blocks_materialized
    @test !atom_growth_summary.pqs_prototype_hamiltonian_data_materialized
    @test !atom_growth_summary.pqs_prototype_artifacts_materialized
    @test atom_growth_report.low_order_pqs_lowering_prototype_available
    @test atom_growth_report.low_order_pqs_transform_prototype_available
    @test atom_growth_report.low_order_pqs_prototype_unit_key ==
          atom_growth_summary.pqs_prototype_unit_key
    @test atom_growth_report.low_order_pqs_prototype_stage == :metadata_only
    @test atom_growth_report.low_order_pqs_prototype_source_cpb_kind ==
          :filled_source_cpb
    @test !atom_growth_report.low_order_pqs_prototype_owned_support_is_cpb
    @test atom_growth_report.low_order_pqs_prototype_intermediate_retained_space ==
          :boundary_comx_product_mode_selection
    @test atom_growth_report.low_order_pqs_prototype_shell_realization ==
          :shell_projection_lowdin_deferred
    @test atom_growth_report.low_order_pqs_prototype_source_count_distinct_from_owned_support_count
    @test !atom_growth_report.low_order_pqs_prototype_coefficient_maps_materialized
    @test !atom_growth_report.low_order_pqs_prototype_source_operator_blocks_materialized
    @test !atom_growth_report.low_order_pqs_prototype_operator_blocks_materialized
    @test !atom_growth_report.low_order_pqs_prototype_pair_operator_blocks_materialized
    @test !atom_growth_report.low_order_pqs_prototype_hamiltonian_data_materialized
    @test !atom_growth_report.low_order_pqs_prototype_artifacts_materialized
    @test hasproperty(atom_growth_report, :route_materializer_payload)
    @test hasproperty(atom_growth_report, :diagnostics)
    @test atom_growth_report.route_skeleton ===
          atom_growth_stages.shells.route_skeleton
    @test atom_growth_report.pair_entries ===
          atom_growth_stages.shells.route_skeleton.pair_entries
    @test isnothing(
        GaussletBases._pqs_source_box_route_driver_terminal_shellification_print_line(
            atom_growth_report,
        ),
    )

    terminal_fixture = _cartesian_report_stage_low_order_policy_fixture(
        probe_parent_axis_construction = false,
        atom_locations = ((-4.0, 0.0, 0.0), (4.0, 0.0, 0.0)),
        parent_axis_counts = (x = 13, y = 7, z = 7),
    )
    terminal_stages = _cartesian_report_stage_low_order_policy_report(
        terminal_fixture;
        policy = :terminal_cartesian_shellification_geometry,
    )
    terminal_report = terminal_stages.report
    terminal_summary = terminal_report.low_order_route_summary
    @test terminal_summary.object_kind ==
          :cartesian_report_stage_low_order_route_summary
    @test terminal_summary.status ==
          :deferred_terminal_shellification_assembly_materialization
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
    @test terminal_summary.assembly_source ==
          :terminal_shellification_pair_terms
    @test terminal_summary.assembly_route_kind ==
          :terminal_shellification_low_order_assembly
    @test terminal_summary.assembly_kind == :terminal_shellification_low_order
    @test terminal_summary.terminal_shellification_selected
    @test !terminal_summary.atom_growth_selected
    @test !terminal_summary.legacy_source_selected
    @test terminal_summary.terminal_shellification_summary_available
    @test terminal_summary.terminal_shellification_scaffold_available
    @test terminal_summary.terminal_shellification_scaffold ===
          terminal_stages.assembly.terminal_shellification_scaffold
    @test terminal_summary.terminal_shellification_region_count ==
          terminal_stages.assembly.terminal_shellification_region_count
    @test terminal_summary.terminal_shellification_unit_inventory_available
    @test terminal_summary.terminal_shellification_unit_inventory ===
          terminal_stages.assembly.terminal_shellification_unit_inventory
    @test terminal_summary.terminal_shellification_unit_count ==
          terminal_stages.assembly.terminal_shellification_unit_count
    @test terminal_summary.terminal_shellification_unit_keys ==
          terminal_stages.assembly.terminal_shellification_unit_keys
    @test terminal_summary.terminal_shellification_unit_roles ==
          terminal_stages.assembly.terminal_shellification_unit_roles
    @test terminal_summary.terminal_shellification_unit_kinds ==
          terminal_stages.assembly.terminal_shellification_unit_kinds
    @test terminal_summary.terminal_shellification_unit_support_counts ==
          terminal_stages.assembly.terminal_shellification_unit_support_counts
    @test !terminal_summary.terminal_shellification_final_retained_unit_inventory_available
    @test !terminal_summary.terminal_shellification_transform_contracts_available
    @test !terminal_summary.terminal_shellification_pair_inventory_available
    @test terminal_summary.terminal_shellification_pair_inventory_status ==
          :deferred_terminal_shellification_pair_inventory
    @test terminal_summary.terminal_shellification_pair_materialization_status ==
          :deferred_terminal_shellification_pair_materialization
    @test terminal_summary.terminal_shellification_assembly_materialization_status ==
          :deferred_terminal_shellification_assembly_materialization
    @test terminal_summary.materialization_required
    @test terminal_summary.materialization_status ==
          :deferred_terminal_shellification_assembly_materialization
    @test terminal_summary.materialization_blocker ==
          :terminal_shellification_pair_blocks_deferred
    @test !terminal_summary.hamiltonian_matrices_materialized
    @test !terminal_summary.operator_matrices_materialized
    @test !terminal_summary.pair_operator_blocks_materialized
    @test terminal_summary.pair_inventory_source ==
          :terminal_shellification_scaffold
    @test !terminal_summary.pair_inventory_known
    @test !terminal_summary.independent_atom_growth_pair_inventory_available
    @test terminal_summary.pair_count == 0
    @test terminal_summary.pair_family_counts == ()
    @test !terminal_summary.route_core_pair_inventory_available
    @test terminal_summary.route_core_pair_inventory_status ==
          :deferred_terminal_shellification_pair_inventory
    @test terminal_summary.route_core_pair_count == 0
    @test terminal_summary.route_core_pair_family_counts == ()
    @test terminal_summary.route_core_summary_status ==
          :deferred_terminal_shellification_pair_inventory
    @test !terminal_summary.route_core_pair_operator_ready
    @test terminal_summary.route_core_pair_operator_readiness_status ==
          :deferred_terminal_shellification_pair_inventory
    @test terminal_summary.route_core_pair_operator_blocker ==
          :deferred_terminal_shellification_pair_inventory
    @test !terminal_summary.route_core_typed_pair_operator_plan_inventory_available
    @test terminal_summary.route_core_typed_pair_operator_plan_inventory_status ==
          :deferred_terminal_shellification_typed_pair_operator_plan_inventory
    @test terminal_summary.route_core_typed_pair_operator_plan_blocker ==
          :deferred_terminal_shellification_pair_inventory
    @test !terminal_summary.route_core_typed_pair_operator_plan_materialized
    @test !terminal_summary.route_core_typed_pair_operator_materialization_ready
    @test terminal_summary.route_core_typed_pair_operator_materialization_readiness_status ==
          :deferred_terminal_shellification_typed_pair_operator_plan_inventory
    @test terminal_summary.terminal_shellification_central_gap_region_count ==
          terminal_stages.assembly.terminal_shellification_central_gap_region_count
    @test terminal_summary.terminal_shellification_central_midpoint_slab_count ==
          terminal_stages.assembly.terminal_shellification_central_midpoint_slab_count
    @test terminal_summary.terminal_shellification_central_distorted_product_box_count ==
          terminal_stages.assembly.terminal_shellification_central_distorted_product_box_count
    @test terminal_summary.terminal_shellification_central_gap_region_count == 3
    @test terminal_summary.terminal_shellification_central_midpoint_slab_count == 3
    @test terminal_summary.terminal_shellification_central_distorted_product_box_count ==
          0
    @test terminal_summary.terminal_shellification_central_distorted_product_box_metadata ==
          terminal_stages.assembly.terminal_shellification_central_distorted_product_box_metadata
    @test terminal_summary.plan_authority
    @test !terminal_summary.active_source_authority
    @test !terminal_summary.legacy_source_authority
    @test all(
        record -> !record.owned_support_is_cpb,
        terminal_summary.terminal_shellification_unit_inventory.terminal_region_units,
    )
    @test all(
        record -> !record.shellification_region_is_cpb,
        terminal_summary.terminal_shellification_unit_inventory.terminal_region_units,
    )
    @test terminal_summary.report_stage_fields_preserved
    @test terminal_summary.summary_only
    @test terminal_report.low_order_terminal_shellification_selected
    @test terminal_report.low_order_terminal_shellification_summary_available
    @test terminal_report.low_order_terminal_shellification_scaffold_available
    @test terminal_report.low_order_terminal_shellification_scaffold ===
          terminal_summary.terminal_shellification_scaffold
    @test terminal_report.low_order_terminal_shellification_region_count ==
          terminal_summary.terminal_shellification_region_count
    @test terminal_report.low_order_terminal_shellification_unit_inventory_available
    @test terminal_report.low_order_terminal_shellification_unit_inventory ===
          terminal_summary.terminal_shellification_unit_inventory
    @test terminal_report.low_order_terminal_shellification_unit_count ==
          terminal_summary.terminal_shellification_unit_count
    @test terminal_report.low_order_terminal_shellification_unit_keys ==
          terminal_summary.terminal_shellification_unit_keys
    @test terminal_report.low_order_terminal_shellification_unit_roles ==
          terminal_summary.terminal_shellification_unit_roles
    @test terminal_report.low_order_terminal_shellification_unit_kinds ==
          terminal_summary.terminal_shellification_unit_kinds
    @test terminal_report.low_order_terminal_shellification_unit_support_counts ==
          terminal_summary.terminal_shellification_unit_support_counts
    @test !terminal_report.low_order_terminal_shellification_final_retained_unit_inventory_available
    @test !terminal_report.low_order_terminal_shellification_transform_contracts_available
    @test !terminal_report.low_order_terminal_shellification_pair_inventory_available
    @test terminal_report.low_order_terminal_shellification_pair_inventory_status ==
          :deferred_terminal_shellification_pair_inventory
    @test terminal_report.low_order_terminal_shellification_pair_materialization_status ==
          :deferred_terminal_shellification_pair_materialization
    @test terminal_report.low_order_terminal_shellification_assembly_materialization_status ==
          :deferred_terminal_shellification_assembly_materialization
    @test terminal_report.low_order_terminal_shellification_central_gap_region_count ==
          terminal_summary.terminal_shellification_central_gap_region_count
    @test terminal_report.low_order_terminal_shellification_central_midpoint_slab_count ==
          terminal_summary.terminal_shellification_central_midpoint_slab_count
    @test terminal_report.low_order_terminal_shellification_central_distorted_product_box_count ==
          terminal_summary.terminal_shellification_central_distorted_product_box_count
    @test terminal_report.low_order_materialization_required
    @test terminal_report.low_order_materialization_status ==
          :deferred_terminal_shellification_assembly_materialization
    @test terminal_report.low_order_materialization_blocker ==
          :terminal_shellification_pair_blocks_deferred
    @test !terminal_report.low_order_hamiltonian_matrices_materialized
    @test !terminal_report.low_order_operator_matrices_materialized
    @test !terminal_report.low_order_pair_operator_blocks_materialized
    @test terminal_report.low_order_pair_inventory_source ==
          :terminal_shellification_scaffold
    @test !terminal_report.low_order_pair_inventory_known
    @test !terminal_report.low_order_independent_atom_growth_pair_inventory_available
    @test terminal_report.low_order_pair_count == 0
    @test terminal_report.low_order_pair_family_counts == ()
    @test !terminal_report.low_order_route_core_pair_inventory_available
    @test terminal_report.low_order_route_core_pair_inventory_status ==
          :deferred_terminal_shellification_pair_inventory
    @test terminal_report.low_order_route_core_pair_count == 0
    @test terminal_report.low_order_route_core_summary_status ==
          :deferred_terminal_shellification_pair_inventory
    @test !terminal_report.low_order_route_core_pair_operator_ready
    @test terminal_report.low_order_route_core_pair_operator_readiness_status ==
          :deferred_terminal_shellification_pair_inventory
    @test !terminal_report.low_order_route_core_typed_pair_operator_plan_inventory_available
    @test terminal_report.low_order_route_core_typed_pair_operator_plan_inventory_status ==
          :deferred_terminal_shellification_typed_pair_operator_plan_inventory
    @test terminal_report.low_order_route_core_typed_pair_operator_plan_blocker ==
          :deferred_terminal_shellification_pair_inventory
    terminal_print_line =
        GaussletBases._pqs_source_box_route_driver_terminal_shellification_print_line(
            terminal_report,
        )
    @test !isnothing(terminal_print_line)
    @test occursin("Terminal shellification: selected", terminal_print_line)
    @test occursin(
        "regions $(terminal_report.low_order_terminal_shellification_region_count)",
        terminal_print_line,
    )
    @test occursin(
        "central gaps $(terminal_report.low_order_terminal_shellification_central_gap_region_count)",
        terminal_print_line,
    )
    @test occursin(
        "midpoint slabs $(terminal_report.low_order_terminal_shellification_central_midpoint_slab_count)",
        terminal_print_line,
    )
    @test occursin(
        "central distorted product boxes $(terminal_report.low_order_terminal_shellification_central_distorted_product_box_count)",
        terminal_print_line,
    )
    @test occursin(
        "pair inventory :deferred_terminal_shellification_pair_inventory",
        terminal_print_line,
    )
    @test occursin(
        "assembly/materialization :deferred_terminal_shellification_assembly_materialization",
        terminal_print_line,
    )
    @test occursin("operators/materialization=no", terminal_print_line)

    atom_growth_materialization = GaussletBases.cartesian_materialization(
        atom_growth_report,
        (;
            low_order_shellization_policy =
                :atom_growth_complete_rectangular,
        ),
    )
    summary_stdout = mktemp() do path, io
        redirect_stdout(io) do
            GaussletBases.cartesian_print_summary(
                atom_growth_report,
                atom_growth_materialization,
            )
        end
        flush(io)
        read(path, String)
    end
    @test occursin(
        "report.low_order_shellization_policy_resolved = :atom_growth_complete_rectangular",
        summary_stdout,
    )
    @test occursin(
        "report.low_order_shellization_policy_source = :explicit_low_order_shellization_policy",
        summary_stdout,
    )
    @test occursin(
        "report.atom_growth_low_order_route_selected = true",
        summary_stdout,
    )
    @test occursin(
        "report.low_order_active_source_authority = false",
        summary_stdout,
    )
    @test occursin(
        "report.low_order_materialization_required = true",
        summary_stdout,
    )
    @test occursin(
        "report.low_order_materialization_status = :deferred_atom_growth_complete_rectangular_pair_block_materialization",
        summary_stdout,
    )
    @test occursin(
        "report.low_order_pair_inventory_source = :atom_growth_unit_inventory",
        summary_stdout,
    )
    @test occursin(
        "report.low_order_pair_count = 36",
        summary_stdout,
    )
    @test occursin(
        "CRC sidecars: final units 8, pairs 36, order match yes",
        summary_stdout,
    )
    @test occursin("CRC pair families: ", summary_stdout)
    @test occursin(
        "white_lindsey_atom_local_child_shellification",
        summary_stdout,
    )
    @test occursin(
        "CRC pair-operator plan: ready metadata plan, final units 8, pairs 36, operator blocks materialized no",
        summary_stdout,
    )
    @test occursin(
        "CRC typed pair-operator inventory: available (:blocked_route_core_pair_operator_plan_inventory), typed plans 36, blocked 15, materialized no, blocker :aggregate_subtree_operator_plan_required",
        summary_stdout,
    )
    @test occursin(
        "report.low_order_hamiltonian_matrices_materialized = false",
        summary_stdout,
    )
    expected_pqs_line = string(
        "CPB/PQS prototype: metadata-only, unit=",
        atom_growth_report.low_order_pqs_prototype_unit_key,
        ", source=filled CPB, owned support is shell, ",
        "retained space=boundary COMX products, ",
        "realization=shell projection + Lowdin deferred, ",
        "operators/materialization=no",
    )
    @test occursin(expected_pqs_line, summary_stdout)
end

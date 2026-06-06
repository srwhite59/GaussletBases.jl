using Test
using GaussletBases

function _cartesian_shell_stage_low_order_policy_fixture(;
    probe_parent_axis_construction = :auto,
    atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
    parent_axis_counts = (x = 9, y = 7, z = 9),
    q::Int = 5,
    n_s::Int = q,
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
        q,
        n_s,
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
    route_inputs = (;
        route_family = :white_lindsey_low_order,
        route_kind = :be2_cartesian_shell_stage_low_order_policy,
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
    return (; parent, spacing_inputs, recipe)
end

@testset "cartesian shell stage carries selected low-order policy" begin
    fixture = _cartesian_shell_stage_low_order_policy_fixture()

    default_shells = GaussletBases.cartesian_shells(
        fixture.parent,
        fixture.spacing_inputs,
        fixture.recipe,
    )
    default_summary = default_shells.low_order_shellization
    @test default_shells.object_kind == :cartesian_shells
    @test default_summary.object_kind ==
          :cartesian_shell_stage_low_order_shellization_summary
    @test default_summary.low_order_shellization_policy_requested === nothing
    @test default_summary.low_order_shellization_policy_resolved ==
          :legacy_diatomic_source
    @test default_summary.low_order_shellization_policy_source ==
          :default_legacy_diatomic_source
    @test default_summary.low_order_shellization_policy_status ==
          :available_low_order_shellization_policy
    @test default_summary.shellization_source ==
          :route_configured_bond_aligned_diatomic_source
    @test default_summary.shellization_kind == :legacy_diatomic_source
    @test default_shells.shellization_stage == :legacy_diatomic_source_summary
    @test default_summary.legacy_source_selected
    @test default_summary.active_source_authority
    @test !default_summary.atom_growth_selected
    @test !default_summary.terminal_shellification_selected
    @test !default_summary.atom_growth_plan_authority
    @test !default_summary.terminal_shellification_authority
    @test !default_summary.atom_growth_plan_summary_available
    @test !default_shells.terminal_shellification_selected
    @test !default_shells.terminal_shellification_plan_available
    @test !default_shells.terminal_shellification_scaffold_available

    atom_growth_shells = GaussletBases.cartesian_shells(
        fixture.parent,
        fixture.spacing_inputs,
        fixture.recipe;
        low_order_shellization_policy = :atom_growth_complete_rectangular,
    )
    atom_growth_summary = atom_growth_shells.low_order_shellization
    @test atom_growth_shells.low_order_shellization_policy_resolved ==
          :atom_growth_complete_rectangular
    @test atom_growth_shells.shellization_source ==
          :bond_aligned_diatomic_atom_growth_construction_plan
    @test atom_growth_shells.shellization_kind ==
          :atom_growth_complete_rectangular
    @test atom_growth_shells.atom_growth_plan_summary_available
    @test atom_growth_shells.atom_growth_plan_available
    @test atom_growth_shells.atom_growth_scaffold_available
    @test atom_growth_shells.atom_growth_shellification_plan_status ==
          :available_atom_growth_shellification_plan
    @test atom_growth_shells.shellization_stage ==
          :available_atom_growth_shellification_plan
    @test atom_growth_shells.atom_growth_plan_authority
    @test !atom_growth_shells.terminal_shellification_selected
    @test !atom_growth_shells.terminal_shellification_plan_available
    @test !atom_growth_shells.terminal_shellification_scaffold_available
    @test !atom_growth_shells.terminal_shellification_authority
    @test !atom_growth_shells.active_source_authority
    @test atom_growth_shells.coverage_complete
    @test atom_growth_shells.materialization_available
    @test atom_growth_shells.full_plan_stored
    @test atom_growth_shells.scaffold_stored
    @test !atom_growth_shells.summary_only
    @test atom_growth_shells.atom_growth_region_count > 0
    @test atom_growth_shells.atom_growth_unsupported_region_count == 0
    @test atom_growth_summary.atom_growth_selected
    @test !atom_growth_summary.legacy_source_selected
    @test atom_growth_summary.status ==
          :available_atom_growth_shellification_plan
    @test atom_growth_summary.low_order_shellization_policy_requested ==
          :atom_growth_complete_rectangular
    @test atom_growth_summary.low_order_shellization_policy_source ==
          :explicit_low_order_shellization_policy
    @test atom_growth_summary.low_order_shellization_policy_status ==
          :available_low_order_shellization_policy
    @test atom_growth_summary.atom_growth_plan_summary_available
    @test atom_growth_summary.atom_growth_plan_available
    @test atom_growth_summary.atom_growth_scaffold_available
    @test atom_growth_summary.atom_growth_shellification_plan_available
    @test atom_growth_summary.atom_growth_shellification_plan_status ==
          :available_atom_growth_shellification_plan
    @test atom_growth_summary.atom_growth_plan_payload.object_kind ==
          :cartesian_shell_stage_atom_growth_plan_payload
    @test atom_growth_summary.atom_growth_plan_payload.status ==
          :available_atom_growth_shellification_plan
    @test atom_growth_summary.atom_growth_construction_plan isa
          GaussletBases._BondAlignedDiatomicAtomGrowthConstructionPlan3D
    @test atom_growth_summary.atom_growth_scaffold.object_kind ==
          :cartesian_atom_growth_shellification_plan3d
    @test atom_growth_summary.atom_growth_scaffold.source_kind ==
          :bond_aligned_diatomic_atom_growth_construction_plan
    @test atom_growth_summary.atom_growth_scaffold.spatial_policy_order ==
          :atom_outward
    @test atom_growth_summary.atom_growth_scaffold.diagnostics.atom_growth_construction_plan_authority
    @test !atom_growth_summary.atom_growth_scaffold.diagnostics.active_source_authority
    @test atom_growth_summary.atom_growth_missing_parent_objects == ()
    @test atom_growth_summary.atom_growth_region_count ==
          atom_growth_summary.atom_growth_scaffold.region_count
    @test atom_growth_summary.atom_growth_unsupported_region_count == 0
    @test atom_growth_summary.atom_growth_spatial_policy_order == :atom_outward
    @test atom_growth_summary.atom_growth_construction_region_order ==
          Tuple(atom_growth_summary.atom_growth_construction_plan.region_order)
    @test atom_growth_summary.atom_growth_materialization_dependency_counts.unsupported_region_count ==
          0
    @test atom_growth_summary.atom_growth_plan_authority
    @test !atom_growth_summary.terminal_shellification_selected
    @test !atom_growth_summary.terminal_shellification_authority
    @test !atom_growth_summary.active_source_authority
    @test atom_growth_summary.coverage_status == :coverage_complete
    @test atom_growth_summary.coverage_complete
    @test atom_growth_summary.materialization_required
    @test atom_growth_summary.materialization_available
    @test atom_growth_summary.materialization_status ==
          :available_atom_growth_shellification_plan
    @test atom_growth_summary.materialization_stage ==
          :atom_growth_shellification_plan_stored_sequence_not_materialized
    @test !atom_growth_summary.materialized_sequence_available
    @test !atom_growth_summary.materialized_operator_matrices_available
    @test !atom_growth_summary.summary_only
    @test atom_growth_summary.full_plan_stored
    @test atom_growth_summary.scaffold_stored

    missing_parent_fixture =
        _cartesian_shell_stage_low_order_policy_fixture(
            probe_parent_axis_construction = false,
        )
    missing_parent_shells = GaussletBases.cartesian_shells(
        missing_parent_fixture.parent,
        missing_parent_fixture.spacing_inputs,
        missing_parent_fixture.recipe;
        low_order_shellization_policy = :atom_growth_complete_rectangular,
    )
    missing_parent_summary = missing_parent_shells.low_order_shellization
    @test missing_parent_shells.shellization_stage ==
          :blocked_atom_growth_missing_parent_objects
    @test missing_parent_summary.status ==
          :blocked_atom_growth_missing_parent_objects
    @test !missing_parent_summary.atom_growth_plan_available
    @test !missing_parent_summary.atom_growth_scaffold_available
    @test !missing_parent_summary.materialization_available
    @test missing_parent_summary.full_plan_stored == false
    @test missing_parent_summary.scaffold_stored == false
    @test missing_parent_summary.summary_only
    @test missing_parent_summary.atom_growth_missing_parent_objects ==
          (:parent_qw_basis_object, :parent_axis_bundle_object)

    terminal_fixture =
        _cartesian_shell_stage_low_order_policy_fixture(
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
    terminal_summary = terminal_shells.low_order_shellization
    @test terminal_shells.low_order_shellization_policy_resolved ==
          :terminal_cartesian_shellification_geometry
    @test terminal_shells.low_order_shellization_policy_status ==
          :available_low_order_shellization_policy
    @test terminal_shells.shellization_source ==
          :terminal_cartesian_shellification_geometry
    @test terminal_shells.shellization_kind ==
          :terminal_cartesian_shellification_geometry
    @test terminal_shells.shellization_stage ==
          :available_terminal_shellification_scaffold
    @test terminal_shells.terminal_shellification_selected
    @test terminal_shells.terminal_shellification_plan_available
    @test terminal_shells.terminal_shellification_scaffold_available
    @test terminal_shells.terminal_shellification_authority
    @test !terminal_shells.atom_growth_plan_authority
    @test !terminal_shells.active_source_authority
    @test terminal_shells.full_plan_stored
    @test terminal_shells.scaffold_stored
    @test !terminal_shells.summary_only
    @test terminal_shells.coverage_complete
    @test !terminal_shells.materialization_available
    @test terminal_shells.terminal_shellification_materialization_status ==
          :available_terminal_shellification_scaffold
    @test terminal_shells.terminal_shellification_spatial_policy_order ==
          :atom_outward
    @test terminal_shells.terminal_shellification_coverage_complete
    @test terminal_shells.terminal_shellification_central_gap_region_count == 3
    @test terminal_shells.terminal_shellification_central_midpoint_slab_count == 3
    @test terminal_shells.terminal_shellification_central_distorted_product_box_count == 0
    @test terminal_shells.terminal_shellification_region_count > 0
    @test terminal_shells.terminal_shellification_scaffold.object_kind ==
          :cartesian_terminal_shellification_scaffold3d

    @test terminal_summary.status == :available_terminal_shellification_scaffold
    @test terminal_summary.terminal_shellification_selected
    @test terminal_summary.terminal_shellification_plan_available
    @test terminal_summary.terminal_shellification_scaffold_available
    @test terminal_summary.terminal_shellification_scaffold.object_kind ==
          :cartesian_terminal_shellification_scaffold3d
    @test terminal_summary.terminal_shellification_region_count ==
          terminal_summary.terminal_shellification_scaffold.region_count
    @test terminal_summary.terminal_shellification_materialization_status ==
          :available_terminal_shellification_scaffold
    @test terminal_summary.terminal_shellification_spatial_policy_order ==
          :atom_outward
    @test terminal_summary.terminal_shellification_coverage_complete
    @test terminal_summary.coverage_status == :coverage_complete
    @test terminal_summary.coverage_complete
    @test terminal_summary.materialization_stage ==
          :deferred_to_route_materialization
    @test !terminal_summary.materialization_available
    @test terminal_summary.materialized_sequence_available == false
    @test terminal_summary.materialized_operator_matrices_available == false
    @test terminal_summary.terminal_shellification_scaffold.diagnostics.terminal_geometry_authority
    @test !terminal_summary.terminal_shellification_scaffold.diagnostics.shellification_regions_are_cpbs
    @test !terminal_summary.terminal_shellification_scaffold.diagnostics.shellification_regions_are_lowering_sources
    @test all(
        !region.shellification_region_is_cpb
        for region in terminal_summary.terminal_shellification_scaffold.regions
    )
    @test all(
        !region.shellification_region_is_lowering_source
        for region in terminal_summary.terminal_shellification_scaffold.regions
    )

    distorted_plan = GaussletBases._cartesian_terminal_shellification_geometry(
        (collect(1:21), collect(1:7), collect(1:7)),
        ((3, 4, 4), (19, 4, 4));
        core_side = 3,
        q = 3,
    )
    distorted_scaffold =
        GaussletBases._cartesian_terminal_shellification_geometry_scaffold(
            distorted_plan,
        )
    @test distorted_scaffold.materialization_status ==
          :deferred_pending_distorted_product_box_lowering
    @test distorted_scaffold.central_distorted_product_box_count == 1
    @test distorted_scaffold.deferred_region_count == 1
    central_region = only(
        region for region in distorted_scaffold.regions
        if region.region_kind == :central_distorted_product_box
    )
    @test !central_region.independently_lowerable
    @test central_region.missing_independent_lowering_reason ==
          :distorted_product_box_lowering_not_implemented
    @test central_region.metadata.q == 3
    @test central_region.metadata.L == 7
    @test central_region.metadata.source_mode_shape == (7, 3, 3)
end

using Test
using GaussletBases

const CRC = GaussletBases.CartesianRouteCore

function _cartesian_route_core_adapter_fixture()
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
        route_kind = :be2_cartesian_route_core_adapter_sidecar,
        route_shape = (:pqs_left, :product, :pqs_right),
        product_body_rule = :centered_single_z_slab,
        pqs_retained_rule = :boundary_comx_product_mode_selection,
        product_retained_rule = :product_doside_retained_unit,
        terms = (:overlap, :position_x, :position_y, :position_z, :kinetic),
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
    shells = GaussletBases.cartesian_shells(
        parent,
        spacing_inputs,
        recipe;
        low_order_shellization_policy = :atom_growth_complete_rectangular,
    )
    units = GaussletBases.cartesian_units(
        parent,
        shells,
        route_probe_inputs,
        recipe,
    )
    return units
end

@testset "CartesianRouteCore staged adapter sidecars" begin
    units = _cartesian_route_core_adapter_fixture()
    plan_inventory = units.plan_unit_inventory
    @test plan_inventory.status == :available_atom_growth_plan_unit_inventory

    lw_unit = first(
        unit for unit in plan_inventory.plan_units
        if unit.lowering_family == :white_lindsey_adaptive_complete_shell
    )
    lw_sidecar = GaussletBases._cartesian_route_core_sidecar(lw_unit)
    @test lw_sidecar.object_kind == :cartesian_route_core_sidecar
    @test lw_sidecar.sidecar_source ==
          :atom_growth_lw_complete_shell_plan_unit
    @test !lw_sidecar.numerical_behavior_changed
    @test !lw_sidecar.materialization_behavior_changed
    @test lw_sidecar.shellification_region isa CRC.ShellificationRegion
    @test lw_sidecar.lowering_source isa CRC.LoweringSource
    @test lw_sidecar.intermediate_retained_space isa CRC.IntermediateRetainedSpace
    @test lw_sidecar.shell_realization isa CRC.ShellRealization
    @test lw_sidecar.final_retained_unit isa CRC.FinalRetainedUnit
    @test CRC.support_count(lw_sidecar.shellification_region) == lw_unit.support_count
    @test CRC.lowering_recipe(lw_sidecar.lowering_source) ==
          :white_lindsey_boundary_strata
    lw_source_cpbs = CRC.source_cpbs(lw_sidecar.lowering_source)
    @test length(lw_source_cpbs) == 26
    @test count(cpb -> cpb.metadata.cpb_family === :facet_cpb, lw_source_cpbs) == 6
    @test count(cpb -> cpb.metadata.cpb_family === :edge_cpb, lw_source_cpbs) == 12
    @test count(cpb -> cpb.metadata.cpb_family === :corner_cpb, lw_source_cpbs) == 8
    @test sum(CRC.support_count, lw_source_cpbs; init = 0) ==
          CRC.support_count(lw_sidecar.shellification_region)
    @test lw_sidecar.intermediate_retained_space.retained_rule ==
          :white_lindsey_boundary_stratum_product
    @test lw_sidecar.shell_realization.realization_kind ==
          :direct_or_trivial_embedding
    @test lw_sidecar.final_retained_unit.unit_key == lw_unit.unit_key

    pqs_prototype = plan_inventory.pqs_lowering_prototype
    pqs_sidecar = GaussletBases._cartesian_route_core_sidecar(pqs_prototype)
    @test pqs_sidecar.object_kind == :cartesian_route_core_sidecar
    @test pqs_sidecar.sidecar_source == :pqs_lowering_metadata_prototype
    @test !pqs_sidecar.numerical_behavior_changed
    @test !pqs_sidecar.materialization_behavior_changed
    @test pqs_sidecar.shellification_region isa CRC.ShellificationRegion
    @test pqs_sidecar.lowering_source isa CRC.LoweringSource
    @test pqs_sidecar.intermediate_retained_space isa CRC.IntermediateRetainedSpace
    @test pqs_sidecar.shell_realization isa CRC.ShellRealization
    @test pqs_sidecar.final_retained_unit isa CRC.FinalRetainedUnit
    @test CRC.lowering_recipe(pqs_sidecar.lowering_source) ==
          :pqs_filled_source_cpb
    pqs_source_cpbs = CRC.source_cpbs(pqs_sidecar.lowering_source)
    @test length(pqs_source_cpbs) == 1
    pqs_source_cpb = only(pqs_source_cpbs)
    @test CRC.codimension(pqs_source_cpb) == 0
    @test CRC.support_count(pqs_source_cpb) ==
          pqs_prototype.source_cpb_support_count
    @test CRC.support_count(pqs_sidecar.shellification_region) ==
          pqs_prototype.owned_support_count
    @test CRC.support_count(pqs_source_cpb) !=
          CRC.support_count(pqs_sidecar.shellification_region)
    @test pqs_sidecar.intermediate_retained_space.retained_rule ==
          :pqs_boundary_comx_product_modes
    @test pqs_sidecar.intermediate_retained_space.dimension ==
          CRC.boundary_product_mode_count(CRC.shape(pqs_source_cpb))
    @test pqs_sidecar.shell_realization.realization_kind ==
          :shell_projection_lowdin
    @test pqs_sidecar.final_retained_unit.unit_key == pqs_prototype.unit_key

    sidecar_inventory =
        GaussletBases._cartesian_route_core_sidecar_inventory(plan_inventory)
    @test sidecar_inventory.object_kind ==
          :cartesian_route_core_sidecar_inventory
    @test sidecar_inventory.status ==
          :blocked_incomplete_route_core_sidecar_inventory
    @test sidecar_inventory.unit_count == plan_inventory.unit_count
    @test sidecar_inventory.supported_unit_count == 6
    @test sidecar_inventory.unsupported_unit_count == 2
    @test sidecar_inventory.final_unit_count ==
          sidecar_inventory.supported_unit_count
    @test all(
        entry -> entry.route_core_sidecar_available,
        sidecar_inventory.supported_entries,
    )
    @test all(
        entry -> entry.route_core_sidecar.object_kind ==
                 :cartesian_route_core_sidecar,
        sidecar_inventory.supported_entries,
    )
    @test all(
        entry -> entry.route_core_sidecar.final_retained_unit isa
                 CRC.FinalRetainedUnit,
        sidecar_inventory.supported_entries,
    )
    @test all(
        unit -> unit isa CRC.FinalRetainedUnit,
        sidecar_inventory.final_units,
    )
    @test sidecar_inventory.unsupported_unit_keys ==
          (:left_atom_box, :right_atom_box)
    @test sidecar_inventory.missing_route_core_sidecar_reasons ==
          (
              :atom_local_child_shellification_not_yet_split_into_crc_child_units,
              :atom_local_child_shellification_not_yet_split_into_crc_child_units,
          )
    @test !sidecar_inventory.crc_pair_inventory_available
    @test sidecar_inventory.crc_pair_inventory_status ==
          :blocked_missing_unit_sidecars
    @test sidecar_inventory.crc_pair_inventory === nothing
    @test isempty(sidecar_inventory.crc_pair_keys)
    @test length(sidecar_inventory.staged_pair_keys) ==
          plan_inventory.unit_count * (plan_inventory.unit_count + 1) ÷ 2
    @test !sidecar_inventory.pair_inventory_order_matches_staged
    @test sidecar_inventory.pqs_prototype_sidecar_available
    @test sidecar_inventory.pqs_prototype_sidecar.final_retained_unit isa
          CRC.FinalRetainedUnit
    @test !sidecar_inventory.numerical_behavior_changed
    @test !sidecar_inventory.materialization_behavior_changed

    direct_entry = only(
        entry for entry in sidecar_inventory.supported_entries
        if entry.unit_key == :contact_cap
    )
    @test direct_entry.route_core_sidecar.sidecar_source ==
          :atom_growth_direct_cpb_plan_unit
    @test CRC.lowering_recipe(direct_entry.route_core_sidecar.lowering_source) ==
          :direct_identity_cpb
    @test direct_entry.route_core_sidecar.intermediate_retained_space.retained_rule ==
          :identity_source_modes

    outer_mismatch_entry = only(
        entry for entry in sidecar_inventory.supported_entries
        if entry.unit_key == :outer_mismatch_shared_molecular_shell
    )
    outer_mismatch_sidecar = outer_mismatch_entry.route_core_sidecar
    @test outer_mismatch_sidecar.sidecar_source ==
          :atom_growth_outer_mismatch_boundary_slab_set_plan_unit
    @test CRC.lowering_recipe(outer_mismatch_sidecar.lowering_source) ==
          :direct_boundary_slab_set
    outer_mismatch_sources =
        CRC.source_cpbs(outer_mismatch_sidecar.lowering_source)
    @test length(outer_mismatch_sources) == 2
    @test all(
        cpb -> cpb isa CRC.CoordinateProductBox,
        outer_mismatch_sources,
    )
    @test sum(CRC.support_count, outer_mismatch_sources; init = 0) ==
          CRC.support_count(outer_mismatch_sidecar.shellification_region)
    outer_mismatch_support =
        CRC.owned_support(outer_mismatch_sidecar.shellification_region)
    @test isempty(outer_mismatch_support.cpbs)
    @test !isnothing(outer_mismatch_support.outer_box)
    @test !isnothing(outer_mismatch_support.inner_exclusion_box)
    @test outer_mismatch_sidecar.intermediate_retained_space.retained_rule ==
          :direct_slab_set_identity_modes
    @test outer_mismatch_sidecar.intermediate_retained_space.dimension ==
          sum(CRC.support_count, outer_mismatch_sources; init = 0)
    @test !outer_mismatch_sidecar.numerical_behavior_changed
    @test !outer_mismatch_sidecar.materialization_behavior_changed

    @test units.route_core_sidecar_inventory_available
    @test units.route_core_sidecar_inventory ===
          units.low_order_units.route_core_sidecar_inventory
    @test units.route_core_sidecar_inventory.status ==
          sidecar_inventory.status
    @test units.route_core_sidecar_inventory.supported_unit_count ==
          sidecar_inventory.supported_unit_count
end

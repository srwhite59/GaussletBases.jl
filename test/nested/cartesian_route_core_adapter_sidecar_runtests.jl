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
    return units.plan_unit_inventory
end

@testset "CartesianRouteCore staged adapter sidecars" begin
    plan_inventory = _cartesian_route_core_adapter_fixture()
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
end

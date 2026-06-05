using Test
using GaussletBases

function _cartesian_atom_growth_checkpoint_fixture()
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
        route_kind = :be2_cartesian_atom_growth_staged_materialized_checkpoint,
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
    transforms = GaussletBases.cartesian_transforms(units, recipe)
    pairs = GaussletBases.cartesian_pair_terms(units, transforms, recipe)
    assembly =
        GaussletBases.cartesian_assembly(parent, shells, units, transforms, pairs, recipe)
    report = GaussletBases.cartesian_report(system, parent, assembly, recipe)
    materialization = GaussletBases._pqs_source_box_route_driver_materialization(
        report;
        materialize_route = true,
        save_basis_artifact = false,
        save_ham_artifact = false,
        low_order_shellization_policy = :atom_growth_complete_rectangular,
        white_lindsey_expansion =
            GaussletBases.coulomb_gaussian_expansion(doacc = false),
    )

    return (; shells, units, transforms, report, materialization)
end

function _cartesian_atom_growth_role_counts(roles)
    return Dict(role => count(==(role), roles) for role in unique(roles))
end

@testset "atom-growth staged facts match materialized checkpoint" begin
    fixture = _cartesian_atom_growth_checkpoint_fixture()
    shells = fixture.shells
    units = fixture.units
    transforms = fixture.transforms
    materialization = fixture.materialization
    materialized_report = materialization.materialized_report

    @test materialization.status ==
          :materialized_route_configured_diatomic_atom_growth_report_available
    @test materialized_report.object_kind ==
          :white_lindsey_low_order_route_configured_diatomic_atom_growth_report
    @test materialization.materialized_report_kind ==
          :white_lindsey_low_order_route_configured_diatomic_atom_growth_report

    shell_summary = shells.low_order_shellization
    @test shells.shellization_source ==
          :bond_aligned_diatomic_atom_growth_construction_plan
    @test materialized_report.shellization_source ==
          :bond_aligned_diatomic_atom_growth_construction_plan
    @test materialized_report.shellization_authority ==
          :bond_aligned_diatomic_atom_growth_construction_plan
    @test shells.atom_growth_plan_authority
    @test shell_summary.atom_growth_plan_authority
    @test materialization.route_configured_diatomic_atom_growth_plan_authority
    @test !shells.active_source_authority
    @test !materialized_report.active_source_authority
    @test !materialization.route_configured_diatomic_atom_growth_active_source_authority
    @test shells.coverage_complete ==
          materialization.route_configured_diatomic_atom_growth_coverage_complete
    @test materialized_report.support_count ==
          materialization.route_configured_diatomic_atom_growth_support_count
    @test shell_summary.atom_growth_scaffold.coverage.expected_support_count ==
          materialized_report.support_count
    @test shell_summary.atom_growth_region_count ==
          shell_summary.atom_growth_scaffold.region_count

    plan_inventory = units.plan_unit_inventory
    route_units = materialized_report.route_units
    materialized_units = route_units.retained_units
    @test plan_inventory.status == :available_atom_growth_plan_unit_inventory
    @test route_units.status == :available_atom_growth_retained_unit_inventory
    @test plan_inventory.unit_count == length(materialized_units)
    @test shell_summary.atom_growth_region_count == plan_inventory.unit_count
    @test Set(plan_inventory.unit_keys) ==
          Set(unit.unit_key for unit in materialized_units)
    @test _cartesian_atom_growth_role_counts(plan_inventory.unit_roles) ==
          _cartesian_atom_growth_role_counts(
              Tuple(unit.unit_role for unit in materialized_units),
          )

    staged_support_counts =
        Dict(unit.unit_key => unit.support_count for unit in plan_inventory.plan_units)
    materialized_support_counts =
        Dict(unit.unit_key => unit.source_dimension for unit in materialized_units)
    @test staged_support_counts == materialized_support_counts
    @test all(unit -> unit.retained_count === nothing, plan_inventory.plan_units)
    @test all(unit -> unit.retained_range === nothing, plan_inventory.plan_units)
    @test all(unit -> !isnothing(unit.retained_range), materialized_units)
    @test route_units.retained_dimension == materialization.retained_dimension

    contract_inventory =
        transforms.low_order_transforms.transform_contract_inventory
    materialized_transform_inventory = materialized_report.transform_inventory
    @test contract_inventory.status ==
          :available_atom_growth_transform_contract_inventory
    @test Set(contract_inventory.unit_keys) ==
          Set(unit.unit_key for unit in materialized_units)
    @test _cartesian_atom_growth_role_counts(contract_inventory.unit_roles) ==
          _cartesian_atom_growth_role_counts(
              Tuple(unit.unit_role for unit in materialized_units),
          )
    @test all(
        contract -> !contract.coefficient_transform_materialized,
        contract_inventory.transform_contracts,
    )
    @test all(
        contract -> !contract.coefficient_map_materialized,
        contract_inventory.transform_contracts,
    )
    @test materialized_transform_inventory.status ==
          :available_atom_growth_transform_inventory
    @test materialized_transform_inventory.retained_dimension ==
          materialization.retained_dimension
    @test materialized_transform_inventory.final_integral_weights_ready
    @test materialized_transform_inventory.final_integral_weight_count ==
          materialization.retained_dimension

    @test materialization.route_configured_materializer_backend_requested ==
          :pgdg_localized_experimental
    @test materialization.route_configured_materializer_backend_source ==
          :system_map_backend
    @test materialization.route_configured_materializer_backend_status ==
          :defaulted_from_system_map_backend
    @test materialization.route_configured_materializer_backend_consumed ==
          materialization.route_configured_materializer_backend_requested
    @test materialization.route_configured_materializer_d_requested == 0.15
    @test materialization.route_configured_materializer_d_source ==
          :recipe_core_spacing
    @test materialization.route_configured_materializer_d_status ==
          :defaulted_from_route_core_spacing
    @test materialization.route_configured_materializer_d_consumed ==
          materialization.route_configured_materializer_d_requested
    @test materialization.route_configured_materializer_nside_requested == 5
    @test materialization.route_configured_materializer_nside_source == :recipe_n_s
    @test materialization.route_configured_materializer_nside_status ==
          :defaulted_from_route_n_s
    @test materialization.route_configured_materializer_nside_consumed ==
          materialization.route_configured_materializer_nside_requested
end

# Integration/slow test. Do not include in default nested runner.

using Test
using GaussletBases

function _terminal_shellification_plan_for_boundary_test()
    policy = GaussletBases.CartesianShellification.AtomOutwardShellification(;
        core_side = 5,
        q = 5,
        bond_axis = :auto,
    )
    return GaussletBases.CartesianShellification.shellify(
        (collect(1:9), collect(1:9), collect(1:9)),
        (5, 5, 5),
        policy,
    )
end

function _terminal_unit_inventory_for_boundary_test(shellification_plan)
    raw_plan = GaussletBases.CartesianShellification.raw_plan(shellification_plan)
    records = Tuple(
        GaussletBases._cartesian_terminal_region_unit_record(region)
        for region in raw_plan.regions
    )
    return (;
        object_kind = :cartesian_terminal_region_unit_inventory,
        status = :available_terminal_region_unit_inventory,
        unit_count = length(records),
        terminal_region_units = records,
        unit_keys = Tuple(record.unit_key for record in records),
        unit_roles = Tuple(record.unit_role for record in records),
        unit_kinds = Tuple(record.unit_kind for record in records),
        terminal_region_roles =
            Tuple(record.terminal_region_role for record in records),
        terminal_region_kinds =
            Tuple(record.terminal_region_kind for record in records),
        support_counts = Tuple(record.support_count for record in records),
        complete_shell_unit_count =
            count(record -> record.terminal_region_kind == :complete_shell, records),
    )
end

@testset "driver helper boundary uses terminal lowering module summary" begin
    shellification_plan = _terminal_shellification_plan_for_boundary_test()
    low_order_shellization = (;
        route_family = :white_lindsey_low_order,
        terminal_shellification_plan = shellification_plan,
    )
    lowering_plan =
        GaussletBases._pqs_source_box_route_driver_terminal_lowering_plan(
            low_order_shellization,
            :white_lindsey_low_order,
        )
    unit_inventory = _terminal_unit_inventory_for_boundary_test(shellification_plan)
    lowering_inventory =
        GaussletBases._pqs_source_box_route_driver_terminal_lowering_contract_inventory_from_plan(
            lowering_plan,
            unit_inventory,
        )
    selected_inventory =
        GaussletBases._pqs_source_box_route_driver_selected_terminal_lowering_contract_inventory_from_plan(
            lowering_plan,
            lowering_inventory,
            :white_lindsey_low_order,
        )

    @test shellification_plan isa
          GaussletBases.CartesianShellification.ShellificationPlan
    @test lowering_plan isa
          GaussletBases.CartesianTerminalLowering.TerminalLoweringPlan
    @test lowering_inventory.inventory_source ==
          :terminal_lowering_plan_compatibility_adapter
    @test selected_inventory.inventory_source ==
          :terminal_lowering_plan_compatibility_adapter

    module_summary = GaussletBases.CartesianTerminalLowering.summary(
        lowering_plan,
    )
    selected_fields =
        GaussletBases._pqs_source_box_route_driver_selected_terminal_lowering_fields(
            nothing,
            :available_terminal_lowering_plan,
            :white_lindsey_low_order,
            lowering_plan,
        )
    selected_inventory_fields =
        GaussletBases._pqs_source_box_route_driver_selected_terminal_lowering_fields(
            selected_inventory,
            selected_inventory.status,
            :white_lindsey_low_order,
        )
    @test module_summary.status == :available_terminal_lowering_plan
    @test module_summary.policy_kind == :white_lindsey_terminal_lowering
    @test module_summary.terminal_region_count ==
          length(GaussletBases.CartesianShellification.terminal_regions(shellification_plan))
    @test module_summary.available_contract_count >
          module_summary.selected_contract_count
    @test module_summary.selected_contract_count ==
          module_summary.terminal_region_count
    @test :white_lindsey_boundary_strata in
          module_summary.selected_contract_kinds
    @test :pqs_filled_source_cpb in module_summary.available_contract_kinds
    @test selected_fields.terminal_shellification_selected_contract_count ==
          module_summary.selected_contract_count
    @test selected_fields.terminal_shellification_selected_contract_kinds ==
          module_summary.selected_contract_kinds
    @test selected_fields.terminal_shellification_unselected_contract_count ==
          module_summary.available_contract_count -
          module_summary.selected_contract_count
    @test lowering_inventory.lowering_contract_count ==
          module_summary.available_contract_count
    @test selected_inventory.selected_contract_count ==
          module_summary.selected_contract_count

    crc_sidecar_summary =
        GaussletBases._pqs_source_box_route_driver_selected_terminal_crc_sidecar_summary(
            selected_inventory,
        )
    terminal_route_state =
        GaussletBases._pqs_source_box_route_driver_terminal_route_state(;
            status = module_summary.status,
            selected = true,
            route_lowering_family = :white_lindsey_low_order,
            shellification_plan,
            unit_inventory,
            lowering_plan,
            lowering_summary = module_summary,
            lowering_contract_inventory = lowering_inventory,
            selected_contract_inventory = selected_inventory,
            selected_crc_sidecar_summary = crc_sidecar_summary,
        )
    terminal_route_summary = terminal_route_state.summary

    @test terminal_route_state.object_kind ==
          :cartesian_driver_terminal_route_state
    @test terminal_route_summary.object_kind ==
          :cartesian_driver_terminal_route_state_summary
    @test terminal_route_summary.status == module_summary.status
    @test terminal_route_summary.unit_count == unit_inventory.unit_count
    @test terminal_route_summary.available_contract_count ==
          module_summary.available_contract_count
    @test terminal_route_summary.selected_contract_count ==
          module_summary.selected_contract_count
    @test terminal_route_summary.selected_contract_kinds ==
          module_summary.selected_contract_kinds
    @test terminal_route_summary.selected_crc_sidecar_status ==
          crc_sidecar_summary.status
    @test terminal_route_summary.selected_crc_sidecar_available_count ==
          crc_sidecar_summary.sidecar_available_count
    @test terminal_route_summary.selected_crc_sidecar_missing_count ==
          crc_sidecar_summary.sidecar_missing_count
    @test !terminal_route_summary.operator_blocks_materialized

    alias_source = merge(
        (;
            terminal_shellification_scaffold_available = true,
            terminal_shellification_scaffold = :synthetic_scaffold,
            terminal_shellification_region_count =
                module_summary.terminal_region_count,
            terminal_shellification_unit_inventory_available = true,
            terminal_shellification_unit_inventory = unit_inventory,
            terminal_shellification_unit_count = unit_inventory.unit_count,
            terminal_shellification_unit_keys = unit_inventory.unit_keys,
            terminal_shellification_unit_roles = unit_inventory.unit_roles,
            terminal_shellification_unit_kinds = unit_inventory.unit_kinds,
            terminal_shellification_unit_support_counts =
                unit_inventory.support_counts,
            terminal_shellification_lowering_contract_inventory_available = true,
            terminal_shellification_lowering_contract_inventory_status =
                lowering_inventory.status,
            terminal_shellification_lowering_contract_inventory =
                lowering_inventory,
            terminal_shellification_lowering_contract_count =
                lowering_inventory.lowering_contract_count,
            terminal_shellification_lowering_contract_kinds =
                lowering_inventory.lowering_contract_kinds,
            terminal_shellification_lowering_contract_kind_counts =
                lowering_inventory.lowering_contract_kind_counts,
            terminal_shellification_selected_crc_sidecar_summary =
                crc_sidecar_summary,
            terminal_shellification_contract_counts_by_unit =
                lowering_inventory.contract_counts_by_unit,
            terminal_shellification_lw_complete_shell_cpb_count =
                lowering_inventory.lw_complete_shell_cpb_count,
            terminal_shellification_lw_complete_shell_cpb_family_counts =
                lowering_inventory.lw_complete_shell_cpb_family_counts,
            terminal_shellification_final_retained_unit_inventory_available =
                lowering_inventory.final_retained_unit_inventory_available,
            terminal_shellification_central_gap_region_count = 0,
            terminal_shellification_central_midpoint_slab_count = 0,
            terminal_shellification_central_distorted_product_box_count = 0,
            terminal_shellification_central_distorted_product_box_metadata = (),
        ),
        selected_inventory_fields,
    )
    alias_fields =
        GaussletBases._pqs_source_box_route_driver_terminal_shellification_alias_fields(
            alias_source,
            true,
        )
    unselected_alias_fields =
        GaussletBases._pqs_source_box_route_driver_terminal_shellification_alias_fields(
            alias_source,
            false;
            include_crc_sidecar_summary = false,
        )

    @test alias_fields.terminal_shellification_unit_count ==
          unit_inventory.unit_count
    @test alias_fields.terminal_shellification_lowering_contract_count ==
          module_summary.available_contract_count
    @test alias_fields.terminal_shellification_selected_contract_count ==
          module_summary.selected_contract_count
    @test alias_fields.terminal_shellification_selected_crc_sidecar_summary.status ==
          crc_sidecar_summary.status
    @test !alias_fields.terminal_shellification_final_retained_unit_inventory_available
    @test unselected_alias_fields.terminal_shellification_unit_count == 0
    @test unselected_alias_fields.terminal_shellification_lowering_contract_kinds ==
          ()
    @test !hasproperty(
        unselected_alias_fields,
        :terminal_shellification_selected_crc_sidecar_summary,
    )

    @test !module_summary.materialized
    @test !module_summary.operator_blocks_materialized
    @test !module_summary.hamiltonian_data_materialized
end

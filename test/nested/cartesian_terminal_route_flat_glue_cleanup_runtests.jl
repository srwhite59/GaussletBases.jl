using Test
using GaussletBases

function _flat_glue_cleanup_low_order_units_fixture()
    terminal_route_state =
        GaussletBases._pqs_source_box_route_driver_terminal_route_state(;
            status = :not_selected,
            selected = false,
            route_lowering_family = nothing,
        )
    return (;
        status = :available_unit_stage_low_order_summary,
        terminal_route_state,
        terminal_route_summary = terminal_route_state.summary,
        low_order_shellization_policy_requested = nothing,
        low_order_shellization_policy_resolved = :not_selected,
        low_order_shellization_policy_source = :not_selected,
        low_order_shellization_policy_status = :not_selected,
        low_order_shellization_policy_blocker = nothing,
        shellization_source = :not_selected,
        shellization_kind = :not_selected,
        unit_route_kind = :not_selected,
        atom_growth_units_selected = false,
        terminal_shellification_units_selected = false,
        legacy_source_units_selected = false,
        materialization_status = :not_selected,
        retained_unit_dimensions_known = false,
        retained_unit_ranges_known = false,
        retained_dimension_known = false,
        retained_dimension = nothing,
        plan_authority = false,
        active_source_authority = false,
        legacy_source_authority = false,
    )
end

@testset "terminal route flat glue transform summary fingerprint" begin
    missing_summary =
        GaussletBases._pqs_source_box_route_driver_transform_stage_low_order_summary((;))
    @test missing_summary.status == :not_available_missing_unit_stage_summary
    @test hasproperty(missing_summary, :terminal_route_summary)
    @test !hasproperty(missing_summary, :terminal_shellification_pair_inventory_available)
    @test !missing_summary.terminal_route_summary.pair_inventory_available
    @test missing_summary.terminal_route_summary.unit_pair_summary.status == :not_selected

    summary =
        GaussletBases._pqs_source_box_route_driver_transform_stage_low_order_summary(
            (; low_order_units = _flat_glue_cleanup_low_order_units_fixture()),
        )
    @test summary.status == :available_transform_stage_low_order_summary
    @test hasproperty(summary, :terminal_route_summary)
    @test !hasproperty(summary, :terminal_shellification_pair_inventory_available)
    @test !summary.terminal_route_summary.pair_inventory_available
    @test summary.terminal_route_summary.unit_pair_summary.status == :not_selected
end

# Runtime role: decomposed White-Lindsey H/H2+ acceptance readiness audit.
#
# The former post-CPB acceptance test used one CPB covering the full parent
# window. That is not a decomposed White-Lindsey acceptance path. Keep this file
# as a small blocker audit until q = 5, ns = 5 scientific H/H2+ acceptance
# assembly consumes real decomposed retained-unit and unit-pair inventories with
# retained column ranges.

using Test
using GaussletBases

const WLAcceptanceReadinessCPBM = GaussletBases.CartesianPairBlockMaterialization

function _wl_decomposed_acceptance_blocker_report()
    adapter = WLAcceptanceReadinessCPBM.white_lindsey_boundary_stratum_one_body_adapter_summary()
    local_terms = adapter.supported_one_body_terms
    route_global_terms = WLAcceptanceReadinessCPBM.route_global_safe_one_body_terms()
    collection = (;
        object_kind = :cartesian_pair_block_local_one_body_block_collection,
        terms = (),
        requested_terms = (),
        materialized_terms = (),
        deferred_terms = (),
        entries = (),
        materialized_entries = (),
        skipped_entries = (),
    )

    placement_plan_blocked = false
    placement_plan_error = nothing
    try
        WLAcceptanceReadinessCPBM.one_body_placement_plan(
            collection;
            term = :electron_nuclear_by_center,
            global_dimension = 1,
        )
    catch error
        placement_plan_blocked = true
        placement_plan_error = error
    end

    return (;
        object_kind = :decomposed_wl_h_h2plus_acceptance_readiness_audit,
        status = :blocked_decomposed_wl_h_h2plus_acceptance,
        blocker = :missing_decomposed_wl_unit_pair_inventory_with_column_ranges,
        q = 5,
        ns = 5,
        n_s = 5,
        decomposed_wl_units_consumed = false,
        full_parent_window_cpb_used = false,
        direct_cartesian_product_assembly_used = false,
        ordinary_cartesian_ida_operators_used = false,
        supported_decomposed_one_body_terms = local_terms,
        route_global_safe_one_body_terms = route_global_terms,
        decomposed_overlap_available = :overlap in local_terms,
        decomposed_kinetic_available = :kinetic in local_terms,
        decomposed_electron_nuclear_by_center_available =
            :electron_nuclear_by_center in local_terms,
        decomposed_electron_nuclear_by_center_selector_available =
            :electron_nuclear_by_center in local_terms,
        route_global_overlap_adapter_available = :overlap in route_global_terms,
        route_global_kinetic_adapter_available = :kinetic in route_global_terms,
        route_global_electron_nuclear_by_center_adapter_available =
            :electron_nuclear_by_center in route_global_terms,
        decomposed_electron_nuclear_by_center_placement_plan_available =
            !placement_plan_blocked,
        decomposed_electron_nuclear_by_center_global_matrix_available =
            !placement_plan_blocked,
        placement_plan_error_type =
            isnothing(placement_plan_error) ? nothing : typeof(placement_plan_error),
        unit_inventory_audit_source =
            :terminal_cartesian_shellification_geometry_route_summary,
        terminal_shellification_unit_inventory_exposed = true,
        terminal_shellification_unit_inventory_granularity =
            :terminal_region_units,
        terminal_shellification_pair_inventory_exposed = false,
        terminal_shellification_pair_inventory_status =
            :deferred_terminal_shellification_pair_inventory,
        terminal_shellification_pair_materialization_status =
            :deferred_terminal_shellification_pair_materialization,
        retained_unit_column_ranges_materialized = false,
        retained_dimension_from_decomposed_unit_inventory_available = false,
        decomposed_unit_pair_column_ranges_available = false,
        route_global_by_center_acceptance_matrix_available = false,
        fixed_block_operator_matrices_available =
            :overlap in route_global_terms && :kinetic in route_global_terms,
        fixed_block_operator_matrices_used = false,
        fixed_block_operator_matrix_source_rejected =
            :nested_fixed_block_is_not_decomposed_wl_acceptance_path,
        acceptance_energy_materialized = false,
        h_atom_acceptance_active = false,
        h2plus_acceptance_active = false,
    )
end

@testset "decomposed WL gausslet-only H/H2+ acceptance readiness" begin
    report = _wl_decomposed_acceptance_blocker_report()
    println("decomposed WL gausslet-only H/H2+ acceptance readiness: ", report)

    @test report.status == :blocked_decomposed_wl_h_h2plus_acceptance
    @test report.blocker ==
          :missing_decomposed_wl_unit_pair_inventory_with_column_ranges
    @test report.q == 5
    @test report.ns == 5
    @test report.n_s == 5
    @test report.decomposed_overlap_available
    @test report.decomposed_kinetic_available
    @test report.decomposed_electron_nuclear_by_center_available
    @test report.decomposed_electron_nuclear_by_center_selector_available
    @test report.route_global_overlap_adapter_available
    @test report.route_global_kinetic_adapter_available
    @test !report.route_global_electron_nuclear_by_center_adapter_available
    @test report.decomposed_electron_nuclear_by_center_placement_plan_available
    @test report.decomposed_electron_nuclear_by_center_global_matrix_available
    @test report.terminal_shellification_unit_inventory_exposed
    @test report.terminal_shellification_unit_inventory_granularity ==
          :terminal_region_units
    @test !report.terminal_shellification_pair_inventory_exposed
    @test report.terminal_shellification_pair_inventory_status ==
          :deferred_terminal_shellification_pair_inventory
    @test !report.retained_unit_column_ranges_materialized
    @test !report.retained_dimension_from_decomposed_unit_inventory_available
    @test !report.decomposed_unit_pair_column_ranges_available
    @test !report.route_global_by_center_acceptance_matrix_available
    @test report.fixed_block_operator_matrices_available
    @test !report.fixed_block_operator_matrices_used
    @test !report.decomposed_wl_units_consumed
    @test !report.full_parent_window_cpb_used
    @test !report.direct_cartesian_product_assembly_used
    @test !report.ordinary_cartesian_ida_operators_used
    @test !report.acceptance_energy_materialized
    @test !report.h_atom_acceptance_active
    @test !report.h2plus_acceptance_active
end

using Test

include("cartesian_pair_stage_fingerprint_helpers.jl")

@testset "cartesian pair-stage compact fingerprint helpers" begin
    synthetic = (;
        terminal_shellification_unit_inventory_available = true,
        terminal_shellification_region_count = 2,
        terminal_shellification_unit_count = 2,
        terminal_shellification_unit_keys = (:left, :right),
        terminal_shellification_unit_roles = (:atom_local, :shared_shell),
        terminal_shellification_unit_kinds = (:direct_core, :complete_shell),
        terminal_shellification_unit_support_counts = (8, 26),
        terminal_shellification_central_gap_region_count = 1,
        terminal_shellification_central_midpoint_slab_count = 1,
        terminal_shellification_central_distorted_product_box_count = 0,
        terminal_shellification_lowering_contract_inventory_available = true,
        terminal_shellification_lowering_contract_inventory_status =
            :available,
        terminal_shellification_lowering_contract_count = 3,
        terminal_shellification_lowering_contract_kinds =
            (:direct_core_identity_cpb, :white_lindsey_boundary_strata),
        terminal_shellification_lowering_contract_kind_counts = (;
            direct_core_identity_cpb_count = 1,
            white_lindsey_boundary_strata_count = 2,
        ),
        terminal_shellification_contract_counts_by_unit = (
            (; unit_key = :left, lowering_contract_count = 1),
        ),
        terminal_shellification_lw_complete_shell_cpb_count = 26,
        terminal_shellification_lw_complete_shell_cpb_family_counts = (;
            facet_cpb_count = 6,
            edge_cpb_count = 12,
            corner_cpb_count = 8,
        ),
        terminal_shellification_selected_lowering_contract_inventory_available =
            true,
        terminal_shellification_selected_lowering_contract_inventory_status =
            :available,
        terminal_shellification_selected_lowering_family =
            :white_lindsey_low_order,
        terminal_shellification_selected_contract_count = 2,
        terminal_shellification_selected_contract_kinds =
            (:direct_core_identity_cpb, :white_lindsey_boundary_strata),
        terminal_shellification_selected_contract_kind_counts = (;
            direct_core_identity_cpb_count = 1,
            white_lindsey_boundary_strata_count = 1,
        ),
        terminal_shellification_selected_contract_counts_by_unit =
            ((; unit_key = :right, selected_contract_count = 1),),
        terminal_shellification_all_units_have_exactly_one_selected_contract =
            true,
        terminal_shellification_unselected_contract_count = 1,
        terminal_shellification_unselected_contract_kinds =
            (:pqs_filled_source_cpb,),
        object_kind =
            :cartesian_unit_stage_selected_terminal_lowering_crc_sidecar_summary,
        status = :available,
        selected_contract_count = 2,
        sidecar_available_count = 2,
        sidecar_missing_count = 0,
        sidecar_inventory_complete = true,
        missing_sidecar_reasons = (),
        missing_sidecar_kinds = (),
        final_retained_unit_inventory_available = false,
        pair_inventory_available = false,
        pair_inventory_status = :deferred,
        operator_blocks_materialized = false,
        pair_operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
        pair_inventory_known = true,
        pair_inventory_source = :final_retained_units,
        pair_entries = ((; pair_key = (:left, :left)),),
        pair_count = 1,
        pair_family_counts = (; direct_pair_count = 1),
        helper_by_pair_family = (; direct_pair = :metadata_only),
        pair_operator_helper_by_family = (; direct_pair = :metadata_only),
        pair_helper_status_by_family = (; direct_pair = :metadata_only),
        operator_pairs_materialized = false,
        route_core_pair_inventory_available = true,
        route_core_pair_inventory_status = :available,
        route_core_pair_count = 1,
        route_core_pair_keys = ((:left, :left),),
        route_core_pair_operator_ready = true,
        route_core_pair_operator_readiness_status = :ready,
        route_core_pair_operator_blocker = nothing,
    )

    unit_fingerprint =
        _cartesian_pair_stage_terminal_unit_inventory_fingerprint(synthetic)
    @test unit_fingerprint.unit_count == 2
    @test unit_fingerprint.unit_keys == (:left, :right)

    lowering_fingerprint =
        _cartesian_pair_stage_lowering_contract_fingerprint(synthetic)
    @test lowering_fingerprint.count == 3
    @test lowering_fingerprint.counts_by_unit == (
        (;
            unit_key = :left,
            lowering_contract_count = 1,
            selected_contract_count = nothing,
        ),
    )

    selected_fingerprint =
        _cartesian_pair_stage_selected_contract_fingerprint(synthetic)
    @test selected_fingerprint.family == :white_lindsey_low_order
    @test selected_fingerprint.unselected_kinds == (:pqs_filled_source_cpb,)
    @test selected_fingerprint.all_units_have_exactly_one_selected_contract

    crc_fingerprint =
        _cartesian_pair_stage_crc_sidecar_summary_fingerprint(synthetic)
    @test crc_fingerprint.sidecar_inventory_complete
    @test !crc_fingerprint.pair_inventory_available

    pair_fingerprint =
        _cartesian_pair_stage_pair_metadata_fingerprint(synthetic)
    @test pair_fingerprint.pair_entry_count == 1
    @test pair_fingerprint.route_core_pair_key_count == 1
    @test pair_fingerprint.route_core_pair_operator_ready
    @test !pair_fingerprint.pair_operator_blocks_materialized
end

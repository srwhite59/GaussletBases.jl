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

function _flat_glue_cleanup_alias_source_without_crc_field()
    terminal_route_state =
        GaussletBases._pqs_source_box_route_driver_terminal_route_state(;
            status = :not_available,
            selected = true,
            route_lowering_family = nothing,
        )
    empty_kind_counts = (
        direct_core_identity_cpb_count = 0,
        direct_slab_identity_cpb_count = 0,
        direct_boundary_slab_identity_cpb_count = 0,
        white_lindsey_boundary_strata_count = 0,
        pqs_filled_source_cpb_count = 0,
        distorted_product_box_comx_count = 0,
    )
    return (;
        terminal_route_state,
        terminal_shellification_selected_lowering_contract_inventory = nothing,
        terminal_shellification_selected_lowering_contract_inventory_status =
            :not_available,
        terminal_shellification_selected_lowering_family = nothing,
        terminal_shellification_scaffold_available = true,
        terminal_shellification_scaffold = nothing,
        terminal_shellification_region_count = 0,
        terminal_shellification_unit_inventory_available = false,
        terminal_shellification_unit_inventory = nothing,
        terminal_shellification_unit_count = 0,
        terminal_shellification_unit_keys = (),
        terminal_shellification_unit_roles = (),
        terminal_shellification_unit_kinds = (),
        terminal_shellification_unit_support_counts = (),
        terminal_shellification_lowering_contract_inventory_available = false,
        terminal_shellification_lowering_contract_inventory_status = :not_available,
        terminal_shellification_lowering_contract_inventory = nothing,
        terminal_shellification_lowering_contract_count = 0,
        terminal_shellification_lowering_contract_kinds = (),
        terminal_shellification_lowering_contract_kind_counts = empty_kind_counts,
        terminal_shellification_contract_counts_by_unit = (),
        terminal_shellification_lw_complete_shell_cpb_count = 0,
        terminal_shellification_lw_complete_shell_cpb_family_counts =
            (facet_cpb = 0, edge_cpb = 0, corner_cpb = 0),
        terminal_shellification_final_retained_unit_inventory_available = false,
        terminal_shellification_central_gap_region_count = 0,
        terminal_shellification_central_midpoint_slab_count = 0,
        terminal_shellification_central_distorted_product_box_count = 0,
        terminal_shellification_central_distorted_product_box_metadata = (),
    )
end

function _flat_glue_cleanup_route_skeleton_fixture()
    return (;
        pair_entries = (),
        pair_family_counts = nothing,
        helper_by_pair_family = nothing,
    )
end

@testset "terminal route flat glue unit summary fingerprint" begin
    missing_summary =
        GaussletBases._pqs_source_box_route_driver_unit_stage_low_order_summary((;))
    @test missing_summary.status == :not_available_missing_shell_stage_summary
    @test hasproperty(missing_summary, :terminal_route_summary)
    @test !hasproperty(
        missing_summary,
        :terminal_shellification_selected_crc_sidecar_summary,
    )
    @test missing_summary.terminal_route_summary.selected_crc_sidecar_status ==
          :not_available_selected_terminal_lowering_contract_inventory
end

@testset "terminal route flat glue transform summary fingerprint" begin
    missing_summary =
        GaussletBases._pqs_source_box_route_driver_transform_stage_low_order_summary((;))
    @test missing_summary.status == :not_available_missing_unit_stage_summary
    @test hasproperty(missing_summary, :terminal_route_summary)
    @test !hasproperty(missing_summary, :terminal_shellification_pair_inventory_available)
    @test !hasproperty(
        missing_summary,
        :terminal_shellification_selected_crc_sidecar_summary,
    )
    @test !missing_summary.terminal_route_summary.pair_inventory_available
    @test missing_summary.terminal_route_summary.unit_pair_summary.status == :not_selected

    summary =
        GaussletBases._pqs_source_box_route_driver_transform_stage_low_order_summary(
            (; low_order_units = _flat_glue_cleanup_low_order_units_fixture()),
        )
    @test summary.status == :available_transform_stage_low_order_summary
    @test hasproperty(summary, :terminal_route_summary)
    @test !hasproperty(summary, :terminal_shellification_pair_inventory_available)
    @test !hasproperty(
        summary,
        :terminal_shellification_selected_crc_sidecar_summary,
    )
    @test !summary.terminal_route_summary.pair_inventory_available
    @test summary.terminal_route_summary.unit_pair_summary.status == :not_selected

    alias_source = _flat_glue_cleanup_alias_source_without_crc_field()
    alias_fields =
        GaussletBases._pqs_source_box_route_driver_terminal_shellification_alias_fields(
            alias_source,
            true,
        )
    @test hasproperty(
        alias_fields,
        :terminal_shellification_selected_crc_sidecar_summary,
    )
    @test alias_fields.terminal_shellification_selected_crc_sidecar_summary.status ==
          alias_source.terminal_route_state.selected_crc_sidecar_summary.status
end

@testset "terminal route flat glue pair summary fingerprint" begin
    missing_summary =
        GaussletBases._pqs_source_box_route_driver_pair_stage_low_order_summary(
            (;),
            (;),
            _flat_glue_cleanup_route_skeleton_fixture(),
        )
    @test missing_summary.status == :not_available_missing_transform_stage_summary
    @test hasproperty(missing_summary, :terminal_route_summary)
    @test !hasproperty(
        missing_summary,
        :terminal_shellification_selected_crc_sidecar_summary,
    )
    @test !missing_summary.terminal_route_summary.pair_inventory_available

    transform_summary =
        GaussletBases._pqs_source_box_route_driver_transform_stage_low_order_summary(
            (; low_order_units = _flat_glue_cleanup_low_order_units_fixture()),
        )
    summary =
        GaussletBases._pqs_source_box_route_driver_pair_stage_low_order_summary(
            (;),
            (; low_order_transforms = transform_summary),
            _flat_glue_cleanup_route_skeleton_fixture(),
        )
    @test summary.status == :available_pair_stage_low_order_summary
    @test hasproperty(summary, :terminal_route_summary)
    @test !hasproperty(
        summary,
        :terminal_shellification_selected_crc_sidecar_summary,
    )
    @test !summary.terminal_route_summary.pair_inventory_available
end

@testset "pair-operator report aliases prefer compact summary" begin
    compact_summary = (;
        object_kind = :cartesian_pair_operator_plan_summary,
        route_core_pair_operator_plan_inventory_available = true,
        route_core_pair_operator_plan_inventory_status =
            :available_route_core_pair_operator_plan_inventory,
        route_core_pair_operator_plan_inventory_blocker = nothing,
        route_core_pair_operator_plan_count = 3,
        route_core_pair_operator_plan_blocked_count = 1,
        materialized = false,
        source_operator_path_counts = (
            (; source_operator_path = :direct_identity_cpb_path, count = 1),
            (; source_operator_path = :pqs_source_cpb_1d_factor_path, count = 2),
        ),
        final_block_path_counts = (
            (; final_block_path = :source_block_direct_to_final_block, count = 1),
            (;
                final_block_path =
                    :source_block_realization_then_final_block,
                count = 2,
            ),
        ),
        materialization_status_counts = (
            (; materialization_status = :metadata_only_not_materialized, count = 2),
            (;
                materialization_status =
                    :blocked_metadata_only_not_materialized,
                count = 1,
            ),
        ),
        blocker_counts = (
            (; blocker = nothing, count = 2),
            (; blocker = :aggregate_subtree_operator_plan_required, count = 1),
        ),
    )
    fallback = (;
        route_core_typed_pair_operator_plan_inventory_available = false,
        route_core_typed_pair_operator_plan_inventory_status =
            :blocked_route_core_pair_operator_plan_inventory,
        route_core_typed_pair_operator_plan_blocker = :fallback_blocker,
        route_core_typed_pair_operator_plan_count = 0,
        route_core_typed_pair_operator_plan_blocked_count = 0,
        route_core_typed_pair_operator_plan_materialized = true,
        route_core_typed_pair_operator_source_path_counts = (),
        route_core_typed_pair_operator_final_block_path_counts = (),
        route_core_typed_pair_operator_materialization_status_counts = (),
        route_core_typed_pair_operator_blocker_counts = (),
    )
    aliases =
        GaussletBases._pqs_source_box_route_driver_pair_operator_report_aliases(
            (; pair_operator_summary = compact_summary),
            fallback,
        )

    @test aliases.route_core_typed_pair_operator_plan_inventory_available
    @test aliases.route_core_typed_pair_operator_plan_inventory_status ==
          fallback.route_core_typed_pair_operator_plan_inventory_status
    @test aliases.route_core_typed_pair_operator_plan_blocker ==
          fallback.route_core_typed_pair_operator_plan_blocker
    @test aliases.route_core_typed_pair_operator_plan_count == 3
    @test aliases.route_core_typed_pair_operator_plan_blocked_count == 1
    @test !aliases.route_core_typed_pair_operator_plan_materialized
    @test aliases.route_core_typed_pair_operator_source_path_counts[1] ==
          (; source_operator_path = :direct_identity_cpb_path, pair_count = 1)
    @test aliases.route_core_typed_pair_operator_final_block_path_counts[2] ==
          (;
              final_block_path = :source_block_realization_then_final_block,
              pair_count = 2,
          )
    @test aliases.route_core_typed_pair_operator_blocker_counts[2] ==
          (;
              blocker = :aggregate_subtree_operator_plan_required,
              pair_count = 1,
          )

    fallback_aliases =
        GaussletBases._pqs_source_box_route_driver_pair_operator_report_aliases(
            (;),
            fallback,
        )
    @test fallback_aliases.route_core_typed_pair_operator_plan_inventory_status ==
          :blocked_route_core_pair_operator_plan_inventory
    @test fallback_aliases.route_core_typed_pair_operator_plan_materialized
end

using Test
using GaussletBases

function _selected_lowering_test_contract(
    contract_index::Int,
    unit_key::Symbol,
    terminal_region_kind::Symbol,
    lowering_contract_kind::Symbol;
    identity_like_source_contract::Bool,
    source_cpb_count::Int = 1,
    source_cpb_family_counts = (facet_cpb = 0, edge_cpb = 0, corner_cpb = 0),
    face_edge_corner_decomposition_required::Bool = false,
    final_unit_count_planned::Int = 1,
)
    return (;
        object_kind = :cartesian_terminal_region_lowering_contract,
        status = :planned_not_materialized,
        contract_index,
        contract_key =
            Symbol(String(unit_key), "_", String(lowering_contract_kind)),
        unit_index = contract_index,
        unit_key,
        terminal_region_kind,
        lowering_contract_kind,
        identity_like_source_contract,
        source_cpb_count,
        source_cpb_family_counts,
        face_edge_corner_decomposition_required,
        final_unit_count_planned,
        coefficient_maps_materialized = false,
        transform_contracts_materialized = false,
        retained_spaces_materialized = false,
        operator_blocks_materialized = false,
        pair_operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
        final_retained_unit_records_materialized = false,
    )
end

function _selected_lowering_test_inventory()
    lowering_contracts = (
        _selected_lowering_test_contract(
            1,
            :direct_core,
            :direct_core,
            :direct_core_identity_cpb;
            identity_like_source_contract = true,
        ),
        _selected_lowering_test_contract(
            2,
            :midpoint_slab,
            :direct_midpoint_slab,
            :direct_slab_identity_cpb;
            identity_like_source_contract = true,
        ),
        _selected_lowering_test_contract(
            3,
            :outer_mismatch,
            :outer_mismatch_slab,
            :direct_boundary_slab_identity_cpb;
            identity_like_source_contract = true,
        ),
        _selected_lowering_test_contract(
            4,
            :complete_shell,
            :complete_shell,
            :white_lindsey_boundary_strata;
            identity_like_source_contract = false,
            source_cpb_count = 26,
            source_cpb_family_counts = (facet_cpb = 6, edge_cpb = 12, corner_cpb = 8),
            face_edge_corner_decomposition_required = true,
            final_unit_count_planned = 26,
        ),
        _selected_lowering_test_contract(
            5,
            :complete_shell,
            :complete_shell,
            :pqs_filled_source_cpb;
            identity_like_source_contract = false,
        ),
        _selected_lowering_test_contract(
            6,
            :central_distorted,
            :central_distorted_product_box,
            :distorted_product_box_comx;
            identity_like_source_contract = false,
        ),
    )
    return (;
        object_kind = :cartesian_terminal_region_lowering_contract_inventory,
        status = :available_terminal_region_lowering_contract_inventory,
        terminal_region_unit_count = 5,
        lowering_contract_count = length(lowering_contracts),
        lowering_contracts,
    )
end

function _selected_lowering_contract_by_unit(inventory, unit_key::Symbol)
    return only(
        contract for contract in inventory.selected_contracts
        if contract.unit_key == unit_key
    )
end

function _assert_selected_lowering_inventory_metadata_only(inventory)
    @test inventory.object_kind ==
          :cartesian_selected_terminal_lowering_contract_inventory
    @test inventory.status ==
          :available_selected_terminal_lowering_contract_inventory
    @test inventory.inventory_source == :terminal_region_lowering_contract_inventory
    @test inventory.source_object_kind ==
          :cartesian_terminal_region_lowering_contract_inventory
    @test inventory.private_development_only
    @test inventory.terminal_region_unit_count == 5
    @test inventory.selected_contract_count == 5
    @test inventory.selected_contract_records === inventory.selected_contracts
    @test inventory.all_units_have_exactly_one_selected_contract
    @test all(
        entry -> entry.selected_contract_count == 1,
        inventory.selected_contract_counts_by_unit,
    )
    @test inventory.unselected_contract_count == 1
    @test !inventory.final_retained_unit_inventory_available
    @test !inventory.pair_inventory_available
    @test inventory.pair_inventory_status ==
          :not_available_selected_lowering_metadata_only
    @test !inventory.coefficient_maps_materialized
    @test !inventory.transform_contracts_materialized
    @test !inventory.retained_spaces_materialized
    @test !inventory.operator_blocks_materialized
    @test !inventory.pair_operator_blocks_materialized
    @test !inventory.hamiltonian_data_materialized
    @test !inventory.artifacts_materialized
    @test inventory.diagnostics.selected_lowering_metadata_only
    @test inventory.diagnostics.all_units_have_exactly_one_selected_contract
    @test !inventory.diagnostics.final_retained_unit_inventory_available
    @test !inventory.diagnostics.pair_inventory_available
    @test !inventory.diagnostics.coefficient_maps_materialized
    @test !inventory.diagnostics.transform_contracts_materialized
    @test !inventory.diagnostics.retained_spaces_materialized
    @test !inventory.diagnostics.operator_blocks_materialized
    @test !inventory.diagnostics.pair_operator_blocks_materialized
    @test !inventory.diagnostics.hamiltonian_data_materialized
    @test !inventory.diagnostics.artifacts_materialized
    @test !inventory.diagnostics.materialization_behavior_changed
    @test !inventory.diagnostics.public_default_behavior_changed
    @test all(
        contract -> !contract.coefficient_maps_materialized &&
                    !contract.transform_contracts_materialized &&
                    !contract.retained_spaces_materialized &&
                    !contract.operator_blocks_materialized &&
                    !contract.pair_operator_blocks_materialized &&
                    !contract.hamiltonian_data_materialized &&
                    !contract.artifacts_materialized &&
                    !contract.final_retained_unit_records_materialized,
        inventory.selected_contracts,
    )
end

@testset "selected terminal lowering contract inventory chooses one route contract per unit" begin
    lowering_inventory = _selected_lowering_test_inventory()

    lw_inventory =
        GaussletBases._cartesian_selected_terminal_lowering_contract_inventory(
            lowering_inventory,
            :white_lindsey_low_order,
        )
    _assert_selected_lowering_inventory_metadata_only(lw_inventory)
    @test lw_inventory.route_lowering_family == :white_lindsey_low_order
    @test _selected_lowering_contract_by_unit(lw_inventory, :direct_core).lowering_contract_kind ==
          :direct_core_identity_cpb
    @test _selected_lowering_contract_by_unit(lw_inventory, :midpoint_slab).lowering_contract_kind ==
          :direct_slab_identity_cpb
    @test _selected_lowering_contract_by_unit(lw_inventory, :outer_mismatch).lowering_contract_kind ==
          :direct_boundary_slab_identity_cpb
    lw_shell = _selected_lowering_contract_by_unit(lw_inventory, :complete_shell)
    @test lw_shell.lowering_contract_kind == :white_lindsey_boundary_strata
    @test lw_shell.face_edge_corner_decomposition_required
    @test lw_shell.source_cpb_count == 26
    @test lw_shell.source_cpb_family_counts ==
          (facet_cpb = 6, edge_cpb = 12, corner_cpb = 8)
    @test only(lw_inventory.unselected_contracts).lowering_contract_kind ==
          :pqs_filled_source_cpb
    @test lw_inventory.selected_contract_kind_counts.white_lindsey_boundary_strata_count ==
          1
    @test lw_inventory.selected_contract_kind_counts.pqs_filled_source_cpb_count ==
          0

    pqs_inventory =
        GaussletBases._cartesian_selected_terminal_lowering_contract_inventory(
            lowering_inventory,
            :pqs,
        )
    _assert_selected_lowering_inventory_metadata_only(pqs_inventory)
    @test pqs_inventory.route_lowering_family == :pqs
    pqs_shell = _selected_lowering_contract_by_unit(pqs_inventory, :complete_shell)
    @test pqs_shell.lowering_contract_kind == :pqs_filled_source_cpb
    @test !pqs_shell.face_edge_corner_decomposition_required
    @test pqs_shell.source_cpb_count == 1
    @test only(pqs_inventory.unselected_contracts).lowering_contract_kind ==
          :white_lindsey_boundary_strata
    @test pqs_inventory.selected_contract_kind_counts.white_lindsey_boundary_strata_count ==
          0
    @test pqs_inventory.selected_contract_kind_counts.pqs_filled_source_cpb_count ==
          1

    for inventory in (lw_inventory, pqs_inventory)
        @test _selected_lowering_contract_by_unit(
            inventory,
            :direct_core,
        ).identity_like_source_contract
        @test _selected_lowering_contract_by_unit(
            inventory,
            :midpoint_slab,
        ).identity_like_source_contract
        @test _selected_lowering_contract_by_unit(
            inventory,
            :outer_mismatch,
        ).identity_like_source_contract
        central =
            _selected_lowering_contract_by_unit(inventory, :central_distorted)
        @test central.lowering_contract_kind == :distorted_product_box_comx
        @test !central.identity_like_source_contract
        @test inventory.selected_contract_kind_counts.direct_core_identity_cpb_count ==
              1
        @test inventory.selected_contract_kind_counts.direct_slab_identity_cpb_count ==
              1
        @test inventory.selected_contract_kind_counts.direct_boundary_slab_identity_cpb_count ==
              1
        @test inventory.selected_contract_kind_counts.distorted_product_box_comx_count ==
              1
    end

    @test_throws ArgumentError GaussletBases._cartesian_selected_terminal_lowering_contract_inventory(
        lowering_inventory,
        :pqs_source_box,
    )
end

# Integration/slow test. Do not include in default nested runner.

using Test
using GaussletBases

const CRCSelectedLowering = GaussletBases.CartesianRouteCore

function _selected_crc_test_cpb(intervals, role::Symbol)
    return (;
        object_kind = :synthetic_selected_lowering_cpb,
        intervals,
        role,
        support_count = prod(length, intervals),
    )
end

function _selected_crc_direct_contract(
    contract_index::Int,
    unit_key::Symbol,
    kind::Symbol,
    intervals;
    real_contract_shape::Bool = false,
)
    source_cpb = _selected_crc_test_cpb(
        intervals,
        Symbol(unit_key, :_source_cpb),
    )
    return (;
        object_kind = :cartesian_terminal_region_lowering_contract,
        status = :planned_not_materialized,
        contract_index,
        contract_key = Symbol(unit_key, :_, kind),
        unit_index = contract_index,
        unit_key,
        unit_role = unit_key,
        terminal_region_kind = unit_key,
        lowering_contract_kind = kind,
        identity_like_source_contract = true,
        support_count = prod(length, intervals),
        source_cpb = real_contract_shape ? nothing : source_cpb,
        source_cpb_plan_box = real_contract_shape ? intervals : nothing,
        source_cpb_plan_kind = kind,
        source_cpb_plan_status = :planned_not_materialized,
        source_cpb_count = 1,
        owned_support = real_contract_shape ? nothing : (;
            object_kind = :synthetic_selected_lowering_owned_support,
            outer_cpb = source_cpb,
            support_kind = :synthetic_direct_owned_cpb,
            support_count = source_cpb.support_count,
        ),
    )
end

function _selected_crc_pqs_contract()
    outer_cpb =
        _selected_crc_test_cpb((1:5, 1:5, 1:5), :pqs_complete_shell_outer_cpb)
    inner_cpb =
        _selected_crc_test_cpb((2:4, 2:4, 2:4), :pqs_complete_shell_inner_cpb)
    return (;
        object_kind = :cartesian_terminal_region_lowering_contract,
        status = :planned_not_materialized,
        contract_index = 4,
        contract_key = :complete_shell_pqs_filled_source_cpb,
        unit_index = 4,
        unit_key = :complete_shell,
        unit_role = :complete_shell,
        terminal_region_kind = :complete_shell,
        lowering_contract_kind = :pqs_filled_source_cpb,
        identity_like_source_contract = false,
        support_count = outer_cpb.support_count - inner_cpb.support_count,
        source_cpb_plan_box = outer_cpb.intervals,
        outer_box = outer_cpb.intervals,
        inner_exclusion_box = inner_cpb.intervals,
        source_cpb_plan_kind = :pqs_filled_source_cpb,
        source_cpb_plan_status = :planned_not_materialized,
        source_cpb_count = 1,
    )
end

function _selected_crc_unsupported_contract(
    contract_index::Int,
    unit_key::Symbol,
    kind::Symbol,
)
    return (;
        object_kind = :cartesian_terminal_region_lowering_contract,
        status = :planned_not_materialized,
        contract_index,
        contract_key = Symbol(unit_key, :_, kind),
        unit_index = contract_index,
        unit_key,
        unit_role = unit_key,
        terminal_region_kind = unit_key,
        lowering_contract_kind = kind,
        identity_like_source_contract = false,
    )
end

function _selected_crc_inventory()
    selected_contracts = (
        _selected_crc_direct_contract(
            1,
            :direct_core,
            :direct_core_identity_cpb,
            (1:2, 1:2, 1:2);
            real_contract_shape = true,
        ),
        _selected_crc_direct_contract(
            2,
            :midpoint_slab,
            :direct_slab_identity_cpb,
            (3:5, 2:2, 1:4),
        ),
        _selected_crc_direct_contract(
            3,
            :boundary_slab,
            :direct_boundary_slab_identity_cpb,
            (1:4, 4:4, 1:3),
        ),
        _selected_crc_pqs_contract(),
        _selected_crc_unsupported_contract(
            5,
            :lw_complete_shell,
            :white_lindsey_boundary_strata,
        ),
        _selected_crc_unsupported_contract(
            6,
            :central_distorted,
            :distorted_product_box_comx,
        ),
    )
    return (;
        object_kind = :cartesian_selected_terminal_lowering_contract_inventory,
        status = :available_selected_terminal_lowering_contract_inventory,
        route_lowering_family = :pqs,
        selected_contract_count = length(selected_contracts),
        selected_contracts,
        selected_contract_records = selected_contracts,
    )
end

function _selected_crc_entry_by_kind(inventory, kind::Symbol)
    return only(
        entry for entry in inventory.sidecar_entries
        if entry.lowering_contract_kind == kind
    )
end

@testset "selected terminal lowering CRC sidecars are metadata-only" begin
    sidecar_inventory =
        GaussletBases._cartesian_route_core_selected_terminal_lowering_sidecar_inventory(
            _selected_crc_inventory(),
        )

    @test sidecar_inventory.object_kind ==
          :cartesian_route_core_selected_terminal_lowering_sidecar_inventory
    @test sidecar_inventory.status ==
          :partial_selected_terminal_lowering_crc_sidecar_inventory
    @test sidecar_inventory.private_development_only
    @test sidecar_inventory.selected_contract_count == 6
    @test sidecar_inventory.sidecar_available_count == 4
    @test sidecar_inventory.sidecar_missing_count == 2
    @test length(sidecar_inventory.sidecar_entries) == 6
    @test length(sidecar_inventory.route_core_sidecars) == 4
    @test !sidecar_inventory.route_core_sidecar_inventory_complete
    @test !sidecar_inventory.final_retained_unit_inventory_available
    @test !sidecar_inventory.pair_inventory_available
    @test sidecar_inventory.pair_inventory_status ==
          :not_available_crc_sidecar_metadata_only
    @test !sidecar_inventory.pair_operator_blocks_materialized
    @test !sidecar_inventory.hamiltonian_data_materialized
    @test !sidecar_inventory.artifacts_materialized

    for kind in (
        :direct_core_identity_cpb,
        :direct_slab_identity_cpb,
        :direct_boundary_slab_identity_cpb,
    )
        entry = _selected_crc_entry_by_kind(sidecar_inventory, kind)
        @test entry.route_core_sidecar_available
        @test isnothing(entry.missing_route_core_sidecar_reason)
        sidecar = entry.route_core_sidecar
        @test sidecar.object_kind == :cartesian_route_core_sidecar
        @test sidecar.sidecar_source ==
              :selected_terminal_direct_identity_lowering_contract
        @test !sidecar.coefficient_maps_materialized
        @test !sidecar.transform_contracts_materialized
        @test !sidecar.retained_spaces_materialized
        @test !sidecar.operator_blocks_materialized
        @test !sidecar.pair_operator_blocks_materialized
        @test !sidecar.hamiltonian_data_materialized
        @test !sidecar.artifacts_materialized
        @test sidecar.shellification_region isa
              CRCSelectedLowering.ShellificationRegion
        @test sidecar.lowering_source isa CRCSelectedLowering.LoweringSource
        @test sidecar.intermediate_retained_space isa
              CRCSelectedLowering.IntermediateRetainedSpace
        @test sidecar.shell_realization isa CRCSelectedLowering.ShellRealization
        @test sidecar.final_retained_unit isa
              CRCSelectedLowering.FinalRetainedUnit
        @test CRCSelectedLowering.lowering_recipe(sidecar.lowering_source) ==
              kind
        @test sidecar.intermediate_retained_space.retained_rule ==
              :identity_source_modes
        @test sidecar.shell_realization.realization_kind ==
              :direct_or_trivial_embedding
        source_cpb = only(CRCSelectedLowering.source_cpbs(sidecar.lowering_source))
        owned_cpb = only(
            CRCSelectedLowering.owned_support(
                sidecar.shellification_region,
            ).cpbs,
        )
        @test CRCSelectedLowering.intervals(source_cpb) ==
              CRCSelectedLowering.intervals(owned_cpb)
        @test sidecar.final_retained_unit.metadata.pair_planning_input == false
    end

    pqs_entry =
        _selected_crc_entry_by_kind(sidecar_inventory, :pqs_filled_source_cpb)
    @test pqs_entry.route_core_sidecar_available
    pqs_sidecar = pqs_entry.route_core_sidecar
    @test pqs_sidecar.sidecar_source ==
          :selected_terminal_pqs_filled_source_contract
    @test CRCSelectedLowering.lowering_recipe(pqs_sidecar.lowering_source) ==
          :pqs_filled_source_cpb
    pqs_source_cpb =
        only(CRCSelectedLowering.source_cpbs(pqs_sidecar.lowering_source))
    @test CRCSelectedLowering.codimension(pqs_source_cpb) == 0
    @test pqs_sidecar.intermediate_retained_space.retained_rule ==
          :pqs_boundary_comx_product_modes
    @test pqs_sidecar.intermediate_retained_space.dimension ==
          CRCSelectedLowering.boundary_product_mode_count(
              CRCSelectedLowering.shape(pqs_source_cpb),
          )
    @test pqs_sidecar.shell_realization.realization_kind ==
          :shell_projection_lowdin
    @test pqs_sidecar.shell_realization.status ==
          :shell_projection_lowdin_planned_not_materialized
    @test pqs_sidecar.shell_realization.metadata.shell_projection_planned
    @test pqs_sidecar.shell_realization.metadata.lowdin_cleanup_planned
    @test CRCSelectedLowering.support_count(pqs_source_cpb) !=
          CRCSelectedLowering.support_count(pqs_sidecar.shellification_region)
    @test !pqs_sidecar.pair_operator_blocks_materialized
    @test !pqs_sidecar.hamiltonian_data_materialized
    @test !pqs_sidecar.artifacts_materialized

    lw_entry =
        _selected_crc_entry_by_kind(sidecar_inventory, :white_lindsey_boundary_strata)
    @test !lw_entry.route_core_sidecar_available
    @test lw_entry.missing_route_core_sidecar_reason ==
          :white_lindsey_leaf_sidecars_not_yet_expanded

    distorted_entry =
        _selected_crc_entry_by_kind(sidecar_inventory, :distorted_product_box_comx)
    @test !distorted_entry.route_core_sidecar_available
    @test distorted_entry.missing_route_core_sidecar_reason ==
          :distorted_product_box_comx_sidecar_not_yet_mapped

    @test_throws ArgumentError GaussletBases._cartesian_route_core_selected_terminal_lowering_sidecar(
        _selected_crc_unsupported_contract(
            7,
            :lw_complete_shell_single,
            :white_lindsey_boundary_strata,
        ),
    )

    bad_support_count_contract = merge(
        _selected_crc_direct_contract(
            8,
            :bad_direct_core,
            :direct_core_identity_cpb,
            (1:2, 1:2, 1:2);
            real_contract_shape = true,
        ),
        (; support_count = 999),
    )
    @test_throws ArgumentError GaussletBases._cartesian_route_core_selected_terminal_lowering_sidecar(
        bad_support_count_contract,
    )
end

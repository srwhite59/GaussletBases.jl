using Test
using GaussletBases

const CTLForTerminalRouteFingerprint = GaussletBases.CartesianTerminalLowering
const CRUForTerminalRouteFingerprint = GaussletBases.CartesianRetainedUnits
const CRTCForTerminalRouteFingerprint =
    GaussletBases.CartesianRetainedUnitTransformContracts
const CUPForTerminalRouteFingerprint = GaussletBases.CartesianUnitPairs
const CPOPForTerminalRouteFingerprint = GaussletBases.CartesianPairOperatorPlans
const CRCForTerminalRouteFingerprint = GaussletBases.CartesianRouteCore
const CPBForTerminalRouteFingerprint = GaussletBases.CartesianCPB

function _terminal_route_fingerprint_direct_contract()
    source = CPBForTerminalRouteFingerprint.filled_cpb(
        1:3,
        1:3,
        1:3;
        role = :terminal_route_fingerprint_direct_core_cpb,
    )
    return CTLForTerminalRouteFingerprint.TerminalLoweringContract(
        :terminal_route_fingerprint_direct_core,
        :terminal_route_fingerprint_direct_core,
        :synthetic_direct_region,
        :direct_core,
        :direct_core_identity_cpb,
        CRCForTerminalRouteFingerprint.owned_cpb(
            source;
            support_kind = :terminal_route_fingerprint_direct_owned_support,
            metadata = (; terminal_region_key = :terminal_route_fingerprint_direct_core),
        ),
        (source,),
        :direct_source_modes,
        :direct_or_trivial_embedding,
        :one_terminal_region,
        false,
        (; identity_like = true),
    )
end

function _terminal_route_fingerprint_pqs_contract()
    outer = CPBForTerminalRouteFingerprint.filled_cpb(
        1:5,
        1:5,
        1:5;
        role = :terminal_route_fingerprint_outer_box,
    )
    inner = CPBForTerminalRouteFingerprint.filled_cpb(
        2:4,
        2:4,
        2:4;
        role = :terminal_route_fingerprint_inner_box,
    )
    support = CRCForTerminalRouteFingerprint.complete_shell_support(
        outer,
        inner;
        metadata = (; terminal_region_key = :terminal_route_fingerprint_pqs_shell),
    )
    return CTLForTerminalRouteFingerprint.TerminalLoweringContract(
        :terminal_route_fingerprint_pqs_shell,
        :terminal_route_fingerprint_pqs_shell,
        :synthetic_complete_shell_region,
        :complete_shell,
        :pqs_filled_source_cpb,
        support,
        (outer,),
        :pqs_boundary_comx_product_modes,
        :shell_projection_lowdin,
        :one_terminal_region,
        false,
        (; q = 3, source_mode_shape = nothing),
    )
end

function _terminal_route_fingerprint_lowering_plan()
    selected = (
        _terminal_route_fingerprint_direct_contract(),
        _terminal_route_fingerprint_pqs_contract(),
    )
    summary = (;
        object_kind = :synthetic_terminal_lowering_plan_summary,
        status = :available_terminal_lowering_plan,
        policy_kind = :pqs_terminal_lowering,
        terminal_region_count = length(selected),
        available_contract_count = length(selected),
        available_contract_kinds = Tuple(contract.lowering_kind for contract in selected),
        selected_contract_count = length(selected),
        selected_contract_kinds = Tuple(contract.lowering_kind for contract in selected),
        all_terminal_regions_have_selected_contract = true,
        materialized = false,
        retained_spaces_materialized = false,
        final_retained_units_materialized = false,
        pair_inventory_materialized = false,
        operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
    )
    return CTLForTerminalRouteFingerprint.TerminalLoweringPlan(
        CTLForTerminalRouteFingerprint.PQSLowering(q = 3),
        selected,
        selected,
        summary,
        (; fixture = :terminal_route_retained_units_fingerprint),
    )
end

function _terminal_route_fingerprint_kind_count(summary, unit_kind::Symbol)
    for entry in summary.unit_kind_counts
        entry.unit_kind == unit_kind && return entry.count
    end
    return 0
end

function _terminal_route_fingerprint_count(counts, field::Symbol, value)
    matches = Tuple(entry for entry in counts if getproperty(entry, field) == value)
    isempty(matches) && return 0
    return only(matches).count
end

@testset "terminal route state metadata-chain fingerprint" begin
    lowering_plan = _terminal_route_fingerprint_lowering_plan()
    lowering_summary = CTLForTerminalRouteFingerprint.summary(lowering_plan)
    terminal_route_state =
        GaussletBases._pqs_source_box_route_driver_terminal_route_state(;
            status = lowering_summary.status,
            selected = true,
            route_lowering_family = :pqs_terminal_route_fingerprint,
            lowering_plan,
            lowering_summary,
            pair_operator_route_core_sidecars = false,
    )
    terminal_route_summary = terminal_route_state.summary
    retained_unit_summary = terminal_route_summary.retained_unit_summary
    transform_contract_summary =
        terminal_route_summary.retained_unit_transform_contract_summary
    unit_pair_summary = terminal_route_summary.unit_pair_summary
    pair_operator_summary = terminal_route_summary.pair_operator_summary

    @test terminal_route_state.retained_unit_plan isa
          CRUForTerminalRouteFingerprint.RetainedUnitPlan
    @test terminal_route_state.retained_unit_transform_contract_plan isa
          CRTCForTerminalRouteFingerprint.RetainedUnitTransformContractPlan
    @test terminal_route_state.unit_pair_plan isa
          CUPForTerminalRouteFingerprint.UnitPairPlan
    @test terminal_route_state.pair_operator_plan isa
          CPOPForTerminalRouteFingerprint.PairOperatorPlan
    @test retained_unit_summary.status == :available_retained_unit_plan
    @test retained_unit_summary.selected_contract_count == 2
    @test retained_unit_summary.retained_unit_count == 2
    @test _terminal_route_fingerprint_kind_count(
        retained_unit_summary,
        :direct_cpb_retained_unit,
    ) == 1
    @test _terminal_route_fingerprint_kind_count(
        retained_unit_summary,
        :pqs_shell_retained_unit,
    ) == 1
    @test retained_unit_summary.route_core_final_unit_available_count ==
          retained_unit_summary.retained_unit_count
    @test retained_unit_summary.route_core_final_unit_blocked_count == 0
    @test !retained_unit_summary.materialized
    @test !retained_unit_summary.transforms_materialized
    @test !retained_unit_summary.coefficient_maps_materialized
    @test !retained_unit_summary.pair_inventory_materialized
    @test !retained_unit_summary.operator_blocks_materialized
    @test !retained_unit_summary.hamiltonian_data_materialized
    @test terminal_route_summary.retained_unit_summary.status ==
          retained_unit_summary.status
    @test terminal_route_summary.retained_unit_transform_contract_summary.status ==
          transform_contract_summary.status
    @test terminal_route_summary.unit_pair_summary.status == unit_pair_summary.status
    @test terminal_route_summary.selected_contract_count == 2
    @test transform_contract_summary.status ==
          :available_retained_unit_transform_contract_plan
    @test transform_contract_summary.retained_unit_count ==
          retained_unit_summary.retained_unit_count
    @test transform_contract_summary.transform_contract_count ==
          retained_unit_summary.retained_unit_count
    @test _terminal_route_fingerprint_count(
        transform_contract_summary.transform_path_counts,
        :transform_path,
        :direct_identity_transform_contract,
    ) == 1
    @test _terminal_route_fingerprint_count(
        transform_contract_summary.transform_path_counts,
        :transform_path,
        :pqs_source_modes_boundary_selection_shell_realization_contract,
    ) == 1
    @test !transform_contract_summary.materialized
    @test !transform_contract_summary.transforms_materialized
    @test !transform_contract_summary.coefficient_maps_materialized
    @test !transform_contract_summary.pair_inventory_materialized
    @test !transform_contract_summary.operator_blocks_materialized
    @test !transform_contract_summary.hamiltonian_data_materialized
    @test !transform_contract_summary.artifacts_materialized
    @test terminal_route_summary.pair_inventory_available
    @test unit_pair_summary.status == :available_unit_pair_plan
    @test unit_pair_summary.retained_unit_count ==
          retained_unit_summary.retained_unit_count
    @test unit_pair_summary.pair_count == 3
    @test unit_pair_summary.expected_upper_triangular_pair_count == 3
    @test unit_pair_summary.route_core_pair_inventory_available
    @test unit_pair_summary.route_core_pair_inventory_status ==
          :available_route_core_pair_inventory
    @test unit_pair_summary.route_core_pair_count == unit_pair_summary.pair_count
    @test !unit_pair_summary.materialized
    @test !unit_pair_summary.pair_inventory_materialized
    @test !unit_pair_summary.source_operator_blocks_materialized
    @test !unit_pair_summary.operator_blocks_materialized
    @test !unit_pair_summary.hamiltonian_data_materialized
    @test !unit_pair_summary.artifacts_materialized
    @test pair_operator_summary.status == :available_pair_operator_plan
    @test pair_operator_summary.retained_unit_count ==
          retained_unit_summary.retained_unit_count
    @test pair_operator_summary.unit_pair_count == unit_pair_summary.pair_count
    @test pair_operator_summary.pair_operator_plan_count ==
          pair_operator_summary.unit_pair_count
    @test pair_operator_summary.expected_pair_operator_plan_count ==
          unit_pair_summary.expected_upper_triangular_pair_count
    @test !pair_operator_summary.route_core_pair_operator_plan_inventory_available
    @test pair_operator_summary.route_core_pair_operator_plan_inventory_status ==
          :not_requested_route_core_pair_operator_plan_inventory
    @test _terminal_route_fingerprint_count(
        pair_operator_summary.realization_path_counts,
        :realization_path,
        (;
            left = :identity_or_trivial_embedding,
            right = :identity_or_trivial_embedding,
        ),
    ) == 1
    @test _terminal_route_fingerprint_count(
        pair_operator_summary.realization_path_counts,
        :realization_path,
        (;
            left = :identity_or_trivial_embedding,
            right = :shell_projection_lowdin_planned,
        ),
    ) == 1
    @test _terminal_route_fingerprint_count(
        pair_operator_summary.realization_path_counts,
        :realization_path,
        (;
            left = :shell_projection_lowdin_planned,
            right = :shell_projection_lowdin_planned,
        ),
    ) == 1
    @test !pair_operator_summary.materialized
    @test !pair_operator_summary.source_operator_blocks_materialized
    @test !pair_operator_summary.final_pair_blocks_materialized
    @test !pair_operator_summary.operator_blocks_materialized
    @test !pair_operator_summary.hamiltonian_data_materialized
    @test !pair_operator_summary.artifacts_materialized
    @test !terminal_route_summary.operator_blocks_materialized
    @test !terminal_route_summary.hamiltonian_data_materialized

    unselected_state =
        GaussletBases._pqs_source_box_route_driver_terminal_route_state(;
            status = :available_terminal_lowering_plan,
            selected = false,
            route_lowering_family = :pqs_terminal_route_fingerprint,
            lowering_plan,
            lowering_summary,
            pair_operator_route_core_sidecars = false,
        )
    unselected_summary = unselected_state.summary.retained_unit_summary
    unselected_transform_summary =
        unselected_state.summary.retained_unit_transform_contract_summary
    unselected_pair_summary = unselected_state.summary.unit_pair_summary
    unselected_pair_operator_summary =
        unselected_state.summary.pair_operator_summary

    @test unselected_state.retained_unit_plan === nothing
    @test unselected_state.retained_unit_transform_contract_plan === nothing
    @test unselected_state.unit_pair_plan === nothing
    @test unselected_state.pair_operator_plan === nothing
    @test unselected_summary.status == :not_selected
    @test unselected_summary.retained_unit_count == 0
    @test unselected_summary.route_core_final_unit_available_count == 0
    @test !unselected_summary.materialized
    @test !unselected_summary.transforms_materialized
    @test !unselected_summary.coefficient_maps_materialized
    @test !unselected_summary.pair_inventory_materialized
    @test !unselected_summary.operator_blocks_materialized
    @test !unselected_summary.hamiltonian_data_materialized
    @test unselected_transform_summary.status == :not_selected
    @test unselected_transform_summary.retained_unit_count == 0
    @test unselected_transform_summary.transform_contract_count == 0
    @test !unselected_transform_summary.materialized
    @test !unselected_transform_summary.transforms_materialized
    @test !unselected_transform_summary.coefficient_maps_materialized
    @test !unselected_transform_summary.pair_inventory_materialized
    @test !unselected_transform_summary.operator_blocks_materialized
    @test !unselected_transform_summary.hamiltonian_data_materialized
    @test !unselected_transform_summary.artifacts_materialized
    @test unselected_pair_summary.status == :not_selected
    @test unselected_pair_summary.pair_count == 0
    @test unselected_pair_summary.expected_upper_triangular_pair_count == 0
    @test !unselected_pair_summary.route_core_pair_inventory_available
    @test !unselected_pair_summary.materialized
    @test !unselected_pair_summary.pair_inventory_materialized
    @test !unselected_pair_summary.source_operator_blocks_materialized
    @test !unselected_pair_summary.operator_blocks_materialized
    @test !unselected_pair_summary.hamiltonian_data_materialized
    @test !unselected_pair_summary.artifacts_materialized
    @test unselected_pair_operator_summary.status == :not_selected
    @test unselected_pair_operator_summary.pair_operator_plan_count == 0
    @test unselected_pair_operator_summary.expected_pair_operator_plan_count == 0
    @test !unselected_pair_operator_summary.route_core_pair_operator_plan_inventory_available
    @test !unselected_pair_operator_summary.materialized
    @test !unselected_pair_operator_summary.source_operator_blocks_materialized
    @test !unselected_pair_operator_summary.final_pair_blocks_materialized
    @test !unselected_pair_operator_summary.operator_blocks_materialized
    @test !unselected_pair_operator_summary.hamiltonian_data_materialized
    @test !unselected_pair_operator_summary.artifacts_materialized
end

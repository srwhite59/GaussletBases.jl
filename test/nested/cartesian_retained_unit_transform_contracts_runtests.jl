using Test
using GaussletBases

const CTLForTransformContracts = GaussletBases.CartesianTerminalLowering
const CRUForTransformContracts = GaussletBases.CartesianRetainedUnits
const CRTC = GaussletBases.CartesianRetainedUnitTransformContracts
const CUPForTransformContracts = GaussletBases.CartesianUnitPairs
const CPOPForTransformContracts = GaussletBases.CartesianPairOperatorPlans
const CRCForTransformContracts = GaussletBases.CartesianRouteCore
const CPBForTransformContracts = GaussletBases.CartesianCPB
const CRPSForTransformContracts = GaussletBases.CartesianRawProductSources

function _transform_contract_count_by_field(summary, field::Symbol, value::Symbol)
    count_field = Symbol(String(field), "_counts")
    for entry in getfield(summary, count_field)
        getfield(entry, field) == value && return entry.count
    end
    return 0
end

function _transform_contract_direct_contract(kind::Symbol, key::Symbol, source_cpb)
    return CTLForTransformContracts.TerminalLoweringContract(
        Symbol(String(key), "_", String(kind)),
        key,
        :synthetic_direct_region,
        kind === :direct_core_identity_cpb ? :direct_core :
        kind === :direct_slab_identity_cpb ? :direct_midpoint_slab :
        :outer_mismatch_slab,
        kind,
        CRCForTransformContracts.owned_cpb(
            source_cpb;
            support_kind = :synthetic_direct_owned_support,
            metadata = (; terminal_region_key = key),
        ),
        (source_cpb,),
        :direct_source_modes,
        :direct_or_trivial_embedding,
        :one_terminal_region,
        false,
        (; identity_like = true),
    )
end

function _transform_contract_complete_shell_contracts()
    outer = CPBForTransformContracts.filled_cpb(
        1:5,
        1:5,
        1:5;
        role = :synthetic_outer_box,
    )
    inner = CPBForTransformContracts.filled_cpb(
        2:4,
        2:4,
        2:4;
        role = :synthetic_inner_box,
    )
    support = CRCForTransformContracts.complete_shell_support(
        outer,
        inner;
        metadata = (; terminal_region_key = :synthetic_complete_shell),
    )
    strata = CPBForTransformContracts.complete_shell_boundary_strata(outer, inner)
    lw_contract = CTLForTransformContracts.TerminalLoweringContract(
        :synthetic_complete_shell_white_lindsey_boundary_strata,
        :synthetic_complete_shell,
        :synthetic_complete_shell_region,
        :complete_shell,
        :white_lindsey_boundary_strata,
        support,
        strata.all_strata,
        :white_lindsey_boundary_stratum_product,
        :direct_or_trivial_embedding,
        :boundary_stratum_children,
        false,
        (;
            facet_count = length(strata.facets),
            edge_count = length(strata.edges),
            corner_count = length(strata.corners),
        ),
    )
    pqs_contract = CTLForTransformContracts.TerminalLoweringContract(
        :synthetic_complete_shell_pqs_filled_source_cpb,
        :synthetic_complete_shell,
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
    return (; lw_contract, pqs_contract)
end

function _transform_contract_distorted_contract()
    source = CPBForTransformContracts.cpb(
        1:7,
        1:3,
        1:3;
        role = :synthetic_distorted_product_source_cpb,
    )
    return CTLForTransformContracts.TerminalLoweringContract(
        :synthetic_distorted_product_box_comx,
        :synthetic_distorted_product_box,
        :synthetic_distorted_region,
        :central_distorted_product_box,
        :distorted_product_box_comx,
        CRCForTransformContracts.owned_cpb(
            source;
            support_kind = :synthetic_distorted_owned_support,
            metadata = (; terminal_region_key = :synthetic_distorted_product_box),
        ),
        (source,),
        :distorted_product_comx_all_axes,
        :distorted_product_realization_planned,
        :one_terminal_region,
        false,
        (; q = 3, L = 7, source_mode_shape = (7, 3, 3), aspect_ratio = 7 / 3),
    )
end

function _transform_contract_pqs_contract(
    source;
    contract_key::Symbol,
    terminal_region_key::Symbol,
    metadata,
)
    return CTLForTransformContracts.TerminalLoweringContract(
        contract_key,
        terminal_region_key,
        :synthetic_pqs_region,
        :complete_shell,
        :pqs_filled_source_cpb,
        CRCForTransformContracts.owned_cpb(
            source;
            support_kind = :synthetic_pqs_owned_support,
            metadata = (; terminal_region_key),
        ),
        (source,),
        :pqs_boundary_comx_product_modes,
        :shell_projection_lowdin,
        :one_terminal_region,
        false,
        metadata,
    )
end

function _transform_contract_lowering_plan(policy, selected, policy_kind::Symbol)
    selected_tuple = Tuple(selected)
    return CTLForTransformContracts.TerminalLoweringPlan(
        policy,
        selected_tuple,
        selected_tuple,
        (;
            object_kind = :synthetic_terminal_lowering_plan_summary,
            status = :available_terminal_lowering_plan,
            policy_kind,
            terminal_region_count = length(selected_tuple),
            available_contract_count = length(selected_tuple),
            available_contract_kinds =
                Tuple(contract.lowering_kind for contract in selected_tuple),
            selected_contract_count = length(selected_tuple),
            selected_contract_kinds =
                Tuple(contract.lowering_kind for contract in selected_tuple),
            all_terminal_regions_have_selected_contract = true,
            materialized = false,
            retained_spaces_materialized = false,
            final_retained_units_materialized = false,
            pair_inventory_materialized = false,
            operator_blocks_materialized = false,
            hamiltonian_data_materialized = false,
        ),
        (; fixture = :cartesian_retained_unit_transform_contracts),
    )
end

function _unknown_transform_contract_retained_plan(source_plan)
    unit = first(CRUForTransformContracts.units(source_plan))
    unknown_unit = CRUForTransformContracts.RetainedUnitRecord(
        unit.unit_key,
        unit.unit_index,
        :unknown_retained_unit,
        unit.source_contract_key,
        unit.terminal_region_key,
        unit.terminal_region_role,
        unit.terminal_region_kind,
        unit.lowering_kind,
        unit.retained_rule,
        unit.realization_rule,
        unit.owned_support,
        unit.source_cpbs,
        unit.source_cpb_index,
        unit.dimension_status,
        unit.dimension,
        unit.column_range_status,
        unit.column_range,
        unit.route_core_final_unit,
        unit.materialized,
        unit.metadata,
    )
    return CRUForTransformContracts.RetainedUnitPlan(
        source_plan.policy,
        source_plan.lowering_plan,
        (unknown_unit,),
        source_plan.summary,
        (; fixture = :unknown_retained_unit_transform_contract),
    )
end

const _TRANSFORM_DIRECT_CORE = _transform_contract_direct_contract(
    :direct_core_identity_cpb,
    :synthetic_direct_core,
    CPBForTransformContracts.filled_cpb(1:3, 1:3, 1:3; role = :synthetic_direct_core_cpb),
)
const _TRANSFORM_DIRECT_SLAB = _transform_contract_direct_contract(
    :direct_slab_identity_cpb,
    :synthetic_direct_slab,
    CPBForTransformContracts.slab_cpb(1:3, 1:3, 4:4; role = :synthetic_direct_slab_cpb),
)
const _TRANSFORM_SHELL_CONTRACTS = _transform_contract_complete_shell_contracts()
const _TRANSFORM_DISTORTED = _transform_contract_distorted_contract()

const _TRANSFORM_WL_LOWERING_PLAN = _transform_contract_lowering_plan(
    CTLForTransformContracts.WhiteLindseyLowering(),
    (
        _TRANSFORM_DIRECT_CORE,
        _TRANSFORM_DIRECT_SLAB,
        _TRANSFORM_SHELL_CONTRACTS.lw_contract,
        _TRANSFORM_DISTORTED,
    ),
    :white_lindsey_terminal_lowering,
)

const _TRANSFORM_PQS_LOWERING_PLAN = _transform_contract_lowering_plan(
    CTLForTransformContracts.PQSLowering(q = 3),
    (
        _TRANSFORM_DIRECT_CORE,
        _TRANSFORM_DIRECT_SLAB,
        _TRANSFORM_SHELL_CONTRACTS.pqs_contract,
        _TRANSFORM_DISTORTED,
    ),
    :pqs_terminal_lowering,
)

const _TRANSFORM_WL_RETAINED_PLAN =
    CRUForTransformContracts.retained_unit_plan(_TRANSFORM_WL_LOWERING_PLAN)
const _TRANSFORM_PQS_RETAINED_PLAN =
    CRUForTransformContracts.retained_unit_plan(_TRANSFORM_PQS_LOWERING_PLAN)

@testset "CartesianRetainedUnitTransformContracts unavailable summary" begin
    contract_summary = CRTC.unavailable_summary(:not_selected, :not_selected_route)

    @test contract_summary.object_kind ==
          :cartesian_retained_unit_transform_contract_plan_summary
    @test contract_summary.status == :not_selected
    @test contract_summary.blocker == :not_selected_route
    @test contract_summary.retained_unit_count == 0
    @test contract_summary.transform_contract_count == 0
    @test contract_summary.transform_path_counts == ()
    @test contract_summary.realization_path_counts == ()
    @test contract_summary.blocked_contract_count == 0
    @test !contract_summary.materialized
    @test !contract_summary.transforms_materialized
    @test !contract_summary.coefficient_maps_materialized
    @test !contract_summary.pair_inventory_materialized
    @test !contract_summary.operator_blocks_materialized
    @test !contract_summary.hamiltonian_data_materialized
    @test !contract_summary.artifacts_materialized
    @test contract_summary.raw_product_source_plan_available_count == 0
    @test contract_summary.raw_product_source_plan_blocked_count == 0
end

@testset "CartesianRetainedUnitTransformContracts path classification" begin
    wl_plan =
        CRTC.retained_unit_transform_contract_plan(_TRANSFORM_WL_RETAINED_PLAN)
    wl_summary = CRTC.summary(wl_plan)

    @test wl_summary.status == :available_retained_unit_transform_contract_plan
    @test wl_summary.retained_unit_count ==
          length(CRUForTransformContracts.units(_TRANSFORM_WL_RETAINED_PLAN))
    @test wl_summary.transform_contract_count == wl_summary.retained_unit_count
    @test wl_summary.blocked_contract_count == 0
    @test _transform_contract_count_by_field(
        wl_summary,
        :transform_path,
        :direct_identity_transform_contract,
    ) == 2
    @test _transform_contract_count_by_field(
        wl_summary,
        :transform_path,
        :white_lindsey_boundary_stratum_product_contract,
    ) == 26
    @test _transform_contract_count_by_field(
        wl_summary,
        :transform_path,
        :distorted_product_comx_contract,
    ) == 1
    @test _transform_contract_count_by_field(
        wl_summary,
        :realization_path,
        :identity_or_trivial_embedding,
    ) == 28
    @test _transform_contract_count_by_field(
        wl_summary,
        :realization_path,
        :distorted_product_realization_planned,
    ) == 1
    @test all(!contract.materialized for contract in CRTC.transform_contracts(wl_plan))
    @test all(
        contract -> !haskey(contract.metadata, :raw_product_source_plan_status),
        CRTC.transform_contracts(wl_plan),
    )
    @test !wl_summary.materialized
    @test !wl_summary.transforms_materialized
    @test !wl_summary.coefficient_maps_materialized
    @test !wl_summary.pair_inventory_materialized
    @test !wl_summary.operator_blocks_materialized
    @test !wl_summary.hamiltonian_data_materialized
    @test !wl_summary.artifacts_materialized

    pqs_plan =
        CRTC.retained_unit_transform_contract_plan(_TRANSFORM_PQS_RETAINED_PLAN)
    pqs_summary = CRTC.summary(pqs_plan)
    pqs_contract = only(
        contract for contract in CRTC.transform_contracts(pqs_plan)
        if contract.unit_kind == :pqs_shell_retained_unit
    )

    @test pqs_summary.status == :available_retained_unit_transform_contract_plan
    @test pqs_summary.transform_contract_count == 4
    @test pqs_summary.raw_product_source_plan_available_count == 1
    @test pqs_summary.raw_product_source_plan_blocked_count == 0
    @test _transform_contract_count_by_field(
        pqs_summary,
        :transform_path,
        :direct_identity_transform_contract,
    ) == 2
    @test _transform_contract_count_by_field(
        pqs_summary,
        :transform_path,
        :pqs_source_modes_boundary_selection_shell_realization_contract,
    ) == 1
    @test _transform_contract_count_by_field(
        pqs_summary,
        :transform_path,
        :distorted_product_comx_contract,
    ) == 1
    @test pqs_contract.retained_rule == :pqs_boundary_comx_product_modes
    @test pqs_contract.realization_rule == :shell_projection_lowdin
    @test pqs_contract.realization_path == :shell_projection_lowdin_planned
    @test pqs_contract.dimension_status == :planned_not_materialized
    @test pqs_contract.column_range_status == :not_materialized
    @test pqs_contract.metadata.raw_product_source_plan_status ==
          :available_raw_product_box_plan
    @test pqs_contract.metadata.raw_product_source_plan isa
          CRPSForTransformContracts.RawProductBoxPlan
    @test pqs_contract.metadata.raw_product_source_plan.source_mode_dims ==
          (3, 3, 3)
    @test pqs_contract.metadata.raw_product_source_plan.source_mode_count == 27
    @test pqs_contract.metadata.raw_product_source_summary.source_mode_dims ==
          (3, 3, 3)
    @test pqs_contract.metadata.raw_product_source_summary.source_mode_count == 27
    @test pqs_contract.metadata.raw_product_source_summary.retained_rule_materialized ==
          false
    @test pqs_contract.metadata.raw_product_source_summary.shell_realization_materialized ==
          false
    @test pqs_contract.metadata.raw_product_source_summary.pair_blocks_materialized ==
          false
    @test pqs_contract.metadata.raw_product_source_summary.hamiltonian_data_materialized ==
          false
    @test pqs_contract.metadata.raw_product_source_summary.artifacts_materialized ==
          false
    @test pqs_contract.metadata.raw_product_source_retained_rule isa
          CRPSForTransformContracts.PQSBoundaryProductModeRetainedRule
    @test pqs_contract.metadata.raw_product_source_retained_rule.source_mode_dims ==
          (3, 3, 3)
    @test pqs_contract.metadata.raw_product_source_retained_rule.retained_count == 26
    @test pqs_contract.metadata.raw_product_source_retained_rule_summary.status ==
          :available_pqs_boundary_product_mode_retained_rule
    @test pqs_contract.metadata.raw_product_source_retained_rule_summary.retained_count ==
          26
    @test pqs_contract.metadata.raw_product_source_retained_rule_summary.transform_kind ==
          :source_mode_column_selector
    @test !pqs_contract.metadata.raw_product_source_retained_rule_summary.shell_realization_materialized
    @test !pqs_contract.metadata.raw_product_source_retained_rule_summary.lowdin_cleanup_used
    @test all(
        fact -> fact.coefficient_status === :not_materialized,
        CRPSForTransformContracts.axis_transform_facts(
            pqs_contract.metadata.raw_product_source_plan,
        ),
    )
    @test all(
        fact -> isnothing(fact.coefficient_matrix),
        CRPSForTransformContracts.axis_transform_facts(
            pqs_contract.metadata.raw_product_source_plan,
        ),
    )
    @test all(
        contract -> contract.unit_kind == :pqs_shell_retained_unit ||
                    !haskey(contract.metadata, :raw_product_source_plan_status),
        CRTC.transform_contracts(pqs_plan),
    )
end

@testset "CartesianRetainedUnitTransformContracts PQS raw source facts" begin
    explicit_source = CPBForTransformContracts.filled_cpb(
        10:11,
        20:22,
        30:33;
        role = :synthetic_explicit_pqs_source_cpb,
    )
    explicit_contract = _transform_contract_pqs_contract(
        explicit_source;
        contract_key = :synthetic_explicit_pqs_source,
        terminal_region_key = :synthetic_explicit_pqs,
        metadata = (; q = 9, source_mode_shape = (5, 4, 3)),
    )
    explicit_lowering_plan = _transform_contract_lowering_plan(
        CTLForTransformContracts.PQSLowering(q = 9),
        (explicit_contract,),
        :pqs_terminal_lowering,
    )
    explicit_retained_plan =
        CRUForTransformContracts.retained_unit_plan(explicit_lowering_plan)
    explicit_transform_plan =
        CRTC.retained_unit_transform_contract_plan(explicit_retained_plan)
    explicit_transform_contract =
        only(CRTC.transform_contracts(explicit_transform_plan))
    explicit_summary = CRTC.summary(explicit_transform_plan)

    @test explicit_summary.raw_product_source_plan_available_count == 1
    @test explicit_summary.raw_product_source_plan_blocked_count == 0
    @test explicit_transform_contract.metadata.raw_product_source_plan_status ==
          :available_raw_product_box_plan
    @test explicit_transform_contract.metadata.raw_product_source_plan.source_mode_dims ==
          (5, 4, 3)
    @test explicit_transform_contract.metadata.raw_product_source_plan.source_mode_count == 60
    @test explicit_transform_contract.metadata.raw_product_source_summary.source_mode_dims ==
          (5, 4, 3)
    @test explicit_transform_contract.metadata.raw_product_source_summary.source_mode_count ==
          60
    @test explicit_transform_contract.metadata.raw_product_source_retained_rule.retained_count ==
          54
    @test explicit_transform_contract.metadata.raw_product_source_retained_rule_summary.retained_count ==
          54
    @test explicit_transform_contract.metadata.raw_product_source_retained_rule_summary.shell_realization_materialized ==
          false
    @test explicit_transform_contract.metadata.raw_product_source_retained_rule_summary.lowdin_cleanup_used ==
          false

    missing_source = CPBForTransformContracts.filled_cpb(
        1:3,
        1:3,
        1:3;
        role = :synthetic_missing_dims_pqs_source_cpb,
    )
    missing_contract = _transform_contract_pqs_contract(
        missing_source;
        contract_key = :synthetic_missing_dims_pqs_source,
        terminal_region_key = :synthetic_missing_dims_pqs,
        metadata = (; q = nothing, source_mode_shape = nothing),
    )
    missing_lowering_plan = _transform_contract_lowering_plan(
        CTLForTransformContracts.PQSLowering(q = 3),
        (missing_contract,),
        :pqs_terminal_lowering,
    )
    missing_retained_plan =
        CRUForTransformContracts.retained_unit_plan(missing_lowering_plan)
    missing_transform_plan =
        CRTC.retained_unit_transform_contract_plan(missing_retained_plan)
    missing_transform_contract =
        only(CRTC.transform_contracts(missing_transform_plan))
    missing_summary = CRTC.summary(missing_transform_plan)

    @test missing_summary.raw_product_source_plan_available_count == 0
    @test missing_summary.raw_product_source_plan_blocked_count == 1
    @test isnothing(missing_transform_contract.metadata.raw_product_source_plan)
    @test missing_transform_contract.metadata.raw_product_source_plan_status ==
          :blocked_missing_source_mode_dims
    @test missing_transform_contract.metadata.raw_product_source_summary.status ==
          :blocked_missing_source_mode_dims
    @test missing_transform_contract.metadata.raw_product_source_summary.blocker ==
          :missing_pqs_source_mode_dims
    @test isnothing(missing_transform_contract.metadata.raw_product_source_retained_rule)
    @test isnothing(
        missing_transform_contract.metadata.raw_product_source_retained_rule_summary,
    )
end

@testset "CartesianRetainedUnitTransformContracts unknown retained unit blocks" begin
    unknown_retained_plan =
        _unknown_transform_contract_retained_plan(_TRANSFORM_PQS_RETAINED_PLAN)
    unknown_plan =
        CRTC.retained_unit_transform_contract_plan(unknown_retained_plan)
    unknown_summary = CRTC.summary(unknown_plan)
    unknown_contract = only(CRTC.transform_contracts(unknown_plan))

    @test unknown_summary.status == :blocked_retained_unit_transform_contract_plan
    @test unknown_summary.retained_unit_count == 1
    @test unknown_summary.transform_contract_count == 1
    @test unknown_summary.blocked_contract_count == 1
    @test unknown_summary.blocker ==
          :unclassified_retained_unit_transform_contract
    @test unknown_contract.transform_path ==
          :pending_retained_unit_transform_contract
    @test unknown_contract.blocker ==
          :unclassified_retained_unit_transform_contract
end

@testset "terminal route state carries retained-unit transform contracts" begin
    lowering_summary = CTLForTransformContracts.summary(_TRANSFORM_PQS_LOWERING_PLAN)
    retained_summary =
        CRUForTransformContracts.summary(_TRANSFORM_PQS_RETAINED_PLAN)
    unit_pair_summary =
        CUPForTransformContracts.unavailable_summary(:not_tested, :test_fixture)
    pair_operator_summary =
        CPOPForTransformContracts.unavailable_summary(:not_tested, :test_fixture)
    terminal_route_state =
        GaussletBases._pqs_source_box_route_driver_terminal_route_state(;
            status = lowering_summary.status,
            selected = true,
            route_lowering_family = :pqs_transform_contract_test,
            lowering_plan = _TRANSFORM_PQS_LOWERING_PLAN,
            lowering_summary,
            retained_unit_plan = _TRANSFORM_PQS_RETAINED_PLAN,
            retained_unit_summary = retained_summary,
            unit_pair_plan = :not_tested,
            unit_pair_summary,
            pair_operator_plan = :not_tested,
            pair_operator_summary,
        )
    terminal_summary = terminal_route_state.summary
    transform_summary =
        terminal_route_state.retained_unit_transform_contract_summary

    @test terminal_route_state.retained_unit_transform_contract_plan isa
          CRTC.RetainedUnitTransformContractPlan
    @test transform_summary.status ==
          :available_retained_unit_transform_contract_plan
    @test terminal_summary.retained_unit_transform_contract_summary.status ==
          transform_summary.status
    @test transform_summary.transform_contract_count ==
          terminal_route_state.retained_unit_summary.retained_unit_count
    @test !transform_summary.transforms_materialized
    @test !transform_summary.coefficient_maps_materialized

    unselected_state =
        GaussletBases._pqs_source_box_route_driver_terminal_route_state(;
            status = lowering_summary.status,
            selected = false,
            route_lowering_family = :pqs_transform_contract_test,
            lowering_plan = _TRANSFORM_PQS_LOWERING_PLAN,
            lowering_summary,
        )
    unselected_transform_summary =
        unselected_state.summary.retained_unit_transform_contract_summary

    @test unselected_state.retained_unit_transform_contract_plan === nothing
    @test unselected_state.retained_unit_transform_contract_summary.status == :not_selected
    @test unselected_transform_summary.transform_contract_count == 0
    @test !unselected_transform_summary.transforms_materialized
end

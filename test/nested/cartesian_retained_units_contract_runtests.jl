using Test
using GaussletBases

const CTLForRetainedUnits = GaussletBases.CartesianTerminalLowering
const CRU = GaussletBases.CartesianRetainedUnits
const CRCForRetainedUnits = GaussletBases.CartesianRouteCore
const CPBForRetainedUnits = GaussletBases.CartesianCPB

function _retained_units_count_by_kind(records, field::Symbol)
    counts = Dict{Symbol,Int}()
    for record in records
        value = getfield(record, field)
        counts[value] = get(counts, value, 0) + 1
    end
    return counts
end

function _synthetic_direct_contract(kind::Symbol, key::Symbol, source_cpb)
    return CTLForRetainedUnits.TerminalLoweringContract(
        Symbol(String(key), "_", String(kind)),
        key,
        :synthetic_direct_region,
        kind === :direct_core_identity_cpb ? :direct_core :
        kind === :direct_slab_identity_cpb ? :direct_midpoint_slab :
        :outer_mismatch_slab,
        kind,
        CRCForRetainedUnits.owned_cpb(
            source_cpb;
            support_kind = :synthetic_direct_owned_support,
            metadata = (; terminal_region_key = key),
        ),
        (source_cpb,),
        :direct_source_modes,
        :direct_or_trivial_embedding,
        :one_terminal_region,
        false,
        (; identity_like = true, source_cpb_equals_owned_support = true),
    )
end

function _synthetic_complete_shell_contracts()
    outer = CPBForRetainedUnits.filled_cpb(
        1:5,
        1:5,
        1:5;
        role = :synthetic_outer_box,
    )
    inner = CPBForRetainedUnits.filled_cpb(
        2:4,
        2:4,
        2:4;
        role = :synthetic_inner_box,
    )
    support = CRCForRetainedUnits.complete_shell_support(
        outer,
        inner;
        metadata = (; terminal_region_key = :synthetic_complete_shell),
    )
    strata = CPBForRetainedUnits.complete_shell_boundary_strata(outer, inner)
    lw_contract = CTLForRetainedUnits.TerminalLoweringContract(
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
            total_source_cpb_count = length(strata.all_strata),
            shell_support_count = strata.shell_support_count,
            stratum_support_count = strata.stratum_support_count,
        ),
    )
    pqs_contract = CTLForRetainedUnits.TerminalLoweringContract(
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
        (;
            q = 3,
            source_mode_shape = nothing,
            face_edge_corner_decomposition_required = false,
        ),
    )
    return (; lw_contract, pqs_contract)
end

function _synthetic_distorted_contract()
    source = CPBForRetainedUnits.cpb(
        1:7,
        1:3,
        1:3;
        role = :synthetic_distorted_product_source_cpb,
    )
    return CTLForRetainedUnits.TerminalLoweringContract(
        :synthetic_distorted_product_box_comx,
        :synthetic_distorted_product_box,
        :synthetic_distorted_region,
        :central_distorted_product_box,
        :distorted_product_box_comx,
        CRCForRetainedUnits.owned_cpb(
            source;
            support_kind = :synthetic_distorted_owned_support,
            metadata = (; terminal_region_key = :synthetic_distorted_product_box),
        ),
        (source,),
        :distorted_product_comx_all_axes,
        :distorted_product_realization_planned,
        :one_terminal_region,
        false,
        (; q = 3, L = 7, source_mode_shape = (7, 3, 3), aspect_ratio = 7 / 3, identity_like = false),
    )
end

function _synthetic_lowering_plan(policy, selected, policy_kind::Symbol)
    selected_tuple = Tuple(selected)
    return CTLForRetainedUnits.TerminalLoweringPlan(
        policy,
        selected_tuple,
        selected_tuple,
        (;
            object_kind = :synthetic_terminal_lowering_plan_summary,
            policy_kind,
            selected_contract_count = length(selected_tuple),
            materialized = false,
            retained_spaces_materialized = false,
            final_retained_units_materialized = false,
            pair_inventory_materialized = false,
            operator_blocks_materialized = false,
            hamiltonian_data_materialized = false,
        ),
        (; fixture = :cartesian_retained_units_contract),
    )
end

const _SYNTHETIC_DIRECT_CORE = _synthetic_direct_contract(
    :direct_core_identity_cpb,
    :synthetic_direct_core,
    CPBForRetainedUnits.filled_cpb(1:3, 1:3, 1:3; role = :synthetic_direct_core_cpb),
)
const _SYNTHETIC_DIRECT_SLAB = _synthetic_direct_contract(
    :direct_slab_identity_cpb,
    :synthetic_direct_slab,
    CPBForRetainedUnits.slab_cpb(1:3, 1:3, 4:4; role = :synthetic_direct_slab_cpb),
)
const _SYNTHETIC_DIRECT_BOUNDARY = _synthetic_direct_contract(
    :direct_boundary_slab_identity_cpb,
    :synthetic_boundary_slab,
    CPBForRetainedUnits.slab_cpb(1:3, 4:4, 1:3; role = :synthetic_boundary_slab_cpb),
)
const _SYNTHETIC_SHELL_CONTRACTS = _synthetic_complete_shell_contracts()
const _SYNTHETIC_DISTORTED = _synthetic_distorted_contract()

const _SYNTHETIC_WL_LOWERING_PLAN = _synthetic_lowering_plan(
    CTLForRetainedUnits.WhiteLindseyLowering(),
    (
        _SYNTHETIC_DIRECT_CORE,
        _SYNTHETIC_DIRECT_SLAB,
        _SYNTHETIC_DIRECT_BOUNDARY,
        _SYNTHETIC_SHELL_CONTRACTS.lw_contract,
        _SYNTHETIC_DISTORTED,
    ),
    :white_lindsey_terminal_lowering,
)

const _SYNTHETIC_PQS_LOWERING_PLAN = _synthetic_lowering_plan(
    CTLForRetainedUnits.PQSLowering(q = 3),
    (
        _SYNTHETIC_DIRECT_CORE,
        _SYNTHETIC_DIRECT_SLAB,
        _SYNTHETIC_DIRECT_BOUNDARY,
        _SYNTHETIC_SHELL_CONTRACTS.pqs_contract,
        _SYNTHETIC_DISTORTED,
    ),
    :pqs_terminal_lowering,
)

const _SYNTHETIC_WL_RETAINED_PLAN = CRU.retained_unit_plan(_SYNTHETIC_WL_LOWERING_PLAN)
const _SYNTHETIC_PQS_RETAINED_PLAN = CRU.retained_unit_plan(_SYNTHETIC_PQS_LOWERING_PLAN)

function _units_for_contract(plan, contract_key::Symbol)
    return Tuple(
        unit for unit in CRU.units(plan)
        if unit.source_contract_key == contract_key
    )
end

@testset "CartesianRetainedUnits White-Lindsey metadata plan" begin
    retained_plan = _SYNTHETIC_WL_RETAINED_PLAN
    retained_summary = CRU.summary(retained_plan)

    @test retained_summary.status == :available_retained_unit_plan
    @test retained_summary.lowering_policy_kind == :white_lindsey_terminal_lowering
    @test retained_summary.selected_contract_count == 5
    @test retained_summary.retained_unit_count == 30
    @test !retained_summary.materialized
    @test !retained_summary.transforms_materialized
    @test !retained_summary.coefficient_maps_materialized
    @test !retained_summary.pair_inventory_materialized
    @test !retained_summary.operator_blocks_materialized
    @test !retained_summary.hamiltonian_data_materialized

    for contract in CTLForRetainedUnits.selected_contracts(_SYNTHETIC_WL_LOWERING_PLAN)
        @test !isempty(_units_for_contract(retained_plan, contract.contract_key))
    end

    shell_units = _units_for_contract(
        retained_plan,
        _SYNTHETIC_SHELL_CONTRACTS.lw_contract.contract_key,
    )
    @test length(shell_units) == 26
    @test count(unit -> unit.metadata.stratum_kind == :facet_cpb, shell_units) == 6
    @test count(unit -> unit.metadata.stratum_kind == :edge_cpb, shell_units) == 12
    @test count(unit -> unit.metadata.stratum_kind == :corner_cpb, shell_units) == 8
    @test all(unit -> unit.source_cpb_index isa Int, shell_units)
    @test all(unit -> length(unit.source_cpbs) == 1, shell_units)
    @test all(
        unit -> unit.unit_kind == :white_lindsey_boundary_stratum_retained_unit,
        shell_units,
    )

    unit_kind_counts = _retained_units_count_by_kind(CRU.units(retained_plan), :unit_kind)
    @test unit_kind_counts[:direct_cpb_retained_unit] == 3
    @test unit_kind_counts[:white_lindsey_boundary_stratum_retained_unit] == 26
    @test unit_kind_counts[:distorted_product_box_retained_unit] == 1
    @test retained_summary.route_core_final_unit_available_count ==
          retained_summary.retained_unit_count
    @test retained_summary.route_core_final_unit_blocked_count == 0
    @test length(CRU.route_core_final_units(retained_plan)) ==
          retained_summary.retained_unit_count
end

@testset "CartesianRetainedUnits PQS metadata plan" begin
    retained_plan = _SYNTHETIC_PQS_RETAINED_PLAN
    retained_summary = CRU.summary(retained_plan)

    @test retained_summary.lowering_policy_kind == :pqs_terminal_lowering
    @test retained_summary.retained_unit_count == 5
    @test count(
        unit -> unit.unit_kind == :white_lindsey_boundary_stratum_retained_unit,
        CRU.units(retained_plan),
    ) == 0

    pqs_units = _units_for_contract(
        retained_plan,
        _SYNTHETIC_SHELL_CONTRACTS.pqs_contract.contract_key,
    )
    @test length(pqs_units) == 1
    unit = only(pqs_units)
    @test unit.unit_kind == :pqs_shell_retained_unit
    @test length(unit.source_cpbs) == 1
    @test CPBForRetainedUnits.codimension(only(unit.source_cpbs)) == 0
    @test unit.retained_rule == :pqs_boundary_comx_product_modes
    @test unit.realization_rule == :shell_projection_lowdin
    @test unit.dimension_status == :planned_not_materialized
    @test isnothing(unit.dimension)
    @test unit.metadata.face_edge_corner_expansion_used == false
    @test !isnothing(unit.route_core_final_unit)
    @test isnothing(unit.route_core_final_unit.intermediate.dimension)
    @test !retained_summary.pair_inventory_materialized
    @test !retained_summary.operator_blocks_materialized
end

@testset "CartesianRetainedUnits distorted product metadata" begin
    retained_plan = _SYNTHETIC_WL_RETAINED_PLAN
    distorted_units = Tuple(
        unit for unit in CRU.units(retained_plan)
        if unit.unit_kind == :distorted_product_box_retained_unit
    )

    @test length(distorted_units) == 1
    unit = only(distorted_units)
    @test unit.retained_rule == :distorted_product_comx_all_axes
    @test unit.realization_rule == :distorted_product_realization_planned
    @test unit.dimension_status == :planned_not_materialized
    @test isnothing(unit.dimension)
    @test unit.metadata.q == 3
    @test unit.metadata.L == 7
    @test unit.metadata.source_mode_shape == (7, 3, 3)
    @test unit.metadata.identity_like == false
end

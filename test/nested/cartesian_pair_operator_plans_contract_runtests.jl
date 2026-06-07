using Test
using GaussletBases

const CPOP = GaussletBases.CartesianPairOperatorPlans
const CUPForPairOps = GaussletBases.CartesianUnitPairs
const CRUForPairOps = GaussletBases.CartesianRetainedUnits
const CTLForPairOps = GaussletBases.CartesianTerminalLowering
const CRCForPairOps = GaussletBases.CartesianRouteCore
const CPBForPairOps = GaussletBases.CartesianCPB

function _pair_ops_direct_contract()
    source = CPBForPairOps.filled_cpb(
        1:3,
        1:3,
        1:3;
        role = :pair_ops_direct_core_cpb,
    )
    return CTLForPairOps.TerminalLoweringContract(
        :pair_ops_direct_core,
        :pair_ops_direct_core,
        :synthetic_direct_region,
        :direct_core,
        :direct_core_identity_cpb,
        CRCForPairOps.owned_cpb(
            source;
            support_kind = :pair_ops_direct_owned_support,
            metadata = (; terminal_region_key = :pair_ops_direct_core),
        ),
        (source,),
        :direct_source_modes,
        :direct_or_trivial_embedding,
        :one_terminal_region,
        false,
        (; identity_like = true),
    )
end

function _pair_ops_pqs_contract()
    outer = CPBForPairOps.filled_cpb(
        1:5,
        1:5,
        1:5;
        role = :pair_ops_pqs_outer_box,
    )
    inner = CPBForPairOps.filled_cpb(
        2:4,
        2:4,
        2:4;
        role = :pair_ops_pqs_inner_box,
    )
    support = CRCForPairOps.complete_shell_support(
        outer,
        inner;
        metadata = (; terminal_region_key = :pair_ops_pqs_shell),
    )
    return CTLForPairOps.TerminalLoweringContract(
        :pair_ops_pqs_shell,
        :pair_ops_pqs_shell,
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

function _pair_ops_white_lindsey_contract()
    outer = CPBForPairOps.filled_cpb(
        1:5,
        1:5,
        1:5;
        role = :pair_ops_lw_outer_box,
    )
    inner = CPBForPairOps.filled_cpb(
        2:4,
        2:4,
        2:4;
        role = :pair_ops_lw_inner_box,
    )
    support = CRCForPairOps.complete_shell_support(
        outer,
        inner;
        metadata = (; terminal_region_key = :pair_ops_lw_shell),
    )
    strata = CPBForPairOps.complete_shell_boundary_strata(outer, inner)
    return CTLForPairOps.TerminalLoweringContract(
        :pair_ops_lw_shell,
        :pair_ops_lw_shell,
        :synthetic_complete_shell_region,
        :complete_shell,
        :white_lindsey_boundary_strata,
        support,
        strata.all_strata,
        :white_lindsey_boundary_stratum_product,
        :direct_or_trivial_embedding,
        :white_lindsey_boundary_strata_children_planned,
        false,
        (; boundary_strata = true),
    )
end

function _pair_ops_lowering_plan()
    selected = (
        _pair_ops_direct_contract(),
        _pair_ops_pqs_contract(),
        _pair_ops_white_lindsey_contract(),
    )
    return CTLForPairOps.TerminalLoweringPlan(
        CTLForPairOps.PQSLowering(q = 3),
        selected,
        selected,
        (;
            object_kind = :synthetic_terminal_lowering_plan_summary,
            status = :available_terminal_lowering_plan,
            policy_kind = :pair_ops_synthetic_lowering,
            terminal_region_count = length(selected),
            available_contract_count = length(selected),
            selected_contract_count = length(selected),
            selected_contract_kinds =
                Tuple(contract.lowering_kind for contract in selected),
            all_terminal_regions_have_selected_contract = true,
            materialized = false,
            retained_spaces_materialized = false,
            final_retained_units_materialized = false,
            pair_inventory_materialized = false,
            operator_blocks_materialized = false,
            hamiltonian_data_materialized = false,
        ),
        (; fixture = :cartesian_pair_operator_plans_contract),
    )
end

function _pair_ops_count(counts, field::Symbol, value)
    matches = Tuple(entry for entry in counts if getproperty(entry, field) == value)
    isempty(matches) && return 0
    return only(matches).count
end

function _pair_ops_small_retained_plan()
    lowering_plan = _pair_ops_lowering_plan()
    direct = _pair_ops_direct_retained_unit()
    pqs = _pair_ops_pqs_retained_unit()
    lw = _pair_ops_white_lindsey_retained_unit()
    return CRUForPairOps.RetainedUnitPlan(
        CRUForPairOps.MetadataOnlyRetainedUnits(),
        lowering_plan,
        (direct, pqs, lw),
        (;
            object_kind = :synthetic_retained_unit_plan_summary,
            status = :available_retained_unit_plan,
            retained_unit_count = 3,
            materialized = false,
            transforms_materialized = false,
            coefficient_maps_materialized = false,
            pair_inventory_materialized = false,
            operator_blocks_materialized = false,
            hamiltonian_data_materialized = false,
        ),
        (; fixture_variant = :small_pair_operator_plan),
    )
end

function _pair_ops_direct_retained_unit()
    source = CPBForPairOps.filled_cpb(
        1:3,
        1:3,
        1:3;
        role = :pair_ops_direct_retained_source,
    )
    support = CRCForPairOps.owned_cpb(source)
    route_core_unit = _pair_ops_route_core_unit(
        :pair_ops_direct_retained_unit,
        :direct_cpb_retained_unit,
        :direct_identity_cpb,
        :identity_source_modes,
        :direct_or_trivial_embedding,
        support,
        (source,),
    )
    return _pair_ops_retained_unit_record(
        :pair_ops_direct_retained_unit,
        1,
        :direct_cpb_retained_unit,
        :direct_core_identity_cpb,
        :identity_source_modes,
        :direct_or_trivial_embedding,
        support,
        (source,),
        route_core_unit,
    )
end

function _pair_ops_pqs_retained_unit()
    outer = CPBForPairOps.filled_cpb(
        1:5,
        1:5,
        1:5;
        role = :pair_ops_pqs_retained_source,
    )
    inner = CPBForPairOps.filled_cpb(
        2:4,
        2:4,
        2:4;
        role = :pair_ops_pqs_retained_inner,
    )
    support = CRCForPairOps.complete_shell_support(outer, inner)
    route_core_unit = _pair_ops_route_core_unit(
        :pair_ops_pqs_retained_unit,
        :pqs_shell_retained_unit,
        :pqs_filled_source_cpb,
        :pqs_boundary_comx_product_modes,
        :shell_projection_lowdin,
        support,
        (outer,),
    )
    return _pair_ops_retained_unit_record(
        :pair_ops_pqs_retained_unit,
        2,
        :pqs_shell_retained_unit,
        :pqs_filled_source_cpb,
        :pqs_boundary_comx_product_modes,
        :shell_projection_lowdin,
        support,
        (outer,),
        route_core_unit,
    )
end

function _pair_ops_white_lindsey_retained_unit()
    source = CPBForPairOps.slab_cpb(
        1:1,
        2:4,
        2:4;
        role = :pair_ops_lw_retained_facet,
    )
    support = CRCForPairOps.owned_cpb(source)
    route_core_unit = _pair_ops_route_core_unit(
        :pair_ops_lw_retained_unit,
        :white_lindsey_boundary_stratum_retained_unit,
        :white_lindsey_boundary_strata,
        :white_lindsey_boundary_stratum_product,
        :direct_or_trivial_embedding,
        support,
        (source,),
    )
    return _pair_ops_retained_unit_record(
        :pair_ops_lw_retained_unit,
        3,
        :white_lindsey_boundary_stratum_retained_unit,
        :white_lindsey_boundary_strata,
        :white_lindsey_boundary_stratum_product,
        :direct_or_trivial_embedding,
        support,
        (source,),
        route_core_unit,
    )
end

function _pair_ops_route_core_unit(
    unit_key,
    unit_kind,
    lowering_recipe,
    retained_rule,
    realization_rule,
    owned_support,
    source_cpbs,
)
    region = CRCForPairOps.shellification_region(unit_kind, owned_support)
    lowering =
        lowering_recipe === :pqs_filled_source_cpb ?
        CRCForPairOps.pqs_filled_source_lowering(region, only(source_cpbs)) :
        CRCForPairOps.lowering_source(lowering_recipe, region, source_cpbs)
    intermediate = CRCForPairOps.intermediate_retained_space(
        lowering;
        retained_rule,
        materialized = false,
    )
    realization =
        realization_rule === :shell_projection_lowdin ?
        CRCForPairOps.pqs_shell_realization(intermediate, region) :
        CRCForPairOps.trivial_shell_realization(intermediate, region)
    return CRCForPairOps.final_retained_unit(
        unit_key,
        unit_kind,
        lowering,
        intermediate,
        realization,
    )
end

function _pair_ops_retained_unit_record(
    unit_key,
    unit_index,
    unit_kind,
    lowering_kind,
    retained_rule,
    realization_rule,
    owned_support,
    source_cpbs,
    route_core_unit,
)
    return CRUForPairOps.RetainedUnitRecord(
        unit_key,
        unit_index,
        unit_kind,
        Symbol(unit_key, "_contract"),
        Symbol(unit_key, "_terminal_region"),
        :synthetic_terminal_region,
        :synthetic_terminal_region,
        lowering_kind,
        retained_rule,
        realization_rule,
        owned_support,
        source_cpbs,
        nothing,
        :not_materialized,
        nothing,
        :not_materialized,
        nothing,
        route_core_unit,
        false,
        (; route_core_sidecar_status = :available),
    )
end

@testset "CartesianPairOperatorPlans unavailable summary" begin
    summary = CPOP.unavailable_summary(:not_selected, :not_selected_route)

    @test summary.object_kind == :cartesian_pair_operator_plan_summary
    @test summary.status == :not_selected
    @test summary.blocker == :not_selected_route
    @test summary.retained_unit_count == 0
    @test summary.unit_pair_count == 0
    @test summary.pair_operator_plan_count == 0
    @test summary.expected_pair_operator_plan_count == 0
    @test summary.pair_family_counts == ()
    @test summary.source_operator_path_counts == ()
    @test !summary.route_core_pair_operator_plan_inventory_available
    @test summary.route_core_pair_operator_plan_inventory_status == :not_available
    @test summary.route_core_pair_operator_plan_inventory_blocker ==
          :not_selected_route
    @test summary.route_core_pair_operator_plan_count == 0
    @test summary.route_core_pair_operator_plan_blocked_count == 0
    @test !summary.materialized
    @test !summary.source_operator_blocks_materialized
    @test !summary.final_pair_blocks_materialized
    @test !summary.operator_blocks_materialized
    @test !summary.hamiltonian_data_materialized
    @test !summary.artifacts_materialized
end

@testset "CartesianPairOperatorPlans metadata inventory" begin
    retained_plan = _pair_ops_small_retained_plan()
    pair_plan = CUPForPairOps.unit_pair_plan(retained_plan)
    operator_plan = CPOP.pair_operator_plan(pair_plan)
    operator_summary = CPOP.summary(operator_plan)

    @test operator_plan isa CPOP.PairOperatorPlan
    @test operator_summary.status == :available_pair_operator_plan
    @test operator_summary.retained_unit_count == 3
    @test operator_summary.unit_pair_count == 6
    @test operator_summary.pair_operator_plan_count ==
          operator_summary.unit_pair_count
    @test operator_summary.expected_pair_operator_plan_count ==
          operator_summary.unit_pair_count
    @test length(CPOP.pair_operator_records(operator_plan)) ==
          operator_summary.pair_operator_plan_count
    @test operator_summary.route_core_pair_operator_plan_inventory_available
    @test operator_summary.route_core_pair_operator_plan_count ==
          operator_summary.pair_operator_plan_count

    @test _pair_ops_count(
        operator_summary.source_operator_path_counts,
        :source_operator_path,
        :direct_identity_cpb_path,
    ) == 1
    @test _pair_ops_count(
        operator_summary.source_operator_path_counts,
        :source_operator_path,
        :pqs_source_cpb_1d_factor_path,
    ) == 3
    @test _pair_ops_count(
        operator_summary.source_operator_path_counts,
        :source_operator_path,
        :white_lindsey_boundary_stratum_adapter_path,
    ) == 2
    @test _pair_ops_count(
        operator_summary.final_block_path_counts,
        :final_block_path,
        :source_block_realization_then_final_block,
    ) == 3
    @test _pair_ops_count(
        operator_summary.final_block_path_counts,
        :final_block_path,
        :source_block_direct_to_final_block,
    ) == 3
    @test operator_summary.blocked_pair_operator_plan_count == 0
    @test _pair_ops_count(
        operator_summary.materialization_status_counts,
        :materialization_status,
        :metadata_only_not_materialized,
    ) == operator_summary.pair_operator_plan_count
    @test !operator_summary.materialized
    @test !operator_summary.source_operator_blocks_materialized
    @test !operator_summary.final_pair_blocks_materialized
    @test !operator_summary.operator_blocks_materialized
    @test !operator_summary.hamiltonian_data_materialized
    @test !operator_summary.artifacts_materialized
end

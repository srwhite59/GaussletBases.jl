using Test
using GaussletBases

const CTLForUnitPairs = GaussletBases.CartesianTerminalLowering
const CRUForUnitPairs = GaussletBases.CartesianRetainedUnits
const CUP = GaussletBases.CartesianUnitPairs
const CRCForUnitPairs = GaussletBases.CartesianRouteCore
const CPBForUnitPairs = GaussletBases.CartesianCPB

function _unit_pairs_direct_contract()
    source = CPBForUnitPairs.filled_cpb(
        1:3,
        1:3,
        1:3;
        role = :unit_pairs_direct_core_cpb,
    )
    return CTLForUnitPairs.TerminalLoweringContract(
        :unit_pairs_direct_core,
        :unit_pairs_direct_core,
        :synthetic_direct_region,
        :direct_core,
        :direct_core_identity_cpb,
        CRCForUnitPairs.owned_cpb(
            source;
            support_kind = :unit_pairs_direct_owned_support,
            metadata = (; terminal_region_key = :unit_pairs_direct_core),
        ),
        (source,),
        :direct_source_modes,
        :direct_or_trivial_embedding,
        :one_terminal_region,
        false,
        (; identity_like = true),
    )
end

function _unit_pairs_pqs_contract()
    outer = CPBForUnitPairs.filled_cpb(
        1:5,
        1:5,
        1:5;
        role = :unit_pairs_pqs_outer_box,
    )
    inner = CPBForUnitPairs.filled_cpb(
        2:4,
        2:4,
        2:4;
        role = :unit_pairs_pqs_inner_box,
    )
    support = CRCForUnitPairs.complete_shell_support(
        outer,
        inner;
        metadata = (; terminal_region_key = :unit_pairs_pqs_shell),
    )
    return CTLForUnitPairs.TerminalLoweringContract(
        :unit_pairs_pqs_shell,
        :unit_pairs_pqs_shell,
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

function _unit_pairs_distorted_contract()
    source = CPBForUnitPairs.cpb(
        1:7,
        1:3,
        1:3;
        role = :unit_pairs_distorted_source_cpb,
    )
    return CTLForUnitPairs.TerminalLoweringContract(
        :unit_pairs_distorted_box,
        :unit_pairs_distorted_box,
        :synthetic_distorted_region,
        :central_distorted_product_box,
        :distorted_product_box_comx,
        CRCForUnitPairs.owned_cpb(
            source;
            support_kind = :unit_pairs_distorted_owned_support,
            metadata = (; terminal_region_key = :unit_pairs_distorted_box),
        ),
        (source,),
        :distorted_product_comx_all_axes,
        :distorted_product_realization_planned,
        :one_terminal_region,
        false,
        (; q = 3, L = 7, source_mode_shape = (7, 3, 3), aspect_ratio = 7 / 3),
    )
end

function _unit_pairs_lowering_plan()
    selected = (
        _unit_pairs_direct_contract(),
        _unit_pairs_pqs_contract(),
        _unit_pairs_distorted_contract(),
    )
    return CTLForUnitPairs.TerminalLoweringPlan(
        CTLForUnitPairs.PQSLowering(q = 3),
        selected,
        selected,
        (;
            object_kind = :synthetic_terminal_lowering_plan_summary,
            status = :available_terminal_lowering_plan,
            policy_kind = :pqs_terminal_lowering,
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
        (; fixture = :cartesian_unit_pairs_contract),
    )
end

function _unit_pairs_blocked_retained_plan(plan)
    units = CRUForUnitPairs.units(plan)
    first_unit = first(units)
    blocked_first_unit = CRUForUnitPairs.RetainedUnitRecord(
        first_unit.unit_key,
        first_unit.unit_index,
        first_unit.unit_kind,
        first_unit.source_contract_key,
        first_unit.terminal_region_key,
        first_unit.terminal_region_role,
        first_unit.terminal_region_kind,
        first_unit.lowering_kind,
        first_unit.retained_rule,
        first_unit.realization_rule,
        first_unit.owned_support,
        first_unit.source_cpbs,
        first_unit.source_cpb_index,
        first_unit.dimension_status,
        first_unit.dimension,
        first_unit.column_range_status,
        first_unit.column_range,
        nothing,
        first_unit.materialized,
        merge(
            first_unit.metadata,
            (;
                route_core_sidecar_status = :blocked,
                route_core_sidecar_blocker = :synthetic_missing_sidecar,
            ),
        ),
    )
    blocked_units = CRUForUnitPairs.RetainedUnitRecord[blocked_first_unit]
    append!(blocked_units, @view units[2:end])
    return CRUForUnitPairs.RetainedUnitPlan(
        plan.policy,
        plan.lowering_plan,
        blocked_units,
        plan.summary,
        merge(plan.metadata, (; fixture_variant = :missing_route_core_sidecar)),
    )
end

function _unit_pair_family_count(summary, pair_family::Symbol)
    for entry in summary.pair_family_counts
        entry.pair_family == pair_family && return entry.count
    end
    return 0
end

@testset "CartesianUnitPairs unavailable summary" begin
    summary = CUP.unavailable_summary(:not_selected, :not_selected_route)

    @test summary.object_kind == :cartesian_unit_pair_plan_summary
    @test summary.status == :not_selected
    @test summary.blocker == :not_selected_route
    @test summary.retained_unit_count == 0
    @test summary.pair_count == 0
    @test summary.expected_upper_triangular_pair_count == 0
    @test summary.pair_family_counts == ()
    @test !summary.route_core_pair_inventory_available
    @test summary.route_core_pair_inventory_status == :not_available
    @test summary.route_core_pair_inventory_blocker == :not_selected_route
    @test summary.route_core_pair_count == 0
    @test summary.route_core_pair_missing_final_unit_indices == ()
    @test !summary.materialized
    @test !summary.pair_inventory_materialized
    @test !summary.source_operator_blocks_materialized
    @test !summary.operator_blocks_materialized
    @test !summary.hamiltonian_data_materialized
    @test !summary.artifacts_materialized
end

@testset "CartesianUnitPairs metadata inventory" begin
    retained_plan = CRUForUnitPairs.retained_unit_plan(_unit_pairs_lowering_plan())
    pair_plan = CUP.unit_pair_plan(retained_plan)
    pair_summary = CUP.summary(pair_plan)

    @test pair_summary.status == :available_unit_pair_plan
    @test pair_summary.retained_unit_count == 3
    @test pair_summary.pair_count == 6
    @test pair_summary.expected_upper_triangular_pair_count == 6
    @test length(CUP.unit_pairs(pair_plan)) == 6
    @test Tuple(pair.pair_index for pair in CUP.unit_pairs(pair_plan)) == (1, 2, 3, 4, 5, 6)
    @test Tuple(pair.left_index => pair.right_index for pair in CUP.unit_pairs(pair_plan)) ==
          (1 => 1, 1 => 2, 1 => 3, 2 => 2, 2 => 3, 3 => 3)
    pair_index_table = CUP.unit_pair_index_table(retained_plan)
    @test length(pair_index_table) == 6
    @test Tuple(pair.left_index => pair.right_index for pair in pair_index_table) ==
          (1 => 1, 1 => 2, 1 => 3, 2 => 2, 2 => 3, 3 => 3)
    @test all(pair -> isnothing(pair.route_core_pair_sidecar), pair_index_table)
    @test pair_summary.route_core_pair_inventory_available
    @test pair_summary.route_core_pair_inventory_status ==
          :available_route_core_pair_inventory
    @test pair_summary.route_core_pair_count == pair_summary.pair_count
    @test pair_summary.pair_storage == :unit_pair_index_table
    @test !pair_summary.rich_unit_pair_records_stored
    @test !isnothing(CUP.route_core_pair_inventory(pair_plan))
    @test all(pair -> isnothing(pair.route_core_pair_sidecar), CUP.unit_pairs(pair_plan))
    @test all(
        pair -> pair.metadata.rich_unit_pair_record_stored === false,
        CUP.unit_pairs(pair_plan),
    )
    @test _unit_pair_family_count(
        pair_summary,
        :direct_cpb_retained_unit__direct_cpb_retained_unit,
    ) == 1
    @test _unit_pair_family_count(
        pair_summary,
        :direct_cpb_retained_unit__pqs_shell_retained_unit,
    ) == 1
    @test _unit_pair_family_count(
        pair_summary,
        :pqs_shell_retained_unit__distorted_product_box_retained_unit,
    ) == 1
    @test !pair_summary.materialized
    @test !pair_summary.pair_inventory_materialized
    @test !pair_summary.source_operator_blocks_materialized
    @test !pair_summary.operator_blocks_materialized
    @test !pair_summary.hamiltonian_data_materialized
    @test !pair_summary.artifacts_materialized
end

@testset "CartesianUnitPairs missing RouteCore sidecars" begin
    retained_plan = CRUForUnitPairs.retained_unit_plan(_unit_pairs_lowering_plan())
    blocked_plan = _unit_pairs_blocked_retained_plan(retained_plan)
    pair_plan = CUP.unit_pair_plan(blocked_plan)
    pair_summary = CUP.summary(pair_plan)

    @test pair_summary.retained_unit_count == 3
    @test pair_summary.pair_count == 6
    @test !pair_summary.route_core_pair_inventory_available
    @test pair_summary.route_core_pair_inventory_status ==
          :blocked_missing_route_core_final_units
    @test pair_summary.route_core_pair_inventory_blocker ==
          :missing_route_core_final_units
    @test pair_summary.route_core_pair_count == 0
    @test pair_summary.route_core_pair_missing_final_unit_indices == (1,)
    @test isnothing(CUP.route_core_pair_inventory(pair_plan))
    @test all(pair -> isnothing(pair.route_core_pair_sidecar), CUP.unit_pairs(pair_plan))
    @test !pair_summary.operator_blocks_materialized
    @test !pair_summary.hamiltonian_data_materialized
end

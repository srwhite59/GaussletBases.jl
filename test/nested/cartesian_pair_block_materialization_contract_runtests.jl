using Test
using GaussletBases

const CPBM = GaussletBases.CartesianPairBlockMaterialization
const CPOPForPairBlocks = GaussletBases.CartesianPairOperatorPlans
const CUPForPairBlocks = GaussletBases.CartesianUnitPairs
const CRTCForPairBlocks = GaussletBases.CartesianRetainedUnitTransformContracts
const CRUForPairBlocks = GaussletBases.CartesianRetainedUnits
const CTLForPairBlocks = GaussletBases.CartesianTerminalLowering

function _pair_block_count(counts, field::Symbol, value)
    matches = Tuple(entry for entry in counts if getproperty(entry, field) == value)
    isempty(matches) && return 0
    return only(matches).count
end

function _pair_block_minimal_lowering_plan()
    return CTLForPairBlocks.TerminalLoweringPlan(
        CTLForPairBlocks.PQSLowering(q = 3),
        (),
        (),
        (;
            object_kind = :synthetic_terminal_lowering_plan_summary,
            status = :available_terminal_lowering_plan,
            policy_kind = :pair_block_materialization_synthetic_lowering,
            terminal_region_count = 0,
            available_contract_count = 0,
            selected_contract_count = 0,
            selected_contract_kinds = (),
            all_terminal_regions_have_selected_contract = true,
            materialized = false,
            retained_spaces_materialized = false,
            final_retained_units_materialized = false,
            pair_inventory_materialized = false,
            operator_blocks_materialized = false,
            hamiltonian_data_materialized = false,
        ),
        (; fixture = :cartesian_pair_block_materialization_contract),
    )
end

function _pair_block_retained_unit(
    unit_key::Symbol,
    unit_index::Int,
    unit_kind::Symbol,
    lowering_kind::Symbol,
    retained_rule::Symbol,
    realization_rule,
)
    return CRUForPairBlocks.RetainedUnitRecord(
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
        nothing,
        (),
        nothing,
        :not_materialized,
        nothing,
        :not_materialized,
        nothing,
        nothing,
        false,
        (; route_core_sidecar_status = :not_materialized),
    )
end

function _pair_block_retained_plan()
    direct = _pair_block_retained_unit(
        :pair_block_direct_unit,
        1,
        :direct_cpb_retained_unit,
        :direct_core_identity_cpb,
        :direct_source_modes,
        :direct_or_trivial_embedding,
    )
    pqs = _pair_block_retained_unit(
        :pair_block_pqs_unit,
        2,
        :pqs_shell_retained_unit,
        :pqs_filled_source_cpb,
        :pqs_boundary_comx_product_modes,
        :shell_projection_lowdin,
    )
    units = (direct, pqs)
    return CRUForPairBlocks.RetainedUnitPlan(
        CRUForPairBlocks.MetadataOnlyRetainedUnits(),
        _pair_block_minimal_lowering_plan(),
        units,
        (;
            object_kind = :synthetic_retained_unit_plan_summary,
            status = :available_retained_unit_plan,
            retained_unit_count = length(units),
            materialized = false,
            transforms_materialized = false,
            coefficient_maps_materialized = false,
            pair_inventory_materialized = false,
            operator_blocks_materialized = false,
            hamiltonian_data_materialized = false,
        ),
        (; fixture = :cartesian_pair_block_materialization_contract),
    )
end

function _pair_block_record_for(plan, left_key::Symbol, right_key::Symbol)
    return only(
        record for record in CPBM.pair_block_materialization_records(plan)
        if record.pair_key == (left_key, right_key)
    )
end

@testset "CartesianPairBlockMaterialization unavailable summary" begin
    summary = CPBM.unavailable_summary(:not_selected, :not_selected_route)

    @test summary.object_kind == :cartesian_pair_block_materialization_plan_summary
    @test summary.status == :not_selected
    @test summary.blocker == :not_selected_route
    @test summary.pair_operator_plan_count == 0
    @test summary.pair_block_record_count == 0
    @test summary.ready_record_count == 0
    @test summary.blocked_record_count == 0
    @test summary.materialization_path_counts == ()
    @test summary.readiness_status_counts == ()
    @test summary.blocker_counts == ()
    @test !summary.materialized
    @test !summary.source_operator_blocks_materialized
    @test !summary.final_pair_blocks_materialized
    @test !summary.operator_blocks_materialized
    @test !summary.hamiltonian_data_materialized
    @test !summary.artifacts_materialized
end

@testset "CartesianPairBlockMaterialization direct/direct preflight" begin
    retained_plan = _pair_block_retained_plan()
    unit_pair_plan = CUPForPairBlocks.unit_pair_plan(retained_plan)
    transform_plan =
        CRTCForPairBlocks.retained_unit_transform_contract_plan(retained_plan)
    pair_operator_plan =
        CPOPForPairBlocks.pair_operator_plan(
            unit_pair_plan,
            transform_plan;
            route_core_sidecars = false,
        )
    materialization_plan =
        CPBM.pair_block_materialization_plan(pair_operator_plan)
    materialization_summary = CPBM.summary(materialization_plan)

    @test materialization_plan isa CPBM.PairBlockMaterializationPlan
    @test length(CPBM.pair_block_materialization_records(materialization_plan)) == 3
    @test materialization_summary.status == :blocked_pair_block_materialization_plan
    @test materialization_summary.pair_operator_plan_count == 3
    @test materialization_summary.pair_block_record_count == 3
    @test materialization_summary.ready_record_count == 1
    @test materialization_summary.blocked_record_count == 2
    @test _pair_block_count(
        materialization_summary.materialization_path_counts,
        :materialization_path,
        :direct_direct_pair_block_materialization_pilot,
    ) == 1
    @test _pair_block_count(
        materialization_summary.materialization_path_counts,
        :materialization_path,
        :deferred_pair_block_materialization_path,
    ) == 2
    @test _pair_block_count(
        materialization_summary.readiness_status_counts,
        :readiness_status,
        :ready_metadata_only_not_materialized,
    ) == 1
    @test _pair_block_count(
        materialization_summary.readiness_status_counts,
        :readiness_status,
        :blocked_pair_block_materialization_not_implemented,
    ) == 2
    @test _pair_block_count(
        materialization_summary.blocker_counts,
        :blocker,
        :non_direct_direct_pair_block_materialization_not_implemented,
    ) == 2

    direct_record = _pair_block_record_for(
        materialization_plan,
        :pair_block_direct_unit,
        :pair_block_direct_unit,
    )
    direct_pqs_record = _pair_block_record_for(
        materialization_plan,
        :pair_block_direct_unit,
        :pair_block_pqs_unit,
    )
    pqs_pqs_record = _pair_block_record_for(
        materialization_plan,
        :pair_block_pqs_unit,
        :pair_block_pqs_unit,
    )

    @test direct_record.materialization_path ==
          :direct_direct_pair_block_materialization_pilot
    @test direct_record.readiness_status == :ready_metadata_only_not_materialized
    @test direct_record.blocker === nothing
    @test !direct_record.materialized
    @test direct_pqs_record.readiness_status ==
          :blocked_pair_block_materialization_not_implemented
    @test direct_pqs_record.blocker ==
          :non_direct_direct_pair_block_materialization_not_implemented
    @test pqs_pqs_record.readiness_status ==
          :blocked_pair_block_materialization_not_implemented
    @test pqs_pqs_record.blocker ==
          :non_direct_direct_pair_block_materialization_not_implemented
    @test !materialization_summary.materialized
    @test !materialization_summary.source_operator_blocks_materialized
    @test !materialization_summary.final_pair_blocks_materialized
    @test !materialization_summary.operator_blocks_materialized
    @test !materialization_summary.hamiltonian_data_materialized
    @test !materialization_summary.artifacts_materialized
end

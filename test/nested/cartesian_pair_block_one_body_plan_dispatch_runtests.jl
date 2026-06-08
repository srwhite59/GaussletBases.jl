using Test
using GaussletBases

const CPBMPlanDispatch = GaussletBases.CartesianPairBlockMaterialization
const CTLPlanDispatch = GaussletBases.CartesianTerminalLowering
const CRUPlanDispatch = GaussletBases.CartesianRetainedUnits
const CRTCPlanDispatch = GaussletBases.CartesianRetainedUnitTransformContracts
const CUPPlanDispatch = GaussletBases.CartesianUnitPairs
const CPOPPlanDispatch = GaussletBases.CartesianPairOperatorPlans

function _plan_dispatch_count(counts, field::Symbol, value)
    matches = Tuple(entry for entry in counts if getproperty(entry, field) == value)
    isempty(matches) && return 0
    return only(matches).count
end

function _plan_dispatch_lowering_plan()
    return CTLPlanDispatch.TerminalLoweringPlan(
        CTLPlanDispatch.PQSLowering(q = 3),
        (),
        (),
        (;
            object_kind = :synthetic_plan_dispatch_lowering_summary,
            status = :available_terminal_lowering_plan,
            materialized = false,
        ),
        (; fixture = :one_body_plan_dispatch),
    )
end

function _plan_dispatch_retained_unit(unit_key::Symbol, unit_index::Int, unit_kind::Symbol)
    return CRUPlanDispatch.RetainedUnitRecord(
        unit_key,
        unit_index,
        unit_kind,
        Symbol(unit_key, "_contract"),
        Symbol(unit_key, "_terminal_region"),
        :synthetic_terminal_region,
        :synthetic_terminal_region,
        :synthetic_lowering,
        :synthetic_retained_rule,
        :synthetic_realization,
        nothing,
        (),
        nothing,
        :not_materialized,
        nothing,
        :not_materialized,
        nothing,
        nothing,
        false,
        (; fixture = :one_body_plan_dispatch),
    )
end

function _plan_dispatch_unit_pair(
    pair_key::Tuple{Symbol,Symbol},
    pair_index::Int,
    pair_family::Symbol,
    units_by_key,
)
    left_unit = units_by_key[pair_key[1]]
    right_unit = units_by_key[pair_key[2]]
    return CUPPlanDispatch.UnitPairRecord(
        pair_key,
        pair_index,
        pair_family,
        left_unit,
        right_unit,
        left_unit.unit_index,
        right_unit.unit_index,
        left_unit.unit_key,
        right_unit.unit_key,
        left_unit.unit_kind,
        right_unit.unit_kind,
        nothing,
        false,
        (; fixture = :one_body_plan_dispatch),
    )
end

function _plan_dispatch_materialization_record(
    pair_key::Tuple{Symbol,Symbol},
    pair_index::Int,
    pair_family::Symbol,
    materialization_path::Symbol;
    readiness_status = :ready_metadata_only_not_materialized,
    blocker = nothing,
)
    return CPBMPlanDispatch.PairBlockMaterializationRecord(
        pair_key,
        pair_index,
        pair_family,
        :synthetic_source_operator_path,
        (; left = :synthetic_left_transform, right = :synthetic_right_transform),
        (; left = :synthetic_left_realization, right = :synthetic_right_realization),
        :synthetic_final_block_path,
        (
            :overlap,
            :position_x,
            :position_y,
            :position_z,
            :x2_x,
            :x2_y,
            :x2_z,
            :kinetic,
        ),
        materialization_path,
        readiness_status,
        blocker,
        false,
        (; fixture = :one_body_plan_dispatch),
    )
end

function _plan_dispatch_materialization_plan(; include_lw_unit_pair::Bool = true)
    units = (
        _plan_dispatch_retained_unit(:direct_left, 1, :direct_unit),
        _plan_dispatch_retained_unit(:direct_right, 2, :direct_unit),
        _plan_dispatch_retained_unit(:pqs_left, 3, :pqs_unit),
        _plan_dispatch_retained_unit(:pqs_right, 4, :pqs_unit),
        _plan_dispatch_retained_unit(:lw_left, 5, :white_lindsey_unit),
        _plan_dispatch_retained_unit(:lw_right, 6, :white_lindsey_unit),
    )
    units_by_key = Dict(unit.unit_key => unit for unit in units)
    retained_plan = CRUPlanDispatch.RetainedUnitPlan(
        CRUPlanDispatch.MetadataOnlyRetainedUnits(),
        _plan_dispatch_lowering_plan(),
        units,
        (;
            object_kind = :synthetic_plan_dispatch_retained_summary,
            status = :available_retained_unit_plan,
            materialized = false,
        ),
        (; fixture = :one_body_plan_dispatch),
    )
    pair_specs = [
        ((:direct_left, :direct_right), :direct_direct),
        ((:pqs_left, :pqs_right), :pqs_pqs),
    ]
    include_lw_unit_pair &&
        push!(pair_specs, ((:lw_left, :lw_right), :white_lindsey_boundary_stratum))
    unit_pairs = Tuple(
        _plan_dispatch_unit_pair(pair_key, index, pair_family, units_by_key)
        for (index, (pair_key, pair_family)) in pairs(pair_specs)
    )
    unit_pair_plan = CUPPlanDispatch.UnitPairPlan(
        CUPPlanDispatch.MetadataOnlyUnitPairs(),
        retained_plan,
        unit_pairs,
        nothing,
        (;
            object_kind = :synthetic_plan_dispatch_unit_pair_summary,
            status = :available_unit_pair_plan,
            materialized = false,
        ),
        (; fixture = :one_body_plan_dispatch),
    )
    transform_plan = CRTCPlanDispatch.RetainedUnitTransformContractPlan(
        CRTCPlanDispatch.MetadataOnlyRetainedUnitTransformContracts(),
        retained_plan,
        (),
        (;
            object_kind = :synthetic_plan_dispatch_transform_summary,
            status = :available_transform_contract_plan,
            materialized = false,
        ),
        (; fixture = :one_body_plan_dispatch),
    )
    pair_operator_plan = CPOPPlanDispatch.PairOperatorPlan(
        CPOPPlanDispatch.MetadataOnlyPairOperatorPlans(),
        unit_pair_plan,
        transform_plan,
        (),
        nothing,
        (;
            object_kind = :synthetic_plan_dispatch_pair_operator_summary,
            status = :available_pair_operator_plan,
            materialized = false,
        ),
        (; fixture = :one_body_plan_dispatch),
    )
    records = (
        _plan_dispatch_materialization_record(
            (:direct_left, :direct_right),
            1,
            :direct_direct,
            :direct_direct_pair_block_materialization_pilot,
        ),
        _plan_dispatch_materialization_record(
            (:pqs_left, :pqs_right),
            2,
            :pqs_pqs,
            :pqs_source_pair_preflight,
        ),
        _plan_dispatch_materialization_record(
            (:lw_left, :lw_right),
            3,
            :white_lindsey_boundary_stratum,
            :white_lindsey_boundary_stratum_adapter_preflight,
        ),
    )
    return CPBMPlanDispatch.PairBlockMaterializationPlan(
        CPBMPlanDispatch.MetadataOnlyPairBlockMaterialization(),
        pair_operator_plan,
        records,
        (;
            object_kind = :synthetic_plan_dispatch_materialization_summary,
            status = :available_pair_block_materialization_plan,
            materialized = false,
        ),
        (; fixture = :one_body_plan_dispatch),
    )
end

function _check_plan_dispatch_nonmaterialization_flags(summary)
    @test !summary.factors_constructed
    @test !summary.numerical_blocks_materialized
    @test !summary.source_operator_blocks_materialized
    @test !summary.final_pair_blocks_materialized
    @test !summary.operator_blocks_materialized
    @test !summary.hamiltonian_data_materialized
    @test !summary.artifacts_materialized
    @test !summary.mixed_dispatcher_materialized
    @test !summary.route_driver_wiring
    @test !summary.pqs_lowdin_materialized
    @test !summary.full_white_lindsey_route_assembled
end

@testset "CartesianPairBlockMaterialization one-body plan dispatch ready summary" begin
    plan = _plan_dispatch_materialization_plan()
    summary = CPBMPlanDispatch._one_body_pair_block_plan_dispatch_summary(
        plan,
        :overlap;
        inputs = (; parent_axis_counts = (3, 4, 5), overlap_1d = :overlap),
    )

    @test summary.object_kind ==
          :cartesian_pair_block_mixed_one_body_plan_dispatch_summary
    @test summary.status == :ready_metadata_only_not_materialized
    @test isnothing(summary.blocker)
    @test summary.requested_term == :overlap
    @test summary.term_descriptor.family == :overlap
    @test summary.record_count == 3
    @test summary.ready_count == 3
    @test summary.blocked_count == 0
    @test length(summary.record_summaries) == 3
    @test _plan_dispatch_count(
        summary.selector_family_counts,
        :selector_family,
        :direct_direct,
    ) == 1
    @test _plan_dispatch_count(
        summary.selector_family_counts,
        :selector_family,
        :pqs_source_pair,
    ) == 1
    @test _plan_dispatch_count(
        summary.selector_family_counts,
        :selector_family,
        :white_lindsey_boundary_stratum,
    ) == 1
    @test _plan_dispatch_count(
        summary.dispatch_status_counts,
        :status,
        :ready_metadata_only_not_materialized,
    ) == 3
    @test summary.blocker_counts == ()
    @test _plan_dispatch_count(
        summary.record_materialization_path_counts,
        :record_materialization_path,
        :white_lindsey_boundary_stratum_adapter_preflight,
    ) == 1
    @test _plan_dispatch_count(
        summary.materialization_path_counts,
        :materialization_path,
        :white_lindsey_boundary_stratum_one_body_selector,
    ) == 1
    @test summary.white_lindsey_unit_pair_required_count == 1
    @test summary.white_lindsey_unit_pair_available_count == 1
    _check_plan_dispatch_nonmaterialization_flags(summary)
end

@testset "CartesianPairBlockMaterialization one-body plan dispatch blockers" begin
    plan = _plan_dispatch_materialization_plan()
    missing_inputs = CPBMPlanDispatch._one_body_pair_block_plan_dispatch_summary(
        plan,
        :kinetic;
        inputs = (; overlap_1d = :overlap),
    )

    @test missing_inputs.status == :blocked_mixed_one_body_plan_dispatch
    @test missing_inputs.blocker == :blocked_mixed_one_body_dispatch_records
    @test missing_inputs.record_count == 3
    @test missing_inputs.ready_count == 0
    @test missing_inputs.blocked_count == 3
    @test _plan_dispatch_count(
        missing_inputs.blocker_counts,
        :blocker,
        :missing_one_body_factor_inputs,
    ) == 3
    @test _plan_dispatch_count(
        missing_inputs.selector_family_counts,
        :selector_family,
        :white_lindsey_boundary_stratum,
    ) == 1
    @test missing_inputs.white_lindsey_unit_pair_required_count == 1
    @test missing_inputs.white_lindsey_unit_pair_available_count == 1
    _check_plan_dispatch_nonmaterialization_flags(missing_inputs)

    missing_lw_pair_plan =
        _plan_dispatch_materialization_plan(; include_lw_unit_pair = false)
    missing_lw_pair =
        CPBMPlanDispatch._one_body_pair_block_plan_dispatch_summary(
            missing_lw_pair_plan,
            :overlap;
            inputs = (; parent_axis_counts = (3, 4, 5), overlap_1d = :overlap),
        )
    @test missing_lw_pair.status == :partially_ready_mixed_one_body_plan_dispatch
    @test missing_lw_pair.record_count == 3
    @test missing_lw_pair.ready_count == 2
    @test missing_lw_pair.blocked_count == 1
    @test _plan_dispatch_count(
        missing_lw_pair.blocker_counts,
        :blocker,
        :missing_white_lindsey_unit_pair,
    ) == 1
    @test missing_lw_pair.white_lindsey_unit_pair_required_count == 1
    @test missing_lw_pair.white_lindsey_unit_pair_available_count == 0
    _check_plan_dispatch_nonmaterialization_flags(missing_lw_pair)

    @test_throws ArgumentError CPBMPlanDispatch._one_body_pair_block_plan_dispatch_summary(
        (;),
        :overlap;
        inputs = (; overlap_1d = :overlap),
    )
end

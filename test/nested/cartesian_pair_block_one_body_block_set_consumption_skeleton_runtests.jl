# Runtime role: tiny contract.
#
# Deferred metadata-only block-set consumption skeleton. This does not call the
# one-term numerical consumer, construct factors, or store batch matrices.

using Test
using GaussletBases

const CPBMSkeleton = GaussletBases.CartesianPairBlockMaterialization
const CTLSkeleton = GaussletBases.CartesianTerminalLowering
const CRUSkeleton = GaussletBases.CartesianRetainedUnits
const CRTCSkeleton = GaussletBases.CartesianRetainedUnitTransformContracts
const CUPSkeleton = GaussletBases.CartesianUnitPairs
const CPOPSkeleton = GaussletBases.CartesianPairOperatorPlans

function _skeleton_pair_operator_plan()
    lowering_plan = CTLSkeleton.TerminalLoweringPlan(
        CTLSkeleton.PQSLowering(q = 2),
        (),
        (),
        (; status = :available_terminal_lowering_plan, materialized = false),
        (; fixture = :one_body_block_set_consumption_skeleton),
    )
    retained_plan = CRUSkeleton.RetainedUnitPlan(
        CRUSkeleton.MetadataOnlyRetainedUnits(),
        lowering_plan,
        (),
        (; status = :available_retained_unit_plan, materialized = false),
        (; fixture = :one_body_block_set_consumption_skeleton),
    )
    unit_pair_plan = CUPSkeleton.UnitPairPlan(
        CUPSkeleton.MetadataOnlyUnitPairs(),
        retained_plan,
        (),
        nothing,
        (; status = :available_unit_pair_plan, materialized = false),
        (; fixture = :one_body_block_set_consumption_skeleton),
    )
    transform_plan = CRTCSkeleton.RetainedUnitTransformContractPlan(
        CRTCSkeleton.MetadataOnlyRetainedUnitTransformContracts(),
        retained_plan,
        (),
        (; status = :available_transform_contract_plan, materialized = false),
        (; fixture = :one_body_block_set_consumption_skeleton),
    )
    return CPOPSkeleton.PairOperatorPlan(
        CPOPSkeleton.MetadataOnlyPairOperatorPlans(),
        unit_pair_plan,
        transform_plan,
        (),
        nothing,
        (; status = :available_pair_operator_plan, materialized = false),
        (; fixture = :one_body_block_set_consumption_skeleton),
    )
end

function _skeleton_record(
    pair_key::Tuple{Symbol,Symbol},
    pair_index::Int,
    pair_family::Symbol,
    materialization_path::Symbol;
    supported_terms = (:overlap, :kinetic),
)
    return CPBMSkeleton.PairBlockMaterializationRecord(
        pair_key,
        pair_index,
        pair_family,
        :synthetic_source_operator_path,
        (; left = :synthetic_left_transform, right = :synthetic_right_transform),
        (; left = :synthetic_left_realization, right = :synthetic_right_realization),
        :synthetic_final_block_path,
        supported_terms,
        materialization_path,
        :ready_metadata_only_not_materialized,
        nothing,
        false,
        (; fixture = :one_body_block_set_consumption_skeleton),
    )
end

function _skeleton_plan()
    records = (
        _skeleton_record(
            (:direct_left, :direct_right),
            1,
            :direct_direct,
            :direct_direct_pair_block_materialization_pilot,
        ),
        _skeleton_record(
            (:pqs_left, :pqs_right),
            2,
            :pqs_pqs,
            :pqs_source_pair_preflight,
        ),
        _skeleton_record(
            (:unsupported_left, :unsupported_right),
            3,
            :unsupported_pair_family,
            :synthetic_unsupported_pair_block_path;
            supported_terms = (:overlap,),
        ),
    )
    return CPBMSkeleton.PairBlockMaterializationPlan(
        CPBMSkeleton.MetadataOnlyPairBlockMaterialization(),
        _skeleton_pair_operator_plan(),
        records,
        (; status = :available_pair_block_materialization_plan),
        (; fixture = :one_body_block_set_consumption_skeleton),
    )
end

@testset "CartesianPairBlockMaterialization one-body block-set consumption skeleton" begin
    plan = _skeleton_plan()
    inputs = (;
        parent_axis_counts = (2, 2, 2),
        overlap_1d = :overlap_factors,
        kinetic_1d = :kinetic_factors,
    )

    consumption = CPBMSkeleton._one_body_pair_block_set_consumption(
        plan;
        terms = (:overlap, :kinetic),
        inputs,
    )
    @test consumption.object_kind ==
          :cartesian_pair_block_mixed_one_body_block_set_consumption
    @test consumption.status ==
          :deferred_metadata_only_mixed_one_body_block_set_consumption
    @test isnothing(consumption.blocker)
    @test consumption.requested_terms == (:overlap, :kinetic)
    @test consumption.term_count == 2
    @test consumption.requested_materialize_terms == ()
    @test consumption.materialized_terms == ()
    @test consumption.deferred_terms == (:overlap, :kinetic)
    @test consumption.preflight_status ==
          :deferred_metadata_only_mixed_one_body_block_set_preflight
    @test isnothing(consumption.preflight_blocker)
    @test consumption.block_set_summary_status ==
          :deferred_metadata_only_mixed_one_body_block_set
    @test consumption.term_statuses == (
        (;
            term = :overlap,
            status = :deferred_metadata_only_mixed_one_body_pair_block_batch,
            batch_result_supplied = false,
        ),
        (;
            term = :kinetic,
            status = :deferred_metadata_only_mixed_one_body_pair_block_batch,
            batch_result_supplied = false,
        ),
    )
    @test consumption.total_materialized_count == 0
    @test consumption.total_skipped_count == 0
    @test !hasproperty(consumption, :term_batch_results)
    @test !hasproperty(consumption.block_set_summary, :term_batch_results)
    @test !hasproperty(consumption.block_set_summary.term_summaries[1], :batch_result)
    @test !hasproperty(consumption.block_set_summary.term_summaries[1], :block)
    @test consumption.result_terms_remain_separated
    @test !consumption.block_set_results_summed
    @test !consumption.term_batch_results_stored
    @test !consumption.one_term_consumer_called
    @test !consumption.factors_constructed
    @test !consumption.numerical_blocks_materialized
    @test !consumption.materialized
    @test !consumption.source_operator_blocks_materialized
    @test !consumption.final_pair_blocks_materialized
    @test !consumption.operator_blocks_materialized
    @test !consumption.hamiltonian_data_materialized
    @test !consumption.artifacts_materialized
    @test !consumption.global_operator_blocks_materialized
    @test !consumption.global_hamiltonian_data_materialized
    @test !consumption.global_artifacts_materialized
    @test !consumption.mixed_dispatcher_materialized
    @test !consumption.route_driver_wiring
    @test !consumption.coulomb_materialized
    @test !consumption.density_density_materialized
    @test !consumption.ida_mwg_data_materialized
    @test !consumption.pqs_lowdin_materialized
    @test !consumption.full_white_lindsey_route_assembled

    missing_input = CPBMSkeleton._one_body_pair_block_set_consumption(
        plan;
        terms = (:overlap, :kinetic),
        inputs = (; parent_axis_counts = (2, 2, 2), overlap_1d = :overlap_factors),
        materialize_terms = (:kinetic,),
    )
    @test missing_input.status == :blocked_mixed_one_body_block_set_consumption
    @test missing_input.blocker == :missing_required_one_body_factors
    @test missing_input.preflight_status ==
          :blocked_mixed_one_body_block_set_preflight
    @test missing_input.preflight_blocker == :missing_required_one_body_factors
    @test missing_input.total_materialized_count == 0
    @test missing_input.total_skipped_count == 0
    @test !missing_input.one_term_consumer_called
    @test !missing_input.numerical_blocks_materialized

    @test_throws ArgumentError CPBMSkeleton._one_body_pair_block_set_consumption(
        plan;
        terms = (:overlap, :coulomb),
        inputs,
    )
    @test_throws ArgumentError CPBMSkeleton._one_body_pair_block_set_consumption((;))
end

# Runtime role: tiny contract / semantic surface.
#
# Metadata-only block-set summary shape for future mixed one-body term-set
# consumption. This does not call the numerical one-term consumer or build
# factors. Use this when the block-set summary contract changes; use the smoke
# test for routine per-pass status/count validation.

using Test
using GaussletBases

const CPBMBlockSet = GaussletBases.CartesianPairBlockMaterialization
const CTLBlockSet = GaussletBases.CartesianTerminalLowering
const CRUBlockSet = GaussletBases.CartesianRetainedUnits
const CRTCBlockSet = GaussletBases.CartesianRetainedUnitTransformContracts
const CUPBlockSet = GaussletBases.CartesianUnitPairs
const CPOPBlockSet = GaussletBases.CartesianPairOperatorPlans

function _block_set_count(counts, field::Symbol, value)
    matches = Tuple(entry for entry in counts if getproperty(entry, field) == value)
    isempty(matches) && return 0
    return only(matches).count
end

function _block_set_pair_operator_plan()
    lowering_plan = CTLBlockSet.TerminalLoweringPlan(
        CTLBlockSet.PQSLowering(q = 2),
        (),
        (),
        (; status = :available_terminal_lowering_plan, materialized = false),
        (; fixture = :one_body_block_set_summary),
    )
    retained_plan = CRUBlockSet.RetainedUnitPlan(
        CRUBlockSet.MetadataOnlyRetainedUnits(),
        lowering_plan,
        (),
        (; status = :available_retained_unit_plan, materialized = false),
        (; fixture = :one_body_block_set_summary),
    )
    unit_pair_plan = CUPBlockSet.UnitPairPlan(
        CUPBlockSet.MetadataOnlyUnitPairs(),
        retained_plan,
        (),
        nothing,
        (; status = :available_unit_pair_plan, materialized = false),
        (; fixture = :one_body_block_set_summary),
    )
    transform_plan = CRTCBlockSet.RetainedUnitTransformContractPlan(
        CRTCBlockSet.MetadataOnlyRetainedUnitTransformContracts(),
        retained_plan,
        (),
        (; status = :available_transform_contract_plan, materialized = false),
        (; fixture = :one_body_block_set_summary),
    )
    return CPOPBlockSet.PairOperatorPlan(
        CPOPBlockSet.MetadataOnlyPairOperatorPlans(),
        unit_pair_plan,
        transform_plan,
        (),
        nothing,
        (; status = :available_pair_operator_plan, materialized = false),
        (; fixture = :one_body_block_set_summary),
    )
end

function _block_set_plan()
    return CPBMBlockSet.PairBlockMaterializationPlan(
        CPBMBlockSet.MetadataOnlyPairBlockMaterialization(),
        _block_set_pair_operator_plan(),
        (),
        (; status = :available_pair_block_materialization_plan),
        (; fixture = :one_body_block_set_summary),
    )
end

function _block_set_batch_result(term::Symbol)
    result = CPBMBlockSet.PairBlockMaterializationResult(
        term,
        (:left, :right),
        reshape([1.0], 1, 1),
        true,
        false,
        true,
        false,
        false,
        false,
        (; selector_family = :direct_direct),
    )
    skipped_records = (
        (;
            selector_family = :pqs_source_pair,
            blocker = :missing_one_body_factor_inputs,
        ),
    )
    return CPBMBlockSet.PairBlockMaterializationBatchResult(
        term,
        (result,),
        skipped_records,
        1,
        1,
        true,
        false,
        true,
        false,
        false,
        false,
        (;
            materialization_path = :mixed_one_body_pair_block_batch_selector,
            mixed_one_body_dispatcher = :direct_pqs_source_lw_plan_dispatcher,
            pair_block_record_count = 2,
            numerical_dispatch_scope =
                :direct_direct_pqs_source_pair_and_white_lindsey_boundary_stratum,
            factors_constructed = false,
            route_driver_wiring = false,
        ),
    )
end

@testset "CartesianPairBlockMaterialization one-body block-set summary shape" begin
    plan = _block_set_plan()
    skeleton = CPBMBlockSet._one_body_pair_block_set_summary(
        plan;
        terms = (:overlap, :kinetic),
    )

    @test skeleton.object_kind ==
          :cartesian_pair_block_mixed_one_body_block_set_summary
    @test skeleton.status == :deferred_metadata_only_mixed_one_body_block_set
    @test skeleton.plan_record_count == 0
    @test skeleton.plan_summary_status == :available_pair_block_materialization_plan
    @test skeleton.requested_terms == (:overlap, :kinetic)
    @test skeleton.term_count == 2
    @test skeleton.supplied_term_count == 0
    @test skeleton.deferred_term_count == 2
    @test skeleton.total_materialized_count == 0
    @test skeleton.total_skipped_count == 0
    @test length(skeleton.term_summaries) == 2
    @test !hasproperty(skeleton, :term_batch_results)
    @test !hasproperty(skeleton.term_summaries[1], :batch_result)
    @test !hasproperty(skeleton.term_summaries[1], :block)
    @test _block_set_count(
        skeleton.term_status_counts,
        :status,
        :deferred_metadata_only_mixed_one_body_pair_block_batch,
    ) == 2
    @test skeleton.term_statuses == (
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
    @test skeleton.result_terms_remain_separated
    @test !skeleton.block_set_results_summed
    @test !skeleton.term_batch_results_stored_in_summary
    @test skeleton.factor_provider_scope == :caller_supplied_or_family_provider
    @test !skeleton.factors_constructed
    @test !skeleton.materialized
    @test !skeleton.source_operator_blocks_materialized
    @test !skeleton.final_pair_blocks_materialized
    @test !skeleton.operator_blocks_materialized
    @test !skeleton.hamiltonian_data_materialized
    @test !skeleton.artifacts_materialized
    @test !skeleton.global_operator_blocks_materialized
    @test !skeleton.global_hamiltonian_data_materialized
    @test !skeleton.global_artifacts_materialized
    @test !skeleton.mixed_dispatcher_materialized
    @test !skeleton.route_driver_wiring
    @test !skeleton.coulomb_materialized
    @test !skeleton.density_density_materialized
    @test !skeleton.ida_mwg_data_materialized
    @test !skeleton.pqs_lowdin_materialized
    @test !skeleton.full_white_lindsey_route_assembled

    supplied = CPBMBlockSet._one_body_pair_block_set_summary(
        plan;
        terms = (:overlap, :kinetic),
        term_batch_results = (; overlap = _block_set_batch_result(:overlap)),
    )
    @test supplied.status == :partially_deferred_mixed_one_body_block_set
    @test supplied.supplied_term_count == 1
    @test supplied.deferred_term_count == 1
    @test supplied.total_materialized_count == 1
    @test supplied.total_skipped_count == 1
    @test supplied.materialized
    @test !supplied.source_operator_blocks_materialized
    @test supplied.final_pair_blocks_materialized
    @test !supplied.operator_blocks_materialized
    @test !supplied.hamiltonian_data_materialized
    @test !supplied.artifacts_materialized
    @test supplied.term_summaries[1].batch_result_supplied
    @test !supplied.term_summaries[2].batch_result_supplied
    @test _block_set_count(
        supplied.materialized_selector_family_counts,
        :selector_family,
        :direct_direct,
    ) == 1
    @test _block_set_count(
        supplied.skipped_selector_family_counts,
        :selector_family,
        :pqs_source_pair,
    ) == 1
    @test _block_set_count(
        supplied.skipped_blocker_counts,
        :blocker,
        :missing_one_body_factor_inputs,
    ) == 1

    tuple_supplied = CPBMBlockSet._one_body_pair_block_set_summary(
        plan;
        terms = (:overlap,),
        term_batch_results = (_block_set_batch_result(:overlap),),
    )
    @test tuple_supplied.status == :partially_materialized_mixed_one_body_block_set
    @test tuple_supplied.total_materialized_count == 1
    @test tuple_supplied.total_skipped_count == 1

    @test_throws ArgumentError CPBMBlockSet._one_body_pair_block_set_summary(
        plan;
        terms = (:overlap, :coulomb),
    )
    @test_throws ArgumentError CPBMBlockSet._one_body_pair_block_set_summary(
        plan;
        terms = (:position_x,),
        term_batch_results = (; position_x = _block_set_batch_result(:overlap)),
    )
    @test_throws ArgumentError CPBMBlockSet._one_body_pair_block_set_summary(
        plan;
        terms = (:overlap,),
        term_batch_results = (; overlap = :not_a_batch_result),
    )
    @test_throws ArgumentError CPBMBlockSet._one_body_pair_block_set_summary((;))
end

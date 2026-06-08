# Runtime role: tiny contract.
#
# Synthetic compact batch-summary contract. This is cheap, but the consumer
# smoke is the preferred single routine per-pass test for mixed-consumer edits.

using Test
using GaussletBases

const CPBMBatchSummary = GaussletBases.CartesianPairBlockMaterialization

function _batch_summary_count(counts, field::Symbol, value)
    matches = Tuple(entry for entry in counts if getproperty(entry, field) == value)
    isempty(matches) && return 0
    return only(matches).count
end

function _batch_summary_result(
    selector_family::Symbol;
    source_operator_blocks_materialized::Bool,
    final_pair_blocks_materialized::Bool,
)
    return CPBMBatchSummary.PairBlockMaterializationResult(
        :overlap,
        (:left, :right),
        reshape([1.0], 1, 1),
        true,
        source_operator_blocks_materialized,
        final_pair_blocks_materialized,
        false,
        false,
        false,
        (; selector_family),
    )
end

function _batch_summary_mixed_batch()
    direct_result = _batch_summary_result(
        :direct_direct;
        source_operator_blocks_materialized = false,
        final_pair_blocks_materialized = true,
    )
    pqs_result = _batch_summary_result(
        :pqs_source_pair;
        source_operator_blocks_materialized = true,
        final_pair_blocks_materialized = false,
    )
    lw_result = _batch_summary_result(
        :white_lindsey_boundary_stratum;
        source_operator_blocks_materialized = false,
        final_pair_blocks_materialized = true,
    )
    skipped_records = (
        (;
            selector_family = :pqs_source_pair,
            blocker = :missing_one_body_factor_inputs,
        ),
        (;
            selector_family = :white_lindsey_boundary_stratum,
            blocker = :missing_white_lindsey_unit_pair,
        ),
    )
    return CPBMBatchSummary.PairBlockMaterializationBatchResult(
        :overlap,
        (direct_result, pqs_result, lw_result),
        skipped_records,
        3,
        2,
        true,
        true,
        true,
        false,
        false,
        false,
        (;
            materialization_path = :mixed_one_body_pair_block_batch_selector,
            mixed_one_body_dispatcher = :direct_pqs_source_lw_plan_dispatcher,
            pair_block_record_count = 5,
            numerical_dispatch_scope =
                :direct_direct_pqs_source_pair_and_white_lindsey_boundary_stratum,
            factors_constructed = false,
            route_driver_wiring = false,
        ),
    )
end

function _batch_summary_all_skipped_batch()
    return CPBMBatchSummary.PairBlockMaterializationBatchResult(
        :position_y,
        (),
        (
            (;
                selector_family = :direct_direct,
                blocker = :missing_one_body_factor_inputs,
            ),
        ),
        0,
        1,
        false,
        false,
        false,
        false,
        false,
        false,
        (;
            materialization_path = :mixed_one_body_pair_block_batch_selector,
            mixed_one_body_dispatcher = :direct_pqs_source_lw_plan_dispatcher,
            pair_block_record_count = 1,
            numerical_dispatch_scope =
                :direct_direct_pqs_source_pair_and_white_lindsey_boundary_stratum,
            factors_constructed = false,
            route_driver_wiring = false,
        ),
    )
end

@testset "CartesianPairBlockMaterialization mixed one-body batch summary" begin
    mixed_summary =
        CPBMBatchSummary._one_body_pair_block_batch_summary(_batch_summary_mixed_batch())

    @test mixed_summary.object_kind ==
          :cartesian_pair_block_mixed_one_body_batch_summary
    @test mixed_summary.status ==
          :partially_materialized_mixed_one_body_pair_block_batch
    @test mixed_summary.term == :overlap
    @test mixed_summary.materialized_count == 3
    @test mixed_summary.skipped_count == 2
    @test mixed_summary.pair_block_record_count == 5
    @test mixed_summary.materialization_path ==
          :mixed_one_body_pair_block_batch_selector
    @test mixed_summary.mixed_one_body_dispatcher ==
          :direct_pqs_source_lw_plan_dispatcher
    @test mixed_summary.numerical_dispatch_scope ==
          :direct_direct_pqs_source_pair_and_white_lindsey_boundary_stratum

    @test _batch_summary_count(
        mixed_summary.materialized_selector_family_counts,
        :selector_family,
        :direct_direct,
    ) == 1
    @test _batch_summary_count(
        mixed_summary.materialized_selector_family_counts,
        :selector_family,
        :pqs_source_pair,
    ) == 1
    @test _batch_summary_count(
        mixed_summary.materialized_selector_family_counts,
        :selector_family,
        :white_lindsey_boundary_stratum,
    ) == 1
    @test _batch_summary_count(
        mixed_summary.skipped_selector_family_counts,
        :selector_family,
        :pqs_source_pair,
    ) == 1
    @test _batch_summary_count(
        mixed_summary.skipped_selector_family_counts,
        :selector_family,
        :white_lindsey_boundary_stratum,
    ) == 1
    @test _batch_summary_count(
        mixed_summary.skipped_blocker_counts,
        :blocker,
        :missing_one_body_factor_inputs,
    ) == 1
    @test _batch_summary_count(
        mixed_summary.skipped_blocker_counts,
        :blocker,
        :missing_white_lindsey_unit_pair,
    ) == 1

    @test mixed_summary.source_space_only_result_count == 1
    @test mixed_summary.final_local_block_result_count == 2
    @test mixed_summary.direct_direct_materialized
    @test mixed_summary.pqs_source_pair_materialized
    @test mixed_summary.white_lindsey_boundary_stratum_materialized
    @test mixed_summary.white_lindsey_materialized
    @test mixed_summary.source_operator_blocks_materialized
    @test mixed_summary.final_pair_blocks_materialized
    @test !mixed_summary.operator_blocks_materialized
    @test !mixed_summary.hamiltonian_data_materialized
    @test !mixed_summary.artifacts_materialized
    @test !mixed_summary.global_operator_blocks_materialized
    @test !mixed_summary.global_hamiltonian_data_materialized
    @test !mixed_summary.global_artifacts_materialized
    @test !mixed_summary.factors_constructed
    @test !mixed_summary.route_driver_wiring

    skipped_summary = CPBMBatchSummary._one_body_pair_block_batch_summary(
        _batch_summary_all_skipped_batch(),
    )
    @test skipped_summary.status == :skipped_mixed_one_body_pair_block_batch
    @test skipped_summary.materialized_count == 0
    @test skipped_summary.skipped_count == 1
    @test skipped_summary.source_space_only_result_count == 0
    @test skipped_summary.final_local_block_result_count == 0
    @test !skipped_summary.direct_direct_materialized
    @test !skipped_summary.pqs_source_pair_materialized
    @test !skipped_summary.white_lindsey_materialized
    @test _batch_summary_count(
        skipped_summary.skipped_selector_family_counts,
        :selector_family,
        :direct_direct,
    ) == 1

    @test_throws ArgumentError CPBMBatchSummary._one_body_pair_block_batch_summary((;))
end

# Runtime role: tiny contract.
#
# Opt-in safe-term materialization for the mixed one-body block-set consumer.
# This delegates only explicit safe one-body terms to the existing one-term
# consumer and leaves other requested terms deferred.

using Test
using GaussletBases

const CPBMBlockOverlap = GaussletBases.CartesianPairBlockMaterialization
const CPBBlockOverlap = GaussletBases.CartesianCPB
const CTLBlockOverlap = GaussletBases.CartesianTerminalLowering
const CRUBlockOverlap = GaussletBases.CartesianRetainedUnits
const CRTCBlockOverlap = GaussletBases.CartesianRetainedUnitTransformContracts
const CUPBlockOverlap = GaussletBases.CartesianUnitPairs
const CPOPBlockOverlap = GaussletBases.CartesianPairOperatorPlans

function _block_overlap_count(counts, field::Symbol, value)
    matches = Tuple(entry for entry in counts if getproperty(entry, field) == value)
    isempty(matches) && return 0
    return only(matches).count
end

function _block_overlap_summary_has_no_matrix_fields(block_set_summary)
    forbidden_fields = (
        :term_batch_results,
        :batch_result,
        :block,
        :blocks,
        :materialized_results,
        :skipped_records,
        :preflight_summary,
        :block_set_summary,
        :term_summaries,
    )
    return !any(field -> hasproperty(block_set_summary, field), forbidden_fields)
end

function _block_overlap_pair_operator_plan()
    lowering_plan = CTLBlockOverlap.TerminalLoweringPlan(
        CTLBlockOverlap.PQSLowering(q = 2),
        (),
        (),
        (; status = :available_terminal_lowering_plan, materialized = false),
        (; fixture = :one_body_block_set_overlap_consumption),
    )
    retained_plan = CRUBlockOverlap.RetainedUnitPlan(
        CRUBlockOverlap.MetadataOnlyRetainedUnits(),
        lowering_plan,
        (),
        (; status = :available_retained_unit_plan, materialized = false),
        (; fixture = :one_body_block_set_overlap_consumption),
    )
    unit_pair_plan = CUPBlockOverlap.UnitPairPlan(
        CUPBlockOverlap.MetadataOnlyUnitPairs(),
        retained_plan,
        (),
        nothing,
        (; status = :available_unit_pair_plan, materialized = false),
        (; fixture = :one_body_block_set_overlap_consumption),
    )
    transform_plan = CRTCBlockOverlap.RetainedUnitTransformContractPlan(
        CRTCBlockOverlap.MetadataOnlyRetainedUnitTransformContracts(),
        retained_plan,
        (),
        (; status = :available_transform_contract_plan, materialized = false),
        (; fixture = :one_body_block_set_overlap_consumption),
    )
    return CPOPBlockOverlap.PairOperatorPlan(
        CPOPBlockOverlap.MetadataOnlyPairOperatorPlans(),
        unit_pair_plan,
        transform_plan,
        (),
        nothing,
        (; status = :available_pair_operator_plan, materialized = false),
        (; fixture = :one_body_block_set_overlap_consumption),
    )
end

function _block_overlap_record(
    pair_key::Tuple{Symbol,Symbol},
    pair_index::Int,
    pair_family::Symbol,
    materialization_path::Symbol;
    metadata = (;),
)
    return CPBMBlockOverlap.PairBlockMaterializationRecord(
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
        :ready_metadata_only_not_materialized,
        nothing,
        false,
        NamedTuple(metadata),
    )
end

function _block_overlap_plan()
    direct_left = CPBBlockOverlap.cpb(
        1:1,
        1:2,
        1:1;
        role = :block_overlap_direct_left_source,
    )
    direct_right = CPBBlockOverlap.cpb(
        2:2,
        1:2,
        1:1;
        role = :block_overlap_direct_right_source,
    )
    records = (
        _block_overlap_record(
            (:direct_left, :direct_right),
            1,
            :direct_direct,
            :direct_direct_pair_block_materialization_pilot;
            metadata = (;
                left_source_cpbs = (direct_left,),
                right_source_cpbs = (direct_right,),
            ),
        ),
        _block_overlap_record(
            (:pqs_left, :pqs_right),
            2,
            :pqs_pqs,
            :pqs_source_pair_preflight;
            metadata = (;
                left_source_mode_dims = (2, 2, 2),
                right_source_mode_dims = (2, 2, 2),
                left_source_mode_count = 8,
                right_source_mode_count = 8,
                source_mode_ordering = :x_major_y_major_z_fast,
                left_source_mode_ordering = :x_major_y_major_z_fast,
                right_source_mode_ordering = :x_major_y_major_z_fast,
            ),
        ),
        _block_overlap_record(
            (:unsupported_left, :unsupported_right),
            3,
            :unsupported_pair_family,
            :synthetic_unsupported_pair_block_path,
        ),
    )
    return CPBMBlockOverlap.PairBlockMaterializationPlan(
        CPBMBlockOverlap.MetadataOnlyPairBlockMaterialization(),
        _block_overlap_pair_operator_plan(),
        records,
        (; status = :available_pair_block_materialization_plan),
        (; fixture = :one_body_block_set_overlap_consumption),
    )
end

function _block_overlap_1d()
    return (;
        x = [1.0 0.2; 0.2 1.1],
        y = [1.2 0.3; 0.3 1.3],
        z = [1.4 0.4; 0.4 1.5],
    )
end

@testset "CartesianPairBlockMaterialization one-body block-set overlap consumption" begin
    plan = _block_overlap_plan()
    inputs = (;
        parent_axis_counts = (2, 2, 2),
        overlap_1d = _block_overlap_1d(),
        position_1d = _block_overlap_1d(),
        x2_1d = _block_overlap_1d(),
        kinetic_1d = _block_overlap_1d(),
    )
    safe_terms = (
        :overlap,
        :position_x,
        :position_y,
        :position_z,
        :x2_x,
        :x2_y,
        :x2_z,
        :kinetic,
    )

    deferred = CPBMBlockOverlap._one_body_pair_block_set_consumption(
        plan;
        terms = (:overlap, :kinetic),
        inputs,
    )
    @test deferred.status ==
          :deferred_metadata_only_mixed_one_body_block_set_consumption
    @test deferred.materialized_terms == ()
    @test deferred.deferred_terms == (:overlap, :kinetic)
    @test !hasproperty(deferred, :term_batch_results)
    @test !deferred.one_term_consumer_called
    @test deferred.total_materialized_count == 0
    @test deferred.total_skipped_count == 0

    no_factor_deferred = CPBMBlockOverlap._one_body_pair_block_set_consumption(
        plan;
        terms = safe_terms,
        inputs = (;),
    )
    @test no_factor_deferred.status ==
          :deferred_metadata_only_mixed_one_body_block_set_consumption
    @test no_factor_deferred.preflight_status ==
          :deferred_metadata_only_mixed_one_body_block_set_preflight
    @test no_factor_deferred.preflight_summary.term_set_input_status ==
          :not_required_deferred_terms_only
    @test no_factor_deferred.materialized_terms == ()
    @test no_factor_deferred.deferred_terms == safe_terms
    @test no_factor_deferred.total_materialized_count == 0
    @test no_factor_deferred.total_skipped_count == 0
    @test !no_factor_deferred.one_term_consumer_called
    @test !hasproperty(no_factor_deferred, :term_batch_results)

    no_factor_deferred_summary =
        CPBMBlockOverlap._one_body_pair_block_set_consumption_summary(
            no_factor_deferred,
        )
    @test no_factor_deferred_summary.object_kind ==
          :cartesian_pair_block_mixed_one_body_block_set_consumption_summary
    @test no_factor_deferred_summary.status ==
          :deferred_metadata_only_mixed_one_body_block_set_consumption
    @test isnothing(no_factor_deferred_summary.blocker)
    @test no_factor_deferred_summary.requested_terms == safe_terms
    @test no_factor_deferred_summary.requested_materialize_terms == ()
    @test no_factor_deferred_summary.materialized_terms == ()
    @test no_factor_deferred_summary.deferred_terms == safe_terms
    @test no_factor_deferred_summary.preflight_status ==
          :deferred_metadata_only_mixed_one_body_block_set_preflight
    @test no_factor_deferred_summary.block_set_summary_status ==
          :deferred_metadata_only_mixed_one_body_block_set
    @test no_factor_deferred_summary.supplied_term_count == 0
    @test no_factor_deferred_summary.deferred_term_count == length(safe_terms)
    @test no_factor_deferred_summary.total_materialized_count == 0
    @test no_factor_deferred_summary.total_skipped_count == 0
    @test !no_factor_deferred_summary.term_batch_results_available_on_consumption
    @test !no_factor_deferred_summary.term_batch_results_stored_in_summary
    @test !no_factor_deferred_summary.materialized
    @test _block_overlap_summary_has_no_matrix_fields(no_factor_deferred_summary)

    overlap = CPBMBlockOverlap._one_body_pair_block_set_consumption(
        plan;
        terms = (:overlap, :kinetic),
        inputs,
        materialize_terms = (:overlap,),
    )
    @test overlap.status ==
          :partially_materialized_mixed_one_body_block_set_consumption
    @test overlap.materialized_terms == (:overlap,)
    @test overlap.deferred_terms == (:kinetic,)
    @test overlap.requested_materialize_terms == (:overlap,)
    @test hasproperty(overlap, :term_batch_results)
    @test overlap.term_batch_results.overlap isa
          CPBMBlockOverlap.PairBlockMaterializationBatchResult
    @test overlap.one_term_consumer_called
    @test overlap.term_batch_results_stored
    @test overlap.total_materialized_count == 2
    @test overlap.total_skipped_count == 1
    @test overlap.materialized
    @test overlap.numerical_blocks_materialized
    @test overlap.source_operator_blocks_materialized
    @test overlap.final_pair_blocks_materialized
    @test !overlap.operator_blocks_materialized
    @test !overlap.hamiltonian_data_materialized
    @test !overlap.artifacts_materialized
    @test !overlap.route_driver_wiring
    @test !overlap.coulomb_materialized
    @test !overlap.ida_mwg_data_materialized
    @test !overlap.pqs_lowdin_materialized
    @test !overlap.full_white_lindsey_route_assembled

    @test overlap.block_set_summary.supplied_term_count == 1
    @test overlap.block_set_summary.deferred_term_count == 1
    @test !hasproperty(overlap.block_set_summary, :term_batch_results)
    @test !hasproperty(overlap.block_set_summary.term_summaries[1], :batch_result)
    @test !hasproperty(overlap.block_set_summary.term_summaries[1], :block)
    @test overlap.term_statuses[1].term == :overlap
    @test overlap.term_statuses[1].status ==
          :partially_materialized_mixed_one_body_pair_block_batch
    @test overlap.term_statuses[1].batch_result_supplied
    @test overlap.term_statuses[2].term == :kinetic
    @test overlap.term_statuses[2].status ==
          :deferred_metadata_only_mixed_one_body_pair_block_batch
    @test !overlap.term_statuses[2].batch_result_supplied
    @test _block_overlap_count(
        overlap.block_set_summary.materialized_selector_family_counts,
        :selector_family,
        :direct_direct,
    ) == 1
    @test _block_overlap_count(
        overlap.block_set_summary.materialized_selector_family_counts,
        :selector_family,
        :pqs_source_pair,
    ) == 1
    @test _block_overlap_count(
        overlap.block_set_summary.skipped_selector_family_counts,
        :selector_family,
        :unsupported,
    ) == 1
    @test _block_overlap_count(
        overlap.block_set_summary.skipped_blocker_counts,
        :blocker,
        :unsupported_pair_block_materialization_path,
    ) == 1

    overlap_from_all_requested = CPBMBlockOverlap._one_body_pair_block_set_consumption(
        plan;
        terms = safe_terms,
        inputs = (;
            parent_axis_counts = (2, 2, 2),
            overlap_1d = _block_overlap_1d(),
        ),
        materialize_terms = (:overlap,),
    )
    @test overlap_from_all_requested.status ==
          :partially_materialized_mixed_one_body_block_set_consumption
    @test overlap_from_all_requested.materialized_terms == (:overlap,)
    @test overlap_from_all_requested.deferred_terms == (
        :position_x,
        :position_y,
        :position_z,
        :x2_x,
        :x2_y,
        :x2_z,
        :kinetic,
    )
    @test overlap_from_all_requested.total_materialized_count == 2
    @test overlap_from_all_requested.total_skipped_count == 1
    @test overlap_from_all_requested.block_set_summary.supplied_term_count == 1
    @test overlap_from_all_requested.block_set_summary.deferred_term_count == 7
    @test overlap_from_all_requested.preflight_summary.required_factor_names ==
          (:overlap_1d,)
    @test overlap_from_all_requested.preflight_summary.missing_factor_names == ()

    overlap_from_all_summary =
        CPBMBlockOverlap._one_body_pair_block_set_consumption_summary(
            overlap_from_all_requested,
        )
    @test overlap_from_all_summary.status ==
          :partially_materialized_mixed_one_body_block_set_consumption
    @test overlap_from_all_summary.requested_terms == safe_terms
    @test overlap_from_all_summary.requested_materialize_terms == (:overlap,)
    @test overlap_from_all_summary.materialized_terms == (:overlap,)
    @test overlap_from_all_summary.deferred_terms == (
        :position_x,
        :position_y,
        :position_z,
        :x2_x,
        :x2_y,
        :x2_z,
        :kinetic,
    )
    @test overlap_from_all_summary.preflight_status ==
          :partially_ready_mixed_one_body_block_set_preflight
    @test overlap_from_all_summary.block_set_summary_status ==
          :partially_deferred_mixed_one_body_block_set
    @test overlap_from_all_summary.supplied_term_count == 1
    @test overlap_from_all_summary.deferred_term_count == 7
    @test overlap_from_all_summary.total_materialized_count == 2
    @test overlap_from_all_summary.total_skipped_count == 1
    @test overlap_from_all_summary.term_batch_results_available_on_consumption
    @test !overlap_from_all_summary.term_batch_results_stored_in_summary
    @test overlap_from_all_summary.source_operator_blocks_materialized
    @test overlap_from_all_summary.final_pair_blocks_materialized
    @test !overlap_from_all_summary.operator_blocks_materialized
    @test _block_overlap_count(
        overlap_from_all_summary.materialized_selector_family_counts,
        :selector_family,
        :direct_direct,
    ) == 1
    @test _block_overlap_count(
        overlap_from_all_summary.materialized_selector_family_counts,
        :selector_family,
        :pqs_source_pair,
    ) == 1
    @test _block_overlap_count(
        overlap_from_all_summary.skipped_blocker_counts,
        :blocker,
        :unsupported_pair_block_materialization_path,
    ) == 1
    @test _block_overlap_summary_has_no_matrix_fields(overlap_from_all_summary)

    position = CPBMBlockOverlap._one_body_pair_block_set_consumption(
        plan;
        terms = (:overlap, :position_y, :kinetic),
        inputs,
        materialize_terms = (:position_y,),
    )
    @test position.status ==
          :partially_materialized_mixed_one_body_block_set_consumption
    @test position.materialized_terms == (:position_y,)
    @test position.deferred_terms == (:overlap, :kinetic)
    @test position.total_materialized_count == 2
    @test position.total_skipped_count == 1
    @test position.source_operator_blocks_materialized
    @test position.final_pair_blocks_materialized
    @test position.term_statuses[1].status ==
          :deferred_metadata_only_mixed_one_body_pair_block_batch
    @test position.term_statuses[2].term == :position_y
    @test position.term_statuses[2].status ==
          :partially_materialized_mixed_one_body_pair_block_batch
    @test position.term_statuses[2].batch_result_supplied
    @test position.term_statuses[3].status ==
          :deferred_metadata_only_mixed_one_body_pair_block_batch
    @test _block_overlap_count(
        position.block_set_summary.materialized_selector_family_counts,
        :selector_family,
        :direct_direct,
    ) == 1
    @test _block_overlap_count(
        position.block_set_summary.materialized_selector_family_counts,
        :selector_family,
        :pqs_source_pair,
    ) == 1

    overlap_position = CPBMBlockOverlap._one_body_pair_block_set_consumption(
        plan;
        terms = (:overlap, :position_x),
        inputs,
        materialize_terms = (:overlap, :position_x),
    )
    @test overlap_position.materialized_terms == (:overlap, :position_x)
    @test overlap_position.deferred_terms == ()
    @test overlap_position.total_materialized_count == 4
    @test overlap_position.total_skipped_count == 2
    @test overlap_position.block_set_summary.supplied_term_count == 2
    @test overlap_position.block_set_summary.deferred_term_count == 0
    @test _block_overlap_count(
        overlap_position.block_set_summary.materialized_selector_family_counts,
        :selector_family,
        :direct_direct,
    ) == 2
    @test _block_overlap_count(
        overlap_position.block_set_summary.materialized_selector_family_counts,
        :selector_family,
        :pqs_source_pair,
    ) == 2

    all_positions = CPBMBlockOverlap._one_body_pair_block_set_consumption(
        plan;
        terms = (:position_x, :position_y, :position_z),
        inputs,
        materialize_terms = (:position_x, :position_y, :position_z),
    )
    @test all_positions.materialized_terms ==
          (:position_x, :position_y, :position_z)
    @test all_positions.deferred_terms == ()
    @test all_positions.total_materialized_count == 6
    @test all_positions.total_skipped_count == 3
    @test all_positions.block_set_summary.supplied_term_count == 3
    @test all_positions.block_set_summary.deferred_term_count == 0
    @test _block_overlap_count(
        all_positions.block_set_summary.materialized_selector_family_counts,
        :selector_family,
        :direct_direct,
    ) == 3
    @test _block_overlap_count(
        all_positions.block_set_summary.materialized_selector_family_counts,
        :selector_family,
        :pqs_source_pair,
    ) == 3
    @test _block_overlap_count(
        all_positions.block_set_summary.skipped_selector_family_counts,
        :selector_family,
        :unsupported,
    ) == 3

    x2 = CPBMBlockOverlap._one_body_pair_block_set_consumption(
        plan;
        terms = (:overlap, :x2_y, :kinetic),
        inputs,
        materialize_terms = (:x2_y,),
    )
    @test x2.status ==
          :partially_materialized_mixed_one_body_block_set_consumption
    @test x2.materialized_terms == (:x2_y,)
    @test x2.deferred_terms == (:overlap, :kinetic)
    @test x2.total_materialized_count == 2
    @test x2.total_skipped_count == 1
    @test x2.source_operator_blocks_materialized
    @test x2.final_pair_blocks_materialized
    @test x2.term_statuses[2].term == :x2_y
    @test x2.term_statuses[2].status ==
          :partially_materialized_mixed_one_body_pair_block_batch
    @test _block_overlap_count(
        x2.block_set_summary.materialized_selector_family_counts,
        :selector_family,
        :direct_direct,
    ) == 1
    @test _block_overlap_count(
        x2.block_set_summary.materialized_selector_family_counts,
        :selector_family,
        :pqs_source_pair,
    ) == 1

    all_x2 = CPBMBlockOverlap._one_body_pair_block_set_consumption(
        plan;
        terms = (:x2_x, :x2_y, :x2_z),
        inputs,
        materialize_terms = (:x2_x, :x2_y, :x2_z),
    )
    @test all_x2.materialized_terms == (:x2_x, :x2_y, :x2_z)
    @test all_x2.deferred_terms == ()
    @test all_x2.total_materialized_count == 6
    @test all_x2.total_skipped_count == 3
    @test all_x2.block_set_summary.supplied_term_count == 3
    @test all_x2.block_set_summary.deferred_term_count == 0
    @test _block_overlap_count(
        all_x2.block_set_summary.materialized_selector_family_counts,
        :selector_family,
        :direct_direct,
    ) == 3
    @test _block_overlap_count(
        all_x2.block_set_summary.materialized_selector_family_counts,
        :selector_family,
        :pqs_source_pair,
    ) == 3
    @test _block_overlap_count(
        all_x2.block_set_summary.skipped_selector_family_counts,
        :selector_family,
        :unsupported,
    ) == 3

    overlap_position_x2 = CPBMBlockOverlap._one_body_pair_block_set_consumption(
        plan;
        terms = (:overlap, :position_x, :x2_z),
        inputs,
        materialize_terms = (:overlap, :position_x, :x2_z),
    )
    @test overlap_position_x2.materialized_terms ==
          (:overlap, :position_x, :x2_z)
    @test overlap_position_x2.deferred_terms == ()
    @test overlap_position_x2.total_materialized_count == 6
    @test overlap_position_x2.total_skipped_count == 3
    @test overlap_position_x2.block_set_summary.supplied_term_count == 3
    @test overlap_position_x2.block_set_summary.deferred_term_count == 0
    @test _block_overlap_count(
        overlap_position_x2.block_set_summary.materialized_selector_family_counts,
        :selector_family,
        :direct_direct,
    ) == 3
    @test _block_overlap_count(
        overlap_position_x2.block_set_summary.materialized_selector_family_counts,
        :selector_family,
        :pqs_source_pair,
    ) == 3

    kinetic = CPBMBlockOverlap._one_body_pair_block_set_consumption(
        plan;
        terms = (:overlap, :kinetic),
        inputs,
        materialize_terms = (:kinetic,),
    )
    @test kinetic.status ==
          :partially_materialized_mixed_one_body_block_set_consumption
    @test kinetic.materialized_terms == (:kinetic,)
    @test kinetic.deferred_terms == (:overlap,)
    @test kinetic.total_materialized_count == 2
    @test kinetic.total_skipped_count == 1
    @test kinetic.source_operator_blocks_materialized
    @test kinetic.final_pair_blocks_materialized
    @test kinetic.term_statuses[2].term == :kinetic
    @test kinetic.term_statuses[2].status ==
          :partially_materialized_mixed_one_body_pair_block_batch

    all_safe = CPBMBlockOverlap._one_body_pair_block_set_consumption(
        plan;
        terms = safe_terms,
        inputs,
        materialize_terms = safe_terms,
    )
    @test all_safe.materialized_terms == safe_terms
    @test all_safe.deferred_terms == ()
    @test all_safe.total_materialized_count == 16
    @test all_safe.total_skipped_count == 8
    @test all_safe.block_set_summary.supplied_term_count == 8
    @test all_safe.block_set_summary.deferred_term_count == 0
    @test all_safe.result_terms_remain_separated
    @test !all_safe.block_set_results_summed
    @test !all_safe.operator_blocks_materialized
    @test !all_safe.hamiltonian_data_materialized
    @test !all_safe.artifacts_materialized
    @test !all_safe.route_driver_wiring
    @test !all_safe.coulomb_materialized
    @test !all_safe.ida_mwg_data_materialized
    @test !all_safe.pqs_lowdin_materialized
    @test !all_safe.full_white_lindsey_route_assembled
    @test !hasproperty(all_safe.block_set_summary, :term_batch_results)
    @test !hasproperty(all_safe.block_set_summary.term_summaries[1], :batch_result)
    @test !hasproperty(all_safe.block_set_summary.term_summaries[1], :block)
    @test _block_overlap_count(
        all_safe.block_set_summary.materialized_selector_family_counts,
        :selector_family,
        :direct_direct,
    ) == 8
    @test _block_overlap_count(
        all_safe.block_set_summary.materialized_selector_family_counts,
        :selector_family,
        :pqs_source_pair,
    ) == 8
    @test _block_overlap_count(
        all_safe.block_set_summary.skipped_selector_family_counts,
        :selector_family,
        :unsupported,
    ) == 8

    all_safe_summary =
        CPBMBlockOverlap._one_body_pair_block_set_consumption_summary(all_safe)
    @test all_safe_summary.status ==
          :partially_materialized_mixed_one_body_block_set_consumption
    @test all_safe_summary.requested_terms == safe_terms
    @test all_safe_summary.requested_materialize_terms == safe_terms
    @test all_safe_summary.materialized_terms == safe_terms
    @test all_safe_summary.deferred_terms == ()
    @test all_safe_summary.block_set_summary_status ==
          :partially_materialized_mixed_one_body_block_set
    @test all_safe_summary.supplied_term_count == 8
    @test all_safe_summary.deferred_term_count == 0
    @test all_safe_summary.total_materialized_count == 16
    @test all_safe_summary.total_skipped_count == 8
    @test all_safe_summary.term_batch_results_available_on_consumption
    @test !all_safe_summary.term_batch_results_stored_in_summary
    @test all_safe_summary.materialized
    @test all_safe_summary.source_operator_blocks_materialized
    @test all_safe_summary.final_pair_blocks_materialized
    @test !all_safe_summary.operator_blocks_materialized
    @test !all_safe_summary.hamiltonian_data_materialized
    @test !all_safe_summary.artifacts_materialized
    @test !all_safe_summary.route_driver_wiring
    @test !all_safe_summary.coulomb_materialized
    @test !all_safe_summary.ida_mwg_data_materialized
    @test !all_safe_summary.pqs_lowdin_materialized
    @test !all_safe_summary.full_white_lindsey_route_assembled
    @test _block_overlap_count(
        all_safe_summary.materialized_selector_family_counts,
        :selector_family,
        :direct_direct,
    ) == 8
    @test _block_overlap_count(
        all_safe_summary.materialized_selector_family_counts,
        :selector_family,
        :pqs_source_pair,
    ) == 8
    @test _block_overlap_count(
        all_safe_summary.skipped_selector_family_counts,
        :selector_family,
        :unsupported,
    ) == 8
    @test _block_overlap_count(
        all_safe_summary.skipped_blocker_counts,
        :blocker,
        :unsupported_pair_block_materialization_path,
    ) == 8
    @test _block_overlap_summary_has_no_matrix_fields(all_safe_summary)
    @test all(
        term_status -> !hasproperty(term_status, :batch_result) &&
                       !hasproperty(term_status, :block),
        all_safe_summary.term_statuses,
    )

    missing_position = CPBMBlockOverlap._one_body_pair_block_set_consumption(
        plan;
        terms = safe_terms,
        inputs = (;
            parent_axis_counts = (2, 2, 2),
            overlap_1d = _block_overlap_1d(),
        ),
        materialize_terms = (:position_x,),
    )
    @test missing_position.status ==
          :blocked_mixed_one_body_block_set_consumption
    @test missing_position.blocker == :missing_required_one_body_factors
    @test missing_position.preflight_status ==
          :blocked_mixed_one_body_block_set_preflight
    @test missing_position.preflight_summary.required_factor_names ==
          (:overlap_1d, :position_1d)
    @test missing_position.preflight_summary.missing_factor_names ==
          (:position_1d,)
    @test missing_position.total_materialized_count == 0
    @test missing_position.total_skipped_count == 0
    @test !missing_position.one_term_consumer_called
    @test !hasproperty(missing_position, :term_batch_results)

    missing_position_summary =
        CPBMBlockOverlap._one_body_pair_block_set_consumption_summary(
            missing_position,
        )
    @test missing_position_summary.status ==
          :blocked_mixed_one_body_block_set_consumption
    @test missing_position_summary.blocker == :missing_required_one_body_factors
    @test missing_position_summary.preflight_status ==
          :blocked_mixed_one_body_block_set_preflight
    @test missing_position_summary.preflight_blocker ==
          :missing_required_one_body_factors
    @test missing_position_summary.block_set_summary_status ==
          :deferred_metadata_only_mixed_one_body_block_set
    @test missing_position_summary.materialized_terms == ()
    @test missing_position_summary.deferred_terms == safe_terms
    @test missing_position_summary.total_materialized_count == 0
    @test missing_position_summary.total_skipped_count == 0
    @test !missing_position_summary.term_batch_results_available_on_consumption
    @test !missing_position_summary.term_batch_results_stored_in_summary
    @test !missing_position_summary.materialized
    @test _block_overlap_summary_has_no_matrix_fields(missing_position_summary)

    @test_throws ArgumentError CPBMBlockOverlap._one_body_pair_block_set_consumption_summary(
        (;),
    )

    @test_throws ArgumentError CPBMBlockOverlap._one_body_pair_block_set_consumption(
        plan;
        terms = (:kinetic,),
        inputs,
        materialize_terms = (:overlap,),
    )
    @test_throws ArgumentError CPBMBlockOverlap._one_body_pair_block_set_consumption(
        plan;
        terms = (:overlap, :coulomb),
        inputs,
        materialize_terms = (:coulomb,),
    )
    @test_throws ArgumentError CPBMBlockOverlap._one_body_pair_block_set_consumption(
        plan;
        terms = (:overlap,),
        inputs,
        materialize_terms = (:overlap, :overlap),
    )
end

# Runtime role: tiny contract.
#
# Opt-in overlap/position/x2 materialization for the mixed one-body block-set
# consumer. This delegates only explicit safe terms enabled so far to the
# existing one-term consumer and leaves other requested terms deferred.

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

    @test_throws ArgumentError CPBMBlockOverlap._one_body_pair_block_set_consumption(
        plan;
        terms = (:overlap, :kinetic),
        inputs,
        materialize_terms = (:kinetic,),
    )
    @test_throws ArgumentError CPBMBlockOverlap._one_body_pair_block_set_consumption(
        plan;
        terms = (:kinetic,),
        inputs,
        materialize_terms = (:overlap,),
    )
    @test_throws ArgumentError CPBMBlockOverlap._one_body_pair_block_set_consumption(
        plan;
        terms = (:overlap,),
        inputs,
        materialize_terms = (:overlap, :overlap),
    )
end

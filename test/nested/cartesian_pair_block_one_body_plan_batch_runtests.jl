using Test
using GaussletBases

const CPBMBatchDispatch = GaussletBases.CartesianPairBlockMaterialization
const CPBBatchDispatch = GaussletBases.CartesianCPB
const CTLBatchDispatch = GaussletBases.CartesianTerminalLowering
const CRUBatchDispatch = GaussletBases.CartesianRetainedUnits
const CRTCBatchDispatch = GaussletBases.CartesianRetainedUnitTransformContracts
const CUPBatchDispatch = GaussletBases.CartesianUnitPairs
const CPOPBatchDispatch = GaussletBases.CartesianPairOperatorPlans

function _batch_dispatch_count(counts, field::Symbol, value)
    matches = Tuple(entry for entry in counts if getproperty(entry, field) == value)
    isempty(matches) && return 0
    return only(matches).count
end

function _batch_dispatch_minimal_pair_operator_plan()
    lowering_plan = CTLBatchDispatch.TerminalLoweringPlan(
        CTLBatchDispatch.PQSLowering(q = 3),
        (),
        (),
        (; status = :available_terminal_lowering_plan, materialized = false),
        (; fixture = :one_body_plan_batch),
    )
    retained_plan = CRUBatchDispatch.RetainedUnitPlan(
        CRUBatchDispatch.MetadataOnlyRetainedUnits(),
        lowering_plan,
        (),
        (; status = :available_retained_unit_plan, materialized = false),
        (; fixture = :one_body_plan_batch),
    )
    unit_pair_plan = CUPBatchDispatch.UnitPairPlan(
        CUPBatchDispatch.MetadataOnlyUnitPairs(),
        retained_plan,
        (),
        nothing,
        (; status = :available_unit_pair_plan, materialized = false),
        (; fixture = :one_body_plan_batch),
    )
    transform_plan = CRTCBatchDispatch.RetainedUnitTransformContractPlan(
        CRTCBatchDispatch.MetadataOnlyRetainedUnitTransformContracts(),
        retained_plan,
        (),
        (; status = :available_transform_contract_plan, materialized = false),
        (; fixture = :one_body_plan_batch),
    )
    return CPOPBatchDispatch.PairOperatorPlan(
        CPOPBatchDispatch.MetadataOnlyPairOperatorPlans(),
        unit_pair_plan,
        transform_plan,
        (),
        nothing,
        (; status = :available_pair_operator_plan, materialized = false),
        (; fixture = :one_body_plan_batch),
    )
end

function _batch_dispatch_record(
    pair_key::Tuple{Symbol,Symbol},
    pair_index::Int,
    pair_family::Symbol,
    materialization_path::Symbol;
    metadata = (;),
)
    return CPBMBatchDispatch.PairBlockMaterializationRecord(
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

function _batch_dispatch_direct_record()
    left_source = CPBBatchDispatch.cpb(
        1:1,
        1:2,
        1:1;
        role = :mixed_one_body_batch_direct_left_source,
    )
    right_source = CPBBatchDispatch.cpb(
        2:3,
        1:1,
        1:1;
        role = :mixed_one_body_batch_direct_right_source,
    )
    return _batch_dispatch_record(
        (:direct_left, :direct_right),
        1,
        :direct_direct,
        :direct_direct_pair_block_materialization_pilot;
        metadata = (;
            left_source_cpbs = (left_source,),
            right_source_cpbs = (right_source,),
        ),
    )
end

function _batch_dispatch_plan()
    records = (
        _batch_dispatch_direct_record(),
        _batch_dispatch_record(
            (:pqs_left, :pqs_right),
            2,
            :pqs_pqs,
            :pqs_source_pair_preflight,
        ),
    )
    return CPBMBatchDispatch.PairBlockMaterializationPlan(
        CPBMBatchDispatch.MetadataOnlyPairBlockMaterialization(),
        _batch_dispatch_minimal_pair_operator_plan(),
        records,
        (; status = :available_pair_block_materialization_plan),
        (; fixture = :one_body_plan_batch),
    )
end

function _batch_dispatch_overlap_1d()
    return (;
        x = [
            1.0 0.1 0.2 0.3
            0.1 1.1 0.4 0.5
            0.2 0.4 1.2 0.6
            0.3 0.5 0.6 1.3
        ],
        y = [
            1.4 0.7 0.8 0.9
            0.7 1.5 1.0 1.1
            0.8 1.0 1.6 1.2
            0.9 1.1 1.2 1.7
        ],
        z = [
            1.8 1.3 1.4 1.5
            1.3 1.9 1.6 1.7
            1.4 1.6 2.0 1.8
            1.5 1.7 1.8 2.1
        ],
    )
end

function _batch_dispatch_position_1d()
    return (;
        x = fill(2.0, 4, 4),
        y = [
            3.0 3.1 3.2 3.3
            3.1 3.4 3.5 3.6
            3.2 3.5 3.7 3.8
            3.3 3.6 3.8 3.9
        ],
        z = fill(4.0, 4, 4),
    )
end

@testset "CartesianPairBlockMaterialization mixed one-body plan batch" begin
    plan = _batch_dispatch_plan()
    parent_axis_counts = (4, 4, 4)
    overlap_1d = _batch_dispatch_overlap_1d()
    position_1d = _batch_dispatch_position_1d()

    overlap_batch = CPBMBatchDispatch._one_body_pair_blocks(
        plan,
        :overlap;
        inputs = (; parent_axis_counts, overlap_1d),
    )
    @test overlap_batch isa CPBMBatchDispatch.PairBlockMaterializationBatchResult
    @test overlap_batch.term == :overlap
    @test overlap_batch.materialized_count == 1
    @test overlap_batch.skipped_count == 1
    @test overlap_batch.materialized
    @test overlap_batch.source_operator_blocks_materialized
    @test overlap_batch.final_pair_blocks_materialized
    @test !overlap_batch.operator_blocks_materialized
    @test !overlap_batch.hamiltonian_data_materialized
    @test !overlap_batch.artifacts_materialized
    @test overlap_batch.metadata.materialization_path ==
          :mixed_one_body_pair_block_batch_selector
    @test overlap_batch.metadata.mixed_one_body_dispatcher ==
          :direct_direct_only_plan_dispatcher
    @test overlap_batch.metadata.pair_block_record_count == 2
    @test _batch_dispatch_count(
        overlap_batch.metadata.materialized_selector_family_counts,
        :selector_family,
        :direct_direct,
    ) == 1
    @test _batch_dispatch_count(
        overlap_batch.metadata.skipped_selector_family_counts,
        :selector_family,
        :pqs_source_pair,
    ) == 1
    @test _batch_dispatch_count(
        overlap_batch.metadata.skipped_blocker_counts,
        :blocker,
        :mixed_one_body_selector_family_not_materialized_yet,
    ) == 1

    overlap_result = only(overlap_batch.materialized_results)
    overlap_skip = only(overlap_batch.skipped_records)
    @test overlap_result.term == :overlap
    @test overlap_result.metadata.selector_family == :direct_direct
    @test overlap_skip.selector_family == :pqs_source_pair
    @test overlap_skip.blocker ==
          :mixed_one_body_selector_family_not_materialized_yet
    @test !overlap_skip.materialized

    position_batch = CPBMBatchDispatch._one_body_pair_blocks(
        plan,
        :position_y;
        inputs = (; parent_axis_counts, overlap_1d, position_1d),
    )
    @test position_batch.term == :position_y
    @test position_batch.materialized_count == 1
    @test position_batch.skipped_count == 1
    @test only(position_batch.materialized_results).term == :position_y
    @test only(position_batch.materialized_results).metadata.position_axis == :y
end

@testset "CartesianPairBlockMaterialization mixed one-body plan batch missing inputs" begin
    plan = _batch_dispatch_plan()
    overlap_1d = _batch_dispatch_overlap_1d()
    batch = CPBMBatchDispatch._one_body_pair_blocks(
        plan,
        :kinetic;
        inputs = (; overlap_1d),
    )

    @test batch.materialized_count == 0
    @test batch.skipped_count == 2
    @test !batch.materialized
    @test !batch.source_operator_blocks_materialized
    @test !batch.final_pair_blocks_materialized
    @test _batch_dispatch_count(
        batch.metadata.skipped_blocker_counts,
        :blocker,
        :missing_one_body_factor_inputs,
    ) == 2
    @test _batch_dispatch_count(
        batch.metadata.skipped_selector_family_counts,
        :selector_family,
        :direct_direct,
    ) == 1
    @test _batch_dispatch_count(
        batch.metadata.skipped_selector_family_counts,
        :selector_family,
        :pqs_source_pair,
    ) == 1
end

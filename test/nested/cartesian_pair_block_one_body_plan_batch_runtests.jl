using Test
using GaussletBases

const CPBMBatchDispatch = GaussletBases.CartesianPairBlockMaterialization
const CPBBatchDispatch = GaussletBases.CartesianCPB
const CTLBatchDispatch = GaussletBases.CartesianTerminalLowering
const CRUBatchDispatch = GaussletBases.CartesianRetainedUnits
const CRTCBatchDispatch = GaussletBases.CartesianRetainedUnitTransformContracts
const CUPBatchDispatch = GaussletBases.CartesianUnitPairs
const CPOPBatchDispatch = GaussletBases.CartesianPairOperatorPlans

const _BATCH_DISPATCH_SAFE_TERMS = (
    :overlap,
    :position_x,
    :position_y,
    :position_z,
    :x2_x,
    :x2_y,
    :x2_z,
    :kinetic,
)

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

function _batch_dispatch_pqs_record(; ready::Bool)
    metadata =
        ready ?
        (;
            left_source_mode_dims = (4, 4, 4),
            right_source_mode_dims = (4, 4, 4),
            left_source_mode_count = 64,
            right_source_mode_count = 64,
            source_mode_ordering = :x_major_y_major_z_fast,
            left_source_mode_ordering = :x_major_y_major_z_fast,
            right_source_mode_ordering = :x_major_y_major_z_fast,
        ) :
        (;)
    return _batch_dispatch_record(
        (:pqs_left, :pqs_right),
        2,
        :pqs_pqs,
        :pqs_source_pair_preflight;
        metadata,
    )
end

function _batch_dispatch_lw_record()
    return _batch_dispatch_record(
        (:lw_left, :lw_right),
        3,
        :white_lindsey_boundary_stratum,
        :white_lindsey_boundary_stratum_adapter_preflight,
    )
end

function _batch_dispatch_plan(; pqs_ready::Bool = false, include_lw::Bool = false)
    records = [
        _batch_dispatch_direct_record(),
        _batch_dispatch_pqs_record(; ready = pqs_ready),
    ]
    include_lw && push!(records, _batch_dispatch_lw_record())
    return CPBMBatchDispatch.PairBlockMaterializationPlan(
        CPBMBatchDispatch.MetadataOnlyPairBlockMaterialization(),
        _batch_dispatch_minimal_pair_operator_plan(),
        Tuple(records),
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

function _batch_dispatch_x2_1d()
    return (;
        x = fill(5.0, 4, 4),
        y = [
            6.0 6.1 6.2 6.3
            6.1 6.4 6.5 6.6
            6.2 6.5 6.7 6.8
            6.3 6.6 6.8 6.9
        ],
        z = fill(7.0, 4, 4),
    )
end

function _batch_dispatch_kinetic_1d()
    return (;
        x = fill(8.0, 4, 4),
        y = [
            9.0 9.1 9.2 9.3
            9.1 9.4 9.5 9.6
            9.2 9.5 9.7 9.8
            9.3 9.6 9.8 9.9
        ],
        z = fill(10.0, 4, 4),
    )
end

function _batch_dispatch_inputs(parent_axis_counts, overlap_1d)
    return (;
        parent_axis_counts,
        overlap_1d,
        position_1d = _batch_dispatch_position_1d(),
        x2_1d = _batch_dispatch_x2_1d(),
        kinetic_1d = _batch_dispatch_kinetic_1d(),
    )
end

function _check_batch_dispatch_global_flags(batch, summary)
    @test batch.materialized
    @test batch.source_operator_blocks_materialized
    @test batch.final_pair_blocks_materialized
    @test !batch.operator_blocks_materialized
    @test !batch.hamiltonian_data_materialized
    @test !batch.artifacts_materialized
    @test summary.materialized
    @test summary.source_operator_blocks_materialized
    @test summary.final_pair_blocks_materialized
    @test !summary.global_operator_blocks_materialized
    @test !summary.global_hamiltonian_data_materialized
    @test !summary.global_artifacts_materialized
    @test !summary.factors_constructed
    @test !summary.route_driver_wiring
end

function _check_batch_dispatch_direct_only(batch, term::Symbol)
    summary = CPBMBatchDispatch._one_body_pair_block_batch_summary(batch)
    @test summary.term == term
    @test summary.materialized_count == 1
    @test summary.skipped_count == 1
    @test summary.mixed_one_body_dispatcher ==
          :direct_direct_only_plan_dispatcher
    @test summary.numerical_dispatch_scope == :direct_direct_only
    @test summary.direct_direct_materialized
    @test !summary.pqs_source_pair_materialized
    @test !summary.white_lindsey_materialized
    @test summary.final_local_block_result_count == 1
    @test summary.source_space_only_result_count == 0
    _check_batch_dispatch_global_flags(batch, summary)
    @test _batch_dispatch_count(
        summary.materialized_selector_family_counts,
        :selector_family,
        :direct_direct,
    ) == 1
    @test _batch_dispatch_count(
        summary.skipped_selector_family_counts,
        :selector_family,
        :pqs_source_pair,
    ) == 1
    @test _batch_dispatch_count(
        summary.skipped_blocker_counts,
        :blocker,
        :mixed_one_body_selector_family_not_materialized_yet,
    ) == 1

    direct_result = only(batch.materialized_results)
    pqs_skip = only(batch.skipped_records)
    @test direct_result.term == term
    @test direct_result.metadata.selector_family == :direct_direct
    @test direct_result.final_pair_blocks_materialized
    @test pqs_skip.selector_family == :pqs_source_pair
    @test pqs_skip.blocker ==
          :mixed_one_body_selector_family_not_materialized_yet
    @test !pqs_skip.materialized
end

function _check_batch_dispatch_direct_pqs(batch, term::Symbol)
    summary = CPBMBatchDispatch._one_body_pair_block_batch_summary(batch)
    @test summary.term == term
    @test summary.materialized_count == 2
    @test summary.skipped_count == 1
    @test summary.mixed_one_body_dispatcher ==
          :direct_pqs_source_lw_plan_dispatcher
    @test summary.numerical_dispatch_scope ==
          :direct_direct_pqs_source_pair_and_white_lindsey_boundary_stratum
    @test summary.direct_direct_materialized
    @test summary.pqs_source_pair_materialized
    @test !summary.white_lindsey_materialized
    @test summary.final_local_block_result_count == 1
    @test summary.source_space_only_result_count == 1
    _check_batch_dispatch_global_flags(batch, summary)
    @test _batch_dispatch_count(
        summary.materialized_selector_family_counts,
        :selector_family,
        :direct_direct,
    ) == 1
    @test _batch_dispatch_count(
        summary.materialized_selector_family_counts,
        :selector_family,
        :pqs_source_pair,
    ) == 1
    @test _batch_dispatch_count(
        summary.skipped_selector_family_counts,
        :selector_family,
        :white_lindsey_boundary_stratum,
    ) == 1

    direct_result = only(
        result for result in batch.materialized_results
        if result.metadata.selector_family == :direct_direct
    )
    pqs_result = only(
        result for result in batch.materialized_results
        if result.metadata.selector_family == :pqs_source_pair
    )
    lw_skip = only(batch.skipped_records)
    @test direct_result.term == term
    @test direct_result.final_pair_blocks_materialized
    @test pqs_result.metadata.block_space == :raw_product_source_modes
    @test pqs_result.source_operator_blocks_materialized
    @test !pqs_result.final_pair_blocks_materialized
    @test !pqs_result.metadata.final_pair_blocks_materialized
    @test !pqs_result.metadata.shell_realization_materialized
    @test !pqs_result.metadata.operator_blocks_materialized
    @test !pqs_result.metadata.hamiltonian_data_materialized
    @test !pqs_result.metadata.artifacts_materialized
    @test lw_skip.selector_family == :white_lindsey_boundary_stratum
    @test !lw_skip.materialized
end

@testset "CartesianPairBlockMaterialization mixed one-body direct-only plan batch" begin
    plan = _batch_dispatch_plan()
    parent_axis_counts = (4, 4, 4)
    overlap_1d = _batch_dispatch_overlap_1d()
    inputs = _batch_dispatch_inputs(parent_axis_counts, overlap_1d)

    for term in _BATCH_DISPATCH_SAFE_TERMS
        batch = CPBMBatchDispatch._one_body_pair_blocks(
            plan,
            term;
            inputs,
            materialize_selector_families = (:direct_direct,),
        )
        _check_batch_dispatch_direct_only(batch, term)
    end
end

@testset "CartesianPairBlockMaterialization mixed one-body direct/PQS plan batch" begin
    plan = _batch_dispatch_plan(; pqs_ready = true, include_lw = true)
    parent_axis_counts = (4, 4, 4)
    overlap_1d = _batch_dispatch_overlap_1d()
    inputs = _batch_dispatch_inputs(parent_axis_counts, overlap_1d)

    for term in _BATCH_DISPATCH_SAFE_TERMS
        batch = CPBMBatchDispatch._one_body_pair_blocks(plan, term; inputs)
        _check_batch_dispatch_direct_pqs(batch, term)
    end
end

@testset "CartesianPairBlockMaterialization mixed one-body consumption wrapper" begin
    plan = _batch_dispatch_plan(; pqs_ready = true, include_lw = true)
    parent_axis_counts = (4, 4, 4)
    overlap_1d = _batch_dispatch_overlap_1d()
    inputs = _batch_dispatch_inputs(parent_axis_counts, overlap_1d)

    consumption = CPBMBatchDispatch._one_body_pair_block_consumption(
        plan,
        :x2_z;
        inputs,
    )

    @test consumption.object_kind ==
          :cartesian_pair_block_mixed_one_body_consumption
    @test consumption.status ==
          :partially_materialized_mixed_one_body_pair_block_batch
    @test consumption.term == :x2_z
    @test consumption.batch_result isa
          CPBMBatchDispatch.PairBlockMaterializationBatchResult
    @test consumption.summary.object_kind ==
          :cartesian_pair_block_mixed_one_body_batch_summary
    @test consumption.materialized_count == 2
    @test consumption.skipped_count == 1
    @test consumption.direct_direct_materialized
    @test consumption.pqs_source_pair_materialized
    @test !consumption.white_lindsey_materialized
    @test consumption.source_space_only_result_count == 1
    @test consumption.final_local_block_result_count == 1
    @test consumption.materialized
    @test consumption.source_operator_blocks_materialized
    @test consumption.final_pair_blocks_materialized
    @test !consumption.operator_blocks_materialized
    @test !consumption.hamiltonian_data_materialized
    @test !consumption.artifacts_materialized
    @test !consumption.route_driver_wiring
    @test !consumption.factors_constructed
    @test consumption.materialization_path ==
          :mixed_one_body_pair_block_batch_selector
    @test consumption.mixed_one_body_dispatcher ==
          :direct_pqs_source_lw_plan_dispatcher
    @test consumption.numerical_dispatch_scope ==
          :direct_direct_pqs_source_pair_and_white_lindsey_boundary_stratum
    @test consumption.materialized_count == consumption.summary.materialized_count
    @test consumption.skipped_count == consumption.summary.skipped_count
    @test consumption.operator_blocks_materialized ==
          consumption.summary.global_operator_blocks_materialized
    @test consumption.hamiltonian_data_materialized ==
          consumption.summary.global_hamiltonian_data_materialized
    @test consumption.artifacts_materialized ==
          consumption.summary.global_artifacts_materialized

    @test_throws ArgumentError CPBMBatchDispatch._one_body_pair_block_consumption(
        (;),
        :overlap,
    )
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

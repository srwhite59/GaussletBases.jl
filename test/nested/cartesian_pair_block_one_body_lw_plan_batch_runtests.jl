using Test

include("cartesian_white_lindsey_adapter_fixture_helpers.jl")

const CPBMMixedLWPlan = GaussletBases.CartesianPairBlockMaterialization
const CTLMixedLWPlan = GaussletBases.CartesianTerminalLowering
const CRUMixedLWPlan = GaussletBases.CartesianRetainedUnits
const CRTCMixedLWPlan = GaussletBases.CartesianRetainedUnitTransformContracts
const CUPMixedLWPlan = GaussletBases.CartesianUnitPairs
const COPMixedLWPlan = GaussletBases.CartesianPairOperatorPlans

function _mixed_lw_plan_record(unit_pair)
    return CPBMMixedLWPlan.PairBlockMaterializationRecord(
        unit_pair.pair_key,
        unit_pair.pair_index,
        unit_pair.pair_family,
        :synthetic_white_lindsey_source_operator_path,
        (; left = :white_lindsey_left_transform, right = :white_lindsey_right_transform),
        (; left = :white_lindsey_left_realization, right = :white_lindsey_right_realization),
        :synthetic_white_lindsey_final_block_path,
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
        :white_lindsey_boundary_stratum_adapter_preflight,
        :ready_metadata_only_not_materialized,
        nothing,
        false,
        (; fixture = :mixed_lw_plan_batch),
    )
end

function _mixed_lw_plan(unit_pair)
    retained_plan = CRUMixedLWPlan.RetainedUnitPlan(
        CRUMixedLWPlan.MetadataOnlyRetainedUnits(),
        CTLMixedLWPlan.TerminalLoweringPlan(
            CTLMixedLWPlan.PQSLowering(q = 3),
            (),
            (),
            (; status = :available_terminal_lowering_plan),
            (; fixture = :mixed_lw_plan_batch),
        ),
        (unit_pair.left_unit, unit_pair.right_unit),
        (; status = :available_retained_unit_plan),
        (; fixture = :mixed_lw_plan_batch),
    )
    unit_pair_plan = CUPMixedLWPlan.UnitPairPlan(
        CUPMixedLWPlan.MetadataOnlyUnitPairs(),
        retained_plan,
        (unit_pair,),
        nothing,
        (; status = :available_unit_pair_plan),
        (; fixture = :mixed_lw_plan_batch),
    )
    transform_plan = CRTCMixedLWPlan.RetainedUnitTransformContractPlan(
        CRTCMixedLWPlan.MetadataOnlyRetainedUnitTransformContracts(),
        retained_plan,
        (),
        (; status = :available_transform_contract_plan),
        (; fixture = :mixed_lw_plan_batch),
    )
    pair_operator_plan = COPMixedLWPlan.PairOperatorPlan(
        COPMixedLWPlan.MetadataOnlyPairOperatorPlans(),
        unit_pair_plan,
        transform_plan,
        (),
        nothing,
        (; status = :available_pair_operator_plan),
        (; fixture = :mixed_lw_plan_batch),
    )
    return CPBMMixedLWPlan.PairBlockMaterializationPlan(
        CPBMMixedLWPlan.MetadataOnlyPairBlockMaterialization(),
        pair_operator_plan,
        (_mixed_lw_plan_record(unit_pair),),
        (; status = :available_pair_block_materialization_plan),
        (; fixture = :mixed_lw_plan_batch),
    )
end

function _check_mixed_lw_batch_result(batch, term::Symbol)
    @test batch isa CPBMMixedLWPlan.PairBlockMaterializationBatchResult
    @test batch.term == term
    @test batch.materialized_count == 1
    @test batch.skipped_count == 0
    @test batch.materialized
    @test batch.source_operator_blocks_materialized
    @test batch.final_pair_blocks_materialized
    @test !batch.operator_blocks_materialized
    @test !batch.hamiltonian_data_materialized
    @test !batch.artifacts_materialized
    @test batch.metadata.mixed_one_body_dispatcher ==
          :direct_pqs_source_lw_plan_dispatcher
    @test batch.metadata.numerical_dispatch_scope ==
          :direct_direct_pqs_source_pair_and_white_lindsey_boundary_stratum
    @test !batch.metadata.pqs_source_pair_materialized
    @test batch.metadata.white_lindsey_materialized

    result = only(batch.materialized_results)
    @test result.metadata.selector_family == :white_lindsey_boundary_stratum
    @test result.metadata.numerical_family_selector ==
          :white_lindsey_boundary_stratum_one_body_block
    @test result.materialized
    @test result.source_operator_blocks_materialized
    @test result.final_pair_blocks_materialized
    @test !result.operator_blocks_materialized
    @test !result.hamiltonian_data_materialized
    @test !result.artifacts_materialized
    @test result.metadata.coefficient_source ==
          :white_lindsey_boundary_stratum_pair_unit_coefficients
    @test !result.metadata.operator_blocks_materialized
    @test !result.metadata.hamiltonian_data_materialized
    @test !result.metadata.artifacts_materialized
    return result
end

@testset "CartesianPairBlockMaterialization mixed White-Lindsey plan batch" begin
    fixture = _lw_adapter_prepared_facet_edge_fixture(;
        prefix = "mixed_lw_plan_batch",
    )
    unit_pair = _lw_adapter_unit_pair(
        fixture.real_units.real_facet_unit,
        fixture.real_units.real_edge_unit,
        1,
    )
    plan = _mixed_lw_plan(unit_pair)
    parent_axis_counts = (7, 7, 7)
    overlap_1d = fixture.factors.overlap_1d
    position_1d = fixture.factors.position_1d

    overlap_batch = CPBMMixedLWPlan._one_body_pair_blocks(
        plan,
        :overlap;
        inputs = (; parent_axis_counts, overlap_1d),
    )
    overlap_result = _check_mixed_lw_batch_result(overlap_batch, :overlap)
    @test overlap_result.metadata.materialization_path ==
          :white_lindsey_boundary_stratum_overlap_adapter

    position_batch = CPBMMixedLWPlan._one_body_pair_blocks(
        plan,
        :position_y;
        inputs = (; parent_axis_counts, overlap_1d, position_1d),
    )
    position_result = _check_mixed_lw_batch_result(position_batch, :position_y)
    @test position_result.metadata.position_axis == :y
    @test position_result.metadata.materialization_path ==
          :white_lindsey_boundary_stratum_position_adapter
end

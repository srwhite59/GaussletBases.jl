using Test

include("cartesian_white_lindsey_adapter_fixture_helpers.jl")

const CPBMMixedLWPlan = GaussletBases.CartesianPairBlockMaterialization
const CTLMixedLWPlan = GaussletBases.CartesianTerminalLowering
const CRUMixedLWPlan = GaussletBases.CartesianRetainedUnits
const CRTCMixedLWPlan = GaussletBases.CartesianRetainedUnitTransformContracts
const CUPMixedLWPlan = GaussletBases.CartesianUnitPairs
const COPMixedLWPlan = GaussletBases.CartesianPairOperatorPlans

const _MIXED_LW_PLAN_SAFE_TERMS = (
    :overlap,
    :position_x,
    :position_y,
    :position_z,
    :x2_x,
    :x2_y,
    :x2_z,
    :kinetic,
)

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

function _mixed_lw_plan_count(counts, field::Symbol, value)
    matches = Tuple(entry for entry in counts if getproperty(entry, field) == value)
    isempty(matches) && return 0
    return only(matches).count
end

function _mixed_lw_plan_axis(term::Symbol)
    term in (:position_x, :x2_x) && return :x
    term in (:position_y, :x2_y) && return :y
    term in (:position_z, :x2_z) && return :z
    return nothing
end

function _check_mixed_lw_plan_summary(batch, term::Symbol)
    summary = CPBMMixedLWPlan._one_body_pair_block_batch_summary(batch)
    @test summary.term == term
    @test summary.materialized_count == 1
    @test summary.skipped_count == 0
    @test !summary.direct_direct_materialized
    @test !summary.pqs_source_pair_materialized
    @test summary.white_lindsey_materialized
    @test summary.source_space_only_result_count == 0
    @test summary.final_local_block_result_count == 1
    @test !summary.global_operator_blocks_materialized
    @test !summary.global_hamiltonian_data_materialized
    @test !summary.global_artifacts_materialized
    @test _mixed_lw_plan_count(
        summary.materialized_selector_family_counts,
        :selector_family,
        :white_lindsey_boundary_stratum,
    ) == 1
end

function _check_mixed_lw_plan_term_metadata(result, term::Symbol)
    if term === :overlap
        @test result.metadata.materialization_path ==
              :white_lindsey_boundary_stratum_overlap_adapter
    elseif term in (:position_x, :position_y, :position_z)
        @test result.metadata.materialization_path ==
              :white_lindsey_boundary_stratum_position_adapter
        @test result.metadata.position_axis == _mixed_lw_plan_axis(term)
    elseif term in (:x2_x, :x2_y, :x2_z)
        @test result.metadata.materialization_path ==
              :white_lindsey_boundary_stratum_x2_adapter
        @test result.metadata.x2_axis == _mixed_lw_plan_axis(term)
    elseif term === :kinetic
        @test result.metadata.materialization_path ==
              :white_lindsey_boundary_stratum_kinetic_adapter
        @test result.metadata.kinetic_component_axes == (:x, :y, :z)
    end
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
    x2_1d = fixture.factors.x2_1d
    kinetic_1d = fixture.factors.kinetic_1d
    inputs = (; parent_axis_counts, overlap_1d, position_1d, x2_1d, kinetic_1d)

    for term in _MIXED_LW_PLAN_SAFE_TERMS
        batch = CPBMMixedLWPlan._one_body_pair_blocks(plan, term; inputs)
        result = _check_mixed_lw_batch_result(batch, term)
        _check_mixed_lw_plan_summary(batch, term)
        _check_mixed_lw_plan_term_metadata(result, term)
    end
end

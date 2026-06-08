# Runtime role: tiny contract.
#
# Metadata-only preflight over a mixed one-body term set. This test does not
# call numerical one-term consumption, construct factors, or build the
# White-Lindsey fixture.

using Test
using GaussletBases

const CPBMPreflight = GaussletBases.CartesianPairBlockMaterialization
const CTLPreflight = GaussletBases.CartesianTerminalLowering
const CRUPreflight = GaussletBases.CartesianRetainedUnits
const CRTCPreflight = GaussletBases.CartesianRetainedUnitTransformContracts
const CUPPreflight = GaussletBases.CartesianUnitPairs
const CPOPPreflight = GaussletBases.CartesianPairOperatorPlans

function _preflight_count(counts, field::Symbol, value)
    matches = Tuple(entry for entry in counts if getproperty(entry, field) == value)
    isempty(matches) && return 0
    return only(matches).count
end

function _preflight_pair_operator_plan()
    lowering_plan = CTLPreflight.TerminalLoweringPlan(
        CTLPreflight.PQSLowering(q = 2),
        (),
        (),
        (; status = :available_terminal_lowering_plan, materialized = false),
        (; fixture = :one_body_block_set_preflight),
    )
    retained_plan = CRUPreflight.RetainedUnitPlan(
        CRUPreflight.MetadataOnlyRetainedUnits(),
        lowering_plan,
        (),
        (; status = :available_retained_unit_plan, materialized = false),
        (; fixture = :one_body_block_set_preflight),
    )
    unit_pair_plan = CUPPreflight.UnitPairPlan(
        CUPPreflight.MetadataOnlyUnitPairs(),
        retained_plan,
        (),
        nothing,
        (; status = :available_unit_pair_plan, materialized = false),
        (; fixture = :one_body_block_set_preflight),
    )
    transform_plan = CRTCPreflight.RetainedUnitTransformContractPlan(
        CRTCPreflight.MetadataOnlyRetainedUnitTransformContracts(),
        retained_plan,
        (),
        (; status = :available_transform_contract_plan, materialized = false),
        (; fixture = :one_body_block_set_preflight),
    )
    return CPOPPreflight.PairOperatorPlan(
        CPOPPreflight.MetadataOnlyPairOperatorPlans(),
        unit_pair_plan,
        transform_plan,
        (),
        nothing,
        (; status = :available_pair_operator_plan, materialized = false),
        (; fixture = :one_body_block_set_preflight),
    )
end

function _preflight_record(
    pair_key::Tuple{Symbol,Symbol},
    pair_index::Int,
    pair_family::Symbol,
    materialization_path::Symbol;
    supported_terms = (:overlap, :kinetic),
)
    return CPBMPreflight.PairBlockMaterializationRecord(
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
        (; fixture = :one_body_block_set_preflight),
    )
end

function _preflight_plan()
    records = (
        _preflight_record(
            (:direct_left, :direct_right),
            1,
            :direct_direct,
            :direct_direct_pair_block_materialization_pilot,
        ),
        _preflight_record(
            (:pqs_left, :pqs_right),
            2,
            :pqs_pqs,
            :pqs_source_pair_preflight,
        ),
        _preflight_record(
            (:unsupported_left, :unsupported_right),
            3,
            :unsupported_pair_family,
            :synthetic_unsupported_pair_block_path;
            supported_terms = (:overlap,),
        ),
    )
    return CPBMPreflight.PairBlockMaterializationPlan(
        CPBMPreflight.MetadataOnlyPairBlockMaterialization(),
        _preflight_pair_operator_plan(),
        records,
        (; status = :available_pair_block_materialization_plan),
        (; fixture = :one_body_block_set_preflight),
    )
end

@testset "CartesianPairBlockMaterialization one-body block-set preflight" begin
    plan = _preflight_plan()
    inputs = (;
        parent_axis_counts = (2, 2, 2),
        overlap_1d = :overlap_factors,
        kinetic_1d = :kinetic_factors,
    )

    summary = CPBMPreflight._one_body_pair_block_set_preflight_summary(
        plan;
        terms = (:overlap, :kinetic),
        inputs,
    )
    @test summary.object_kind ==
          :cartesian_pair_block_mixed_one_body_block_set_preflight_summary
    @test summary.status == :partially_ready_mixed_one_body_block_set_preflight
    @test summary.blocker == :blocked_mixed_one_body_dispatch_records
    @test summary.requested_terms == (:overlap, :kinetic)
    @test summary.term_count == 2
    @test summary.plan_record_count == 3
    @test summary.term_set_input_status ==
          :available_one_body_term_set_factor_inputs
    @test isnothing(summary.term_set_input_blocker)
    @test summary.term_set_input_blockers == ()
    @test summary.input_source == :named_tuple
    @test summary.factor_provider_scope == :caller_supplied_or_family_provider
    @test summary.parent_axis_counts_required
    @test summary.parent_axis_counts_status == :available_parent_axis_counts
    @test summary.required_factor_names == (:overlap_1d, :kinetic_1d)
    @test summary.present_factor_names == (:overlap_1d, :kinetic_1d)
    @test summary.missing_factor_names == ()
    @test summary.total_ready_record_count == 4
    @test summary.total_blocked_record_count == 2
    @test _preflight_count(
        summary.selector_family_counts,
        :selector_family,
        :direct_direct,
    ) == 2
    @test _preflight_count(
        summary.selector_family_counts,
        :selector_family,
        :pqs_source_pair,
    ) == 2
    @test _preflight_count(
        summary.selector_family_counts,
        :selector_family,
        :unsupported,
    ) == 2
    @test _preflight_count(
        summary.materialization_path_counts,
        :materialization_path,
        :direct_direct_one_body_selector,
    ) == 2
    @test _preflight_count(
        summary.materialization_path_counts,
        :materialization_path,
        :pqs_source_pair_one_body_selector,
    ) == 2
    @test _preflight_count(
        summary.blocker_counts,
        :blocker,
        :unsupported_pair_block_materialization_path,
    ) == 2
    @test _preflight_count(
        summary.blocker_counts,
        :blocker,
        :unsupported_one_body_term_for_record,
    ) == 1
    @test summary.term_statuses == (
        (;
            term = :overlap,
            status = :partially_ready_mixed_one_body_plan_dispatch,
            ready_count = 2,
            blocked_count = 1,
        ),
        (;
            term = :kinetic,
            status = :partially_ready_mixed_one_body_plan_dispatch,
            ready_count = 2,
            blocked_count = 1,
        ),
    )
    @test !hasproperty(summary, :plan_dispatch_summaries)
    @test !hasproperty(summary.term_summaries[1], :record_summaries)
    @test !hasproperty(summary, :factor_values)
    @test !summary.factor_values_stored
    @test !summary.factors_constructed
    @test !summary.numerical_blocks_materialized
    @test !summary.materialized
    @test !summary.source_operator_blocks_materialized
    @test !summary.final_pair_blocks_materialized
    @test !summary.operator_blocks_materialized
    @test !summary.hamiltonian_data_materialized
    @test !summary.artifacts_materialized
    @test !summary.global_operator_blocks_materialized
    @test !summary.global_hamiltonian_data_materialized
    @test !summary.global_artifacts_materialized
    @test !summary.mixed_dispatcher_materialized
    @test !summary.route_driver_wiring
    @test !summary.coulomb_materialized
    @test !summary.density_density_materialized
    @test !summary.ida_mwg_data_materialized
    @test !summary.pqs_lowdin_materialized
    @test !summary.full_white_lindsey_route_assembled

    provider_values = Dict{Symbol, Any}(
        :parent_axis_counts => (2, 2, 2),
        :overlap_1d => :overlap_factors,
        :kinetic_1d => :kinetic_factors,
    )
    provider = name -> get(provider_values, name, nothing)
    provider_summary = CPBMPreflight._one_body_pair_block_set_preflight_summary(
        plan;
        terms = (:overlap, :kinetic),
        provider,
    )
    @test provider_summary.input_source == :provider_callback
    @test provider_summary.term_set_input_status ==
          :available_one_body_term_set_factor_inputs
    @test provider_summary.total_ready_record_count == 4
    @test provider_summary.total_blocked_record_count == 2

    missing_factor = CPBMPreflight._one_body_pair_block_set_preflight_summary(
        plan;
        terms = (:overlap, :kinetic),
        inputs = (; parent_axis_counts = (2, 2, 2), overlap_1d = :overlap_factors),
    )
    @test missing_factor.status == :blocked_mixed_one_body_block_set_preflight
    @test missing_factor.blocker == :missing_required_one_body_factors
    @test missing_factor.term_set_input_status ==
          :blocked_missing_one_body_term_set_inputs
    @test missing_factor.term_set_input_blockers ==
          (:missing_required_one_body_factors,)
    @test missing_factor.missing_factor_names == (:kinetic_1d,)

    missing_parent = CPBMPreflight._one_body_pair_block_set_preflight_summary(
        plan;
        terms = (:overlap,),
        inputs = (; overlap_1d = :overlap_factors),
    )
    @test missing_parent.status == :blocked_mixed_one_body_block_set_preflight
    @test missing_parent.blocker == :missing_parent_axis_counts
    @test missing_parent.term_set_input_blockers == (:missing_parent_axis_counts,)
    @test missing_parent.parent_axis_counts_required
    @test missing_parent.parent_axis_counts_status == :missing_parent_axis_counts

    @test_throws ArgumentError CPBMPreflight._one_body_pair_block_set_preflight_summary(
        plan;
        terms = (:overlap, :coulomb),
        inputs,
    )
    @test_throws ArgumentError CPBMPreflight._one_body_pair_block_set_preflight_summary(
        (;);
        terms = (:overlap,),
        inputs,
    )
end

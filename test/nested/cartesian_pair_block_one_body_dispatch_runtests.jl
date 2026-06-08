# Runtime role: contract.
#
# Metadata-only mixed dispatch classifier contract. Use for classifier
# semantics; the tiny consumer smoke is preferred for routine mixed-consumer
# edits.

using Test
using GaussletBases

const CPBMOneBodyDispatch = GaussletBases.CartesianPairBlockMaterialization

function _one_body_dispatch_test_record(
    materialization_path::Symbol;
    pair_key = (:left_unit, :right_unit),
    pair_index = 1,
    pair_family = :synthetic_pair_family,
    readiness_status = :ready_metadata_only_not_materialized,
    blocker = nothing,
    supported_terms = (
        :overlap,
        :position_x,
        :position_y,
        :position_z,
        :x2_x,
        :x2_y,
        :x2_z,
        :kinetic,
    ),
)
    return CPBMOneBodyDispatch.PairBlockMaterializationRecord(
        pair_key,
        pair_index,
        pair_family,
        :synthetic_source_operator_path,
        (; left = :synthetic_left_transform, right = :synthetic_right_transform),
        (; left = :synthetic_left_realization, right = :synthetic_right_realization),
        :synthetic_final_block_path,
        supported_terms,
        materialization_path,
        readiness_status,
        blocker,
        false,
        (; fixture = :one_body_dispatch_skeleton_test),
    )
end

function _check_dispatch_nonmaterialization_flags(summary)
    @test !summary.factors_constructed
    @test !summary.numerical_blocks_materialized
    @test !summary.source_operator_blocks_materialized
    @test !summary.final_pair_blocks_materialized
    @test !summary.operator_blocks_materialized
    @test !summary.hamiltonian_data_materialized
    @test !summary.artifacts_materialized
    @test !summary.mixed_dispatcher_materialized
    @test !summary.route_driver_wiring
    @test !summary.pqs_lowdin_materialized
    @test !summary.full_white_lindsey_route_assembled
end

@testset "CartesianPairBlockMaterialization direct one-body dispatch metadata" begin
    record = _one_body_dispatch_test_record(
        :direct_direct_pair_block_materialization_pilot;
        pair_family = :direct_direct,
    )
    summary = CPBMOneBodyDispatch._one_body_pair_block_dispatch_summary(
        record,
        :overlap;
        inputs = (; parent_axis_counts = (3, 4, 5), overlap_1d = :overlap),
    )

    @test summary.object_kind ==
          :cartesian_pair_block_mixed_one_body_dispatch_summary
    @test summary.status == :ready_metadata_only_not_materialized
    @test isnothing(summary.blocker)
    @test summary.blockers == ()
    @test summary.pair_key == (:left_unit, :right_unit)
    @test summary.pair_index == 1
    @test summary.pair_family == :direct_direct
    @test summary.selector_family == :direct_direct
    @test summary.selector_family_status == :available
    @test summary.requested_term == :overlap
    @test summary.term_descriptor.family == :overlap
    @test summary.term_supported_by_record
    @test summary.record_materialization_path ==
          :direct_direct_pair_block_materialization_pilot
    @test summary.materialization_path == :direct_direct_one_body_selector
    @test summary.record_readiness_status ==
          :ready_metadata_only_not_materialized
    @test summary.record_ready
    @test summary.factor_input_status == :available_one_body_factor_inputs
    @test summary.factor_input_blockers == ()
    @test summary.parent_axis_counts_required
    @test summary.parent_axis_counts_status == :available_parent_axis_counts
    @test summary.required_factor_names == (:overlap_1d,)
    @test summary.present_factor_names == (:overlap_1d,)
    @test summary.missing_factor_names == ()
    @test summary.required_factors_available
    @test !summary.unit_pair_required
    @test summary.unit_pair_requirement_status == :not_required
    _check_dispatch_nonmaterialization_flags(summary)
end

@testset "CartesianPairBlockMaterialization PQS source one-body dispatch metadata" begin
    record = _one_body_dispatch_test_record(
        :pqs_source_pair_preflight;
        pair_family = :pqs_pqs,
    )
    summary = CPBMOneBodyDispatch._one_body_pair_block_dispatch_summary(
        record,
        :x2_z;
        inputs = (; overlap_1d = :overlap, x2_1d = :x2),
    )

    @test summary.status == :ready_metadata_only_not_materialized
    @test summary.selector_family == :pqs_source_pair
    @test summary.materialization_path == :pqs_source_pair_one_body_selector
    @test summary.term_descriptor.family == :x2
    @test summary.term_descriptor.axis == :z
    @test summary.factor_input_status == :available_one_body_factor_inputs
    @test !summary.parent_axis_counts_required
    @test summary.parent_axis_counts_status == :not_required
    @test summary.required_factor_names == (:overlap_1d, :x2_1d)
    @test summary.present_factor_names == (:overlap_1d, :x2_1d)
    @test !summary.unit_pair_required
    _check_dispatch_nonmaterialization_flags(summary)
end

@testset "CartesianPairBlockMaterialization White-Lindsey unit-pair requirement" begin
    record = _one_body_dispatch_test_record(
        :white_lindsey_boundary_stratum_adapter_preflight;
        pair_key = (:lw_left, :lw_right),
        pair_family = :white_lindsey_boundary_stratum,
    )
    inputs = (; parent_axis_counts = (5, 5, 5), overlap_1d = :overlap)
    missing_unit_pair = CPBMOneBodyDispatch._one_body_pair_block_dispatch_summary(
        record,
        :overlap;
        inputs,
    )

    @test missing_unit_pair.status == :blocked_mixed_one_body_dispatch
    @test missing_unit_pair.blocker == :missing_white_lindsey_unit_pair
    @test missing_unit_pair.blockers == (:missing_white_lindsey_unit_pair,)
    @test missing_unit_pair.selector_family == :white_lindsey_boundary_stratum
    @test missing_unit_pair.materialization_path ==
          :white_lindsey_boundary_stratum_one_body_selector
    @test missing_unit_pair.factor_input_status ==
          :available_one_body_factor_inputs
    @test missing_unit_pair.parent_axis_counts_required
    @test missing_unit_pair.unit_pair_required
    @test !missing_unit_pair.unit_pair_available
    @test missing_unit_pair.unit_pair_requirement_status ==
          :missing_white_lindsey_unit_pair
    _check_dispatch_nonmaterialization_flags(missing_unit_pair)

    aligned_unit_pair = (; pair_key = (:lw_left, :lw_right))
    aligned = CPBMOneBodyDispatch._one_body_pair_block_dispatch_summary(
        record,
        :overlap;
        inputs,
        unit_pair = aligned_unit_pair,
    )
    @test aligned.status == :ready_metadata_only_not_materialized
    @test aligned.blockers == ()
    @test aligned.unit_pair_available
    @test aligned.unit_pair_key == (:lw_left, :lw_right)
    @test aligned.unit_pair_requirement_status ==
          :available_aligned_white_lindsey_unit_pair

    mismatched = CPBMOneBodyDispatch._one_body_pair_block_dispatch_summary(
        record,
        :overlap;
        inputs,
        unit_pair = (; pair_key = (:lw_left, :other_right)),
    )
    @test mismatched.status == :blocked_mixed_one_body_dispatch
    @test mismatched.blocker == :mismatched_white_lindsey_unit_pair
    @test mismatched.unit_pair_requirement_status ==
          :mismatched_white_lindsey_unit_pair
end

@testset "CartesianPairBlockMaterialization one-body dispatch blockers" begin
    missing_inputs_record = _one_body_dispatch_test_record(
        :direct_direct_pair_block_materialization_pilot;
        pair_family = :direct_direct,
    )
    missing_inputs = CPBMOneBodyDispatch._one_body_pair_block_dispatch_summary(
        missing_inputs_record,
        :kinetic;
        inputs = (; overlap_1d = :overlap),
    )
    @test missing_inputs.status == :blocked_mixed_one_body_dispatch
    @test missing_inputs.blocker == :missing_one_body_factor_inputs
    @test missing_inputs.factor_input_blockers ==
          (:missing_required_one_body_factors, :missing_parent_axis_counts)
    @test missing_inputs.missing_factor_names == (:kinetic_1d,)
    @test missing_inputs.parent_axis_counts_status ==
          :missing_parent_axis_counts
    _check_dispatch_nonmaterialization_flags(missing_inputs)

    unsupported_term_record = _one_body_dispatch_test_record(
        :direct_direct_pair_block_materialization_pilot;
        pair_family = :direct_direct,
        supported_terms = (:overlap,),
    )
    unsupported_term = CPBMOneBodyDispatch._one_body_pair_block_dispatch_summary(
        unsupported_term_record,
        :kinetic;
        inputs = (;
            parent_axis_counts = (2, 2, 2),
            overlap_1d = :overlap,
            kinetic_1d = :kinetic,
        ),
    )
    @test unsupported_term.status == :blocked_mixed_one_body_dispatch
    @test unsupported_term.blocker == :unsupported_one_body_term_for_record
    @test !unsupported_term.term_supported_by_record
    @test unsupported_term.factor_input_status == :available_one_body_factor_inputs

    blocked_record = _one_body_dispatch_test_record(
        :pqs_source_pair_preflight;
        pair_family = :pqs_pqs,
        readiness_status = :blocked_missing_raw_product_source_plan,
        blocker = :missing_left_raw_product_source_plan,
    )
    blocked_summary = CPBMOneBodyDispatch._one_body_pair_block_dispatch_summary(
        blocked_record,
        :overlap;
        inputs = (; overlap_1d = :overlap),
    )
    @test blocked_summary.status == :blocked_mixed_one_body_dispatch
    @test blocked_summary.blocker == :pair_block_record_not_ready
    @test !blocked_summary.record_ready
    @test blocked_summary.record_readiness_status ==
          :blocked_missing_raw_product_source_plan
    @test blocked_summary.record_blocker ==
          :missing_left_raw_product_source_plan

    unsupported_record = _one_body_dispatch_test_record(
        :deferred_pair_block_materialization_path;
        pair_family = :unsupported_family,
    )
    unsupported_summary = CPBMOneBodyDispatch._one_body_pair_block_dispatch_summary(
        unsupported_record,
        :overlap;
        inputs = (; overlap_1d = :overlap),
    )
    @test unsupported_summary.status == :blocked_mixed_one_body_dispatch
    @test unsupported_summary.blocker ==
          :unsupported_pair_block_materialization_path
    @test unsupported_summary.selector_family == :unsupported
    @test unsupported_summary.materialization_path === nothing
    @test unsupported_summary.materialization_path_status == :unknown

    @test_throws ArgumentError CPBMOneBodyDispatch._one_body_pair_block_dispatch_summary(
        (;),
        :overlap;
        inputs = (; overlap_1d = :overlap),
    )
end

# Runtime role: contract.
#
# Focused caller-supplied factor input contract. Use for input convention
# changes; the tiny consumer smoke is preferred for routine mixed-consumer
# edits.

using Test
using GaussletBases

const CPBMOneBodyInputs = GaussletBases.CartesianPairBlockMaterialization

@testset "CartesianPairBlockMaterialization one-body named-tuple factor inputs" begin
    inputs = (;
        parent_axis_counts = (3, 4, 5),
        overlap_1d = :overlap_factors,
        position_1d = :position_factors,
        x2_1d = :x2_factors,
        kinetic_1d = :kinetic_factors,
    )
    summary = CPBMOneBodyInputs._one_body_factor_input_summary(
        :position_y;
        inputs,
        selector_family = :direct_direct,
        parent_axis_counts_required = true,
    )

    @test summary.object_kind ==
          :cartesian_pair_block_one_body_factor_input_summary
    @test summary.status == :available_one_body_factor_inputs
    @test isnothing(summary.blocker)
    @test summary.blockers == ()
    @test summary.requested_term == :position_y
    @test summary.term_descriptor.requested_term == :position_y
    @test summary.term_descriptor.family == :position
    @test summary.term_descriptor.axis == :y
    @test summary.selector_family == :direct_direct
    @test summary.selector_family_status == :specified
    @test summary.materialization_path == :direct_direct_one_body_selector
    @test summary.materialization_path_status == :available
    @test summary.input_source == :named_tuple
    @test summary.parent_axis_counts_required
    @test summary.parent_axis_counts_status == :available_parent_axis_counts
    @test summary.parent_axis_counts == (3, 4, 5)
    @test summary.required_factor_names == (:overlap_1d, :position_1d)
    @test summary.present_factor_names == (:overlap_1d, :position_1d)
    @test summary.missing_factor_names == ()
    @test summary.required_factors_available
    @test hasproperty(summary.factor_values, :overlap_1d)
    @test hasproperty(summary.factor_values, :position_1d)
    @test summary.factor_values.overlap_1d == :overlap_factors
    @test summary.factor_values.position_1d == :position_factors
    @test !summary.factors_constructed
    @test !summary.numerical_blocks_materialized
    @test !summary.mixed_dispatcher_materialized
    @test !summary.route_driver_wiring
    @test !summary.hamiltonian_data_materialized
    @test !summary.artifacts_materialized
end

@testset "CartesianPairBlockMaterialization one-body provider factor inputs" begin
    values = Dict{Symbol, Any}(
        :overlap_1d => :overlap_factors,
        :x2_1d => :x2_factors,
    )
    provider = name -> get(values, name, nothing)
    summary = CPBMOneBodyInputs._one_body_factor_input_summary(
        :x2_z;
        provider,
        selector_family = :pqs_source_pair,
        parent_axis_counts_required = false,
    )

    @test summary.status == :available_one_body_factor_inputs
    @test summary.input_source == :provider_callback
    @test summary.parent_axis_counts_status == :not_required
    @test summary.parent_axis_counts === nothing
    @test summary.materialization_path == :pqs_source_pair_one_body_selector
    @test summary.required_factor_names == (:overlap_1d, :x2_1d)
    @test summary.present_factor_names == (:overlap_1d, :x2_1d)
    @test summary.missing_factor_names == ()
    @test summary.required_factors_available
    @test hasproperty(summary.factor_values, :overlap_1d)
    @test hasproperty(summary.factor_values, :x2_1d)
    @test !summary.factors_constructed
    @test !summary.numerical_blocks_materialized
end

@testset "CartesianPairBlockMaterialization one-body term-set factor inputs" begin
    complete_inputs = (;
        parent_axis_counts = (3, 4, 5),
        overlap_1d = :overlap_factors,
        position_1d = :position_factors,
        x2_1d = :x2_factors,
        kinetic_1d = :kinetic_factors,
    )

    default_summary = CPBMOneBodyInputs._one_body_term_set_factor_input_summary(
        inputs = complete_inputs,
        parent_axis_counts_required = true,
    )
    @test default_summary.object_kind ==
          :cartesian_pair_block_one_body_term_set_factor_input_summary
    @test default_summary.status == :available_one_body_term_set_factor_inputs
    @test isnothing(default_summary.blocker)
    @test default_summary.blockers == ()
    @test default_summary.requested_terms == (
        :overlap,
        :position_x,
        :position_y,
        :position_z,
        :x2_x,
        :x2_y,
        :x2_z,
        :kinetic,
    )
    @test default_summary.term_count == 8
    @test default_summary.required_factor_names ==
          (:overlap_1d, :position_1d, :x2_1d, :kinetic_1d)
    @test default_summary.present_factor_names ==
          (:overlap_1d, :position_1d, :x2_1d, :kinetic_1d)
    @test default_summary.missing_factor_names == ()
    @test default_summary.required_factors_available
    @test default_summary.input_source == :named_tuple
    @test default_summary.parent_axis_counts_required
    @test default_summary.parent_axis_counts_status ==
          :available_parent_axis_counts
    @test default_summary.parent_axis_counts == (3, 4, 5)
    @test default_summary.factor_provider_scope ==
          :caller_supplied_or_family_provider
    @test !hasproperty(default_summary, :factor_values)
    @test !default_summary.factor_values_stored
    @test !default_summary.factors_constructed
    @test !default_summary.numerical_blocks_materialized
    @test !default_summary.mixed_dispatcher_materialized
    @test !default_summary.route_driver_wiring
    @test !default_summary.global_operator_blocks_materialized
    @test !default_summary.hamiltonian_data_materialized
    @test !default_summary.artifacts_materialized
    @test !default_summary.coulomb_materialized
    @test !default_summary.density_density_materialized
    @test !default_summary.ida_mwg_data_materialized
    @test !default_summary.pqs_lowdin_materialized
    @test !default_summary.full_white_lindsey_route_assembled

    custom_summary = CPBMOneBodyInputs._one_body_term_set_factor_input_summary(
        (:kinetic, :overlap, :position_y);
        inputs = complete_inputs,
    )
    @test custom_summary.requested_terms == (:kinetic, :overlap, :position_y)
    @test custom_summary.term_count == 3
    @test custom_summary.required_factor_names ==
          (:overlap_1d, :kinetic_1d, :position_1d)
    @test custom_summary.present_factor_names ==
          (:overlap_1d, :kinetic_1d, :position_1d)
    @test custom_summary.parent_axis_counts_status ==
          :available_parent_axis_counts_not_required

    values = Dict{Symbol, Any}(
        :overlap_1d => :overlap_factors,
        :x2_1d => :x2_factors,
        :parent_axis_counts => (2, 2, 2),
    )
    provider = name -> get(values, name, nothing)
    provider_summary = CPBMOneBodyInputs._one_body_term_set_factor_input_summary(
        (:x2_z, :overlap);
        provider,
        parent_axis_counts_required = true,
    )
    @test provider_summary.status == :available_one_body_term_set_factor_inputs
    @test provider_summary.input_source == :provider_callback
    @test provider_summary.required_factor_names == (:overlap_1d, :x2_1d)
    @test provider_summary.present_factor_names == (:overlap_1d, :x2_1d)
    @test provider_summary.missing_factor_names == ()
    @test provider_summary.parent_axis_counts_status ==
          :available_parent_axis_counts
    @test provider_summary.parent_axis_counts == (2, 2, 2)

    missing_factor = CPBMOneBodyInputs._one_body_term_set_factor_input_summary(
        (:overlap, :kinetic);
        inputs = (; overlap_1d = :overlap_factors),
    )
    @test missing_factor.status == :blocked_missing_one_body_term_set_inputs
    @test missing_factor.blocker == :missing_required_one_body_factors
    @test missing_factor.blockers == (:missing_required_one_body_factors,)
    @test missing_factor.required_factor_names == (:overlap_1d, :kinetic_1d)
    @test missing_factor.present_factor_names == (:overlap_1d,)
    @test missing_factor.missing_factor_names == (:kinetic_1d,)
    @test !missing_factor.required_factors_available
    @test missing_factor.parent_axis_counts_status == :not_required

    missing_parent = CPBMOneBodyInputs._one_body_term_set_factor_input_summary(
        (:overlap,);
        inputs = (; overlap_1d = :overlap_factors),
        parent_axis_counts_required = true,
    )
    @test missing_parent.status == :blocked_missing_one_body_term_set_inputs
    @test missing_parent.blocker == :missing_parent_axis_counts
    @test missing_parent.blockers == (:missing_parent_axis_counts,)
    @test missing_parent.required_factors_available
    @test missing_parent.parent_axis_counts_status == :missing_parent_axis_counts

    @test_throws ArgumentError CPBMOneBodyInputs._one_body_term_set_factor_input_summary(
        (:overlap, :coulomb);
        inputs = complete_inputs,
    )
    @test_throws ArgumentError CPBMOneBodyInputs._one_body_term_set_factor_input_summary(
        (:overlap,);
        inputs = (; overlap_1d = :overlap_factors),
        provider = name -> nothing,
    )
end

@testset "CartesianPairBlockMaterialization one-body missing inputs" begin
    missing_factors = CPBMOneBodyInputs._one_body_factor_input_summary(
        :kinetic;
        inputs = (; overlap_1d = :overlap_factors),
        selector_family = :white_lindsey_boundary_stratum,
        parent_axis_counts_required = true,
    )

    @test missing_factors.status == :blocked_missing_one_body_inputs
    @test missing_factors.blocker == :missing_required_one_body_factors
    @test missing_factors.blockers ==
          (:missing_required_one_body_factors, :missing_parent_axis_counts)
    @test missing_factors.materialization_path ==
          :white_lindsey_boundary_stratum_one_body_selector
    @test missing_factors.parent_axis_counts_status ==
          :missing_parent_axis_counts
    @test missing_factors.required_factor_names == (:overlap_1d, :kinetic_1d)
    @test missing_factors.present_factor_names == (:overlap_1d,)
    @test missing_factors.missing_factor_names == (:kinetic_1d,)
    @test !missing_factors.required_factors_available
    @test hasproperty(missing_factors.factor_values, :overlap_1d)
    @test !hasproperty(missing_factors.factor_values, :kinetic_1d)
    @test !missing_factors.factors_constructed
    @test !missing_factors.numerical_blocks_materialized

    @test_throws ArgumentError CPBMOneBodyInputs._one_body_factor_input_summary(
        :overlap;
        inputs = (; overlap_1d = :overlap_factors),
        provider = name -> nothing,
    )
    @test_throws ArgumentError CPBMOneBodyInputs._one_body_factor_input_summary(
        "overlap";
        inputs = (; overlap_1d = :overlap_factors),
    )
end

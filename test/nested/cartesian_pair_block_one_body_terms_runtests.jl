# Runtime role: contract.
#
# Focused metadata contract for supported one-body term descriptors. Use after
# selector-surface changes; the tiny consumer smoke is preferred for routine
# mixed-consumer edits.

using Test
using GaussletBases

const CPBMOneBodyTerms = GaussletBases.CartesianPairBlockMaterialization

@testset "CartesianPairBlockMaterialization one-body term descriptors" begin
    expected = (;
        overlap = (;
            family = :overlap,
            term_kind = :one_body_overlap,
            axis = nothing,
            axis_index = nothing,
            required_factor_roles = (:overlap,),
            required_factor_names = (:overlap_1d,),
        ),
        position_x = (;
            family = :position,
            term_kind = :one_body_axis_position,
            axis = :x,
            axis_index = 1,
            required_factor_roles = (:overlap, :position),
            required_factor_names = (:overlap_1d, :position_1d),
        ),
        position_y = (;
            family = :position,
            term_kind = :one_body_axis_position,
            axis = :y,
            axis_index = 2,
            required_factor_roles = (:overlap, :position),
            required_factor_names = (:overlap_1d, :position_1d),
        ),
        position_z = (;
            family = :position,
            term_kind = :one_body_axis_position,
            axis = :z,
            axis_index = 3,
            required_factor_roles = (:overlap, :position),
            required_factor_names = (:overlap_1d, :position_1d),
        ),
        x2_x = (;
            family = :x2,
            term_kind = :one_body_axis_x2,
            axis = :x,
            axis_index = 1,
            required_factor_roles = (:overlap, :x2),
            required_factor_names = (:overlap_1d, :x2_1d),
        ),
        x2_y = (;
            family = :x2,
            term_kind = :one_body_axis_x2,
            axis = :y,
            axis_index = 2,
            required_factor_roles = (:overlap, :x2),
            required_factor_names = (:overlap_1d, :x2_1d),
        ),
        x2_z = (;
            family = :x2,
            term_kind = :one_body_axis_x2,
            axis = :z,
            axis_index = 3,
            required_factor_roles = (:overlap, :x2),
            required_factor_names = (:overlap_1d, :x2_1d),
        ),
        kinetic = (;
            family = :kinetic,
            term_kind = :one_body_cartesian_kinetic_sum,
            axis = nothing,
            axis_index = nothing,
            required_factor_roles = (:overlap, :kinetic),
            required_factor_names = (:overlap_1d, :kinetic_1d),
        ),
    )
    supported_terms = (
        :overlap,
        :position_x,
        :position_y,
        :position_z,
        :x2_x,
        :x2_y,
        :x2_z,
        :kinetic,
    )

    for term in supported_terms
        descriptor = CPBMOneBodyTerms._one_body_term_descriptor(term)
        expected_term = getproperty(expected, term)
        @test descriptor.object_kind ==
              :cartesian_pair_block_one_body_term_descriptor
        @test descriptor.status == :available_one_body_term_descriptor
        @test descriptor.requested_term == term
        @test descriptor.term == term
        @test descriptor.family == expected_term.family
        @test descriptor.term_family == expected_term.family
        @test descriptor.term_kind == expected_term.term_kind
        @test descriptor.axis == expected_term.axis
        @test descriptor.axis_index == expected_term.axis_index
        @test descriptor.required_factor_roles ==
              expected_term.required_factor_roles
        @test descriptor.required_factor_names ==
              expected_term.required_factor_names
        @test descriptor.factor_provider_scope ==
              :caller_supplied_or_family_provider
        @test !descriptor.factors_constructed
        @test !descriptor.numerical_blocks_materialized
        @test !descriptor.mixed_dispatcher_materialized
        @test !descriptor.route_driver_wiring
        @test !descriptor.hamiltonian_data_materialized
        @test !descriptor.artifacts_materialized
    end

    @test_throws ArgumentError CPBMOneBodyTerms._one_body_term_descriptor(:coulomb)
    @test_throws ArgumentError CPBMOneBodyTerms._one_body_term_descriptor("overlap")
end

@testset "CartesianPairBlockMaterialization one-body term-set descriptor" begin
    default_terms = (
        :overlap,
        :position_x,
        :position_y,
        :position_z,
        :x2_x,
        :x2_y,
        :x2_z,
        :kinetic,
    )

    descriptor = CPBMOneBodyTerms._one_body_term_set_descriptor()
    @test descriptor.object_kind ==
          :cartesian_pair_block_one_body_term_set_descriptor
    @test descriptor.status == :available_one_body_term_set_descriptor
    @test descriptor.requested_terms == default_terms
    @test descriptor.terms == default_terms
    @test descriptor.term_count == 8
    @test length(descriptor.term_descriptors) == 8
    @test getproperty.(descriptor.term_descriptors, :requested_term) ==
          default_terms
    @test descriptor.term_families == (
        :overlap,
        :position,
        :position,
        :position,
        :x2,
        :x2,
        :x2,
        :kinetic,
    )
    @test descriptor.required_factor_names ==
          (:overlap_1d, :position_1d, :x2_1d, :kinetic_1d)
    @test descriptor.result_terms_remain_separated
    @test !descriptor.block_set_results_summed
    @test descriptor.factor_provider_scope == :caller_supplied_or_family_provider
    @test !descriptor.factors_constructed
    @test !descriptor.numerical_blocks_materialized
    @test !descriptor.mixed_dispatcher_materialized
    @test !descriptor.route_driver_wiring
    @test !descriptor.global_operator_blocks_materialized
    @test !descriptor.hamiltonian_data_materialized
    @test !descriptor.artifacts_materialized
    @test !descriptor.coulomb_materialized
    @test !descriptor.density_density_materialized
    @test !descriptor.ida_mwg_data_materialized
    @test !descriptor.pqs_lowdin_materialized
    @test !descriptor.full_white_lindsey_route_assembled

    custom = CPBMOneBodyTerms._one_body_term_set_descriptor(
        (:kinetic, :overlap, :position_y),
    )
    @test custom.terms == (:kinetic, :overlap, :position_y)
    @test custom.term_count == 3
    @test getproperty.(custom.term_descriptors, :requested_term) ==
          (:kinetic, :overlap, :position_y)
    @test custom.required_factor_names ==
          (:overlap_1d, :kinetic_1d, :position_1d)

    singleton = CPBMOneBodyTerms._one_body_term_set_descriptor(:overlap)
    @test singleton.terms == (:overlap,)
    @test singleton.term_count == 1

    vector_terms = CPBMOneBodyTerms._one_body_term_set_descriptor([:x2_z, :x2_x])
    @test vector_terms.terms == (:x2_z, :x2_x)

    @test_throws ArgumentError CPBMOneBodyTerms._one_body_term_set_descriptor(())
    @test_throws ArgumentError CPBMOneBodyTerms._one_body_term_set_descriptor(
        (:overlap, :coulomb),
    )
    @test_throws ArgumentError CPBMOneBodyTerms._one_body_term_set_descriptor(
        (:density_density,),
    )
    @test_throws ArgumentError CPBMOneBodyTerms._one_body_term_set_descriptor(
        (:nuclear_attraction,),
    )
    @test_throws ArgumentError CPBMOneBodyTerms._one_body_term_set_descriptor(
        (:gaussian_local,),
    )
    @test_throws ArgumentError CPBMOneBodyTerms._one_body_term_set_descriptor(
        "overlap",
    )
end

@testset "CartesianPairBlockMaterialization one-body selector surface audit" begin
    summary = CPBMOneBodyTerms._one_body_selector_surface_summary()
    @test summary.object_kind ==
          :cartesian_pair_block_one_body_selector_surface_summary
    @test summary.status ==
          :available_internal_one_body_selector_surface_audit
    @test summary.supported_terms == (
        :overlap,
        :position_x,
        :position_y,
        :position_z,
        :x2_x,
        :x2_y,
        :x2_z,
        :kinetic,
    )
    @test summary.selector_family_count == 3
    @test summary.selector_families ==
          (:direct_direct, :pqs_source_pair, :white_lindsey_boundary_stratum)
    @test summary.record_selectors.direct_direct ==
          :direct_direct_one_body_block
    @test summary.record_selectors.pqs_source_pair ==
          :pqs_source_pair_one_body_block
    @test summary.record_selectors.white_lindsey_boundary_stratum ==
          :white_lindsey_boundary_stratum_one_body_block
    @test summary.batch_selectors.direct_direct ==
          :direct_direct_one_body_blocks
    @test summary.batch_selectors.pqs_source_pair ==
          :pqs_source_pair_one_body_blocks
    @test summary.batch_selectors.white_lindsey_boundary_stratum ==
          :white_lindsey_boundary_stratum_one_body_blocks
    @test summary.factor_provider_scope == :caller_supplied_or_family_provider
    @test !summary.factors_constructed
    @test !summary.numerical_blocks_materialized
    @test !summary.mixed_dispatcher_materialized
    @test !summary.route_driver_wiring
    @test !summary.hamiltonian_data_materialized
    @test !summary.artifacts_materialized
    @test !summary.coulomb_materialized
    @test !summary.ida_mwg_data_materialized
    @test !summary.pqs_lowdin_materialized
    @test !summary.full_white_lindsey_route_assembled
end

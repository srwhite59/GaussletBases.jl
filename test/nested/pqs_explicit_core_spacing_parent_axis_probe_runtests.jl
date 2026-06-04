using Test
using GaussletBases

@testset "PQS explicit-core-spacing parent-axis probe" begin
    metrics_module = GaussletBases.CartesianContractedParentMetrics

    unavailable_setup =
        metrics_module._pqs_standard_source_box_route_setup(
            nuclear_charges = (4, 4),
            atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
            q = 5,
            radius = 3.0,
            q_to_core_spacing_rule = :explicit_core_spacing_only,
        )
    unavailable_probe =
        metrics_module._pqs_explicit_core_spacing_parent_axis_probe(unavailable_setup)
    @test unavailable_probe.object_kind == :pqs_explicit_core_spacing_parent_axis_probe
    @test unavailable_probe.status == :not_constructed_pending_facts
    @test unavailable_probe.axis_lengths === nothing
    @test !unavailable_probe.parent_axis_metadata_constructed
    @test !unavailable_probe.carry_objects_requested
    @test !unavailable_probe.basis_object_available
    @test !unavailable_probe.axis_bundle_object_available
    @test unavailable_probe.basis_object === nothing
    @test unavailable_probe.axis_bundle_object === nothing
    @test :explicit_core_spacing in unavailable_probe.pending_facts
    @test !unavailable_probe.diagnostics.explicit_spacing_probe_only
    @test !unavailable_probe.diagnostics.default_standard_rule
    @test unavailable_probe.core_spacing_source == :unavailable

    standard_setup =
        metrics_module._pqs_standard_source_box_route_setup(
            nuclear_charges = (4, 4),
            atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
            q = 5,
            radius = 3.0,
        )
    standard_probe =
        metrics_module._pqs_explicit_core_spacing_parent_axis_probe(standard_setup)
    @test standard_probe.status ==
          :constructed_explicit_core_spacing_parent_axis_metadata
    @test standard_probe.parent_axis_metadata_constructed
    @test !standard_probe.carry_objects_requested
    @test !standard_probe.basis_object_available
    @test !standard_probe.axis_bundle_object_available
    @test standard_probe.basis_object === nothing
    @test standard_probe.axis_bundle_object === nothing
    @test standard_probe.core_spacing == 0.15
    @test standard_probe.core_spacing_source == :standard_n_s_default
    @test !standard_probe.explicit_spacing_probe_only
    @test standard_probe.default_standard_rule

    explicit_setup =
        metrics_module._pqs_standard_source_box_route_setup(
            nuclear_charges = (4, 4),
            atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
            q = 5,
            radius = 3.0,
            core_spacing = 0.15,
        )
    explicit_probe =
        metrics_module._pqs_explicit_core_spacing_parent_axis_probe(explicit_setup)

    @test explicit_probe.status ==
          :constructed_explicit_core_spacing_parent_axis_metadata
    @test explicit_probe.parent_axis_metadata_constructed
    @test explicit_probe.axis_bundle_metadata.status == :constructed
    @test explicit_probe.axis_bundle_metadata.object_kind == :_CartesianNestedAxisBundles3D
    @test explicit_probe.basis_metadata.object_kind == :BondAlignedDiatomicQWBasis3D
    @test !explicit_probe.carry_objects_requested
    @test !explicit_probe.basis_object_available
    @test !explicit_probe.axis_bundle_object_available
    @test explicit_probe.basis_object === nothing
    @test explicit_probe.axis_bundle_object === nothing
    @test explicit_probe.basis_metadata.constructor == :bond_aligned_homonuclear_qw_basis
    @test explicit_probe.basis_metadata.bond_axis == :x
    @test explicit_probe.basis_metadata.bond_length == 4.0
    @test explicit_probe.physical_extent_inputs.xmax_parallel == 5.0
    @test explicit_probe.physical_extent_inputs.xmax_transverse == 3.0
    @test explicit_probe.core_spacing == 0.15
    @test explicit_probe.reference_spacing == 1.0
    @test explicit_probe.tail_spacing == 10.0
    @test explicit_probe.gausslet_backend == :pgdg_localized_experimental
    @test explicit_probe.gausslet_backend_role == :default_development_pgdg_localized
    @test explicit_probe.expansion_source ==
          :default_coulomb_gaussian_expansion_doacc_false
    @test explicit_probe.axis_lengths isa NTuple{3,Int}
    @test all(length_value -> length_value > 0, explicit_probe.axis_lengths)
    @test explicit_probe.pending_facts == ()
    @test explicit_probe.diagnostics.axis_lengths == explicit_probe.axis_lengths
    @test explicit_probe.explicit_spacing_probe_only
    @test !explicit_probe.default_standard_rule
    @test explicit_probe.core_spacing_source == :explicit_core_spacing_override

    carried_probe =
        metrics_module._pqs_explicit_core_spacing_parent_axis_probe(
            explicit_setup;
            carry_objects = true,
        )
    @test carried_probe.status ==
          :constructed_explicit_core_spacing_parent_axis_metadata
    @test carried_probe.carry_objects_requested
    @test carried_probe.basis_object_available
    @test carried_probe.axis_bundle_object_available
    @test carried_probe.basis_object !== nothing
    @test carried_probe.axis_bundle_object !== nothing
    @test carried_probe.basis_object_type_label == "BondAlignedDiatomicQWBasis3D"
    @test carried_probe.axis_bundle_object_type_label ==
          "_CartesianNestedAxisBundles3D"

    numerical_probe =
        metrics_module._pqs_explicit_core_spacing_parent_axis_probe(
            explicit_setup;
            gausslet_backend = :numerical_reference,
        )
    @test numerical_probe.parent_axis_metadata_constructed
    @test numerical_probe.gausslet_backend == :numerical_reference
    @test numerical_probe.gausslet_backend_role ==
          :explicit_non_default_numerical_reference

    disabled_probe =
        metrics_module._pqs_explicit_core_spacing_parent_axis_probe(
            explicit_setup;
            construct_axis_bundles = false,
        )
    @test disabled_probe.status == :not_constructed_pending_facts
    @test !disabled_probe.parent_axis_metadata_constructed
    @test :probe_parent_axis_construction_flag in disabled_probe.pending_facts

    heteronuclear_setup =
        metrics_module._pqs_standard_source_box_route_setup(
            nuclear_charges = (4, 3),
            atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
            q = 5,
            radius = 3.0,
            core_spacing = 0.15,
        )
    heteronuclear_probe =
        metrics_module._pqs_explicit_core_spacing_parent_axis_probe(heteronuclear_setup)
    @test heteronuclear_probe.status == :not_constructed_pending_facts
    @test !heteronuclear_probe.parent_axis_metadata_constructed
    @test :homonuclear_setup in heteronuclear_probe.pending_facts

    for diagnostics in (
        unavailable_probe.diagnostics,
        standard_probe.diagnostics,
        explicit_probe.diagnostics,
        carried_probe.diagnostics,
        disabled_probe.diagnostics,
        heteronuclear_probe.diagnostics,
    )
        @test diagnostics.private_development_only
        @test !diagnostics.production_route
        @test !diagnostics.public_default_consumes
        @test !diagnostics.packet_adoption
        @test !diagnostics.fixed_block_routing
        @test !diagnostics.qwhamiltonian_consumes
        @test !diagnostics.hamiltonian_matrix_built
        @test !diagnostics.shell_projection_used
        @test !diagnostics.lowdin_cleanup_used
        @test !diagnostics.support_local_shell_row_algorithm
        @test !diagnostics.support_coefficient_matrix_used
        @test !diagnostics.retained_pqs_weights_used
        @test !diagnostics.retained_weight_division_allowed
        @test !diagnostics.repo_side_ray_id
        @test !diagnostics.mwg_ida_semantics_changed
        @test !diagnostics.ecp_terms_implemented
        @test !diagnostics.cr2_science_status_changed
    end
end

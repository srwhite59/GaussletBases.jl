using LinearAlgebra
using Test
using GaussletBases
using JLD2

@testset "White-Lindsey materialized complete-shell seed fixture" begin
    density_expansion = coulomb_gaussian_expansion(doacc = false)
    fixture = GaussletBases._white_lindsey_low_order_materialized_seed_fixture(
        expansion = density_expansion,
    )
    sequence = fixture.sequence
    fixed_block = fixture.fixed_block
    structure = fixture.structure
    fixture_shellization_summary = fixture.shellization_summary
    inventory = fixture.inventory
    route_units = GaussletBases._white_lindsey_low_order_materialized_seed_route_units(fixture)
    operator_inventory =
        GaussletBases._white_lindsey_low_order_materialized_seed_operator_inventory(fixture)

    shell_retained_count = 98
    core_retained_count = 5^3
    total_retained_dimension = core_retained_count + shell_retained_count

    @test fixture.packet_kernel == :factorized_direct
    @test sequence isa GaussletBases._CartesianNestedShellSequence3D
    @test length(sequence.shell_layers) == 1
    shell = only(sequence.shell_layers)
    @test shell isa GaussletBases._CartesianNestedCompleteShell3D

    shellization_summary = GaussletBases._cartesian_shellization_route_summary(sequence)
    @test shellization_summary.object_kind == :cartesian_shellization_route_summary
    @test shellization_summary.status == :private_development_summary
    @test shellization_summary.private_development_only
    @test shellization_summary.route_family == :one_center_full_parent
    @test shellization_summary.source_kind == :one_center_shell_sequence
    @test shellization_summary.shellization_stage == :route_neutral_spatial_planning
    @test shellization_summary.lowering_stage == :not_lowered_by_shellization_summary
    @test shellization_summary.parent_box == sequence.working_box
    @test shellization_summary.working_box == sequence.working_box
    @test shellization_summary.core_column_range == sequence.core_column_range
    @test shellization_summary.core_retained_count == length(sequence.core_column_range)
    @test shellization_summary.shell_layer_count == length(sequence.shell_layers)
    @test shellization_summary.shell_layer_kinds == (:_CartesianNestedCompleteShell3D,)
    @test shellization_summary.shell_layer_column_ranges == Tuple(sequence.layer_column_ranges)
    @test shellization_summary.retained_dimension == size(sequence.coefficient_matrix, 2)
    @test shellization_summary.support_count == length(sequence.support_indices)
    @test shellization_summary.split_status == :not_applicable
    @test shellization_summary.bond_axis === nothing
    @test !shellization_summary.midpoint_slab_present
    @test shellization_summary.child_sequence_count == 0
    @test shellization_summary.shared_shell_layer_count == length(sequence.shell_layers)
    @test isempty(shellization_summary.child_column_ranges)
    @test shellization_summary.contact_or_merge_status == :not_applicable
    @test shellization_summary.diagnostics.route_neutral_spatial_planning
    @test !shellization_summary.diagnostics.lowering_applied_by_summary
    @test !shellization_summary.diagnostics.white_lindsey_lowering_adopted_by_summary
    @test !shellization_summary.diagnostics.pqs_lowering_adopted_by_summary
    @test !shellization_summary.diagnostics.public_default_behavior_changed
    @test !shellization_summary.diagnostics.hamiltonian_schema_changed
    @test !shellization_summary.diagnostics.gto_supplement_semantics_changed
    @test !shellization_summary.diagnostics.raw_or_diagnostic_weights_promoted

    @test fixture_shellization_summary.object_kind == :cartesian_shellization_route_summary
    @test fixture_shellization_summary.route_family == :white_lindsey_low_order
    @test fixture_shellization_summary.source_kind == :white_lindsey_one_center_seed
    @test fixture_shellization_summary.shellization_role ==
          :seed_one_center_full_parent_shellization
    @test fixture_shellization_summary.shellization_stage == :route_neutral_spatial_planning
    @test fixture_shellization_summary.lowering_stage ==
          :not_lowered_by_shellization_summary
    @test fixture_shellization_summary.retained_dimension == total_retained_dimension

    @test length(shell.faces) == 6
    @test length(shell.edges) == 12
    @test length(shell.corners) == 8
    @test all(face -> face.side_first.retained_count == 3, shell.faces)
    @test all(face -> face.side_second.retained_count == 3, shell.faces)
    @test all(edge -> edge.side.retained_count == 3, shell.edges)

    @test length(shell.support_indices) == 7^3 - 5^3
    @test isempty(intersect(sequence.core_indices, shell.support_indices))
    @test sequence.working_box == (1:7, 1:7, 1:7)
    @test shell.provenance.source_box == (1:7, 1:7, 1:7)
    @test shell.provenance.next_inner_box == (2:6, 2:6, 2:6)
    @test shell.provenance.source_point_count == 7^3 - 5^3
    @test shell.provenance.retained_fixed_count == shell_retained_count

    @test sum(length, shell.face_column_ranges) == 54
    @test sum(length, shell.edge_column_ranges) == 36
    @test sum(length, shell.corner_column_ranges) == 8
    @test length(sequence.core_column_range) == core_retained_count
    @test length(only(sequence.layer_column_ranges)) == shell_retained_count
    @test first(sequence.core_column_range) == 1
    @test last(sequence.core_column_range) == core_retained_count
    @test only(sequence.layer_column_ranges) ==
          (core_retained_count + 1):total_retained_dimension

    @test fixed_block isa GaussletBases._NestedFixedBlock3D
    @test fixed_block.shell === sequence
    @test size(fixed_block.coefficient_matrix, 2) == total_retained_dimension
    @test size(fixed_block.overlap) == (total_retained_dimension, total_retained_dimension)
    @test norm(fixed_block.overlap - I, Inf) < 1.0e-10
    # These are retained-basis integral weights carried by the nested fixed block.
    @test length(fixed_block.weights) == total_retained_dimension
    @test all(isfinite, fixed_block.weights)
    @test minimum(fixed_block.weights) > 0.0

    @test structure.parent_side_count == 7
    @test structure.working_box_side_count == 7
    @test structure.nside == 5
    @test structure.core_side_count == 5
    @test structure.shell_layer_count == 1
    @test structure.expected_shell_increment == shell_retained_count
    @test structure.expected_face_retained_count == 54
    @test structure.expected_edge_retained_count == 36
    @test structure.expected_corner_retained_count == 8
    @test structure.total_face_retained_count == 54
    @test structure.total_edge_retained_count == 36
    @test structure.total_corner_retained_count == 8
    @test structure.total_expected_gausslet_count == total_retained_dimension
    @test structure.total_actual_gausslet_count == total_retained_dimension
    @test structure.layers_match_expected

    @test inventory.object_kind == :white_lindsey_low_order_materialized_seed_inventory
    @test inventory.route_family == :white_lindsey_low_order
    @test inventory.status == :private_development_seed
    @test inventory.private_development_only
    @test inventory.packet_kernel == :factorized_direct
    @test inventory.parent_side_count == 7
    @test inventory.source_side_count == 7
    @test inventory.nside == 5
    @test inventory.piece_counts == (core = 1, faces = 6, edges = 12, corners = 8)
    @test inventory.support_counts == (core = 125, shell = 218, total_source = 343)
    @test inventory.retained_counts == (
        core = 125,
        faces = 54,
        edges = 36,
        corners = 8,
        shell = 98,
        total = 223,
    )
    @test inventory.retained_ranges.core == 1:125
    @test inventory.retained_ranges.shell == 126:223
    @test inventory.retained_ranges.faces ==
          Tuple((first(range) + 125):(last(range) + 125) for range in shell.face_column_ranges)
    @test inventory.retained_ranges.edges ==
          Tuple((first(range) + 125):(last(range) + 125) for range in shell.edge_column_ranges)
    @test inventory.retained_ranges.corners ==
          Tuple((first(range) + 125):(last(range) + 125) for range in shell.corner_column_ranges)
    @test inventory.materialized_shell_local_ranges.faces == Tuple(shell.face_column_ranges)
    @test inventory.materialized_shell_local_ranges.edges == Tuple(shell.edge_column_ranges)
    @test inventory.materialized_shell_local_ranges.corners == Tuple(shell.corner_column_ranges)
    @test inventory.fixed_block_ready
    @test inventory.overlap_ready
    @test inventory.retained_basis_integral_weights_ready
    @test inventory.weight_semantics == :retained_basis_integral_weights

    @test route_units.object_kind == :white_lindsey_low_order_materialized_seed_route_units
    @test route_units.route_family == :white_lindsey_low_order
    @test route_units.status == :private_development_seed
    @test route_units.packet_kernel == :factorized_direct
    @test route_units.retained_dimension == 223
    @test !route_units.operator_pairs_materialized
    @test isempty(route_units.pair_entries)
    @test route_units.pair_family_counts == (white_lindsey_low_order = 0,)
    @test route_units.weight_semantics == :retained_basis_integral_weights

    @test Tuple(unit.unit_key for unit in route_units.retained_units) == (
        :low_order_core_direct,
        :low_order_face_interiors,
        :low_order_edges,
        :low_order_corners,
    )
    @test Tuple(unit.retained_unit_kind for unit in route_units.retained_units) == (
        :white_lindsey_direct_core,
        :white_lindsey_face_interior_2d_products,
        :white_lindsey_edge_1d_side_functions,
        :white_lindsey_corner_direct_single_sites,
    )
    @test Tuple(unit.retained_rule_kind for unit in route_units.retained_units) == (
        :direct_parent_sites,
        :face_interior_2d_products_of_1d_retained_side_functions,
        :edge_1d_retained_side_functions,
        :corner_direct_single_site_pieces,
    )
    @test all(
        unit -> unit.weight_semantics == :retained_basis_integral_weights,
        route_units.retained_units,
    )
    @test route_units.unit_inventory.retained_counts == (
        low_order_core_direct = 125,
        low_order_face_interiors = 54,
        low_order_edges = 36,
        low_order_corners = 8,
    )
    @test route_units.unit_inventory.ranges == (
        low_order_core_direct = 1:125,
        low_order_face_interiors = 126:179,
        low_order_edges = 180:215,
        low_order_corners = 216:223,
    )
    @test route_units.standard_unit_inventory.retained_counts_materialized
    @test route_units.standard_unit_inventory.retained_ranges_materialized
    @test route_units.standard_unit_inventory.retained_dimension == 223
    @test route_units.standard_unit_inventory.pair_count == 0
    @test route_units.standard_unit_inventory.pair_family_counts ==
          (white_lindsey_low_order = 0,)

    @test operator_inventory.object_kind ==
          :white_lindsey_low_order_materialized_seed_operator_inventory
    @test operator_inventory.route_family == :white_lindsey_low_order
    @test operator_inventory.status == :private_development_seed
    @test operator_inventory.packet_kernel == :factorized_direct
    @test operator_inventory.operator_source == :nested_fixed_block
    @test operator_inventory.retained_dimension == 223
    @test operator_inventory.terms == (
        :overlap,
        :position_x,
        :position_y,
        :position_z,
        :x2_x,
        :x2_y,
        :x2_z,
        :kinetic,
    )
    @test all(==((223, 223)), values(operator_inventory.matrix_sizes))
    @test all(values(operator_inventory.finite_ready))
    @test operator_inventory.all_finite
    @test operator_inventory.overlap_identity_error < 1.0e-10
    @test operator_inventory.overlap_identity_ready
    @test operator_inventory.symmetric_ready.overlap
    @test operator_inventory.symmetric_ready.x2_x
    @test operator_inventory.symmetric_ready.x2_y
    @test operator_inventory.symmetric_ready.x2_z
    @test operator_inventory.symmetric_ready.kinetic
    @test !operator_inventory.operator_pairs_materialized
    @test !operator_inventory.electron_electron_materialized

    report = GaussletBases._white_lindsey_low_order_materialized_seed_report()
    @test report.object_kind == :white_lindsey_low_order_materialized_seed_report
    @test report.route_family == :white_lindsey_low_order
    @test report.status == :private_development_seed
    @test report.private_development_only
    @test report.packet_kernel == :factorized_direct
    @test report.fixture.packet_kernel == :factorized_direct
    @test report.inventory === report.fixture.inventory
    @test report.route_units.retained_dimension == route_units.retained_dimension
    @test report.retained_dimension == route_units.retained_dimension
    @test report.operator_inventory.retained_dimension == operator_inventory.retained_dimension
    @test report.shellization_summary === report.fixture.shellization_summary
    @test report.shellization_summary_available
    @test report.shellization_source == :white_lindsey_one_center_seed
    @test !report.route_configured_shellization_consumed
    @test report.materialized_shellization_stage == :route_neutral_spatial_planning
    @test report.seed_materialization_status == :seed_based_private_materialization
    @test Tuple(unit.unit_key for unit in report.route_units.retained_units) ==
          Tuple(unit.unit_key for unit in route_units.retained_units)
    @test report.operator_inventory.terms == operator_inventory.terms
    @test !report.operator_pairs_materialized
    @test !report.electron_electron_materialized
    @test report.weight_semantics == :retained_basis_integral_weights

    direct_inventory = GaussletBases._white_lindsey_low_order_materialized_seed_inventory(
        sequence,
        fixed_block,
        structure,
    )
    direct_route_units =
        GaussletBases._white_lindsey_low_order_materialized_seed_route_units(direct_inventory)
    direct_operator_inventory =
        GaussletBases._white_lindsey_low_order_materialized_seed_operator_inventory(fixed_block)
    @test direct_inventory.packet_kernel === nothing
    @test direct_route_units.packet_kernel === nothing
    @test direct_operator_inventory.packet_kernel === nothing

    one_center_report = (
        route_family = :white_lindsey_low_order,
        system_metadata = (
            atom_symbols = ("Be",),
            nuclear_charges = (4,),
            atom_locations = ((0.0, 0.0, 0.0),),
            parent_axis_counts = (x = 7, y = 7, z = 7),
            parent_axis_counts_source = :manual_fixture,
            parent_box = (x = -3.0:3.0, y = -3.0:3.0, z = -3.0:3.0),
        ),
        recipe_metadata = (
            route_kind = :one_center_low_order_probe,
            route_shape = (:low_order_units,),
        ),
    )
    one_center_request =
        GaussletBases._cartesian_shellization_route_configured_request(one_center_report)
    one_center_plan =
        GaussletBases._cartesian_shellization_route_planning_stub(one_center_request)
    one_center_helper_map =
        GaussletBases._cartesian_shellization_route_planning_helper_map(one_center_plan)
    one_center_readiness =
        GaussletBases._cartesian_shellization_route_materializer_input_readiness(
            one_center_request,
            one_center_plan,
            one_center_helper_map,
        )
    one_center_config =
        GaussletBases._cartesian_shellization_route_materializer_config(
            one_center_request,
            one_center_plan,
            one_center_helper_map,
            one_center_readiness,
        )
    one_center_materialization =
        GaussletBases._cartesian_shellization_route_materialize_one_center_low_order(
            one_center_config;
            expansion = density_expansion,
        )

    @test one_center_materialization.object_kind ==
          :cartesian_shellization_route_one_center_materialization
    @test one_center_materialization.status ==
          :materialized_route_configured_one_center_low_order
    @test one_center_materialization.private_development_only
    @test one_center_materialization.route_family == :white_lindsey_low_order
    @test one_center_materialization.route_kind == :one_center_low_order_probe
    @test one_center_materialization.planning_family == :one_center_atomic_shellization
    @test one_center_materialization.route_configured_shellization_consumed
    @test one_center_materialization.calls_white_lindsey_seed_fixture
    @test !one_center_materialization.calls_lower_level_one_center_helpers_directly
    @test !one_center_materialization.public_default_behavior_changed
    @test one_center_materialization.materializer_options.parent_side_count == 7
    @test one_center_materialization.materializer_options.nside == 5
    @test one_center_materialization.materializer_options.Z == 4.0
    @test one_center_materialization.fixture.parent_side_count == 7
    @test one_center_materialization.fixture.nside == 5
    @test one_center_materialization.sequence === one_center_materialization.fixture.sequence
    @test one_center_materialization.fixed_block ===
          one_center_materialization.fixture.fixed_block
    @test one_center_materialization.retained_dimension == total_retained_dimension
    @test one_center_materialization.shellization_summary.source_kind ==
          :route_configured_one_center_low_order
    @test one_center_materialization.shellization_summary.shellization_role ==
          :route_configured_one_center_full_parent_shellization
    @test one_center_materialization.shellization_summary.retained_dimension ==
          total_retained_dimension

    support_fixture =
        GaussletBases._white_lindsey_low_order_materialized_seed_fixture(
            expansion = density_expansion,
            packet_kernel = :support_reference,
        )
    support_inventory = support_fixture.inventory
    support_route_units =
        GaussletBases._white_lindsey_low_order_materialized_seed_route_units(support_fixture)
    support_operator_inventory =
        GaussletBases._white_lindsey_low_order_materialized_seed_operator_inventory(
            support_fixture,
        )

    @test support_fixture.packet_kernel == :support_reference
    @test support_inventory.packet_kernel == :support_reference
    @test support_route_units.packet_kernel == :support_reference
    @test support_operator_inventory.packet_kernel == :support_reference
    @test support_inventory.retained_counts == inventory.retained_counts
    @test support_inventory.retained_ranges == inventory.retained_ranges
    @test support_inventory.retained_basis_integral_weights_ready
    @test support_route_units.retained_dimension == route_units.retained_dimension
    @test support_route_units.unit_inventory.retained_counts ==
          route_units.unit_inventory.retained_counts
    @test support_route_units.unit_inventory.ranges == route_units.unit_inventory.ranges
    @test Tuple(unit.unit_key for unit in support_route_units.retained_units) ==
          Tuple(unit.unit_key for unit in route_units.retained_units)
    @test Tuple(unit.retained_unit_kind for unit in support_route_units.retained_units) ==
          Tuple(unit.retained_unit_kind for unit in route_units.retained_units)
    @test support_operator_inventory.terms == operator_inventory.terms

    support_report =
        GaussletBases._white_lindsey_low_order_materialized_seed_report(
            packet_kernel = :support_reference,
        )
    @test support_report.packet_kernel == :support_reference
    @test support_report.fixture.packet_kernel == :support_reference
    @test support_report.inventory.packet_kernel == :support_reference
    @test support_report.route_units.packet_kernel == :support_reference
    @test support_report.operator_inventory.packet_kernel == :support_reference
    @test support_report.retained_dimension == support_route_units.retained_dimension
    @test Tuple(unit.unit_key for unit in support_report.route_units.retained_units) ==
          Tuple(unit.unit_key for unit in support_route_units.retained_units)
    @test support_report.operator_inventory.terms == support_operator_inventory.terms
    @test !support_report.operator_pairs_materialized
    @test !support_report.electron_electron_materialized

    @test !isnothing(fixed_block.pair_sum)
    @test !isnothing(support_fixture.fixed_block.pair_sum)
    direct_interaction =
        GaussletBases._qwrg_fixed_block_interaction_matrix(fixed_block, density_expansion)
    support_interaction = GaussletBases._qwrg_fixed_block_interaction_matrix(
        support_fixture.fixed_block,
        density_expansion,
    )
    @test size(direct_interaction) == (total_retained_dimension, total_retained_dimension)
    @test size(support_interaction) == size(direct_interaction)
    @test size(fixed_block.coefficient_matrix, 2) == size(direct_interaction, 1)
    @test size(support_fixture.fixed_block.coefficient_matrix, 2) ==
          size(support_interaction, 1)
    @test length(fixed_block.weights) == size(direct_interaction, 1)
    @test length(support_fixture.fixed_block.weights) == size(support_interaction, 1)
    @test all(isfinite, fixed_block.weights)
    @test all(isfinite, support_fixture.fixed_block.weights)
    @test minimum(fixed_block.weights) > 0.0
    @test minimum(support_fixture.fixed_block.weights) > 0.0
    @test all(isfinite, direct_interaction)
    @test all(isfinite, support_interaction)
    direct_interaction_symmetry_error = norm(
        direct_interaction - transpose(direct_interaction),
        Inf,
    )
    support_interaction_symmetry_error = norm(
        support_interaction - transpose(support_interaction),
        Inf,
    )
    @test direct_interaction_symmetry_error <= 1.0e-12
    @test support_interaction_symmetry_error <= 1.0e-12
    density_interaction_comparison_error = norm(
        direct_interaction - support_interaction,
        Inf,
    )
    @test density_interaction_comparison_error <= 1.0e-10

    ham_candidate =
        GaussletBases._white_lindsey_low_order_materialized_seed_ham_payload_candidate(
            fixture;
            expansion = density_expansion,
            Z = 2.0,
        )
    expected_one_body = GaussletBases._qwrg_fixed_block_one_body_matrix(
        fixed_block,
        density_expansion;
        Z = 2.0,
    )
    @test ham_candidate.object_kind == :white_lindsey_low_order_ham_payload_candidate
    @test ham_candidate.route_family == :white_lindsey_low_order
    @test ham_candidate.status == :private_payload_candidate_not_writer_adapted
    @test ham_candidate.private_development_only
    @test !ham_candidate.public_api
    @test !ham_candidate.writer_ready
    @test !ham_candidate.export_ready
    @test ham_candidate.export_status == :private_payload_candidate_not_writer_adapted
    @test ham_candidate.ham_bundle_export_status ==
          :blocked_private_payload_candidate_not_writer_adapted
    @test ham_candidate.packet_kernel == :factorized_direct
    @test ham_candidate.retained_dimension == total_retained_dimension
    @test ham_candidate.weight_semantics == :retained_basis_integral_weights
    @test ham_candidate.one_body_source == :fixed_block_kinetic_minus_Z_gaussian_sum
    @test ham_candidate.interaction_source == :fixed_block_pair_sum
    @test ham_candidate.expansion_term_count == length(density_expansion.exponents)
    @test ham_candidate.nuclear_charge == 2.0
    @test ham_candidate.overlap == Matrix{Float64}(fixed_block.overlap)
    @test ham_candidate.one_body_hamiltonian == expected_one_body
    @test ham_candidate.interaction_matrix == direct_interaction
    @test ham_candidate.final_integral_weights == Float64[Float64(w) for w in fixed_block.weights]
    @test size(ham_candidate.overlap) == (total_retained_dimension, total_retained_dimension)
    @test size(ham_candidate.one_body_hamiltonian) ==
          (total_retained_dimension, total_retained_dimension)
    @test size(ham_candidate.interaction_matrix) ==
          (total_retained_dimension, total_retained_dimension)
    @test length(ham_candidate.final_integral_weights) == total_retained_dimension
    @test ham_candidate.checks.matrix_sizes == (
        overlap = (223, 223),
        one_body_hamiltonian = (223, 223),
        interaction_matrix = (223, 223),
    )
    @test ham_candidate.checks.expected_matrix_size == (223, 223)
    @test ham_candidate.checks.matrix_size_ready
    @test all(values(ham_candidate.checks.finite_ready))
    @test ham_candidate.checks.all_finite
    @test all(values(ham_candidate.checks.symmetric_ready))
    @test ham_candidate.checks.all_symmetric
    @test ham_candidate.checks.symmetry_errors.overlap <= 1.0e-12
    @test ham_candidate.checks.symmetry_errors.one_body_hamiltonian <= 1.0e-12
    @test ham_candidate.checks.symmetry_errors.interaction_matrix <= 1.0e-12
    @test ham_candidate.checks.weight_length_ready
    @test ham_candidate.checks.weights_finite
    @test ham_candidate.checks.weights_ready
    @test ham_candidate.checks.gaussian_sum_available
    @test ham_candidate.checks.pair_sum_available
    @test !ham_candidate.checks.writer_ready
    @test !ham_candidate.checks.export_ready
    @test ham_candidate.checks.export_status ==
          :private_payload_candidate_not_writer_adapted

    fixed_block_candidate =
        GaussletBases._white_lindsey_low_order_materialized_seed_ham_payload_candidate(
            fixed_block;
            expansion = density_expansion,
            Z = 2.0,
        )
    @test fixed_block_candidate.packet_kernel === nothing
    @test fixed_block_candidate.one_body_hamiltonian == expected_one_body
    @test fixed_block_candidate.interaction_matrix == direct_interaction
    @test !fixed_block_candidate.export_ready

    ham_adapter =
        GaussletBases._white_lindsey_low_order_materialized_seed_ham_bundle_adapter(
            fixture;
            expansion = density_expansion,
            Z = 2.0,
        )
    @test ham_adapter isa GaussletBases._WhiteLindseyLowOrderHamBundleAdapter
    @test ham_adapter.fixed_block === fixed_block
    @test ham_adapter.candidate == ham_candidate
    @test ham_adapter.expansion === density_expansion

    adapter_bundle = cartesian_basis_bundle_payload(ham_adapter)
    @test adapter_bundle.basis["basis_kind"] == "nested_fixed_block"
    @test adapter_bundle.ham !== nothing
    @test adapter_bundle.ham["format"] == "cartesian_hamiltonian_bundle_v1"
    @test adapter_bundle.ham["model_kind"] == "white_lindsey_low_order"
    @test adapter_bundle.ham["route_family"] == "white_lindsey_low_order"
    @test adapter_bundle.ham["interaction_model"] == "density_density"
    @test adapter_bundle.ham["interaction_treatment"] == "nested_fixed_block_pair_sum"
    @test adapter_bundle.ham["payload_candidate_status"] ==
          "private_payload_candidate_not_writer_adapted"
    @test size(adapter_bundle.ham["overlap"]) ==
          (total_retained_dimension, total_retained_dimension)
    @test size(adapter_bundle.ham["one_body_hamiltonian"]) ==
          (total_retained_dimension, total_retained_dimension)
    @test size(adapter_bundle.ham["interaction_matrix"]) ==
          (total_retained_dimension, total_retained_dimension)
    @test adapter_bundle.ham["basis_integral_weights"] ==
          adapter_bundle.basis["final_integral_weights"]
    @test adapter_bundle.ham["basis_integral_weights"] == ham_candidate.final_integral_weights
    @test adapter_bundle.ham["nuclear_charge"] == 2.0
    @test adapter_bundle.meta["has_ham"]

    mktempdir() do dir
        hamfile = joinpath(dir, "white_lindsey_adapter_ham.jld2")
        @test write_cartesian_basis_bundle_jld2(hamfile, ham_adapter) == hamfile
        @test isfile(hamfile)
        jldopen(hamfile, "r") do file
            top_keys = Set(
                key isa AbstractVector ? join(string.(key), "/") : string(key) for key in keys(file)
            )
            @test "basis" in top_keys
            @test "ham" in top_keys
            @test "meta" in top_keys
            @test String(file["ham/format"]) == "cartesian_hamiltonian_bundle_v1"
            @test String(file["ham/model_kind"]) == "white_lindsey_low_order"
            @test String(file["ham/route_family"]) == "white_lindsey_low_order"
            @test String(file["ham/interaction_model"]) == "density_density"
            @test String(file["ham/interaction_treatment"]) == "nested_fixed_block_pair_sum"
            @test size(file["ham/overlap"]) ==
                  (total_retained_dimension, total_retained_dimension)
            @test size(file["ham/one_body_hamiltonian"]) ==
                  (total_retained_dimension, total_retained_dimension)
            @test size(file["ham/interaction_matrix"]) ==
                  (total_retained_dimension, total_retained_dimension)
            @test file["ham/basis_integral_weights"] == file["basis/final_integral_weights"]
            @test length(file["ham/basis_integral_weights"]) == total_retained_dimension
            @test Bool(file["meta/has_ham"])
        end
    end

    direct_matrices =
        GaussletBases._white_lindsey_low_order_materialized_seed_operator_matrices(fixed_block)
    support_matrices =
        GaussletBases._white_lindsey_low_order_materialized_seed_operator_matrices(
            support_fixture.fixed_block,
        )
    matrix_comparison_tolerance = 1.0e-10
    for term in operator_inventory.terms
        @test norm(
            getproperty(direct_matrices, term) - getproperty(support_matrices, term),
            Inf,
        ) <= matrix_comparison_tolerance
    end
end

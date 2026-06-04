using Test
using GaussletBases
using JLD2

const _ROUTE_DRIVER_STANDARD_UNIT_INVENTORY_KEYS = (
    :unit_count,
    :unit_keys,
    :retained_unit_kinds,
    :source_families,
    :source_dimensions,
    :retained_counts,
    :retained_dimension,
    :retained_counts_materialized,
    :retained_ranges_materialized,
    :pair_count,
    :pair_family_counts,
    :pair_families,
    :output_representations,
)
const _ROUTE_DRIVER_PARENT_CONTRACT_KEYS = (
    :object_kind,
    :status,
    :atom_count,
    :center_count,
    :atom_symbols,
    :nuclear_charges,
    :atom_locations,
    :center_table,
    :center_axis_metadata,
    :system_classification,
    :system_classification_status,
    :bond_axis,
    :chain_axis,
    :parent_axis_counts,
    :parent_axis_counts_source,
    :parent_axis_counts_status,
    :parent_box,
    :parent_box_rule,
    :parent_materialization_plan,
    :parent_materialization_plan_status,
    :parent_materialization_planning_family,
    :parent_materialization_blocker,
    :parent_basis_object_available,
    :parent_qw_basis_object_available,
    :parent_axis_bundle_object_available,
    :parent_basis_object_type_label,
    :parent_qw_basis_object_type_label,
    :parent_axis_bundle_object_type_label,
    :parent_basis_materialization_status,
    :parent_basis_materialization,
    :parent_basis_materialized,
    :parent_axis_metadata_constructed,
    :axis_bundle_materialized,
    :diagnostics,
)
const _PQS_ROUTE_DRIVER_TEST_ROUTE_KIND = :be2_cartesian_nesting_route_driver_spine

function _pqs_route_driver_report_for_test(;
    route_family = :pqs_source_box,
    route_kind = _PQS_ROUTE_DRIVER_TEST_ROUTE_KIND,
    atom_symbols = ("Be", "Be"),
    nuclear_charges = (4, 4),
    atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
    parent_axis_counts = (x = 9, y = 7, z = 9),
)
    return GaussletBases._pqs_source_box_route_driver_dry_run(
        ;
        route_family,
        route_kind,
        atom_symbols,
        nuclear_charges,
        atom_locations,
        radius = 15.0,
        parent_axis_counts,
        map_backend = :pgdg_localized_experimental,
        q = 5,
        n_s = 5,
        reference_spacing = 1.0,
        tail_spacing = 10.0,
        q_to_core_spacing_rule = :standard_pqs_ns_equals_q,
        core_spacing = nothing,
        probe_parent_axis_construction = false,
        parent_axis_probe_backend = :pgdg_localized_experimental,
        parent_axis_probe_family = :G10,
        probe_raw_product_box_plans = false,
        raw_product_box_probe_backend = :pgdg_localized_experimental,
        route_shape = (:pqs_left, :product, :pqs_right),
        product_body_rule = :centered_single_z_slab,
        pqs_retained_rule = :boundary_comx_product_mode_selection,
        product_retained_rule = :product_doside_retained_unit,
        terms = (
            :overlap,
            :position_x,
            :position_y,
            :position_z,
            :x2_x,
            :x2_y,
            :x2_z,
            :kinetic,
        ),
        pair_factor_normalization = :density_normalized,
        support_dense_direct_allowed = false,
        reference_only_authorities = (:support_row_oracle, :dense_parent_projection),
    )
end

function _pqs_route_driver_one_center_report_for_test()
    return _pqs_route_driver_report_for_test(
        ;
        route_family = :white_lindsey_low_order,
        route_kind = :one_center_parent_contract_report_probe,
        atom_symbols = ("Be",),
        nuclear_charges = (4,),
        atom_locations = ((0.0, 0.0, 0.0),),
        parent_axis_counts = (x = 7, y = 7, z = 7),
    )
end

function _pqs_route_driver_one_center_white_lindsey_report_for_test()
    return (;
        route_family = :white_lindsey_low_order,
        retained_dimension = 223,
        system_metadata = (;
            atom_symbols = ("Be",),
            nuclear_charges = (4,),
            atom_locations = ((0.0, 0.0, 0.0),),
            parent_axis_counts = (x = 7, y = 7, z = 7),
            parent_axis_counts_source = :manual_fixture,
            parent_box = (x = -3.0:3.0, y = -3.0:3.0, z = -3.0:3.0),
            map_backend = :pgdg_localized_experimental,
        ),
        recipe_metadata = (;
            route_kind = :one_center_low_order_probe,
            route_shape = (:standard_cartesian_units, :low_order_comx_coarsening),
            benchmark_role = :published_cartesian_baseline_for_pqs_comparison,
            n_s = 5,
            core_spacing = 0.15,
            reference_spacing = 1.0,
            tail_spacing = 10.0,
            parent_axis_probe_backend = :pgdg_localized_experimental,
        ),
    )
end

function _pqs_route_driver_check_standard_unit_inventory(
    report;
    retained_dimension,
    unit_keys,
    pair_count,
    pair_family_counts,
    retained_counts_materialized,
    retained_ranges_materialized,
    pair_families,
    output_representations,
)
    @test report.retained_dimension == retained_dimension
    @test report.recipe_metadata.route_kind == _PQS_ROUTE_DRIVER_TEST_ROUTE_KIND
    @test Tuple(unit.unit_key for unit in report.retained_units) == unit_keys
    @test length(report.pair_entries) == pair_count
    @test report.pair_family_counts == pair_family_counts

    @test hasproperty(report, :standard_unit_inventory)
    summary = report.standard_unit_inventory
    @test Tuple(keys(summary)) == _ROUTE_DRIVER_STANDARD_UNIT_INVENTORY_KEYS
    @test summary.unit_count == length(unit_keys)
    @test summary.unit_keys == unit_keys
    @test summary.retained_unit_kinds ==
          Tuple(unit.retained_unit_kind for unit in report.retained_units)
    @test summary.source_families ==
          Tuple(unit.source_family for unit in report.retained_units)
    @test summary.source_dimensions == report.source_dimensions
    @test summary.retained_counts == report.retained_counts
    @test summary.retained_dimension == retained_dimension
    @test summary.retained_counts_materialized == retained_counts_materialized
    @test summary.retained_ranges_materialized == retained_ranges_materialized
    @test summary.pair_count == pair_count
    @test summary.pair_family_counts == pair_family_counts
    @test summary.pair_families == pair_families
    @test summary.output_representations == output_representations
    return nothing
end

function _pqs_route_driver_check_parent_contract(
    report;
    atom_count,
    system_classification,
    system_classification_status,
    bond_axis,
    chain_axis,
    parent_materialization_planning_family,
    constructs_basis_now = false,
    constructs_axis_bundle_now = false,
    parent_basis_materialization_status = :metadata_only_not_materialized,
)
    @test hasproperty(report, :parent_contract)
    parent_contract = report.parent_contract
    @test Tuple(keys(parent_contract)) == _ROUTE_DRIVER_PARENT_CONTRACT_KEYS
    @test parent_contract.object_kind == :cartesian_route_parent_contract
    @test parent_contract.atom_count == atom_count
    @test parent_contract.center_count == atom_count
    @test parent_contract.system_classification == system_classification
    @test parent_contract.system_classification_status == system_classification_status
    @test parent_contract.bond_axis == bond_axis
    @test parent_contract.chain_axis == chain_axis
    @test parent_contract.parent_materialization_plan.object_kind ==
          :cartesian_parent_materialization_plan
    @test parent_contract.parent_materialization_planning_family ==
          parent_materialization_planning_family
    @test parent_contract.parent_materialization_plan_status ==
          parent_contract.parent_materialization_plan.status
    @test parent_contract.parent_materialization_blocker ==
          parent_contract.parent_materialization_plan.blocker
    @test !hasproperty(parent_contract, :parent_basis_object)
    @test !hasproperty(parent_contract, :parent_qw_basis_object)
    @test !hasproperty(parent_contract, :parent_axis_bundle_object)
    @test parent_contract.parent_basis_object_available ==
          parent_contract.parent_materialization_plan.constructs_basis_now
    @test parent_contract.parent_axis_bundle_object_available ==
          parent_contract.parent_materialization_plan.constructs_axis_bundle_now
    @test parent_contract.parent_materialization_plan.constructs_basis_now ==
          constructs_basis_now
    @test parent_contract.parent_materialization_plan.constructs_axis_bundle_now ==
          constructs_axis_bundle_now
    @test parent_contract.parent_axis_counts == report.system_metadata.parent_axis_counts
    @test parent_contract.parent_axis_counts_source ==
          report.system_metadata.parent_axis_counts_source
    @test parent_contract.parent_box == report.system_metadata.parent_box
    @test parent_contract.parent_basis_materialization_status ==
          parent_basis_materialization_status
    @test parent_contract.parent_basis_materialized == constructs_basis_now
    @test parent_contract.axis_bundle_materialized == constructs_axis_bundle_now
    @test parent_contract.diagnostics.parent_contract_driven_downstream_metadata
    @test parent_contract.diagnostics.public_default_behavior_changed == false

    @test report.parent_description.parent_contract_status == parent_contract.status
    @test report.parent_description.system_classification == system_classification
    @test report.parent_description.system_classification_status ==
          system_classification_status
    @test report.parent_description.bond_axis == bond_axis
    @test report.parent_description.chain_axis == chain_axis
    @test report.parent_description.parent_materialization_plan_status ==
          parent_contract.parent_materialization_plan_status
    @test report.parent_description.parent_materialization_planning_family ==
          parent_materialization_planning_family
    @test report.parent_description.parent_materialization_blocker ==
          parent_contract.parent_materialization_blocker
    @test report.parent_description.parent_basis_materialization_status ==
          parent_contract.parent_basis_materialization_status
    return nothing
end

function _pqs_route_driver_check_report_output_sections(report)
    tmpdir = mktempdir()
    try
        textfile = joinpath(tmpdir, "report.txt")
        open(textfile, "w") do io
            redirect_stdout(io) do
                GaussletBases._pqs_source_box_route_driver_print_details(report)
            end
        end
        text = read(textfile, String)
        @test occursin("[standard_unit_inventory]", text)
        @test occursin("[parent_contract]", text)
        @test occursin("[retained_units]", text)
        @test occursin("[pair_inventory]", text)

        jld2file = joinpath(tmpdir, "report.jld2")
        open(joinpath(tmpdir, "jld2_stdout.txt"), "w") do io
            redirect_stdout(io) do
                GaussletBases._pqs_source_box_route_driver_save(
                    report;
                    save_artifact = true,
                    save_tsv = false,
                    outfile = jld2file,
                    tsvfile = joinpath(tmpdir, "unused.tsv"),
                )
            end
        end
        jldopen(jld2file, "r") do file
            top_keys = Set(
                key isa AbstractVector ? join(string.(key), "/") : string(key) for key in keys(file)
            )
            @test "report" in top_keys
            @test !("materialization" in top_keys)
        end

        tsvfile = joinpath(tmpdir, "report.tsv")
        stdout_file = joinpath(tmpdir, "save_stdout.txt")
        open(stdout_file, "w") do io
            redirect_stdout(io) do
                GaussletBases._pqs_source_box_route_driver_save(
                    report;
                    save_artifact = false,
                    save_tsv = true,
                    outfile = joinpath(tmpdir, "unused.jld2"),
                    tsvfile,
                )
            end
        end
        tsv = read(tsvfile, String)
        @test occursin("parent_contract\tsystem_classification", tsv)
        @test occursin("parent_contract\tparent_materialization_plan_status", tsv)
        @test occursin("parent_contract\tparent_basis_materialization_status", tsv)
        @test occursin("standard_unit_inventory\tunit_count", tsv)
        @test occursin("retained_unit\t", tsv)
        @test occursin("pair_entry\t", tsv)
    finally
        rm(tmpdir; recursive = true, force = true)
    end
    @test !isfile("cartesian_ham_builder_report.jld2")
    @test !isfile("cartesian_ham_builder_report.tsv")
    return nothing
end

function _pqs_route_driver_check_be2_shellization_request(
    materialization;
    current_materialization_source,
)
    @test materialization.route_configured_shellization_request_available
    @test materialization.route_configured_shellization_request_status ==
          :metadata_only_pending_materializer
    @test materialization.route_configured_system_classification ==
          :bond_aligned_diatomic
    @test materialization.route_configured_system_classification_status ==
          :explicit_two_atom_single_axis_separation
    @test materialization.route_configured_bond_axis == :x
    @test materialization.route_configured_shellization_plan_available
    @test materialization.route_configured_shellization_plan_status ==
          :metadata_only_pending_materializer
    @test materialization.route_configured_shellization_planning_status ==
          :planned_metadata_only
    @test materialization.route_configured_shellization_planning_family ==
          :bond_aligned_diatomic_shellization
    @test materialization.route_configured_midpoint_slab_status ==
          :conditional_diatomic_midpoint_slab
    @test materialization.route_configured_shellization_helper_map_available
    @test materialization.route_configured_shellization_helper_map_status ==
          :metadata_only_pending_materializer_inputs
    @test materialization.route_configured_primary_planned_helper ==
          :_nested_bond_aligned_diatomic_source
    @test materialization.route_configured_missing_input_count == 8
    @test materialization.route_configured_helper_map_blocker ==
          :pending_route_configured_bond_aligned_diatomic_materializer_inputs
    @test materialization.route_configured_input_readiness_available
    @test materialization.route_configured_input_readiness_status ==
          :blocked_missing_materializer_inputs
    @test materialization.route_configured_available_fact_count == 10
    @test materialization.route_configured_materializer_missing_input_count == 8
    @test materialization.route_configured_input_readiness_blocker ==
          :pending_route_configured_bond_aligned_diatomic_materializer_inputs
    @test materialization.route_configured_materializer_config_available
    @test materialization.route_configured_materializer_config_status ==
          :blocked_missing_materializer_inputs
    @test materialization.route_configured_materializer_config_planning_family ==
          :bond_aligned_diatomic_shellization
    @test materialization.route_configured_materializer_config_pending_input_count == 8

    request = materialization.route_configured_shellization_request
    @test request.object_kind == :cartesian_shellization_route_configured_request
    @test request.status == materialization.route_configured_shellization_request_status
    @test request.private_development_only
    @test request.route_family == materialization.route_family
    @test request.route_kind == _PQS_ROUTE_DRIVER_TEST_ROUTE_KIND
    @test request.atom_count == 2
    @test request.atom_symbols == ("Be", "Be")
    @test request.system_classification == materialization.route_configured_system_classification
    @test request.system_classification_status ==
          materialization.route_configured_system_classification_status
    @test request.bond_axis == materialization.route_configured_bond_axis
    @test request.classification_source == :parent_contract
    @test request.parent_contract_available
    @test request.parent_contract_status == :available
    @test request.parent_basis_materialization_status ==
          :metadata_only_not_materialized
    @test request.requested_shellization_stage == :route_neutral_spatial_planning
    @test request.expected_next_materializer_status ==
          :pending_route_configured_shellization_materializer
    @test request.current_materialization_source == current_materialization_source
    @test !request.route_configured_shellization_consumed
    @test !request.constructs_basis
    @test !request.constructs_shell_sequence
    @test !request.constructs_fixed_block

    plan = materialization.route_configured_shellization_plan
    @test plan.object_kind == :cartesian_shellization_route_planning_stub
    @test plan.status == materialization.route_configured_shellization_plan_status
    @test plan.planning_status ==
          materialization.route_configured_shellization_planning_status
    @test plan.private_development_only
    @test plan.request_object_kind == request.object_kind
    @test plan.route_family == request.route_family
    @test plan.route_kind == request.route_kind
    @test plan.system_classification == request.system_classification
    @test plan.system_classification_status == request.system_classification_status
    @test plan.bond_axis == request.bond_axis
    @test plan.planning_family ==
          materialization.route_configured_shellization_planning_family
    @test plan.spatial_stage_order == (
        :atom_local_uncontracted_cores,
        :atom_local_shells,
        :contact_merge,
        :optional_midpoint_slab,
        :outer_rectangular_shell_boxes,
        :final_boundary_edge_adjustment,
    )
    @test plan.atom_local_core_status == :planned_two_atom_uncontracted_cores
    @test plan.atom_local_shell_status == :planned_two_atom_local_shells
    @test plan.contact_merge_status == :planned_contact_or_merge_decision
    @test plan.midpoint_slab_status ==
          materialization.route_configured_midpoint_slab_status
    @test plan.outer_rectangular_shell_status == :planned_after_contact_merge_region
    @test plan.boundary_edge_adjustment_status == :planned_final_adjustment
    @test plan.expected_next_materializer_status ==
          :pending_route_configured_shellization_materializer
    @test !plan.constructs_basis
    @test !plan.constructs_shell_sequence
    @test !plan.constructs_fixed_block
    @test !plan.route_configured_shellization_consumed

    helper_map = materialization.route_configured_shellization_helper_map
    @test helper_map.object_kind == :cartesian_shellization_route_planning_helper_map
    @test helper_map.status ==
          materialization.route_configured_shellization_helper_map_status
    @test helper_map.private_development_only
    @test helper_map.planning_family == plan.planning_family
    @test helper_map.primary_planned_helper ==
          materialization.route_configured_primary_planned_helper
    @test helper_map.helper_chain == (
        :_nested_bond_aligned_diatomic_source,
        :_nested_bond_aligned_diatomic_split_geometry,
        :_nested_bond_aligned_diatomic_sequence_for_box,
        :_nested_shell_sequence_from_core_block,
        :_nested_fixed_block,
    )
    @test helper_map.missing_input_count ==
          materialization.route_configured_missing_input_count
    @test helper_map.blocker == materialization.route_configured_helper_map_blocker
    @test !helper_map.route_configured_shellization_consumed
    @test !helper_map.calls_mapped_helpers

    readiness = materialization.route_configured_input_readiness
    @test readiness.object_kind ==
          :cartesian_shellization_route_materializer_input_readiness
    @test readiness.status == materialization.route_configured_input_readiness_status
    @test readiness.private_development_only
    @test readiness.route_family == request.route_family
    @test readiness.route_kind == request.route_kind
    @test readiness.planning_family == plan.planning_family
    @test readiness.available_fact_count ==
          materialization.route_configured_available_fact_count
    @test readiness.missing_inputs == helper_map.missing_inputs
    @test readiness.missing_input_count ==
          materialization.route_configured_materializer_missing_input_count
    @test readiness.blocker == materialization.route_configured_input_readiness_blocker
    @test !in(:nside, readiness.driver_defaults_not_materializer_contract)
    @test :refinement_options in readiness.driver_defaults_not_materializer_contract
    @test !readiness.materializer_ready
    @test !readiness.route_configured_shellization_consumed
    @test !readiness.constructs_basis
    @test !readiness.constructs_axis_bundles
    @test !readiness.constructs_shell_sequence
    @test !readiness.constructs_fixed_block

    config = materialization.route_configured_materializer_config
    @test config.object_kind == :cartesian_shellization_route_materializer_config
    @test config.status == materialization.route_configured_materializer_config_status
    @test config.private_development_only
    @test config.route_family == request.route_family
    @test config.route_kind == request.route_kind
    @test config.system_classification == request.system_classification
    @test config.bond_axis == request.bond_axis
    @test config.atom_symbols == request.atom_symbols
    @test config.atom_locations == request.atom_locations
    @test config.nuclear_charges == request.nuclear_charges
    @test config.parent_axis_counts == request.parent_axis_counts
    @test config.parent_box == request.parent_box
    @test config.planning_family ==
          materialization.route_configured_materializer_config_planning_family
    @test config.primary_planned_helper == helper_map.primary_planned_helper
    @test config.helper_chain == helper_map.helper_chain
    @test config.pending_inputs == readiness.missing_inputs
    @test config.pending_input_count ==
          materialization.route_configured_materializer_config_pending_input_count
    @test config.blocker == readiness.blocker
    @test config.driver_defaults_not_materializer_contract ==
          readiness.driver_defaults_not_materializer_contract
    @test !config.materializer_ready
    @test !config.route_configured_shellization_consumed
    @test !config.constructs_basis
    @test !config.constructs_axis_bundles
    @test !config.constructs_shell_sequence
    @test !config.constructs_fixed_block
    return nothing
end

function _pqs_route_driver_check_materialization_status(pqs_report, white_lindsey_report)
    density_expansion = coulomb_gaussian_expansion(doacc = false)
    default_status = GaussletBases._pqs_source_box_route_driver_materialization(
        pqs_report;
        materialize_route = false,
        save_basis_artifact = false,
        save_ham_artifact = false,
        basisfile = "unused_basis.jld2",
        hamfile = "unused.jld2",
    )
    @test default_status.status == :not_requested_metadata_only
    @test default_status.materialized_report === nothing
    @test default_status.shellization_summary === nothing
    @test !default_status.shellization_summary_available
    @test !default_status.route_configured_shellization_consumed
    @test default_status.materialized_shellization_stage == :not_checked_metadata_only
    _pqs_route_driver_check_be2_shellization_request(
        default_status;
        current_materialization_source = :pending_source_box_route_materializer,
    )
    @test !default_status.route_configured_one_center_materializer_probe_requested
    @test default_status.route_configured_one_center_materializer_probe_status ==
          :not_requested
    @test !default_status.route_configured_one_center_materializer_probe_materialized
    @test !default_status.route_configured_one_center_materializer_probe_consumed
    @test default_status.route_configured_one_center_materializer_probe_blocker === nothing
    @test default_status.route_configured_one_center_materializer_probe.materialization === nothing
    @test default_status.ham_bundle_export_status == :not_requested
    @test default_status.ham_preflight_status == :not_checked_metadata_only
    @test default_status.ham_missing_builder === nothing
    @test default_status.ham_operator_payload_status == :not_checked_metadata_only
    @test default_status.ham_interaction_status == :not_checked_metadata_only
    @test default_status.ham_preflight === nothing
    @test default_status.basis_artifact_status == :not_requested
    @test !default_status.basis_artifact_written
    @test !default_status.ham_artifact_written

    pqs_status = GaussletBases._pqs_source_box_route_driver_materialization(
        pqs_report;
        materialize_route = true,
        save_basis_artifact = false,
        save_ham_artifact = false,
        basisfile = "unused_basis.jld2",
        hamfile = "unused.jld2",
    )
    @test pqs_status.status == :pending_source_box_retained_route
    @test pqs_status.materialized_report === nothing
    @test pqs_status.shellization_summary === nothing
    @test !pqs_status.shellization_summary_available
    @test pqs_status.shellization_source == :pending_source_box_route_shellization
    @test !pqs_status.route_configured_shellization_consumed
    @test pqs_status.materialized_shellization_stage ==
          :pending_source_box_retained_route
    @test pqs_status.seed_materialization_status == :not_applicable
    _pqs_route_driver_check_be2_shellization_request(
        pqs_status;
        current_materialization_source = :pending_source_box_route_materializer,
    )
    @test pqs_status.final_integral_weights_status == :pending_final_ida_weights
    @test pqs_status.basis_bundle_export_status == :pending_final_retained_basis
    @test pqs_status.basis_artifact_status == :not_requested
    @test pqs_status.ham_preflight_status == :not_applicable_to_pqs_source_box_route
    @test pqs_status.ham_missing_builder == :pending_source_box_retained_route
    @test pqs_status.ham_operator_payload_status ==
          :pending_source_box_retained_operator_payload
    @test pqs_status.ham_interaction_status ==
          :pending_source_box_retained_density_density_blocks
    @test pqs_status.ham_bundle_export_status == :pending_source_box_retained_route
    @test pqs_status.ham_preflight === nothing
    @test pqs_status.pqs_materialization_status == :pending_source_box_retained_route

    white_lindsey_probe_blocked =
        GaussletBases._pqs_source_box_route_driver_materialization(
            white_lindsey_report;
            materialize_route = true,
            probe_route_configured_one_center_materializer = true,
            save_basis_artifact = false,
            save_ham_artifact = false,
            basisfile = "unused_basis.jld2",
            hamfile = "unused.jld2",
            white_lindsey_expansion = density_expansion,
        )
    @test white_lindsey_probe_blocked.status == :materialized_seed_report_available
    @test white_lindsey_probe_blocked.shellization_source ==
          :white_lindsey_one_center_seed
    @test !white_lindsey_probe_blocked.route_configured_shellization_consumed
    @test white_lindsey_probe_blocked.route_configured_one_center_materializer_probe_requested
    @test white_lindsey_probe_blocked.route_configured_one_center_materializer_probe_status ==
          :blocked_not_one_center
    @test !white_lindsey_probe_blocked.route_configured_one_center_materializer_probe_materialized
    @test !white_lindsey_probe_blocked.route_configured_one_center_materializer_probe_consumed
    @test white_lindsey_probe_blocked.route_configured_one_center_materializer_probe_blocker ==
          :route_config_not_one_center
    @test white_lindsey_probe_blocked.route_configured_one_center_materializer_probe.materialization ===
          nothing

    one_center_probe =
        GaussletBases._pqs_source_box_route_driver_materialization(
            _pqs_route_driver_one_center_white_lindsey_report_for_test();
            materialize_route = false,
            probe_route_configured_one_center_materializer = true,
            save_basis_artifact = false,
            save_ham_artifact = false,
            basisfile = "unused_basis.jld2",
            hamfile = "unused.jld2",
            white_lindsey_expansion = density_expansion,
        )
    @test one_center_probe.status == :not_requested_metadata_only
    @test one_center_probe.route_configured_system_classification == :one_center
    @test one_center_probe.route_configured_shellization_request.classification_source ==
          :system_metadata_fallback
    @test one_center_probe.route_configured_one_center_materializer_probe_requested
    @test one_center_probe.route_configured_one_center_materializer_probe_status ==
          :materialized_route_configured_one_center_low_order
    @test one_center_probe.route_configured_one_center_materializer_probe_materialized
    @test one_center_probe.route_configured_one_center_materializer_probe_consumed
    @test one_center_probe.route_configured_one_center_materializer_probe_blocker === nothing
    route_configured_materialization =
        one_center_probe.route_configured_one_center_materializer_probe.materialization
    @test route_configured_materialization.object_kind ==
          :cartesian_shellization_route_one_center_materialization
    @test route_configured_materialization.retained_dimension == 223
    @test route_configured_materialization.route_configured_shellization_consumed
    @test !route_configured_materialization.calls_white_lindsey_seed_fixture
    @test route_configured_materialization.calls_lower_level_one_center_helpers_directly
    @test route_configured_materialization.shellization_summary.source_kind ==
          :route_configured_one_center_low_order
    @test route_configured_materialization.materializer_options.gausslet_backend ==
          :pgdg_localized_experimental
    @test route_configured_materialization.materializer_options.d == 0.15
    @test route_configured_materialization.materializer_options.nside == 5
    @test one_center_probe.route_configured_materializer_backend_requested ==
          :pgdg_localized_experimental
    @test one_center_probe.route_configured_materializer_backend_consumed ==
          :pgdg_localized_experimental
    @test one_center_probe.route_configured_materializer_d_requested == 0.15
    @test one_center_probe.route_configured_materializer_d_consumed == 0.15
    @test one_center_probe.route_configured_materializer_nside_requested == 5
    @test one_center_probe.route_configured_materializer_nside_consumed == 5
    @test !route_configured_materialization.public_default_behavior_changed

    one_center_materialized =
        GaussletBases._pqs_source_box_route_driver_materialization(
            _pqs_route_driver_one_center_white_lindsey_report_for_test();
            materialize_route = true,
            save_basis_artifact = false,
            save_ham_artifact = false,
            basisfile = "unused_basis.jld2",
            hamfile = "unused.jld2",
            white_lindsey_expansion = density_expansion,
        )
    @test one_center_materialized.status ==
          :materialized_route_configured_one_center_report_available
    @test one_center_materialized.materialized_report.object_kind ==
          :white_lindsey_low_order_route_configured_one_center_report
    @test one_center_materialized.materialized_report_kind ==
          :white_lindsey_low_order_route_configured_one_center_report
    @test one_center_materialized.shellization_source ==
          :route_configured_one_center_low_order
    @test one_center_materialized.route_configured_shellization_consumed
    @test one_center_materialized.seed_materialization_status ==
          :not_seed_route_configured_materialization
    @test one_center_materialized.retained_dimension == 223
    @test one_center_materialized.route_configured_one_center_materializer_probe_requested
    @test one_center_materialized.route_configured_one_center_materializer_probe_materialized
    @test one_center_materialized.route_configured_one_center_materializer_probe_consumed
    @test one_center_materialized.route_configured_materializer_backend_requested ==
          :pgdg_localized_experimental
    @test one_center_materialized.route_configured_materializer_backend_consumed ==
          :pgdg_localized_experimental
    @test one_center_materialized.route_configured_materializer_d_requested == 0.15
    @test one_center_materialized.route_configured_materializer_d_consumed == 0.15
    @test one_center_materialized.route_configured_materializer_nside_requested == 5
    @test one_center_materialized.route_configured_materializer_nside_consumed == 5
    @test one_center_materialized.materialized_report.shellization_summary.source_kind ==
          :route_configured_one_center_low_order
    @test one_center_materialized.materialized_report.weight_semantics ==
          :retained_basis_integral_weights
    @test one_center_materialized.basis_bundle_export_status ==
          :supported_route_configured_one_center_basis_only_fixed_block
    @test one_center_materialized.basis_artifact_status == :not_requested
    @test one_center_materialized.ham_artifact_status == :not_requested

    mktempdir() do dir
        basisfile = joinpath(dir, "route_configured_one_center_basis.jld2")
        hamfile = joinpath(dir, "unused_ham.jld2")
        basis_status = GaussletBases._pqs_source_box_route_driver_materialization(
            _pqs_route_driver_one_center_white_lindsey_report_for_test();
            materialize_route = true,
            save_basis_artifact = true,
            save_ham_artifact = false,
            basisfile,
            hamfile,
            white_lindsey_expansion = density_expansion,
        )
        @test basis_status.status ==
              :materialized_route_configured_one_center_report_available
        @test basis_status.basis_bundle_export_status ==
              :supported_route_configured_one_center_basis_only_fixed_block
        @test basis_status.basis_artifact_status ==
              :written_route_configured_one_center_basis_only_bundle
        @test basis_status.basis_export_blocker === nothing
        @test basis_status.basis_artifact_written
        @test basis_status.basis_artifact_path == basisfile
        @test basis_status.ham_artifact_status == :not_requested
        @test !basis_status.ham_artifact_written
        @test isfile(basisfile)
        @test !isfile(hamfile)
        jldopen(basisfile, "r") do file
            top_keys = Set(
                key isa AbstractVector ? join(string.(key), "/") : string(key) for key in keys(file)
            )
            @test "basis" in top_keys
            @test "meta" in top_keys
            @test !("ham" in top_keys)
            @test String(file["basis/format"]) == "cartesian_basis_bundle_v1"
            @test String(file["basis/basis_kind"]) == "nested_fixed_block"
            @test length(file["basis/final_integral_weights"]) == 223
            @test Bool(file["meta/has_ham"]) == false
            @test String(file["meta/materialized_report_kind"]) ==
                  "white_lindsey_low_order_route_configured_one_center_report"
            @test String(file["meta/shellization_source"]) ==
                  "route_configured_one_center_low_order"
            @test Bool(file["meta/route_configured_shellization_consumed"])
            @test String(file["meta/seed_materialization_status"]) ==
                  "not_seed_route_configured_materialization"
            @test String(file["meta/export_status"]) == "basis_only"
            @test String(file["meta/basis_export_status"]) ==
                  "supported_route_configured_one_center_basis_only_fixed_block"
            @test String(file["meta/route_configured_materializer_backend_requested"]) ==
                  "pgdg_localized_experimental"
            @test String(file["meta/route_configured_materializer_backend_consumed"]) ==
                  "pgdg_localized_experimental"
            @test file["meta/route_configured_materializer_d_requested"] == 0.15
            @test file["meta/route_configured_materializer_d_consumed"] == 0.15
            @test file["meta/route_configured_materializer_nside_requested"] == 5
            @test file["meta/route_configured_materializer_nside_consumed"] == 5
        end
    end

    mktempdir() do dir
        hamfile = joinpath(dir, "route_configured_one_center_ham.jld2")
        ham_status = GaussletBases._pqs_source_box_route_driver_materialization(
            _pqs_route_driver_one_center_white_lindsey_report_for_test();
            materialize_route = true,
            save_basis_artifact = false,
            save_ham_artifact = true,
            basisfile = joinpath(dir, "unused_basis.jld2"),
            hamfile,
            white_lindsey_expansion = density_expansion,
            white_lindsey_Z = 2.0,
        )
        @test ham_status.status ==
              :materialized_route_configured_one_center_report_available
        @test ham_status.basis_bundle_export_status ==
              :supported_route_configured_one_center_basis_only_fixed_block
        @test ham_status.basis_artifact_status == :not_requested
        @test ham_status.basis_export_blocker === nothing
        @test ham_status.ham_preflight_status ==
              :available_private_low_order_ham_bundle_adapter
        @test ham_status.ham_missing_builder === nothing
        @test ham_status.ham_operator_payload_status ==
              :available_low_order_operator_payload
        @test ham_status.ham_interaction_status ==
              :available_low_order_density_density_interaction_matrix
        @test ham_status.ham_bundle_export_status ==
              :available_low_order_ham_bundle_payload
        @test ham_status.ham_export_blocker === nothing
        @test ham_status.ham_artifact_status ==
              :written_route_configured_one_center_ham_bundle
        @test ham_status.ham_artifact_written
        @test isfile(hamfile)
        @test !isfile(joinpath(dir, "unused_basis.jld2"))
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
            @test size(file["ham/overlap"]) == (223, 223)
            @test size(file["ham/one_body_hamiltonian"]) == (223, 223)
            @test size(file["ham/interaction_matrix"]) == (223, 223)
            @test length(file["ham/basis_integral_weights"]) == 223
            @test file["ham/basis_integral_weights"] == file["basis/final_integral_weights"]
            @test Bool(file["meta/has_ham"])
            @test String(file["meta/materialized_report_kind"]) ==
                  "white_lindsey_low_order_route_configured_one_center_report"
            @test String(file["meta/shellization_source"]) ==
                  "route_configured_one_center_low_order"
            @test Bool(file["meta/route_configured_shellization_consumed"])
            @test String(file["meta/seed_materialization_status"]) ==
                  "not_seed_route_configured_materialization"
            @test String(file["meta/export_status"]) == "basis_and_ham"
            @test String(file["meta/ham_preflight_status"]) ==
                  "available_private_low_order_ham_bundle_adapter"
            @test Bool(file["meta/ham_missing_builder/is_nothing"])
            @test String(file["meta/ham_operator_payload_status"]) ==
                  "available_low_order_operator_payload"
            @test String(file["meta/ham_interaction_status"]) ==
                  "available_low_order_density_density_interaction_matrix"
            @test String(file["meta/ham_export_status"]) ==
                  "available_low_order_ham_bundle_payload"
            @test Bool(file["meta/ham_export_blocker/is_nothing"])
            @test String(file["meta/route_configured_materializer_backend_requested"]) ==
                  "pgdg_localized_experimental"
            @test String(file["meta/route_configured_materializer_backend_consumed"]) ==
                  "pgdg_localized_experimental"
            @test file["meta/route_configured_materializer_d_requested"] == 0.15
            @test file["meta/route_configured_materializer_d_consumed"] == 0.15
            @test file["meta/route_configured_materializer_nside_requested"] == 5
            @test file["meta/route_configured_materializer_nside_consumed"] == 5
        end
    end

    mktempdir() do dir
        @test_throws ArgumentError GaussletBases._pqs_source_box_route_driver_materialization(
            _pqs_route_driver_one_center_white_lindsey_report_for_test();
            materialize_route = true,
            save_basis_artifact = false,
            save_ham_artifact = true,
            basisfile = joinpath(dir, "unused_basis.jld2"),
            hamfile = joinpath(dir, "missing_expansion_ham.jld2"),
        )
        @test_throws ArgumentError GaussletBases._pqs_source_box_route_driver_materialization(
            _pqs_route_driver_one_center_white_lindsey_report_for_test();
            materialize_route = true,
            save_basis_artifact = false,
            save_ham_artifact = true,
            basisfile = joinpath(dir, "unused_basis.jld2"),
            hamfile = joinpath(dir, "missing_z_ham.jld2"),
            white_lindsey_expansion = density_expansion,
        )
        @test !isfile(joinpath(dir, "missing_expansion_ham.jld2"))
        @test !isfile(joinpath(dir, "missing_z_ham.jld2"))
    end

    mktempdir() do dir
        pqs_hamfile = joinpath(dir, "pqs_pending_ham.jld2")
        pqs_ham_status = GaussletBases._pqs_source_box_route_driver_materialization(
            pqs_report;
            materialize_route = true,
            save_basis_artifact = false,
            save_ham_artifact = true,
            basisfile = joinpath(dir, "unused_basis.jld2"),
            hamfile = pqs_hamfile,
        )
        @test pqs_ham_status.ham_bundle_export_status ==
              :pending_source_box_retained_route
        @test pqs_ham_status.ham_artifact_status ==
              :not_written_pending_source_box_retained_route
        @test pqs_ham_status.ham_export_blocker == :pending_source_box_retained_route
        @test !pqs_ham_status.ham_artifact_written
        @test !isfile(pqs_hamfile)
    end

    white_lindsey_status = GaussletBases._pqs_source_box_route_driver_materialization(
        white_lindsey_report;
        materialize_route = true,
        save_basis_artifact = false,
        save_ham_artifact = false,
        basisfile = "unused_basis.jld2",
        hamfile = "unused.jld2",
    )
    @test white_lindsey_status.status == :materialized_seed_report_available
    @test white_lindsey_status.materialized_report.object_kind ==
          :white_lindsey_low_order_materialized_seed_report
    @test white_lindsey_status.materialized_report_kind ==
          :white_lindsey_low_order_materialized_seed_report
    @test white_lindsey_status.shellization_summary ===
          white_lindsey_status.materialized_report.shellization_summary
    @test white_lindsey_status.shellization_summary_available
    @test white_lindsey_status.shellization_source == :white_lindsey_one_center_seed
    @test !white_lindsey_status.route_configured_shellization_consumed
    @test white_lindsey_status.materialized_shellization_stage ==
          :route_neutral_spatial_planning
    @test white_lindsey_status.seed_materialization_status ==
          :seed_based_private_materialization
    _pqs_route_driver_check_be2_shellization_request(
        white_lindsey_status;
        current_materialization_source = :white_lindsey_one_center_seed,
    )
    @test white_lindsey_status.retained_dimension == 223
    @test white_lindsey_status.final_integral_weights_status ==
          :available_retained_basis_integral_weights
    @test white_lindsey_status.one_body_operator_status ==
          :materialized_finite_one_body_inventory
    @test white_lindsey_status.basis_bundle_export_status ==
          :supported_basis_only_fixed_block
    @test white_lindsey_status.basis_artifact_status == :not_requested
    @test !white_lindsey_status.basis_artifact_written
    @test white_lindsey_status.ham_preflight.object_kind ==
          :white_lindsey_low_order_ham_preflight
    @test white_lindsey_status.ham_preflight_status ==
          :blocked_missing_pure_low_order_interaction_builder
    @test white_lindsey_status.ham_missing_builder ==
          :missing_pure_low_order_fixed_block_density_density_interaction_builder
    @test white_lindsey_status.ham_operator_payload_status ==
          :pending_low_order_operator_payload
    @test white_lindsey_status.ham_interaction_status ==
          :pending_low_order_density_density_interaction_matrix
    @test white_lindsey_status.ham_bundle_export_status ==
          :pending_low_order_density_density_interaction_matrix
    @test white_lindsey_status.ham_artifact_status == :not_requested
    @test white_lindsey_status.ham_export_blocker ==
          :missing_pure_low_order_fixed_block_density_density_interaction_builder
    @test !white_lindsey_status.ham_artifact_written

    mktempdir() do dir
        hamfile = joinpath(dir, "white_lindsey_ham.jld2")
        save_status = GaussletBases._pqs_source_box_route_driver_materialization(
            white_lindsey_report;
            materialize_route = true,
            save_basis_artifact = false,
            save_ham_artifact = true,
            basisfile = joinpath(dir, "unused_basis.jld2"),
            hamfile,
            white_lindsey_expansion = density_expansion,
            white_lindsey_Z = 2.0,
        )
        @test save_status.ham_preflight_status ==
              :available_private_low_order_ham_bundle_adapter
        @test save_status.ham_missing_builder === nothing
        @test save_status.ham_operator_payload_status ==
              :available_low_order_operator_payload
        @test save_status.ham_interaction_status ==
              :available_low_order_density_density_interaction_matrix
        @test save_status.ham_bundle_export_status ==
              :available_low_order_ham_bundle_payload
        @test save_status.ham_artifact_status ==
              :written_white_lindsey_low_order_ham_bundle
        @test save_status.shellization_summary_available
        @test save_status.shellization_source == :white_lindsey_one_center_seed
        @test !save_status.route_configured_shellization_consumed
        @test save_status.materialized_shellization_stage ==
              :route_neutral_spatial_planning
        _pqs_route_driver_check_be2_shellization_request(
            save_status;
            current_materialization_source = :white_lindsey_one_center_seed,
        )
        @test save_status.ham_export_blocker === nothing
        @test save_status.ham_artifact_written
        @test isfile(hamfile)
        @test !isfile(joinpath(dir, "unused_basis.jld2"))
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
            @test size(file["ham/overlap"]) == (223, 223)
            @test size(file["ham/one_body_hamiltonian"]) == (223, 223)
            @test size(file["ham/interaction_matrix"]) == (223, 223)
            @test length(file["ham/basis_integral_weights"]) == 223
            @test file["ham/basis_integral_weights"] == file["basis/final_integral_weights"]
            @test Bool(file["meta/has_ham"])
            @test String(file["meta/route_family"]) == "white_lindsey_low_order"
            @test String(file["meta/export_status"]) == "basis_and_ham"
            @test String(file["meta/route_configured_shellization_request_status"]) ==
                  "metadata_only_pending_materializer"
            @test String(file["meta/route_configured_system_classification"]) ==
                  "bond_aligned_diatomic"
            @test String(file["meta/route_configured_system_classification_status"]) ==
                  "explicit_two_atom_single_axis_separation"
            @test String(file["meta/route_configured_bond_axis"]) == "x"
            @test String(file["meta/route_configured_shellization_plan_status"]) ==
                  "metadata_only_pending_materializer"
            @test String(file["meta/route_configured_shellization_planning_status"]) ==
                  "planned_metadata_only"
            @test String(file["meta/route_configured_shellization_planning_family"]) ==
                  "bond_aligned_diatomic_shellization"
            @test String(file["meta/route_configured_midpoint_slab_status"]) ==
                  "conditional_diatomic_midpoint_slab"
            @test String(file["meta/route_configured_shellization_helper_map_status"]) ==
                  "metadata_only_pending_materializer_inputs"
            @test String(file["meta/route_configured_primary_planned_helper"]) ==
                  "_nested_bond_aligned_diatomic_source"
            @test file["meta/route_configured_missing_input_count"] == 8
            @test String(file["meta/route_configured_helper_map_blocker"]) ==
                  "pending_route_configured_bond_aligned_diatomic_materializer_inputs"
            @test String(file["meta/route_configured_input_readiness_status"]) ==
                  "blocked_missing_materializer_inputs"
            @test file["meta/route_configured_available_fact_count"] == 10
            @test file["meta/route_configured_materializer_missing_input_count"] == 8
            @test String(file["meta/route_configured_input_readiness_blocker"]) ==
                  "pending_route_configured_bond_aligned_diatomic_materializer_inputs"
            @test String(file["meta/route_configured_materializer_config_status"]) ==
                  "blocked_missing_materializer_inputs"
            @test String(file["meta/route_configured_materializer_config_planning_family"]) ==
                  "bond_aligned_diatomic_shellization"
            @test file["meta/route_configured_materializer_config_pending_input_count"] == 8
            @test Bool(file["meta/shellization_summary_available"])
            @test String(file["meta/shellization_source"]) ==
                  "white_lindsey_one_center_seed"
            @test Bool(file["meta/route_configured_shellization_consumed"]) == false
            @test String(file["meta/materialized_shellization_stage"]) ==
                  "route_neutral_spatial_planning"
            @test String(file["meta/seed_materialization_status"]) ==
                  "seed_based_private_materialization"
            @test String(file["meta/ham_preflight_status"]) ==
                  "available_private_low_order_ham_bundle_adapter"
            @test Bool(file["meta/ham_missing_builder/is_nothing"])
            @test String(file["meta/ham_operator_payload_status"]) ==
                  "available_low_order_operator_payload"
            @test String(file["meta/ham_interaction_status"]) ==
                  "available_low_order_density_density_interaction_matrix"
            @test String(file["meta/ham_export_status"]) ==
                  "available_low_order_ham_bundle_payload"
            @test Bool(file["meta/ham_export_blocker/is_nothing"])
        end
    end

    mktempdir() do dir
        basisfile = joinpath(dir, "white_lindsey_basis_bundle.jld2")
        basis_status = GaussletBases._pqs_source_box_route_driver_materialization(
            white_lindsey_report;
            materialize_route = true,
            save_basis_artifact = true,
            save_ham_artifact = false,
            basisfile,
            hamfile = joinpath(dir, "unused_ham.jld2"),
        )
        @test basis_status.basis_artifact_status == :written_basis_only_bundle
        @test basis_status.basis_artifact_written
        @test basis_status.basis_artifact_path == basisfile
        @test isfile(basisfile)
        @test !isfile(joinpath(dir, "unused_ham.jld2"))
        jldopen(basisfile, "r") do file
            top_keys = Set(
                key isa AbstractVector ? join(string.(key), "/") : string(key) for key in keys(file)
            )
            @test "basis" in top_keys
            @test "meta" in top_keys
            @test !("ham" in top_keys)
            @test String(file["basis/format"]) == "cartesian_basis_bundle_v1"
            @test String(file["basis/basis_kind"]) == "nested_fixed_block"
            @test length(file["basis/final_integral_weights"]) == 223
            @test Bool(file["meta/has_ham"]) == false
            @test String(file["meta/route_family"]) == "white_lindsey_low_order"
            @test String(file["meta/benchmark_role"]) ==
                  "published_cartesian_baseline_for_pqs_comparison"
            @test String(file["meta/materialized_report_kind"]) ==
                  "white_lindsey_low_order_materialized_seed_report"
            @test String(file["meta/route_configured_shellization_request_status"]) ==
                  "metadata_only_pending_materializer"
            @test String(file["meta/route_configured_system_classification"]) ==
                  "bond_aligned_diatomic"
            @test String(file["meta/route_configured_system_classification_status"]) ==
                  "explicit_two_atom_single_axis_separation"
            @test String(file["meta/route_configured_bond_axis"]) == "x"
            @test String(file["meta/route_configured_shellization_plan_status"]) ==
                  "metadata_only_pending_materializer"
            @test String(file["meta/route_configured_shellization_planning_status"]) ==
                  "planned_metadata_only"
            @test String(file["meta/route_configured_shellization_planning_family"]) ==
                  "bond_aligned_diatomic_shellization"
            @test String(file["meta/route_configured_midpoint_slab_status"]) ==
                  "conditional_diatomic_midpoint_slab"
            @test String(file["meta/route_configured_shellization_helper_map_status"]) ==
                  "metadata_only_pending_materializer_inputs"
            @test String(file["meta/route_configured_primary_planned_helper"]) ==
                  "_nested_bond_aligned_diatomic_source"
            @test file["meta/route_configured_missing_input_count"] == 8
            @test String(file["meta/route_configured_helper_map_blocker"]) ==
                  "pending_route_configured_bond_aligned_diatomic_materializer_inputs"
            @test String(file["meta/route_configured_input_readiness_status"]) ==
                  "blocked_missing_materializer_inputs"
            @test file["meta/route_configured_available_fact_count"] == 10
            @test file["meta/route_configured_materializer_missing_input_count"] == 8
            @test String(file["meta/route_configured_input_readiness_blocker"]) ==
                  "pending_route_configured_bond_aligned_diatomic_materializer_inputs"
            @test String(file["meta/route_configured_materializer_config_status"]) ==
                  "blocked_missing_materializer_inputs"
            @test String(file["meta/route_configured_materializer_config_planning_family"]) ==
                  "bond_aligned_diatomic_shellization"
            @test file["meta/route_configured_materializer_config_pending_input_count"] == 8
            @test Bool(file["meta/shellization_summary_available"])
            @test String(file["meta/shellization_source"]) ==
                  "white_lindsey_one_center_seed"
            @test Bool(file["meta/route_configured_shellization_consumed"]) == false
            @test String(file["meta/materialized_shellization_stage"]) ==
                  "route_neutral_spatial_planning"
            @test String(file["meta/seed_materialization_status"]) ==
                  "seed_based_private_materialization"
            @test String(file["meta/export_status"]) == "basis_only"
            @test String(file["meta/basis_export_status"]) ==
                  "supported_basis_only_fixed_block"
            @test String(file["meta/ham_preflight_status"]) ==
                  "blocked_missing_pure_low_order_interaction_builder"
            @test String(file["meta/ham_missing_builder"]) ==
                  "missing_pure_low_order_fixed_block_density_density_interaction_builder"
            @test String(file["meta/ham_operator_payload_status"]) ==
                  "pending_low_order_operator_payload"
            @test String(file["meta/ham_interaction_status"]) ==
                  "pending_low_order_density_density_interaction_matrix"
            @test String(file["meta/ham_export_status"]) ==
                  "pending_low_order_density_density_interaction_matrix"
            @test String(file["meta/ham_export_blocker"]) ==
                  "missing_pure_low_order_fixed_block_density_density_interaction_builder"
        end
    end
    return nothing
end

function _pqs_route_driver_check_white_lindsey_ham_preflight()
    seed_report = GaussletBases._white_lindsey_low_order_materialized_seed_report()
    preflight = GaussletBases._pqs_source_box_route_driver_white_lindsey_ham_preflight(
        seed_report,
    )
    fixed_block_preflight =
        GaussletBases._pqs_source_box_route_driver_white_lindsey_ham_preflight(
            seed_report.fixture.fixed_block,
        )

    @test preflight == fixed_block_preflight
    @test preflight.object_kind == :white_lindsey_low_order_ham_preflight
    @test preflight.route_family == :white_lindsey_low_order
    @test occursin("_NestedFixedBlock3D", preflight.fixed_block_type_label)
    @test occursin("MappedUniformBasis", preflight.parent_basis_type_label)
    @test !preflight.ordinary_qw_fixed_block_applicable
    @test !preflight.nested_cartesian_fixed_block_applicable
    @test preflight.ordinary_cartesian_ida_builder_name_defined
    @test !preflight.ordinary_cartesian_ida_fixed_block_applicable
    @test preflight.basis_bundle_include_ham_checked
    @test !preflight.basis_bundle_ham_payload_available
    @test preflight.basis_bundle_ham_payload_status == :absent_for_fixed_block
    @test !preflight.pure_operator_payload_available
    @test preflight.status == :blocked_missing_pure_low_order_interaction_builder
    @test preflight.required_builder_contract ==
          :white_lindsey_low_order_fixed_block_density_density_builder
    @test preflight.ham_operator_payload_status == :pending_low_order_operator_payload
    @test preflight.ham_interaction_status ==
          :pending_low_order_density_density_interaction_matrix
    @test preflight.ham_bundle_export_status ==
          :pending_low_order_density_density_interaction_matrix
    @test preflight.missing_builder ==
          :missing_pure_low_order_fixed_block_density_density_interaction_builder
    @test preflight.supplement_required_paths_policy ==
          :diagnostic_only_not_benchmark_route
    @test !preflight.full_ham_export_ready

    density_expansion = coulomb_gaussian_expansion(doacc = false)
    ham_adapter =
        GaussletBases._white_lindsey_low_order_materialized_seed_ham_bundle_adapter(
            seed_report;
            expansion = density_expansion,
            Z = 2.0,
        )
    adapter_preflight =
        GaussletBases._pqs_source_box_route_driver_white_lindsey_ham_preflight(
            seed_report;
            ham_bundle_adapter = ham_adapter,
        )
    @test adapter_preflight.object_kind == :white_lindsey_low_order_ham_preflight
    @test adapter_preflight.private_writer_adapter_used
    @test adapter_preflight.private_payload_candidate_status ==
          :private_payload_candidate_not_writer_adapted
    @test adapter_preflight.basis_bundle_ham_payload_available
    @test adapter_preflight.basis_bundle_ham_payload_status ==
          :available_private_writer_adapter
    @test adapter_preflight.pure_operator_payload_available
    @test adapter_preflight.status == :available_private_low_order_ham_bundle_adapter
    @test adapter_preflight.ham_operator_payload_status ==
          :available_low_order_operator_payload
    @test adapter_preflight.ham_interaction_status ==
          :available_low_order_density_density_interaction_matrix
    @test adapter_preflight.ham_bundle_export_status ==
          :available_low_order_ham_bundle_payload
    @test adapter_preflight.missing_builder === nothing
    @test adapter_preflight.full_ham_export_ready
    return nothing
end

function _pqs_route_driver_check_materialization_report_artifacts(white_lindsey_report)
    mktempdir() do dir
        basisfile = joinpath(dir, "white_lindsey_basis_bundle.jld2")
        materialization = GaussletBases._pqs_source_box_route_driver_materialization(
            white_lindsey_report;
            materialize_route = true,
            save_basis_artifact = true,
            save_ham_artifact = false,
            basisfile,
            hamfile = joinpath(dir, "unused_ham.jld2"),
        )
        reportfile = joinpath(dir, "route_report.jld2")
        tsvfile = joinpath(dir, "route_report.tsv")
        stdout_file = joinpath(dir, "save_stdout.txt")
        open(stdout_file, "w") do io
            redirect_stdout(io) do
                GaussletBases._pqs_source_box_route_driver_save(
                    white_lindsey_report;
                    save_artifact = true,
                    save_tsv = true,
                    outfile = reportfile,
                    tsvfile,
                    materialization,
                )
            end
        end

        @test isfile(reportfile)
        @test isfile(tsvfile)
        @test isfile(basisfile)
        @test !isfile(joinpath(dir, "unused_ham.jld2"))
        jldopen(reportfile, "r") do file
            top_keys = Set(
                key isa AbstractVector ? join(string.(key), "/") : string(key) for key in keys(file)
            )
            @test "report" in top_keys
            @test "materialization" in top_keys
            saved_materialization = file["materialization"]
            @test saved_materialization.status == :materialized_seed_report_available
            @test saved_materialization.basis_artifact_status == :written_basis_only_bundle
            @test saved_materialization.shellization_summary_available
            @test saved_materialization.shellization_source == :white_lindsey_one_center_seed
            @test !saved_materialization.route_configured_shellization_consumed
            @test saved_materialization.materialized_shellization_stage ==
                  :route_neutral_spatial_planning
            @test saved_materialization.seed_materialization_status ==
                  :seed_based_private_materialization
            _pqs_route_driver_check_be2_shellization_request(
                saved_materialization;
                current_materialization_source = :white_lindsey_one_center_seed,
            )
            @test saved_materialization.ham_preflight.status ==
                  :blocked_missing_pure_low_order_interaction_builder
            @test saved_materialization.ham_preflight.required_builder_contract ==
                  :white_lindsey_low_order_fixed_block_density_density_builder
            @test saved_materialization.ham_preflight.missing_builder ==
                  :missing_pure_low_order_fixed_block_density_density_interaction_builder
            @test saved_materialization.ham_preflight_status ==
                  :blocked_missing_pure_low_order_interaction_builder
            @test saved_materialization.ham_missing_builder ==
                  :missing_pure_low_order_fixed_block_density_density_interaction_builder
            @test saved_materialization.ham_operator_payload_status ==
                  :pending_low_order_operator_payload
            @test saved_materialization.ham_interaction_status ==
                  :pending_low_order_density_density_interaction_matrix
            @test saved_materialization.ham_bundle_export_status ==
                  :pending_low_order_density_density_interaction_matrix
            @test saved_materialization.ham_export_blocker ==
                  :missing_pure_low_order_fixed_block_density_density_interaction_builder
        end

        tsv = read(tsvfile, String)
        @test occursin("route_materialization\tstatus\t:materialized_seed_report_available", tsv)
        @test occursin(
            "route_materialization\troute_configured_shellization_request_available\ttrue",
            tsv,
        )
        @test occursin(
            "route_materialization\troute_configured_shellization_request_status\t:metadata_only_pending_materializer",
            tsv,
        )
        @test occursin(
            "route_materialization\troute_configured_system_classification\t:bond_aligned_diatomic",
            tsv,
        )
        @test occursin(
            "route_materialization\troute_configured_system_classification_status\t:explicit_two_atom_single_axis_separation",
            tsv,
        )
        @test occursin(
            "route_materialization\troute_configured_bond_axis\t:x",
            tsv,
        )
        @test occursin(
            "route_materialization\troute_configured_shellization_plan_available\ttrue",
            tsv,
        )
        @test occursin(
            "route_materialization\troute_configured_shellization_plan_status\t:metadata_only_pending_materializer",
            tsv,
        )
        @test occursin(
            "route_materialization\troute_configured_shellization_planning_status\t:planned_metadata_only",
            tsv,
        )
        @test occursin(
            "route_materialization\troute_configured_shellization_planning_family\t:bond_aligned_diatomic_shellization",
            tsv,
        )
        @test occursin(
            "route_materialization\troute_configured_midpoint_slab_status\t:conditional_diatomic_midpoint_slab",
            tsv,
        )
        @test occursin(
            "route_materialization\troute_configured_shellization_helper_map_available\ttrue",
            tsv,
        )
        @test occursin(
            "route_materialization\troute_configured_shellization_helper_map_status\t:metadata_only_pending_materializer_inputs",
            tsv,
        )
        @test occursin(
            "route_materialization\troute_configured_primary_planned_helper\t:_nested_bond_aligned_diatomic_source",
            tsv,
        )
        @test occursin(
            "route_materialization\troute_configured_missing_input_count\t8",
            tsv,
        )
        @test occursin(
            "route_materialization\troute_configured_helper_map_blocker\t:pending_route_configured_bond_aligned_diatomic_materializer_inputs",
            tsv,
        )
        @test occursin(
            "route_materialization\troute_configured_input_readiness_available\ttrue",
            tsv,
        )
        @test occursin(
            "route_materialization\troute_configured_input_readiness_status\t:blocked_missing_materializer_inputs",
            tsv,
        )
        @test occursin(
            "route_materialization\troute_configured_available_fact_count\t10",
            tsv,
        )
        @test occursin(
            "route_materialization\troute_configured_materializer_missing_input_count\t8",
            tsv,
        )
        @test occursin(
            "route_materialization\troute_configured_input_readiness_blocker\t:pending_route_configured_bond_aligned_diatomic_materializer_inputs",
            tsv,
        )
        @test occursin(
            "route_materialization\troute_configured_materializer_config_available\ttrue",
            tsv,
        )
        @test occursin(
            "route_materialization\troute_configured_materializer_config_status\t:blocked_missing_materializer_inputs",
            tsv,
        )
        @test occursin(
            "route_materialization\troute_configured_materializer_config_planning_family\t:bond_aligned_diatomic_shellization",
            tsv,
        )
        @test occursin(
            "route_materialization\troute_configured_materializer_config_pending_input_count\t8",
            tsv,
        )
        @test occursin(
            "route_materialization\tshellization_summary_available\ttrue",
            tsv,
        )
        @test occursin(
            "route_materialization\tshellization_source\t:white_lindsey_one_center_seed",
            tsv,
        )
        @test occursin(
            "route_materialization\troute_configured_shellization_consumed\tfalse",
            tsv,
        )
        @test occursin(
            "route_materialization\tmaterialized_shellization_stage\t:route_neutral_spatial_planning",
            tsv,
        )
        @test occursin(
            "route_materialization\tseed_materialization_status\t:seed_based_private_materialization",
            tsv,
        )
        @test occursin("route_materialization\tbasis_artifact_status\t:written_basis_only_bundle", tsv)
        @test occursin(
            "route_materialization\tham_preflight_status\t:blocked_missing_pure_low_order_interaction_builder",
            tsv,
        )
        @test occursin(
            "route_materialization\tham_missing_builder\t:missing_pure_low_order_fixed_block_density_density_interaction_builder",
            tsv,
        )
        @test occursin(
            "route_materialization\tham_operator_payload_status\t:pending_low_order_operator_payload",
            tsv,
        )
        @test occursin(
            "route_materialization\tham_interaction_status\t:pending_low_order_density_density_interaction_matrix",
            tsv,
        )
        @test occursin(
            "route_materialization\tham_bundle_export_status\t:pending_low_order_density_density_interaction_matrix",
            tsv,
        )
        @test occursin(
            "route_materialization\tham_export_blocker\t:missing_pure_low_order_fixed_block_density_density_interaction_builder",
            tsv,
        )
        @test !occursin("route_materialization\tmaterialized_report\t", tsv)
    end
    return nothing
end

@testset "Route-driver standard unit inventory report" begin
    pqs_report = _pqs_route_driver_report_for_test()
    _pqs_route_driver_check_parent_contract(
        pqs_report;
        atom_count = 2,
        system_classification = :bond_aligned_diatomic,
        system_classification_status =
            :explicit_two_atom_single_axis_separation,
        bond_axis = :x,
        chain_axis = :x,
        parent_materialization_planning_family =
            :bond_aligned_diatomic_parent_lattice,
    )
    _pqs_route_driver_check_standard_unit_inventory(
        pqs_report;
        retained_dimension = 221,
        unit_keys = (:pqs_left, :pqs_right, :product),
        pair_count = 6,
        pair_family_counts =
            (pqs_pqs = 3, pqs_product = 2, product_pqs = 0, product_product = 1),
        retained_counts_materialized = true,
        retained_ranges_materialized = true,
        pair_families = (:pqs_pqs, :pqs_product, :product_product),
        output_representations = (:retained_two_index_density_density,),
    )
    _pqs_route_driver_check_report_output_sections(pqs_report)

    white_lindsey_report =
        _pqs_route_driver_report_for_test(route_family = :white_lindsey_low_order)
    _pqs_route_driver_check_parent_contract(
        white_lindsey_report;
        atom_count = 2,
        system_classification = :bond_aligned_diatomic,
        system_classification_status =
            :explicit_two_atom_single_axis_separation,
        bond_axis = :x,
        chain_axis = :x,
        parent_materialization_planning_family =
            :bond_aligned_diatomic_parent_lattice,
    )
    _pqs_route_driver_check_standard_unit_inventory(
        white_lindsey_report;
        retained_dimension = nothing,
        unit_keys = (:low_order_units,),
        pair_count = 1,
        pair_family_counts = (
            pqs_pqs = 0,
            pqs_product = 0,
            product_pqs = 0,
            product_product = 0,
            white_lindsey_low_order = 1,
        ),
        retained_counts_materialized = false,
        retained_ranges_materialized = false,
        pair_families = (:white_lindsey_low_order,),
        output_representations = (:low_order_nested_cartesian_basis,),
    )
    _pqs_route_driver_check_report_output_sections(white_lindsey_report)
    one_center_report = _pqs_route_driver_one_center_report_for_test()
    _pqs_route_driver_check_parent_contract(
        one_center_report;
        atom_count = 1,
        system_classification = :one_center,
        system_classification_status = :explicit_atom_count_one,
        bond_axis = nothing,
        chain_axis = nothing,
        parent_materialization_planning_family = :one_center_parent_lattice,
        constructs_basis_now = true,
        constructs_axis_bundle_now = true,
        parent_basis_materialization_status =
            :materialized_parent_objects_available,
    )
    one_center_request =
        GaussletBases._cartesian_shellization_route_configured_request(one_center_report)
    @test one_center_request.system_classification ==
          one_center_report.parent_contract.system_classification
    @test one_center_request.classification_source == :parent_contract
    _pqs_route_driver_check_white_lindsey_ham_preflight()
    _pqs_route_driver_check_materialization_status(pqs_report, white_lindsey_report)
    _pqs_route_driver_check_materialization_report_artifacts(white_lindsey_report)
end

# Integration/slow test. Do not include in default nested runner.

using Test
using GaussletBases
using JLD2

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
)
    @test report.retained_dimension == retained_dimension
    @test report.recipe_metadata.route_kind == _PQS_ROUTE_DRIVER_TEST_ROUTE_KIND
    @test Tuple(unit.unit_key for unit in report.retained_units) == unit_keys
    @test length(report.pair_entries) == pair_count
    @test report.pair_family_counts == pair_family_counts

    @test hasproperty(report, :standard_unit_inventory)
    summary = report.standard_unit_inventory
    @test summary.unit_count == length(unit_keys)
    @test summary.unit_keys == unit_keys
    @test summary.retained_counts == report.retained_counts
    @test summary.retained_dimension == retained_dimension
    @test summary.retained_counts_materialized == retained_counts_materialized
    @test summary.retained_ranges_materialized == retained_ranges_materialized
    @test summary.pair_count == pair_count
    @test summary.pair_family_counts == pair_family_counts
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
    @test materialization.route_configured_bond_axis == :x
    @test materialization.route_configured_shellization_plan_available
    @test materialization.route_configured_shellization_plan_status ==
          :metadata_only_pending_materializer
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
    @test materialization.route_configured_materializer_missing_input_count == 8
    @test materialization.route_configured_materializer_config_available
    @test materialization.route_configured_materializer_config_status ==
          :blocked_missing_materializer_inputs

    request = materialization.route_configured_shellization_request
    @test request.object_kind == :cartesian_shellization_route_configured_request
    @test request.status == materialization.route_configured_shellization_request_status
    @test request.route_kind == _PQS_ROUTE_DRIVER_TEST_ROUTE_KIND
    @test request.atom_count == 2
    @test request.system_classification == materialization.route_configured_system_classification
    @test request.bond_axis == materialization.route_configured_bond_axis
    @test request.current_materialization_source == current_materialization_source
    @test !request.route_configured_shellization_consumed
    @test !request.constructs_basis

    plan = materialization.route_configured_shellization_plan
    @test plan.object_kind == :cartesian_shellization_route_planning_stub
    @test plan.status == materialization.route_configured_shellization_plan_status
    @test plan.request_object_kind == request.object_kind
    @test plan.route_kind == request.route_kind
    @test plan.system_classification == request.system_classification
    @test plan.bond_axis == request.bond_axis
    @test !plan.constructs_basis
    @test !plan.route_configured_shellization_consumed

    helper_map = materialization.route_configured_shellization_helper_map
    @test helper_map.object_kind == :cartesian_shellization_route_planning_helper_map
    @test helper_map.status ==
          materialization.route_configured_shellization_helper_map_status
    @test helper_map.primary_planned_helper ==
          materialization.route_configured_primary_planned_helper
    @test helper_map.missing_input_count ==
          materialization.route_configured_missing_input_count
    @test helper_map.blocker == materialization.route_configured_helper_map_blocker
    @test !helper_map.route_configured_shellization_consumed

    readiness = materialization.route_configured_input_readiness
    @test readiness.object_kind ==
          :cartesian_shellization_route_materializer_input_readiness
    @test readiness.status == materialization.route_configured_input_readiness_status
    @test readiness.missing_inputs == helper_map.missing_inputs
    @test readiness.missing_input_count ==
          materialization.route_configured_materializer_missing_input_count
    @test readiness.blocker == materialization.route_configured_input_readiness_blocker
    @test !readiness.materializer_ready
    @test !readiness.route_configured_shellization_consumed
    @test !readiness.constructs_basis

    config = materialization.route_configured_materializer_config
    @test config.object_kind == :cartesian_shellization_route_materializer_config
    @test config.status == materialization.route_configured_materializer_config_status
    @test config.route_kind == request.route_kind
    @test config.system_classification == request.system_classification
    @test config.bond_axis == request.bond_axis
    @test config.primary_planned_helper == helper_map.primary_planned_helper
    @test config.pending_inputs == readiness.missing_inputs
    @test config.pending_input_count ==
          materialization.route_configured_materializer_config_pending_input_count
    @test config.blocker == readiness.blocker
    @test !config.materializer_ready
    @test !config.route_configured_shellization_consumed
    @test !config.constructs_basis
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
end

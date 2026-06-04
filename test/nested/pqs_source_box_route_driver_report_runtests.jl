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
const _PQS_ROUTE_DRIVER_TEST_ROUTE_KIND = :be2_cartesian_nesting_route_driver_spine

function _pqs_route_driver_report_for_test(; route_family = :pqs_source_box)
    return GaussletBases._pqs_source_box_route_driver_dry_run(
        ;
        route_family,
        route_kind = _PQS_ROUTE_DRIVER_TEST_ROUTE_KIND,
        atom_symbols = ("Be", "Be"),
        nuclear_charges = (4, 4),
        atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
        radius = 15.0,
        parent_axis_counts = (x = 9, y = 7, z = 9),
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
        @test occursin("standard_unit_inventory\tunit_count", tsv)
        @test occursin("retained_unit\t", tsv)
        @test occursin("pair_entry\t", tsv)
    finally
        rm(tmpdir; recursive = true, force = true)
    end
    @test !isfile("pqs_source_box_route_driver_report.jld2")
    @test !isfile("pqs_source_box_route_driver_report.tsv")
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
    _pqs_route_driver_check_white_lindsey_ham_preflight()
    _pqs_route_driver_check_materialization_status(pqs_report, white_lindsey_report)
    _pqs_route_driver_check_materialization_report_artifacts(white_lindsey_report)
end

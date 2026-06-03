using Test

@testset "Be2 strict-PQS q5 source metadata acceptance opt-in" begin
    artifact_dir = strip(get(ENV, "BE2_PQS_Q5_ARTIFACT_DIR", ""))
    if isempty(artifact_dir)
        @test_skip "set BE2_PQS_Q5_ARTIFACT_DIR to run the real-artifact source metadata acceptance probe"
    else
        probe_path = normpath(
            joinpath(
                @__DIR__,
                "..",
                "..",
                "tmp",
                "work",
                "be2_pqs_source_metadata_acceptance_probe.jl",
            ),
        )
        @test isfile(probe_path)
        @test isdir(artifact_dir)
        if isfile(probe_path) && isdir(artifact_dir)
            probe_module =
                Module(:Be2PQSQ5SourceMetadataAcceptanceOptInProbe)
            Core.eval(
                probe_module,
                :(include(path::AbstractString) = Base.include(@__MODULE__, path)),
            )
            Base.include(probe_module, probe_path)
            probe_result = Core.eval(
                probe_module,
                quote
                    checks = run_acceptance($artifact_dir)
                    write_report(ACCEPTANCE_REPORT_PATH, checks)
                    write_summary(ACCEPTANCE_SUMMARY_PATH, checks)
                    (failures = copy(checks.failures), rows = Dict(checks.rows))
                end,
            )
            @test isempty(probe_result.failures)
            rows = probe_result.rows
            required_true_checks = (
                "fixed_dimension_is_1483",
                "source_shells_table_present",
                "source_modes_table_present",
                "source_shell_count_is_8",
                "source_mode_count_is_2299",
                "product_doside_source_shells_is_3",
                "product_doside_source_modes_is_531",
                "support_dense_source_shells_is_2",
                "support_dense_source_modes_is_1458",
                "shell_realized_pqs_source_shells_is_3",
                "shell_realized_pqs_source_modes_is_310",
                "unavailable_source_units_is_0",
                "unavailable_source_columns_is_0",
                "source_inventory_no_center_inference",
                "source_inventory_no_nearest_grid_inference",
                "source_inventory_no_support_order_inference",
                "source_inventory_no_support_index_inference",
                "source_inventory_no_raw_to_final_inference",
                "source_shell_shell_statuses_unavailable",
                "source_shell_ray_statuses_unavailable",
                "source_shell_radial_statuses_unavailable",
                "source_mode_shell_statuses_unavailable",
                "source_mode_ray_statuses_unavailable",
                "source_mode_radial_statuses_unavailable",
                "product_doside_relations_status",
                "product_doside_relation_kinds_only",
                "product_axis_tuples_not_rays",
                "product_axis_tuples_are_not_ray_labels",
                "fixed_column_source_relations_product_doside_only",
            )
            @test all(
                check_name -> get(rows, check_name, false) === true,
                required_true_checks,
            )
        end
    end
end

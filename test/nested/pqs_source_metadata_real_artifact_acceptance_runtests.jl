using Test

include("pqs_source_metadata_real_artifact_acceptance_support.jl")

@testset "Be2 strict-PQS q5 source metadata acceptance opt-in" begin
    artifact_dir = strip(get(ENV, "BE2_PQS_Q5_ARTIFACT_DIR", ""))
    if isempty(artifact_dir)
        @test_skip "set BE2_PQS_Q5_ARTIFACT_DIR to run the real-artifact source metadata acceptance probe"
    else
        @test isdir(artifact_dir)
        if isdir(artifact_dir)
            probe_result =
                _be2_pqs_q5_source_metadata_acceptance(artifact_dir)
            @test isempty(probe_result.failures)
            rows = probe_result.row_dict
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
                "source_mode_header_has_explicit_parent_lattice_axes",
                "source_mode_local_axes_in_contracted_range",
                "source_mode_parent_lattice_axis_statuses_explicit",
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

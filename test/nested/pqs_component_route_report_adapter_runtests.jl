using Test
using GaussletBases

@testset "PQS component-route report adapter" begin
    metrics_module = GaussletBases.CartesianContractedParentMetrics

    function _synthetic_no_go(; overrides = (;))
        return merge(
            (
                source_box_first = true,
                source_box_algorithmic_path_true_for_every_pair = true,
                shell_projection_used = false,
                lowdin_cleanup_used = false,
                support_local_oracle_used = false,
                support_local_pqs_oracle_used = false,
                support_local_shell_row_algorithm = false,
                support_coefficient_matrix_used = false,
                shell_row_algorithm = false,
                retained_pqs_weights_used = false,
                retained_pqs_weights_positive_checked = false,
                retained_weight_division_allowed = false,
                retained_pqs_weight_division_allowed = false,
                ida_weight_division_allowed = false,
                packet_adoption = false,
                fixed_block_routing = false,
                qwhamiltonian_consumes = false,
                hamiltonian_matrix_built = false,
                public_default_consumes = false,
                mwg_supplement_residual_path = false,
                mwg_supplement_residual_provenance_adapted = false,
                ecp_terms_implemented = false,
                cr2_science_status_changed = false,
                dense_parent_projection_algorithmic = false,
            ),
            overrides,
        )
    end

    function _synthetic_summary(mode::Symbol)
        dense_mode = mode == :density_normalized
        return (
            object_kind = :pqs_pqs_product_component_route_smoke_summary,
            pair_factor_normalization = mode,
            retained_dimension = 221,
            nuclear_pair_count = 6,
            electron_electron_pair_count = 6,
            ida_term_count = 45,
            finite_checks = (
                output_finite = true,
                nuclear_output_finite = true,
                electron_electron_output_finite = true,
            ),
            symmetry_errors = (
                nuclear = 1.0e-16,
                electron_electron = 0.0,
            ),
            nuclear_total_from_center_error = 0.0,
            dense_parent_ida_authority = (
                available = dense_mode,
                max_error = dense_mode ? 1.33e-15 : nothing,
                within_tolerance = dense_mode ? true : nothing,
                skip_reason =
                    dense_mode ? nothing : :density_normalized_authority_only,
                validation_only = dense_mode,
            ),
            source_weight_division_owner =
                dense_mode ? :caller_supplied_density_normalized_pair_factors :
                :pgdg_auxiliary_source_weights,
            source_weight_division_applied_by_helper = !dense_mode,
            no_go_diagnostics = _synthetic_no_go(),
        )
    end

    function _synthetic_mwg_facts(; overrides = Dict{String,String}())
        facts = Dict(
            "route" =>
                "ordinary_cartesian_qiu_white_operators synthetic Be2 path",
            "residual_nucleus_indices" => "[1, 1, 1, 2, 2, 2]",
            "residual_owner_counts" => "{1=>3, 2=>3}",
            "owner_metadata_source" => "synthetic explicit owners",
            "owner_semantics_inferred_from_raw_to_final_support" => "false",
            "component_helper" => "_qwrg_final_residual_mwg_component_blocks",
            "max_authority_error" => "0",
            "fixed_fixed_shape" => "(125, 125)",
            "fixed_residual_shape" => "(125, 6)",
            "residual_residual_shape" => "(6, 6)",
            "final_interaction_shape" => "(131, 131)",
            "diagnostics.raw_gto_rows_role" =>
                "residual_construction_inputs_only",
            "diagnostics.raw_gto_gto_mwg_interaction_blocks_used" => "false",
            "diagnostics.fixed_raw_gto_mwg_interaction_blocks_used" => "false",
        )
        merge!(facts, overrides)
        return facts
    end

    summaries = (
        _synthetic_summary(:density_normalized),
        _synthetic_summary(:raw_weighted),
    )
    facts = _synthetic_mwg_facts()
    report = metrics_module._pqs_pqs_product_component_route_smoke_report_adapter(
        summaries,
        facts;
        generated_at = "synthetic",
        source_report = "synthetic_mwg_report.txt",
        parent_dims = (5, 5, 7),
        source_mode_dims = (5, 5, 5),
        left_source_box = (1:5, 1:5, 1:5),
        right_source_box = (1:5, 1:5, 3:7),
        product_source_box = (1:5, 1:5, 4:4),
        nuclear_centers = ((0.0, 0.0, -2.5), (0.0, 0.0, 2.5)),
        nuclear_charges = (4.0, 4.0),
    )

    @test report.object_kind ==
          :pqs_pqs_product_component_route_smoke_report_adapter
    @test report.status == :private_component_route_smoke_report
    @test report.source_box_pqs_ida_component_smoke.rows[1].retained_dimension ==
          221
    @test report.source_box_pqs_ida_component_smoke.rows[1].no_go_clear
    @test report.source_box_pqs_ida_component_smoke.rows[2].no_go_clear
    @test report.source_box_pqs_ida_component_smoke.rows[1].dense_parent_ida_authority_max_error ==
          1.33e-15
    @test report.final_residual_mwg_supplement_component_facts.residual_nucleus_indices ==
          "[1, 1, 1, 2, 2, 2]"
    @test report.final_residual_mwg_supplement_component_facts.max_authority_error ==
          "0"
    @test !report.lane_boundaries.pqs_source_box_ida_and_mwg_residual_same_algorithm
    @test report.lane_boundaries.no_owner_inference_from_raw_to_final_support
    @test report.lane_boundaries.no_raw_gto_gto_mwg_blocks
    @test report.lane_boundaries.no_fixed_raw_gto_mwg_blocks
    @test report.lane_boundaries.no_retained_weight_or_ida_division
    @test report.lane_boundaries.no_packet_fixed_block_qw_hamiltonian_public_adoption
    @test report.lane_boundaries.no_ecp_scf_hf_cr2_science_claim
    @test report.lane_boundaries.no_mwg_ida_semantic_change
    @test report.diagnostics.source_box_rows_all_finite_and_no_go_clear
    @test report.diagnostics.final_residual_mwg_authority_error_zero
    @test report.diagnostics.final_residual_mwg_owner_vector_available
    @test !report.diagnostics.packet_adoption
    @test !report.diagnostics.fixed_block_routing
    @test !report.diagnostics.qwhamiltonian_consumes
    @test !report.diagnostics.hamiltonian_matrix_built
    @test !report.diagnostics.public_default_consumes
    @test !report.diagnostics.ecp_terms_implemented
    @test !report.diagnostics.scf_hf_validation_claim
    @test !report.diagnostics.cr2_science_status_changed
    @test !report.diagnostics.science_route_adoption
    @test !report.diagnostics.mwg_ida_semantics_changed

    io = IOBuffer()
    metrics_module._write_pqs_pqs_product_component_route_smoke_report(io, report)
    text = String(take!(io))
    @test occursin("[source_box_pqs_ida_component_smoke]", text)
    @test occursin(
        "[source_box_pqs_ida_component_smoke.density_normalized]",
        text,
    )
    @test occursin("[source_box_pqs_ida_component_smoke.raw_weighted]", text)
    @test occursin("[final_residual_mwg_supplement_component_facts]", text)
    @test occursin("[lane_boundaries]", text)
    @test occursin("status\tprivate_summary_only", text)
    @test occursin("residual_nucleus_indices\t[1, 1, 1, 2, 2, 2]", text)

    missing_owner_facts = _synthetic_mwg_facts()
    delete!(missing_owner_facts, "residual_nucleus_indices")
    @test_throws ArgumentError metrics_module._pqs_pqs_product_component_route_smoke_report_adapter(
        summaries,
        missing_owner_facts,
    )
    @test_throws ArgumentError metrics_module._pqs_pqs_product_component_route_smoke_report_adapter(
        summaries,
        _synthetic_mwg_facts(; overrides = Dict("max_authority_error" => "1.0e-3")),
    )
    @test_throws ArgumentError metrics_module._pqs_pqs_product_component_route_smoke_report_adapter(
        summaries,
        _synthetic_mwg_facts(;
            overrides = Dict(
                "owner_semantics_inferred_from_raw_to_final_support" => "true",
            ),
        ),
    )
    @test_throws ArgumentError metrics_module._pqs_pqs_product_component_route_smoke_report_adapter(
        summaries,
        _synthetic_mwg_facts(;
            overrides = Dict(
                "diagnostics.raw_gto_gto_mwg_interaction_blocks_used" => "true",
            ),
        ),
    )
    @test_throws ArgumentError metrics_module._pqs_pqs_product_component_route_smoke_report_adapter(
        summaries,
        _synthetic_mwg_facts(;
            overrides = Dict(
                "diagnostics.fixed_raw_gto_mwg_interaction_blocks_used" =>
                    "true",
            ),
        ),
    )
end

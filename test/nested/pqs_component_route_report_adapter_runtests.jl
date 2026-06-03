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

    function _synthetic_retained_units()
        return (
            (
                unit_key = :pqs_left,
                retained_unit_kind = :pqs,
                source_family = :mode_selected_raw_product_box,
                retained_rule_kind = :boundary_comx_product_mode_selection,
                retained_range = 1:98,
                source_dimensions = (5, 5, 5),
                source_dimension = 125,
                axis_intervals = (1:5, 1:5, 1:5),
                source_mode_ordering = :z_fast,
                boundary_mode_count = 98,
                retained_count = 98,
                supported_safe_terms = (:overlap, :kinetic),
            ),
            (
                unit_key = :pqs_right,
                retained_unit_kind = :pqs,
                source_family = :mode_selected_raw_product_box,
                retained_rule_kind = :boundary_comx_product_mode_selection,
                retained_range = 99:196,
                source_dimensions = (5, 5, 5),
                source_dimension = 125,
                axis_intervals = (1:5, 1:5, 3:7),
                source_mode_ordering = :z_fast,
                boundary_mode_count = 98,
                retained_count = 98,
                supported_safe_terms = (:overlap, :kinetic),
            ),
            (
                unit_key = :product,
                retained_unit_kind = :product_doside,
                source_family = :product_doside,
                retained_rule_kind = :product_doside,
                retained_range = 197:221,
                source_dimensions = (5, 5, 1),
                source_dimension = 25,
                axis_intervals = (1:5, 1:5, 4:4),
                retained_count = 25,
                supported_safe_terms = (:overlap, :kinetic),
            ),
        )
    end

    function _synthetic_summary(mode::Symbol)
        dense_mode = mode == :density_normalized
        return (
            object_kind = :pqs_pqs_product_component_route_smoke_summary,
            pair_factor_normalization = mode,
            route_shape = (:pqs_left, :pqs_right, :product),
            retained_units = _synthetic_retained_units(),
            ranges = (pqs_left = 1:98, pqs_right = 99:196, product = 197:221),
            retained_dimension = 221,
            nuclear_pair_count = 6,
            nuclear_pair_family_counts =
                (pqs_pqs = 3, pqs_product = 2, product_product = 1),
            electron_electron_pair_count = 6,
            electron_electron_pair_family_counts =
                (pqs_pqs = 3, pqs_product = 2, product_product = 1),
            helper_used_for_nuclear_pair_families = (
                pqs_pqs = :_pqs_pqs_source_box_nuclear_attraction_by_center,
                pqs_product =
                    :_pqs_product_source_box_nuclear_attraction_by_center,
                product_product =
                    :_product_doside_source_box_nuclear_attraction_by_center,
            ),
            helper_used_for_electron_electron_pair_families = (
                pqs_pqs =
                    :_pqs_pqs_source_box_density_density_interaction_block,
                pqs_product =
                    :_pqs_product_source_box_density_density_interaction_block,
                product_product =
                    :_product_doside_source_box_density_density_interaction_block,
            ),
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
            source_weight_division_shape =
                dense_mode ? nothing : :axis_pair_weight_outer,
            no_go_diagnostics = _synthetic_no_go(),
            performance = (
                nuclear = (
                    elapsed_seconds = dense_mode ? 0.011 : 0.012,
                    allocated_bytes = dense_mode ? 1100 : 1200,
                    gc_time_seconds = 0.0,
                ),
                electron_electron = (
                    elapsed_seconds = dense_mode ? 0.021 : 0.022,
                    allocated_bytes = dense_mode ? 2100 : 2200,
                    gc_time_seconds = 0.0,
                ),
            ),
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

    function _synthetic_fixed_side_inventory_unit(
        role::Symbol,
        category::Symbol,
        kind::Symbol,
        column_range::UnitRange{Int};
        support_count::Int,
        primitive_family = nothing,
        representation_kind = category,
        staged_unit = nothing,
        raw_box_auxiliary_metadata = nothing,
        shell_transform = nothing,
        shell_row_oracle_only::Bool = false,
        support_local_oracle_used::Bool = false,
        source_axis_center_vectors = nothing,
        support_states = nothing,
    )
        source_center_convention =
            isnothing(source_axis_center_vectors) ? :unavailable : :comx_construction
        source_center_status =
            isnothing(source_axis_center_vectors) ? :unavailable : :native_representative
        return (
            role = role,
            category = category,
            kind = kind,
            column_range = column_range,
            retained_count = length(column_range),
            support_count = support_count,
            support_states = isnothing(support_states) ?
                NTuple{3,Int}[] :
                NTuple{3,Int}[state for state in support_states],
            staged_unit = staged_unit,
            support_source_semantics = category == :support_dense ?
                :support_local_direct_rows : :synthetic_source_semantics,
            safe_term_capability = category == :shell_realized_pqs_fixture ?
                :support_local_oracle_for_shell_realization :
                :synthetic_safe_terms,
            active_representation_stage = representation_kind,
            raw_box_auxiliary_metadata = raw_box_auxiliary_metadata,
            shell_realization_transform_fact = shell_transform,
            source_axis_center_vectors = source_axis_center_vectors,
            source_center_convention = source_center_convention,
            source_center_status = source_center_status,
            raw_product_box_operator_contract = category == :product_doside,
            route_descriptor_emitted = false,
            construction_mutated = false,
            sidecar_installation = false,
            packet_adoption = false,
            provenance = (
                source = :synthetic_fixed_side_inventory_unit,
                primitive_family = primitive_family,
            ),
            diagnostics = (
                retained_weight_semantics = :not_positive_quadrature_weights,
                ida_weight_division_allowed = false,
                shell_row_oracle_only = shell_row_oracle_only,
                support_local_oracle_used = support_local_oracle_used,
            ),
        )
    end

    function _synthetic_support_dense_source_states()
        return NTuple{3,Int}[
            (ix, iy, iz) for iz in 2:4 for iy in 5:7 for ix in 4:6
        ]
    end

    function _synthetic_shell_realized_boundary_modes()
        return NTuple{3,Int}[
            (1, 1, 1),
            (1, 1, 5),
            (1, 5, 1),
            (5, 1, 1),
            (5, 5, 1),
            (5, 1, 5),
            (1, 5, 5),
        ]
    end

    function _synthetic_product_axis(kind::Symbol, fixed_index, interval)
        return (
            kind = kind,
            fixed_index = fixed_index,
            interval = interval,
        )
    end

    function _synthetic_product_staged_unit()
        return (
            axes = (
                _synthetic_product_axis(:fixed, 1, nothing),
                _synthetic_product_axis(:active, nothing, 1:2),
                _synthetic_product_axis(:active, nothing, 3:4),
            ),
            axis_function_indices = ((1, 1, 1), (1, 2, 1), (1, 1, 2), (1, 2, 2)),
        )
    end

    function _synthetic_fixed_side_inventory(;
        shell_source_box_operator_ready::Bool = false,
        product_source_axis_center_vectors = nothing,
    )
        units = (
            _synthetic_fixed_side_inventory_unit(
                :outer_mismatch_z_low_slab,
                :product_doside,
                :product_doside,
                1:4;
                support_count = 4,
                primitive_family = :outer_mismatch_boundary_slab_set,
                representation_kind = :product_doside_bridge,
                staged_unit = _synthetic_product_staged_unit(),
                source_axis_center_vectors = product_source_axis_center_vectors,
            ),
            _synthetic_fixed_side_inventory_unit(
                :left_atom_box,
                :support_dense,
                :atom_core_cube,
                5:8;
                support_count = 27,
                primitive_family = :atom_local_complete_shell_sequence,
                representation_kind = :support_dense_direct_support,
                support_states = _synthetic_support_dense_source_states(),
            ),
            _synthetic_fixed_side_inventory_unit(
                :regular_shared_molecular_shell_1,
                :shell_realized_pqs_fixture,
                :projected_q_shell,
                9:15;
                support_count = 50,
                primitive_family = :projected_q_shell,
                representation_kind = :shell_realized_pqs_fixture,
                raw_box_auxiliary_metadata = (
                    available = true,
                    source_mode_dims = (5, 5, 5),
                    axis_intervals = (10:24, 20:34, 30:54),
                    boundary_mode_indices =
                        _synthetic_shell_realized_boundary_modes(),
                    selection_rule = :synthetic_boundary_fixture,
                    reference_only = true,
                    active_current_route_contract = false,
                ),
                shell_transform = (
                    source_box_operator_application_ready =
                        shell_source_box_operator_ready,
                    compact_source_space_transform = (available = false,),
                ),
                shell_row_oracle_only = true,
                support_local_oracle_used = true,
            ),
        )
        return (
            object_kind = :pqs_current_route_retained_unit_inventory_fixture,
            status = :private_diagnostic_only,
            units = units,
            coverage = (
                first_column = 1,
                last_column = 15,
                represented_count = 15,
                covers_every_column_once = true,
            ),
            diagnostics = (
                fixed_dimension = 15,
                coverage_complete = true,
            ),
        )
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
    @test report.source_box_pqs_ida_component_smoke.retained_unit_count == 3
    @test report.source_box_pqs_ida_component_smoke.source_unit_label_status ==
          :explicit_route_descriptor_unit_keys
    @test report.source_box_pqs_ida_component_smoke.source_unit_labels ==
          (:pqs_left, :pqs_right, :product)
    @test map(
        unit -> unit.retained_range,
        report.source_box_pqs_ida_component_smoke.retained_units,
    ) == (1:98, 99:196, 197:221)
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

    sidecar =
        metrics_module._pqs_pqs_product_component_route_smoke_cr2_sidecar_schema(
            report;
            provenance = (source = :synthetic_sidecar_test,),
        )
    @test sidecar.object_kind ==
          :pqs_pqs_product_component_route_smoke_cr2_sidecar_schema
    @test sidecar.status == :private_cr2_sidecar_schema
    @test sidecar.lanes.source_box_pqs_ida_fixed_side.retained_dimension == 221
    @test sidecar.lanes.source_box_pqs_ida_fixed_side.source_unit_label_status ==
          :explicit_route_descriptor_unit_keys
    @test sidecar.lanes.source_box_pqs_ida_fixed_side.source_unit_labels ==
          (:pqs_left, :pqs_right, :product)
    @test map(
        unit -> unit.retained_range,
        sidecar.lanes.source_box_pqs_ida_fixed_side.retained_units,
    ) == (1:98, 99:196, 197:221)
    @test length(sidecar.lanes.source_box_pqs_ida_fixed_side.components) == 2
    @test sidecar.lanes.source_box_pqs_ida_fixed_side.components[1].electron_electron_density_density.representation ==
          :retained_two_index_density_density
    @test sidecar.lanes.source_box_pqs_ida_fixed_side.components[1].electron_electron_density_density.not_four_index_galerkin_tensor
    @test sidecar.lanes.source_box_pqs_ida_fixed_side.components[1].nuclear_attraction_by_center.helper_by_family.pqs_pqs ==
          :_pqs_pqs_source_box_nuclear_attraction_by_center
    @test isnothing(sidecar.fixed_side_retained_unit_metadata)
    @test isnothing(sidecar.source_shell_mode_inventory)
    @test sidecar.lanes.final_residual_mwg_supplement.fixed_dimension == 125
    @test sidecar.lanes.final_residual_mwg_supplement.residual_dimension == 6
    @test sidecar.lanes.final_residual_mwg_supplement.final_dimension == 131
    @test sidecar.lanes.final_residual_mwg_supplement.fixed_column_range == 1:125
    @test sidecar.lanes.final_residual_mwg_supplement.residual_column_range ==
          126:131
    @test sidecar.lanes.final_residual_mwg_supplement.final_column_range == 1:131
    @test sidecar.lanes.final_residual_mwg_supplement.residual_owner_metadata.residual_nucleus_indices ==
          (1, 1, 1, 2, 2, 2)
    @test sidecar.lanes.final_residual_mwg_supplement.residual_owner_metadata.owner_count_matches_residual_rows
    @test !sidecar.lanes.final_residual_mwg_supplement.residual_owner_metadata.owner_semantics_inferred_from_raw_to_final_support
    @test sidecar.lanes.final_residual_mwg_supplement.components.fixed_residual.shape ==
          (125, 6)
    @test sidecar.labels.source_unit_label_status ==
          :explicit_route_descriptor_unit_keys
    @test sidecar.labels.shell_label_status == :unavailable
    @test !sidecar.labels.label_reconstruction_from_centers
    @test !sidecar.labels.nearest_grid_or_center_label_heuristic
    @test sidecar.absences_by_contract.raw_gto_gto_mwg_interaction_blocks
    @test sidecar.absences_by_contract.fixed_raw_gto_mwg_interaction_blocks
    @test sidecar.absences_by_contract.owner_inference_from_raw_to_final_support
    @test sidecar.absences_by_contract.retained_weight_ida_division
    @test sidecar.absences_by_contract.packet_fixed_block_qw_hamiltonian_adoption
    @test sidecar.diagnostics.lanes_remain_separate
    @test sidecar.diagnostics.source_box_unit_records_available
    @test !sidecar.diagnostics.fixed_side_retained_unit_metadata_available
    @test !sidecar.diagnostics.source_shell_mode_inventory_available
    @test sidecar.diagnostics.source_shell_count == 0
    @test sidecar.diagnostics.source_mode_count == 0
    @test sidecar.diagnostics.source_shell_mode_center_status == :unavailable
    @test !sidecar.diagnostics.source_shell_mode_product_doside_only
    @test sidecar.diagnostics.final_residual_shape_consistent
    @test sidecar.diagnostics.residual_owner_rows_match
    @test !sidecar.diagnostics.packet_adoption
    @test !sidecar.diagnostics.fixed_block_routing
    @test !sidecar.diagnostics.qwhamiltonian_consumes
    @test !sidecar.diagnostics.hamiltonian_matrix_built
    @test !sidecar.diagnostics.public_default_consumes
    @test !sidecar.diagnostics.ecp_terms_implemented
    @test !sidecar.diagnostics.scf_hf_validation_claim
    @test !sidecar.diagnostics.cr2_science_status_changed
    @test !sidecar.diagnostics.mwg_ida_semantics_changed

    readiness =
        metrics_module._pqs_pqs_product_private_source_box_route_adapter_readiness_summary(
            report,
        )
    @test readiness.object_kind ==
          :pqs_pqs_product_private_source_box_route_adapter_readiness_summary
    @test readiness.status == :private_route_adapter_inputs_ready
    @test readiness.route_shape ==
          "left raw-box PQS | product/doside slab | right raw-box PQS"
    @test readiness.retained_dimension == 221
    @test readiness.retained_unit_count == 3
    @test readiness.pair_factor_normalization_modes ==
          (:density_normalized, :raw_weighted)
    @test isempty(readiness.missing_required_pieces)
    @test :cr2_sidecar_schema in readiness.missing_optional_pieces
    @test readiness.required_input_availability.source_box_fixed_side_facts
    @test readiness.required_input_availability.by_center_nuclear_attraction
    @test readiness.required_input_availability.ida_source_box_electron_electron
    @test readiness.required_input_availability.source_box_algorithmic_path
    @test readiness.required_input_availability.final_residual_mwg_component_facts
    @test readiness.required_input_availability.final_residual_mwg_owner_metadata
    @test readiness.required_input_availability.final_residual_mwg_authority
    @test readiness.required_input_availability.authority_comparison_accounted_for
    @test readiness.optional_input_availability.timing_allocation_fields
    @test readiness.optional_input_availability.dense_parent_ida_authority_comparison
    @test readiness.ida_source_box_electron_electron.authority_comparison.available_modes ==
          (:density_normalized,)
    @test readiness.ida_source_box_electron_electron.authority_comparison.skip_reasons ==
          (:density_normalized_authority_only,)
    @test readiness.timing_allocation.available
    @test readiness.timing_allocation.rows[1].nuclear_allocated_bytes == 1100
    @test isempty(readiness.no_go_violations)
    @test !readiness.no_go_flags.public_default_behavior
    @test !readiness.no_go_flags.packet_fixed_block_qw_hamiltonian_adoption
    @test !readiness.no_go_flags.mwg_ida_semantic_change
    @test !readiness.no_go_flags.retained_weight_division
    @test !readiness.no_go_flags.raw_gto_gto_mwg_blocks
    @test !readiness.no_go_flags.fixed_raw_gto_mwg_blocks
    @test !readiness.no_go_flags.owner_shell_ray_inference
    @test !readiness.no_go_flags.cr2_science_status_changed
    @test readiness.ready_for_next_private_adapter_pass
    @test !readiness.diagnostics.builds_route_matrices
    @test !readiness.diagnostics.construction_behavior_changed
    @test readiness.diagnostics.lanes_remain_separate

    readiness_with_sidecar =
        metrics_module._pqs_pqs_product_private_source_box_route_adapter_readiness_summary(
            report;
            cr2_sidecar = sidecar,
    )
    @test readiness_with_sidecar.cr2_sidecar.available
    @test readiness_with_sidecar.optional_input_availability.cr2_sidecar_schema
    @test !in(
        :cr2_sidecar_schema,
        readiness_with_sidecar.missing_optional_pieces,
    )
    @test isempty(readiness_with_sidecar.no_go_violations)
    @test readiness_with_sidecar.ready_for_next_private_adapter_pass

    sidecar_io = IOBuffer()
    metrics_module._write_pqs_pqs_product_component_route_smoke_cr2_sidecar_schema_report(
        sidecar_io,
        sidecar,
    )
    sidecar_text = String(take!(sidecar_io))
    @test occursin("[source_box_pqs_ida_fixed_side]", sidecar_text)
    @test occursin("[final_residual_mwg_supplement]", sidecar_text)
    @test occursin("[labels]", sidecar_text)
    @test occursin("[boundaries]", sidecar_text)
    @test occursin(
        "source_unit_label_status\texplicit_route_descriptor_unit_keys",
        sidecar_text,
    )
    @test occursin("source_unit_labels\t(:pqs_left, :pqs_right, :product)", sidecar_text)
    @test occursin("shell_label_status\tunavailable", sidecar_text)
    @test occursin("label_reconstruction_from_centers\tfalse", sidecar_text)
    @test occursin("nearest_grid_or_center_label_heuristic\tfalse", sidecar_text)
    @test occursin("fixed_dimension\t125", sidecar_text)
    @test occursin("residual_dimension\t6", sidecar_text)
    @test occursin("final_dimension\t131", sidecar_text)
    @test occursin("fixed_column_range\t1:125", sidecar_text)
    @test occursin("residual_column_range\t126:131", sidecar_text)
    @test occursin("residual_nucleus_indices\t(1, 1, 1, 2, 2, 2)", sidecar_text)
    @test occursin("owner_count_matches_residual_rows\ttrue", sidecar_text)
    @test occursin("component.fixed_residual.shape\t(125, 6)", sidecar_text)
    @test occursin(
        "raw_gto_gto_mwg_interaction_blocks_absent_by_contract\ttrue",
        sidecar_text,
    )
    @test occursin(
        "owner_inference_from_raw_to_final_support_absent_by_contract\ttrue",
        sidecar_text,
    )
    @test occursin("packet_adoption\tfalse", sidecar_text)
    @test occursin("qwhamiltonian_consumes\tfalse", sidecar_text)
    @test occursin("cr2_science_status_changed\tfalse", sidecar_text)
    @test !occursin("[fixed_side_retained_unit_metadata]", sidecar_text)
    @test !occursin("[source_shell_mode_inventory]", sidecar_text)

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
    @test occursin("source_unit_label_status\texplicit_route_descriptor_unit_keys", text)
    @test occursin("residual_nucleus_indices\t[1, 1, 1, 2, 2, 2]", text)

    fixed_side_metadata =
        metrics_module._pqs_current_route_fixed_side_retained_unit_metadata(
            _synthetic_fixed_side_inventory();
            provenance = (source = :synthetic_fixed_side_metadata_test,),
        )
    @test fixed_side_metadata.object_kind ==
          :pqs_current_route_fixed_side_retained_unit_metadata
    @test fixed_side_metadata.status == :private_fixed_side_retained_unit_metadata
    @test fixed_side_metadata.fixed_dimension == 15
    @test fixed_side_metadata.unit_count == 3
    @test fixed_side_metadata.labels.source_unit_label_status ==
          :explicit_inventory_unit_keys
    @test fixed_side_metadata.labels.source_unit_labels == (
        :outer_mismatch_z_low_slab,
        :left_atom_box,
        :regular_shared_molecular_shell_1,
    )
    @test fixed_side_metadata.labels.shell_label_status == :unavailable
    @test isempty(fixed_side_metadata.labels.shell_labels)
    @test !fixed_side_metadata.labels.label_reconstruction_from_centers
    @test !fixed_side_metadata.labels.nearest_grid_or_center_label_heuristic
    @test map(record -> record.retained_range, fixed_side_metadata.retained_units) ==
          (1:4, 5:8, 9:15)
    @test map(record -> record.retained_count, fixed_side_metadata.retained_units) ==
          (4, 4, 7)
    @test fixed_side_metadata.retained_units[1].is_product_doside
    @test fixed_side_metadata.retained_units[2].is_support_dense_direct_support
    @test fixed_side_metadata.retained_units[3].is_shell_realized_pqs_fixture
    @test fixed_side_metadata.retained_units[3].source_mode_dims == (5, 5, 5)
    @test fixed_side_metadata.retained_units[3].source_dimension == 125
    @test fixed_side_metadata.retained_units[3].primitive_family ==
          :projected_q_shell
    @test fixed_side_metadata.retained_units[3].representation_kind ==
          :shell_realized_pqs_fixture
    @test fixed_side_metadata.retained_units[3].shell_realized_pqs_metadata_oracle_fixture
    @test !fixed_side_metadata.retained_units[3].shell_realized_pqs_source_box_operator_ready
    @test !fixed_side_metadata.retained_units[3].compact_source_space_transform_available
    @test fixed_side_metadata.retained_units[3].shell_row_oracle_only
    @test fixed_side_metadata.retained_units[3].support_local_oracle_used
    @test !fixed_side_metadata.retained_units[3].ida_weight_division_allowed
    @test fixed_side_metadata.absences_by_contract.packet_fixed_block_qw_hamiltonian_adoption
    @test fixed_side_metadata.absences_by_contract.retained_weight_ida_division
    @test fixed_side_metadata.absences_by_contract.shell_label_reconstruction_from_centers
    @test fixed_side_metadata.diagnostics.fixed_side_records_are_explicit_metadata
    @test fixed_side_metadata.diagnostics.coverage_complete
    @test fixed_side_metadata.diagnostics.shell_realized_pqs_fixture_count == 1
    @test fixed_side_metadata.diagnostics.shell_realized_pqs_source_box_operator_ready_count ==
          0
    @test fixed_side_metadata.diagnostics.shell_realized_pqs_fixtures_are_metadata_oracle_only
    @test !fixed_side_metadata.diagnostics.route_construction_changed
    @test !fixed_side_metadata.diagnostics.packet_adoption
    @test !fixed_side_metadata.diagnostics.fixed_block_routing
    @test !fixed_side_metadata.diagnostics.qwhamiltonian_changed
    @test !fixed_side_metadata.diagnostics.hamiltonian_matrix_built
    @test !fixed_side_metadata.diagnostics.public_default_consumes
    @test !fixed_side_metadata.diagnostics.mwg_ida_semantics_changed
    @test !fixed_side_metadata.diagnostics.ecp_terms_implemented
    @test !fixed_side_metadata.diagnostics.scf_hf_validation_claim
    @test !fixed_side_metadata.diagnostics.cr2_science_status_changed
    @test !fixed_side_metadata.diagnostics.retained_weight_or_ida_division

    fixed_column_label_inventory =
        metrics_module._pqs_current_route_fixed_column_label_inventory(
            _synthetic_fixed_side_inventory();
            provenance = (source = :synthetic_fixed_column_label_inventory_test,),
        )
    @test fixed_column_label_inventory.object_kind ==
          :pqs_current_route_fixed_column_label_inventory
    @test fixed_column_label_inventory.status ==
          :private_fixed_column_label_inventory
    @test fixed_column_label_inventory.schema_version ==
          :pqs_fixed_column_labels_private_v1
    @test fixed_column_label_inventory.fixed_dimension == 15
    @test fixed_column_label_inventory.row_count == 15
    @test fixed_column_label_inventory.fixed_cols == collect(1:15)
    @test fixed_column_label_inventory.unit_indices ==
          [fill(1, 4); fill(2, 4); fill(3, 7)]
    @test fixed_column_label_inventory.unit_labels ==
          [
              fill(:outer_mismatch_z_low_slab, 4);
              fill(:left_atom_box, 4);
              fill(:regular_shared_molecular_shell_1, 7)
          ]
    @test fixed_column_label_inventory.unit_categories ==
          [
              fill(:product_doside, 4);
              fill(:support_dense, 4);
              fill(:shell_realized_pqs_fixture, 7)
          ]
    @test fixed_column_label_inventory.unit_kinds ==
          [
              fill(:product_doside, 4);
              fill(:atom_core_cube, 4);
              fill(:projected_q_shell, 7)
          ]
    @test fixed_column_label_inventory.unit_retained_starts ==
          [fill(1, 4); fill(5, 4); fill(9, 7)]
    @test fixed_column_label_inventory.unit_retained_stops ==
          [fill(4, 4); fill(8, 4); fill(15, 7)]
    @test fixed_column_label_inventory.source_region_labels ==
          fixed_column_label_inventory.unit_labels
    @test all(
        ==(:retained_unit_region_label),
        fixed_column_label_inventory.source_region_label_statuses,
    )
    @test all(isnothing, fixed_column_label_inventory.source_box_labels)
    @test all(==(:unavailable), fixed_column_label_inventory.source_box_label_statuses)
    @test all(isnothing, fixed_column_label_inventory.owner_labels)
    @test all(==(:unavailable), fixed_column_label_inventory.owner_label_statuses)
    @test all(==(:unavailable), fixed_column_label_inventory.shell_label_statuses)
    @test all(==(0), fixed_column_label_inventory.shell_indices)
    @test all(==(:unavailable), fixed_column_label_inventory.ray_label_statuses)
    @test all(isnothing, fixed_column_label_inventory.ray_ids)
    @test all(isnothing, fixed_column_label_inventory.ray_family_labels)
    @test all(==(:unavailable), fixed_column_label_inventory.radial_order_statuses)
    @test all(==(0), fixed_column_label_inventory.radial_orders)
    @test !any(fixed_column_label_inventory.inferred_from_centers)
    @test !any(fixed_column_label_inventory.inferred_from_nearest_grid)
    @test !any(fixed_column_label_inventory.inferred_from_support_order)
    @test !any(fixed_column_label_inventory.inferred_from_support_indices)
    @test !any(fixed_column_label_inventory.inferred_from_raw_to_final_support)
    @test fixed_column_label_inventory.label_status.source_region_label_status ==
          :retained_unit_region_label
    @test fixed_column_label_inventory.label_status.source_box_label_status ==
          :unavailable
    @test fixed_column_label_inventory.label_status.owner_label_status ==
          :unavailable
    @test fixed_column_label_inventory.label_status.shell_label_status ==
          :unavailable
    @test fixed_column_label_inventory.label_status.ray_label_status ==
          :unavailable
    @test fixed_column_label_inventory.label_status.radial_order_status ==
          :unavailable
    @test fixed_column_label_inventory.absences_by_contract.source_box_labels
    @test fixed_column_label_inventory.absences_by_contract.owner_labels
    @test fixed_column_label_inventory.absences_by_contract.shell_labels
    @test fixed_column_label_inventory.absences_by_contract.ray_labels
    @test fixed_column_label_inventory.absences_by_contract.radial_order_labels
    @test fixed_column_label_inventory.absences_by_contract.coordinate_or_nearest_grid_reconstruction
    @test fixed_column_label_inventory.absences_by_contract.support_row_order_or_support_index_inference
    @test fixed_column_label_inventory.absences_by_contract.raw_to_final_support_inference
    @test fixed_column_label_inventory.diagnostics.fixed_cols_cover_1_to_fixed_dimension
    @test fixed_column_label_inventory.diagnostics.unit_ranges_match_inventory
    @test fixed_column_label_inventory.diagnostics.source_region_labels_match_unit_labels
    @test fixed_column_label_inventory.diagnostics.product_doside_row_count == 4
    @test fixed_column_label_inventory.diagnostics.support_dense_row_count == 4
    @test fixed_column_label_inventory.diagnostics.shell_realized_pqs_row_count == 7
    @test fixed_column_label_inventory.diagnostics.shell_realized_pqs_shell_ray_radial_unavailable
    @test !fixed_column_label_inventory.diagnostics.inferred_from_centers
    @test !fixed_column_label_inventory.diagnostics.inferred_from_nearest_grid
    @test !fixed_column_label_inventory.diagnostics.inferred_from_support_order
    @test !fixed_column_label_inventory.diagnostics.inferred_from_support_indices
    @test !fixed_column_label_inventory.diagnostics.inferred_from_raw_to_final_support
    @test !fixed_column_label_inventory.diagnostics.retained_weight_or_ida_division
    @test !fixed_column_label_inventory.diagnostics.route_construction_changed
    @test !fixed_column_label_inventory.diagnostics.packet_adoption
    @test !fixed_column_label_inventory.diagnostics.qwhamiltonian_changed
    @test !fixed_column_label_inventory.diagnostics.public_default_consumes

    source_shell_mode_inventory =
        metrics_module._pqs_current_route_source_shell_mode_inventory(
            _synthetic_fixed_side_inventory();
            provenance = (source = :synthetic_source_shell_mode_inventory_test,),
        )
    @test source_shell_mode_inventory.object_kind ==
          :pqs_current_route_source_shell_mode_inventory
    @test source_shell_mode_inventory.status ==
          :native_current_route_source_shell_modes
    @test source_shell_mode_inventory.schema_version ==
          :pqs_source_shell_modes_private_v1
    @test source_shell_mode_inventory.fixed_dimension == 15
    @test source_shell_mode_inventory.source_shell_count == 3
    @test source_shell_mode_inventory.source_mode_count == 38
    @test source_shell_mode_inventory.source_shells.source_shell_ids == [1, 2, 3]
    @test source_shell_mode_inventory.source_shells.unit_indices == [1, 2, 3]
    @test source_shell_mode_inventory.source_shells.unit_labels ==
          [
              :outer_mismatch_z_low_slab,
              :left_atom_box,
              :regular_shared_molecular_shell_1,
          ]
    @test source_shell_mode_inventory.source_shells.unit_categories ==
          [:product_doside, :support_dense, :shell_realized_pqs_fixture]
    @test source_shell_mode_inventory.source_shells.unit_kinds ==
          [:product_doside, :atom_core_cube, :projected_q_shell]
    @test source_shell_mode_inventory.source_shells.retained_starts == [1, 5, 9]
    @test source_shell_mode_inventory.source_shells.retained_stops == [4, 8, 15]
    @test source_shell_mode_inventory.source_shells.source_shell_labels ==
          [
              :outer_mismatch_z_low_slab,
              :left_atom_box,
              :regular_shared_molecular_shell_1,
          ]
    @test source_shell_mode_inventory.source_shells.source_shell_statuses ==
          [
              :native_product_doside_source_box,
              :native_support_dense_source_support_states,
              :native_shell_realized_boundary_source_box,
          ]
    @test source_shell_mode_inventory.source_shells.construction_kinds ==
          [
              :product_doside,
              :support_dense_direct_support,
              :shell_realized_pqs_fixture,
          ]
    @test source_shell_mode_inventory.source_shells.axis_kinds ==
          [
              :fixed :active :active
              :parent_lattice_support_state :parent_lattice_support_state :parent_lattice_support_state
              :raw_box_axis :raw_box_axis :raw_box_axis
          ]
    @test source_shell_mode_inventory.source_shells.axis_starts ==
          [
              1 1 3
              4 5 2
              10 20 30
          ]
    @test source_shell_mode_inventory.source_shells.axis_stops ==
          [
              1 2 4
              6 7 4
              24 34 54
          ]
    @test source_shell_mode_inventory.source_shells.axis_starts[3, :] == [10, 20, 30]
    @test source_shell_mode_inventory.source_shells.axis_stops[3, :] == [24, 34, 54]
    @test source_shell_mode_inventory.source_shells.fixed_axis_indices ==
          [
              1 0 0
              0 0 0
              0 0 0
          ]
    @test source_shell_mode_inventory.source_shells.contracted_dims ==
          [
              1 2 2
              3 3 3
              5 5 5
          ]
    @test source_shell_mode_inventory.source_shells.source_mode_counts == [4, 27, 7]
    @test source_shell_mode_inventory.source_shells.source_mode_orderings ==
          [
              :axis_function_indices_order,
              :construction_support_state_order,
              :boundary_mode_indices_order,
          ]
    @test all(==(:unavailable), source_shell_mode_inventory.source_shells.center_definitions)
    @test all(==(:unavailable), source_shell_mode_inventory.source_shells.center_statuses)
    @test !any(source_shell_mode_inventory.source_shells.lowdin_correction_applied)
    @test all(==(:unavailable), source_shell_mode_inventory.source_shells.shell_label_statuses)
    @test all(==(:unavailable), source_shell_mode_inventory.source_shells.ray_label_statuses)
    @test all(==(:unavailable), source_shell_mode_inventory.source_shells.radial_order_statuses)
    @test source_shell_mode_inventory.source_modes.source_shell_ids ==
          [fill(1, 4); fill(2, 27); fill(3, 7)]
    @test source_shell_mode_inventory.source_modes.mode_indices ==
          [collect(1:4); collect(1:27); collect(1:7)]
    @test source_shell_mode_inventory.source_modes.unit_labels ==
          [
              fill(:outer_mismatch_z_low_slab, 4);
              fill(:left_atom_box, 27);
              fill(:regular_shared_molecular_shell_1, 7)
          ]
    @test source_shell_mode_inventory.source_modes.native_source_id_labels[1:4] == [
        "source_mode:1:1,1,1",
        "source_mode:1:1,2,1",
        "source_mode:1:1,1,2",
        "source_mode:1:1,2,2",
    ]
    @test source_shell_mode_inventory.source_modes.native_source_id_labels[5] ==
          "source_mode:2:1,1,1"
    @test source_shell_mode_inventory.source_modes.native_source_id_labels[31] ==
          "source_mode:2:3,3,3"
    @test source_shell_mode_inventory.source_modes.native_source_id_labels[32:38] == [
        "source_mode:3:1,1,1",
        "source_mode:3:1,1,5",
        "source_mode:3:1,5,1",
        "source_mode:3:5,1,1",
        "source_mode:3:5,5,1",
        "source_mode:3:5,1,5",
        "source_mode:3:1,5,5",
    ]
    @test source_shell_mode_inventory.source_modes.local_axis_function_indices[1:4, :] == [
        1 1 1
        1 2 1
        1 1 2
        1 2 2
    ]
    @test source_shell_mode_inventory.source_modes.local_axis_function_indices[5, :] ==
          [1, 1, 1]
    @test source_shell_mode_inventory.source_modes.local_axis_function_indices[31, :] ==
          [3, 3, 3]
    @test source_shell_mode_inventory.source_modes.local_axis_function_indices[32:38, :] == [
        1 1 1
        1 1 5
        1 5 1
        5 1 1
        5 5 1
        5 1 5
        1 5 5
    ]
    @test all(1:source_shell_mode_inventory.source_mode_count) do row
        shell_id = source_shell_mode_inventory.source_modes.source_shell_ids[row]
        all(1:3) do axis
            local_axis =
                source_shell_mode_inventory.source_modes.local_axis_function_indices[
                    row,
                    axis,
                ]
            local_axis >= 1 &&
                local_axis <= source_shell_mode_inventory.source_shells.contracted_dims[
                    shell_id,
                    axis,
                ]
        end
    end
    @test source_shell_mode_inventory.source_modes.source_axis_indices[1:4, :] == [
        1 1 3
        1 2 3
        1 1 4
        1 2 4
    ]
    @test source_shell_mode_inventory.source_modes.source_axis_indices[5, :] ==
          [4, 5, 2]
    @test source_shell_mode_inventory.source_modes.source_axis_indices[31, :] ==
          [6, 7, 4]
    @test source_shell_mode_inventory.source_modes.source_axis_indices[32:38, :] ==
          source_shell_mode_inventory.source_modes.local_axis_function_indices[32:38, :]
    @test source_shell_mode_inventory.source_modes.source_axis_indices[32, :] ==
          [1, 1, 1]
    @test source_shell_mode_inventory.source_modes.source_axis_indices[33, :] ==
          [1, 1, 5]
    @test source_shell_mode_inventory.source_modes.source_axis_indices[32, :] !=
          source_shell_mode_inventory.source_shells.axis_starts[3, :]
    @test source_shell_mode_inventory.source_modes.parent_lattice_axis_indices[1:4, :] ==
          source_shell_mode_inventory.source_modes.source_axis_indices[1:4, :]
    @test source_shell_mode_inventory.source_modes.parent_lattice_axis_indices[5, :] ==
          [4, 5, 2]
    @test source_shell_mode_inventory.source_modes.parent_lattice_axis_indices[31, :] ==
          [6, 7, 4]
    @test all(
        ==(0),
        source_shell_mode_inventory.source_modes.parent_lattice_axis_indices[32:38, :],
    )
    @test all(
        ==(:native_product_doside_source_mode),
        source_shell_mode_inventory.source_modes.source_mode_statuses[1:4],
    )
    @test all(
        ==(:native_support_dense_source_support_state),
        source_shell_mode_inventory.source_modes.source_mode_statuses[5:31],
    )
    @test all(
        ==(:native_shell_realized_boundary_source_mode),
        source_shell_mode_inventory.source_modes.source_mode_statuses[32:38],
    )
    @test all(
        ==(:native_product_axis_tuple),
        source_shell_mode_inventory.source_modes.source_axis_tuple_statuses[1:4],
    )
    @test all(
        ==(:native_parent_lattice_support_state),
        source_shell_mode_inventory.source_modes.source_axis_tuple_statuses[5:31],
    )
    @test all(
        ==(:native_boundary_source_mode_tuple),
        source_shell_mode_inventory.source_modes.source_axis_tuple_statuses[32:38],
    )
    @test all(
        ==(:native_product_parent_lattice_axis_tuple),
        source_shell_mode_inventory.source_modes.parent_lattice_axis_statuses[1:4],
    )
    @test all(
        ==(:native_parent_lattice_support_state),
        source_shell_mode_inventory.source_modes.parent_lattice_axis_statuses[5:31],
    )
    @test all(
        ==(:unavailable),
        source_shell_mode_inventory.source_modes.parent_lattice_axis_statuses[32:38],
    )
    @test all(isnan, source_shell_mode_inventory.source_modes.center_coordinates)
    @test all(==(:unavailable), source_shell_mode_inventory.source_modes.center_definitions)
    @test all(==(:unavailable), source_shell_mode_inventory.source_modes.center_statuses)
    @test !any(source_shell_mode_inventory.source_modes.lowdin_correction_applied)
    @test all(==(:unavailable), source_shell_mode_inventory.source_modes.shell_label_statuses)
    @test all(==(:unavailable), source_shell_mode_inventory.source_modes.ray_label_statuses)
    @test all(==(:unavailable), source_shell_mode_inventory.source_modes.radial_order_statuses)
    @test !any(source_shell_mode_inventory.source_modes.inferred_from_centers)
    @test !any(source_shell_mode_inventory.source_modes.inferred_from_nearest_grid)
    @test !any(source_shell_mode_inventory.source_modes.inferred_from_support_order)
    @test !any(source_shell_mode_inventory.source_modes.inferred_from_support_indices)
    @test !any(source_shell_mode_inventory.source_modes.inferred_from_raw_to_final_support)
    @test source_shell_mode_inventory.center_status ==
          :unavailable_missing_native_comx_center_facts
    @test source_shell_mode_inventory.covered_unit_categories ==
          (:product_doside, :support_dense, :shell_realized_pqs_fixture)
    @test source_shell_mode_inventory.non_product_source_mode_status ==
          :native_non_product_source_shell_mode_labels
    @test source_shell_mode_inventory.source_mode_label_status ==
          :native_source_mode_tuple_labels_shell_ray_radial_unavailable
    @test source_shell_mode_inventory.absences_by_contract.repo_ray_grouping_policy
    @test source_shell_mode_inventory.absences_by_contract.product_axis_tuples_not_interpreted_as_ray_labels
    @test source_shell_mode_inventory.absences_by_contract.representative_centers_as_identity_labels
    @test source_shell_mode_inventory.absences_by_contract.native_comx_centers
    @test !source_shell_mode_inventory.absences_by_contract.support_dense_source_shell_modes
    @test source_shell_mode_inventory.absences_by_contract.shell_realized_pqs_source_relations
    @test source_shell_mode_inventory.absences_by_contract.lowdin_mixture_weights_or_spans
    @test source_shell_mode_inventory.diagnostics.source_shell_count == 3
    @test source_shell_mode_inventory.diagnostics.source_mode_count == 38
    @test source_shell_mode_inventory.diagnostics.product_doside_source_shell_count == 1
    @test source_shell_mode_inventory.diagnostics.product_doside_source_mode_count == 4
    @test source_shell_mode_inventory.diagnostics.support_dense_source_shell_count == 1
    @test source_shell_mode_inventory.diagnostics.support_dense_source_mode_count == 27
    @test source_shell_mode_inventory.diagnostics.shell_realized_pqs_source_shell_count ==
          1
    @test source_shell_mode_inventory.diagnostics.shell_realized_pqs_source_mode_count ==
          7
    @test source_shell_mode_inventory.diagnostics.support_dense_unavailable_unit_count == 0
    @test source_shell_mode_inventory.diagnostics.support_dense_unavailable_column_count == 0
    @test source_shell_mode_inventory.diagnostics.shell_realized_pqs_unavailable_unit_count == 0
    @test source_shell_mode_inventory.diagnostics.shell_realized_pqs_unavailable_column_count ==
          0
    @test source_shell_mode_inventory.diagnostics.total_unavailable_unit_count == 0
    @test source_shell_mode_inventory.diagnostics.total_unavailable_column_count == 0
    @test !source_shell_mode_inventory.diagnostics.product_doside_source_shell_modes_only
    @test source_shell_mode_inventory.diagnostics.repo_exports_native_facts_not_ray_policy
    @test source_shell_mode_inventory.diagnostics.product_axis_tuples_not_interpreted_as_ray_labels
    @test !source_shell_mode_inventory.diagnostics.product_axis_tuples_are_ray_labels
    @test !source_shell_mode_inventory.diagnostics.representative_centers_are_identity_labels
    @test !source_shell_mode_inventory.diagnostics.inferred_from_centers
    @test !source_shell_mode_inventory.diagnostics.inferred_from_nearest_grid
    @test !source_shell_mode_inventory.diagnostics.inferred_from_support_order
    @test !source_shell_mode_inventory.diagnostics.inferred_from_support_indices
    @test !source_shell_mode_inventory.diagnostics.inferred_from_raw_to_final_support
    @test !source_shell_mode_inventory.diagnostics.retained_weight_or_ida_division
    @test !source_shell_mode_inventory.diagnostics.route_construction_changed
    @test !source_shell_mode_inventory.diagnostics.packet_adoption
    @test !source_shell_mode_inventory.diagnostics.qwhamiltonian_changed
    @test !source_shell_mode_inventory.diagnostics.public_default_consumes

    source_metadata_contract =
        metrics_module._pqs_source_metadata_export_contract()
    @test source_metadata_contract.schema_version ==
          :pqs_source_shell_modes_private_v1
    @test source_metadata_contract.label_semantics ==
          :construction_native_identifiers_not_relations
    @test source_metadata_contract.source_mode_local_axis_semantics ==
          :normalized_source_shell_local_coordinates
    @test source_metadata_contract.parent_lattice_axis_coordinate_policy ==
          :explicit_columns_when_native_available
    @test source_metadata_contract.shell_realized_pqs_source_axis_indices ==
          :local_native_source_mode_coordinates
    @test source_metadata_contract.repo_ray_id_policy == :not_exported
    @test source_metadata_contract.relation_weight_span_policy ==
          :not_in_source_metadata_tables
    @test !("ray_id" in source_metadata_contract.source_shells_header)
    @test !("ray_id" in source_metadata_contract.source_modes_header)

    source_shell_table_io = IOBuffer()
    metrics_module._write_pqs_source_shells_table(
        source_shell_table_io,
        source_shell_mode_inventory,
    )
    source_shell_table_lines =
        split(chomp(String(take!(source_shell_table_io))), '\n')
    @test length(source_shell_table_lines) ==
          source_shell_mode_inventory.source_shell_count + 1
    @test split(source_shell_table_lines[1], '\t') ==
          collect(source_metadata_contract.source_shells_header)
    @test !occursin("ray_id", source_shell_table_lines[1])
    shell_realized_source_shell_fields =
        split(source_shell_table_lines[4], '\t')
    @test shell_realized_source_shell_fields[4] == "shell_realized_pqs_fixture"
    @test shell_realized_source_shell_fields[11:13] ==
          ["raw_box_axis", "raw_box_axis", "raw_box_axis"]
    @test shell_realized_source_shell_fields[14:16] == ["10", "20", "30"]
    @test shell_realized_source_shell_fields[17:19] == ["24", "34", "54"]
    @test shell_realized_source_shell_fields[23:25] == ["5", "5", "5"]
    @test shell_realized_source_shell_fields[31:33] ==
          ["unavailable", "unavailable", "unavailable"]

    source_mode_table_io = IOBuffer()
    metrics_module._write_pqs_source_modes_table(
        source_mode_table_io,
        source_shell_mode_inventory,
    )
    source_mode_table_lines =
        split(chomp(String(take!(source_mode_table_io))), '\n')
    @test length(source_mode_table_lines) ==
          source_shell_mode_inventory.source_mode_count + 1
    @test split(source_mode_table_lines[1], '\t') ==
          collect(source_metadata_contract.source_modes_header)
    @test !occursin("ray_id", source_mode_table_lines[1])
    @test occursin("parent_lattice_axis_x", source_mode_table_lines[1])
    shell_realized_source_mode_fields =
        split(source_mode_table_lines[33], '\t')
    @test shell_realized_source_mode_fields[1:4] ==
          ["3", "1", "regular_shared_molecular_shell_1", "source_mode:3:1,1,1"]
    @test shell_realized_source_mode_fields[5:7] == ["1", "1", "1"]
    @test shell_realized_source_mode_fields[8:10] == ["1", "1", "1"]
    @test shell_realized_source_mode_fields[11:13] == ["0", "0", "0"]
    @test shell_realized_source_mode_fields[8:13] != ["10", "20", "30"]
    @test shell_realized_source_mode_fields[14] ==
          "native_shell_realized_boundary_source_mode"
    @test shell_realized_source_mode_fields[15] ==
          "native_boundary_source_mode_tuple"
    @test shell_realized_source_mode_fields[16] == "unavailable"
    @test shell_realized_source_mode_fields[23:25] ==
          ["unavailable", "unavailable", "unavailable"]
    @test shell_realized_source_mode_fields[26:30] ==
          ["false", "false", "false", "false", "false"]

    centered_source_shell_mode_inventory =
        metrics_module._pqs_current_route_source_shell_mode_inventory(
            _synthetic_fixed_side_inventory(;
                product_source_axis_center_vectors = (
                    [10.0],
                    [20.0, 21.0],
                    [30.0, 31.0],
                ),
            );
            provenance = (source = :synthetic_native_source_center_test,),
        )
    @test centered_source_shell_mode_inventory.center_status ==
          :partial_native_representative_product_doside
    @test centered_source_shell_mode_inventory.source_shells.center_definitions ==
          [:comx_construction, :unavailable, :unavailable]
    @test centered_source_shell_mode_inventory.source_shells.center_statuses ==
          [:native_representative, :unavailable, :unavailable]
    @test centered_source_shell_mode_inventory.source_modes.center_coordinates[1:4, :] == [
        10.0 20.0 30.0
        10.0 21.0 30.0
        10.0 20.0 31.0
        10.0 21.0 31.0
    ]
    @test all(
        isnan,
        centered_source_shell_mode_inventory.source_modes.center_coordinates[5:38, :],
    )
    @test all(
        ==(:comx_construction),
        centered_source_shell_mode_inventory.source_modes.center_definitions[1:4],
    )
    @test all(
        ==(:unavailable),
        centered_source_shell_mode_inventory.source_modes.center_definitions[5:38],
    )
    @test all(
        ==(:native_representative),
        centered_source_shell_mode_inventory.source_modes.center_statuses[1:4],
    )
    @test all(
        ==(:unavailable),
        centered_source_shell_mode_inventory.source_modes.center_statuses[5:38],
    )
    @test centered_source_shell_mode_inventory.diagnostics.center_status ==
          :partial_native_representative_product_doside
    @test centered_source_shell_mode_inventory.diagnostics.center_definition ==
          :comx_construction
    @test centered_source_shell_mode_inventory.diagnostics.native_center_shell_count == 1
    @test centered_source_shell_mode_inventory.diagnostics.native_center_mode_count == 4
    @test centered_source_shell_mode_inventory.absences_by_contract.native_comx_centers
    @test !any(centered_source_shell_mode_inventory.source_modes.inferred_from_centers)
    @test !centered_source_shell_mode_inventory.diagnostics.retained_weight_or_ida_division
    @test_throws DimensionMismatch metrics_module._pqs_current_route_source_shell_mode_inventory(
        _synthetic_fixed_side_inventory(;
            product_source_axis_center_vectors = ([10.0], [20.0], [30.0, 31.0]),
        ),
    )

    source_relation_inventory =
        metrics_module._pqs_current_route_fixed_column_source_relation_inventory(
            _synthetic_fixed_side_inventory();
            provenance = (source = :synthetic_source_relation_inventory_test,),
        )
    @test source_relation_inventory.object_kind ==
          :pqs_current_route_fixed_column_source_relation_inventory
    @test source_relation_inventory.status ==
          :product_doside_axis_tuple_relations_only
    @test source_relation_inventory.schema_version ==
          :pqs_fixed_column_source_relations_private_v1
    @test source_relation_inventory.fixed_dimension == 15
    @test source_relation_inventory.row_count == 4
    @test source_relation_inventory.relation_rows_available
    @test source_relation_inventory.fixed_cols == [1, 2, 3, 4]
    @test source_relation_inventory.relation_indices == ones(Int, 4)
    @test source_relation_inventory.relation_kinds ==
          fill(:product_axis_tuple, 4)
    @test source_relation_inventory.source_unit_labels ==
          fill(:outer_mismatch_z_low_slab, 4)
    @test source_relation_inventory.source_mode_labels == [
        "product_axis_tuple:1,1,3",
        "product_axis_tuple:1,2,3",
        "product_axis_tuple:1,1,4",
        "product_axis_tuple:1,2,4",
    ]
    @test source_relation_inventory.source_axis_indices == [
        1 1 3
        1 2 3
        1 1 4
        1 2 4
    ]
    @test source_relation_inventory.local_axis_function_indices == [
        1 1 1
        1 2 1
        1 1 2
        1 2 2
    ]
    @test all(==(:native_product_axis_tuple), source_relation_inventory.relation_statuses)
    @test all(==(:unavailable), source_relation_inventory.shell_label_statuses)
    @test all(==(:unavailable), source_relation_inventory.ray_label_statuses)
    @test all(==(:unavailable), source_relation_inventory.radial_order_statuses)
    @test all(==(:unavailable), source_relation_inventory.coefficient_statuses)
    @test all(==(:unavailable), source_relation_inventory.weight_statuses)
    @test all(==(:unavailable), source_relation_inventory.span_statuses)
    @test !any(source_relation_inventory.inferred_from_centers)
    @test !any(source_relation_inventory.inferred_from_nearest_grid)
    @test !any(source_relation_inventory.inferred_from_support_order)
    @test !any(source_relation_inventory.inferred_from_support_indices)
    @test !any(source_relation_inventory.inferred_from_raw_to_final_support)
    @test source_relation_inventory.covered_unit_categories ==
          (:product_doside,)
    @test source_relation_inventory.non_product_relation_status ==
          :unavailable_missing_native_non_product_relation_producer
    @test source_relation_inventory.relation_label_status ==
          :native_product_axis_tuple_only_shell_ray_radial_unavailable
    @test source_relation_inventory.relation_weight_status == :unavailable
    @test source_relation_inventory.relation_span_status == :unavailable
    @test source_relation_inventory.missing_producer ==
          :construction_native_non_product_shell_ray_relation_producer
    @test source_relation_inventory.absences_by_contract.product_axis_tuples_not_interpreted_as_ray_labels
    @test source_relation_inventory.absences_by_contract.support_dense_relation_rows
    @test source_relation_inventory.absences_by_contract.shell_realized_pqs_relation_rows
    @test source_relation_inventory.absences_by_contract.relation_weights_or_spans
    @test source_relation_inventory.diagnostics.product_doside_row_count == 4
    @test source_relation_inventory.diagnostics.support_dense_unavailable_column_count == 4
    @test source_relation_inventory.diagnostics.shell_realized_pqs_unavailable_column_count ==
          7
    @test source_relation_inventory.diagnostics.total_unavailable_column_count == 11
    @test source_relation_inventory.diagnostics.product_axis_tuples_not_interpreted_as_ray_labels
    @test !source_relation_inventory.diagnostics.product_axis_tuples_are_ray_labels
    @test !source_relation_inventory.diagnostics.inferred_from_centers
    @test !source_relation_inventory.diagnostics.inferred_from_nearest_grid
    @test !source_relation_inventory.diagnostics.inferred_from_support_order
    @test !source_relation_inventory.diagnostics.inferred_from_support_indices
    @test !source_relation_inventory.diagnostics.inferred_from_raw_to_final_support
    @test !source_relation_inventory.diagnostics.retained_weight_or_ida_division
    @test !source_relation_inventory.diagnostics.route_construction_changed
    @test !source_relation_inventory.diagnostics.packet_adoption
    @test !source_relation_inventory.diagnostics.qwhamiltonian_changed
    @test !source_relation_inventory.diagnostics.public_default_consumes

    source_shell_report =
        metrics_module._pqs_pqs_product_component_route_smoke_report_adapter(
            summaries,
            _synthetic_mwg_facts(;
                overrides = Dict(
                    "fixed_fixed_shape" => "(15, 15)",
                    "fixed_residual_shape" => "(15, 6)",
                    "final_interaction_shape" => "(21, 21)",
                ),
            ),
        )
    sidecar_with_source_shell_modes =
        metrics_module._pqs_pqs_product_component_route_smoke_cr2_sidecar_schema(
            source_shell_report;
            source_shell_mode_inventory = source_shell_mode_inventory,
            provenance = (source = :synthetic_sidecar_source_shell_mode_test,),
        )
    @test sidecar_with_source_shell_modes.source_shell_mode_inventory ===
          source_shell_mode_inventory
    @test sidecar_with_source_shell_modes.diagnostics.source_shell_mode_inventory_available
    @test sidecar_with_source_shell_modes.diagnostics.source_shell_count == 3
    @test sidecar_with_source_shell_modes.diagnostics.source_mode_count == 38
    @test sidecar_with_source_shell_modes.diagnostics.source_shell_mode_center_status ==
          :unavailable_missing_native_comx_center_facts
    @test !sidecar_with_source_shell_modes.diagnostics.source_shell_mode_product_doside_only
    @test sidecar_with_source_shell_modes.diagnostics.source_shell_mode_support_dense_unavailable_column_count ==
          0
    @test sidecar_with_source_shell_modes.diagnostics.source_shell_mode_shell_realized_pqs_unavailable_column_count ==
          0
    @test sidecar_with_source_shell_modes.source_shell_mode_inventory.source_shells.unit_categories ==
          [:product_doside, :support_dense, :shell_realized_pqs_fixture]
    @test sidecar_with_source_shell_modes.source_shell_mode_inventory.covered_unit_categories ==
          (:product_doside, :support_dense, :shell_realized_pqs_fixture)
    @test sidecar_with_source_shell_modes.source_shell_mode_inventory.diagnostics.product_doside_source_mode_count ==
          4
    @test sidecar_with_source_shell_modes.source_shell_mode_inventory.diagnostics.support_dense_source_mode_count ==
          27
    @test sidecar_with_source_shell_modes.source_shell_mode_inventory.diagnostics.shell_realized_pqs_source_mode_count ==
          7
    @test sidecar_with_source_shell_modes.source_shell_mode_inventory.diagnostics.total_unavailable_column_count ==
          0

    source_shell_mode_io = IOBuffer()
    metrics_module._write_pqs_pqs_product_component_route_smoke_cr2_sidecar_schema_report(
        source_shell_mode_io,
        sidecar_with_source_shell_modes,
    )
    source_shell_mode_text = String(take!(source_shell_mode_io))
    @test occursin("[source_shell_mode_inventory]", source_shell_mode_text)
    @test occursin("schema_version\tpqs_source_shell_modes_private_v1", source_shell_mode_text)
    @test occursin("status\tnative_current_route_source_shell_modes", source_shell_mode_text)
    @test occursin("source_shell_count\t3", source_shell_mode_text)
    @test occursin("source_mode_count\t38", source_shell_mode_text)
    @test occursin(
        "center_status\tunavailable_missing_native_comx_center_facts",
        source_shell_mode_text,
    )
    @test occursin(
        "covered_unit_categories\t(:product_doside, :support_dense, :shell_realized_pqs_fixture)",
        source_shell_mode_text,
    )
    @test occursin(
        "non_product_source_mode_status\tnative_non_product_source_shell_mode_labels",
        source_shell_mode_text,
    )
    @test occursin(
        "source_mode_label_status\tnative_source_mode_tuple_labels_shell_ray_radial_unavailable",
        source_shell_mode_text,
    )
    @test occursin("product_doside_only\tfalse", source_shell_mode_text)
    @test occursin("product_doside_source_shell_count\t1", source_shell_mode_text)
    @test occursin("product_doside_source_mode_count\t4", source_shell_mode_text)
    @test occursin("support_dense_source_shell_count\t1", source_shell_mode_text)
    @test occursin("support_dense_source_mode_count\t27", source_shell_mode_text)
    @test occursin("support_dense_unavailable_column_count\t0", source_shell_mode_text)
    @test occursin("shell_realized_pqs_source_shell_count\t1", source_shell_mode_text)
    @test occursin("shell_realized_pqs_source_mode_count\t7", source_shell_mode_text)
    @test occursin(
        "shell_realized_pqs_unavailable_column_count\t0",
        source_shell_mode_text,
    )
    @test occursin("total_unavailable_unit_count\t0", source_shell_mode_text)
    @test occursin("total_unavailable_column_count\t0", source_shell_mode_text)
    @test occursin("inferred_from_centers\tfalse", source_shell_mode_text)
    @test occursin("inferred_from_nearest_grid\tfalse", source_shell_mode_text)
    @test occursin("inferred_from_support_order\tfalse", source_shell_mode_text)
    @test occursin("inferred_from_support_indices\tfalse", source_shell_mode_text)
    @test occursin("inferred_from_raw_to_final_support\tfalse", source_shell_mode_text)

    @test_throws ArgumentError metrics_module._pqs_pqs_product_component_route_smoke_cr2_sidecar_schema(
        source_shell_report;
        source_shell_mode_inventory =
            merge(source_shell_mode_inventory, (object_kind = :not_source_shell_modes,)),
    )
    @test_throws DimensionMismatch metrics_module._pqs_pqs_product_component_route_smoke_cr2_sidecar_schema(
        source_shell_report;
        source_shell_mode_inventory =
            merge(source_shell_mode_inventory, (fixed_dimension = 14,)),
    )
    @test_throws ArgumentError metrics_module._pqs_pqs_product_component_route_smoke_cr2_sidecar_schema(
        source_shell_report;
        source_shell_mode_inventory = merge(
            source_shell_mode_inventory,
            (
                diagnostics = merge(
                    source_shell_mode_inventory.diagnostics,
                    (inferred_from_centers = true,),
                ),
            ),
        ),
    )

    sidecar_with_fixed =
        metrics_module._pqs_pqs_product_component_route_smoke_cr2_sidecar_schema(
            report;
            fixed_side_retained_unit_metadata = fixed_side_metadata,
            provenance = (source = :synthetic_sidecar_fixed_metadata_test,),
        )
    @test sidecar_with_fixed.fixed_side_retained_unit_metadata ===
          fixed_side_metadata
    @test sidecar_with_fixed.diagnostics.fixed_side_retained_unit_metadata_available
    @test sidecar_with_fixed.fixed_side_retained_unit_metadata.labels.shell_label_status ==
          :unavailable
    @test !sidecar_with_fixed.fixed_side_retained_unit_metadata.labels.label_reconstruction_from_centers
    @test !sidecar_with_fixed.fixed_side_retained_unit_metadata.labels.nearest_grid_or_center_label_heuristic
    @test !sidecar_with_fixed.fixed_side_retained_unit_metadata.diagnostics.retained_weight_or_ida_division

    sidecar_with_fixed_io = IOBuffer()
    metrics_module._write_pqs_pqs_product_component_route_smoke_cr2_sidecar_schema_report(
        sidecar_with_fixed_io,
        sidecar_with_fixed,
    )
    sidecar_with_fixed_text = String(take!(sidecar_with_fixed_io))
    @test occursin("[fixed_side_retained_unit_metadata]", sidecar_with_fixed_text)
    @test occursin("fixed_dimension\t15", sidecar_with_fixed_text)
    @test occursin("unit_count\t3", sidecar_with_fixed_text)
    @test occursin(
        "source_unit_label_status\texplicit_inventory_unit_keys",
        sidecar_with_fixed_text,
    )
    @test occursin(
        "source_unit_labels\t(:outer_mismatch_z_low_slab, :left_atom_box, :regular_shared_molecular_shell_1)",
        sidecar_with_fixed_text,
    )
    @test occursin("shell_label_status\tunavailable", sidecar_with_fixed_text)
    @test occursin(
        "label_reconstruction_from_centers\tfalse",
        sidecar_with_fixed_text,
    )
    @test occursin(
        "nearest_grid_or_center_label_heuristic\tfalse",
        sidecar_with_fixed_text,
    )
    @test occursin(
        "unit.regular_shared_molecular_shell_1.retained_range\t9:15",
        sidecar_with_fixed_text,
    )
    @test occursin(
        "unit.regular_shared_molecular_shell_1.shell_realized_pqs_metadata_oracle_fixture\ttrue",
        sidecar_with_fixed_text,
    )
    @test occursin(
        "unit.regular_shared_molecular_shell_1.shell_realized_pqs_source_box_operator_ready\tfalse",
        sidecar_with_fixed_text,
    )

    bad_inventory = merge(
        _synthetic_fixed_side_inventory(),
        (object_kind = :not_current_route_inventory,),
    )
    @test_throws ArgumentError metrics_module._pqs_current_route_fixed_side_retained_unit_metadata(
        bad_inventory,
    )
    @test_throws ArgumentError metrics_module._pqs_current_route_fixed_column_label_inventory(
        bad_inventory,
    )
    @test_throws ArgumentError metrics_module._pqs_current_route_source_shell_mode_inventory(
        bad_inventory,
    )
    @test_throws ArgumentError metrics_module._pqs_current_route_fixed_column_source_relation_inventory(
        bad_inventory,
    )
    bad_coverage = merge(
        _synthetic_fixed_side_inventory(),
        (
            coverage = merge(
                _synthetic_fixed_side_inventory().coverage,
                (covers_every_column_once = false,),
            ),
        ),
    )
    @test_throws ArgumentError metrics_module._pqs_current_route_fixed_side_retained_unit_metadata(
        bad_coverage,
    )
    @test_throws ArgumentError metrics_module._pqs_current_route_fixed_column_label_inventory(
        bad_coverage,
    )
    @test_throws ArgumentError metrics_module._pqs_current_route_source_shell_mode_inventory(
        bad_coverage,
    )
    @test_throws ArgumentError metrics_module._pqs_current_route_fixed_column_source_relation_inventory(
        bad_coverage,
    )
    bad_ready_inventory = _synthetic_fixed_side_inventory(;
        shell_source_box_operator_ready = true,
    )
    @test_throws ArgumentError metrics_module._pqs_current_route_fixed_side_retained_unit_metadata(
        bad_ready_inventory,
    )
    bad_ready_metadata =
        metrics_module._pqs_current_route_fixed_side_retained_unit_metadata(
            bad_ready_inventory;
            strict = false,
        )
    @test bad_ready_metadata.diagnostics.shell_realized_pqs_source_box_operator_ready_count ==
          1
    @test !bad_ready_metadata.diagnostics.shell_realized_pqs_fixtures_are_metadata_oracle_only
    @test_throws ArgumentError metrics_module._pqs_pqs_product_component_route_smoke_cr2_sidecar_schema(
        report;
        fixed_side_retained_unit_metadata = bad_ready_metadata,
    )

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
    @test_throws DimensionMismatch metrics_module._pqs_pqs_product_component_route_smoke_cr2_sidecar_schema(
        metrics_module._pqs_pqs_product_component_route_smoke_report_adapter(
            summaries,
            _synthetic_mwg_facts(;
                overrides = Dict(
                    "residual_nucleus_indices" => "[1, 2]",
                ),
            ),
        ),
    )
end

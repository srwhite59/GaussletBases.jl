# Text report and optional artifact helpers for `bin/cartesian_ham_builder.jl`.
#
# Keep these after report construction so the core route helper file can focus
# on setup, skeletons, diagnostics, and dry-run report assembly.


# Text report helpers. These are intentionally simple print utilities, not a
# general logging framework.

function _pqs_route_driver_print_section(title)
    println()
    println("[", title, "]")
    return nothing
end

function _pqs_route_driver_print_kv(key, value)
    println(rpad(String(key), 42), "  ", value)
    return nothing
end

function _pqs_route_driver_print_named_tuple(title, values)
    _pqs_route_driver_print_section(title)
    for field in keys(values)
        _pqs_route_driver_print_kv(field, getproperty(values, field))
    end
    return nothing
end

function _pqs_source_box_route_driver_print_details(report)
    _pqs_route_driver_print_named_tuple("system_metadata", report.system_metadata)
    _pqs_route_driver_print_named_tuple("recipe_metadata", report.recipe_metadata)

    _pqs_route_driver_print_section("standard_setup")
    for field in (
        :status,
        :n_s,
        :n_s_source,
        :core_cube_side,
        :parent_box,
        :core_spacing,
        :mapping_s,
        :mapping_s_by_atom,
        :q_to_core_spacing_rule,
    )
        _pqs_route_driver_print_kv(field, getproperty(report.standard_setup, field))
    end
    for field in (
        :q_to_core_spacing_rule_status,
        :q_to_core_spacing_provenance,
        :core_spacing_source,
        :core_spacing_default_formula,
        :q_to_core_spacing_non_optimality_claim,
    )
        _pqs_route_driver_print_kv(
            field,
            getproperty(report.standard_setup.diagnostics, field),
        )
    end

    _pqs_route_driver_print_section("parent_axis_readiness")
    for field in (
        :status,
        :core_spacing_available,
        :white_lindsey_spacing_facts_available,
        :charge_family,
        :geometry,
        :extent_candidates,
        :parent_axis_counts,
        :parent_axis_counts_status,
        :parent_axis_counts_manual_fixture,
        :parent_axis_counts_derived,
        :existing_parent_api_appears_applicable,
        :standard_parent_axis_rule_ready,
        :parent_axis_metadata_constructed,
        :construction_decision,
        :pending_facts,
    )
        _pqs_route_driver_print_kv(field, getproperty(report.parent_axis_readiness, field))
    end

    if !isnothing(report.parent_axis_probe)
        _pqs_route_driver_print_section("parent_axis_probe")
        for field in (
            :status,
            :basis_metadata,
            :axis_bundle_metadata,
            :axis_lengths,
            :physical_extent_inputs,
            :core_spacing,
            :reference_spacing,
            :tail_spacing,
            :gausslet_backend,
            :gausslet_backend_role,
            :expansion_source,
            :explicit_spacing_probe_only,
            :default_standard_rule,
            :parent_axis_metadata_constructed,
            :pending_facts,
        )
            _pqs_route_driver_print_kv(field, getproperty(report.parent_axis_probe, field))
        end
    end

    _pqs_route_driver_print_section("route_axis_counts")
    for field in (
        :status,
        :parent_axis_counts,
        :parent_axis_counts_source,
        :parent_axis_counts_derived,
        :parent_axis_counts_manual_fixture,
        :parent_axis_probe_status,
        :parent_axis_readiness_status,
        :q_minimum_satisfied,
        :pending_facts,
    )
        _pqs_route_driver_print_kv(field, getproperty(report.route_axis_counts, field))
    end

    if !isnothing(report.raw_product_box_probe)
        _pqs_route_driver_print_section("raw_product_box_probe")
        for field in (
            :status,
            :raw_product_box_plan_count,
            :all_pgdg_exact,
            :any_numerical_reference_fallback,
            :max_axis_overlap_error,
            :gausslet_backend,
            :gausslet_backend_role,
            :pending_facts,
        )
            _pqs_route_driver_print_kv(field, getproperty(report.raw_product_box_probe, field))
        end
        for metadata in report.raw_product_box_probe.unit_plan_metadata
            println(
                metadata.unit_key, '\t', metadata.source_box, '\t',
                metadata.source_mode_dims, '\t', metadata.source_mode_count, '\t',
                metadata.integration_contract, '\t',
                metadata.numerical_reference_fallback, '\t',
                metadata.max_axis_overlap_error,
            )
        end
    end

    _pqs_route_driver_print_named_tuple("parent_contract", report.parent_contract)
    _pqs_route_driver_print_named_tuple("parent_description", report.parent_description)
    _pqs_route_driver_print_named_tuple("source_boxes", report.source_boxes)
    _pqs_route_driver_print_named_tuple(
        "standard_unit_inventory", report.standard_unit_inventory)

    _pqs_route_driver_print_section("retained_units")
    for unit in report.retained_units
        println(
            unit.unit_key, '\t', unit.unit_role, '\t',
            unit.retained_count, '\t', unit.retained_range, '\t',
            unit.source_box, '\t',
            unit.retained_rule_derivation,
        )
    end

    _pqs_route_driver_print_section("pair_inventory")
    pair_entries = report.pair_entries
    pair_family_counts = report.pair_family_counts
    @show length(pair_entries)
    @show pair_family_counts
    for entry in pair_entries
        println(
            entry.pair_key, '\t', entry.pair_family, '\t',
            entry.density_density_helper, '\t',
            entry.transpose_policy,
        )
    end

    _pqs_route_driver_print_named_tuple("linear_algebra_plan", report.linear_algebra_plan)
    _pqs_route_driver_print_named_tuple("diagnostics", report.diagnostics)
    return nothing
end

function _pqs_source_box_route_driver_materialization_status_fields()
    return (
        :route_family,
        :materialize_route_requested,
        :save_basis_artifact_requested,
        :save_ham_artifact_requested,
        :route_configured_diatomic_ham_interaction_treatment_requested,
        :route_configured_diatomic_ham_interaction_treatment_consumed,
        :route_configured_diatomic_ham_interaction_treatment_status,
        :status,
        :materialized_report_kind,
        :route_configured_shellization_request_available,
        :route_configured_shellization_request_status,
        :route_configured_system_classification,
        :route_configured_system_classification_status,
        :route_configured_bond_axis,
        :route_configured_shellization_plan_available,
        :route_configured_shellization_plan_status,
        :route_configured_shellization_planning_status,
        :route_configured_shellization_planning_family,
        :route_configured_midpoint_slab_status,
        :route_configured_shellization_helper_map_available,
        :route_configured_shellization_helper_map_status,
        :route_configured_primary_planned_helper,
        :route_configured_missing_input_count,
        :route_configured_helper_map_blocker,
        :route_configured_input_readiness_available,
        :route_configured_input_readiness_status,
        :route_configured_available_fact_count,
        :route_configured_materializer_missing_input_count,
        :route_configured_input_readiness_blocker,
        :route_configured_materializer_config_available,
        :route_configured_materializer_config_status,
        :route_configured_materializer_config_planning_family,
        :route_configured_materializer_config_pending_input_count,
        :route_configured_materializer_backend_requested,
        :route_configured_materializer_backend_source,
        :route_configured_materializer_backend_status,
        :route_configured_materializer_backend_consumed,
        :route_configured_materializer_d_requested,
        :route_configured_materializer_d_source,
        :route_configured_materializer_d_status,
        :route_configured_materializer_d_consumed,
        :route_configured_materializer_nside_requested,
        :route_configured_materializer_nside_source,
        :route_configured_materializer_nside_status,
        :route_configured_materializer_nside_consumed,
        :route_configured_materializer_reference_spacing_requested,
        :route_configured_materializer_reference_spacing_source,
        :route_configured_materializer_reference_spacing_status,
        :route_configured_materializer_reference_spacing_consumed,
        :route_configured_materializer_tail_spacing_requested,
        :route_configured_materializer_tail_spacing_source,
        :route_configured_materializer_tail_spacing_status,
        :route_configured_materializer_tail_spacing_consumed,
        :route_configured_materializer_options_ready,
        :route_configured_materializer_missing_options,
        :route_configured_materializer_option_blocker,
        :route_configured_one_center_materializer_probe_requested,
        :route_configured_one_center_materializer_probe_status,
        :route_configured_one_center_materializer_probe_materialized,
        :route_configured_one_center_materializer_probe_consumed,
        :route_configured_one_center_materializer_probe_blocker,
        :route_configured_diatomic_materializer_probe_requested,
        :route_configured_diatomic_materializer_probe_status,
        :route_configured_diatomic_materializer_probe_materialized,
        :route_configured_diatomic_materializer_probe_consumed,
        :route_configured_diatomic_materializer_probe_blocker,
        :route_configured_diatomic_materializer_missing_contract,
        :route_configured_diatomic_materializer_payload_available,
        :route_configured_diatomic_parent_qw_basis_object_handoff_available,
        :route_configured_diatomic_parent_axis_bundle_object_handoff_available,
        :route_configured_diatomic_axis_bundle_backend_handoff_available,
        :route_configured_diatomic_axis_bundle_backend_handoff,
        :route_configured_diatomic_shared_shell_layer_policy,
        :route_configured_diatomic_packet_kernel,
        :route_configured_diatomic_policy_source,
        :route_configured_diatomic_seed_fallback,
        :route_configured_diatomic_atom_growth_materializer_probe_requested,
        :route_configured_diatomic_atom_growth_materializer_probe_status,
        :route_configured_diatomic_atom_growth_materializer_probe_materialized,
        :route_configured_diatomic_atom_growth_materializer_probe_consumed,
        :route_configured_diatomic_atom_growth_materializer_probe_blocker,
        :route_configured_diatomic_atom_growth_materializer_missing_contract,
        :route_configured_diatomic_atom_growth_sequence_available,
        :route_configured_diatomic_atom_growth_retained_dimension,
        :route_configured_diatomic_atom_growth_support_count,
        :route_configured_diatomic_atom_growth_coverage_complete,
        :route_configured_diatomic_atom_growth_active_source_authority,
        :route_configured_diatomic_atom_growth_plan_authority,
        :route_configured_diatomic_atom_growth_calls_shellification_plan_materializer,
        :route_configured_diatomic_atom_growth_fixed_block_available,
        :route_configured_diatomic_atom_growth_fixed_block_status,
        :route_configured_diatomic_atom_growth_basis_adapter_status,
        :route_configured_diatomic_atom_growth_basis_adapter_blocker,
        :route_configured_diatomic_atom_growth_final_integral_weights_status,
        :route_configured_diatomic_atom_growth_ham_adapter_status,
        :route_configured_diatomic_atom_growth_ham_adapter_blocker,
        :shellization_summary_available,
        :shellization_source,
        :route_configured_shellization_consumed,
        :materialized_shellization_stage,
        :seed_materialization_status,
        :retained_dimension,
        :final_integral_weights_status,
        :one_body_operator_status,
        :basis_bundle_export_status,
        :basis_artifact_status,
        :basis_artifact_written,
        :basisfile,
        :basis_artifact_path,
        :basis_export_blocker,
        :ham_preflight_status,
        :ham_missing_builder,
        :ham_operator_payload_status,
        :ham_interaction_status,
        :ham_bundle_export_status,
        :ham_artifact_status,
        :ham_artifact_written,
        :hamfile,
        :ham_export_blocker,
        :pqs_materialization_status,
    )
end

function _pqs_source_box_route_driver_print_materialization(materialization)
    _pqs_route_driver_print_section("route_materialization")
    for field in _pqs_source_box_route_driver_materialization_status_fields()
        _pqs_route_driver_print_kv(field, getproperty(materialization, field))
    end
    return nothing
end


# Optional artifact helpers for the private driver.

function _pqs_route_driver_write_tsv_row(io, section, key, value)
    println(io, section, '\t', key, '\t', repr(value))
    return nothing
end

function _pqs_route_driver_write_named_tuple_tsv(io, section, values)
    for field in keys(values)
        _pqs_route_driver_write_tsv_row(io, section, field, getproperty(values, field))
    end
    return nothing
end


function _pqs_source_box_route_driver_durable_materializer_payload(payload)
    isnothing(payload) && return nothing
    return (;
        object_kind = payload.object_kind,
        private_development_only = payload.private_development_only,
        transient_only = payload.transient_only,
        durable_report_serialization = :sanitized_before_save,
        source = payload.source,
        parent_qw_basis_object_available = payload.parent_qw_basis_object_available,
        parent_axis_bundle_object_available =
            payload.parent_axis_bundle_object_available,
        parent_qw_basis_object_type_label = payload.parent_qw_basis_object_type_label,
        parent_axis_bundle_object_type_label =
            payload.parent_axis_bundle_object_type_label,
        axis_bundle_backend = payload.axis_bundle_backend,
        axis_bundle_backend_available = payload.axis_bundle_backend_available,
        heavy_objects_elided = true,
    )
end

function _pqs_source_box_route_driver_durable_report(report)
    hasproperty(report, :route_materializer_payload) || return report
    return (;
        (
            field => field == :route_materializer_payload ?
                     _pqs_source_box_route_driver_durable_materializer_payload(
                         report.route_materializer_payload,
                     ) :
                     getproperty(report, field) for field in keys(report)
        )...,
    )
end

function _pqs_source_box_route_driver_durable_diatomic_materialization(materialization)
    isnothing(materialization) && return nothing
    hasproperty(materialization, :source) || return materialization
    sanitized = (;
        (
            field => field == :source ?
                     nothing :
                     getproperty(materialization, field) for
            field in keys(materialization)
        )...,
    )
    return merge(sanitized, (; source_elided_for_durable_report = true))
end

function _pqs_source_box_route_driver_durable_diatomic_probe(probe)
    isnothing(probe) && return nothing
    hasproperty(probe, :materialization) || return probe
    return (;
        (
            field => field == :materialization ?
                     _pqs_source_box_route_driver_durable_diatomic_materialization(
                         probe.materialization,
                     ) :
                     getproperty(probe, field) for field in keys(probe)
        )...,
    )
end

function _pqs_source_box_route_driver_durable_materialization(materialization)
    isnothing(materialization) && return nothing
    hasproperty(materialization, :route_configured_diatomic_materializer_probe) ||
        return materialization
    return (;
        (
            field => field == :route_configured_diatomic_materializer_probe ?
                     _pqs_source_box_route_driver_durable_diatomic_probe(
                         materialization.route_configured_diatomic_materializer_probe,
                     ) :
                     getproperty(materialization, field) for field in keys(materialization)
        )...,
    )
end

function _pqs_source_box_route_driver_save(
    report;
    save_artifact, save_tsv, outfile, tsvfile, materialization = nothing,
)
    durable_report = _pqs_source_box_route_driver_durable_report(report)
    durable_materialization =
        _pqs_source_box_route_driver_durable_materialization(materialization)

    if save_artifact
        println("saving JLD2 report ", outfile)
        if isnothing(durable_materialization)
            jldsave(outfile; report = durable_report)
        else
            jldsave(outfile; report = durable_report, materialization = durable_materialization)
        end
    end

    if save_tsv
        println("saving TSV report ", tsvfile)
        open(tsvfile, "w") do io
            println(io, "section\tkey\tvalue")
            _pqs_route_driver_write_named_tuple_tsv(io, "system_metadata", report.system_metadata)
            _pqs_route_driver_write_named_tuple_tsv(io, "recipe_metadata", report.recipe_metadata)
            _pqs_route_driver_write_named_tuple_tsv(
                io, "parent_contract", report.parent_contract)
            _pqs_route_driver_write_named_tuple_tsv(
                io, "standard_setup_diagnostics", report.standard_setup.diagnostics,
            )
            _pqs_route_driver_write_named_tuple_tsv(
                io, "parent_axis_readiness_diagnostics",
                report.parent_axis_readiness.diagnostics,
            )
            if !isnothing(report.parent_axis_probe)
                _pqs_route_driver_write_named_tuple_tsv(
                    io, "parent_axis_probe_diagnostics",
                    report.parent_axis_probe.diagnostics,
                )
            end
            _pqs_route_driver_write_named_tuple_tsv(
                io, "route_axis_counts_diagnostics", report.route_axis_counts.diagnostics,
            )
            if !isnothing(report.raw_product_box_probe)
                _pqs_route_driver_write_named_tuple_tsv(
                    io, "raw_product_box_probe_diagnostics",
                    report.raw_product_box_probe.diagnostics,
                )
            end
            _pqs_route_driver_write_named_tuple_tsv(
                io, "standard_unit_inventory", report.standard_unit_inventory)
            for unit in report.retained_units
                _pqs_route_driver_write_tsv_row(io, "retained_unit", unit.unit_key, unit)
            end
            for entry in report.pair_entries
                _pqs_route_driver_write_tsv_row(io, "pair_entry", entry.pair_key, entry)
            end
            if !isnothing(materialization)
                for field in _pqs_source_box_route_driver_materialization_status_fields()
                    _pqs_route_driver_write_tsv_row(
                        io,
                        "route_materialization",
                        field,
                        getproperty(materialization, field),
                    )
                end
            end
            _pqs_route_driver_write_named_tuple_tsv(io, "diagnostics", report.diagnostics)
        end
    end
    return nothing
end

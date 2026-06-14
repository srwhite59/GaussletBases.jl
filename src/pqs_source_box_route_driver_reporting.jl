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
        :low_order_shellization_policy_requested,
        :low_order_shellization_policy_resolved,
        :low_order_shellization_policy_source,
        :low_order_shellization_policy_status,
        :low_order_shellization_policy_blocker,
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
        :route_configured_diatomic_atom_growth_shellification_consumed,
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
        :route_configured_legacy_diatomic_source_consumed,
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

function _pqs_source_box_route_driver_write_group!(file, group, values)
    for key in keys(values)
        file[string(group, "/", key)] = getproperty(values, key)
    end
end

function _pqs_source_box_route_driver_write_present_group!(file, group, values)
    for key in keys(values)
        value = getproperty(values, key)
        isnothing(value) || (file[string(group, "/", key)] = value)
    end
end

function _pqs_source_box_route_driver_ordered_target_counts(counts, order)
    isnothing(counts) && return nothing
    isempty(order) && return ()
    return Tuple(getproperty(counts, key) for key in order)
end

function _pqs_source_box_route_driver_write_pqs_he_artifact!(file, report, input_path)
    hasproperty(report, :complete_core_shell_h1_j_driver_route_materialized) &&
        report.complete_core_shell_h1_j_driver_route_materialized || return nothing
    recipe = report.recipe_metadata
    h1 = report.complete_core_shell_h1_j_h1_energy
    j = report.complete_core_shell_h1_j_self_coulomb
    wl_h1 = recipe.wl_h1_lowest
    wl_j = recipe.wl_h1_self_coulomb
    for (group, values) in (
        ("config", (; input_path = isnothing(input_path) ? "" : String(input_path), route_family = recipe.route_family, route_kind = recipe.route_kind, parent_mapping_rule = recipe.parent_mapping_rule, parent_mapping_Z = recipe.parent_mapping_Z, parent_mapping_d = recipe.parent_mapping_d, tail_spacing = recipe.tail_spacing, q = recipe.q, n_s = recipe.n_s)),
        ("basis", (; final_dimension = report.complete_core_shell_h1_j_final_dimension, core_support_count = report.complete_core_shell_core_support_count, shell_support_count = report.complete_core_shell_shell_support_count, shell_layer_count = report.complete_core_shell_shell_layer_count, retained_per_shell = report.complete_core_shell_retained_per_shell, shell_final_retained_count = report.complete_core_shell_shell_final_retained_count, final_overlap_identity_error = report.complete_core_shell_final_overlap_identity_error)),
        ("physics", (; h1_lowest = h1, h1_j_self_coulomb = j)),
        ("density_interaction", (; status = report.complete_core_shell_h1_j_diagnostic_status, density_gauge = report.complete_core_shell_h1_j_density_gauge, raw_pair_factor_convention = report.complete_core_shell_raw_pair_factor_convention)),
        ("comparison", (; reference_label = something(recipe.comparison_reference_label, ""), wl_h1_lowest = wl_h1, wl_h1_self_coulomb = wl_j, delta_h1 = h1 - wl_h1, delta_h1_j = j - wl_j)),
    )
        _pqs_source_box_route_driver_write_group!(file, group, values)
    end
    if get(report, :private_rhf_requested, false)
        _pqs_source_box_route_driver_write_present_group!(
            file,
            "private_rhf",
            (;
                status = report.private_rhf_status,
                blocker = report.private_rhf_blocker,
                total_energy = report.private_rhf_total_energy,
                iteration_count = report.private_rhf_iteration_count,
                converged = report.private_rhf_converged,
                residual = report.private_rhf_residual,
                mixing_kind = report.private_rhf_mixing_kind,
            ),
        )
        _pqs_source_box_route_driver_write_present_group!(
            file,
            "comparison",
            (;
                wl_rhf_total = report.private_rhf_wl_total,
                delta_rhf = report.private_rhf_delta,
            ),
        )
    end
    return nothing
end

function _pqs_source_box_route_driver_bond_length(atom_locations)
    length(atom_locations) == 2 || return nothing
    a = atom_locations[1]
    b = atom_locations[2]
    length(a) == 3 && length(b) == 3 || return nothing
    return sqrt(sum((Float64(a[axis]) - Float64(b[axis]))^2 for axis in 1:3))
end

function _pqs_source_box_route_driver_write_pqs_diatomic_readiness_artifact!(
    file,
    report,
    input_path,
)
    recipe = report.recipe_metadata
    system = report.system_metadata
    recipe.route_family === :pqs_source_box || return nothing
    length(system.atom_symbols) == 2 || return nothing
    get(recipe, :supplement_policy, nothing) === :none || return nothing

    atom_locations = Tuple(Tuple(location) for location in system.atom_locations)
    bond_axis =
        !isnothing(get(system, :bond_axis, nothing)) ?
        system.bond_axis :
        get(report.parent_contract, :bond_axis, nothing)
    bond_length =
        !isnothing(get(system, :bond_length, nothing)) ?
        system.bond_length :
        _pqs_source_box_route_driver_bond_length(atom_locations)
    comparison_ready = get(recipe, :comparison_ready, false)
    comparison_blocker =
        get(recipe, :comparison_blocker, nothing)
    readiness = get(
        report,
        :diatomic_complete_core_shell_readiness_summary,
        (;
            status = :not_available_missing_diatomic_complete_core_shell_readiness,
            blocker = :missing_diatomic_complete_core_shell_readiness,
        ),
    )
    target = get(
        report,
        :physical_gausslet_target_summary,
        (;
            status = :not_available_missing_physical_gausslet_target_payload,
            blocker = :missing_physical_gausslet_target_payload,
            support_units = (),
            retained_units = (),
            retained_order = (),
            support_counts = nothing,
            retained_counts = nothing,
            expected_final_dimension = nothing,
            retained_atom_core_interiors = nothing,
            source_plan_role = nothing,
            source_plan_status = :not_available,
            source_plan_blocker = nothing,
            supplement_policy = nothing,
        ),
    )
    physical_target_artifact =
        get(recipe, :artifact_role, nothing) === :physical_gausslet_endpoint_target
    route_source_plan_status =
        physical_target_artifact ?
        get(target, :source_plan_status, get(readiness, :source_plan_status, :not_available)) :
        get(readiness, :source_plan_status, :not_available)
    physics_endpoint_blocker =
        physical_target_artifact ?
        get(target, :source_plan_blocker, get(recipe, :physics_endpoint_blocker, nothing)) :
        get(recipe, :physics_endpoint_blocker, nothing)

    for (group, values) in (
        ("system", (; atom_symbols = system.atom_symbols, nuclear_charges = system.nuclear_charges, atom_locations, bond_axis, bond_length)),
        ("config", (; input_path = isnothing(input_path) ? "" : String(input_path), route_family = recipe.route_family, route_kind = recipe.route_kind, q = recipe.q, n_s = recipe.n_s, core_spacing = recipe.core_spacing, xmax_parallel = get(recipe, :xmax_parallel, nothing), xmax_transverse = get(recipe, :xmax_transverse, nothing), supplement_policy = recipe.supplement_policy, comparison_ready, run_final_basis = get(recipe, :run_final_basis, false))),
        ("comparison", (; ready = comparison_ready, blocker = comparison_blocker, reference_label = something(recipe.comparison_reference_label, ""))),
        ("parent", (; parent_axis_counts = report.parent_contract.parent_axis_counts, parent_axis_counts_source = report.parent_contract.parent_axis_counts_source, parent_materialization_blocker = report.parent_contract.parent_materialization_blocker, parent_basis_object_available = report.parent_contract.parent_basis_object_available, parent_qw_basis_object_available = report.parent_contract.parent_qw_basis_object_available, parent_axis_bundle_object_available = report.parent_contract.parent_axis_bundle_object_available, parent_basis_object_type_label = report.parent_contract.parent_basis_object_type_label, parent_qw_basis_object_type_label = report.parent_contract.parent_qw_basis_object_type_label, parent_axis_bundle_object_type_label = report.parent_contract.parent_axis_bundle_object_type_label)),
        ("route", (; artifact_role = get(recipe, :artifact_role, nothing), readiness_status = get(readiness, :status, :not_available), readiness_blocker = get(readiness, :blocker, nothing), source_plan_status = route_source_plan_status, final_basis_status = get(readiness, :final_basis_status, :not_available), h1_status = get(readiness, :h1_status, :not_available), h1_materialized = get(readiness, :h1_materialized, false), h1_j_materialized = get(readiness, :h1_j_materialized, false), ham_input_status = get(readiness, :ham_input_payload_status, :not_available), hamiltonian_handoff_status = get(readiness, :hamiltonian_handoff_payload_status, :not_available), private_rhf_materialized = get(readiness, :rhf_materialized, false), public_api = get(readiness, :public_api, false), exports_materialized = get(readiness, :exports_materialized, false), artifacts_materialized = get(readiness, :artifacts_materialized, false))),
    )
        _pqs_source_box_route_driver_write_group!(file, group, values)
    end
    _pqs_source_box_route_driver_write_present_group!(
        file,
        "basis",
        (;
            final_dimension = get(readiness, :final_dimension, nothing),
            final_overlap_identity_error =
                get(readiness, :final_overlap_identity_error, nothing),
            retained_atom_core_interiors =
                get(recipe, :retained_atom_core_interiors, nothing),
            source_plan_role = get(recipe, :source_plan_role, nothing),
        ),
    )
    _pqs_source_box_route_driver_write_group!(
        file,
        "target",
        (;
            status = get(target, :status, nothing),
            blocker = get(target, :blocker, nothing),
            support_units = get(target, :support_units, ()),
            support_counts =
                _pqs_source_box_route_driver_ordered_target_counts(
                    get(target, :support_counts, nothing),
                    get(target, :support_units, ()),
                ),
            retained_units = get(target, :retained_units, ()),
            retained_counts =
                _pqs_source_box_route_driver_ordered_target_counts(
                    get(target, :retained_counts, nothing),
                    get(target, :retained_order, ()),
                ),
            retained_order = get(target, :retained_order, ()),
            expected_final_dimension =
                get(target, :expected_final_dimension, nothing),
            retained_atom_core_interiors =
                get(target, :retained_atom_core_interiors, nothing),
            source_plan_role = get(target, :source_plan_role, nothing),
            source_plan_status = get(target, :source_plan_status, nothing),
            source_plan_blocker = get(target, :source_plan_blocker, nothing),
            supplement_policy = get(target, :supplement_policy, nothing),
        ),
    )
    _pqs_source_box_route_driver_write_present_group!(
        file,
        "physics",
        (;
            endpoint_ready = get(recipe, :physics_endpoint_ready, nothing),
            endpoint_blocker = physics_endpoint_blocker,
            h1_lowest = get(readiness, :h1_lowest_energy, nothing),
            h1_hamiltonian_matrix_finite =
                get(readiness, :h1_hamiltonian_matrix_finite, nothing),
            h1_hamiltonian_symmetry_error =
                get(readiness, :h1_hamiltonian_symmetry_error, nothing),
        ),
    )
    _pqs_source_box_route_driver_write_present_group!(
        file,
        "density_interaction",
        (;
            density_gauge = get(readiness, :density_gauge, nothing),
            raw_pair_factor_convention =
                get(readiness, :raw_pair_factor_convention, nothing),
        ),
    )
    _pqs_source_box_route_driver_write_group!(
        file,
        "private_rhf",
        (;
            requested = get(recipe, :run_private_rhf, false),
            materialized = get(readiness, :rhf_materialized, false),
        ),
    )
    return nothing
end

function _pqs_source_box_route_driver_save(
    report;
    save_artifact, save_tsv, outfile, tsvfile, materialization = nothing,
    input_path = nothing,
)
    durable_report = _pqs_source_box_route_driver_durable_report(report)
    durable_materialization =
        _pqs_source_box_route_driver_durable_materialization(materialization)

    if save_artifact
        println("saving JLD2 report ", outfile)
        jldopen(outfile, "w") do file
            file["report"] = durable_report
            isnothing(durable_materialization) ||
                (file["materialization"] = durable_materialization)
            _pqs_source_box_route_driver_write_pqs_he_artifact!(
                file,
                durable_report,
                input_path,
            )
            _pqs_source_box_route_driver_write_pqs_diatomic_readiness_artifact!(
                file,
                durable_report,
                input_path,
            )
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

function _pqs_route_driver_probe_requested(value, core_spacing)
    value == :auto && return !isnothing(core_spacing)
    value isa Bool && return value
    throw(ArgumentError("probe_parent_axis_construction must be :auto, true, or false"))
end

function _pqs_route_driver_raw_box_probe_requested(value, parent_axis_probe, route_axis_counts)
    if value == :auto
        return !isnothing(parent_axis_probe) &&
               parent_axis_probe.parent_axis_metadata_constructed &&
               route_axis_counts.parent_axis_counts_source == :constructed_parent_axis_probe
    elseif value isa Bool
        return value
    end
    throw(ArgumentError("probe_raw_product_box_plans must be :auto, true, or false"))
end

function _pqs_source_box_route_driver_dry_run(;
    route_kind,
    atom_symbols,
    nuclear_charges,
    atom_locations,
    radius,
    parent_axis_counts,
    map_backend,
    q,
    n_s,
    reference_spacing,
    tail_spacing,
    q_to_core_spacing_rule,
    core_spacing,
    probe_parent_axis_construction,
    parent_axis_probe_backend,
    parent_axis_probe_family,
    probe_raw_product_box_plans,
    raw_product_box_probe_backend,
    route_shape,
    product_body_rule,
    pqs_retained_rule,
    product_retained_rule,
    terms,
    pair_factor_normalization,
    support_dense_direct_allowed,
    reference_only_authorities,
)
    pair_factor_normalization in (:density_normalized, :raw_weighted) || throw(
        ArgumentError(
            "pair_factor_normalization must be :density_normalized or :raw_weighted",
        ),
    )

    metrics = CartesianContractedParentMetrics
    standard_setup =
        metrics._pqs_standard_source_box_route_setup(
            ;
            nuclear_charges,
            atom_locations,
            q,
            radius,
            reference_spacing,
            tail_spacing,
            q_to_core_spacing_rule,
            core_spacing,
            n_s,
        )

    parent_axis_readiness =
        metrics._pqs_standard_parent_axis_construction_readiness(
            standard_setup;
            parent_axis_counts,
        )

    parent_axis_probe_requested =
        _pqs_route_driver_probe_requested(
            probe_parent_axis_construction,
            standard_setup.core_spacing,
        )
    parent_axis_probe = parent_axis_probe_requested ?
        metrics._pqs_explicit_core_spacing_parent_axis_probe(
            standard_setup;
            gausslet_backend = parent_axis_probe_backend,
            family = parent_axis_probe_family,
            construct_axis_bundles = true,
        ) : nothing
    parent_axis_probe_status =
        isnothing(parent_axis_probe) ? :not_requested : parent_axis_probe.status
    parent_axis_probe_constructed =
        isnothing(parent_axis_probe) ? false : parent_axis_probe.parent_axis_metadata_constructed
    parent_axis_probe_pending_facts =
        isnothing(parent_axis_probe) ? () : parent_axis_probe.pending_facts
    route_axis_counts =
        metrics._pqs_source_box_route_parent_axis_counts_for_skeleton(
            standard_setup,
            parent_axis_readiness,
            parent_axis_probe;
            manual_parent_axis_counts = parent_axis_counts,
        )

    system_metadata = (;
        atom_symbols,
        nuclear_charges = standard_setup.nuclear_charges,
        atom_locations = standard_setup.atom_locations,
        radius = standard_setup.radius,
        manual_parent_axis_counts = parent_axis_counts,
        parent_axis_counts = route_axis_counts.parent_axis_counts,
        parent_axis_counts_status = route_axis_counts.status,
        parent_axis_counts_source = route_axis_counts.parent_axis_counts_source,
        parent_box = standard_setup.parent_box,
        parent_box_rule = standard_setup.parent_box_rule,
        map_backend,
    )

    route_skeleton =
        metrics._pqs_pqs_product_source_box_route_skeleton(
            ;
            q,
            parent_axis_counts = route_axis_counts.parent_axis_counts,
            route_shape,
            product_body_rule,
            pqs_retained_rule,
            product_retained_rule,
            pair_factor_normalization,
        )
    raw_product_box_probe_requested =
        _pqs_route_driver_raw_box_probe_requested(
            probe_raw_product_box_plans,
            parent_axis_probe,
            route_axis_counts,
        )
    raw_product_box_probe = raw_product_box_probe_requested ?
        metrics._pqs_explicit_core_spacing_route_raw_product_box_plan_probe(
            standard_setup,
            route_skeleton;
            gausslet_backend = raw_product_box_probe_backend,
        ) : nothing
    raw_product_box_probe_status =
        isnothing(raw_product_box_probe) ? :not_requested : raw_product_box_probe.status
    raw_product_box_probe_pending_facts =
        isnothing(raw_product_box_probe) ? () : raw_product_box_probe.pending_facts

    recipe_metadata = (;
        route_kind,
        q,
        n_s = standard_setup.n_s,
        n_s_source = standard_setup.n_s_source,
        core_cube_side = standard_setup.core_cube_side,
        reference_spacing = standard_setup.reference_spacing,
        tail_spacing = standard_setup.tail_spacing,
        q_to_core_spacing_rule = standard_setup.q_to_core_spacing_rule,
        core_spacing = standard_setup.core_spacing,
        core_spacing_source = standard_setup.spacing.core_spacing_source,
        q_to_core_spacing_rule_status =
            standard_setup.spacing.q_to_core_spacing_rule_status,
        probe_parent_axis_construction,
        parent_axis_probe_requested,
        parent_axis_probe_backend,
        route_axis_counts_source = route_axis_counts.parent_axis_counts_source,
        probe_raw_product_box_plans,
        raw_product_box_probe_requested,
        raw_product_box_probe_backend,
        route_shape,
        product_body_rule,
        pqs_source_box_rule = :mode_selected_raw_box_pqs,
        pqs_retained_rule,
        product_retained_rule,
        terms,
        pair_factor_normalization,
        support_dense_direct_allowed,
        reference_only_authorities,
    )

    parent_description = (;
        status = :described_not_constructed,
        standard_setup,
        parent_axis_readiness,
        parent_axis_probe,
        route_axis_counts,
        raw_product_box_probe,
        physical_parent_box = standard_setup.parent_box,
        physical_parent_box_rule = standard_setup.parent_box_rule,
        axis_transform_status = parent_axis_readiness.status,
        one_dimensional_transforms = (:x_axis_transform, :y_axis_transform, :z_axis_transform),
        parent_lattice = :raw_product_box_parent_lattice,
        parent_axis_counts = route_skeleton.parent_axis_counts,
        parent_axis_counts_source = route_axis_counts.parent_axis_counts_source,
        source_boxes = route_skeleton.source_boxes,
        raw_product_box_plan_status = raw_product_box_probe_status,
        pending_facts = (
            route_skeleton.pending_facts...,
            :parent_axis_counts_from_standard_parent_constructor,
            parent_axis_readiness.pending_facts...,
            route_axis_counts.pending_facts...,
            raw_product_box_probe_pending_facts...,
        ),
    )

    source_boxes = route_skeleton.source_boxes
    source_dimensions = route_skeleton.source_dimensions
    retained_units = route_skeleton.retained_units
    retained_counts = route_skeleton.retained_counts
    ranges = route_skeleton.ranges
    retained_dimension = route_skeleton.retained_dimension
    pair_entries = route_skeleton.pair_entries
    pair_family_counts = route_skeleton.pair_family_counts
    helper_by_pair_family = route_skeleton.helper_by_pair_family

    linear_algebra_plan = (;
        retained_block_formula = "O_final[i,j] = T_i' * O_source_box_pair * T_j",
        assemble_complete_retained_matrix = true,
        dense_parent_matrix_required = false,
        product_pqs_policy = :transpose_of_pqs_product_only_after_symmetric_pair_factor_check,
        product_pqs_lower_blocks = ((:product, :pqs_left), (:product, :pqs_right)),
        finite_output_check = :required_when_operator_blocks_are_materialized,
        symmetry_error_check = :required_for_symmetric_same_route_input,
    )
    no_go_flags = (
        public_default_behavior = false,
        packet_fixed_block_qw_hamiltonian_adoption = false,
        mwg_ida_semantic_change = false,
        retained_weight_division = false,
        retained_pqs_weights_used = false,
        repo_side_ray_id = false,
        ecp_behavior = false,
        cr2_science_claim = false,
        shell_projection = false,
        lowdin_cleanup = false,
        support_local_shell_row_algorithm = false,
        support_coefficient_matrix = false,
    )

    stage_table = (
        (stage = 1, name = :collect_system_metadata, status = :represented),
        (stage = 2, name = :collect_recipe_metadata, status = :represented),
        (stage = 3, name = :construct_parent_object, status = :described_not_constructed),
        (stage = 4, name = :split_parent_into_product_type_units, status = :derived_by_helper),
        (stage = 5, name = :define_each_unit, status = :derived_by_helper),
        (stage = 6, name = :loop_over_unit_pairs, status = :derived_by_helper),
        (stage = 7, name = :apply_final_linear_algebra, status = :plan_reported),
        (stage = 8, name = :validate_report_save, status = :metadata_dry_run),
    )
    dry_run_validation = (;
        builds_real_hamiltonian = false,
        builds_route_matrices = false,
        finite_output = :not_run_metadata_only,
        symmetry_error = :not_run_metadata_only,
        reference_error = :unavailable_metadata_only,
        timing_allocation = :placeholder_only,
    )
    diagnostics = merge(
        route_skeleton.diagnostics,
        (
            source = :pqs_source_box_route_driver_skeleton,
            standard_setup_helper = :_pqs_standard_source_box_route_setup,
            standard_setup_status = standard_setup.status,
            standard_setup_diagnostics = standard_setup.diagnostics,
            parent_axis_readiness_helper =
                :_pqs_standard_parent_axis_construction_readiness,
            parent_axis_readiness_status = parent_axis_readiness.status,
            parent_axis_readiness_diagnostics = parent_axis_readiness.diagnostics,
            route_skeleton_helper = :_pqs_pqs_product_source_box_route_skeleton,
            n_s = standard_setup.n_s,
            n_s_source = standard_setup.n_s_source,
            core_cube_side = standard_setup.core_cube_side,
            core_cube_side_rule = standard_setup.core_cube_side_rule,
            parent_box_rule = standard_setup.parent_box_rule,
            parent_box = standard_setup.parent_box,
            core_spacing = standard_setup.core_spacing,
            mapping_s = standard_setup.mapping_s,
            q_to_core_spacing_rule = standard_setup.q_to_core_spacing_rule,
            q_to_core_spacing_rule_status =
                standard_setup.spacing.q_to_core_spacing_rule_status,
            q_to_core_spacing_provenance = standard_setup.spacing.provenance,
            core_spacing_source = standard_setup.spacing.core_spacing_source,
            core_spacing_default_formula =
                standard_setup.spacing.core_spacing_default_formula,
            q_to_core_spacing_non_optimality_claim =
                standard_setup.spacing.non_optimality_claim,
            parent_axis_counts_status =
                parent_axis_readiness.parent_axis_counts_status,
            parent_axis_counts_manual_fixture =
                parent_axis_readiness.parent_axis_counts_manual_fixture,
            parent_axis_counts_derived =
                parent_axis_readiness.parent_axis_counts_derived,
            existing_parent_api_appears_applicable =
                parent_axis_readiness.existing_parent_api_appears_applicable,
            standard_parent_axis_rule_ready =
                parent_axis_readiness.standard_parent_axis_rule_ready,
            route_axis_counts_helper =
                :_pqs_source_box_route_parent_axis_counts_for_skeleton,
            route_axis_counts_status = route_axis_counts.status,
            route_axis_counts_source = route_axis_counts.parent_axis_counts_source,
            route_axis_counts_derived = route_axis_counts.parent_axis_counts_derived,
            route_axis_counts_manual_fixture =
                route_axis_counts.parent_axis_counts_manual_fixture,
            route_axis_counts_diagnostics = route_axis_counts.diagnostics,
            parent_axis_probe_requested,
            parent_axis_probe_status,
            parent_axis_metadata_constructed = parent_axis_probe_constructed,
            parent_axis_probe_pending_facts,
            raw_product_box_probe_requested,
            raw_product_box_probe_status,
            raw_product_box_probe_pending_facts,
            raw_product_box_plan_count =
                isnothing(raw_product_box_probe) ?
                0 :
                raw_product_box_probe.raw_product_box_plan_count,
            raw_product_box_all_pgdg_exact =
                isnothing(raw_product_box_probe) ?
                false :
                raw_product_box_probe.all_pgdg_exact,
            raw_product_box_any_numerical_reference_fallback =
                isnothing(raw_product_box_probe) ?
                false :
                raw_product_box_probe.any_numerical_reference_fallback,
            parent_axis_pending_facts = parent_axis_readiness.pending_facts,
            output_representation = :retained_two_index_density_density,
            no_go_flags,
            driver_builds_real_hamiltonian = false,
            driver_builds_route_matrices = false,
        ),
    )

    return (;
        object_kind = :pqs_source_box_route_driver_skeleton_report,
        generated_at = string(time()),
        standard_setup,
        parent_axis_readiness,
        parent_axis_probe,
        route_axis_counts,
        raw_product_box_probe,
        system_metadata,
        recipe_metadata,
        parent_description,
        route_skeleton,
        route_shape = route_skeleton.route_shape,
        retained_unit_order = route_skeleton.retained_unit_order,
        source_boxes,
        source_dimensions,
        retained_units,
        retained_counts,
        ranges,
        retained_dimension,
        pair_entries,
        pair_family_counts,
        helper_by_pair_family,
        linear_algebra_plan,
        stage_table,
        dry_run_validation,
        diagnostics,
    )
end

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
                metadata.unit_key,
                '\t',
                metadata.source_box,
                '\t',
                metadata.source_mode_dims,
                '\t',
                metadata.source_mode_count,
                '\t',
                metadata.integration_contract,
                '\t',
                metadata.numerical_reference_fallback,
                '\t',
                metadata.max_axis_overlap_error,
            )
        end
    end

    _pqs_route_driver_print_named_tuple("parent_description", report.parent_description)
    _pqs_route_driver_print_named_tuple("source_boxes", report.source_boxes)

    _pqs_route_driver_print_section("retained_units")
    for unit in report.retained_units
        println(
            unit.unit_key,
            '\t',
            unit.unit_role,
            '\t',
            unit.retained_count,
            '\t',
            unit.retained_range,
            '\t',
            unit.source_box,
            '\t',
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
            entry.pair_key,
            '\t',
            entry.pair_family,
            '\t',
            entry.density_density_helper,
            '\t',
            entry.transpose_policy,
        )
    end

    _pqs_route_driver_print_named_tuple("linear_algebra_plan", report.linear_algebra_plan)
    _pqs_route_driver_print_named_tuple("diagnostics", report.diagnostics)
    return nothing
end

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

function _pqs_source_box_route_driver_save(
    report;
    save_artifact,
    save_tsv,
    outfile,
    tsvfile,
)
    if save_artifact
        println("saving JLD2 report ", outfile)
        jldsave(outfile; report)
    end

    if save_tsv
        println("saving TSV report ", tsvfile)
        open(tsvfile, "w") do io
            println(io, "section\tkey\tvalue")
            _pqs_route_driver_write_named_tuple_tsv(io, "system_metadata", report.system_metadata)
            _pqs_route_driver_write_named_tuple_tsv(io, "recipe_metadata", report.recipe_metadata)
            _pqs_route_driver_write_named_tuple_tsv(
                io,
                "standard_setup_diagnostics",
                report.standard_setup.diagnostics,
            )
            _pqs_route_driver_write_named_tuple_tsv(
                io,
                "parent_axis_readiness_diagnostics",
                report.parent_axis_readiness.diagnostics,
            )
            if !isnothing(report.parent_axis_probe)
                _pqs_route_driver_write_named_tuple_tsv(
                    io,
                    "parent_axis_probe_diagnostics",
                    report.parent_axis_probe.diagnostics,
                )
            end
            _pqs_route_driver_write_named_tuple_tsv(
                io,
                "route_axis_counts_diagnostics",
                report.route_axis_counts.diagnostics,
            )
            if !isnothing(report.raw_product_box_probe)
                _pqs_route_driver_write_named_tuple_tsv(
                    io,
                    "raw_product_box_probe_diagnostics",
                    report.raw_product_box_probe.diagnostics,
                )
            end
            for unit in report.retained_units
                _pqs_route_driver_write_tsv_row(io, "retained_unit", unit.unit_key, unit)
            end
            for entry in report.pair_entries
                _pqs_route_driver_write_tsv_row(io, "pair_entry", entry.pair_key, entry)
            end
            _pqs_route_driver_write_named_tuple_tsv(io, "diagnostics", report.diagnostics)
        end
    end
    return nothing
end

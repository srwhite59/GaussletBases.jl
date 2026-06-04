# Private support for `bin/pqs_source_box_route_driver.jl`.
#
# Keep route bookkeeping here so the executable driver can stay human-facing:
# editable defaults, overrides, visible stages, print, and save.


# Small local utility helpers.

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

function _pqs_route_driver_axis_counts(parent_axis_counts)
    isnothing(parent_axis_counts) && return nothing
    if hasproperty(parent_axis_counts, :x)
        return (
            x = Int(parent_axis_counts.x),
            y = Int(parent_axis_counts.y),
            z = Int(parent_axis_counts.z),
        )
    end
    length(parent_axis_counts) == 3 || throw(
        ArgumentError("route driver parent_axis_counts must have three axes"),
    )
    return (
        x = Int(parent_axis_counts[1]),
        y = Int(parent_axis_counts[2]),
        z = Int(parent_axis_counts[3]),
    )
end

function _pqs_route_driver_axis_count_tuple(counts)
    isnothing(counts) && return nothing
    return (counts.x, counts.y, counts.z)
end

function _pqs_route_driver_parent_source_box(counts)
    isnothing(counts) && return nothing
    return (x = 1:counts.x, y = 1:counts.y, z = 1:counts.z)
end


# Input checks and standard setup.

function _pqs_source_box_route_driver_check_inputs(route_recipe)
    route_recipe.route_family in (:pqs_source_box, :white_lindsey_low_order) || throw(
        ArgumentError("route_family must be :pqs_source_box or :white_lindsey_low_order"),
    )
    route_recipe.pair_factor_normalization in (:density_normalized, :raw_weighted) || throw(
        ArgumentError(
            "pair_factor_normalization must be :density_normalized or :raw_weighted",
        ),
    )
    return nothing
end

function _pqs_source_box_route_driver_standard_setup(system_inputs, spacing_inputs)
    metrics = CartesianContractedParentMetrics
    return metrics._pqs_standard_source_box_route_setup(
        ;
        nuclear_charges = system_inputs.nuclear_charges,
        atom_locations = system_inputs.atom_locations,
        q = spacing_inputs.q,
        radius = system_inputs.radius,
        reference_spacing = spacing_inputs.reference_spacing,
        tail_spacing = spacing_inputs.tail_spacing,
        q_to_core_spacing_rule = spacing_inputs.q_to_core_spacing_rule,
        core_spacing = spacing_inputs.core_spacing,
        n_s = spacing_inputs.n_s,
    )
end


# Parent-axis readiness/probe.

function _pqs_source_box_route_driver_parent_axis(
    standard_setup,
    system_inputs,
    probe_inputs,
)
    metrics = CartesianContractedParentMetrics
    parent_axis_readiness =
        metrics._pqs_standard_parent_axis_construction_readiness(
            standard_setup;
            parent_axis_counts = system_inputs.parent_axis_counts,
        )

    parent_axis_probe_requested =
        _pqs_route_driver_probe_requested(
            probe_inputs.probe_parent_axis_construction,
            standard_setup.core_spacing,
        )
    parent_axis_probe = parent_axis_probe_requested ?
        metrics._pqs_explicit_core_spacing_parent_axis_probe(
            standard_setup;
            gausslet_backend = probe_inputs.parent_axis_probe_backend,
            family = probe_inputs.parent_axis_probe_family,
            construct_axis_bundles = true,
        ) : nothing
    parent_axis_probe_status =
        isnothing(parent_axis_probe) ? :not_requested : parent_axis_probe.status
    parent_axis_probe_constructed =
        isnothing(parent_axis_probe) ? false : parent_axis_probe.parent_axis_metadata_constructed
    parent_axis_probe_pending_facts =
        isnothing(parent_axis_probe) ? () : parent_axis_probe.pending_facts

    return (;
        parent_axis_readiness,
        parent_axis_probe_requested,
        parent_axis_probe,
        parent_axis_probe_status,
        parent_axis_probe_constructed,
        parent_axis_probe_pending_facts,
    )
end


# Route axis counts.

function _pqs_source_box_route_driver_route_axis_counts(
    standard_setup,
    parent_axis,
    system_inputs,
    route_recipe,
)
    metrics = CartesianContractedParentMetrics
    if route_recipe.route_family == :pqs_source_box
        return metrics._pqs_source_box_route_parent_axis_counts_for_skeleton(
            standard_setup, parent_axis.parent_axis_readiness, parent_axis.parent_axis_probe;
            manual_parent_axis_counts = system_inputs.parent_axis_counts,
        )
    end

    parent_axis_probe = parent_axis.parent_axis_probe
    probe_constructed =
        !isnothing(parent_axis_probe) &&
        hasproperty(parent_axis_probe, :parent_axis_metadata_constructed) &&
        parent_axis_probe.parent_axis_metadata_constructed
    counts =
        probe_constructed ?
        _pqs_route_driver_axis_counts(parent_axis_probe.axis_lengths) :
        _pqs_route_driver_axis_counts(system_inputs.parent_axis_counts)
    counts_source =
        probe_constructed ? :constructed_parent_axis_probe :
        isnothing(system_inputs.parent_axis_counts) ? :unavailable : :manual_fixture
    pending_facts = isnothing(counts) ? (:manual_parent_axis_counts_or_constructed_parent_axis_probe,) : ()

    return (;
        object_kind = :cartesian_nesting_route_parent_axis_counts_for_skeleton,
        status = isnothing(counts) ? :not_available_pending_facts : :available,
        parent_axis_counts = counts,
        parent_axis_counts_source = counts_source,
        parent_axis_counts_derived = probe_constructed,
        parent_axis_counts_manual_fixture =
            counts_source == :manual_fixture,
        parent_axis_probe_status = parent_axis.parent_axis_probe_status,
        parent_axis_readiness_status = parent_axis.parent_axis_readiness.status,
        setup_object_kind = standard_setup.object_kind,
        q = standard_setup.q,
        q_minimum_satisfied =
            !isnothing(counts) &&
            counts.x >= standard_setup.q &&
            counts.y >= standard_setup.q &&
            counts.z >= standard_setup.q,
        pending_facts,
        diagnostics = (
            source = :cartesian_nesting_route_parent_axis_counts_for_skeleton,
            route_family = route_recipe.route_family,
            private_development_only = true,
            production_route = false,
            parent_axis_counts_source = counts_source,
            parent_axis_counts_derived = probe_constructed,
            parent_axis_counts_manual_fixture = counts_source == :manual_fixture,
            published_benchmark_route = true,
            source_box_route = false,
            public_default_consumes = false,
            packet_adoption = false,
            fixed_block_routing = false,
            qwhamiltonian_consumes = false,
            hamiltonian_matrix_built = false,
            shell_projection_used = false,
            lowdin_cleanup_used = false,
            support_local_shell_row_algorithm = false,
            support_coefficient_matrix_used = false,
            retained_pqs_weights_used = false,
            retained_weight_division_allowed = false,
            repo_side_ray_id = false,
            mwg_ida_semantics_changed = false,
            ecp_terms_implemented = false,
            cr2_science_status_changed = false,
        ),
    )
end


# Compatibility route-axis count helper for older source-box-only callers.

function _pqs_source_box_route_driver_route_axis_counts(
    standard_setup,
    parent_axis,
    system_inputs,
)
    metrics = CartesianContractedParentMetrics
    return metrics._pqs_source_box_route_parent_axis_counts_for_skeleton(
        standard_setup, parent_axis.parent_axis_readiness, parent_axis.parent_axis_probe;
        manual_parent_axis_counts = system_inputs.parent_axis_counts,
    )
end


# System metadata.

function _pqs_source_box_route_driver_system_metadata(
    standard_setup,
    route_axis_counts,
    system_inputs,
)
    return (;
        atom_symbols = system_inputs.atom_symbols,
        nuclear_charges = standard_setup.nuclear_charges,
        atom_locations = standard_setup.atom_locations,
        radius = standard_setup.radius,
        manual_parent_axis_counts = system_inputs.parent_axis_counts,
        parent_axis_counts = route_axis_counts.parent_axis_counts,
        parent_axis_counts_status = route_axis_counts.status,
        parent_axis_counts_source = route_axis_counts.parent_axis_counts_source,
        parent_box = standard_setup.parent_box,
        parent_box_rule = standard_setup.parent_box_rule,
        map_backend = system_inputs.map_backend,
    )
end

# Route-specific probes.

function _pqs_source_box_route_driver_raw_box_probe(
    standard_setup,
    route_skeleton,
    parent_axis,
    route_axis_counts,
    probe_inputs,
    route_recipe,
)
    metrics = CartesianContractedParentMetrics
    if route_recipe.route_family != :pqs_source_box
        return (;
            raw_product_box_probe_requested = false,
            raw_product_box_probe = nothing,
            raw_product_box_probe_status = :not_applicable_to_route_family,
            raw_product_box_probe_pending_facts = (),
        )
    end

    raw_product_box_probe_requested =
        _pqs_route_driver_raw_box_probe_requested(
            probe_inputs.probe_raw_product_box_plans,
            parent_axis.parent_axis_probe,
            route_axis_counts,
        )
    raw_product_box_probe = raw_product_box_probe_requested ?
        metrics._pqs_explicit_core_spacing_route_raw_product_box_plan_probe(
            standard_setup, route_skeleton;
            gausslet_backend = probe_inputs.raw_product_box_probe_backend,
        ) : nothing
    raw_product_box_probe_status =
        isnothing(raw_product_box_probe) ? :not_requested : raw_product_box_probe.status
    raw_product_box_probe_pending_facts =
        isnothing(raw_product_box_probe) ? () : raw_product_box_probe.pending_facts

    return (;
        raw_product_box_probe_requested,
        raw_product_box_probe,
        raw_product_box_probe_status,
        raw_product_box_probe_pending_facts,
    )
end


# Metadata and parent description.

function _pqs_source_box_route_driver_recipe_metadata(
    standard_setup,
    route_axis_counts,
    parent_axis,
    raw_box,
    spacing_inputs, probe_inputs, route_recipe,
)
    common_metadata = (;
        route_family = route_recipe.route_family,
        route_kind = route_recipe.route_kind,
        q = spacing_inputs.q,
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
        probe_parent_axis_construction = probe_inputs.probe_parent_axis_construction,
        parent_axis_probe_requested = parent_axis.parent_axis_probe_requested,
        parent_axis_probe_backend = probe_inputs.parent_axis_probe_backend,
        route_axis_counts_source = route_axis_counts.parent_axis_counts_source,
        probe_raw_product_box_plans = probe_inputs.probe_raw_product_box_plans,
        raw_product_box_probe_requested = raw_box.raw_product_box_probe_requested,
        raw_product_box_probe_backend = probe_inputs.raw_product_box_probe_backend,
        terms = route_recipe.terms,
        pair_factor_normalization = route_recipe.pair_factor_normalization,
    )

    if route_recipe.route_family == :pqs_source_box
        source_box_recipe = route_recipe.source_box
        return merge(
            common_metadata,
            (;
                route_shape = source_box_recipe.route_shape,
                product_body_rule = source_box_recipe.product_body_rule,
                pqs_source_box_rule = :mode_selected_raw_box_pqs,
                pqs_retained_rule = source_box_recipe.pqs_retained_rule,
                product_retained_rule = source_box_recipe.product_retained_rule,
                support_dense_direct_allowed =
                    source_box_recipe.support_dense_direct_allowed,
                reference_only_authorities =
                    source_box_recipe.reference_only_authorities,
            ),
        )
    end

    low_order_recipe = route_recipe.white_lindsey
    return merge(
        common_metadata,
        (;
            route_shape = low_order_recipe.route_shape,
            white_lindsey_mapping_rule = low_order_recipe.mapping_rule,
            white_lindsey_nesting_rule = low_order_recipe.nesting_rule,
            white_lindsey_retained_rule = low_order_recipe.retained_rule,
            white_lindsey_operator_rule = low_order_recipe.operator_rule,
            benchmark_role = low_order_recipe.benchmark_role,
            pqs_source_box_rule = :not_applicable,
            pqs_retained_rule = :not_applicable,
            product_retained_rule = :not_applicable,
            support_dense_direct_allowed = false,
            reference_only_authorities = (),
        ),
    )
end

function _pqs_source_box_route_driver_parent_description(
    standard_setup,
    parent_axis,
    route_axis_counts,
    route_skeleton,
    raw_box,
)
    source_box_route = route_skeleton.route_family == :pqs_source_box
    return (;
        status = :described_not_constructed,
        route_family = route_skeleton.route_family,
        standard_setup,
        parent_axis_readiness = parent_axis.parent_axis_readiness,
        parent_axis_probe = parent_axis.parent_axis_probe,
        route_axis_counts,
        raw_product_box_probe = raw_box.raw_product_box_probe,
        physical_parent_box = standard_setup.parent_box,
        physical_parent_box_rule = standard_setup.parent_box_rule,
        axis_transform_status = parent_axis.parent_axis_readiness.status,
        one_dimensional_transforms = (:x_axis_transform, :y_axis_transform, :z_axis_transform),
        parent_lattice =
            source_box_route ?
            :raw_product_box_parent_lattice :
            :white_lindsey_nested_cartesian_parent_lattice,
        parent_axis_counts = route_skeleton.parent_axis_counts,
        parent_axis_counts_source = route_axis_counts.parent_axis_counts_source,
        source_boxes = route_skeleton.source_boxes,
        raw_product_box_plan_status = raw_box.raw_product_box_probe_status,
        pending_facts = (
            route_skeleton.pending_facts...,
            :parent_axis_counts_from_standard_parent_constructor,
            parent_axis.parent_axis_readiness.pending_facts...,
            route_axis_counts.pending_facts...,
            raw_box.raw_product_box_probe_pending_facts...,
        ),
    )
end


# Route facts, contracts, diagnostics, and report assembly.

function _pqs_source_box_route_driver_shared_unit_fields()
    return (
        :unit_key,
        :unit_role,
        :retained_unit_kind,
        :source_family,
        :source_box,
        :source_dimensions,
        :source_dimension,
        :retained_rule_kind,
        :retained_rule_derivation,
        :retained_range,
        :retained_count,
        :provenance_label,
        :weight_semantics,
    )
end

function _pqs_source_box_route_driver_pair_entry_fields()
    return (
        :pair_key,
        :pair_family,
        :pair_kind,
        :density_density_helper,
        :source_box_algorithmic_path,
        :fallback_oracle_path,
        :transpose_policy,
        :output_representation,
    )
end

function _pqs_source_box_route_driver_check_record_fields(record, expected_fields, label)
    actual_fields = Tuple(keys(record))
    actual_fields == expected_fields || throw(
        ArgumentError("$(label) fields $(actual_fields) do not match $(expected_fields)"),
    )
    return record
end

function _pqs_source_box_route_driver_unit_record(;
    unit_key,
    unit_role,
    retained_unit_kind,
    source_family,
    source_box,
    source_dimensions,
    source_dimension = isnothing(source_dimensions) ? nothing : prod(source_dimensions),
    retained_rule_kind,
    retained_rule_derivation,
    retained_range,
    retained_count,
    provenance_label,
    weight_semantics,
)
    if !isnothing(source_dimensions) && !isnothing(source_dimension)
        derived_source_dimension = prod(source_dimensions)
        source_dimension == derived_source_dimension || throw(
            ArgumentError("source_dimension $(source_dimension) does not match $(derived_source_dimension)"),
        )
    end

    record = (;
        unit_key,
        unit_role,
        retained_unit_kind,
        source_family,
        source_box,
        source_dimensions,
        source_dimension,
        retained_rule_kind,
        retained_rule_derivation,
        retained_range,
        retained_count,
        provenance_label,
        weight_semantics,
    )
    return _pqs_source_box_route_driver_check_record_fields(
        record,
        _pqs_source_box_route_driver_shared_unit_fields(),
        :retained_unit,
    )
end

function _pqs_source_box_route_driver_named_tuple_from_units(retained_units, field)
    unit_keys = Tuple(unit.unit_key for unit in retained_units)
    values = Tuple(getproperty(unit, field) for unit in retained_units)
    return NamedTuple{unit_keys}(values)
end

function _pqs_source_box_route_driver_inventory_retained_dimension(retained_units)
    isempty(retained_units) && return 0
    any(unit -> isnothing(unit.retained_range), retained_units) && return nothing
    return maximum(last(unit.retained_range) for unit in retained_units)
end

function _pqs_source_box_route_driver_unit_inventory(retained_units)
    expected_fields = _pqs_source_box_route_driver_shared_unit_fields()
    for unit in retained_units
        _pqs_source_box_route_driver_check_record_fields(
            unit,
            expected_fields,
            :retained_unit,
        )
    end

    return (;
        retained_units,
        source_boxes =
            _pqs_source_box_route_driver_named_tuple_from_units(retained_units, :source_box),
        source_dimensions =
            _pqs_source_box_route_driver_named_tuple_from_units(
                retained_units, :source_dimensions),
        retained_counts =
            _pqs_source_box_route_driver_named_tuple_from_units(
                retained_units, :retained_count),
        ranges =
            _pqs_source_box_route_driver_named_tuple_from_units(
                retained_units, :retained_range),
        retained_dimension =
            _pqs_source_box_route_driver_inventory_retained_dimension(retained_units),
    )
end

function _pqs_source_box_route_driver_inventory_map_matches(existing, derived)
    length(keys(existing)) == length(keys(derived)) || return false
    for key in keys(derived)
        hasproperty(existing, key) || return false
        getproperty(existing, key) == getproperty(derived, key) || return false
    end
    return true
end

function _pqs_source_box_route_driver_pair_inventory(
    pair_entries;
    expected_families = nothing,
)
    expected_fields = _pqs_source_box_route_driver_pair_entry_fields()
    for entry in pair_entries
        _pqs_source_box_route_driver_check_record_fields(
            entry,
            expected_fields,
            :pair_entry,
        )
    end
    families = isnothing(expected_families) ?
        Tuple(unique(entry.pair_family for entry in pair_entries)) :
        Tuple(expected_families)
    pair_family_counts = NamedTuple{families}(
        Tuple(count(entry -> entry.pair_family == family, pair_entries) for family in families),
    )
    return (; pair_entries, pair_family_counts)
end

function _pqs_source_box_route_driver_standard_unit_inventory_summary(route_facts)
    retained_units = route_facts.retained_units
    pair_entries = route_facts.pair_entries
    retained_counts = route_facts.retained_counts
    retained_ranges = route_facts.ranges
    return (;
        unit_count = length(retained_units),
        unit_keys = Tuple(unit.unit_key for unit in retained_units),
        retained_unit_kinds = Tuple(unit.retained_unit_kind for unit in retained_units),
        source_families = Tuple(unit.source_family for unit in retained_units),
        source_dimensions = route_facts.source_dimensions,
        retained_counts = retained_counts,
        retained_dimension = route_facts.retained_dimension,
        retained_counts_materialized =
            all(!isnothing(getproperty(retained_counts, key)) for key in keys(retained_counts)),
        retained_ranges_materialized =
            all(!isnothing(getproperty(retained_ranges, key)) for key in keys(retained_ranges)),
        pair_count = length(pair_entries),
        pair_family_counts = route_facts.pair_family_counts,
        pair_families = Tuple(unique(entry.pair_family for entry in pair_entries)),
        output_representations =
            Tuple(unique(entry.output_representation for entry in pair_entries)),
    )
end

function _pqs_source_box_route_driver_route_facts(route_skeleton)
    unit_inventory =
        _pqs_source_box_route_driver_unit_inventory(route_skeleton.retained_units)
    _pqs_source_box_route_driver_inventory_map_matches(
        route_skeleton.source_boxes, unit_inventory.source_boxes) || throw(
        ArgumentError("route skeleton source_boxes do not match retained unit records"),
    )
    _pqs_source_box_route_driver_inventory_map_matches(
        route_skeleton.source_dimensions, unit_inventory.source_dimensions) || throw(
        ArgumentError("route skeleton source_dimensions do not match retained unit records"),
    )
    _pqs_source_box_route_driver_inventory_map_matches(
        route_skeleton.retained_counts, unit_inventory.retained_counts) || throw(
        ArgumentError("route skeleton retained_counts do not match retained unit records"),
    )
    _pqs_source_box_route_driver_inventory_map_matches(
        route_skeleton.ranges, unit_inventory.ranges) || throw(
        ArgumentError("route skeleton ranges do not match retained unit records"),
    )
    route_skeleton.retained_dimension == unit_inventory.retained_dimension || throw(
        ArgumentError("route skeleton retained_dimension does not match retained unit records"),
    )
    pair_inventory = _pqs_source_box_route_driver_pair_inventory(
        route_skeleton.pair_entries;
        expected_families = keys(route_skeleton.pair_family_counts),
    )
    route_skeleton.pair_family_counts == pair_inventory.pair_family_counts || throw(
        ArgumentError("route skeleton pair_family_counts do not match pair entries"),
    )

    return (;
        source_boxes = route_skeleton.source_boxes,
        source_dimensions = route_skeleton.source_dimensions,
        retained_units = route_skeleton.retained_units,
        retained_counts = route_skeleton.retained_counts,
        ranges = route_skeleton.ranges,
        retained_dimension = route_skeleton.retained_dimension,
        pair_entries = route_skeleton.pair_entries,
        pair_family_counts = route_skeleton.pair_family_counts,
        helper_by_pair_family = route_skeleton.helper_by_pair_family,
    )
end

function _pqs_source_box_route_driver_contract_metadata(route_recipe)
    if route_recipe.route_family == :white_lindsey_low_order
        linear_algebra_plan = (;
            retained_block_formula =
                "low-order nested Cartesian operator blocks in the White-Lindsey retained basis",
            assemble_complete_retained_matrix = :planned_benchmark_step,
            dense_parent_matrix_required = :implementation_dependent,
            product_pqs_policy = :not_applicable,
            product_pqs_lower_blocks = (),
            finite_output_check = :required_when_operator_blocks_are_materialized,
            symmetry_error_check = :required_for_symmetric_operator_blocks,
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
            pqs_source_box_algorithm_claim = false,
        )
        stage_table = (
            (stage = 1, name = :collect_system_metadata, status = :represented),
            (stage = 2, name = :select_route_family, status = :represented),
            (stage = 3, name = :construct_standard_unit_backbone_maps, status = :contract_only),
            (stage = 4, name = :apply_low_order_comx_to_unit_boxes, status = :contract_only),
            (stage = 5, name = :form_low_order_unit_basis, status = :pending_implementation),
            (stage = 6, name = :assemble_low_order_unit_operator_blocks, status = :pending_implementation),
            (stage = 7, name = :compare_against_pqs_source_box_route, status = :planned_benchmark),
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
        return (; linear_algebra_plan, no_go_flags, stage_table, dry_run_validation)
    end

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

    return (; linear_algebra_plan, no_go_flags, stage_table, dry_run_validation)
end

_pqs_source_box_route_driver_contract_metadata() =
    _pqs_source_box_route_driver_contract_metadata((; route_family = :pqs_source_box))

function _pqs_source_box_route_driver_diagnostics(
    standard_setup,
    parent_axis,
    route_axis_counts,
    route_skeleton,
    raw_box,
    contract,
)
    parent_axis_readiness = parent_axis.parent_axis_readiness
    raw_product_box_probe = raw_box.raw_product_box_probe
    source_box_route = route_skeleton.route_family == :pqs_source_box
    route_skeleton_helper =
        source_box_route ?
        :_pqs_pqs_product_source_box_route_skeleton :
        :_pqs_source_box_route_driver_white_lindsey_low_order_skeleton
    route_axis_counts_helper =
        source_box_route ?
        :_pqs_source_box_route_parent_axis_counts_for_skeleton :
        :_pqs_source_box_route_driver_route_axis_counts
    output_representation =
        hasproperty(route_skeleton.diagnostics, :output_representation) ?
        route_skeleton.diagnostics.output_representation :
        :retained_two_index_density_density
    diagnostics = merge(
        route_skeleton.diagnostics,
        (
            source = :cartesian_nesting_route_driver_skeleton,
            route_family = route_skeleton.route_family,
            standard_setup_helper = :_pqs_standard_source_box_route_setup,
            standard_setup_status = standard_setup.status,
            standard_setup_diagnostics = standard_setup.diagnostics,
            parent_axis_readiness_helper =
                :_pqs_standard_parent_axis_construction_readiness,
            parent_axis_readiness_status = parent_axis_readiness.status,
            parent_axis_readiness_diagnostics = parent_axis_readiness.diagnostics,
            route_skeleton_helper = route_skeleton_helper,
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
            route_axis_counts_helper = route_axis_counts_helper,
            route_axis_counts_status = route_axis_counts.status,
            route_axis_counts_source = route_axis_counts.parent_axis_counts_source,
            route_axis_counts_derived = route_axis_counts.parent_axis_counts_derived,
            route_axis_counts_manual_fixture =
                route_axis_counts.parent_axis_counts_manual_fixture,
            route_axis_counts_diagnostics = route_axis_counts.diagnostics,
            parent_axis_probe_requested = parent_axis.parent_axis_probe_requested,
            parent_axis_probe_status = parent_axis.parent_axis_probe_status,
            parent_axis_metadata_constructed = parent_axis.parent_axis_probe_constructed,
            parent_axis_probe_pending_facts = parent_axis.parent_axis_probe_pending_facts,
            raw_product_box_probe_requested = raw_box.raw_product_box_probe_requested,
            raw_product_box_probe_status = raw_box.raw_product_box_probe_status,
            raw_product_box_probe_pending_facts = raw_box.raw_product_box_probe_pending_facts,
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
            output_representation = output_representation,
            no_go_flags = contract.no_go_flags,
            driver_builds_real_hamiltonian = false,
            driver_builds_route_matrices = false,
        ),
    )
    return diagnostics
end

function _pqs_source_box_route_driver_report(
    standard_setup,
    parent_axis,
    route_axis_counts,
    raw_box,
    system_metadata,
    recipe_metadata,
    parent_description,
    route_skeleton,
    route_facts,
    contract,
    diagnostics,
)
    return (;
        object_kind = :cartesian_nesting_route_driver_skeleton_report,
        generated_at = Base.Libc.strftime("%Y-%m-%dT%H:%M:%S", time()),
        route_family = route_skeleton.route_family,
        standard_setup,
        parent_axis_readiness = parent_axis.parent_axis_readiness,
        parent_axis_probe = parent_axis.parent_axis_probe,
        route_axis_counts,
        raw_product_box_probe = raw_box.raw_product_box_probe,
        system_metadata,
        recipe_metadata,
        parent_description,
        route_skeleton,
        route_shape = route_skeleton.route_shape,
        retained_unit_order = route_skeleton.retained_unit_order,
        source_boxes = route_facts.source_boxes,
        source_dimensions = route_facts.source_dimensions,
        retained_units = route_facts.retained_units,
        retained_counts = route_facts.retained_counts,
        ranges = route_facts.ranges,
        retained_dimension = route_facts.retained_dimension,
        pair_entries = route_facts.pair_entries,
        pair_family_counts = route_facts.pair_family_counts,
        helper_by_pair_family = route_facts.helper_by_pair_family,
        standard_unit_inventory =
            _pqs_source_box_route_driver_standard_unit_inventory_summary(route_facts),
        linear_algebra_plan = contract.linear_algebra_plan,
        stage_table = contract.stage_table,
        dry_run_validation = contract.dry_run_validation,
        diagnostics,
    )
end

function _pqs_source_box_route_driver_white_lindsey_preflight_fixed_block(seed_or_fixed_block)
    if hasproperty(seed_or_fixed_block, :fixture) &&
       hasproperty(seed_or_fixed_block.fixture, :fixed_block)
        return seed_or_fixed_block.fixture.fixed_block
    elseif hasproperty(seed_or_fixed_block, :fixed_block)
        return seed_or_fixed_block.fixed_block
    end
    return seed_or_fixed_block
end

function _pqs_source_box_route_driver_white_lindsey_ham_preflight(
    seed_or_fixed_block;
    ham_bundle_adapter = nothing,
)
    fixed_block =
        isnothing(ham_bundle_adapter) ?
        _pqs_source_box_route_driver_white_lindsey_preflight_fixed_block(
            seed_or_fixed_block,
        ) :
        ham_bundle_adapter.fixed_block
    ordinary_qw_fixed_block_applicable =
        applicable(ordinary_cartesian_qiu_white_operators, fixed_block)
    nested_cartesian_fixed_block_applicable =
        applicable(nested_cartesian_operators, fixed_block)
    ida_builder_name_defined = isdefined(@__MODULE__, :ordinary_cartesian_ida_operators)
    ordinary_cartesian_ida_fixed_block_applicable =
        ida_builder_name_defined &&
        applicable(getfield(@__MODULE__, :ordinary_cartesian_ida_operators), fixed_block)
    bundle_object = isnothing(ham_bundle_adapter) ? fixed_block : ham_bundle_adapter
    basis_bundle_payload = cartesian_basis_bundle_payload(bundle_object; include_ham = true)
    basis_bundle_ham_payload_available = !isnothing(basis_bundle_payload.ham)
    pure_operator_payload_available =
        ordinary_qw_fixed_block_applicable ||
        nested_cartesian_fixed_block_applicable ||
        ordinary_cartesian_ida_fixed_block_applicable ||
        basis_bundle_ham_payload_available
    missing_builder =
        pure_operator_payload_available ?
        nothing :
        :missing_pure_low_order_fixed_block_density_density_interaction_builder
    status =
        basis_bundle_ham_payload_available && !isnothing(ham_bundle_adapter) ?
        :available_private_low_order_ham_bundle_adapter :
        pure_operator_payload_available ?
        :available_pure_low_order_operator_payload :
        :blocked_missing_pure_low_order_interaction_builder
    ham_operator_payload_status =
        pure_operator_payload_available ?
        :available_low_order_operator_payload :
        :pending_low_order_operator_payload
    ham_interaction_status =
        basis_bundle_ham_payload_available ?
        :available_low_order_density_density_interaction_matrix :
        :pending_low_order_density_density_interaction_matrix
    ham_bundle_export_status =
        basis_bundle_ham_payload_available ?
        :available_low_order_ham_bundle_payload :
        :pending_low_order_density_density_interaction_matrix

    return (;
        object_kind = :white_lindsey_low_order_ham_preflight,
        route_family = :white_lindsey_low_order,
        fixed_block_type_label = string(typeof(fixed_block)),
        parent_basis_type_label =
            hasproperty(fixed_block, :parent_basis) ?
            string(typeof(fixed_block.parent_basis)) :
            "unavailable",
        ordinary_qw_fixed_block_applicable,
        nested_cartesian_fixed_block_applicable,
        ordinary_cartesian_ida_builder_name_defined = ida_builder_name_defined,
        ordinary_cartesian_ida_fixed_block_applicable,
        basis_bundle_include_ham_checked = true,
        basis_bundle_ham_payload_available,
        basis_bundle_ham_payload_status =
            basis_bundle_ham_payload_available ?
            (
                isnothing(ham_bundle_adapter) ?
                :available :
                :available_private_writer_adapter
            ) :
            :absent_for_fixed_block,
        pure_operator_payload_available,
        status,
        required_builder_contract =
            :white_lindsey_low_order_fixed_block_density_density_builder,
        ham_operator_payload_status,
        ham_interaction_status,
        ham_bundle_export_status,
        missing_builder,
        supplement_required_paths_policy = :diagnostic_only_not_benchmark_route,
        full_ham_export_ready = basis_bundle_ham_payload_available,
        private_writer_adapter_used = !isnothing(ham_bundle_adapter),
        private_payload_candidate_status =
            isnothing(ham_bundle_adapter) ? nothing : ham_bundle_adapter.candidate.status,
    )
end

function _pqs_source_box_route_driver_materialization(
    report;
    materialize_route::Bool = false,
    save_basis_artifact::Bool = false,
    save_ham_artifact::Bool = false,
    basisfile::AbstractString = "cartesian_nesting_route_driver_basis_bundle.jld2",
    hamfile::AbstractString = "cartesian_nesting_route_driver_ham_bundle.jld2",
    white_lindsey_expansion = nothing,
    white_lindsey_Z = nothing,
)
    route_family = report.route_family
    route_configured_shellization_request =
        _cartesian_shellization_route_configured_request(report)
    route_configured_shellization_request_status =
        route_configured_shellization_request.status
    route_configured_system_classification =
        route_configured_shellization_request.system_classification
    route_configured_system_classification_status =
        route_configured_shellization_request.system_classification_status
    route_configured_bond_axis = route_configured_shellization_request.bond_axis
    route_configured_shellization_plan =
        _cartesian_shellization_route_planning_stub(route_configured_shellization_request)
    route_configured_shellization_plan_status =
        route_configured_shellization_plan.status
    route_configured_shellization_planning_status =
        route_configured_shellization_plan.planning_status
    route_configured_shellization_planning_family =
        route_configured_shellization_plan.planning_family
    route_configured_midpoint_slab_status =
        route_configured_shellization_plan.midpoint_slab_status

    if !materialize_route
        return (;
            object_kind = :cartesian_nesting_route_driver_materialization,
            route_family,
            private_development_only = true,
            materialize_route_requested = false,
            save_basis_artifact_requested = save_basis_artifact,
            save_ham_artifact_requested = save_ham_artifact,
            status = :not_requested_metadata_only,
            materialized_report = nothing,
            materialized_report_kind = nothing,
            route_configured_shellization_request,
            route_configured_shellization_request_available = true,
            route_configured_shellization_request_status,
            route_configured_system_classification,
            route_configured_system_classification_status,
            route_configured_bond_axis,
            route_configured_shellization_plan,
            route_configured_shellization_plan_available = true,
            route_configured_shellization_plan_status,
            route_configured_shellization_planning_status,
            route_configured_shellization_planning_family,
            route_configured_midpoint_slab_status,
            shellization_summary = nothing,
            shellization_summary_available = false,
            shellization_source =
                route_family == :white_lindsey_low_order ?
                :white_lindsey_one_center_seed_not_materialized :
                nothing,
            route_configured_shellization_consumed = false,
            materialized_shellization_stage = :not_checked_metadata_only,
            seed_materialization_status =
                route_family == :white_lindsey_low_order ?
                :not_requested_seed_materialization :
                :not_applicable,
            retained_dimension = report.retained_dimension,
            final_integral_weights_status = :not_checked_metadata_only,
            one_body_operator_status = :not_checked_metadata_only,
            basis_bundle_export_status = :not_requested,
            basis_artifact_status =
                save_basis_artifact ?
                :not_written_materialization_not_requested :
                :not_requested,
            basis_artifact_written = false,
            basisfile,
            basis_artifact_path = nothing,
            basis_export_blocker =
                save_basis_artifact ? :materialize_route_false : nothing,
            ham_preflight_status = :not_checked_metadata_only,
            ham_missing_builder = nothing,
            ham_operator_payload_status = :not_checked_metadata_only,
            ham_interaction_status = :not_checked_metadata_only,
            ham_bundle_export_status = :not_requested,
            ham_artifact_status =
                save_ham_artifact ?
                :not_written_materialization_not_requested :
                :not_requested,
            ham_artifact_written = false,
            hamfile,
            ham_export_blocker =
                save_ham_artifact ? :materialize_route_false : nothing,
            ham_preflight = nothing,
            pqs_materialization_status =
                route_family == :pqs_source_box ?
                :pending_source_box_retained_route :
                :not_applicable,
        )
    end

    if route_family == :white_lindsey_low_order
        materialized_report = _white_lindsey_low_order_materialized_seed_report()
        basis_export_status = :supported_basis_only_fixed_block
        shellization_summary = materialized_report.shellization_summary
        shellization_summary_available = materialized_report.shellization_summary_available
        shellization_source = materialized_report.shellization_source
        route_configured_shellization_consumed =
            materialized_report.route_configured_shellization_consumed
        materialized_shellization_stage = materialized_report.materialized_shellization_stage
        seed_materialization_status = materialized_report.seed_materialization_status
        ham_bundle_adapter = nothing
        if save_ham_artifact
            isnothing(white_lindsey_expansion) && throw(
                ArgumentError(
                    "White-Lindsey Ham artifact export requires explicit white_lindsey_expansion",
                ),
            )
            isnothing(white_lindsey_Z) && throw(
                ArgumentError("White-Lindsey Ham artifact export requires explicit white_lindsey_Z"),
            )
            ham_bundle_adapter =
                _white_lindsey_low_order_materialized_seed_ham_bundle_adapter(
                    materialized_report;
                    expansion = white_lindsey_expansion,
                    Z = white_lindsey_Z,
                )
        end
        ham_preflight =
            isnothing(ham_bundle_adapter) ?
            _pqs_source_box_route_driver_white_lindsey_ham_preflight(materialized_report) :
            _pqs_source_box_route_driver_white_lindsey_ham_preflight(
                materialized_report;
                ham_bundle_adapter = ham_bundle_adapter,
            )
        ham_operator_payload_status = ham_preflight.ham_operator_payload_status
        ham_interaction_status = ham_preflight.ham_interaction_status
        ham_bundle_export_status = ham_preflight.ham_bundle_export_status
        ham_export_blocker = ham_preflight.missing_builder
        basis_artifact_written = false
        basis_artifact_status =
            save_basis_artifact ? :written_basis_only_bundle : :not_requested
        if save_basis_artifact
            write_cartesian_basis_bundle_jld2(
                basisfile,
                materialized_report.fixture.fixed_block;
                include_ham = false,
                meta = (;
                    route_family,
                    route_kind = report.recipe_metadata.route_kind,
                    benchmark_role = report.recipe_metadata.benchmark_role,
                    materialized_report_kind = materialized_report.object_kind,
                    route_configured_shellization_request_status,
                    route_configured_system_classification,
                    route_configured_system_classification_status,
                    route_configured_bond_axis,
                    route_configured_shellization_plan_status,
                    route_configured_shellization_planning_status,
                    route_configured_shellization_planning_family,
                    route_configured_midpoint_slab_status,
                    shellization_summary_available,
                    shellization_source,
                    route_configured_shellization_consumed,
                    materialized_shellization_stage,
                    seed_materialization_status,
                    export_status = :basis_only,
                    basis_export_status,
                    ham_preflight_status = ham_preflight.status,
                    ham_missing_builder = ham_preflight.missing_builder,
                    ham_operator_payload_status,
                    ham_interaction_status,
                    ham_export_status = ham_bundle_export_status,
                    ham_export_blocker,
                    private_development_only = true,
                ),
            )
            basis_artifact_written = true
        end
        ham_artifact_written = false
        ham_artifact_status =
            save_ham_artifact ?
            :not_written_private_white_lindsey_ham_adapter_not_ready :
            :not_requested
        if save_ham_artifact && ham_preflight.full_ham_export_ready
            write_cartesian_basis_bundle_jld2(
                hamfile,
                ham_bundle_adapter;
                include_ham = true,
                meta = (;
                    route_family,
                    route_kind = report.recipe_metadata.route_kind,
                    benchmark_role = report.recipe_metadata.benchmark_role,
                    materialized_report_kind = materialized_report.object_kind,
                    route_configured_shellization_request_status,
                    route_configured_system_classification,
                    route_configured_system_classification_status,
                    route_configured_bond_axis,
                    route_configured_shellization_plan_status,
                    route_configured_shellization_planning_status,
                    route_configured_shellization_planning_family,
                    route_configured_midpoint_slab_status,
                    shellization_summary_available,
                    shellization_source,
                    route_configured_shellization_consumed,
                    materialized_shellization_stage,
                    seed_materialization_status,
                    export_status = :basis_and_ham,
                    basis_export_status,
                    ham_preflight_status = ham_preflight.status,
                    ham_missing_builder = ham_preflight.missing_builder,
                    ham_operator_payload_status,
                    ham_interaction_status,
                    ham_export_status = ham_bundle_export_status,
                    ham_export_blocker,
                    private_development_only = true,
                    private_writer_adapter =
                        :_WhiteLindseyLowOrderHamBundleAdapter,
                    ham_payload_candidate_status =
                        ham_bundle_adapter.candidate.status,
                ),
            )
            ham_artifact_written = true
            ham_artifact_status = :written_white_lindsey_low_order_ham_bundle
        end
        return (;
            object_kind = :cartesian_nesting_route_driver_materialization,
            route_family,
            private_development_only = true,
            materialize_route_requested = true,
            save_basis_artifact_requested = save_basis_artifact,
            save_ham_artifact_requested = save_ham_artifact,
            status = :materialized_seed_report_available,
            materialized_report,
            materialized_report_kind = materialized_report.object_kind,
            route_configured_shellization_request,
            route_configured_shellization_request_available = true,
            route_configured_shellization_request_status,
            route_configured_system_classification,
            route_configured_system_classification_status,
            route_configured_bond_axis,
            route_configured_shellization_plan,
            route_configured_shellization_plan_available = true,
            route_configured_shellization_plan_status,
            route_configured_shellization_planning_status,
            route_configured_shellization_planning_family,
            route_configured_midpoint_slab_status,
            shellization_summary,
            shellization_summary_available,
            shellization_source,
            route_configured_shellization_consumed,
            materialized_shellization_stage,
            seed_materialization_status,
            retained_dimension = materialized_report.retained_dimension,
            final_integral_weights_status =
                materialized_report.inventory.retained_basis_integral_weights_ready ?
                :available_retained_basis_integral_weights :
                :not_ready,
            one_body_operator_status =
                materialized_report.operator_inventory.all_finite ?
                :materialized_finite_one_body_inventory :
                :not_ready,
            basis_bundle_export_status = basis_export_status,
            basis_artifact_status,
            basis_artifact_written,
            basisfile,
            basis_artifact_path = basis_artifact_written ? basisfile : nothing,
            basis_export_blocker = nothing,
            ham_preflight_status = ham_preflight.status,
            ham_missing_builder = ham_preflight.missing_builder,
            ham_operator_payload_status,
            ham_interaction_status,
            ham_bundle_export_status,
            ham_artifact_status,
            ham_artifact_written,
            hamfile,
            ham_export_blocker,
            ham_preflight,
            pqs_materialization_status = :not_applicable,
        )
    end

    ham_export_blocker = :pending_source_box_retained_route
    return (;
        object_kind = :cartesian_nesting_route_driver_materialization,
        route_family,
        private_development_only = true,
        materialize_route_requested = true,
        save_basis_artifact_requested = save_basis_artifact,
        save_ham_artifact_requested = save_ham_artifact,
        status = :pending_source_box_retained_route,
        materialized_report = nothing,
        materialized_report_kind = nothing,
        route_configured_shellization_request,
        route_configured_shellization_request_available = true,
        route_configured_shellization_request_status,
        route_configured_system_classification,
        route_configured_system_classification_status,
        route_configured_bond_axis,
        route_configured_shellization_plan,
        route_configured_shellization_plan_available = true,
        route_configured_shellization_plan_status,
        route_configured_shellization_planning_status,
        route_configured_shellization_planning_family,
        route_configured_midpoint_slab_status,
        shellization_summary = nothing,
        shellization_summary_available = false,
        shellization_source = :pending_source_box_route_shellization,
        route_configured_shellization_consumed = false,
        materialized_shellization_stage = :pending_source_box_retained_route,
        seed_materialization_status = :not_applicable,
        retained_dimension = report.retained_dimension,
        final_integral_weights_status = :pending_final_ida_weights,
        one_body_operator_status = :pending_source_box_retained_blocks,
        basis_bundle_export_status = :pending_final_retained_basis,
        basis_artifact_status =
            save_basis_artifact ?
            :not_written_pending_final_retained_basis :
            :not_requested,
        basis_artifact_written = false,
        basisfile,
        basis_artifact_path = nothing,
        basis_export_blocker =
            save_basis_artifact ? :pending_final_retained_basis : nothing,
        ham_preflight_status = :not_applicable_to_pqs_source_box_route,
        ham_missing_builder = :pending_source_box_retained_route,
        ham_operator_payload_status = :pending_source_box_retained_operator_payload,
        ham_interaction_status = :pending_source_box_retained_density_density_blocks,
        ham_bundle_export_status = :pending_source_box_retained_route,
        ham_artifact_status =
            save_ham_artifact ?
            :not_written_pending_source_box_retained_route :
            :not_requested,
        ham_artifact_written = false,
        hamfile,
        ham_export_blocker,
        ham_preflight = nothing,
        pqs_materialization_status = :pending_source_box_retained_route,
    )
end


# Compatibility dry-run wrapper. This mirrors the executable driver stages,
# but returns a report directly for focused validation and tests.

function _pqs_source_box_route_driver_dry_run(;
    route_family = :pqs_source_box,
    route_kind, atom_symbols, nuclear_charges, atom_locations,
    radius, parent_axis_counts, map_backend,
    q, n_s, reference_spacing, tail_spacing, q_to_core_spacing_rule, core_spacing,
    probe_parent_axis_construction, parent_axis_probe_backend, parent_axis_probe_family,
    probe_raw_product_box_plans, raw_product_box_probe_backend,
    route_shape, product_body_rule, pqs_retained_rule, product_retained_rule,
    terms, pair_factor_normalization,
    support_dense_direct_allowed, reference_only_authorities,
    white_lindsey_route_shape = (:standard_cartesian_units, :low_order_comx_coarsening,),
    white_lindsey_mapping_rule = :standard_unit_backbone_mapping_family,
    white_lindsey_nesting_rule = :unit_box_low_order_comx_coarsening,
    white_lindsey_retained_rule = :low_order_unit_comx_retained_basis,
    white_lindsey_operator_rule = :low_order_unit_operator_blocks,
    white_lindsey_benchmark_role = :published_cartesian_baseline_for_pqs_comparison,
)
    system_inputs = (; atom_symbols, nuclear_charges, atom_locations,
        radius, parent_axis_counts, map_backend,)
    spacing_inputs = (; q, n_s, reference_spacing, tail_spacing,
        q_to_core_spacing_rule, core_spacing,)
    probe_inputs = (; probe_parent_axis_construction, parent_axis_probe_backend,
        parent_axis_probe_family, probe_raw_product_box_plans, raw_product_box_probe_backend,)
    source_box_recipe = (; route_shape, product_body_rule,
        pqs_retained_rule, product_retained_rule,
        support_dense_direct_allowed, reference_only_authorities,)
    white_lindsey_recipe = (; route_shape = white_lindsey_route_shape,
        mapping_rule = white_lindsey_mapping_rule, nesting_rule = white_lindsey_nesting_rule,
        retained_rule = white_lindsey_retained_rule, operator_rule = white_lindsey_operator_rule,
        benchmark_role = white_lindsey_benchmark_role,)
    route_recipe = (; route_family, route_kind, terms, pair_factor_normalization,
        source_box = source_box_recipe, white_lindsey = white_lindsey_recipe,)

    _pqs_source_box_route_driver_check_inputs(route_recipe)
    standard_setup = _pqs_source_box_route_driver_standard_setup(system_inputs, spacing_inputs)
    parent_axis = _pqs_source_box_route_driver_parent_axis(
        standard_setup, system_inputs, probe_inputs)
    route_axis_counts = _pqs_source_box_route_driver_route_axis_counts(
        standard_setup, parent_axis, system_inputs, route_recipe)
    system_metadata = _pqs_source_box_route_driver_system_metadata(
        standard_setup, route_axis_counts, system_inputs)
    route_skeleton = _pqs_source_box_route_driver_route_skeleton(
        route_axis_counts, spacing_inputs, route_recipe)
    raw_box = _pqs_source_box_route_driver_raw_box_probe(
        standard_setup, route_skeleton, parent_axis, route_axis_counts,
        probe_inputs, route_recipe)
    recipe_metadata = _pqs_source_box_route_driver_recipe_metadata(
        standard_setup, route_axis_counts, parent_axis, raw_box,
        spacing_inputs, probe_inputs, route_recipe)
    parent_description = _pqs_source_box_route_driver_parent_description(
        standard_setup, parent_axis, route_axis_counts, route_skeleton, raw_box)
    route_facts = _pqs_source_box_route_driver_route_facts(route_skeleton)
    contract = _pqs_source_box_route_driver_contract_metadata(route_recipe)
    diagnostics = _pqs_source_box_route_driver_diagnostics(
        standard_setup, parent_axis, route_axis_counts, route_skeleton, raw_box, contract)
    return _pqs_source_box_route_driver_report(
        standard_setup, parent_axis, route_axis_counts, raw_box,
        system_metadata, recipe_metadata, parent_description,
        route_skeleton, route_facts, contract, diagnostics)
end

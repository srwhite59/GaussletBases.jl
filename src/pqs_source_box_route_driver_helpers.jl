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


# Route skeletons.
#
# The default PQS/source-box skeleton delegates to the source-box route helper.
# The White-Lindsey route remains a low-order benchmark metadata skeleton.

function _pqs_source_box_route_driver_route_skeleton(
    route_axis_counts,
    spacing_inputs,
    route_recipe,
)
    metrics = CartesianContractedParentMetrics
    if route_recipe.route_family == :white_lindsey_low_order
        return _pqs_source_box_route_driver_white_lindsey_low_order_skeleton(
            route_axis_counts, spacing_inputs, route_recipe)
    end

    source_box_recipe = route_recipe.source_box
    skeleton = metrics._pqs_pqs_product_source_box_route_skeleton(
        ;
        q = spacing_inputs.q,
        parent_axis_counts = route_axis_counts.parent_axis_counts,
        route_shape = source_box_recipe.route_shape,
        product_body_rule = source_box_recipe.product_body_rule,
        pqs_retained_rule = source_box_recipe.pqs_retained_rule,
        product_retained_rule = source_box_recipe.product_retained_rule,
        pair_factor_normalization = route_recipe.pair_factor_normalization,
    )
    return merge(skeleton, (; route_family = route_recipe.route_family))
end

function _pqs_source_box_route_driver_white_lindsey_low_order_skeleton(
    route_axis_counts,
    spacing_inputs,
    route_recipe,
)
    low_order_recipe = route_recipe.white_lindsey
    counts = route_axis_counts.parent_axis_counts
    source_box = _pqs_route_driver_parent_source_box(counts)
    source_dimensions = _pqs_route_driver_axis_count_tuple(counts)
    source_dimension = isnothing(source_dimensions) ? nothing : prod(source_dimensions)

    retained_units = ((
        unit_key = :low_order_units,
        unit_role = :standard_box_units_with_white_lindsey_low_order_transform,
        retained_unit_kind = :low_order_unit_partition,
        source_family = :standard_cartesian_unit_boxes,
        source_box = source_box,
        source_dimensions = source_dimensions,
        source_dimension = source_dimension,
        retained_rule_kind = low_order_recipe.retained_rule,
        retained_rule_derivation = :unit_box_comx_coarsening_not_historical_split_rule,
        retained_range = nothing,
        retained_count = nothing,
        provenance_label = :white_lindsey_low_order_units,
        weight_semantics = :not_pqs_retained_weights,
    ),)
    pair_entries = ((
        pair_key = (:low_order_units, :low_order_units),
        pair_family = :white_lindsey_low_order,
        pair_kind = :low_order_unit_operator_pair,
        density_density_helper = :not_applicable_low_order_route,
        source_box_algorithmic_path = false,
        fallback_oracle_path = false,
        transpose_policy = :none,
        output_representation = :low_order_nested_cartesian_basis,
    ),)
    pair_family_counts = (
        pqs_pqs = 0,
        pqs_product = 0,
        product_pqs = 0,
        product_product = 0,
        white_lindsey_low_order = length(pair_entries),
    )
    helper_by_pair_family = (
        white_lindsey_low_order = :pending_white_lindsey_low_order_operator_builder,
    )
    pending_facts = (
        :materialized_standard_unit_partition,
        :low_order_comx_transform_per_unit,
        :low_order_operator_block_builder,
        :benchmark_validation_against_pqs_source_box_route,
    )

    return (;
        object_kind = :white_lindsey_low_order_route_skeleton,
        status = :published_benchmark_metadata_skeleton,
        route_family = route_recipe.route_family,
        route_shape = low_order_recipe.route_shape,
        retained_unit_order = (:low_order_units,),
        q = spacing_inputs.q,
        parent_axis_counts = counts,
        source_boxes = (low_order_units = source_box,),
        source_dimensions = (low_order_units = source_dimensions,),
        retained_units,
        retained_counts = (low_order_units = nothing,),
        ranges = (low_order_units = nothing,),
        retained_dimension = nothing,
        pair_entries,
        pair_family_counts,
        helper_by_pair_family,
        low_order_recipe,
        pending_facts,
        diagnostics = (
            source = :white_lindsey_low_order_route_skeleton,
            route_family = route_recipe.route_family,
            route_shape = low_order_recipe.route_shape,
            private_development_only = true,
            production_route = false,
            published_benchmark_route = true,
            literature_reference = :white_lindsey_nested_gausslet_basis_sets_jcp_2023,
            author_spelling = :lindsey,
            mapping_rule = low_order_recipe.mapping_rule,
            nesting_rule = low_order_recipe.nesting_rule,
            retained_rule = low_order_recipe.retained_rule,
            operator_rule = low_order_recipe.operator_rule,
            benchmark_role = low_order_recipe.benchmark_role,
            standard_unit_box_organization = true,
            atom_count_special_case_required = false,
            historical_split_rule_preserved = false,
            comx_coarsening_reference = true,
            one_dimensional_diagonal_basis_coarsening = true,
            source_box_first = false,
            source_box_algorithmic_path_true_for_every_pair = false,
            derived_retained_counts = false,
            derived_retained_ranges = false,
            derived_pair_inventory = true,
            pending_facts = pending_facts,
            pair_factor_normalization = route_recipe.pair_factor_normalization,
            pair_count = length(pair_entries),
            pair_family_counts = pair_family_counts,
            retained_dimension = nothing,
            retained_unit_count = length(retained_units),
            output_representation = :low_order_nested_cartesian_basis,
            four_index_galerkin_tensor = false,
            retained_pqs_weights_used = false,
            retained_weight_division_allowed = false,
            packet_adoption = false,
            fixed_block_routing = false,
            qwhamiltonian_consumes = false,
            public_default_consumes = false,
            mwg_ida_semantics_changed = false,
            repo_side_ray_id = false,
            ecp_terms_implemented = false,
            cr2_science_status_changed = false,
            shell_projection_used = false,
            lowdin_cleanup_used = false,
            support_local_shell_row_algorithm = false,
            support_coefficient_matrix_used = false,
        ),
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

function _pqs_source_box_route_driver_route_facts(route_skeleton)
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
        linear_algebra_plan = contract.linear_algebra_plan,
        stage_table = contract.stage_table,
        dry_run_validation = contract.dry_run_validation,
        diagnostics,
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

    _pqs_route_driver_print_named_tuple("parent_description", report.parent_description)
    _pqs_route_driver_print_named_tuple("source_boxes", report.source_boxes)

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

function _pqs_source_box_route_driver_save(
    report;
    save_artifact, save_tsv, outfile, tsvfile,
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

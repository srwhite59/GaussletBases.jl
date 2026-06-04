# Route-family skeleton helpers for `bin/cartesian_ham_builder.jl`.
#
# These functions keep route-family-specific skeleton construction out of the
# main driver helper file. They consume setup/axis-count utilities defined in
# `pqs_source_box_route_driver_helpers.jl` and feed report assembly there.


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

    retained_units = (
        _pqs_source_box_route_driver_unit_record(
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
        ),
    )
    unit_inventory = _pqs_source_box_route_driver_unit_inventory(retained_units)
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
        source_boxes = unit_inventory.source_boxes,
        source_dimensions = unit_inventory.source_dimensions,
        retained_units,
        retained_counts = unit_inventory.retained_counts,
        ranges = unit_inventory.ranges,
        retained_dimension = unit_inventory.retained_dimension,
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

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
    if route_recipe.route_kind in (
        :bond_aligned_diatomic_physical_gausslet_core_shell_pqs,
        :bond_aligned_diatomic_independent_pqs_source_box_core_shell,
    )
        return _pqs_source_box_route_driver_physical_gausslet_core_shell_skeleton(
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

function _pqs_source_box_route_driver_physical_gausslet_core_shell_skeleton(
    route_axis_counts,
    spacing_inputs,
    route_recipe,
)
    independent_target = route_recipe.route_kind ===
                         :bond_aligned_diatomic_independent_pqs_source_box_core_shell
    support_units = (:atom_contact_core, :shared_shell_1, :shared_shell_2)
    support_counts = (;
        atom_contact_core = 275,
        shared_shell_1 = 578,
        shared_shell_2 = 362,
    )
    retained_units =
        independent_target ? () :
        (
            _pqs_source_box_route_driver_unit_record(
            unit_key = :atom_contact_core,
            unit_role = :full_atom_contact_core_interiors,
            retained_unit_kind = :atom_contact_core,
            source_family = :reviewed_wl_qw_atom_contact_core_inventory,
            source_box = nothing,
            source_dimensions = nothing,
            source_dimension = 275,
            retained_rule_kind = :full_atom_contact_core_interiors,
            retained_rule_derivation = :pass_200_reviewed_wl_qw_inventory,
            retained_range = 1:251,
            retained_count = 251,
            provenance_label = :wl_qw_h2_r4_gausslet_only_core_inventory,
            weight_semantics = :not_ida_weights,
        ),
        _pqs_source_box_route_driver_unit_record(
            unit_key = :shared_shell_1,
            unit_role = :first_shared_molecular_shell,
            retained_unit_kind = :pqs_shared_shell_layer,
            source_family = :pqs_shared_molecular_shell_target,
            source_box = nothing,
            source_dimensions = nothing,
            source_dimension = 578,
            retained_rule_kind = :pqs_shared_shell_boundary_target,
            retained_rule_derivation = :pass_200_reviewed_wl_qw_inventory,
            retained_range = 252:349,
            retained_count = 98,
            provenance_label = :wl_qw_h2_r4_gausslet_only_shared_shell_1,
            weight_semantics = :not_ida_weights,
        ),
        _pqs_source_box_route_driver_unit_record(
            unit_key = :shared_shell_2,
            unit_role = :second_shared_molecular_shell,
            retained_unit_kind = :pqs_shared_shell_layer,
            source_family = :pqs_shared_molecular_shell_target,
            source_box = nothing,
            source_dimensions = nothing,
            source_dimension = 362,
            retained_rule_kind = :pqs_shared_shell_boundary_target,
            retained_rule_derivation = :pass_200_reviewed_wl_qw_inventory,
            retained_range = 350:463,
            retained_count = 114,
            provenance_label = :wl_qw_h2_r4_gausslet_only_shared_shell_2,
            weight_semantics = :not_ida_weights,
        ),
    )
    unit_inventory = _pqs_source_box_route_driver_unit_inventory(retained_units)
    retained_counts = unit_inventory.retained_counts
    target_status =
        independent_target ?
        :blocked_independent_pqs_source_box_target_readiness :
        :available_physical_gausslet_core_shell_target_inventory
    target_blocker =
        independent_target ?
        :missing_independent_pqs_atom_contact_core_retained_rule :
        nothing
    target_inventory = (;
        status = target_status,
        blocker = target_blocker,
        support_units,
        support_counts,
        retained_units = independent_target ? () : support_units,
        retained_counts,
        retained_order = independent_target ? () : support_units,
        expected_final_dimension = independent_target ? nothing : unit_inventory.retained_dimension,
        retained_atom_core_interiors = !independent_target,
        source_plan_role =
            independent_target ?
            :independent_pqs_source_box_construction :
            :atom_contact_core_plus_pqs_shared_shells,
        supplement_policy = :none,
        provenance =
            independent_target ?
            :pass_230_independent_pqs_recovery_audit :
            :pass_200_reviewed_wl_qw_inventory,
        source_backed_fixed_source_oracle_used = false,
        retained_transform_authority = :pqs_source_box_construction,
        primary_blocker = target_blocker,
        secondary_blocker =
            independent_target ?
            :missing_independent_pqs_shared_shell_2_retained_rule :
            nothing,
        source_plan_blocker =
            independent_target ?
            :missing_independent_pqs_physical_source_plan_materializer :
            nothing,
        support_plan =
            independent_target ?
            (;
                status = :blocked_independent_pqs_support_region_plan,
                blocker = :missing_independent_pqs_support_region_materializer,
                authority = :pqs_source_box_route_geometry_pending_materializer,
                counts_generated = false,
                counts_source = :target_constants_pending_support_region_materializer,
            ) :
            nothing,
    )
    pair_entries = ()
    pair_family_counts =
        independent_target ?
        (independent_pqs_target_readiness = 0,) :
        (physical_gausslet_target = 0,)
    pending_facts = (
        independent_target ?
        :independent_pqs_atom_contact_core_retained_rule :
        :physical_gausslet_source_plan_producer,
        independent_target ?
        :independent_pqs_shared_shell_2_retained_rule :
        :physical_gausslet_final_basis_builder,
        independent_target ?
        :independent_pqs_physical_source_plan_materializer :
        :physical_gausslet_h1_builder,
    )

    return (;
        object_kind =
            independent_target ?
            :pqs_diatomic_independent_source_box_core_shell_target_skeleton :
            :pqs_diatomic_physical_gausslet_core_shell_target_skeleton,
        status = target_status,
        route_family = route_recipe.route_family,
        route_kind = route_recipe.route_kind,
        route_shape = support_units,
        retained_unit_order = support_units,
        q = spacing_inputs.q,
        parent_axis_counts = route_axis_counts.parent_axis_counts,
        source_boxes = unit_inventory.source_boxes,
        source_dimensions = unit_inventory.source_dimensions,
        retained_units,
        retained_counts,
        ranges = unit_inventory.ranges,
        retained_dimension = unit_inventory.retained_dimension,
        pair_entries,
        pair_family_counts,
        helper_by_pair_family =
            independent_target ?
            (independent_pqs_target_readiness = :target_readiness_only_not_materialized,) :
            (physical_gausslet_target = :target_inventory_only_not_materialized,),
        physical_target_inventory = target_inventory,
        pending_facts,
        diagnostics = (
            source = :pqs_source_box_route_driver_physical_gausslet_core_shell_skeleton,
            route_family = route_recipe.route_family,
            route_kind = route_recipe.route_kind,
            private_development_only = true,
            production_route = false,
            source_box_first = true,
            target_inventory_only = true,
            source_plan_materialized = false,
            final_basis_materialized = false,
            h1_materialized = false,
            h1_j_materialized = false,
            rhf_materialized = false,
            retained_atom_core_interiors = !independent_target,
            supplement_policy = :none,
            expected_final_dimension =
                independent_target ? nothing : unit_inventory.retained_dimension,
            source_backed_fixed_source_oracle_used = false,
            retained_transform_authority = :pqs_source_box_construction,
        ),
    )
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

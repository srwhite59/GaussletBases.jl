# Route-family skeleton helpers for `bin/cartesian_ham_builder.jl`.
#
# These functions keep route-family-specific skeleton construction out of the
# main driver helper file. They consume setup/axis-count utilities defined in
# `pqs_source_box_route_driver_helpers.jl` and feed report assembly there.


# Route skeletons.
#
# The default PQS/source-box skeleton delegates to the source-box route helper.
# The White-Lindsey route remains a low-order benchmark metadata skeleton.

_pqs_source_box_route_driver_skeleton_source_dimension(box) =
    prod(length(getproperty(box, axis)) for axis in (:x, :y, :z))

_pqs_source_box_route_driver_skeleton_range(offset::Int, count::Int) =
    (offset + 1):(offset + count)

function _pqs_source_box_route_driver_skeleton_pair_family(left_kind::Symbol, right_kind::Symbol)
    left_kind == :pqs && right_kind == :pqs && return :pqs_pqs
    left_kind == :pqs && right_kind == :product_doside && return :pqs_product
    left_kind == :product_doside && right_kind == :pqs && return :product_pqs
    left_kind == :product_doside && right_kind == :product_doside && return :product_product
    throw(ArgumentError("unsupported source-box route unit pair $(left_kind), $(right_kind)"))
end

function _pqs_source_box_route_driver_skeleton_density_density_helper(
    pair_family::Symbol,
    pair_factor_normalization::Symbol,
)
    raw = pair_factor_normalization == :raw_weighted
    pair_family == :pqs_pqs && return raw ?
        :_pqs_pqs_source_box_raw_weighted_density_density_interaction_block :
        :_pqs_pqs_source_box_density_density_interaction_block
    pair_family in (:pqs_product, :product_pqs) && return raw ?
        :_pqs_product_source_box_raw_weighted_density_density_interaction_block :
        :_pqs_product_source_box_density_density_interaction_block
    pair_family == :product_product && return raw ?
        :_product_doside_source_box_raw_weighted_density_density_interaction_block :
        :_product_doside_source_box_density_density_interaction_block
    throw(ArgumentError("unsupported density-density pair family $(pair_family)"))
end

function _pqs_source_box_route_driver_skeleton_unit(;
    unit_key,
    unit_role,
    retained_unit_kind,
    source_family,
    source_box,
    source_dimensions,
    retained_rule_kind,
    retained_rule_derivation,
    retained_range,
    retained_count,
    provenance_label,
    weight_semantics,
)
    return _pqs_source_box_route_driver_unit_record(
        ;
        unit_key,
        unit_role,
        retained_unit_kind,
        source_family,
        source_box,
        source_dimensions,
        source_dimension =
            _pqs_source_box_route_driver_skeleton_source_dimension(source_box),
        retained_rule_kind,
        retained_rule_derivation,
        retained_range,
        retained_count,
        provenance_label,
        weight_semantics,
    )
end

function _pqs_source_box_route_driver_generic_source_box_skeleton(
    route_axis_counts,
    spacing_inputs,
    route_recipe,
)
    source_box_recipe = route_recipe.source_box
    q = Int(spacing_inputs.q)
    q >= 2 || throw(ArgumentError("PQS source-box route skeleton requires q >= 2"))
    route_shape = source_box_recipe.route_shape
    route_shape == (:pqs_left, :product, :pqs_right) ||
        throw(ArgumentError("generic PQS source-box skeleton only supports (:pqs_left, :product, :pqs_right)"))
    source_box_recipe.pqs_retained_rule == :boundary_comx_product_mode_selection ||
        throw(ArgumentError("generic PQS source-box skeleton requires boundary COMX product-mode selection"))
    source_box_recipe.product_retained_rule == :product_doside_retained_unit ||
        throw(ArgumentError("generic PQS source-box skeleton requires product/doside retained units"))
    route_recipe.pair_factor_normalization in (:density_normalized, :raw_weighted) ||
        throw(ArgumentError("generic PQS source-box skeleton requires a reviewed pair-factor normalization"))

    raw_counts = route_axis_counts.parent_axis_counts
    counts = (; x = Int(raw_counts.x), y = Int(raw_counts.y), z = Int(raw_counts.z))
    counts.x >= q && counts.y >= q && counts.z >= q ||
        throw(ArgumentError("PQS source-box route skeleton parent axis counts must be >= q"))
    source_box_recipe.product_body_rule == :centered_single_z_slab ||
        throw(ArgumentError("generic PQS source-box skeleton only keeps centered single-slab metadata"))

    product_z = div(counts.z - 1, 2) + 1
    source_boxes = (;
        pqs_left = (x = 1:q, y = 1:q, z = 1:q),
        pqs_right = (x = 1:q, y = 1:q, z = (counts.z - q + 1):counts.z),
        product = (x = 1:q, y = 1:q, z = product_z:product_z),
    )
    source_dimensions = (; pqs_left = (q, q, q), pqs_right = (q, q, q), product = (q, q, 1))
    retained_counts = (; pqs_left = q^3 - (q - 2)^3, pqs_right = q^3 - (q - 2)^3, product = q^2)
    ranges = (;
        pqs_left =
            _pqs_source_box_route_driver_skeleton_range(0, retained_counts.pqs_left),
        pqs_right =
            _pqs_source_box_route_driver_skeleton_range(retained_counts.pqs_left, retained_counts.pqs_right),
        product =
            _pqs_source_box_route_driver_skeleton_range(
                retained_counts.pqs_left + retained_counts.pqs_right,
                retained_counts.product,
            ),
    )
    retained_units = (
        _pqs_source_box_route_driver_skeleton_unit(
            unit_key = :pqs_left,
            unit_role = :left_mode_selected_raw_box_pqs_unit,
            retained_unit_kind = :pqs,
            source_family = :mode_selected_raw_product_box,
            source_box = source_boxes.pqs_left,
            source_dimensions = source_dimensions.pqs_left,
            retained_rule_kind = source_box_recipe.pqs_retained_rule,
            retained_rule_derivation = :boundary_comx_product_mode_selector,
            retained_range = ranges.pqs_left,
            retained_count = retained_counts.pqs_left,
            provenance_label = :pqs_left_source_modes,
            weight_semantics = :retained_columns_not_positive_quadrature_weights,
        ),
        _pqs_source_box_route_driver_skeleton_unit(
            unit_key = :pqs_right,
            unit_role = :right_mode_selected_raw_box_pqs_unit,
            retained_unit_kind = :pqs,
            source_family = :mode_selected_raw_product_box,
            source_box = source_boxes.pqs_right,
            source_dimensions = source_dimensions.pqs_right,
            retained_rule_kind = source_box_recipe.pqs_retained_rule,
            retained_rule_derivation = :boundary_comx_product_mode_selector,
            retained_range = ranges.pqs_right,
            retained_count = retained_counts.pqs_right,
            provenance_label = :pqs_right_source_modes,
            weight_semantics = :retained_columns_not_positive_quadrature_weights,
        ),
        _pqs_source_box_route_driver_skeleton_unit(
            unit_key = :product,
            unit_role = :middle_product_doside_slab_unit,
            retained_unit_kind = :product_doside,
            source_family = :product_doside,
            source_box = source_boxes.product,
            source_dimensions = source_dimensions.product,
            retained_rule_kind = source_box_recipe.product_retained_rule,
            retained_rule_derivation = :product_source_dimension,
            retained_range = ranges.product,
            retained_count = retained_counts.product,
            provenance_label = :product_doside_source_modes,
            weight_semantics = :product_source_weights_owned_by_source_box_helpers,
        ),
    )
    retained_unit_order = (:pqs_left, :pqs_right, :product)
    unit_by_key = Dict(unit.unit_key => unit for unit in retained_units)
    pair_entries = Tuple(begin
        left_key = retained_unit_order[i]
        right_key = retained_unit_order[j]
        pair_family = _pqs_source_box_route_driver_skeleton_pair_family(
            unit_by_key[left_key].retained_unit_kind,
            unit_by_key[right_key].retained_unit_kind,
        )
        (;
            pair_key = (left_key, right_key),
            pair_family,
            pair_kind =
                pair_family == :pqs_pqs ?
                :pqs_pqs_source_box_density_density_pair :
                pair_family == :pqs_product ?
                :pqs_product_source_box_density_density_pair :
                :product_doside_source_box_density_density_pair,
            density_density_helper =
                _pqs_source_box_route_driver_skeleton_density_density_helper(
                    pair_family,
                    route_recipe.pair_factor_normalization,
                ),
            source_box_algorithmic_path = true,
            fallback_oracle_path = false,
            transpose_policy = left_key == right_key ? :none :
                :lower_block_uses_transpose_when_pair_factors_are_symmetric,
            output_representation = :retained_two_index_density_density,
        )
    end for i in eachindex(retained_unit_order) for j in i:length(retained_unit_order))
    pair_inventory = _pqs_source_box_route_driver_pair_inventory(
        pair_entries;
        expected_families = (:pqs_pqs, :pqs_product, :product_pqs, :product_product),
    )
    unit_inventory = _pqs_source_box_route_driver_unit_inventory(retained_units)
    helper_by_pair_family = (;
        pqs_pqs =
            _pqs_source_box_route_driver_skeleton_density_density_helper(
                :pqs_pqs, route_recipe.pair_factor_normalization),
        pqs_product =
            _pqs_source_box_route_driver_skeleton_density_density_helper(
                :pqs_product, route_recipe.pair_factor_normalization),
        product_pqs = :transpose_of_pqs_product_helper_for_lower_blocks_only,
        product_product =
            _pqs_source_box_route_driver_skeleton_density_density_helper(
                :product_product, route_recipe.pair_factor_normalization),
    )
    return (;
        route_family = route_recipe.route_family,
        route_kind = route_recipe.route_kind,
        route_shape,
        retained_unit_order,
        q,
        parent_axis_counts = counts,
        source_boxes = unit_inventory.source_boxes,
        source_dimensions = unit_inventory.source_dimensions,
        retained_units,
        retained_counts = unit_inventory.retained_counts,
        ranges = unit_inventory.ranges,
        retained_dimension = unit_inventory.retained_dimension,
        pair_entries = pair_inventory.pair_entries,
        pair_family_counts = pair_inventory.pair_family_counts,
        helper_by_pair_family,
    )
end

function _pqs_source_box_route_driver_route_skeleton(
    route_axis_counts,
    spacing_inputs,
    route_recipe,
)
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

    return _pqs_source_box_route_driver_generic_source_box_skeleton(
        route_axis_counts, spacing_inputs, route_recipe)
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
        :missing_independent_pqs_physical_source_plan_materializer :
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
            nothing,
        source_plan_blocker =
            independent_target ?
            :missing_independent_pqs_physical_source_plan_materializer :
            nothing,
    )
    pair_entries = ()
    pair_family_counts =
        independent_target ?
        (independent_pqs_target_readiness = 0,) :
        (physical_gausslet_target = 0,)
    return (;
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
    return (;
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
    )
end

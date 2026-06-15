# Private support for `bin/cartesian_ham_builder.jl`.
#
# Keep route bookkeeping here so the executable driver can stay human-facing:
# editable defaults, overrides, visible stages, print, and save.


# Small local utility helpers.

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

function _pqs_route_driver_type_label(object)
    isnothing(object) && return "unavailable"
    return string(nameof(typeof(object)))
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
    q = Int(spacing_inputs.q)
    n_s = Int(spacing_inputs.n_s)
    core_cube_side = 2 * n_s + 1
    core_spacing = spacing_inputs.core_spacing
    parent_radius = max(Int(ceil(Float64(system_inputs.radius))), n_s)
    parent_axis = -parent_radius:parent_radius
    return (;
        nuclear_charges = Tuple(system_inputs.nuclear_charges),
        atom_locations = Tuple(system_inputs.atom_locations),
        q,
        radius = system_inputs.radius,
        reference_spacing = spacing_inputs.reference_spacing,
        tail_spacing = spacing_inputs.tail_spacing,
        q_to_core_spacing_rule = spacing_inputs.q_to_core_spacing_rule,
        core_spacing,
        xmax_parallel = get(spacing_inputs, :xmax_parallel, nothing),
        xmax_transverse = get(spacing_inputs, :xmax_transverse, nothing),
        n_s,
        core_cube_side,
        core_cube_side_rule = :two_n_s_plus_one,
        parent_box_rule = :radius_index_box,
        parent_box = (x = parent_axis, y = parent_axis, z = parent_axis),
        mapping_s = core_spacing,
        mapping_s_by_atom = ntuple(_ -> core_spacing, length(system_inputs.nuclear_charges)),
    )
end

function _cartesian_center_list_axis_records(
    center_table,
    axis_index::Int,
    core_spacing::Real,
)
    records = Dict{Float64,Tuple{Float64,Float64}}()
    for center in center_table
        coordinate = Float64(center.location[axis_index])
        charge = Float64(center.nuclear_charge)
        charge > 0.0 || throw(ArgumentError("Cartesian parent centers require positive nuclear charges"))
        core_range = sqrt(Float64(core_spacing) / charge)
        target_spacing = Float64(core_spacing)
        records[coordinate] =
            haskey(records, coordinate) ?
            (min(records[coordinate][1], core_range), min(records[coordinate][2], target_spacing)) :
            (core_range, target_spacing)
    end
    coordinates = sort(collect(keys(records)))
    return coordinates,
           [records[coordinate][1] for coordinate in coordinates],
           [records[coordinate][2] for coordinate in coordinates]
end

function _cartesian_center_list_axis_extent(
    center_table,
    axis_index::Int,
    standard_setup,
)
    coordinates = Float64[Float64(center.location[axis_index]) for center in center_table]
    active_axis = maximum(coordinates) - minimum(coordinates) > 1.0e-12
    requested =
        active_axis ?
        get(standard_setup, :xmax_parallel, nothing) :
        get(standard_setup, :xmax_transverse, nothing)
    return Float64(something(requested, standard_setup.radius))
end

function _cartesian_center_list_parent_axis(
    center_table,
    standard_setup,
    parent_inputs,
    axis_index::Int,
)
    isnothing(standard_setup.core_spacing) &&
        throw(ArgumentError("center-list Cartesian parent construction requires core_spacing"))
    coordinates, core_ranges, target_spacings =
        _cartesian_center_list_axis_records(
            center_table,
            axis_index,
            standard_setup.core_spacing,
        )
    mapping = fit_combined_invsqrt_mapping(
        centers = coordinates,
        core_ranges = core_ranges,
        target_spacings = target_spacings,
        tail_spacing = standard_setup.tail_spacing,
    )
    extent = _cartesian_center_list_axis_extent(center_table, axis_index, standard_setup)
    count = _qwrg_mapped_odd_count_for_extent(
        mapping,
        extent;
        reference_spacing = standard_setup.reference_spacing,
    )
    family = get(parent_inputs, :parent_axis_family,
        :G10)
    axis = build_basis(MappedUniformBasisSpec(
        family;
        count,
        mapping,
        reference_spacing = standard_setup.reference_spacing,
    ))
    return axis, count
end

function _cartesian_center_list_parent_basis_object(
    center_table,
    classification,
    standard_setup,
    parent_inputs,
)
    classification.system_classification == :one_center && return nothing
    axes_and_counts = ntuple(
        axis_index ->
            _cartesian_center_list_parent_axis(
                center_table,
                standard_setup,
                parent_inputs,
                axis_index,
            ),
        3,
    )
    axes = (x = axes_and_counts[1][1], y = axes_and_counts[2][1],
        z = axes_and_counts[3][1])
    counts = (x = axes_and_counts[1][2], y = axes_and_counts[2][2],
        z = axes_and_counts[3][2])
    return CartesianParentGaussletBases.CartesianParentGaussletBasis3D(
        axes;
        metadata = (
            basis_family = :center_list_mapped_uniform_cartesian_product,
            axis_basis_helper = :_cartesian_center_list_parent_basis_object,
            reference_spacing = standard_setup.reference_spacing,
            core_spacing = standard_setup.core_spacing,
            atom_symbols = Tuple(center.atom_symbol for center in center_table),
            nuclear_charges = Tuple(center.nuclear_charge for center in center_table),
            atom_locations = Tuple(center.location for center in center_table),
            axis_counts = counts,
        ),
    )
end

function _cartesian_axis_bundle_from_parent_basis_object(parent_basis_object, parent_inputs)
    axes = CartesianParentGaussletBases.parent_axes(parent_basis_object)
    expansion = coulomb_gaussian_expansion(doacc = false)
    backend =
        hasproperty(parent_inputs, :parent_axis_bundle_backend) ?
        parent_inputs.parent_axis_bundle_backend :
        :pgdg_localized_experimental
    function _axis_bundle(axis)
        return _mapped_ordinary_gausslet_1d_bundle(
            axis;
            exponents = expansion.exponents,
            center = 0.0,
            backend,
        )
    end
    return _CartesianNestedAxisBundles3D(
        _axis_bundle(axes.x),
        _axis_bundle(axes.y),
        _axis_bundle(axes.z),
    )
end

function _cartesian_parent_symbol_tuple(atom_symbols)
    atom_symbols isa AbstractString && return (atom_symbols,)
    atom_symbols isa Symbol && return (atom_symbols,)
    return Tuple(atom_symbols)
end

function _cartesian_parent_location_tuple(location)
    length(location) == 3 || throw(
        ArgumentError("cartesian parent center locations must be three-dimensional"),
    )
    return (Float64(location[1]), Float64(location[2]), Float64(location[3]))
end

function _cartesian_parent_center_table(system, standard_setup)
    atom_symbols = _cartesian_parent_symbol_tuple(system.atom_symbols)
    nuclear_charges = Tuple(standard_setup.nuclear_charges)
    atom_locations = Tuple(standard_setup.atom_locations)
    length(atom_symbols) == length(nuclear_charges) == length(atom_locations) || throw(
        ArgumentError("cartesian parent requires matching atom symbols, charges, and locations"),
    )
    return Tuple(
        (
            center_index = index,
            center_key = Symbol(:center_, index),
            atom_symbol = atom_symbols[index],
            nuclear_charge = nuclear_charges[index],
            location = _cartesian_parent_location_tuple(atom_locations[index]),
        ) for index in eachindex(atom_symbols)
    )
end

function _cartesian_parent_axis_metadata(center_table; atol::Float64 = 1.0e-12)
    isempty(center_table) && throw(
        ArgumentError("cartesian parent requires at least one center"),
    )
    reference_location = first(center_table).location
    active_axes = Tuple(
        axis for (axis, axis_index) in zip((:x, :y, :z), 1:3)
        if any(
            center ->
                abs(center.location[axis_index] - reference_location[axis_index]) > atol,
            center_table,
        )
    )
    axis_aligned = length(active_axes) <= 1
    chain_axis = axis_aligned && !isempty(active_axes) ? only(active_axes) : nothing
    axis_index_by_symbol = (x = 1, y = 2, z = 3)
    center_axis_coordinates =
        isnothing(chain_axis) ?
        () :
        Tuple(center.location[getproperty(axis_index_by_symbol, chain_axis)] for center in center_table)

    return (
        active_axes,
        axis_aligned,
        bond_axis = length(center_table) == 2 ? chain_axis : nothing,
        chain_axis,
        center_axis_coordinates,
        center_axis_order = isempty(center_axis_coordinates) ?
            () :
            Tuple(sortperm(collect(center_axis_coordinates))),
    )
end

function _cartesian_parent_system_classification(center_table, axis_metadata)
    atom_count = length(center_table)
    if atom_count == 1
        return (
            system_classification = :one_center,
            bond_axis = nothing,
            chain_axis = nothing,
        )
    elseif atom_count == 2 && !isnothing(axis_metadata.bond_axis)
        return (
            system_classification = :bond_aligned_diatomic,
            bond_axis = axis_metadata.bond_axis,
            chain_axis = axis_metadata.chain_axis,
        )
    elseif atom_count == 2
        return (
            system_classification = :non_axis_aligned_diatomic,
            bond_axis = nothing,
            chain_axis = nothing,
        )
    elseif !isnothing(axis_metadata.chain_axis)
        return (
            system_classification = :axis_aligned_chain_metadata_only,
            bond_axis = nothing,
            chain_axis = axis_metadata.chain_axis,
        )
    end
    return (
        system_classification = :unsupported_general_multi_atom,
        bond_axis = nothing,
        chain_axis = nothing,
    )
end

function _cartesian_one_center_parent_mapping(center, standard_setup, parent_inputs)
    rule = get(parent_inputs, :parent_mapping_rule, :identity_mapping)
    rule === :identity_mapping && return IdentityMapping(), :IdentityMapping
    rule === :white_lindsey_atomic_mapping || throw(
        ArgumentError("unsupported one-center parent_mapping_rule $(repr(rule))"),
    )
    d = get(parent_inputs, :parent_mapping_d, nothing)
    isnothing(d) && throw(
        ArgumentError("white_lindsey_atomic_mapping parent mapping requires parent_mapping_d"),
    )
    Z = something(get(parent_inputs, :parent_mapping_Z, nothing), center.nuclear_charge)
    tail_spacing =
        something(get(parent_inputs, :parent_mapping_tail_spacing, nothing),
            standard_setup.tail_spacing)
    return white_lindsey_atomic_mapping(;
        Z = Float64(Z),
        d = Float64(d),
        tail_spacing = Float64(tail_spacing),
    ), :white_lindsey_atomic_mapping
end

function _cartesian_one_center_parent_basis_object(
    center_table,
    classification,
    standard_setup,
    route_axis_counts,
    parent_inputs,
)
    classification.system_classification == :one_center || return nothing

    counts = route_axis_counts.parent_axis_counts
    isnothing(counts) &&
        throw(ArgumentError("one-center Cartesian parent construction requires axis counts"))

    center = first(center_table)
    origin_centered = all(abs.(center.location) .<= 1.0e-12)
    origin_centered || throw(
        ArgumentError("one-center Cartesian parent construction expects an origin-centered atom"),
    )

    family = get(parent_inputs, :parent_axis_family, :G10)
    mapping, mapping_label =
        _cartesian_one_center_parent_mapping(center, standard_setup, parent_inputs)
    axis_cache = Dict{Int, Any}()
    function _axis_for_count(count)
        axis_count = Int(count)
        get!(axis_cache, axis_count) do
            build_basis(MappedUniformBasisSpec(
                family;
                count = axis_count,
                mapping,
                reference_spacing = standard_setup.reference_spacing,
            ))
        end
    end

    return CartesianParentGaussletBases.CartesianParentGaussletBasis3D(
        (;
            x = _axis_for_count(counts.x),
            y = _axis_for_count(counts.y),
            z = _axis_for_count(counts.z),
        );
        metadata = (
            basis_family = :one_center_mapped_uniform_cartesian_product,
            axis_basis_helper = :_cartesian_one_center_parent_basis_object,
            axis_family = Symbol(family),
            mapping = mapping_label,
            reference_spacing = standard_setup.reference_spacing,
            atom_symbol = center.atom_symbol,
            nuclear_charge = center.nuclear_charge,
            atom_location = center.location,
            axis_counts = counts,
        ),
    )
end

function _cartesian_parent_object_carry(
    center_table,
    classification,
    standard_setup,
    parent_axis,
    route_axis_counts,
    parent_inputs,
)
    one_center_parent_basis =
        _cartesian_one_center_parent_basis_object(
            center_table, classification, standard_setup,
            route_axis_counts, parent_inputs)
    center_list_parent_basis =
        isnothing(one_center_parent_basis) ?
        _cartesian_center_list_parent_basis_object(
            center_table,
            classification,
            standard_setup,
            parent_inputs,
        ) :
        nothing
    parent_basis_object =
        isnothing(one_center_parent_basis) ?
        center_list_parent_basis :
        one_center_parent_basis
    isnothing(parent_basis_object) &&
        throw(ArgumentError("Cartesian parent construction did not produce a parent basis object"))

    axis_bundle_object =
        _cartesian_axis_bundle_from_parent_basis_object(
            parent_basis_object,
            parent_inputs,
        )

    return (;
        parent_basis_object,
        parent_axis_bundle_object = axis_bundle_object,
        parent_axis_counts =
            _pqs_route_driver_axis_counts(
                CartesianParentGaussletBases.parent_axis_counts(parent_basis_object),
            ),
    )
end

# Parent-axis construction.

function _pqs_source_box_route_driver_parent_axis(
    standard_setup,
    system_inputs,
)
    requested_counts = _pqs_route_driver_axis_counts(system_inputs.parent_axis_counts)
    axis_counts = isnothing(requested_counts) ? (;
        x = length(standard_setup.parent_box.x),
        y = length(standard_setup.parent_box.y),
        z = length(standard_setup.parent_box.z),
    ) : requested_counts
    geometry = (;
        bond_axis = hasproperty(system_inputs, :bond_axis) ? system_inputs.bond_axis : nothing,
        bond_length =
            hasproperty(system_inputs, :bond_length) ? system_inputs.bond_length : nothing,
    )

    return (;
        parent_axis_counts = axis_counts,
        geometry,
    )
end


# Route axis counts.

function _pqs_source_box_route_driver_route_axis_counts(
    standard_setup,
    parent_axis,
)
    counts = parent_axis.parent_axis_counts
    counts.x >= standard_setup.q &&
        counts.y >= standard_setup.q &&
        counts.z >= standard_setup.q ||
        throw(ArgumentError("parent axis counts must be at least q on every axis"))

    return (;
        parent_axis_counts = counts,
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
        bond_axis = get(system_inputs, :bond_axis, nothing),
        bond_length = get(system_inputs, :bond_length, nothing),
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

# Metadata and parent description.

function _pqs_source_box_route_driver_parent_contract(parent)
    return (;
        object_kind = :cartesian_route_parent_contract,
        status = parent.status,
        atom_count = parent.atom_count,
        center_count = parent.center_count,
        atom_symbols = parent.atom_symbols,
        nuclear_charges = parent.nuclear_charges,
        atom_locations = parent.atom_locations,
        center_table = parent.center_table,
        center_axis_metadata = parent.center_axis_metadata,
        system_classification = parent.system_classification,
        system_classification_status = parent.system_classification_status,
        bond_axis = parent.bond_axis,
        chain_axis = parent.chain_axis,
        parent_axis_counts = parent.axis_counts,
        parent_axis_counts_source = parent.axis_counts_source,
        parent_axis_counts_status = parent.axis_counts_status,
        parent_box = parent.physical_box,
        parent_box_rule = parent.physical_box_rule,
        parent_materialization_plan = parent.parent_materialization_plan,
        parent_materialization_plan_status =
            parent.parent_materialization_plan_status,
        parent_materialization_planning_family =
            parent.parent_materialization_planning_family,
        parent_materialization_blocker = parent.parent_materialization_blocker,
        parent_basis_object_available = parent.parent_basis_object_available,
        parent_qw_basis_object_available =
            parent.parent_qw_basis_object_available,
        parent_axis_bundle_object_available =
            parent.parent_axis_bundle_object_available,
        parent_basis_object_type_label = parent.parent_basis_object_type_label,
        parent_qw_basis_object_type_label =
            parent.parent_qw_basis_object_type_label,
        parent_axis_bundle_object_type_label =
            parent.parent_axis_bundle_object_type_label,
        parent_basis_materialization_status =
            parent.parent_basis_materialization_status,
        parent_basis_materialization = parent.parent_basis_materialization,
        parent_basis_materialized = parent.parent_basis_materialized,
        parent_axis_metadata_constructed =
            parent.parent_axis_metadata_constructed,
        axis_bundle_materialized = parent.axis_bundle_materialized,
        diagnostics = (
            source = :cartesian_parent,
            private_development_only = true,
            report_parent_contract = true,
            parent_contract_driven_downstream_metadata = true,
            parent_basis_materialized = parent.parent_basis_materialized,
            axis_bundle_materialized = parent.axis_bundle_materialized,
            public_default_behavior_changed = false,
        ),
    )
end

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
        xmax_parallel = get(spacing_inputs, :xmax_parallel, nothing),
        xmax_transverse = get(spacing_inputs, :xmax_transverse, nothing),
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
        parent_mapping_rule = get(probe_inputs, :parent_mapping_rule, :identity_mapping),
        parent_mapping_Z = get(probe_inputs, :parent_mapping_Z, nothing),
        parent_mapping_d = get(probe_inputs, :parent_mapping_d, nothing),
        parent_mapping_tail_spacing =
            get(probe_inputs, :parent_mapping_tail_spacing, nothing),
        retained_atom_core_interiors =
            get(route_recipe, :retained_atom_core_interiors, nothing),
        supplement_policy = get(route_recipe, :supplement_policy, nothing),
        run_final_basis = get(route_recipe, :run_final_basis, true),
        run_h1 = get(route_recipe, :run_h1, true),
        run_h1_j = get(route_recipe, :run_h1_j, true),
        run_private_rhf =
            get(get(route_recipe, :private_rhf_inputs, (;)), :run_private_rhf, false),
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
    parent_contract,
    route_skeleton,
    raw_box,
)
    source_box_route = route_skeleton.route_family == :pqs_source_box
    return (;
        status = :described_not_constructed,
        route_family = route_skeleton.route_family,
        standard_setup,
        route_axis_counts,
        parent_contract_status = parent_contract.status,
        parent_contract_object_kind = parent_contract.object_kind,
        center_count = parent_contract.center_count,
        system_classification = parent_contract.system_classification,
        system_classification_status =
            parent_contract.system_classification_status,
        bond_axis = parent_contract.bond_axis,
        chain_axis = parent_contract.chain_axis,
        parent_materialization_plan_status =
            parent_contract.parent_materialization_plan_status,
        parent_materialization_planning_family =
            parent_contract.parent_materialization_planning_family,
        parent_materialization_blocker =
            parent_contract.parent_materialization_blocker,
        parent_basis_object_available =
            parent_contract.parent_basis_object_available,
        parent_axis_bundle_object_available =
            parent_contract.parent_axis_bundle_object_available,
        parent_basis_materialization_status =
            parent_contract.parent_basis_materialization_status,
        parent_basis_materialized = parent_contract.parent_basis_materialized,
        axis_bundle_materialized = parent_contract.axis_bundle_materialized,
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
        :_pqs_source_box_route_driver_generic_source_box_skeleton :
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
    parent,
    parent_axis,
    route_axis_counts,
    raw_box,
    system_metadata,
    recipe_metadata,
    parent_contract,
    parent_description,
    route_skeleton,
    route_facts,
    contract,
    diagnostics,
    low_order_route_summary,
)
    route_materializer_payload =
        _pqs_source_box_route_driver_materializer_payload(parent)
    return (;
        object_kind = :cartesian_nesting_route_driver_skeleton_report,
        generated_at = Base.Libc.strftime("%Y-%m-%dT%H:%M:%S", time()),
        route_family = route_skeleton.route_family,
        standard_setup,
        route_axis_counts,
        raw_product_box_probe = raw_box.raw_product_box_probe,
        system_metadata,
        recipe_metadata,
        parent_contract,
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
        low_order_route_summary,
        diagnostics,
        route_materializer_payload,
    )
end

function _pqs_source_box_route_driver_materializer_payload(parent)
    parent_axis_probe =
        hasproperty(parent, :parent_axis_probe) ? parent.parent_axis_probe : nothing
    parent_basis_object =
        hasproperty(parent, :parent_basis_object) ? parent.parent_basis_object : nothing
    parent_qw_basis_object =
        hasproperty(parent, :parent_qw_basis_object) ?
        parent.parent_qw_basis_object :
        nothing
    parent_axis_bundle_object =
        hasproperty(parent, :parent_axis_bundle_object) ?
        parent.parent_axis_bundle_object :
        nothing
    axis_bundle_backend =
        !isnothing(parent_axis_probe) &&
        hasproperty(parent_axis_probe, :gausslet_backend) ?
        parent_axis_probe.gausslet_backend :
        nothing

    return (;
        object_kind = :cartesian_route_materializer_transient_payload,
        private_development_only = true,
        transient_only = true,
        durable_report_serialization = :sanitize_before_save,
        source = :parent_axis_probe_object_carry,
        parent_basis_object,
        parent_qw_basis_object,
        parent_axis_bundle_object,
        parent_basis_object_available = !isnothing(parent_basis_object),
        parent_qw_basis_object_available = !isnothing(parent_qw_basis_object),
        parent_axis_bundle_object_available = !isnothing(parent_axis_bundle_object),
        parent_basis_object_type_label =
            _pqs_route_driver_type_label(parent_basis_object),
        parent_qw_basis_object_type_label =
            _pqs_route_driver_type_label(parent_qw_basis_object),
        parent_axis_bundle_object_type_label =
            _pqs_route_driver_type_label(parent_axis_bundle_object),
        axis_bundle_backend,
        axis_bundle_backend_available = !isnothing(axis_bundle_backend),
    )
end

# High-level driver facade.

function cartesian_system(system_inputs)
    return system_inputs
end

function cartesian_recipe(route_inputs)
    if hasproperty(route_inputs, :source_box) &&
       hasproperty(route_inputs, :white_lindsey)
        route_recipe = route_inputs
    else
        run_h1 = get(route_inputs, :run_h1, true)
        run_h1_j = get(route_inputs, :run_h1_j, true)
        private_rhf_inputs =
            get(route_inputs, :private_rhf_inputs, (; run_private_rhf = false))
        run_private_rhf = get(private_rhf_inputs, :run_private_rhf, false)
        run_final_basis = get(route_inputs, :run_final_basis, nothing)
        source_box_recipe = (;
            route_shape = route_inputs.route_shape,
            product_body_rule = route_inputs.product_body_rule,
            pqs_retained_rule = route_inputs.pqs_retained_rule,
            product_retained_rule = route_inputs.product_retained_rule,
        )
        white_lindsey_recipe = (;
            route_shape = route_inputs.white_lindsey_route_shape,
            mapping_rule = route_inputs.white_lindsey_mapping_rule,
            nesting_rule = route_inputs.white_lindsey_nesting_rule,
            retained_rule = route_inputs.white_lindsey_retained_rule,
            operator_rule = route_inputs.white_lindsey_operator_rule,
        )
        route_recipe = (;
            route_family = route_inputs.route_family,
            route_kind = route_inputs.route_kind,
            terms = route_inputs.terms,
            pair_factor_normalization = route_inputs.pair_factor_normalization,
            source_box = source_box_recipe,
            white_lindsey = white_lindsey_recipe,
            supplement_policy = get(route_inputs, :supplement_policy, nothing),
            run_final_basis =
                isnothing(run_final_basis) ?
                (run_h1 || run_h1_j || run_private_rhf) :
                run_final_basis,
            run_h1,
            run_h1_j,
            private_rhf_inputs,
        )
    end

    _pqs_source_box_route_driver_check_inputs(route_recipe)
    return route_recipe
end

function cartesian_parent(system, spacing_inputs, parent_inputs, recipe)
    standard_setup =
        _pqs_source_box_route_driver_standard_setup(system, spacing_inputs)
    center_table = _cartesian_parent_center_table(system, standard_setup)
    center_axis_metadata = _cartesian_parent_axis_metadata(center_table)
    classification =
        _cartesian_parent_system_classification(center_table, center_axis_metadata)
    parent_axis =
        _pqs_source_box_route_driver_parent_axis(
            standard_setup, system)
    route_axis_counts =
        _pqs_source_box_route_driver_route_axis_counts(
            standard_setup, parent_axis)
    object_carry = _cartesian_parent_object_carry(
        center_table, classification, standard_setup,
        parent_axis, route_axis_counts, parent_inputs)
    parent_axis_counts =
        isnothing(object_carry.parent_axis_counts) ?
        route_axis_counts.parent_axis_counts :
        object_carry.parent_axis_counts

    return (;
        system,
        spacing_inputs,
        parent_inputs,
        standard_setup,
        atom_count = length(center_table),
        atom_symbols = Tuple(center.atom_symbol for center in center_table),
        nuclear_charges = Tuple(center.nuclear_charge for center in center_table),
        atom_locations = Tuple(center.location for center in center_table),
        center_table,
        system_classification = classification.system_classification,
        bond_axis = classification.bond_axis,
        axis_counts = parent_axis_counts,
        physical_box = standard_setup.parent_box,
        parent_basis_object = object_carry.parent_basis_object,
        parent_axis_bundle_object = object_carry.parent_axis_bundle_object,
    )
end

function _pqs_source_box_route_driver_terminal_parent_axes(parent)
    counts = _pqs_route_driver_axis_count_tuple(parent.axis_counts)
    isnothing(counts) &&
        throw(ArgumentError("terminal shellification requires parent axis counts"))
    all(>(0), counts) ||
        throw(ArgumentError("terminal shellification requires positive parent axis counts"))
    return (collect(1:counts[1]), collect(1:counts[2]), collect(1:counts[3]))
end

function _pqs_source_box_route_driver_axis_index(axis::Symbol)
    axis == :x && return 1
    axis == :y && return 2
    axis == :z && return 3
    throw(ArgumentError("terminal shellification requires bond_axis = :x, :y, or :z"))
end

function _pqs_source_box_route_driver_terminal_center_index(count::Int)
    count > 0 || throw(ArgumentError("terminal shellification axis count must be positive"))
    return cld(count, 2)
end

function _pqs_source_box_route_driver_terminal_nuclear_index_positions(parent)
    counts = _pqs_route_driver_axis_count_tuple(parent.axis_counts)
    isnothing(counts) &&
        throw(ArgumentError("terminal shellification requires parent axis counts"))
    core_side = parent.standard_setup.core_cube_side
    radius = div(core_side - 1, 2)
    center = ntuple(
        axis -> _pqs_source_box_route_driver_terminal_center_index(counts[axis]),
        3,
    )
    if parent.system_classification == :one_center
        return (center,)
    end
    parent.system_classification == :bond_aligned_diatomic || throw(
        ArgumentError("terminal shellification supports one-center or bond-aligned diatomic parents"),
    )
    bond_axis_index = _pqs_source_box_route_driver_axis_index(parent.bond_axis)
    low_center = 1 + radius
    high_center = counts[bond_axis_index] - radius
    low_center > 0 && high_center <= counts[bond_axis_index] || throw(
        ArgumentError("terminal shellification atom core centers fall outside parent axis"),
    )
    atom_axis_coordinates =
        Tuple(location[bond_axis_index] for location in parent.atom_locations)
    atom_order = Tuple(sortperm(collect(atom_axis_coordinates)))
    length(atom_order) == 2 || throw(
        ArgumentError("terminal shellification diatomic path requires two atom locations"),
    )
    positions = Vector{NTuple{3,Int}}(undef, 2)
    for (atom_index, axis_center) in zip(atom_order, (low_center, high_center))
        positions[atom_index] =
            ntuple(axis -> axis == bond_axis_index ? axis_center : center[axis], 3)
    end
    return Tuple(positions)
end

function _pqs_source_box_route_driver_shell_stage_terminal_shellification(parent)
    parent_axes = _pqs_source_box_route_driver_terminal_parent_axes(parent)
    nuclear_positions =
        _pqs_source_box_route_driver_terminal_nuclear_index_positions(parent)
    policy =
        parent.system_classification == :one_center ?
        CartesianShellification.OneCenterShellification(;
            core_side = parent.standard_setup.core_cube_side,
            q = parent.standard_setup.q,
        ) :
        CartesianShellification.AtomOutwardShellification(;
            core_side = parent.standard_setup.core_cube_side,
            q = parent.standard_setup.q,
            bond_axis = parent.bond_axis,
        )
    plan = CartesianShellification.shellify(parent_axes, nuclear_positions, policy)
    scaffold = CartesianShellification.scaffold(
        plan;
        route_family = :white_lindsey_low_order,
    )
    return (;
        policy,
        plan,
        raw_plan = CartesianShellification.raw_plan(plan),
        scaffold,
        region_count = scaffold.region_count,
        ordered_region_roles = scaffold.ordered_region_roles,
        spatial_policy_order = scaffold.spatial_policy_order,
        coverage_complete = scaffold.coverage.coverage_complete,
        central_gap_region_count = scaffold.central_gap_region_count,
        central_midpoint_slab_count = scaffold.central_midpoint_slab_count,
        central_distorted_product_box_count =
            scaffold.central_distorted_product_box_count,
        central_distorted_product_box_metadata =
            scaffold.central_distorted_product_box_metadata,
    )
end

function _pqs_source_box_route_driver_shell_stage_low_order_shellization(
    parent,
    recipe;
)
    uses_terminal_shellification =
        (
            recipe.route_family == :white_lindsey_low_order &&
            parent.system_classification == :one_center
        ) ||
        (
            recipe.route_family == :pqs_source_box &&
            parent.system_classification == :one_center
        )
    uses_terminal_shellification ||
        return nothing

    terminal = _pqs_source_box_route_driver_shell_stage_terminal_shellification(parent)
    return (;
        route_family = recipe.route_family,
        shellification_kind = :terminal_cartesian_shellification_geometry,
        terminal_shellification = terminal,
        shellification_plan = terminal.plan,
        raw_shellification_plan = terminal.raw_plan,
        shellification_scaffold = terminal.scaffold,
        region_count = terminal.region_count,
        ordered_region_roles = terminal.ordered_region_roles,
        spatial_policy_order = terminal.spatial_policy_order,
        coverage_complete = terminal.coverage_complete,
        central_gap_region_count = terminal.central_gap_region_count,
        central_midpoint_slab_count = terminal.central_midpoint_slab_count,
        central_distorted_product_box_count =
            terminal.central_distorted_product_box_count,
        central_distorted_product_box_metadata =
            terminal.central_distorted_product_box_metadata,
    )
end

function cartesian_shells(
    parent,
    spacing_inputs,
    recipe;
)
    route_skeleton =
        _pqs_source_box_route_driver_route_skeleton(
            parent.axis_counts, spacing_inputs, recipe)
    shellification =
        _pqs_source_box_route_driver_shell_stage_low_order_shellization(
            parent,
            recipe;
        )

    return (;
        spacing_inputs,
        route_skeleton,
        shellification_plan =
            isnothing(shellification) ?
            nothing :
            shellification.shellification_plan,
        shellification_scaffold =
            isnothing(shellification) ?
            nothing :
            shellification.shellification_scaffold,
        shellification_kind =
            isnothing(shellification) ?
            nothing :
            shellification.shellification_kind,
        route_shape = route_skeleton.route_shape,
        source_boxes = route_skeleton.source_boxes,
    )
end

function _pqs_source_box_route_driver_terminal_lowering_family(route_family)
    route_family == :white_lindsey_low_order &&
        return :white_lindsey_low_order
    route_family == :pqs_source_box && return :pqs
    return nothing
end

function _pqs_source_box_route_driver_empty_terminal_lowering_contract_kind_counts()
    return _cartesian_terminal_region_lowering_contract_kind_counts(())
end

function _pqs_source_box_route_driver_terminal_lowering_policy(
    route_lowering_family,
    shellification_plan,
)
    route_lowering_family == :white_lindsey_low_order &&
        return CartesianTerminalLowering.WhiteLindseyLowering()
    if route_lowering_family == :pqs
        raw_plan = CartesianShellification.raw_plan(shellification_plan)
        q = !isnothing(raw_plan) && hasproperty(raw_plan, :q) ? raw_plan.q : nothing
        isnothing(q) &&
            throw(ArgumentError("PQS terminal lowering requires shellification q"))
        return CartesianTerminalLowering.PQSLowering(q = q)
    end
    return nothing
end

function _pqs_source_box_route_driver_terminal_lowering_plan(
    shellification_plan,
    route_lowering_family,
)
    shellification_plan isa CartesianShellification.ShellificationPlan || return nothing

    policy =
        _pqs_source_box_route_driver_terminal_lowering_policy(
            route_lowering_family,
            shellification_plan,
        )
    isnothing(policy) && return nothing
    return CartesianTerminalLowering.lower_terminal_regions(shellification_plan, policy)
end

function _pqs_source_box_route_driver_terminal_lowering_kind_counts(contracts)
    return (;
        direct_core_identity_cpb_count =
            count(
                contract ->
                    CartesianTerminalLowering.lowering_kind(contract) ==
                    :direct_core_identity_cpb,
                contracts,
            ),
        direct_slab_identity_cpb_count =
            count(
                contract ->
                    CartesianTerminalLowering.lowering_kind(contract) ==
                    :direct_slab_identity_cpb,
                contracts,
            ),
        direct_boundary_slab_identity_cpb_count =
            count(
                contract ->
                    CartesianTerminalLowering.lowering_kind(contract) ==
                    :direct_boundary_slab_identity_cpb,
                contracts,
            ),
        white_lindsey_boundary_strata_count =
            count(
                contract ->
                    CartesianTerminalLowering.lowering_kind(contract) ==
                    :white_lindsey_boundary_strata,
                contracts,
            ),
        pqs_filled_source_cpb_count =
            count(
                contract ->
                    CartesianTerminalLowering.lowering_kind(contract) ==
                    :pqs_filled_source_cpb,
                contracts,
            ),
        distorted_product_box_comx_count =
            count(
                contract ->
                    CartesianTerminalLowering.lowering_kind(contract) ==
                    :distorted_product_box_comx,
                contracts,
            ),
    )
end

function _pqs_source_box_route_driver_terminal_lowering_metadata(
    contract,
    field::Symbol,
    default = nothing,
)
    return hasproperty(contract.metadata, field) ?
           getproperty(contract.metadata, field) :
           default
end

function _pqs_source_box_route_driver_terminal_region_key(unit)
    return Symbol("terminal_region_", string(unit.terminal_region_order_index))
end

function _pqs_source_box_route_driver_unit_for_terminal_lowering_contract(
    unit_inventory,
    contract,
)
    matches = Tuple(
        unit for unit in unit_inventory.terminal_region_units
        if _pqs_source_box_route_driver_terminal_region_key(unit) ==
           contract.terminal_region_key
    )
    length(matches) == 1 ||
        throw(
            ArgumentError(
                "expected exactly one terminal-region unit for $(contract.terminal_region_key)",
            ),
        )
    return only(matches)
end

function _pqs_source_box_route_driver_terminal_lowering_source_cpb_plan_kind(
    contract,
)
    kind = CartesianTerminalLowering.lowering_kind(contract)
    kind in (
        :direct_core_identity_cpb,
        :direct_slab_identity_cpb,
        :direct_boundary_slab_identity_cpb,
    ) && return :direct_coordinate_product_source
    kind == :white_lindsey_boundary_strata &&
        return :complete_shell_boundary_strata
    kind == :pqs_filled_source_cpb && return :filled_source_cpb
    kind == :distorted_product_box_comx &&
        return :central_distorted_product_box_source
    return :terminal_lowering_source_cpb
end

function _pqs_source_box_route_driver_terminal_lowering_source_cpb_plan_box(
    contract,
)
    source_cpbs = CartesianTerminalLowering.source_cpbs(contract)
    length(source_cpbs) == 1 || return nothing
    return CartesianCPB.intervals(only(source_cpbs))
end

function _pqs_source_box_route_driver_terminal_lowering_source_family_counts(
    contract,
)
    return (;
        facet_cpb =
            _pqs_source_box_route_driver_terminal_lowering_metadata(
                contract,
                :facet_count,
                0,
            ),
        edge_cpb =
            _pqs_source_box_route_driver_terminal_lowering_metadata(
                contract,
                :edge_count,
                0,
            ),
        corner_cpb =
            _pqs_source_box_route_driver_terminal_lowering_metadata(
                contract,
                :corner_count,
                0,
            ),
    )
end

function _pqs_source_box_route_driver_terminal_lowering_retained_rule(contract)
    kind = CartesianTerminalLowering.lowering_kind(contract)
    kind == :white_lindsey_boundary_strata &&
        return :white_lindsey_boundary_strata_children
    kind == :pqs_filled_source_cpb &&
        return :boundary_comx_product_mode_selection
    return contract.retained_rule
end

function _pqs_source_box_route_driver_terminal_lowering_shell_realization_status(
    contract,
)
    kind = CartesianTerminalLowering.lowering_kind(contract)
    kind == :pqs_filled_source_cpb &&
        return :projection_lowdin_planned_not_materialized
    return :not_materialized
end

function _pqs_source_box_route_driver_terminal_lowering_final_granularity(
    contract,
)
    kind = CartesianTerminalLowering.lowering_kind(contract)
    kind == :white_lindsey_boundary_strata &&
        return :white_lindsey_boundary_strata_children_planned
    kind == :pqs_filled_source_cpb &&
        return :single_pqs_shell_realized_unit_planned
    kind == :distorted_product_box_comx &&
        return :single_distorted_product_box_unit_planned
    return :single_direct_source_cpb
end

function _pqs_source_box_route_driver_terminal_lowering_contract_record(
    contract,
    unit,
    contract_index::Int,
)
    lowering_kind = CartesianTerminalLowering.lowering_kind(contract)
    source_cpbs = CartesianTerminalLowering.source_cpbs(contract)
    source_cpb_plan_box =
        _pqs_source_box_route_driver_terminal_lowering_source_cpb_plan_box(
            contract,
        )
    source_cpb_count = length(source_cpbs)
    identity_like =
        _pqs_source_box_route_driver_terminal_lowering_metadata(
            contract,
            :identity_like,
            false,
        )
    source_cpb_plan_equals_owned_support =
        unit.terminal_region_kind != :complete_shell
    source_mode_shape =
        _pqs_source_box_route_driver_terminal_lowering_metadata(
            contract,
            :source_mode_shape,
            nothing,
        )
    return (;
        object_kind = :cartesian_terminal_region_lowering_contract,
        status = :planned_not_materialized,
        contract_index,
        contract_key =
            Symbol(String(unit.unit_key), "_", String(lowering_kind)),
        terminal_region_key = contract.terminal_region_key,
        unit_index = unit.unit_index,
        unit_key = unit.unit_key,
        unit_role = unit.unit_role,
        unit_kind = unit.unit_kind,
        terminal_region_order_index = unit.terminal_region_order_index,
        terminal_region_role = unit.terminal_region_role,
        terminal_region_kind = unit.terminal_region_kind,
        lowering_contract_kind = lowering_kind,
        outer_box = unit.outer_box,
        box = unit.outer_box,
        inner_exclusion_box = unit.inner_exclusion_box,
        support_count = unit.support_count,
        owned_support_status = unit.owned_support_status,
        owned_support_is_cpb = unit.owned_support_is_cpb,
        shellification_region_is_cpb = unit.shellification_region_is_cpb,
        shellification_region_is_lowering_source =
            unit.shellification_region_is_lowering_source,
        source_cpb = source_cpb_count == 1 ? only(source_cpbs) : nothing,
        source_cpbs,
        source_cpb_plan_status = :planned_not_materialized,
        source_cpb_plan_box,
        source_cpb_plan_kind =
            _pqs_source_box_route_driver_terminal_lowering_source_cpb_plan_kind(
                contract,
            ),
        source_cpb_plan_equals_owned_support,
        source_cpb_count,
        source_cpb_family_counts =
            _pqs_source_box_route_driver_terminal_lowering_source_family_counts(
                contract,
            ),
        source_cpbs_materialized = false,
        identity_like_source_contract = identity_like,
        retained_rule =
            _pqs_source_box_route_driver_terminal_lowering_retained_rule(
                contract,
            ),
        source_mode_shape,
        q =
            _pqs_source_box_route_driver_terminal_lowering_metadata(
                contract,
                :q,
                nothing,
            ),
        L =
            _pqs_source_box_route_driver_terminal_lowering_metadata(
                contract,
                :L,
                nothing,
            ),
        aspect_ratio =
            _pqs_source_box_route_driver_terminal_lowering_metadata(
                contract,
                :aspect_ratio,
                nothing,
            ),
        intermediate_retained_space_status =
            lowering_kind in (:pqs_filled_source_cpb, :distorted_product_box_comx) ?
            :planned_not_materialized :
            :not_materialized,
        shell_realization_status =
            _pqs_source_box_route_driver_terminal_lowering_shell_realization_status(
                contract,
            ),
        face_edge_corner_decomposition_required =
            lowering_kind == :white_lindsey_boundary_strata,
        coefficient_maps_materialized = false,
        transform_contracts_materialized = false,
        retained_spaces_materialized = false,
        operator_blocks_materialized = false,
        pair_operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
        final_retained_unit_status = :not_materialized,
        final_retained_unit_records_materialized = false,
        final_unit_count_planned =
            lowering_kind == :white_lindsey_boundary_strata ? source_cpb_count : 1,
        final_unit_granularity =
            _pqs_source_box_route_driver_terminal_lowering_final_granularity(
                contract,
            ),
        metadata = merge(
            unit.metadata,
            (;
                terminal_lowering_contract_key = contract.contract_key,
                terminal_lowering_plan_contract = true,
            ),
        ),
        provenance = (;
            source = :terminal_lowering_plan_compatibility_adapter,
            unit_inventory_source =
                :_cartesian_terminal_shellification_region_unit_inventory,
            unit_key = unit.unit_key,
            terminal_region_key = contract.terminal_region_key,
            terminal_region_role = unit.terminal_region_role,
            terminal_region_kind = unit.terminal_region_kind,
        ),
    )
end

function _pqs_source_box_route_driver_terminal_lowering_contract_inventory_from_plan(
    terminal_lowering_plan,
    unit_inventory,
)
    terminal_lowering_plan isa CartesianTerminalLowering.TerminalLoweringPlan ||
        throw(
            ArgumentError(
                "terminal lowering compatibility inventory requires a TerminalLoweringPlan",
            ),
        )
    unit_inventory.object_kind == :cartesian_terminal_region_unit_inventory ||
        throw(
            ArgumentError(
                "terminal lowering compatibility inventory requires a terminal-region unit inventory",
            ),
        )

    typed_contracts =
        CartesianTerminalLowering.available_contracts(terminal_lowering_plan)
    lowering_contracts = Tuple(
        _pqs_source_box_route_driver_terminal_lowering_contract_record(
            contract,
            _pqs_source_box_route_driver_unit_for_terminal_lowering_contract(
                unit_inventory,
                contract,
            ),
            contract_index,
        )
        for (contract_index, contract) in enumerate(typed_contracts)
    )
    contract_counts_by_unit = Tuple(
        (;
            unit_key = unit.unit_key,
            lowering_contract_count =
                count(
                    contract -> contract.unit_key == unit.unit_key,
                    lowering_contracts,
                ),
        ) for unit in unit_inventory.terminal_region_units
    )
    complete_shell_contracts = Tuple(
        contract for contract in lowering_contracts
        if contract.terminal_region_kind == :complete_shell
    )
    lw_contracts = Tuple(
        contract for contract in lowering_contracts
        if contract.lowering_contract_kind == :white_lindsey_boundary_strata
    )

    return (;
        object_kind = :cartesian_terminal_region_lowering_contract_inventory,
        status = :available_terminal_region_lowering_contract_inventory,
        inventory_source = :terminal_lowering_plan_compatibility_adapter,
        source_object_kind = :cartesian_terminal_lowering_plan,
        private_development_only = true,
        terminal_region_unit_count = unit_inventory.unit_count,
        lowering_contract_count = length(lowering_contracts),
        lowering_contracts,
        contract_records = lowering_contracts,
        unit_keys = unit_inventory.unit_keys,
        unit_roles = unit_inventory.unit_roles,
        unit_kinds = unit_inventory.unit_kinds,
        terminal_region_roles = unit_inventory.terminal_region_roles,
        terminal_region_kinds = unit_inventory.terminal_region_kinds,
        support_counts = unit_inventory.support_counts,
        contract_counts_by_unit,
        lowering_contract_kinds =
            Tuple(contract.lowering_contract_kind for contract in lowering_contracts),
        lowering_contract_kind_counts =
            _cartesian_terminal_region_lowering_contract_kind_counts(
                lowering_contracts,
            ),
        complete_shell_unit_count = unit_inventory.complete_shell_unit_count,
        complete_shell_lowering_contract_count = length(complete_shell_contracts),
        lw_complete_shell_lowering_contract_count = length(lw_contracts),
        lw_complete_shell_cpb_count =
            sum(contract -> contract.source_cpb_count, lw_contracts; init = 0),
        lw_complete_shell_cpb_family_counts =
            (
                facet_cpb =
                    sum(
                        contract -> contract.source_cpb_family_counts.facet_cpb,
                        lw_contracts;
                        init = 0,
                    ),
                edge_cpb =
                    sum(
                        contract -> contract.source_cpb_family_counts.edge_cpb,
                        lw_contracts;
                        init = 0,
                    ),
                corner_cpb =
                    sum(
                        contract -> contract.source_cpb_family_counts.corner_cpb,
                        lw_contracts;
                        init = 0,
                    ),
            ),
        final_retained_unit_inventory_available = false,
        pair_inventory_available = false,
        pair_inventory_status = :not_available_lowering_contract_metadata_only,
        coefficient_maps_materialized = false,
        transform_contracts_materialized = false,
        retained_spaces_materialized = false,
        operator_blocks_materialized = false,
        pair_operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
        diagnostics = (;
            source = :terminal_lowering_plan_compatibility_adapter,
            private_development_only = true,
            terminal_region_metadata_only = true,
            lowering_contracts_metadata_only = true,
            all_units_have_lowering_contracts =
                all(entry -> entry.lowering_contract_count >= 1, contract_counts_by_unit),
            final_retained_unit_inventory_available = false,
            pair_inventory_available = false,
            coefficient_maps_materialized = false,
            transform_contracts_materialized = false,
            retained_spaces_materialized = false,
            operator_blocks_materialized = false,
            pair_operator_blocks_materialized = false,
            hamiltonian_data_materialized = false,
            artifacts_materialized = false,
            materialization_behavior_changed = false,
            public_default_behavior_changed = false,
        ),
    )
end

function _pqs_source_box_route_driver_compat_contract_for_terminal_contract(
    lowering_contract_inventory,
    terminal_contract,
)
    matches = Tuple(
        contract for contract in lowering_contract_inventory.lowering_contracts
        if contract.terminal_region_key == terminal_contract.terminal_region_key &&
           contract.lowering_contract_kind ==
           CartesianTerminalLowering.lowering_kind(terminal_contract)
    )
    length(matches) == 1 ||
        throw(
            ArgumentError(
                "expected exactly one compatibility contract for $(terminal_contract.contract_key)",
            ),
        )
    return only(matches)
end

function _pqs_source_box_route_driver_selected_terminal_lowering_contract_inventory_from_plan(
    terminal_lowering_plan,
    lowering_contract_inventory,
    route_lowering_family::Symbol,
)
    terminal_lowering_plan isa CartesianTerminalLowering.TerminalLoweringPlan ||
        throw(
            ArgumentError(
                "selected terminal lowering compatibility inventory requires a TerminalLoweringPlan",
            ),
        )
    selected_contracts = Tuple(
        _pqs_source_box_route_driver_compat_contract_for_terminal_contract(
            lowering_contract_inventory,
            contract,
        )
        for contract in CartesianTerminalLowering.selected_contracts(
            terminal_lowering_plan,
        )
    )
    selected_contract_keys =
        Tuple(contract.contract_key for contract in selected_contracts)
    unselected_contracts = Tuple(
        contract for contract in lowering_contract_inventory.lowering_contracts
        if !(contract.contract_key in selected_contract_keys)
    )
    unit_keys = lowering_contract_inventory.unit_keys
    selected_contract_counts_by_unit = Tuple(
        (;
            unit_key,
            selected_contract_count =
                count(
                    contract -> contract.unit_key == unit_key,
                    selected_contracts,
                ),
        ) for unit_key in unit_keys
    )
    all_units_have_exactly_one_selected_contract =
        all(entry -> entry.selected_contract_count == 1, selected_contract_counts_by_unit)

    return (;
        object_kind = :cartesian_selected_terminal_lowering_contract_inventory,
        status = :available_selected_terminal_lowering_contract_inventory,
        inventory_source = :terminal_lowering_plan_compatibility_adapter,
        source_object_kind = lowering_contract_inventory.object_kind,
        private_development_only = true,
        route_lowering_family,
        terminal_region_unit_count = length(unit_keys),
        selected_contract_count = length(selected_contracts),
        selected_contracts,
        selected_contract_records = selected_contracts,
        selected_contract_kinds =
            Tuple(contract.lowering_contract_kind for contract in selected_contracts),
        selected_contract_kind_counts =
            _cartesian_terminal_region_lowering_contract_kind_counts(
                selected_contracts,
            ),
        selected_contract_counts_by_unit,
        all_units_have_exactly_one_selected_contract,
        unselected_contract_count = length(unselected_contracts),
        unselected_contracts,
        unselected_contract_kinds =
            Tuple(contract.lowering_contract_kind for contract in unselected_contracts),
        final_retained_unit_inventory_available = false,
        pair_inventory_available = false,
        pair_inventory_status = :not_available_selected_lowering_metadata_only,
        coefficient_maps_materialized = false,
        transform_contracts_materialized = false,
        retained_spaces_materialized = false,
        operator_blocks_materialized = false,
        pair_operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
        diagnostics = (;
            source = :terminal_lowering_plan_compatibility_adapter,
            private_development_only = true,
            selected_lowering_metadata_only = true,
            all_units_have_exactly_one_selected_contract,
            final_retained_unit_inventory_available = false,
            pair_inventory_available = false,
            coefficient_maps_materialized = false,
            transform_contracts_materialized = false,
            retained_spaces_materialized = false,
            operator_blocks_materialized = false,
            pair_operator_blocks_materialized = false,
            hamiltonian_data_materialized = false,
            artifacts_materialized = false,
            materialization_behavior_changed = false,
            public_default_behavior_changed = false,
        ),
    )
end

function _pqs_source_box_route_driver_selected_terminal_lowering_fields(
    selected_terminal_lowering_contract_inventory,
    status::Symbol,
    route_lowering_family,
    terminal_lowering_plan = nothing,
)
    selected_inventory_available =
        !isnothing(selected_terminal_lowering_contract_inventory) &&
        selected_terminal_lowering_contract_inventory.status ==
        :available_selected_terminal_lowering_contract_inventory
    terminal_lowering_plan_available =
        terminal_lowering_plan isa CartesianTerminalLowering.TerminalLoweringPlan
    terminal_lowering_summary =
        terminal_lowering_plan_available ?
        CartesianTerminalLowering.summary(terminal_lowering_plan) :
        nothing
    selected_contracts =
        terminal_lowering_plan_available ?
        CartesianTerminalLowering.selected_contracts(terminal_lowering_plan) :
        ()
    available_contracts =
        terminal_lowering_plan_available ?
        CartesianTerminalLowering.available_contracts(terminal_lowering_plan) :
        ()
    selected_contract_keys =
        Tuple(contract.contract_key for contract in selected_contracts)
    unselected_contracts =
        terminal_lowering_plan_available ?
        Tuple(
            contract for contract in available_contracts
            if !(contract.contract_key in selected_contract_keys)
        ) :
        ()
    return (;
        terminal_shellification_selected_lowering_contract_inventory_available =
            selected_inventory_available,
        terminal_shellification_selected_lowering_contract_inventory_status =
            status,
        terminal_shellification_selected_lowering_contract_inventory =
            selected_terminal_lowering_contract_inventory,
        terminal_shellification_selected_lowering_family =
            route_lowering_family,
        terminal_shellification_selected_contract_count =
            terminal_lowering_plan_available ?
            terminal_lowering_summary.selected_contract_count :
            selected_inventory_available ?
            selected_terminal_lowering_contract_inventory.selected_contract_count :
            0,
        terminal_shellification_selected_contract_kinds =
            terminal_lowering_plan_available ?
            terminal_lowering_summary.selected_contract_kinds :
            selected_inventory_available ?
            selected_terminal_lowering_contract_inventory.selected_contract_kinds :
            (),
        terminal_shellification_selected_contract_kind_counts =
            terminal_lowering_plan_available ?
            _pqs_source_box_route_driver_terminal_lowering_kind_counts(
                selected_contracts,
            ) :
            selected_inventory_available ?
            selected_terminal_lowering_contract_inventory.selected_contract_kind_counts :
            _pqs_source_box_route_driver_empty_terminal_lowering_contract_kind_counts(),
        terminal_shellification_selected_contract_counts_by_unit =
            selected_inventory_available ?
            selected_terminal_lowering_contract_inventory.selected_contract_counts_by_unit :
            (),
        terminal_shellification_all_units_have_exactly_one_selected_contract =
            terminal_lowering_plan_available ?
            terminal_lowering_summary.all_terminal_regions_have_selected_contract :
            selected_inventory_available &&
            selected_terminal_lowering_contract_inventory.all_units_have_exactly_one_selected_contract,
        terminal_shellification_unselected_contract_count =
            terminal_lowering_plan_available ?
            length(unselected_contracts) :
            selected_inventory_available ?
            selected_terminal_lowering_contract_inventory.unselected_contract_count :
            0,
        terminal_shellification_unselected_contract_kinds =
            terminal_lowering_plan_available ?
            Tuple(contract.lowering_kind for contract in unselected_contracts) :
            selected_inventory_available ?
            selected_terminal_lowering_contract_inventory.unselected_contract_kinds :
            (),
    )
end

function _pqs_source_box_route_driver_terminal_route_state_summary(;
    status,
    selected::Bool,
    route_lowering_family,
    shellification_plan_available::Bool,
    unit_inventory_available::Bool,
    unit_inventory,
    lowering_plan_available::Bool,
    lowering_summary,
    selected_contract_inventory_available::Bool,
    selected_contract_inventory,
    selected_crc_sidecar_summary,
    retained_unit_summary,
    retained_unit_transform_contract_summary,
    unit_pair_summary,
    pair_operator_summary,
    pair_block_materialization_summary,
)
    return (;
        object_kind = :cartesian_driver_terminal_route_state_summary,
        status,
        selected,
        route_lowering_family,
        shellification_plan_available,
        unit_inventory_available,
        unit_count =
            unit_inventory_available && !isnothing(unit_inventory) ?
            unit_inventory.unit_count :
            0,
        unit_keys =
            unit_inventory_available && !isnothing(unit_inventory) ?
            unit_inventory.unit_keys :
            (),
        unit_kinds =
            unit_inventory_available && !isnothing(unit_inventory) ?
            unit_inventory.unit_kinds :
            (),
        lowering_plan_available,
        lowering_plan_status =
            lowering_plan_available && !isnothing(lowering_summary) ?
            lowering_summary.status :
            :not_available,
        available_contract_count =
            lowering_plan_available && !isnothing(lowering_summary) ?
            lowering_summary.available_contract_count :
            0,
        available_contract_kinds =
            lowering_plan_available && !isnothing(lowering_summary) ?
            lowering_summary.available_contract_kinds :
            (),
        selected_contract_inventory_available,
        selected_contract_count =
            lowering_plan_available && !isnothing(lowering_summary) ?
            lowering_summary.selected_contract_count :
            selected_contract_inventory_available &&
            !isnothing(selected_contract_inventory) ?
            selected_contract_inventory.selected_contract_count :
            0,
        selected_contract_kinds =
            lowering_plan_available && !isnothing(lowering_summary) ?
            lowering_summary.selected_contract_kinds :
            selected_contract_inventory_available &&
            !isnothing(selected_contract_inventory) ?
            selected_contract_inventory.selected_contract_kinds :
            (),
        all_units_have_selected_contract =
            lowering_plan_available && !isnothing(lowering_summary) ?
            lowering_summary.all_terminal_regions_have_selected_contract :
            selected_contract_inventory_available &&
            !isnothing(selected_contract_inventory) &&
            selected_contract_inventory.all_units_have_exactly_one_selected_contract,
        selected_crc_sidecar_status =
            !isnothing(selected_crc_sidecar_summary) ?
            selected_crc_sidecar_summary.status :
            :not_available,
        selected_crc_sidecar_complete =
            !isnothing(selected_crc_sidecar_summary) ?
            selected_crc_sidecar_summary.sidecar_inventory_complete :
            false,
        selected_crc_sidecar_available_count =
            !isnothing(selected_crc_sidecar_summary) ?
            selected_crc_sidecar_summary.sidecar_available_count :
            0,
        selected_crc_sidecar_missing_count =
            !isnothing(selected_crc_sidecar_summary) ?
            selected_crc_sidecar_summary.sidecar_missing_count :
            0,
        retained_unit_summary,
        retained_unit_transform_contract_summary,
        unit_pair_summary,
        pair_operator_summary,
        pair_block_materialization_summary,
        final_retained_unit_inventory_available =
            !isnothing(selected_crc_sidecar_summary) &&
            selected_crc_sidecar_summary.final_retained_unit_inventory_available,
        pair_inventory_available =
            !isnothing(unit_pair_summary) &&
            unit_pair_summary.status == :available_unit_pair_plan,
        operator_blocks_materialized = false,
        pair_operator_blocks_materialized =
            !isnothing(selected_crc_sidecar_summary) &&
            selected_crc_sidecar_summary.pair_operator_blocks_materialized,
        hamiltonian_data_materialized =
            !isnothing(selected_crc_sidecar_summary) &&
            selected_crc_sidecar_summary.hamiltonian_data_materialized,
        artifacts_materialized =
            !isnothing(selected_crc_sidecar_summary) &&
            selected_crc_sidecar_summary.artifacts_materialized,
    )
end

function _pqs_source_box_route_driver_terminal_route_state(;
    status,
    selected::Bool,
    route_lowering_family = nothing,
    shellification_plan = nothing,
    unit_inventory = nothing,
    lowering_plan = nothing,
    lowering_summary = nothing,
    lowering_contract_inventory = nothing,
    selected_contract_inventory = nothing,
    selected_crc_sidecar_summary = nothing,
    retained_unit_plan = nothing,
    retained_unit_summary = nothing,
    retained_unit_transform_contract_plan = nothing,
    retained_unit_transform_contract_summary = nothing,
    unit_pair_plan = nothing,
    unit_pair_summary = nothing,
    pair_operator_plan = nothing,
    pair_operator_summary = nothing,
    pair_block_materialization_plan = nothing,
    pair_block_materialization_summary = nothing,
    pair_operator_route_core_sidecars::Bool = true,
    blocker = nothing,
)
    shellification_plan_available =
        shellification_plan isa CartesianShellification.ShellificationPlan
    unit_inventory_available =
        !isnothing(unit_inventory) &&
        hasproperty(unit_inventory, :status) &&
        unit_inventory.status == :available_terminal_region_unit_inventory
    lowering_plan_available =
        lowering_plan isa CartesianTerminalLowering.TerminalLoweringPlan
    lowering_contract_inventory_available =
        !isnothing(lowering_contract_inventory) &&
        hasproperty(lowering_contract_inventory, :status) &&
        lowering_contract_inventory.status ==
        :available_terminal_region_lowering_contract_inventory
    selected_contract_inventory_available =
        !isnothing(selected_contract_inventory) &&
        hasproperty(selected_contract_inventory, :status) &&
        selected_contract_inventory.status ==
        :available_selected_terminal_lowering_contract_inventory
    selected_crc_sidecar_summary = nothing
    retained_unit_plan =
        isnothing(retained_unit_plan) && selected && lowering_plan_available ?
        CartesianRetainedUnits.retained_unit_plan(lowering_plan) :
        retained_unit_plan
    retained_unit_summary =
        isnothing(retained_unit_summary) && retained_unit_plan isa
        CartesianRetainedUnits.RetainedUnitPlan ?
        CartesianRetainedUnits.summary(retained_unit_plan) :
        isnothing(retained_unit_summary) ?
        CartesianRetainedUnits.unavailable_summary(
            selected ? status : :not_selected,
            blocker,
        ) :
        retained_unit_summary
    retained_unit_transform_contract_plan =
        isnothing(retained_unit_transform_contract_plan) &&
        retained_unit_plan isa CartesianRetainedUnits.RetainedUnitPlan ?
        CartesianRetainedUnitTransformContracts.retained_unit_transform_contract_plan(
            retained_unit_plan,
        ) :
        retained_unit_transform_contract_plan
    retained_unit_transform_contract_summary =
        isnothing(retained_unit_transform_contract_summary) &&
        retained_unit_transform_contract_plan isa
        CartesianRetainedUnitTransformContracts.RetainedUnitTransformContractPlan ?
        CartesianRetainedUnitTransformContracts.summary(
            retained_unit_transform_contract_plan,
        ) :
        isnothing(retained_unit_transform_contract_summary) ?
        CartesianRetainedUnitTransformContracts.unavailable_summary(
            selected ? status : :not_selected,
            blocker,
        ) :
        retained_unit_transform_contract_summary
    unit_pair_plan =
        isnothing(unit_pair_plan) && retained_unit_plan isa
        CartesianRetainedUnits.RetainedUnitPlan ?
        CartesianUnitPairs.unit_pair_plan(retained_unit_plan) :
        unit_pair_plan
    unit_pair_summary =
        isnothing(unit_pair_summary) && unit_pair_plan isa
        CartesianUnitPairs.UnitPairPlan ?
        CartesianUnitPairs.summary(unit_pair_plan) :
        isnothing(unit_pair_summary) ?
        CartesianUnitPairs.unavailable_summary(
            selected ? status : :not_selected,
            blocker,
        ) :
        unit_pair_summary
    pair_operator_plan =
        isnothing(pair_operator_plan) &&
        unit_pair_plan isa CartesianUnitPairs.UnitPairPlan &&
        retained_unit_transform_contract_plan isa
        CartesianRetainedUnitTransformContracts.RetainedUnitTransformContractPlan ?
        CartesianPairOperatorPlans.pair_operator_plan(
            unit_pair_plan,
            retained_unit_transform_contract_plan,
            route_core_sidecars = pair_operator_route_core_sidecars,
        ) :
        isnothing(pair_operator_plan) && unit_pair_plan isa
        CartesianUnitPairs.UnitPairPlan ?
        CartesianPairOperatorPlans.pair_operator_plan(
            unit_pair_plan;
            route_core_sidecars = pair_operator_route_core_sidecars,
        ) :
        pair_operator_plan
    pair_operator_summary =
        isnothing(pair_operator_summary) && pair_operator_plan isa
        CartesianPairOperatorPlans.PairOperatorPlan ?
        CartesianPairOperatorPlans.summary(pair_operator_plan) :
        isnothing(pair_operator_summary) ?
        CartesianPairOperatorPlans.unavailable_summary(
            selected ? status : :not_selected,
            blocker,
        ) :
        pair_operator_summary
    pair_block_materialization_plan =
        isnothing(pair_block_materialization_plan) &&
        pair_operator_plan isa CartesianPairOperatorPlans.PairOperatorPlan ?
        CartesianPairBlockMaterialization.pair_block_materialization_plan(
            pair_operator_plan,
        ) :
        pair_block_materialization_plan
    pair_block_materialization_summary =
        isnothing(pair_block_materialization_summary) &&
        pair_block_materialization_plan isa
        CartesianPairBlockMaterialization.PairBlockMaterializationPlan ?
        CartesianPairBlockMaterialization.summary(pair_block_materialization_plan) :
        isnothing(pair_block_materialization_summary) ?
        CartesianPairBlockMaterialization.unavailable_summary(
            selected ? status : :not_selected,
            blocker,
        ) :
        pair_block_materialization_summary
    summary =
        _pqs_source_box_route_driver_terminal_route_state_summary(;
            status,
            selected,
            route_lowering_family,
            shellification_plan_available,
            unit_inventory_available,
            unit_inventory,
            lowering_plan_available,
            lowering_summary,
            selected_contract_inventory_available,
            selected_contract_inventory,
            selected_crc_sidecar_summary,
            retained_unit_summary,
            retained_unit_transform_contract_summary,
            unit_pair_summary,
            pair_operator_summary,
            pair_block_materialization_summary,
        )

    return (;
        object_kind = :cartesian_driver_terminal_route_state,
        status,
        selected,
        available = summary.shellification_plan_available,
        private_development_only = true,
        route_lowering_family,
        shellification_plan_available,
        shellification_plan,
        unit_inventory_available,
        unit_inventory,
        lowering_plan_available,
        lowering_plan,
        lowering_summary,
        lowering_contract_inventory_available,
        lowering_contract_inventory,
        selected_contract_inventory_available,
        selected_contract_inventory,
        selected_crc_sidecar_summary,
        retained_unit_plan,
        retained_unit_summary,
        retained_unit_transform_contract_plan,
        retained_unit_transform_contract_summary,
        unit_pair_plan,
        unit_pair_summary,
        pair_operator_plan,
        pair_operator_summary,
        pair_block_materialization_plan,
        pair_block_materialization_summary,
        summary,
        blocker,
        compatibility_alias_status =
            :legacy_scalar_aliases_derived_from_terminal_route_state,
    )
end

function _pqs_source_box_route_driver_terminal_route_state_unavailable(
    status::Symbol,
    blocker = nothing,
)
    return _pqs_source_box_route_driver_terminal_route_state(;
        status,
        selected = false,
        blocker,
    )
end

function _pqs_source_box_route_driver_unit_stage_low_order_summary(shells)
    shellification_plan = get(shells, :shellification_plan, nothing)
    shellification_scaffold = get(shells, :shellification_scaffold, nothing)
    shellification_kind = get(shells, :shellification_kind, nothing)
    (isnothing(shellification_plan) || isnothing(shellification_scaffold)) &&
        return nothing

    unit_inventory =
        _cartesian_terminal_shellification_region_unit_inventory(
            shellification_scaffold,
        )
    route_lowering_family =
        _pqs_source_box_route_driver_terminal_lowering_family(
            shells.route_skeleton.route_family)
    lowering_plan =
        _pqs_source_box_route_driver_terminal_lowering_plan(
            shellification_plan,
            route_lowering_family,
        )
    lowering_contract_inventory =
        lowering_plan isa CartesianTerminalLowering.TerminalLoweringPlan ?
        _pqs_source_box_route_driver_terminal_lowering_contract_inventory_from_plan(
            lowering_plan,
            unit_inventory,
        ) :
        nothing

    return (;
        shellification_kind,
        shellification_plan,
        shellification_scaffold,
        unit_inventory,
        route_lowering_family,
        lowering_plan,
        lowering_contract_inventory,
        region_count = shellification_scaffold.region_count,
        ordered_region_roles = shellification_scaffold.ordered_region_roles,
        central_gap_region_count =
            shellification_scaffold.central_gap_region_count,
        central_midpoint_slab_count =
            shellification_scaffold.central_midpoint_slab_count,
        central_distorted_product_box_count =
            shellification_scaffold.central_distorted_product_box_count,
        central_distorted_product_box_metadata =
            shellification_scaffold.central_distorted_product_box_metadata,
    )
end

function cartesian_units(parent, shells, recipe)
    low_order_units =
        _pqs_source_box_route_driver_unit_stage_low_order_summary(shells)

    return (;
        parent,
        route_family = recipe.route_family,
        route_kind = recipe.route_kind,
        route_skeleton = shells.route_skeleton,
        low_order_units,
        shellification_plan = shells.shellification_plan,
        shellification_scaffold = shells.shellification_scaffold,
        shellification_kind = shells.shellification_kind,
        route_shape = shells.route_shape,
        source_boxes = shells.source_boxes,
        source_dimensions = shells.route_skeleton.source_dimensions,
        retained_units = shells.route_skeleton.retained_units,
        retained_unit_order = shells.route_skeleton.retained_unit_order,
        retained_counts = shells.route_skeleton.retained_counts,
        ranges = shells.route_skeleton.ranges,
        retained_dimension = shells.route_skeleton.retained_dimension,
        pair_entries = shells.route_skeleton.pair_entries,
        pair_family_counts = shells.route_skeleton.pair_family_counts,
        helper_by_pair_family = shells.route_skeleton.helper_by_pair_family,
    )
end

function _pqs_source_box_route_driver_transform_stage_low_order_summary(units)
    low_order_units = get(units, :low_order_units, nothing)
    isnothing(low_order_units) && return nothing

    retained_units = get(low_order_units, :retained_units, ())
    return (;
        shellification_kind = low_order_units.shellification_kind,
        shellification_plan = low_order_units.shellification_plan,
        shellification_scaffold = low_order_units.shellification_scaffold,
        unit_inventory = low_order_units.unit_inventory,
        route_lowering_family = low_order_units.route_lowering_family,
        lowering_plan = low_order_units.lowering_plan,
        lowering_contract_inventory = low_order_units.lowering_contract_inventory,
        retained_units,
        retained_counts =
            _pqs_source_box_route_driver_named_tuple_from_units(
                retained_units, :retained_count),
        ranges =
            _pqs_source_box_route_driver_named_tuple_from_units(
                retained_units, :retained_range),
        retained_dimension = get(low_order_units, :retained_dimension, nothing),
    )
end

function cartesian_transforms(units, recipe)
    retained_units = units.retained_units
    low_order_transforms =
        _pqs_source_box_route_driver_transform_stage_low_order_summary(units)
    return (;
        route_family = recipe.route_family,
        route_kind = recipe.route_kind,
        retained_units,
        low_order_transforms,
        shellification_plan = units.shellification_plan,
        shellification_scaffold = units.shellification_scaffold,
        shellification_kind = units.shellification_kind,
        route_skeleton = units.route_skeleton,
        pair_entries = units.pair_entries,
        pair_family_counts = units.pair_family_counts,
        helper_by_pair_family = units.helper_by_pair_family,
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

function _pqs_source_box_route_driver_pair_keys_from_entries(pair_entries)
    return Tuple(pair.pair_key for pair in pair_entries)
end

function _pqs_source_box_route_driver_pair_stage_low_order_summary(
    transforms,
)
    low_order_transforms = get(transforms, :low_order_transforms, nothing)
    isnothing(low_order_transforms) && return nothing

    return (;
        shellification_kind = low_order_transforms.shellification_kind,
        shellification_plan = low_order_transforms.shellification_plan,
        shellification_scaffold = low_order_transforms.shellification_scaffold,
        unit_inventory = low_order_transforms.unit_inventory,
        route_lowering_family = low_order_transforms.route_lowering_family,
        lowering_plan = low_order_transforms.lowering_plan,
        lowering_contract_inventory =
            low_order_transforms.lowering_contract_inventory,
        pair_entries = transforms.pair_entries,
        pair_keys = _pqs_source_box_route_driver_pair_keys_from_entries(
            transforms.pair_entries),
        pair_family_counts = transforms.pair_family_counts,
        helper_by_pair_family = transforms.helper_by_pair_family,
    )
end

function cartesian_pair_terms(units, transforms, recipe)
    low_order_pairs =
        _pqs_source_box_route_driver_pair_stage_low_order_summary(
            transforms,
        )

    return (;
        route_family = recipe.route_family,
        route_kind = recipe.route_kind,
        retained_dimension = transforms.retained_dimension,
        low_order_pairs,
        pair_entries = transforms.pair_entries,
        pair_keys = _pqs_source_box_route_driver_pair_keys_from_entries(
            transforms.pair_entries),
        pair_family_counts = transforms.pair_family_counts,
        helper_by_pair_family = transforms.helper_by_pair_family,
    )
end

function _pqs_source_box_route_driver_assembly_stage_low_order_summary(pairs)
    low_order_pairs =
        hasproperty(pairs, :low_order_pairs) ? pairs.low_order_pairs : nothing
    isnothing(low_order_pairs) && return nothing

    pair_entries = hasproperty(low_order_pairs, :pair_entries) ?
                   low_order_pairs.pair_entries :
                   pairs.pair_entries
    return (;
        shellification_kind = low_order_pairs.shellification_kind,
        shellification_plan = low_order_pairs.shellification_plan,
        shellification_scaffold = low_order_pairs.shellification_scaffold,
        unit_inventory = low_order_pairs.unit_inventory,
        route_lowering_family = low_order_pairs.route_lowering_family,
        lowering_plan = low_order_pairs.lowering_plan,
        lowering_contract_inventory =
            low_order_pairs.lowering_contract_inventory,
        pair_entries,
        pair_keys = _pqs_source_box_route_driver_pair_keys_from_entries(
            pair_entries),
        pair_family_counts = hasproperty(low_order_pairs, :pair_family_counts) ?
                             low_order_pairs.pair_family_counts :
                             pairs.pair_family_counts,
        helper_by_pair_family =
            hasproperty(low_order_pairs, :helper_by_pair_family) ?
            low_order_pairs.helper_by_pair_family :
            pairs.helper_by_pair_family,
    )
end

function _pqs_source_box_route_driver_complete_core_shell_h1_j_missing_inputs(;
    region_plan,
    source_plan,
    final_basis,
    h1_payload,
    axis_weights,
    raw_pair_factor_terms,
    coulomb_expansion,
)
    missing = Symbol[]
    isnothing(region_plan) && push!(missing, :pqs_multilayer_shell_region_plan)
    isnothing(source_plan) && push!(missing, :pqs_multilayer_shell_source_plan)
    isnothing(final_basis) &&
        push!(missing, :pqs_multilayer_complete_core_shell_final_basis)
    isnothing(h1_payload) &&
        push!(missing, :pqs_multilayer_complete_core_shell_h1_payload)
    isnothing(axis_weights) && push!(missing, :axis_weights)
    isnothing(raw_pair_factor_terms) && push!(missing, :raw_pair_factor_terms)
    isnothing(coulomb_expansion) && push!(missing, :coulomb_expansion)
    return Tuple(missing)
end

function _pqs_source_box_route_driver_complete_core_shell_source_plan_payload(
    parent,
    transforms,
    recipe,
)
    recipe.route_family === :pqs_source_box || return (;
        status = :not_applicable_complete_core_shell_source_plan_non_pqs_route,
        blocker = nothing,
        region_plan = nothing,
        source_plan = nothing,
        coulomb_expansion = nothing,
        missing_inputs = (),
        summary = (;
            status = :not_applicable_complete_core_shell_source_plan_non_pqs_route,
            blocker = nothing,
            region_plan_status = :not_applicable,
            source_plan_status = :not_applicable,
        ),
    )

    terminal_route_state =
        hasproperty(transforms, :terminal_route_state) ?
        transforms.terminal_route_state :
        nothing
    shellification_plan =
        !isnothing(terminal_route_state) &&
        hasproperty(terminal_route_state, :shellification_plan) ?
        terminal_route_state.shellification_plan :
        nothing
    lowering_plan =
        !isnothing(terminal_route_state) &&
        hasproperty(terminal_route_state, :lowering_plan) ?
        terminal_route_state.lowering_plan :
        nothing
    parent_axis_bundle_object =
        hasproperty(parent, :parent_axis_bundle_object) ?
        parent.parent_axis_bundle_object :
        nothing

    missing_inputs = Symbol[]
    shellification_plan isa CartesianShellification.ShellificationPlan ||
        push!(missing_inputs, :terminal_shellification_plan)
    lowering_plan isa CartesianTerminalLowering.TerminalLoweringPlan ||
        push!(missing_inputs, :terminal_lowering_plan)
    isnothing(parent_axis_bundle_object) &&
        push!(missing_inputs, :parent_axis_bundle_object)
    if !isempty(missing_inputs)
        status = :blocked_missing_complete_core_shell_source_plan_inputs
        return (;
            status,
            blocker = :missing_complete_core_shell_source_plan_inputs,
            region_plan = nothing,
            source_plan = nothing,
            coulomb_expansion = nothing,
            missing_inputs = Tuple(missing_inputs),
            summary = (;
                status,
                blocker = :missing_complete_core_shell_source_plan_inputs,
                region_plan_status = :not_available,
                source_plan_status = :not_available,
            ),
        )
    end

    expansion = coulomb_gaussian_expansion(doacc = false)
    bond_axis =
        hasproperty(parent, :bond_axis) && parent.bond_axis in (:x, :y, :z) ?
        parent.bond_axis :
        :z
    metadata = (;
        source = :pqs_source_box_route_driver_complete_core_shell_source_plan_payload,
        route_kind = recipe.route_kind,
        route_family = recipe.route_family,
        shellification_backed_geometry = true,
        explicit_box_bridge = false,
    )
    try
        region_plan = pqs_multilayer_shell_region_plan(
            shellification_plan,
            lowering_plan;
            metadata,
        )
        if region_plan.status !== :available_pqs_multilayer_shell_region_plan
            status = :blocked_complete_core_shell_region_plan
            return (;
                status,
                blocker = region_plan.blocker,
                region_plan = nothing,
                source_plan = nothing,
                coulomb_expansion = expansion,
                missing_inputs = (:pqs_multilayer_shell_region_plan,),
                summary = (;
                    status,
                    blocker = region_plan.blocker,
                    region_plan_status = region_plan.status,
                    source_plan_status = :not_available,
                ),
            )
        end

        source_plan = pqs_multilayer_shell_source_plan(
            parent_axis_bundle_object,
            region_plan;
            bond_axis,
            term_coefficients = Float64.(expansion.coefficients),
            metadata,
        )
        source_plan.status === :available_pqs_multilayer_shell_source_plan ||
            return (;
                status = :blocked_complete_core_shell_source_plan,
                blocker = source_plan.blocker,
                region_plan,
                source_plan = nothing,
                coulomb_expansion = expansion,
                missing_inputs = (:pqs_multilayer_shell_source_plan,),
                summary = (;
                    status = :blocked_complete_core_shell_source_plan,
                    blocker = source_plan.blocker,
                    region_plan_status = region_plan.status,
                    source_plan_status = source_plan.status,
                ),
            )

        return (;
            status = :available_complete_core_shell_source_plan_payload,
            blocker = nothing,
            region_plan,
            source_plan,
            coulomb_expansion = expansion,
            missing_inputs = (),
            summary = (;
                status = :available_complete_core_shell_source_plan_payload,
                blocker = nothing,
                region_plan_status = region_plan.status,
                source_plan_status = source_plan.status,
                source_kind = source_plan.source_kind,
                layer_count = source_plan.layer_count,
                core_support_count = length(source_plan.core_support_indices),
                shell_support_count = length(source_plan.shell_support_indices),
            ),
        )
    catch error
        error isa ArgumentError || rethrow()
        status = :blocked_complete_core_shell_source_plan_error
        return (;
            status,
            blocker = :complete_core_shell_source_plan_error,
            region_plan = nothing,
            source_plan = nothing,
            coulomb_expansion = expansion,
            missing_inputs = (:pqs_multilayer_shell_region_plan,
                              :pqs_multilayer_shell_source_plan),
            summary = (;
                status,
                blocker = :complete_core_shell_source_plan_error,
                region_plan_status = :blocked_argument_error,
                source_plan_status = :blocked_argument_error,
                error_message = sprint(showerror, error),
            ),
        )
    end
end

function _pqs_source_box_route_driver_descriptor_property(
    object,
    key::Symbol,
    default = nothing,
)
    object isa NamedTuple && haskey(object, key) && return getfield(object, key)
    hasproperty(object, key) && return getproperty(object, key)
    return default
end

function _pqs_source_box_route_driver_complete_core_shell_center_records(parent)
    center_table =
        hasproperty(parent, :center_table) ? parent.center_table : nothing
    isnothing(center_table) && return nothing, (:center_table,)
    centers = Tuple(center_table)
    isempty(centers) && return nothing, (:center_table,)
    records = Tuple(
        (;
            center_key = center.center_key,
            center_index = center.center_index,
            location = center.location,
            charge = Float64(center.nuclear_charge),
            nuclear_charge = Float64(center.nuclear_charge),
        ) for center in centers
    )
    return records, ()
end

function _pqs_source_box_route_driver_complete_core_shell_axis_layers(
    parent_axis_bundle_object,
)
    missing = Symbol[]
    layers = Dict{Symbol,Any}()
    for axis in (:x, :y, :z)
        pgdg = _nested_axis_pgdg(parent_axis_bundle_object, axis)
        if isnothing(pgdg)
            push!(missing, Symbol(axis, "_pgdg_intermediate"))
            continue
        end
        layer =
            _pqs_source_box_route_driver_descriptor_property(
                pgdg,
                :auxiliary_layer,
            )
        if isnothing(layer)
            push!(missing, Symbol(axis, "_auxiliary_layer"))
            continue
        end
        layers[axis] = layer
    end
    isempty(missing) || return nothing, Tuple(missing)
    return (x = layers[:x], y = layers[:y], z = layers[:z]), ()
end

function _pqs_source_box_route_driver_complete_core_shell_h1_payload(
    parent,
    source_plan_payload,
    recipe,
)
    source_plan = source_plan_payload.source_plan
    if isnothing(source_plan)
        status = :blocked_missing_complete_core_shell_source_plan_for_h1
        return (;
            status,
            blocker = :missing_complete_core_shell_source_plan,
            final_basis = nothing,
            h1_payload = nothing,
            missing_inputs = (:pqs_multilayer_shell_source_plan,),
            summary = (;
                status,
                blocker = :missing_complete_core_shell_source_plan,
                final_basis_status = :not_available,
                h1_payload_status = :not_available,
            ),
        )
    end

    center_records, missing_centers =
        _pqs_source_box_route_driver_complete_core_shell_center_records(parent)
    axis_layers, missing_axis_layers =
        _pqs_source_box_route_driver_complete_core_shell_axis_layers(
            parent.parent_axis_bundle_object,
        )
    missing_inputs = (missing_centers..., missing_axis_layers...)
    if !isempty(missing_inputs)
        status = :blocked_missing_complete_core_shell_h1_inputs
        return (;
            status,
            blocker = :missing_complete_core_shell_h1_inputs,
            final_basis = nothing,
            h1_payload = nothing,
            missing_inputs,
            summary = (;
                status,
                blocker = :missing_complete_core_shell_h1_inputs,
                final_basis_status = :not_available,
                h1_payload_status = :not_available,
            ),
        )
    end

    metadata = (;
        source = :pqs_source_box_route_driver_complete_core_shell_h1_payload,
        route_kind = recipe.route_kind,
        route_family = recipe.route_family,
        shellification_backed_geometry = true,
        explicit_box_bridge = false,
    )
    try
        final_basis = pqs_multilayer_complete_core_shell_final_basis(
            source_plan;
            metadata,
        )
        final_basis_status = get(final_basis, :status, nothing)
        if final_basis_status !== :available_pqs_complete_core_shell_final_basis
            status = :blocked_complete_core_shell_final_basis
            return (;
                status,
                blocker = get(final_basis, :blocker, :complete_core_shell_final_basis_blocked),
                final_basis = nothing,
                h1_payload = nothing,
                missing_inputs = (:pqs_multilayer_complete_core_shell_final_basis,),
                summary = (;
                    status,
                    blocker = get(final_basis, :blocker, :complete_core_shell_final_basis_blocked),
                    final_basis_status,
                    h1_payload_status = :not_available,
                ),
            )
        end

        h1_payload = pqs_multilayer_complete_core_shell_h1_payload(
            source_plan;
            final_basis,
            coulomb_expansion = source_plan_payload.coulomb_expansion,
            center_records,
            axis_layers,
            metadata,
        )
        h1_payload_status = get(h1_payload, :status, nothing)
        h1_payload_status ===
        :materialized_pqs_multilayer_complete_core_shell_h1_payload || return (;
            status = :blocked_complete_core_shell_h1_payload,
            blocker = get(h1_payload, :blocker, :complete_core_shell_h1_payload_blocked),
            final_basis,
            h1_payload = nothing,
            missing_inputs = (:pqs_multilayer_complete_core_shell_h1_payload,),
            summary = (;
                status = :blocked_complete_core_shell_h1_payload,
                blocker = get(h1_payload, :blocker, :complete_core_shell_h1_payload_blocked),
                final_basis_status,
                h1_payload_status,
            ),
        )

        return (;
            status = :available_complete_core_shell_h1_payload,
            blocker = nothing,
            final_basis,
            h1_payload,
            missing_inputs = (),
            summary = (;
                status = :available_complete_core_shell_h1_payload,
                blocker = nothing,
                final_basis_status,
                h1_payload_status,
                final_dimension = h1_payload.summary.final_dimension,
                h1_energy = h1_payload.h1.lowest_energy,
                center_count = length(center_records),
            ),
        )
    catch error
        error isa ArgumentError || rethrow()
        status = :blocked_complete_core_shell_h1_payload_error
        return (;
            status,
            blocker = :complete_core_shell_h1_payload_error,
            final_basis = nothing,
            h1_payload = nothing,
            missing_inputs = (:pqs_multilayer_complete_core_shell_final_basis,
                              :pqs_multilayer_complete_core_shell_h1_payload),
            summary = (;
                status,
                blocker = :complete_core_shell_h1_payload_error,
                final_basis_status = :blocked_argument_error,
                h1_payload_status = :blocked_argument_error,
                error_message = sprint(showerror, error),
            ),
        )
    end
end

function _pqs_source_box_route_driver_complete_core_shell_density_inputs(
    source_plan,
    coulomb_expansion;
    metadata = (;),
)
    source_plan_kind =
        isnothing(source_plan) ?
        nothing :
        _pqs_source_box_route_driver_descriptor_property(
            source_plan,
            :object_kind,
        )
    source_plan_status =
        isnothing(source_plan) ?
        nothing :
        _pqs_source_box_route_driver_descriptor_property(
            source_plan,
            :status,
        )
    source_plan_blocker =
        isnothing(source_plan) ?
        nothing :
        _pqs_source_box_route_driver_descriptor_property(
            source_plan,
            :blocker,
        )
    bundles =
        !isnothing(source_plan) && hasproperty(source_plan, :bundles) ?
        source_plan.bundles :
        nothing

    missing_inputs = Symbol[]
    if source_plan_kind !== :pqs_multilayer_shell_source_plan ||
       source_plan_status !== :available_pqs_multilayer_shell_source_plan
        push!(missing_inputs, :pqs_multilayer_shell_source_plan)
    end
    isnothing(bundles) &&
        push!(missing_inputs, :pqs_multilayer_shell_source_plan_bundles)
    isnothing(coulomb_expansion) &&
        push!(missing_inputs, :coulomb_expansion)
    if !isnothing(coulomb_expansion) && !hasproperty(coulomb_expansion, :coefficients)
        push!(missing_inputs, :coulomb_expansion_coefficients)
    end
    if !isempty(missing_inputs)
        status = :blocked_missing_complete_core_shell_density_inputs
        return (;
            status,
            blocker = :missing_complete_core_shell_density_inputs,
            axis_weights = nothing,
            raw_pair_factor_terms = nothing,
            missing_inputs = Tuple(missing_inputs),
            summary = (;
                status,
                blocker = :missing_complete_core_shell_density_inputs,
                source_plan_status,
                source_plan_blocker,
                source_plan_kind,
                axis_weights_available = false,
                raw_pair_factor_terms_available = false,
            ),
            metadata = merge(
                NamedTuple(metadata),
                (;
                    source =
                        :pqs_source_box_route_driver_complete_core_shell_density_inputs,
                    source_plan_status,
                ),
            ),
        )
    end

    try
        expected_term_count = length(coulomb_expansion.coefficients)
        provenance =
            CartesianContractedParentMetrics._pqs_source_box_ida_factor_provenance(
                bundles;
                expected_term_count,
            )
        return (;
            status = :available_complete_core_shell_density_inputs,
            blocker = nothing,
            axis_weights = provenance.axis_weights,
            raw_pair_factor_terms = provenance.raw_axis_pair_factor_terms,
            missing_inputs = (),
            summary = (;
                status = :available_complete_core_shell_density_inputs,
                blocker = nothing,
                source_plan_status,
                source_plan_kind,
                provenance_status = get(
                    provenance,
                    :object_kind,
                    :pqs_source_box_ida_factor_provenance,
                ),
                term_count = provenance.term_count,
                factor_dimensions = provenance.factor_dimensions,
                axis_weights_available = true,
                raw_pair_factor_terms_available = true,
                private_diagnostic_only = true,
                retained_pqs_weights_used = false,
                density_normalized_pair_terms_used_as_authority = false,
            ),
            metadata = merge(
                NamedTuple(metadata),
                (;
                    source =
                        :pqs_source_box_route_driver_complete_core_shell_density_inputs,
                    provenance_source = :pqs_source_box_ida_factor_provenance,
                    expected_term_count,
                    source_plan_status,
                ),
            ),
        )
    catch error
        error isa ArgumentError || error isa DimensionMismatch || rethrow()
        status = :blocked_complete_core_shell_density_inputs_error
        return (;
            status,
            blocker = :complete_core_shell_density_inputs_error,
            axis_weights = nothing,
            raw_pair_factor_terms = nothing,
            missing_inputs = (:axis_weights, :raw_pair_factor_terms),
            summary = (;
                status,
                blocker = :complete_core_shell_density_inputs_error,
                source_plan_status,
                source_plan_blocker,
                source_plan_kind,
                axis_weights_available = false,
                raw_pair_factor_terms_available = false,
                error_message = sprint(showerror, error),
            ),
            metadata = merge(
                NamedTuple(metadata),
                (;
                    source =
                        :pqs_source_box_route_driver_complete_core_shell_density_inputs,
                    provenance_source = :pqs_source_box_ida_factor_provenance,
                    source_plan_status,
                ),
            ),
        )
    end
end

function _pqs_source_box_route_driver_blocked_complete_core_shell_h1_j_payload(;
    route_family,
    status,
    blocker,
    missing_inputs,
    final_basis = nothing,
    h1_payload = nothing,
    metadata = (;),
)
    final_dimension =
        !isnothing(h1_payload) && hasproperty(h1_payload, :summary) &&
        hasproperty(h1_payload.summary, :final_dimension) ?
        h1_payload.summary.final_dimension :
        nothing
    h1_energy =
        !isnothing(h1_payload) && hasproperty(h1_payload, :h1) &&
        hasproperty(h1_payload.h1, :lowest_energy) ?
        h1_payload.h1.lowest_energy :
        nothing
    summary = (;
        status,
        blocker,
        final_dimension,
        h1_energy,
        self_coulomb = nothing,
        density_gauge = nothing,
        missing_inputs,
        support_density_input_source = :not_materialized,
        h1_orbital_source = :not_materialized,
        signed_final_weight_division_used = false,
        raw_no_division_used = false,
        density_normalized_pair_terms_used_as_authority = false,
        driver_route_materialized = false,
        ida_data_materialized = false,
        density_density_materialized = false,
        rhf_materialized = false,
        gto_materialized = false,
        exports_materialized = false,
        artifacts_materialized = false,
    )
    return (;
        object_kind = :pqs_complete_core_shell_h1_j_diagnostic_route_payload,
        status,
        blocker,
        route_family,
        materialized = false,
        complete_core_shell_route_payload_available =
            !isnothing(final_basis) && !isnothing(h1_payload),
        region_plan = nothing,
        source_plan = nothing,
        final_basis = nothing,
        h1_payload = nothing,
        h1_j_payload = nothing,
        missing_inputs,
        summary,
        metadata = merge(
            NamedTuple(metadata),
            (;
                source =
                    :pqs_source_box_route_driver_complete_core_shell_h1_j_diagnostic_payload,
                driver_route_materialized = false,
                report_placeholder = false,
            ),
        ),
    )
end

function _pqs_source_box_route_driver_complete_core_shell_h1_j_diagnostic_payload(;
    route_family,
    region_plan = nothing,
    source_plan = nothing,
    final_basis = nothing,
    h1_payload = nothing,
    axis_weights = nothing,
    raw_pair_factor_terms = nothing,
    coulomb_expansion = nothing,
    metadata = (;),
)
    if route_family !== :pqs_source_box
        return _pqs_source_box_route_driver_blocked_complete_core_shell_h1_j_payload(
            ;
            route_family,
            status =
                :not_applicable_complete_core_shell_h1_j_non_pqs_source_box_route,
            blocker = nothing,
            missing_inputs = (),
            metadata,
        )
    end

    missing_inputs =
        _pqs_source_box_route_driver_complete_core_shell_h1_j_missing_inputs(
            ;
            region_plan,
            source_plan,
            final_basis,
            h1_payload,
            axis_weights,
            raw_pair_factor_terms,
            coulomb_expansion,
        )
    if !isempty(missing_inputs)
        return _pqs_source_box_route_driver_blocked_complete_core_shell_h1_j_payload(
            ;
            route_family,
            status = :blocked_missing_complete_core_shell_h1_j_route_inputs,
            blocker = :missing_complete_core_shell_h1_j_route_inputs,
            missing_inputs,
            final_basis,
            h1_payload,
            metadata,
        )
    end

    h1_j_payload = pqs_multilayer_complete_core_shell_h1_j_payload(
        source_plan;
        final_basis,
        h1_payload,
        axis_weights,
        raw_pair_factor_terms,
        coulomb_expansion,
        metadata = merge(
            NamedTuple(metadata),
            (;
                source =
                    :pqs_source_box_route_driver_complete_core_shell_h1_j_diagnostic_payload,
            ),
        ),
    )
    h1_j_summary = h1_j_payload.summary
    summary = merge(
        h1_j_summary,
        (;
            driver_route_materialized =
                h1_j_summary.status ===
                :materialized_pqs_multilayer_complete_core_shell_h1_j_payload,
            report_placeholder = false,
        ),
    )
    return (;
        object_kind = :pqs_complete_core_shell_h1_j_diagnostic_route_payload,
        status = summary.status,
        blocker = summary.blocker,
        route_family,
        materialized = summary.driver_route_materialized,
        complete_core_shell_route_payload_available = true,
        region_plan,
        source_plan,
        final_basis,
        h1_payload,
        h1_j_payload,
        missing_inputs = (),
        summary,
        metadata = merge(
            NamedTuple(metadata),
            (;
                source =
                    :pqs_source_box_route_driver_complete_core_shell_h1_j_diagnostic_payload,
                driver_route_materialized = summary.driver_route_materialized,
                report_placeholder = false,
            ),
        ),
    )
end

function _pqs_source_box_route_driver_complete_core_shell_ham_missing_blocker(
    missing_inputs,
)
    isempty(missing_inputs) && return nothing
    length(missing_inputs) > 1 && return :missing_complete_core_shell_ham_inputs
    only_missing = only(missing_inputs)
    only_missing === :pqs_multilayer_complete_core_shell_final_basis &&
        return :missing_complete_core_shell_final_basis
    only_missing === :pqs_multilayer_complete_core_shell_h1_payload &&
        return :missing_complete_core_shell_h1_payload
    only_missing === :pqs_complete_core_shell_final_one_electron_hamiltonian &&
        return :missing_complete_core_shell_final_hamiltonian
    only_missing === :complete_core_shell_density_inputs &&
        return :missing_complete_core_shell_density_inputs
    only_missing === :complete_core_shell_h1_j_diagnostic_payload &&
        return :missing_complete_core_shell_h1_j_diagnostic_payload
    only_missing === :pqs_complete_core_shell_pre_final_density_interaction &&
        return :missing_complete_core_shell_density_interaction
    return :missing_complete_core_shell_ham_inputs
end

function _pqs_source_box_route_driver_blocked_complete_core_shell_ham_payload(;
    route_family,
    status,
    blocker,
    missing_inputs,
    source_payload = nothing,
    final_basis = nothing,
    h1_payload = nothing,
    density_inputs = nothing,
    h1_j_payload = nothing,
    metadata = (;),
)
    final_dimension =
        !isnothing(final_basis) && hasproperty(final_basis, :final_retained_count) ?
        final_basis.final_retained_count :
        nothing
    summary = (;
        status,
        blocker,
        route_family,
        final_dimension,
        electron_electron_representation = :not_materialized,
        density_gauge = nothing,
        raw_pair_factor_convention = nothing,
        support_row_order = nothing,
        signed_final_weight_division_used = false,
        raw_no_division_used = false,
        density_normalized_pair_terms_used_as_authority = false,
        public_api = false,
        exports_materialized = false,
        artifacts_materialized = false,
        rhf_product_surface = false,
        serious_hf_claim = false,
        missing_inputs,
    )
    return (;
        object_kind = :pqs_source_box_complete_core_shell_ham_payload,
        status,
        blocker,
        route_family,
        source_payload,
        source_plan_summary =
            !isnothing(source_payload) && hasproperty(source_payload, :summary) ?
            source_payload.summary :
            nothing,
        final_basis,
        final_basis_summary = nothing,
        h1_payload,
        one_body_hamiltonian = nothing,
        density_inputs,
        density_input_summary =
            !isnothing(density_inputs) && hasproperty(density_inputs, :summary) ?
            density_inputs.summary :
            nothing,
        h1_j_payload,
        density_interaction = nothing,
        electron_electron_representation = :not_materialized,
        coulomb_expansion =
            !isnothing(source_payload) && hasproperty(source_payload, :coulomb_expansion) ?
            source_payload.coulomb_expansion :
            nothing,
        coulomb_expansion_summary = nothing,
        dimension_summary = (; final_dimension,),
        ordering_summary = (; support_row_order = nothing,),
        convention_labels = (;
            electron_electron_representation = :not_materialized,
            density_gauge = nothing,
            raw_pair_factor_convention = nothing,
            support_row_order = nothing,
            signed_final_weight_division_used = false,
            raw_no_division_used = false,
            density_normalized_pair_terms_used_as_authority = false,
            fixed_block_pair_data_authority_used = false,
        ),
        missing_inputs,
        summary,
        metadata = merge(
            NamedTuple(metadata),
            (;
                source =
                    :pqs_source_box_route_driver_complete_core_shell_ham_payload,
                public_api = false,
                exports_materialized = false,
                artifacts_materialized = false,
                rhf_product_surface = false,
                serious_hf_claim = false,
            ),
        ),
    )
end

function _pqs_source_box_route_driver_complete_core_shell_ham_payload(;
    route_family,
    source_payload = nothing,
    final_basis = nothing,
    h1_payload = nothing,
    density_inputs = nothing,
    h1_j_payload = nothing,
    metadata = (;),
)
    if route_family !== :pqs_source_box
        return _pqs_source_box_route_driver_blocked_complete_core_shell_ham_payload(
            ;
            route_family,
            status = :not_applicable_complete_core_shell_ham_non_pqs_source_box_route,
            blocker = nothing,
            missing_inputs = (),
            source_payload,
            final_basis,
            h1_payload,
            density_inputs,
            h1_j_payload,
            metadata,
        )
    end

    one_body_hamiltonian =
        _pqs_source_box_route_driver_descriptor_property(
            h1_payload,
            :final_hamiltonian,
        )
    inner_h1_j_payload =
        _pqs_source_box_route_driver_descriptor_property(
            h1_j_payload,
            :h1_j_payload,
        )
    density_interaction =
        _pqs_source_box_route_driver_descriptor_property(
            h1_j_payload,
            :density_interaction,
        )
    isnothing(density_interaction) &&
        (density_interaction =
            _pqs_source_box_route_driver_descriptor_property(
                inner_h1_j_payload,
                :density_interaction,
            ))

    missing_inputs = Symbol[]
    (
        _pqs_source_box_route_driver_descriptor_property(final_basis, :status) ===
        :available_pqs_complete_core_shell_final_basis
    ) || push!(
        missing_inputs,
        :pqs_multilayer_complete_core_shell_final_basis,
    )
    (
        _pqs_source_box_route_driver_descriptor_property(h1_payload, :status) ===
        :materialized_pqs_multilayer_complete_core_shell_h1_payload
    ) || push!(
        missing_inputs,
        :pqs_multilayer_complete_core_shell_h1_payload,
    )
    (
        _pqs_source_box_route_driver_descriptor_property(
            one_body_hamiltonian,
            :status,
        ) === :materialized_pqs_complete_core_shell_final_one_electron_hamiltonian
    ) || push!(
        missing_inputs,
        :pqs_complete_core_shell_final_one_electron_hamiltonian,
    )
    (
        _pqs_source_box_route_driver_descriptor_property(density_inputs, :status) ===
        :available_complete_core_shell_density_inputs
    ) || push!(missing_inputs, :complete_core_shell_density_inputs)
    (
        _pqs_source_box_route_driver_descriptor_property(h1_j_payload, :status) ===
        :materialized_pqs_multilayer_complete_core_shell_h1_j_payload
    ) || push!(missing_inputs, :complete_core_shell_h1_j_diagnostic_payload)
    (
        _pqs_source_box_route_driver_descriptor_property(
            density_interaction,
            :status,
        ) === :materialized_pqs_complete_core_shell_pre_final_density_interaction
    ) || push!(
        missing_inputs,
        :pqs_complete_core_shell_pre_final_density_interaction,
    )

    if !isempty(missing_inputs)
        blocker =
            _pqs_source_box_route_driver_complete_core_shell_ham_missing_blocker(
                Tuple(missing_inputs),
            )
        return _pqs_source_box_route_driver_blocked_complete_core_shell_ham_payload(
            ;
            route_family,
            status = :blocked_complete_core_shell_ham_payload,
            blocker,
            missing_inputs = Tuple(missing_inputs),
            source_payload,
            final_basis,
            h1_payload,
            density_inputs,
            h1_j_payload,
            metadata,
        )
    end

    support_row_order = final_basis.support_row_order
    density_gauge = density_interaction.density_gauge
    raw_pair_factor_convention =
        hasproperty(density_interaction, :metadata) &&
        hasproperty(density_interaction.metadata, :raw_pair_factor_convention) ?
        density_interaction.metadata.raw_pair_factor_convention :
        :raw_numerator
    convention_labels = (;
        electron_electron_representation = :pre_final_density_interaction,
        interaction_model = :density_density_pre_final_gauge,
        density_gauge,
        raw_pair_factor_convention,
        support_row_order,
        weight_application_stage =
            hasproperty(density_interaction, :metadata) &&
            hasproperty(density_interaction.metadata, :weight_application_stage) ?
            density_interaction.metadata.weight_application_stage :
            :pre_final_density_interaction_boundary,
        final_orbital_consumption_rule =
            hasproperty(density_interaction, :metadata) &&
            hasproperty(density_interaction.metadata, :final_orbital_consumption_rule) ?
            density_interaction.metadata.final_orbital_consumption_rule :
            :combined_lowdin_cleanup_times_final_coefficients,
        signed_final_weight_division_used =
            density_interaction.signed_final_weight_division_used,
        raw_no_division_used = density_interaction.raw_no_division_used,
        density_normalized_pair_terms_used_as_authority = false,
        fixed_block_pair_data_authority_used =
            density_interaction.fixed_block_pair_data_authority_used,
        public_api = false,
        exports_materialized = false,
        artifacts_materialized = false,
        rhf_product_surface = false,
        serious_hf_claim = false,
    )
    final_basis_summary = (;
        status = final_basis.status,
        blocker = final_basis.blocker,
        final_retained_count = final_basis.final_retained_count,
        core_support_count = final_basis.core_support_count,
        shell_support_count = final_basis.shell_support_count,
        shell_final_retained_count = final_basis.shell_final_retained_count,
        support_row_order,
        final_overlap_identity_error = final_basis.final_overlap_identity_error,
        final_overlap_is_identity = final_basis.final_overlap_is_identity,
    )
    density_input_summary = density_inputs.summary
    coulomb_expansion = source_payload.coulomb_expansion
    coulomb_expansion_summary = (;
        coefficient_count = length(coulomb_expansion.coefficients),
        coefficients_positive = all(>(0.0), coulomb_expansion.coefficients),
    )
    dimension_summary = (;
        final_dimension = final_basis.final_retained_count,
        hamiltonian_shape = one_body_hamiltonian.hamiltonian_matrix_shape,
        pre_final_dimension = density_interaction.pre_final_weight_count,
        pre_final_pair_matrix_shape =
            density_interaction.pre_final_pair_matrix_shape,
    )
    ordering_summary = (;
        final_basis_support_row_order = support_row_order,
        density_interaction_support_row_order =
            density_interaction.support_row_order,
        final_matrix_order = :final_basis_column_order,
        pre_final_density_order = :pre_final_complete_core_shell_order,
    )
    center_summaries =
        hasproperty(one_body_hamiltonian, :center_summaries) ?
        one_body_hamiltonian.center_summaries :
        ()
    summary = (;
        status = :materialized_pqs_source_box_complete_core_shell_ham_payload,
        blocker = nothing,
        route_family,
        final_dimension = final_basis.final_retained_count,
        one_body_hamiltonian_status = one_body_hamiltonian.status,
        density_interaction_status = density_interaction.status,
        electron_electron_representation = :pre_final_density_interaction,
        density_gauge,
        raw_pair_factor_convention,
        support_row_order,
        center_count = length(center_summaries),
        signed_final_weight_division_used =
            convention_labels.signed_final_weight_division_used,
        raw_no_division_used = convention_labels.raw_no_division_used,
        density_normalized_pair_terms_used_as_authority = false,
        public_api = false,
        exports_materialized = false,
        artifacts_materialized = false,
        rhf_product_surface = false,
        serious_hf_claim = false,
    )
    return (;
        object_kind = :pqs_source_box_complete_core_shell_ham_payload,
        status = summary.status,
        blocker = nothing,
        route_family,
        source_payload,
        source_plan_summary = source_payload.summary,
        final_basis,
        final_basis_summary,
        h1_payload,
        one_body_hamiltonian,
        density_inputs,
        density_input_summary,
        h1_j_payload,
        density_interaction,
        electron_electron_representation = :pre_final_density_interaction,
        coulomb_expansion,
        coulomb_expansion_summary,
        center_summaries,
        dimension_summary,
        ordering_summary,
        convention_labels,
        missing_inputs = (),
        summary,
        metadata = merge(
            NamedTuple(metadata),
            (;
                source =
                    :pqs_source_box_route_driver_complete_core_shell_ham_payload,
                electron_electron_representation =
                    :pre_final_density_interaction,
                density_gauge,
                raw_pair_factor_convention,
                support_row_order,
                public_api = false,
                exports_materialized = false,
                artifacts_materialized = false,
                rhf_product_surface = false,
                serious_hf_claim = false,
            ),
        ),
    )
end

struct _PQSCompleteCoreShellDiagnosticRoutePayload
    status::Symbol
    blocker
    route_family::Symbol
    source_payload
    final_basis
    h1_payload
    density_inputs
    h1_j_payload
    complete_core_shell_ham_payload
    missing_inputs::Tuple
    summary
    metadata
end

function _pqs_source_box_route_driver_complete_core_shell_diagnostic_route_payload(
    parent,
    transforms,
    recipe,
    route_skeleton,
)
    source_payload =
        _pqs_source_box_route_driver_complete_core_shell_source_plan_payload(
            parent,
            transforms,
            recipe,
        )
    h1_payload =
        _pqs_source_box_route_driver_complete_core_shell_h1_payload(
            parent,
            source_payload,
            recipe,
        )
    density_inputs =
        _pqs_source_box_route_driver_complete_core_shell_density_inputs(
            source_payload.source_plan,
            source_payload.coulomb_expansion;
            metadata = (;
                source = :cartesian_assembly,
                route_kind = recipe.route_kind,
                complete_core_shell_source_plan_status = source_payload.status,
            ),
        )
    h1_j_payload =
        _pqs_source_box_route_driver_complete_core_shell_h1_j_diagnostic_payload(
            ;
            route_family = route_skeleton.route_family,
            region_plan = source_payload.region_plan,
            source_plan = source_payload.source_plan,
            final_basis = h1_payload.final_basis,
            h1_payload = h1_payload.h1_payload,
            axis_weights = density_inputs.axis_weights,
            raw_pair_factor_terms = density_inputs.raw_pair_factor_terms,
            coulomb_expansion = source_payload.coulomb_expansion,
            metadata = (;
                source = :cartesian_assembly,
                route_kind = recipe.route_kind,
                complete_core_shell_source_plan_status = source_payload.status,
                complete_core_shell_h1_payload_status = h1_payload.status,
                complete_core_shell_density_inputs_status = density_inputs.status,
            ),
        )
    complete_core_shell_ham_payload =
        _pqs_source_box_route_driver_complete_core_shell_ham_payload(
            ;
            route_family = route_skeleton.route_family,
            source_payload,
            final_basis = h1_payload.final_basis,
            h1_payload = h1_payload.h1_payload,
            density_inputs,
            h1_j_payload,
            metadata = (;
                source = :cartesian_assembly,
                route_kind = recipe.route_kind,
                complete_core_shell_source_plan_status = source_payload.status,
                complete_core_shell_h1_payload_status = h1_payload.status,
                complete_core_shell_density_inputs_status = density_inputs.status,
                complete_core_shell_h1_j_diagnostic_status = h1_j_payload.status,
            ),
        )
    metadata = (;
        source = :pqs_source_box_route_driver_complete_core_shell_diagnostic_route_payload,
        route_kind = recipe.route_kind,
        route_family = route_skeleton.route_family,
        complete_core_shell_source_plan_status = source_payload.status,
        complete_core_shell_h1_payload_status = h1_payload.status,
        complete_core_shell_density_inputs_status = density_inputs.status,
        complete_core_shell_h1_j_diagnostic_status = h1_j_payload.status,
        complete_core_shell_ham_payload_status =
            complete_core_shell_ham_payload.status,
        driver_route_materialized = h1_j_payload.summary.driver_route_materialized,
        report_placeholder = false,
    )
    return _PQSCompleteCoreShellDiagnosticRoutePayload(
        h1_j_payload.status,
        h1_j_payload.blocker,
        route_skeleton.route_family,
        source_payload,
        h1_payload.final_basis,
        h1_payload.h1_payload,
        density_inputs,
        h1_j_payload,
        complete_core_shell_ham_payload,
        h1_j_payload.missing_inputs,
        h1_j_payload.summary,
        metadata,
    )
end

function _pqs_source_box_route_driver_private_rhf_electron_count(parent, requested)
    !isnothing(requested) && return requested, :explicit
    if parent.system_classification == :one_center && length(parent.nuclear_charges) == 1
        charge = only(parent.nuclear_charges)
        if charge isa Real && isfinite(Float64(charge)) && isinteger(Float64(charge))
            return Int(charge), :one_center_nuclear_charge
        end
        return nothing, :invalid_one_center_nuclear_charge
    end
    return nothing, :not_inferred_non_one_center
end

function _pqs_source_box_route_driver_complete_core_shell_private_rhf_payload(
    parent,
    recipe,
    diagnostic_route_payload,
)
    controls = get(recipe, :private_rhf_inputs, (; run_private_rhf = false))
    get(controls, :run_private_rhf, false) || return nothing
    electron_count, electron_count_source =
        _pqs_source_box_route_driver_private_rhf_electron_count(
            parent,
            get(controls, :private_rhf_electron_count, nothing),
        )
    source_payload = diagnostic_route_payload.source_payload
    input_contract = _pqs_multilayer_complete_core_shell_rhf_input_contract(
        ;
        source_plan = source_payload.source_plan,
        final_basis = diagnostic_route_payload.final_basis,
        h1_payload = diagnostic_route_payload.h1_payload,
        density_inputs = diagnostic_route_payload.density_inputs,
        coulomb_expansion = source_payload.coulomb_expansion,
        electron_count,
        fixture_role = get(controls, :private_rhf_fixture_role, :route_smoke),
        metadata = (; source = :cartesian_assembly, electron_count_source),
    )
    return _pqs_multilayer_complete_core_shell_rhf_scf_payload(
        ;
        input_contract,
        h1_payload = diagnostic_route_payload.h1_payload,
        h1_j_payload = diagnostic_route_payload.h1_j_payload,
        density_interaction =
            diagnostic_route_payload.complete_core_shell_ham_payload.density_interaction,
        mixing_kind = get(controls, :private_rhf_mixing_kind, :fock_diis),
        max_iterations = get(controls, :private_rhf_max_iterations, 25),
        density_atol = get(controls, :private_rhf_density_atol, 1.0e-8),
        energy_atol = get(controls, :private_rhf_energy_atol, 1.0e-10),
        residual_atol = get(controls, :private_rhf_residual_atol, 1.0e-8),
        trace_atol = get(controls, :private_rhf_trace_atol, 1.0e-8),
        idempotency_atol = get(controls, :private_rhf_idempotency_atol, 1.0e-8),
        max_history = get(controls, :private_rhf_max_history, nothing),
        diis_start_iteration =
            get(controls, :private_rhf_diis_start_iteration, 2),
        diis_regularization =
            get(controls, :private_rhf_diis_regularization, 1.0e-12),
        diis_coefficient_max_abs =
            get(controls, :private_rhf_diis_coefficient_max_abs, 25.0),
        metadata = (;
            source = :cartesian_assembly,
            electron_count_source,
            private_diagnostic_only = true,
        )
    )
end

function cartesian_assembly(parent, shells, units, transforms, pairs, recipe)
    route_skeleton = shells.route_skeleton
    route_facts = _pqs_source_box_route_driver_route_facts(route_skeleton)
    contract = _pqs_source_box_route_driver_contract_metadata(recipe)
    low_order_assembly =
        _pqs_source_box_route_driver_assembly_stage_low_order_summary(pairs)
    diatomic_physical_gausslet_target_payload =
        _pqs_source_box_route_driver_diatomic_physical_gausslet_target_payload(
            parent,
            route_skeleton,
            recipe,
        )
    diatomic_physical_gausslet_supplement_request_payload =
        _pqs_source_box_route_driver_diatomic_physical_gausslet_supplement_request_payload(
            parent,
            diatomic_physical_gausslet_target_payload,
        )
    diatomic_physical_gausslet_supplement_representation_payload =
        _pqs_source_box_route_driver_diatomic_physical_gausslet_supplement_representation_payload(
            diatomic_physical_gausslet_supplement_request_payload,
        )
    diatomic_physical_gausslet_supplement_preflight_payload =
        _pqs_source_box_route_driver_diatomic_physical_gausslet_supplement_preflight_payload(
            diatomic_physical_gausslet_target_payload,
            diatomic_physical_gausslet_supplement_request_payload,
            diatomic_physical_gausslet_supplement_representation_payload,
        )
    independent_physical_source_plan_route =
        recipe.route_kind ===
        :bond_aligned_diatomic_independent_pqs_source_box_core_shell
    diatomic_physical_gausslet_source_plan_candidate_payload =
        independent_physical_source_plan_route ?
        nothing :
        _pqs_source_box_route_driver_diatomic_physical_gausslet_source_plan_candidate_payload(
            parent,
            route_skeleton,
            recipe,
            diatomic_physical_gausslet_target_payload,
        )
    diatomic_physical_gausslet_source_plan_payload =
        _pqs_source_box_route_driver_diatomic_physical_gausslet_source_plan_payload(
            parent,
            diatomic_physical_gausslet_target_payload,
            diatomic_physical_gausslet_source_plan_candidate_payload,
        )
    independent_h2_pqs_supplement_support_partition_payload =
        _pqs_source_box_route_driver_independent_h2_pqs_supplement_support_partition_payload(
            diatomic_physical_gausslet_target_payload,
            diatomic_physical_gausslet_source_plan_payload,
        )
    diatomic_physical_gausslet_final_basis_payload =
        _pqs_source_box_route_driver_diatomic_physical_gausslet_final_basis_payload(
            route_skeleton,
            recipe,
            diatomic_physical_gausslet_source_plan_payload,
        )
    diatomic_physical_gausslet_h1_payload =
        _pqs_source_box_route_driver_diatomic_physical_gausslet_h1_payload(
            parent,
            route_skeleton,
            recipe,
            diatomic_physical_gausslet_source_plan_payload,
            diatomic_physical_gausslet_final_basis_payload,
        )
    diatomic_physical_gausslet_h1_j_payload =
        _pqs_source_box_route_driver_diatomic_physical_gausslet_h1_j_payload(
            route_skeleton,
            recipe,
            diatomic_physical_gausslet_source_plan_payload,
            diatomic_physical_gausslet_final_basis_payload,
            diatomic_physical_gausslet_h1_payload,
        )
    run_private_rhf =
        get(get(recipe, :private_rhf_inputs, (;)), :run_private_rhf, false)
    diatomic_physical_gausslet_rhf_input_contract = run_private_rhf ?
        _pqs_source_box_route_driver_diatomic_physical_gausslet_rhf_input_contract(
            parent,
            route_skeleton,
            recipe,
            diatomic_physical_gausslet_source_plan_payload,
            diatomic_physical_gausslet_final_basis_payload,
            diatomic_physical_gausslet_h1_payload,
            diatomic_physical_gausslet_h1_j_payload,
        ) : nothing
    diatomic_physical_gausslet_rhf_execution_payload =
        isnothing(diatomic_physical_gausslet_rhf_input_contract) ?
        nothing :
        _pqs_source_box_route_driver_diatomic_physical_gausslet_rhf_execution_payload(
            diatomic_physical_gausslet_rhf_input_contract,
            diatomic_physical_gausslet_h1_payload,
            diatomic_physical_gausslet_h1_j_payload,
        )
    h2_wl_gausslet_only_reference_candidate =
        _pqs_source_box_route_driver_h2_wl_gausslet_only_reference_candidate(
            parent,
            route_skeleton,
            recipe,
            diatomic_physical_gausslet_target_payload,
            diatomic_physical_gausslet_final_basis_payload,
        )
    complete_core_shell_diagnostic_route_payload =
        _pqs_source_box_route_driver_complete_core_shell_diagnostic_route_payload(
            parent,
            transforms,
            recipe,
            route_skeleton,
        )
    complete_core_shell_h1_j_diagnostic_payload =
        complete_core_shell_diagnostic_route_payload.h1_j_payload
    complete_core_shell_private_rhf_payload =
        _pqs_source_box_route_driver_complete_core_shell_private_rhf_payload(
            parent,
            recipe,
            complete_core_shell_diagnostic_route_payload,
        )
    diatomic_complete_core_shell_support_window_payload =
        _pqs_source_box_route_driver_diatomic_complete_core_shell_support_window_payload(
            parent,
            route_skeleton,
            recipe,
        )
    diatomic_raw_box_route_payload =
        _pqs_source_box_route_driver_diatomic_raw_box_route_payload(
            parent,
            route_skeleton,
            recipe,
            diatomic_complete_core_shell_support_window_payload,
        )
    diatomic_complete_core_shell_source_realization_payload =
        _pqs_source_box_route_driver_diatomic_complete_core_shell_source_realization_payload(
            parent,
            route_skeleton,
            recipe,
            diatomic_complete_core_shell_support_window_payload,
            diatomic_raw_box_route_payload,
        )
    diatomic_complete_core_shell_source_plan_payload =
        _pqs_source_box_route_driver_diatomic_complete_core_shell_source_plan_payload(
            parent,
            route_skeleton,
            recipe,
            diatomic_complete_core_shell_support_window_payload,
            diatomic_raw_box_route_payload,
            diatomic_complete_core_shell_source_realization_payload,
        )
    diatomic_complete_core_shell_final_basis_payload =
        _pqs_source_box_route_driver_diatomic_complete_core_shell_final_basis_payload(
            route_skeleton,
            recipe,
            diatomic_complete_core_shell_support_window_payload,
            diatomic_raw_box_route_payload,
            diatomic_complete_core_shell_source_realization_payload,
            diatomic_complete_core_shell_source_plan_payload,
        )
    diatomic_complete_core_shell_h1_payload =
        _pqs_source_box_route_driver_diatomic_complete_core_shell_h1_payload(
            parent,
            route_skeleton,
            recipe,
            diatomic_complete_core_shell_source_plan_payload,
            diatomic_complete_core_shell_final_basis_payload,
        )
    diatomic_complete_core_shell_ham_input_payload =
        _pqs_source_box_route_driver_diatomic_complete_core_shell_ham_input_payload(
            route_skeleton,
            recipe,
            diatomic_complete_core_shell_source_plan_payload,
            diatomic_complete_core_shell_final_basis_payload,
            diatomic_complete_core_shell_h1_payload,
        )
    diatomic_complete_core_shell_hamiltonian_handoff_payload =
        _pqs_source_box_route_driver_diatomic_complete_core_shell_hamiltonian_handoff_payload(
            parent,
            route_skeleton,
            recipe,
            diatomic_complete_core_shell_source_plan_payload,
            diatomic_complete_core_shell_final_basis_payload,
            diatomic_complete_core_shell_h1_payload,
            diatomic_complete_core_shell_ham_input_payload,
        )
    diatomic_complete_core_shell_hamiltonian_consumer_contract_payload =
        _pqs_source_box_route_driver_diatomic_complete_core_shell_hamiltonian_consumer_contract_payload(
            route_skeleton,
            recipe,
            diatomic_complete_core_shell_hamiltonian_handoff_payload,
        )
    diatomic_complete_core_shell_ham_readiness_payload =
        _pqs_source_box_route_driver_diatomic_complete_core_shell_ham_readiness_payload(
            parent,
            route_skeleton,
            recipe,
            diatomic_complete_core_shell_source_plan_payload,
            diatomic_complete_core_shell_final_basis_payload,
            diatomic_complete_core_shell_h1_payload,
            diatomic_complete_core_shell_ham_input_payload,
            diatomic_complete_core_shell_hamiltonian_handoff_payload,
            diatomic_complete_core_shell_hamiltonian_consumer_contract_payload,
        )

    return (;
        spacing_inputs = shells.spacing_inputs,
        route_skeleton,
        route_facts,
        contract,
        shells,
        units,
        transforms,
        pairs,
        low_order_assembly,
        diatomic_physical_gausslet_target_payload,
        diatomic_physical_gausslet_supplement_request_payload,
        diatomic_physical_gausslet_supplement_representation_payload,
        diatomic_physical_gausslet_supplement_preflight_payload,
        diatomic_physical_gausslet_source_plan_candidate_payload,
        diatomic_physical_gausslet_source_plan_payload,
        independent_h2_pqs_supplement_support_partition_payload,
        diatomic_physical_gausslet_final_basis_payload,
        diatomic_physical_gausslet_h1_payload,
        diatomic_physical_gausslet_h1_j_payload,
        diatomic_physical_gausslet_rhf_input_contract,
        diatomic_physical_gausslet_rhf_execution_payload,
        h2_wl_gausslet_only_reference_candidate,
        complete_core_shell_diagnostic_route_payload,
        diatomic_complete_core_shell_support_window_payload,
        diatomic_raw_box_route_payload,
        diatomic_complete_core_shell_source_realization_payload,
        diatomic_complete_core_shell_source_plan_payload,
        diatomic_complete_core_shell_final_basis_payload,
        diatomic_complete_core_shell_h1_payload,
        diatomic_complete_core_shell_ham_input_payload,
        diatomic_complete_core_shell_hamiltonian_handoff_payload,
        diatomic_complete_core_shell_hamiltonian_consumer_contract_payload,
        diatomic_complete_core_shell_ham_readiness_payload,
        complete_core_shell_h1_j_diagnostic_payload,
        complete_core_shell_private_rhf_payload,
    )
end

function _pqs_source_box_route_driver_report_stage_low_order_route_summary(assembly)
    low_order_assembly =
        hasproperty(assembly, :low_order_assembly) ?
        assembly.low_order_assembly :
        nothing
    isnothing(low_order_assembly) && return (;
        object_kind = :cartesian_report_stage_low_order_route_summary,
        status = :not_available_missing_assembly_stage_summary,
        materialization_status = :blocked_missing_assembly_stage_summary,
        materialization_blocker = :missing_assembly_stage_low_order_summary,
        plan_authority = false,
        active_source_authority = false,
        legacy_source_authority = false,
        summary_only = true,
    )
    return (;
        object_kind = :cartesian_report_stage_low_order_route_summary,
        status = get(low_order_assembly, :status, :available_report_stage_low_order_route_summary),
        materialization_status = get(low_order_assembly, :assembly_materialization_status, :not_available),
        materialization_blocker = get(low_order_assembly, :assembly_blocker, nothing),
        plan_authority = get(low_order_assembly, :plan_authority, false),
        active_source_authority = get(low_order_assembly, :active_source_authority, false),
        legacy_source_authority = get(low_order_assembly, :legacy_source_authority, false),
        summary_only = true,
    )
end

_pqs_source_box_route_driver_payload_summary(payload) =
    isnothing(payload) ? nothing :
    hasproperty(payload, :summary) ? payload.summary :
    payload

function _pqs_source_box_route_driver_complete_core_shell_h1_j_report_fields(
    assembly,
)
    payload =
        hasproperty(assembly, :complete_core_shell_h1_j_diagnostic_payload) ?
        assembly.complete_core_shell_h1_j_diagnostic_payload :
        nothing
    summary = _pqs_source_box_route_driver_payload_summary(payload)
    return (;
        complete_core_shell_h1_j_diagnostic_summary = summary,
        complete_core_shell_h1_j_diagnostic_status =
            isnothing(summary) ? :not_available : get(summary, :status, :available),
        complete_core_shell_h1_j_diagnostic_blocker =
            isnothing(summary) ? :missing_complete_core_shell_h1_j_payload :
            get(summary, :blocker, nothing),
    )
end

function _pqs_source_box_route_driver_complete_core_shell_private_rhf_report_fields(
    assembly,
    recipe,
)
    payload =
        hasproperty(assembly, :complete_core_shell_private_rhf_payload) ?
        assembly.complete_core_shell_private_rhf_payload :
        nothing
    summary = _pqs_source_box_route_driver_payload_summary(payload)
    return (;
        private_rhf_summary = summary,
        private_rhf_status =
            isnothing(summary) ? :not_requested : get(summary, :status, :available),
        private_rhf_blocker =
            isnothing(summary) ? nothing : get(summary, :blocker, nothing),
    )
end

function _pqs_source_box_route_driver_diatomic_complete_core_shell_report_fields(
    assembly,
)
    readiness =
        hasproperty(assembly, :diatomic_complete_core_shell_ham_readiness_payload) ?
        assembly.diatomic_complete_core_shell_ham_readiness_payload :
        nothing
    summary = _pqs_source_box_route_driver_payload_summary(readiness)
    return (;
        diatomic_complete_core_shell_readiness_summary = summary,
    )
end

function _pqs_source_box_route_driver_physical_gausslet_target_report_fields(
    assembly,
)
    target =
        hasproperty(assembly, :diatomic_physical_gausslet_target_payload) ?
        assembly.diatomic_physical_gausslet_target_payload :
        nothing
    source_plan =
        hasproperty(assembly, :diatomic_physical_gausslet_source_plan_payload) ?
        assembly.diatomic_physical_gausslet_source_plan_payload :
        nothing
    final_basis =
        hasproperty(assembly, :diatomic_physical_gausslet_final_basis_payload) ?
        assembly.diatomic_physical_gausslet_final_basis_payload :
        nothing
    h1 =
        hasproperty(assembly, :diatomic_physical_gausslet_h1_payload) ?
        assembly.diatomic_physical_gausslet_h1_payload :
        nothing
    h1_j =
        hasproperty(assembly, :diatomic_physical_gausslet_h1_j_payload) ?
        assembly.diatomic_physical_gausslet_h1_j_payload :
        nothing
    rhf =
        hasproperty(assembly, :diatomic_physical_gausslet_rhf_execution_payload) ?
        assembly.diatomic_physical_gausslet_rhf_execution_payload :
        nothing
    supplement =
        hasproperty(assembly, :diatomic_physical_gausslet_supplement_preflight_payload) ?
        assembly.diatomic_physical_gausslet_supplement_preflight_payload :
        nothing
    target_summary = _pqs_source_box_route_driver_payload_summary(target)
    return (;
        physical_gausslet_target_summary = target_summary,
        physical_gausslet_target_status =
            isnothing(target_summary) ? :not_available :
            get(target_summary, :status, :available),
        physical_gausslet_target_blocker =
            isnothing(target_summary) ? :missing_physical_gausslet_target :
            get(target_summary, :blocker, nothing),
        physical_gausslet_source_plan_summary =
            _pqs_source_box_route_driver_payload_summary(source_plan),
        physical_gausslet_final_basis_summary =
            _pqs_source_box_route_driver_payload_summary(final_basis),
        physical_gausslet_h1_summary =
            _pqs_source_box_route_driver_payload_summary(h1),
        physical_gausslet_h1_j_summary =
            _pqs_source_box_route_driver_payload_summary(h1_j),
        physical_gausslet_private_rhf_summary =
            _pqs_source_box_route_driver_payload_summary(rhf),
        physical_gausslet_supplement_preflight_summary =
            _pqs_source_box_route_driver_payload_summary(supplement),
    )
end

function cartesian_report(system, parent, assembly, recipe)
    standard_setup = parent.standard_setup
    parent_axis = parent.parent_axis
    route_axis_counts = parent.route_axis_counts
    route_skeleton = assembly.route_skeleton
    raw_box = assembly.raw_box
    route_facts = assembly.route_facts
    contract = assembly.contract
    probe_inputs = merge(parent.parent_inputs, assembly.route_inputs)

    system_metadata =
        _pqs_source_box_route_driver_system_metadata(
            standard_setup, route_axis_counts, system)
    recipe_metadata =
        _pqs_source_box_route_driver_recipe_metadata(
            standard_setup, route_axis_counts, parent_axis, raw_box,
            assembly.spacing_inputs, probe_inputs, recipe)
    parent_contract = _pqs_source_box_route_driver_parent_contract(parent)
    parent_description =
        _pqs_source_box_route_driver_parent_description(
            standard_setup, parent_axis, route_axis_counts, parent_contract,
            route_skeleton, raw_box)
    diagnostics =
        _pqs_source_box_route_driver_diagnostics(
            standard_setup, parent_axis, route_axis_counts,
            route_skeleton, raw_box, contract)
    low_order_route_summary =
        _pqs_source_box_route_driver_report_stage_low_order_route_summary(
            assembly,
        )
    complete_core_shell_h1_j_report_fields =
        _pqs_source_box_route_driver_complete_core_shell_h1_j_report_fields(
            assembly,
        )
    complete_core_shell_private_rhf_report_fields =
        _pqs_source_box_route_driver_complete_core_shell_private_rhf_report_fields(
            assembly,
            recipe,
        )
    diatomic_complete_core_shell_report_fields =
        _pqs_source_box_route_driver_diatomic_complete_core_shell_report_fields(
            assembly,
        )
    physical_gausslet_target_report_fields =
        _pqs_source_box_route_driver_physical_gausslet_target_report_fields(
            assembly,
        )

    report = _pqs_source_box_route_driver_report(
        standard_setup, parent, parent_axis, route_axis_counts, raw_box,
        system_metadata, recipe_metadata, parent_contract, parent_description,
        route_skeleton, route_facts, contract, diagnostics,
        low_order_route_summary)
    return merge(
        report,
        complete_core_shell_h1_j_report_fields,
        complete_core_shell_private_rhf_report_fields,
        diatomic_complete_core_shell_report_fields,
        physical_gausslet_target_report_fields,
    )
end

function cartesian_materialization(report, materialization_inputs)
    return _pqs_source_box_route_driver_materialization(
        report;
        materialization_inputs...,
    )
end

function cartesian_print_summary(report, materialization)
    recipe = report.recipe_metadata
    setup = report.standard_setup
    readiness = report.parent_axis_readiness
    route_axis_counts = report.route_axis_counts
    diagnostics = report.diagnostics
    retained_counts = report.retained_counts
    retained_dimension = report.retained_dimension

    route_family = report.route_family
    route_kind = recipe.route_kind
    q = recipe.q
    radius = report.system_metadata.radius

    println("Cartesian nesting route driver")
    @show route_family route_kind q radius
    if route_family == :pqs_source_box
        route_shape = recipe.route_shape
        product_body_rule = recipe.product_body_rule
        pair_factor_normalization = recipe.pair_factor_normalization
        @show route_shape product_body_rule pair_factor_normalization
    else
        white_lindsey_route_shape = recipe.route_shape
        white_lindsey_mapping_rule = recipe.white_lindsey_mapping_rule
        white_lindsey_nesting_rule = recipe.white_lindsey_nesting_rule
        @show white_lindsey_route_shape white_lindsey_mapping_rule white_lindsey_nesting_rule
    end

    @show setup.n_s setup.core_cube_side setup.core_spacing
    @show setup.spacing.q_to_core_spacing_rule_status
    @show readiness.status readiness.parent_axis_counts_status
    @show route_axis_counts.parent_axis_counts_source route_axis_counts.parent_axis_counts
    @show diagnostics.parent_axis_probe_requested diagnostics.parent_axis_probe_status
    @show diagnostics.raw_product_box_probe_requested diagnostics.raw_product_box_probe_status
    @show retained_counts retained_dimension
    @show get(report.low_order_route_summary, :materialization_status, nothing)
    @show materialization.status materialization.pqs_materialization_status
    @show materialization.materialized_report_kind
    return nothing
end

function cartesian_print_details(report, materialization)
    _pqs_source_box_route_driver_print_materialization(materialization)
    return nothing
end

function cartesian_save(report, save_inputs, materialization)
    return _pqs_source_box_route_driver_save(
        report;
        save_inputs...,
        materialization,
    )
end

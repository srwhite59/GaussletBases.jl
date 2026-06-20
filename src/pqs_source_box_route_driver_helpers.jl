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

function _pqs_source_box_route_driver_positive_integer(value, name::Symbol)
    if value isa Integer
        result = Int(value)
    elseif value isa Real && isinteger(value)
        result = Int(value)
    else
        throw(ArgumentError("$(name) must be a positive integer"))
    end
    result > 0 || throw(ArgumentError("$(name) must be positive"))
    return result
end

function _pqs_source_box_route_driver_odd_q_core_side(q)
    q_int = _pqs_source_box_route_driver_positive_integer(q, :q)
    core_cube_side = isodd(q_int) ? q_int : q_int + 1
    core_cube_side > 0 && isodd(core_cube_side) && core_cube_side >= q_int ||
        throw(ArgumentError("PQS core side must be positive, odd, and at least q"))
    return core_cube_side
end

function _pqs_source_box_route_driver_standard_setup(system_inputs, spacing_inputs)
    q = _pqs_source_box_route_driver_positive_integer(spacing_inputs.q, :q)
    n_s = _pqs_source_box_route_driver_positive_integer(spacing_inputs.n_s, :n_s)
    core_cube_side = _pqs_source_box_route_driver_odd_q_core_side(q)
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
        core_cube_side_rule = :odd_q_core_side,
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
            supplement_basis = get(route_inputs, :supplement_basis, nothing),
            supplement_lmax = get(route_inputs, :supplement_lmax, nothing),
            supplement_uncontracted =
                get(route_inputs, :supplement_uncontracted, nothing),
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
    parent_basis_object =
        hasproperty(parent, :parent_basis_object) ? parent.parent_basis_object : nothing
    if !isnothing(parent_basis_object)
        axes = CartesianParentGaussletBases.parent_axes(parent_basis_object)
        return (
            collect(axes.x.center_data),
            collect(axes.y.center_data),
            collect(axes.z.center_data),
        )
    end
    counts = _pqs_route_driver_axis_count_tuple(parent.axis_counts)
    isnothing(counts) &&
        throw(ArgumentError("terminal shellification requires parent axis counts"))
    all(>(0), counts) ||
        throw(ArgumentError("terminal shellification requires positive parent axis counts"))
    return (collect(1:counts[1]), collect(1:counts[2]), collect(1:counts[3]))
end

function _pqs_source_box_route_driver_terminal_center_index(count::Int)
    count > 0 || throw(ArgumentError("terminal shellification axis count must be positive"))
    return cld(count, 2)
end

function _pqs_source_box_route_driver_terminal_nuclear_positions(parent)
    parent_basis_object =
        hasproperty(parent, :parent_basis_object) ? parent.parent_basis_object : nothing
    if !isnothing(parent_basis_object)
        return Tuple(
            ntuple(axis -> Float64(location[axis]), 3)
            for location in parent.atom_locations
        )
    end

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

function _pqs_source_box_route_driver_shell_stage_terminal_shellification(
    parent,
    recipe,
)
    parent_axes = _pqs_source_box_route_driver_terminal_parent_axes(parent)
    nuclear_positions =
        _pqs_source_box_route_driver_terminal_nuclear_positions(parent)
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
        route_family = recipe.route_family,
    )
    return (;
        status = :available_terminal_cartesian_shellification_geometry,
        blocker = nothing,
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

function _pqs_source_box_route_driver_blocked_terminal_shellification(
    parent,
    recipe,
    error,
)
    return (;
        status = :blocked_terminal_cartesian_shellification_geometry,
        blocker = Symbol(nameof(typeof(error))),
        blocker_message = sprint(showerror, error),
        policy = nothing,
        plan = nothing,
        raw_plan = nothing,
        scaffold = nothing,
        route_family = recipe.route_family,
        system_classification = parent.system_classification,
        bond_axis = parent.bond_axis,
        core_side = parent.standard_setup.core_cube_side,
        q = parent.standard_setup.q,
        parent_axis_counts = parent.axis_counts,
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
            parent.system_classification in (:one_center, :bond_aligned_diatomic)
        )
    uses_terminal_shellification ||
        return nothing

    terminal = try
        _pqs_source_box_route_driver_shell_stage_terminal_shellification(
            parent,
            recipe,
        )
    catch error
        (error isa ArgumentError || error isa DimensionMismatch) || rethrow()
        _pqs_source_box_route_driver_blocked_terminal_shellification(
            parent,
            recipe,
            error,
        )
    end
    return (;
        status = terminal.status,
        blocker = terminal.blocker,
        blocker_message = get(terminal, :blocker_message, nothing),
        route_family = recipe.route_family,
        shellification_kind = :terminal_cartesian_shellification_geometry,
        terminal_shellification = terminal,
        shellification_plan = get(terminal, :plan, nothing),
        raw_shellification_plan = get(terminal, :raw_plan, nothing),
        shellification_scaffold = get(terminal, :scaffold, nothing),
        region_count = get(terminal, :region_count, nothing),
        ordered_region_roles = get(terminal, :ordered_region_roles, ()),
        spatial_policy_order = get(terminal, :spatial_policy_order, nothing),
        coverage_complete = get(terminal, :coverage_complete, false),
        central_gap_region_count = get(terminal, :central_gap_region_count, 0),
        central_midpoint_slab_count =
            get(terminal, :central_midpoint_slab_count, 0),
        central_distorted_product_box_count =
            get(terminal, :central_distorted_product_box_count, 0),
        central_distorted_product_box_metadata =
            get(terminal, :central_distorted_product_box_metadata, ()),
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
        low_order_shellification = shellification,
        shellification_status =
            isnothing(shellification) ?
            :not_requested :
            shellification.status,
        shellification_blocker =
            isnothing(shellification) ?
            nothing :
            shellification.blocker,
        shellification_blocker_message =
            isnothing(shellification) ?
            nothing :
            get(shellification, :blocker_message, nothing),
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
        all_units_have_lowering_contracts =
            all(entry -> entry.lowering_contract_count >= 1, contract_counts_by_unit),
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
    retained_unit_plan =
        lowering_plan isa CartesianTerminalLowering.TerminalLoweringPlan ?
        CartesianRetainedUnits.retained_unit_plan(lowering_plan) :
        nothing

    return (;
        shellification_kind,
        shellification_plan,
        shellification_scaffold,
        unit_inventory,
        route_lowering_family,
        lowering_plan,
        lowering_contract_inventory,
        retained_unit_plan,
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
        low_order_shellification = get(shells, :low_order_shellification, nothing),
        shellification_status =
            get(shells, :shellification_status, :not_requested),
        shellification_blocker =
            get(shells, :shellification_blocker, nothing),
        shellification_blocker_message =
            get(shells, :shellification_blocker_message, nothing),
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
    retained_unit_plan = low_order_units.retained_unit_plan
    retained_unit_transform_contract_plan =
        retained_unit_plan isa CartesianRetainedUnits.RetainedUnitPlan ?
        CartesianRetainedUnitTransformContracts.retained_unit_transform_contract_plan(
            retained_unit_plan,
        ) :
        nothing
    terminal_retained_rule_plan =
        _pqs_source_box_route_driver_terminal_retained_rule_plan(
            units.parent,
            low_order_units,
            retained_unit_plan,
            retained_unit_transform_contract_plan,
        )
    return (;
        shellification_kind = low_order_units.shellification_kind,
        shellification_plan = low_order_units.shellification_plan,
        shellification_scaffold = low_order_units.shellification_scaffold,
        unit_inventory = low_order_units.unit_inventory,
        route_lowering_family = low_order_units.route_lowering_family,
        lowering_plan = low_order_units.lowering_plan,
        lowering_contract_inventory = low_order_units.lowering_contract_inventory,
        retained_unit_plan,
        retained_unit_transform_contract_plan,
        terminal_retained_rule_plan,
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
        low_order_shellification =
            get(units, :low_order_shellification, nothing),
        shellification_status =
            get(units, :shellification_status, :not_requested),
        shellification_blocker =
            get(units, :shellification_blocker, nothing),
        shellification_blocker_message =
            get(units, :shellification_blocker_message, nothing),
        low_order_transforms,
        terminal_retained_rule_plan =
            isnothing(low_order_transforms) ?
            nothing :
            low_order_transforms.terminal_retained_rule_plan,
        retained_unit_plan =
            isnothing(low_order_transforms) ?
            nothing :
            low_order_transforms.retained_unit_plan,
        retained_unit_transform_contract_plan =
            isnothing(low_order_transforms) ?
            nothing :
            low_order_transforms.retained_unit_transform_contract_plan,
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
        retained_unit_plan = low_order_transforms.retained_unit_plan,
        retained_unit_transform_contract_plan =
            low_order_transforms.retained_unit_transform_contract_plan,
        terminal_retained_rule_plan =
            low_order_transforms.terminal_retained_rule_plan,
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
        retained_unit_plan = low_order_pairs.retained_unit_plan,
        retained_unit_transform_contract_plan =
            low_order_pairs.retained_unit_transform_contract_plan,
        terminal_retained_rule_plan =
            low_order_pairs.terminal_retained_rule_plan,
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

function cartesian_assembly(parent, shells, units, transforms, pairs, recipe)
    route_skeleton = shells.route_skeleton
    route_facts = _pqs_source_box_route_driver_route_facts(route_skeleton)
    low_order_assembly =
        _pqs_source_box_route_driver_assembly_stage_low_order_summary(pairs)
    diatomic_physical_gausslet_target_payload =
        _pqs_source_box_route_driver_diatomic_physical_gausslet_target_payload(
            parent,
            route_skeleton,
            recipe,
            low_order_assembly,
        )
    diatomic_physical_gausslet_supplement_request_payload =
        _pqs_source_box_route_driver_diatomic_physical_gausslet_supplement_request_payload(
            parent,
            diatomic_physical_gausslet_target_payload,
            recipe,
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
            low_order_assembly,
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
    return (;
        spacing_inputs = shells.spacing_inputs,
        route_skeleton,
        route_facts,
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
    )
end

_pqs_source_box_route_driver_payload_summary(payload) =
    isnothing(payload) ? nothing :
    hasproperty(payload, :summary) ? payload.summary :
    payload

function _pqs_source_box_route_driver_pqs_gto_sidecar_inputs(assembly)
    source_plan_payload =
        hasproperty(assembly, :diatomic_physical_gausslet_source_plan_payload) ?
        assembly.diatomic_physical_gausslet_source_plan_payload :
        nothing
    final_basis_payload =
        hasproperty(assembly, :diatomic_physical_gausslet_final_basis_payload) ?
        assembly.diatomic_physical_gausslet_final_basis_payload :
        nothing
    h1_payload =
        hasproperty(assembly, :diatomic_physical_gausslet_h1_payload) ?
        assembly.diatomic_physical_gausslet_h1_payload :
        nothing
    h1_j_payload =
        hasproperty(assembly, :diatomic_physical_gausslet_h1_j_payload) ?
        assembly.diatomic_physical_gausslet_h1_j_payload :
        nothing
    supplement_payload =
        hasproperty(assembly, :diatomic_physical_gausslet_supplement_representation_payload) ?
        assembly.diatomic_physical_gausslet_supplement_representation_payload :
        nothing

    source_plan =
        !isnothing(source_plan_payload) &&
        hasproperty(source_plan_payload, :source_plan) ?
        source_plan_payload.source_plan :
        nothing
    final_basis =
        !isnothing(final_basis_payload) &&
        hasproperty(final_basis_payload, :final_basis) ?
        final_basis_payload.final_basis :
        nothing
    h1 =
        !isnothing(h1_payload) &&
        hasproperty(h1_payload, :h1) ?
        h1_payload.h1 :
        nothing
    h1_hamiltonian =
        !isnothing(h1_payload) &&
        hasproperty(h1_payload, :final_hamiltonian) ?
        h1_payload.final_hamiltonian :
        nothing
    final_kinetic =
        !isnothing(h1_payload) &&
        hasproperty(h1_payload, :final_kinetic) ?
        h1_payload.final_kinetic :
        nothing
    final_nuclear_by_center =
        !isnothing(h1_payload) &&
        hasproperty(h1_payload, :final_electron_nuclear_by_center) ?
        h1_payload.final_electron_nuclear_by_center :
        nothing
    density_interaction =
        !isnothing(h1_j_payload) &&
        hasproperty(h1_j_payload, :density_interaction) ?
        h1_j_payload.density_interaction :
        nothing
    h1_j =
        !isnothing(h1_j_payload) &&
        hasproperty(h1_j_payload, :h1_j_diagnostic) ?
        h1_j_payload.h1_j_diagnostic :
        nothing
    supplement_representation =
        !isnothing(supplement_payload) &&
        hasproperty(supplement_payload, :representation) ?
        supplement_payload.representation :
        nothing

    any(
        isnothing,
        (
            source_plan,
            final_basis,
            h1,
            h1_hamiltonian,
            final_kinetic,
            final_nuclear_by_center,
            density_interaction,
            h1_j,
            supplement_representation,
        ),
    ) && return nothing

    return (;
        source_plan,
        final_basis,
        h1,
        h1_hamiltonian,
        final_kinetic,
        final_nuclear_by_center,
        density_interaction,
        h1_j,
        supplement_representation,
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
    route_skeleton = assembly.route_skeleton
    route_facts = assembly.route_facts

    source_recipe =
        recipe.route_family == :pqs_source_box ? recipe.source_box : recipe.white_lindsey
    private_rhf_inputs = get(recipe, :private_rhf_inputs, (;))
    system_metadata = (;
        atom_symbols = parent.atom_symbols,
        nuclear_charges = parent.nuclear_charges,
        atom_locations = parent.atom_locations,
        nup = get(system, :nup, nothing),
        ndn = get(system, :ndn, nothing),
        bond_axis = parent.bond_axis,
        bond_length = get(system, :bond_length, nothing),
        radius = get(system, :radius, nothing),
        map_backend = get(system, :map_backend, nothing),
    )
    recipe_metadata = (;
        route_family = recipe.route_family,
        route_kind = recipe.route_kind,
        route_shape = get(source_recipe, :route_shape, nothing),
        q = assembly.spacing_inputs.q,
        n_s = standard_setup.n_s,
        core_cube_side = standard_setup.core_cube_side,
        reference_spacing = standard_setup.reference_spacing,
        tail_spacing = standard_setup.tail_spacing,
        q_to_core_spacing_rule = standard_setup.q_to_core_spacing_rule,
        core_spacing = standard_setup.core_spacing,
        xmax_parallel = get(assembly.spacing_inputs, :xmax_parallel, nothing),
        xmax_transverse = get(assembly.spacing_inputs, :xmax_transverse, nothing),
        terms = recipe.terms,
        pair_factor_normalization = recipe.pair_factor_normalization,
        supplement_policy = get(recipe, :supplement_policy, nothing),
        run_final_basis = get(recipe, :run_final_basis, false),
        run_h1 = get(recipe, :run_h1, false),
        run_h1_j = get(recipe, :run_h1_j, false),
        run_private_rhf = get(private_rhf_inputs, :run_private_rhf, false),
        product_body_rule = get(source_recipe, :product_body_rule, nothing),
        pqs_retained_rule = get(source_recipe, :pqs_retained_rule, nothing),
        product_retained_rule = get(source_recipe, :product_retained_rule, nothing),
        white_lindsey_mapping_rule = get(source_recipe, :mapping_rule, nothing),
        white_lindsey_nesting_rule = get(source_recipe, :nesting_rule, nothing),
        white_lindsey_retained_rule = get(source_recipe, :retained_rule, nothing),
        white_lindsey_operator_rule = get(source_recipe, :operator_rule, nothing),
    )
    parent_summary = (;
        atom_count = parent.atom_count,
        system_classification = parent.system_classification,
        center_table = parent.center_table,
        axis_counts = parent.axis_counts,
        physical_box = parent.physical_box,
    )
    route_summary = (;
        route_shape = get(route_skeleton, :route_shape, nothing),
        retained_unit_order = get(route_skeleton, :retained_unit_order, ()),
        source_boxes = route_facts.source_boxes,
        source_dimensions = route_facts.source_dimensions,
        retained_counts = route_facts.retained_counts,
        retained_dimension = route_facts.retained_dimension,
        shellification_kind = get(assembly.shells, :shellification_kind, nothing),
    )
    pair_summary = (;
        pair_entries = assembly.pairs.pair_entries,
        pair_family_counts = assembly.pairs.pair_family_counts,
        helper_by_pair_family = assembly.pairs.helper_by_pair_family,
    )

    return merge(
        (;
            generated_at = Base.Libc.strftime("%Y-%m-%dT%H:%M:%S", time()),
            route_family = route_skeleton.route_family,
            route_kind = get(route_skeleton, :route_kind, recipe.route_kind),
            standard_setup,
            system_metadata,
            recipe_metadata,
            parent_summary,
            route_summary,
            pair_summary,
            source_boxes = route_facts.source_boxes,
            source_dimensions = route_facts.source_dimensions,
            retained_units = route_facts.retained_units,
            retained_counts = route_facts.retained_counts,
            ranges = route_facts.ranges,
            retained_dimension = route_facts.retained_dimension,
            pair_entries = assembly.pairs.pair_entries,
            pair_family_counts = assembly.pairs.pair_family_counts,
            helper_by_pair_family = assembly.pairs.helper_by_pair_family,
            parent_basis_object = parent.parent_basis_object,
            parent_axis_bundle_object = parent.parent_axis_bundle_object,
            axis_bundle_backend = parent.parent_inputs.parent_axis_bundle_backend,
            pqs_gto_sidecar_inputs =
                _pqs_source_box_route_driver_pqs_gto_sidecar_inputs(assembly),
        ),
        _pqs_source_box_route_driver_physical_gausslet_target_report_fields(
            assembly,
        ),
    )
end

function cartesian_materialization(report, materialization_inputs)
    return _pqs_source_box_route_driver_materialization(
        report;
        materialization_inputs...,
    )
end

function cartesian_print_summary(report, materialization)
    maybe_get(obj, key) = isnothing(obj) ? nothing : get(obj, key, nothing)
    recipe = report.recipe_metadata
    setup = report.standard_setup
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
    @show report.parent_summary.axis_counts
    @show report.parent_summary.system_classification
    @show retained_counts retained_dimension
    @show report.pair_family_counts
    @show report.route_summary.shellification_kind
    final_basis = report.physical_gausslet_final_basis_summary
    h1 = report.physical_gausslet_h1_summary
    h1_j = report.physical_gausslet_h1_j_summary
    rhf = report.physical_gausslet_private_rhf_summary
    supplement = report.physical_gausslet_supplement_preflight_summary
    @show maybe_get(final_basis, :final_dimension)
    @show maybe_get(final_basis, :final_overlap_identity_error)
    @show maybe_get(h1, :lowest_energy)
    @show maybe_get(h1, :h1_hamiltonian_symmetry_error)
    @show maybe_get(h1_j, :self_coulomb)
    @show maybe_get(h1_j, :electron_electron_ida_symmetry_error)
    @show maybe_get(rhf, :total_with_nuclear_repulsion)
    @show maybe_get(rhf, :commutator_residual)
    @show maybe_get(supplement, :supplement_policy)
    @show maybe_get(supplement, :provider_block_count)
    @show materialization.result_kind materialization.materialized
    @show get(materialization, :h1_lowest, nothing)
    @show get(materialization, :overlap_identity_error, nothing)
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

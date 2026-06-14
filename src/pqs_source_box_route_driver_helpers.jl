# Private support for `bin/cartesian_ham_builder.jl`.
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

function _cartesian_parent_symbol_tuple(atom_symbols)
    atom_symbols isa AbstractString && return (atom_symbols,)
    atom_symbols isa Symbol && return (atom_symbols,)
    return Tuple(atom_symbols)
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
            location = _cartesian_shellization_route_location_tuple(atom_locations[index]),
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
        object_kind = :cartesian_parent_axis_metadata,
        status =
            axis_aligned ?
            :axis_aligned_or_one_center :
            :not_axis_aligned,
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
            system_classification_status = :explicit_atom_count_one,
            bond_axis = nothing,
            chain_axis = nothing,
        )
    elseif atom_count == 2 && !isnothing(axis_metadata.bond_axis)
        return (
            system_classification = :bond_aligned_diatomic,
            system_classification_status = :explicit_two_atom_single_axis_separation,
            bond_axis = axis_metadata.bond_axis,
            chain_axis = axis_metadata.chain_axis,
        )
    elseif atom_count == 2
        return (
            system_classification = :pending_system_classification,
            system_classification_status = :diatomic_not_axis_aligned_by_metadata,
            bond_axis = nothing,
            chain_axis = nothing,
        )
    elseif !isnothing(axis_metadata.chain_axis)
        return (
            system_classification = :axis_aligned_chain_metadata_only,
            system_classification_status = :multi_center_single_axis_chain_not_materialized,
            bond_axis = nothing,
            chain_axis = axis_metadata.chain_axis,
        )
    end
    return (
        system_classification = :unsupported_general_multi_atom,
        system_classification_status = :general_multi_atom_parent_materializer_not_planned,
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
    classification.system_classification == :one_center || return (;
        parent_basis_object = nothing,
        parent_basis_metadata_available = false,
        parent_basis_object_source = :not_applicable,
        status = :not_applicable_to_system_classification,
        pending_facts = (),
    )

    counts = route_axis_counts.parent_axis_counts
    isnothing(counts) && return (;
        parent_basis_object = nothing,
        parent_basis_metadata_available = false,
        parent_basis_object_source = :unavailable_parent_axis_counts,
        status = :pending_parent_axis_counts,
        pending_facts = (:parent_axis_counts,),
    )

    center = first(center_table)
    origin_centered = all(abs.(center.location) .<= 1.0e-12)
    origin_centered || return (;
        parent_basis_object = nothing,
        parent_basis_metadata_available = false,
        parent_basis_object_source = :pending_non_origin_one_center_parent_mapping,
        status = :pending_origin_centered_or_translated_one_center_parent_mapping,
        pending_facts = (:origin_centered_one_center_parent_or_translated_mapping,),
    )

    family = get(parent_inputs, :parent_axis_family,
        get(parent_inputs, :parent_axis_probe_family, :G10))
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

    parent_basis_object =
        CartesianParentGaussletBases.CartesianParentGaussletBasis3D(
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

    return (;
        parent_basis_object,
        parent_basis_metadata_available = true,
        parent_basis_object_source = :one_center_mapped_uniform_axes,
        status = :constructed_one_center_parent_basis_object,
        pending_facts = (:one_center_axis_bundle_builder,),
    )
end

function _cartesian_one_center_parent_axis_bundle_object(
    classification,
    parent_basis_object,
    one_center_parent_basis,
    parent_inputs,
)
    classification.system_classification == :one_center || return (;
        axis_bundle_object = nothing,
        axis_bundle_metadata_available = false,
        axis_bundle_object_source = :not_applicable,
        status = :not_applicable_to_system_classification,
        pending_facts = (),
    )
    isnothing(parent_basis_object) && return (;
        axis_bundle_object = nothing,
        axis_bundle_metadata_available = false,
        axis_bundle_object_source =
            isnothing(one_center_parent_basis) ?
            :unavailable_parent_basis_object :
            one_center_parent_basis.parent_basis_object_source,
        status = :pending_one_center_parent_basis_object,
        pending_facts = (:one_center_parent_basis_object,),
    )
    isnothing(one_center_parent_basis) && return (;
        axis_bundle_object = nothing,
        axis_bundle_metadata_available = false,
        axis_bundle_object_source = :not_applicable_to_non_one_center_parent_basis,
        status = :not_applicable_to_non_one_center_parent_basis,
        pending_facts = (),
    )
    one_center_parent_basis.parent_basis_object_source ==
        :one_center_mapped_uniform_axes || return (;
            axis_bundle_object = nothing,
            axis_bundle_metadata_available = false,
            axis_bundle_object_source =
                one_center_parent_basis.parent_basis_object_source,
            status = :pending_reviewed_one_center_axis_bundle_contract,
            pending_facts = one_center_parent_basis.pending_facts,
        )

    axes = CartesianParentGaussletBases.parent_axes(parent_basis_object)
    expansion = coulomb_gaussian_expansion(doacc = false)
    backend =
        hasproperty(parent_inputs, :parent_axis_probe_backend) ?
        parent_inputs.parent_axis_probe_backend :
        :pgdg_localized_experimental
    function _axis_bundle(axis)
        return _mapped_ordinary_gausslet_1d_bundle(
            axis;
            exponents = expansion.exponents,
            center = 0.0,
            backend,
        )
    end
    axis_bundle_object =
        _CartesianNestedAxisBundles3D(
            _axis_bundle(axes.x),
            _axis_bundle(axes.y),
            _axis_bundle(axes.z),
        )

    return (;
        axis_bundle_object,
        axis_bundle_metadata_available = true,
        axis_bundle_object_source = :one_center_mapped_ordinary_axis_bundles,
        status = :constructed_one_center_parent_axis_bundle_object,
        pending_facts = (),
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
    parent_axis_probe = parent_axis.parent_axis_probe
    qw_basis_object =
        !isnothing(parent_axis_probe) &&
        hasproperty(parent_axis_probe, :basis_object) ?
        parent_axis_probe.basis_object :
        nothing
    axis_bundle_object =
        !isnothing(parent_axis_probe) &&
        hasproperty(parent_axis_probe, :axis_bundle_object) ?
        parent_axis_probe.axis_bundle_object :
        nothing
    parent_axis_probe_bundle_object_available = !isnothing(axis_bundle_object)
    parent_basis_object =
        isnothing(qw_basis_object) ?
        nothing :
        CartesianParentGaussletBases.cartesian_parent_gausslet_basis(qw_basis_object)
    one_center_parent_basis =
        isnothing(parent_basis_object) ?
        _cartesian_one_center_parent_basis_object(
            center_table, classification, standard_setup,
            route_axis_counts, parent_inputs) :
        nothing
    parent_basis_object =
        isnothing(parent_basis_object) && !isnothing(one_center_parent_basis) ?
        one_center_parent_basis.parent_basis_object :
        parent_basis_object
    parent_basis_metadata_available =
        !isnothing(qw_basis_object) ||
        (
            !isnothing(one_center_parent_basis) &&
            one_center_parent_basis.parent_basis_metadata_available
        )
    parent_basis_object_source =
        !isnothing(qw_basis_object) ?
        :parent_axis_probe_qw_basis :
        isnothing(one_center_parent_basis) ?
        :unavailable :
        one_center_parent_basis.parent_basis_object_source
    one_center_axis_bundle =
        isnothing(axis_bundle_object) ?
        _cartesian_one_center_parent_axis_bundle_object(
            classification, parent_basis_object, one_center_parent_basis,
            parent_inputs) :
        nothing
    axis_bundle_object =
        isnothing(axis_bundle_object) && !isnothing(one_center_axis_bundle) ?
        one_center_axis_bundle.axis_bundle_object :
        axis_bundle_object
    axis_bundle_metadata_available =
        (
            !isnothing(parent_axis_probe) &&
            parent_axis_probe.axis_bundle_metadata.status == :constructed
        ) ||
        (
            !isnothing(one_center_axis_bundle) &&
            one_center_axis_bundle.axis_bundle_metadata_available
        )
    axis_bundle_object_source =
        parent_axis_probe_bundle_object_available ?
        :parent_axis_probe_axis_bundle :
        isnothing(one_center_axis_bundle) ?
        :unavailable :
        one_center_axis_bundle.axis_bundle_object_source

    return (;
        object_kind = :cartesian_parent_object_carry,
        status =
            !isnothing(parent_basis_object) && !isnothing(axis_bundle_object) ?
            :materialized_parent_objects_available :
            !isnothing(parent_basis_object) ?
            :parent_basis_object_available_axis_bundle_pending :
            :not_materialized,
        parent_basis_object,
        parent_qw_basis_object = qw_basis_object,
        parent_axis_bundle_object = axis_bundle_object,
        parent_basis_object_available = !isnothing(parent_basis_object),
        parent_qw_basis_object_available = !isnothing(qw_basis_object),
        parent_axis_bundle_object_available = !isnothing(axis_bundle_object),
        parent_basis_metadata_available,
        parent_axis_bundle_metadata_available = axis_bundle_metadata_available,
        parent_basis_object_source,
        parent_axis_bundle_object_source = axis_bundle_object_source,
        parent_basis_object_type_label = _pqs_route_driver_type_label(parent_basis_object),
        parent_qw_basis_object_type_label = _pqs_route_driver_type_label(qw_basis_object),
        parent_axis_bundle_object_type_label =
            _pqs_route_driver_type_label(axis_bundle_object),
        parent_axis_counts =
            isnothing(parent_basis_object) ?
            nothing :
            CartesianParentGaussletBases.parent_axis_counts(parent_basis_object),
        object_carry_scope = :transient_cartesian_parent_only,
        report_must_sanitize = true,
    )
end

function _cartesian_parent_materialization_status(parent_axis, object_carry)
    parent_axis_probe = parent_axis.parent_axis_probe
    parent_axis_metadata_constructed = parent_axis.parent_axis_probe_constructed
    axis_bundle_status =
        object_carry.parent_axis_bundle_object_available ?
        :constructed :
        !isnothing(parent_axis_probe) ?
        parent_axis_probe.axis_bundle_metadata.status :
        object_carry.parent_basis_object_available ?
        :pending_one_center_axis_bundle_builder :
        :not_requested
    parent_basis_object_available = object_carry.parent_basis_object_available
    axis_bundle_object_available = object_carry.parent_axis_bundle_object_available

    return (
        object_kind = :cartesian_parent_materialization_status,
        status =
            parent_basis_object_available && axis_bundle_object_available ?
            :materialized_parent_objects_available :
            parent_basis_object_available ?
            :parent_basis_object_available_axis_bundle_pending :
            parent_axis_metadata_constructed ?
            :metadata_constructed_probe_only :
            :metadata_only_not_materialized,
        parent_axis_metadata_constructed,
        parent_basis_materialized = parent_basis_object_available,
        parent_basis_object_available,
        parent_basis_object_type_label = object_carry.parent_basis_object_type_label,
        parent_basis_metadata_available =
            object_carry.parent_basis_metadata_available,
        axis_bundle_materialized = axis_bundle_object_available,
        axis_bundle_object_available,
        axis_bundle_object_type_label =
            object_carry.parent_axis_bundle_object_type_label,
        axis_bundle_metadata_status = axis_bundle_status,
        axis_bundle_metadata_available =
            object_carry.parent_axis_bundle_metadata_available,
        materialization_scope =
            parent_basis_object_available && axis_bundle_object_available ?
            :transient_parent_objects_carried_report_metadata_only :
            parent_basis_object_available ?
            :transient_parent_basis_carried_axis_bundle_pending_report_metadata_only :
            :metadata_only_parent_contract,
    )
end

function _cartesian_parent_materialization_plan_available_inputs(
    standard_setup,
    parent_axis,
    route_axis_counts,
)
    inputs = Symbol[
        :center_table,
        :center_axis_metadata,
        :system_classification,
        :standard_setup,
        :parent_axis_readiness,
        :route_axis_counts,
    ]
    isnothing(standard_setup.core_spacing) || push!(inputs, :core_spacing)
    isnothing(route_axis_counts.parent_axis_counts) ||
        push!(inputs, :parent_axis_counts)
    parent_axis.parent_axis_probe_constructed &&
        push!(inputs, :parent_axis_probe_metadata)
    return Tuple(inputs)
end

function _cartesian_parent_materialization_plan(
    center_table,
    center_axis_metadata,
    classification,
    standard_setup,
    parent_axis,
    route_axis_counts,
    object_carry,
)
    atom_count = length(center_table)
    basis_metadata_available =
        !isnothing(parent_axis.parent_axis_probe) &&
        !isnothing(parent_axis.parent_axis_probe.basis_metadata) ||
        object_carry.parent_basis_metadata_available
    axis_bundle_metadata_available =
        !isnothing(parent_axis.parent_axis_probe) &&
        parent_axis.parent_axis_probe.axis_bundle_metadata.status == :constructed
    basis_object_available = object_carry.parent_basis_object_available
    axis_bundle_object_available = object_carry.parent_axis_bundle_object_available
    available_inputs =
        _cartesian_parent_materialization_plan_available_inputs(
            standard_setup, parent_axis, route_axis_counts)
    route_axis_pending = route_axis_counts.pending_facts

    family = classification.system_classification
    if family == :one_center
        carried_objects_available =
            basis_object_available && axis_bundle_object_available
        missing_inputs =
            carried_objects_available ?
            () :
            basis_object_available ?
            (:one_center_axis_bundle_builder,) :
            (
                route_axis_pending...,
                :reviewed_one_center_parent_axis_builder,
                :one_center_parent_basis_or_axis_bundle_metadata_probe,
            )
        return (;
            object_kind = :cartesian_parent_materialization_plan,
            status =
                carried_objects_available ?
                :materialized_parent_objects_available :
                basis_object_available ?
                :one_center_parent_basis_carried_axis_bundle_pending :
                :metadata_only_pending_one_center_parent_axis_builder,
            planning_family = :one_center_parent_lattice,
            system_classification = family,
            system_classification_status =
                classification.system_classification_status,
            atom_count,
            center_count = atom_count,
            bond_axis = nothing,
            chain_axis = nothing,
            intended_parent_constructor = :CartesianParentGaussletBasis3D,
            intended_axis_basis_helper = :pending_one_center_mapped_axis_builder,
            intended_axis_bundle_helper =
                carried_objects_available ?
                :_cartesian_one_center_parent_axis_bundle_object :
                :pending_one_center_axis_bundle_builder,
            route_neutral_parent_wrapper = :cartesian_parent_gausslet_basis,
            one_center_compatible = true,
            bond_aligned_diatomic_compatible = false,
            axis_aligned_chain_compatible = false,
            metadata_only = !carried_objects_available,
            constructs_basis_now = basis_object_available,
            constructs_axis_bundle_now = axis_bundle_object_available,
            basis_metadata_available,
            axis_bundle_metadata_available,
            basis_object_available,
            axis_bundle_object_available,
            required_inputs_available = available_inputs,
            available_input_count = length(available_inputs),
            missing_inputs,
            missing_input_count = length(missing_inputs),
            pending_facts = missing_inputs,
            blocked = !carried_objects_available,
            blocker =
                carried_objects_available ?
                nothing :
                basis_object_available ?
                :pending_one_center_axis_bundle_builder :
                :pending_one_center_parent_axis_builder,
            diagnostics = (
                source = :cartesian_parent_materialization_plan,
                private_development_only = true,
                materialization_scope =
                    carried_objects_available ?
                    :transient_parent_objects_available :
                    basis_object_available ?
                    :transient_one_center_parent_basis_object_available :
                    :metadata_only_parent_planning,
                public_default_behavior_changed = false,
                shellification_restructured = false,
                hamiltonian_export = false,
            ),
        )
    elseif family == :bond_aligned_diatomic
        homonuclear = parent_axis.parent_axis_readiness.homonuclear
        constructor =
            homonuclear ?
            :bond_aligned_homonuclear_qw_basis :
            :bond_aligned_heteronuclear_qw_basis
        probe_constructed = parent_axis.parent_axis_probe_constructed
        carried_objects_available =
            basis_object_available && axis_bundle_object_available
        missing_inputs =
            carried_objects_available ?
            () :
            probe_constructed ?
            (:real_parent_basis_object_carry_or_reviewed_materializer_handoff,) :
            (
                route_axis_pending...,
                parent_axis.parent_axis_readiness.pending_facts...,
                :parent_axis_metadata_probe_or_real_parent_basis_object,
            )
        return (;
            object_kind = :cartesian_parent_materialization_plan,
            status =
                carried_objects_available ?
                :materialized_parent_objects_available :
                probe_constructed ?
                :metadata_probe_available_pending_parent_basis_carry :
                :metadata_only_diatomic_parent_api_candidate,
            planning_family = :bond_aligned_diatomic_parent_lattice,
            system_classification = family,
            system_classification_status =
                classification.system_classification_status,
            atom_count,
            center_count = atom_count,
            bond_axis = classification.bond_axis,
            chain_axis = classification.chain_axis,
            intended_parent_constructor = constructor,
            intended_axis_basis_helper = constructor,
            intended_axis_bundle_helper = :_qwrg_bond_aligned_axis_bundles,
            route_neutral_parent_wrapper = :cartesian_parent_gausslet_basis,
            one_center_compatible = false,
            bond_aligned_diatomic_compatible = true,
            axis_aligned_chain_compatible = false,
            metadata_only = !carried_objects_available,
            constructs_basis_now = carried_objects_available,
            constructs_axis_bundle_now = carried_objects_available,
            basis_metadata_available,
            axis_bundle_metadata_available,
            basis_object_available,
            axis_bundle_object_available,
            required_inputs_available = available_inputs,
            available_input_count = length(available_inputs),
            missing_inputs,
            missing_input_count = length(missing_inputs),
            pending_facts = missing_inputs,
            blocked = !carried_objects_available,
            blocker =
                carried_objects_available ?
                nothing :
                probe_constructed ?
                :pending_real_parent_basis_carry :
                :pending_reviewed_diatomic_parent_materializer,
            diagnostics = (
                source = :cartesian_parent_materialization_plan,
                private_development_only = true,
                materialization_scope =
                    carried_objects_available ?
                    :transient_parent_objects_available :
                    :metadata_only_parent_planning,
                public_default_behavior_changed = false,
                shellification_restructured = false,
                hamiltonian_export = false,
            ),
        )
    elseif family == :axis_aligned_chain_metadata_only
        homonuclear = parent_axis.parent_axis_readiness.homonuclear
        constructor =
            homonuclear ? :bond_aligned_homonuclear_chain_qw_basis : nothing
        missing_inputs = (
            route_axis_pending...,
            homonuclear ?
            :reviewed_chain_parent_axis_constructor_call :
            :heteronuclear_chain_parent_materializer_design,
            homonuclear ?
            :chain_parent_axis_bundle_metadata_probe :
            :heteronuclear_chain_axis_bundle_design,
        )
        return (;
            object_kind = :cartesian_parent_materialization_plan,
            status =
                homonuclear ?
                :metadata_only_chain_parent_constructor_candidate :
                :blocked_unsupported_heteronuclear_chain_parent_materializer,
            planning_family = :axis_aligned_chain_parent_lattice,
            system_classification = family,
            system_classification_status =
                classification.system_classification_status,
            atom_count,
            center_count = atom_count,
            bond_axis = nothing,
            chain_axis = classification.chain_axis,
            intended_parent_constructor = constructor,
            intended_axis_basis_helper = constructor,
            intended_axis_bundle_helper =
                homonuclear ?
                :pending_homonuclear_chain_axis_bundle_builder :
                nothing,
            route_neutral_parent_wrapper = :cartesian_parent_gausslet_basis,
            one_center_compatible = false,
            bond_aligned_diatomic_compatible = false,
            axis_aligned_chain_compatible = true,
            metadata_only = true,
            constructs_basis_now = false,
            constructs_axis_bundle_now = false,
            basis_metadata_available,
            axis_bundle_metadata_available,
            basis_object_available,
            axis_bundle_object_available,
            required_inputs_available = available_inputs,
            available_input_count = length(available_inputs),
            missing_inputs,
            missing_input_count = length(missing_inputs),
            pending_facts = missing_inputs,
            blocked = true,
            blocker =
                homonuclear ?
                :chain_parent_materializer_not_connected :
                :heteronuclear_chain_parent_materializer_not_available,
            diagnostics = (
                source = :cartesian_parent_materialization_plan,
                private_development_only = true,
                materialization_scope = :metadata_only_parent_planning,
                public_default_behavior_changed = false,
                shellification_restructured = false,
                hamiltonian_export = false,
            ),
        )
    elseif family == :pending_system_classification
        missing_inputs = (:axis_aligned_or_supported_parent_geometry_classification,)
        return (;
            object_kind = :cartesian_parent_materialization_plan,
            status = :blocked_pending_system_classification,
            planning_family = :pending_system_classification_parent_lattice,
            system_classification = family,
            system_classification_status =
                classification.system_classification_status,
            atom_count,
            center_count = atom_count,
            bond_axis = nothing,
            chain_axis = nothing,
            intended_parent_constructor = nothing,
            intended_axis_basis_helper = nothing,
            intended_axis_bundle_helper = nothing,
            route_neutral_parent_wrapper = :cartesian_parent_gausslet_basis,
            one_center_compatible = false,
            bond_aligned_diatomic_compatible = false,
            axis_aligned_chain_compatible = false,
            metadata_only = true,
            constructs_basis_now = false,
            constructs_axis_bundle_now = false,
            basis_metadata_available,
            axis_bundle_metadata_available,
            basis_object_available,
            axis_bundle_object_available,
            required_inputs_available = available_inputs,
            available_input_count = length(available_inputs),
            missing_inputs,
            missing_input_count = length(missing_inputs),
            pending_facts = missing_inputs,
            blocked = true,
            blocker = :pending_system_classification,
            diagnostics = (
                source = :cartesian_parent_materialization_plan,
                private_development_only = true,
                materialization_scope = :metadata_only_parent_planning,
                public_default_behavior_changed = false,
                shellification_restructured = false,
                hamiltonian_export = false,
            ),
        )
    end

    missing_inputs = (:general_multi_atom_parent_materializer_design,)
    return (;
        object_kind = :cartesian_parent_materialization_plan,
        status = :blocked_unsupported_parent_materializer,
        planning_family = :unsupported_parent_lattice,
        system_classification = family,
        system_classification_status = classification.system_classification_status,
        atom_count,
        center_count = atom_count,
        bond_axis = nothing,
        chain_axis = nothing,
        intended_parent_constructor = nothing,
        intended_axis_basis_helper = nothing,
        intended_axis_bundle_helper = nothing,
        route_neutral_parent_wrapper = :cartesian_parent_gausslet_basis,
        one_center_compatible = false,
        bond_aligned_diatomic_compatible = false,
        axis_aligned_chain_compatible = false,
        metadata_only = true,
        constructs_basis_now = false,
        constructs_axis_bundle_now = false,
        basis_metadata_available,
        axis_bundle_metadata_available,
        basis_object_available,
        axis_bundle_object_available,
        required_inputs_available = available_inputs,
        available_input_count = length(available_inputs),
        missing_inputs,
        missing_input_count = length(missing_inputs),
        pending_facts = missing_inputs,
        blocked = true,
        blocker = :unsupported_general_multi_atom_parent_materializer,
        diagnostics = (
            source = :cartesian_parent_materialization_plan,
            private_development_only = true,
            materialization_scope = :metadata_only_parent_planning,
            public_default_behavior_changed = false,
            shellification_restructured = false,
            hamiltonian_export = false,
        ),
    )
end

function _cartesian_parent_axis_probe_report_summary(parent_axis_probe)
    isnothing(parent_axis_probe) && return nothing
    carry_objects_requested =
        hasproperty(parent_axis_probe, :carry_objects_requested) &&
        parent_axis_probe.carry_objects_requested
    basis_object_available =
        hasproperty(parent_axis_probe, :basis_object_available) &&
        parent_axis_probe.basis_object_available
    axis_bundle_object_available =
        hasproperty(parent_axis_probe, :axis_bundle_object_available) &&
        parent_axis_probe.axis_bundle_object_available
    basis_object_type_label =
        hasproperty(parent_axis_probe, :basis_object_type_label) ?
        parent_axis_probe.basis_object_type_label :
        "unavailable"
    axis_bundle_object_type_label =
        hasproperty(parent_axis_probe, :axis_bundle_object_type_label) ?
        parent_axis_probe.axis_bundle_object_type_label :
        "unavailable"

    return (
        object_kind = parent_axis_probe.object_kind,
        status = parent_axis_probe.status,
        readiness = parent_axis_probe.readiness,
        basis_metadata = parent_axis_probe.basis_metadata,
        axis_bundle_metadata = parent_axis_probe.axis_bundle_metadata,
        axis_lengths = parent_axis_probe.axis_lengths,
        physical_extent_inputs = parent_axis_probe.physical_extent_inputs,
        core_spacing = parent_axis_probe.core_spacing,
        reference_spacing = parent_axis_probe.reference_spacing,
        tail_spacing = parent_axis_probe.tail_spacing,
        gausslet_backend = parent_axis_probe.gausslet_backend,
        gausslet_backend_role = parent_axis_probe.gausslet_backend_role,
        expansion_source = parent_axis_probe.expansion_source,
        explicit_spacing_probe_only = parent_axis_probe.explicit_spacing_probe_only,
        default_standard_rule = parent_axis_probe.default_standard_rule,
        core_spacing_source = parent_axis_probe.core_spacing_source,
        parent_axis_metadata_constructed =
            parent_axis_probe.parent_axis_metadata_constructed,
        carry_objects_requested,
        basis_object_available,
        axis_bundle_object_available,
        basis_object_type_label,
        axis_bundle_object_type_label,
        heavy_objects_sanitized = true,
        pending_facts = parent_axis_probe.pending_facts,
        diagnostics = merge(
            parent_axis_probe.diagnostics,
            (
                carry_objects_requested = carry_objects_requested,
                basis_object_available = basis_object_available,
                axis_bundle_object_available = axis_bundle_object_available,
                heavy_objects_sanitized = true,
            ),
        ),
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
    parent_axis_probe_carry_objects_requested =
        parent_axis_probe_requested &&
        parent_axis_readiness.homonuclear &&
        parent_axis_readiness.geometry.existing_bond_aligned_api_geometry_ready
    parent_axis_probe = parent_axis_probe_requested ?
        metrics._pqs_explicit_core_spacing_parent_axis_probe(
            standard_setup;
            gausslet_backend = probe_inputs.parent_axis_probe_backend,
            family = probe_inputs.parent_axis_probe_family,
            construct_axis_bundles = true,
            carry_objects = parent_axis_probe_carry_objects_requested,
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
        parent_axis_probe_carry_objects_requested,
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
        comparison_reference_label =
            get(route_recipe, :comparison_reference_label, nothing),
        comparison_ready = get(route_recipe, :comparison_ready, true),
        comparison_blocker = get(route_recipe, :comparison_blocker, nothing),
        artifact_role = get(route_recipe, :artifact_role, nothing),
        physics_endpoint_ready =
            get(route_recipe, :physics_endpoint_ready, nothing),
        physics_endpoint_blocker =
            get(route_recipe, :physics_endpoint_blocker, nothing),
        retained_atom_core_interiors =
            get(route_recipe, :retained_atom_core_interiors, nothing),
        source_plan_role = get(route_recipe, :source_plan_role, nothing),
        supplement_policy = get(route_recipe, :supplement_policy, nothing),
        wl_h1_lowest = get(route_recipe, :wl_h1_lowest, nothing),
        wl_h1_self_coulomb = get(route_recipe, :wl_h1_self_coulomb, nothing),
        wl_rhf_one_electron_energy =
            get(route_recipe, :wl_rhf_one_electron_energy, nothing),
        wl_rhf_electron_electron_energy =
            get(route_recipe, :wl_rhf_electron_electron_energy, nothing),
        wl_rhf_electronic_energy =
            get(route_recipe, :wl_rhf_electronic_energy, nothing),
        wl_rhf_nuclear_repulsion =
            get(route_recipe, :wl_rhf_nuclear_repulsion, nothing),
        wl_rhf_total_with_nuclear_repulsion =
            get(route_recipe, :wl_rhf_total_with_nuclear_repulsion, nothing),
        run_final_basis = get(route_recipe, :run_final_basis, true),
        run_h1 = get(route_recipe, :run_h1, true),
        run_h1_j = get(route_recipe, :run_h1_j, true),
        run_private_rhf =
            get(get(route_recipe, :private_rhf_inputs, (;)), :run_private_rhf, false),
        wl_rhf_total = get(route_recipe, :wl_rhf_total, nothing),
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
        parent_axis_readiness = parent_axis.parent_axis_readiness,
        parent_axis_probe =
            _cartesian_parent_axis_probe_report_summary(parent_axis.parent_axis_probe),
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

function _pqs_source_box_route_driver_report_terminal_route_aliases(
    low_order_route_summary,
)
    return (;
        low_order_terminal_shellification_selected =
            low_order_route_summary.terminal_shellification_selected,
        low_order_terminal_shellification_summary_available =
            low_order_route_summary.terminal_shellification_summary_available,
        low_order_terminal_shellification_scaffold_available =
            low_order_route_summary.terminal_shellification_scaffold_available,
        low_order_terminal_shellification_scaffold =
            low_order_route_summary.terminal_shellification_scaffold,
        low_order_terminal_shellification_region_count =
            low_order_route_summary.terminal_shellification_region_count,
        low_order_terminal_shellification_unit_inventory_available =
            low_order_route_summary.terminal_shellification_unit_inventory_available,
        low_order_terminal_shellification_unit_inventory =
            low_order_route_summary.terminal_shellification_unit_inventory,
        low_order_terminal_shellification_unit_count =
            low_order_route_summary.terminal_shellification_unit_count,
        low_order_terminal_shellification_unit_keys =
            low_order_route_summary.terminal_shellification_unit_keys,
        low_order_terminal_shellification_unit_roles =
            low_order_route_summary.terminal_shellification_unit_roles,
        low_order_terminal_shellification_unit_kinds =
            low_order_route_summary.terminal_shellification_unit_kinds,
        low_order_terminal_shellification_unit_support_counts =
            low_order_route_summary.terminal_shellification_unit_support_counts,
        low_order_terminal_shellification_lowering_contract_inventory_available =
            low_order_route_summary.terminal_shellification_lowering_contract_inventory_available,
        low_order_terminal_shellification_lowering_contract_inventory_status =
            low_order_route_summary.terminal_shellification_lowering_contract_inventory_status,
        low_order_terminal_shellification_lowering_contract_inventory =
            low_order_route_summary.terminal_shellification_lowering_contract_inventory,
        low_order_terminal_shellification_lowering_plan_available =
            low_order_route_summary.terminal_shellification_lowering_plan_available,
        low_order_terminal_shellification_lowering_plan_status =
            low_order_route_summary.terminal_shellification_lowering_plan_status,
        low_order_terminal_shellification_lowering_plan =
            low_order_route_summary.terminal_shellification_lowering_plan,
        low_order_terminal_shellification_lowering_summary =
            low_order_route_summary.terminal_shellification_lowering_summary,
        low_order_terminal_shellification_lowering_contract_count =
            low_order_route_summary.terminal_shellification_lowering_contract_count,
        low_order_terminal_shellification_lowering_contract_kinds =
            low_order_route_summary.terminal_shellification_lowering_contract_kinds,
        low_order_terminal_shellification_lowering_contract_kind_counts =
            low_order_route_summary.terminal_shellification_lowering_contract_kind_counts,
        low_order_terminal_shellification_contract_counts_by_unit =
            low_order_route_summary.terminal_shellification_contract_counts_by_unit,
        _pqs_source_box_route_driver_report_selected_terminal_lowering_fields(
            low_order_route_summary,
        )...,
        low_order_terminal_shellification_lw_complete_shell_cpb_count =
            low_order_route_summary.terminal_shellification_lw_complete_shell_cpb_count,
        low_order_terminal_shellification_lw_complete_shell_cpb_family_counts =
            low_order_route_summary.terminal_shellification_lw_complete_shell_cpb_family_counts,
        low_order_terminal_shellification_final_retained_unit_inventory_available =
            low_order_route_summary.terminal_shellification_final_retained_unit_inventory_available,
        low_order_terminal_shellification_transform_contracts_available =
            low_order_route_summary.terminal_shellification_transform_contracts_available,
        low_order_terminal_shellification_pair_inventory_available =
            low_order_route_summary.terminal_shellification_pair_inventory_available,
        low_order_terminal_shellification_pair_inventory_status =
            low_order_route_summary.terminal_shellification_pair_inventory_status,
        low_order_terminal_shellification_pair_materialization_status =
            low_order_route_summary.terminal_shellification_pair_materialization_status,
        low_order_terminal_shellification_assembly_materialization_status =
            low_order_route_summary.terminal_shellification_assembly_materialization_status,
        low_order_terminal_shellification_central_gap_region_count =
            low_order_route_summary.terminal_shellification_central_gap_region_count,
        low_order_terminal_shellification_central_midpoint_slab_count =
            low_order_route_summary.terminal_shellification_central_midpoint_slab_count,
        low_order_terminal_shellification_central_distorted_product_box_count =
            low_order_route_summary.terminal_shellification_central_distorted_product_box_count,
        low_order_terminal_shellification_central_distorted_product_box_metadata =
            low_order_route_summary.terminal_shellification_central_distorted_product_box_metadata,
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
        parent_axis_readiness = parent_axis.parent_axis_readiness,
        parent_axis_probe =
            _cartesian_parent_axis_probe_report_summary(parent_axis.parent_axis_probe),
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
        low_order_shellization_policy_requested =
            low_order_route_summary.low_order_shellization_policy_requested,
        low_order_shellization_policy_resolved =
            low_order_route_summary.low_order_shellization_policy_resolved,
        low_order_shellization_policy_source =
            low_order_route_summary.low_order_shellization_policy_source,
        low_order_shellization_policy_status =
            low_order_route_summary.low_order_shellization_policy_status,
        low_order_shellization_policy_blocker =
            low_order_route_summary.low_order_shellization_policy_blocker,
        low_order_shellization_source =
            low_order_route_summary.shellization_source,
        low_order_shellization_kind = low_order_route_summary.shellization_kind,
        low_order_unit_route_kind = low_order_route_summary.unit_route_kind,
        low_order_transform_route_kind =
            low_order_route_summary.transform_route_kind,
        low_order_pair_route_kind = low_order_route_summary.pair_route_kind,
        low_order_assembly_source = low_order_route_summary.assembly_source,
        low_order_assembly_route_kind =
            low_order_route_summary.assembly_route_kind,
        low_order_assembly_kind = low_order_route_summary.assembly_kind,
        atom_growth_low_order_route_selected =
            low_order_route_summary.atom_growth_selected,
        _pqs_source_box_route_driver_report_terminal_route_aliases(
            low_order_route_summary,
        )...,
        legacy_source_low_order_route_selected =
            low_order_route_summary.legacy_source_selected,
        low_order_plan_authority = low_order_route_summary.plan_authority,
        low_order_active_source_authority =
            low_order_route_summary.active_source_authority,
        low_order_legacy_source_authority =
            low_order_route_summary.legacy_source_authority,
        low_order_materialization_required =
            low_order_route_summary.materialization_required,
        low_order_materialization_status =
            low_order_route_summary.materialization_status,
        low_order_materialization_blocker =
            low_order_route_summary.materialization_blocker,
        low_order_hamiltonian_matrices_materialized =
            low_order_route_summary.hamiltonian_matrices_materialized,
        low_order_operator_matrices_materialized =
            low_order_route_summary.operator_matrices_materialized,
        low_order_pair_operator_blocks_materialized =
            low_order_route_summary.pair_operator_blocks_materialized,
        low_order_pair_inventory_source =
            low_order_route_summary.pair_inventory_source,
        low_order_pair_inventory_known =
            low_order_route_summary.pair_inventory_known,
        low_order_independent_atom_growth_pair_inventory_available =
            low_order_route_summary.independent_atom_growth_pair_inventory_available,
        low_order_pair_count = low_order_route_summary.pair_count,
        low_order_pair_family_counts = low_order_route_summary.pair_family_counts,
        low_order_route_core_final_unit_count =
            low_order_route_summary.route_core_final_unit_count,
        low_order_route_core_pair_inventory_available =
            low_order_route_summary.route_core_pair_inventory_available,
        low_order_route_core_pair_inventory_status =
            low_order_route_summary.route_core_pair_inventory_status,
        low_order_route_core_pair_count =
            low_order_route_summary.route_core_pair_count,
        low_order_route_core_pair_order_matches_staged =
            low_order_route_summary.route_core_pair_order_matches_staged,
        low_order_route_core_pair_order_comparison_source =
            low_order_route_summary.route_core_pair_order_comparison_source,
        low_order_route_core_pair_family_counts =
            low_order_route_summary.route_core_pair_family_counts,
        low_order_route_core_summary_status =
            low_order_route_summary.route_core_summary_status,
        low_order_route_core_pair_operator_ready =
            low_order_route_summary.route_core_pair_operator_ready,
        low_order_route_core_pair_operator_readiness_status =
            low_order_route_summary.route_core_pair_operator_readiness_status,
        low_order_route_core_pair_operator_blocker =
            low_order_route_summary.route_core_pair_operator_blocker,
        low_order_route_core_pair_operator_readiness_requirements =
            low_order_route_summary.route_core_pair_operator_readiness_requirements,
        low_order_route_core_pair_operator_preflight_available =
            low_order_route_summary.route_core_pair_operator_preflight_available,
        low_order_route_core_pair_operator_preflight_status =
            low_order_route_summary.route_core_pair_operator_preflight_status,
        low_order_route_core_pair_operator_preflight =
            low_order_route_summary.route_core_pair_operator_preflight,
        low_order_route_core_pair_operator_preflight_blocker =
            low_order_route_summary.route_core_pair_operator_preflight_blocker,
        low_order_route_core_pair_operator_plan_available =
            low_order_route_summary.route_core_pair_operator_plan_available,
        low_order_route_core_pair_operator_plan_status =
            low_order_route_summary.route_core_pair_operator_plan_status,
        low_order_route_core_pair_operator_plan =
            low_order_route_summary.route_core_pair_operator_plan,
        low_order_route_core_pair_operator_plan_blocker =
            low_order_route_summary.route_core_pair_operator_plan_blocker,
        low_order_route_core_typed_pair_operator_plan_inventory_available =
            low_order_route_summary.route_core_typed_pair_operator_plan_inventory_available,
        low_order_route_core_typed_pair_operator_plan_inventory_status =
            low_order_route_summary.route_core_typed_pair_operator_plan_inventory_status,
        low_order_route_core_typed_pair_operator_plan_blocker =
            low_order_route_summary.route_core_typed_pair_operator_plan_blocker,
        low_order_route_core_typed_pair_operator_plan_count =
            low_order_route_summary.route_core_typed_pair_operator_plan_count,
        low_order_route_core_typed_pair_operator_plan_blocked_count =
            low_order_route_summary.route_core_typed_pair_operator_plan_blocked_count,
        low_order_route_core_typed_pair_operator_plan_materialized =
            low_order_route_summary.route_core_typed_pair_operator_plan_materialized,
        low_order_route_core_typed_pair_operator_source_path_counts =
            low_order_route_summary.route_core_typed_pair_operator_source_path_counts,
        low_order_route_core_typed_pair_operator_final_block_path_counts =
            low_order_route_summary.route_core_typed_pair_operator_final_block_path_counts,
        low_order_route_core_typed_pair_operator_materialization_status_counts =
            low_order_route_summary.route_core_typed_pair_operator_materialization_status_counts,
        low_order_route_core_typed_pair_operator_blocker_counts =
            low_order_route_summary.route_core_typed_pair_operator_blocker_counts,
        low_order_route_core_typed_pair_operator_plan_family_counts =
            low_order_route_summary.route_core_typed_pair_operator_plan_family_counts,
        low_order_route_core_typed_pair_operator_materialization_ready =
            low_order_route_summary.route_core_typed_pair_operator_materialization_ready,
        low_order_route_core_typed_pair_operator_materialization_readiness_status =
            low_order_route_summary.route_core_typed_pair_operator_materialization_readiness_status,
        low_order_route_core_typed_pair_operator_materialization_readiness_blocker =
            low_order_route_summary.route_core_typed_pair_operator_materialization_readiness_blocker,
        low_order_route_core_typed_pair_operator_materialization_readiness_requirements =
            low_order_route_summary.route_core_typed_pair_operator_materialization_readiness_requirements,
        low_order_route_core_typed_pair_operator_materialization_readiness_plan_count =
            low_order_route_summary.route_core_typed_pair_operator_materialization_readiness_plan_count,
        low_order_route_core_typed_pair_operator_materialization_readiness_blocked_count =
            low_order_route_summary.route_core_typed_pair_operator_materialization_readiness_blocked_count,
        low_order_route_core_typed_pair_operator_materialization_readiness_materialized_count =
            low_order_route_summary.route_core_typed_pair_operator_materialization_readiness_materialized_count,
        low_order_lw_complete_shell_cpb_enumeration_available =
            low_order_route_summary.lw_complete_shell_cpb_enumeration_available,
        low_order_lw_complete_shell_region_count =
            low_order_route_summary.lw_complete_shell_region_count,
        low_order_lw_complete_shell_cpb_count =
            low_order_route_summary.lw_complete_shell_cpb_count,
        low_order_lw_complete_shell_cpb_family_counts =
            low_order_route_summary.lw_complete_shell_cpb_family_counts,
        low_order_lw_complete_shell_enumeration_policy =
            low_order_route_summary.lw_complete_shell_enumeration_policy,
        low_order_lw_complete_shell_coefficient_maps_materialized =
            low_order_route_summary.lw_complete_shell_coefficient_maps_materialized,
        low_order_lw_complete_shell_operator_blocks_materialized =
            low_order_route_summary.lw_complete_shell_operator_blocks_materialized,
        low_order_lw_complete_shell_pair_operator_blocks_materialized =
            low_order_route_summary.lw_complete_shell_pair_operator_blocks_materialized,
        low_order_lw_complete_shell_hamiltonian_data_materialized =
            low_order_route_summary.lw_complete_shell_hamiltonian_data_materialized,
        low_order_pqs_lowering_prototype_available =
            low_order_route_summary.pqs_lowering_prototype_available,
        low_order_pqs_transform_prototype_available =
            low_order_route_summary.pqs_transform_prototype_available,
        low_order_pqs_prototype_unit_key =
            low_order_route_summary.pqs_prototype_unit_key,
        low_order_pqs_prototype_stage =
            low_order_route_summary.pqs_prototype_stage,
        low_order_pqs_prototype_source_cpb_kind =
            low_order_route_summary.pqs_prototype_source_cpb_kind,
        low_order_pqs_prototype_owned_support_is_cpb =
            low_order_route_summary.pqs_prototype_owned_support_is_cpb,
        low_order_pqs_prototype_intermediate_retained_space =
            low_order_route_summary.pqs_prototype_intermediate_retained_space,
        low_order_pqs_prototype_shell_realization =
            low_order_route_summary.pqs_prototype_shell_realization,
        low_order_pqs_prototype_source_count_distinct_from_owned_support_count =
            low_order_route_summary.pqs_prototype_source_count_distinct_from_owned_support_count,
        low_order_pqs_prototype_coefficient_maps_materialized =
            low_order_route_summary.pqs_prototype_coefficient_maps_materialized,
        low_order_pqs_prototype_source_operator_blocks_materialized =
            low_order_route_summary.pqs_prototype_source_operator_blocks_materialized,
        low_order_pqs_prototype_operator_blocks_materialized =
            low_order_route_summary.pqs_prototype_operator_blocks_materialized,
        low_order_pqs_prototype_pair_operator_blocks_materialized =
            low_order_route_summary.pqs_prototype_pair_operator_blocks_materialized,
        low_order_pqs_prototype_hamiltonian_data_materialized =
            low_order_route_summary.pqs_prototype_hamiltonian_data_materialized,
        low_order_pqs_prototype_artifacts_materialized =
            low_order_route_summary.pqs_prototype_artifacts_materialized,
        diagnostics,
        route_materializer_payload,
    )
end

function _pqs_source_box_route_driver_materializer_payload(parent)
    parent_axis_probe =
        hasproperty(parent, :parent_axis_probe) ? parent.parent_axis_probe : nothing
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
        parent_qw_basis_object,
        parent_axis_bundle_object,
        parent_qw_basis_object_available = !isnothing(parent_qw_basis_object),
        parent_axis_bundle_object_available = !isnothing(parent_axis_bundle_object),
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
    return merge(
        (;
            object_kind = :cartesian_route_system,
            status = :initialized_driver_system,
        ),
        system_inputs,
    )
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
            support_dense_direct_allowed = route_inputs.support_dense_direct_allowed,
            reference_only_authorities = route_inputs.reference_only_authorities,
        )
        white_lindsey_recipe = (;
            route_shape = route_inputs.white_lindsey_route_shape,
            mapping_rule = route_inputs.white_lindsey_mapping_rule,
            nesting_rule = route_inputs.white_lindsey_nesting_rule,
            retained_rule = route_inputs.white_lindsey_retained_rule,
            operator_rule = route_inputs.white_lindsey_operator_rule,
            benchmark_role = route_inputs.white_lindsey_benchmark_role,
        )
        route_recipe = (;
            route_family = route_inputs.route_family,
            route_kind = route_inputs.route_kind,
            terms = route_inputs.terms,
            pair_factor_normalization = route_inputs.pair_factor_normalization,
            source_box = source_box_recipe,
            white_lindsey = white_lindsey_recipe,
            comparison_reference_label =
                get(route_inputs, :comparison_reference_label, nothing),
            comparison_ready = get(route_inputs, :comparison_ready, true),
            comparison_blocker = get(route_inputs, :comparison_blocker, nothing),
            artifact_role = get(route_inputs, :artifact_role, nothing),
            physics_endpoint_ready =
                get(route_inputs, :physics_endpoint_ready, nothing),
            physics_endpoint_blocker =
                get(route_inputs, :physics_endpoint_blocker, nothing),
            retained_atom_core_interiors =
                get(route_inputs, :retained_atom_core_interiors, nothing),
            source_plan_role = get(route_inputs, :source_plan_role, nothing),
            supplement_policy = get(route_inputs, :supplement_policy, nothing),
            wl_h1_lowest = get(route_inputs, :wl_h1_lowest, nothing),
            wl_h1_self_coulomb = get(route_inputs, :wl_h1_self_coulomb, nothing),
            wl_rhf_one_electron_energy =
                get(route_inputs, :wl_rhf_one_electron_energy, nothing),
            wl_rhf_electron_electron_energy =
                get(route_inputs, :wl_rhf_electron_electron_energy, nothing),
            wl_rhf_electronic_energy =
                get(route_inputs, :wl_rhf_electronic_energy, nothing),
            wl_rhf_nuclear_repulsion =
                get(route_inputs, :wl_rhf_nuclear_repulsion, nothing),
            wl_rhf_total_with_nuclear_repulsion =
                get(route_inputs, :wl_rhf_total_with_nuclear_repulsion, nothing),
            run_final_basis =
                isnothing(run_final_basis) ?
                (run_h1 || run_h1_j || run_private_rhf) :
                run_final_basis,
            run_h1,
            run_h1_j,
            private_rhf_inputs,
            wl_rhf_total = get(route_inputs, :wl_rhf_total, nothing),
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
            standard_setup, system, parent_inputs)
    route_axis_counts =
        _pqs_source_box_route_driver_route_axis_counts(
            standard_setup, parent_axis, system, recipe)
    object_carry = _cartesian_parent_object_carry(
        center_table, classification, standard_setup,
        parent_axis, route_axis_counts, parent_inputs)
    materialization_plan =
        _cartesian_parent_materialization_plan(
            center_table, center_axis_metadata, classification,
            standard_setup, parent_axis, route_axis_counts, object_carry)
    materialization_status =
        _cartesian_parent_materialization_status(parent_axis, object_carry)

    return (;
        object_kind = :cartesian_route_parent,
        status = route_axis_counts.status,
        system,
        spacing_inputs,
        parent_inputs,
        standard_setup,
        parent_axis,
        parent_axis_readiness = parent_axis.parent_axis_readiness,
        parent_axis_probe = parent_axis.parent_axis_probe,
        route_axis_counts,
        atom_count = length(center_table),
        atom_symbols = Tuple(center.atom_symbol for center in center_table),
        nuclear_charges = Tuple(center.nuclear_charge for center in center_table),
        atom_locations = Tuple(center.location for center in center_table),
        center_table,
        center_count = length(center_table),
        center_axis_metadata,
        system_classification = classification.system_classification,
        system_classification_status = classification.system_classification_status,
        bond_axis = classification.bond_axis,
        chain_axis = classification.chain_axis,
        axis_counts = route_axis_counts.parent_axis_counts,
        axis_counts_source = route_axis_counts.parent_axis_counts_source,
        axis_counts_status = route_axis_counts.status,
        physical_box = standard_setup.parent_box,
        physical_box_rule = standard_setup.parent_box_rule,
        parent_object_carry = object_carry,
        parent_basis_object = object_carry.parent_basis_object,
        parent_qw_basis_object = object_carry.parent_qw_basis_object,
        parent_axis_bundle_object = object_carry.parent_axis_bundle_object,
        parent_basis_object_available = object_carry.parent_basis_object_available,
        parent_qw_basis_object_available =
            object_carry.parent_qw_basis_object_available,
        parent_axis_bundle_object_available =
            object_carry.parent_axis_bundle_object_available,
        parent_basis_object_type_label = object_carry.parent_basis_object_type_label,
        parent_qw_basis_object_type_label =
            object_carry.parent_qw_basis_object_type_label,
        parent_axis_bundle_object_type_label =
            object_carry.parent_axis_bundle_object_type_label,
        parent_materialization_plan = materialization_plan,
        parent_materialization_plan_status = materialization_plan.status,
        parent_materialization_planning_family =
            materialization_plan.planning_family,
        parent_materialization_blocker = materialization_plan.blocker,
        parent_basis_materialization = materialization_status,
        parent_basis_materialization_status = materialization_status.status,
        parent_basis_materialized = materialization_status.parent_basis_materialized,
        parent_axis_metadata_constructed =
            materialization_status.parent_axis_metadata_constructed,
        axis_bundle_materialized = materialization_status.axis_bundle_materialized,
    )
end

function _pqs_source_box_route_driver_shell_stage_atom_growth_plan(parent)
    parent_qw_basis_object =
        hasproperty(parent, :parent_qw_basis_object) ?
        parent.parent_qw_basis_object :
        nothing
    parent_axis_bundle_object =
        hasproperty(parent, :parent_axis_bundle_object) ?
        parent.parent_axis_bundle_object :
        nothing
    missing_parent_objects = Symbol[]
    isnothing(parent_qw_basis_object) &&
        push!(missing_parent_objects, :parent_qw_basis_object)
    isnothing(parent_axis_bundle_object) &&
        push!(missing_parent_objects, :parent_axis_bundle_object)
    missing_parent_objects = Tuple(missing_parent_objects)
    if !isempty(missing_parent_objects)
        return (;
            object_kind = :cartesian_shell_stage_atom_growth_plan_payload,
            status = :blocked_atom_growth_missing_parent_objects,
            construction_plan = nothing,
            scaffold = nothing,
            construction_plan_available = false,
            scaffold_available = false,
            shellification_plan_materialization_available = false,
            missing_parent_objects,
            blocker = :blocked_atom_growth_missing_parent_objects,
            error_message = nothing,
            region_count = 0,
            unsupported_region_count = nothing,
            materialization_dependency_counts = nothing,
            construction_region_order = (),
            ordered_region_roles = (),
            spatial_policy_order = :atom_outward,
            coverage_status = :not_checked_missing_parent_objects,
            coverage_complete = nothing,
            parent_qw_basis_object_available = !isnothing(parent_qw_basis_object),
            parent_axis_bundle_object_available =
                !isnothing(parent_axis_bundle_object),
        )
    end

    try
        nside = parent.standard_setup.n_s
        anatomy =
            _nested_bond_aligned_diatomic_atom_growth_anatomy(
                parent_qw_basis_object,
                parent_axis_bundle_object;
                bond_axis = parent.bond_axis,
                protected_atom_side_count = nside,
            )
        construction_plan =
            _nested_bond_aligned_diatomic_atom_growth_construction_plan(
                anatomy,
            )
        retention = _nested_resolve_complete_shell_retention(nside)
        protect_rows =
            _nested_diatomic_resolve_core_near_nucleus_protect_rows(
                :auto,
                nside,
            )
        scaffold =
            _cartesian_shellification_plan_atom_growth_complete_rectangular_low_order(
                construction_plan,
                parent_axis_bundle_object;
                nside,
                child_retention_policy = retention,
                shared_retention_policy = retention,
                reference_fudge_factor = 1.2,
                core_near_nucleus_protect_rows = protect_rows,
                shared_shell_angular_resolution_scale = 1.4,
                route_family = :white_lindsey_low_order,
            )
        unsupported_region_count = scaffold.unsupported_region_count
        blocked_on_unsupported = unsupported_region_count > 0

        return (;
            object_kind = :cartesian_shell_stage_atom_growth_plan_payload,
            status =
                blocked_on_unsupported ?
                :blocked_atom_growth_unsupported_regions :
                :available_atom_growth_shellification_plan,
            construction_plan,
            scaffold,
            construction_plan_available = true,
            scaffold_available = true,
            shellification_plan_materialization_available =
                !blocked_on_unsupported,
            missing_parent_objects = (),
            blocker =
                blocked_on_unsupported ?
                :blocked_atom_growth_unsupported_regions :
                nothing,
            error_message = nothing,
            region_count = scaffold.region_count,
            unsupported_region_count,
            materialization_dependency_counts =
                scaffold.materialization_dependency_counts,
            construction_region_order = scaffold.construction_region_order,
            ordered_region_roles = scaffold.ordered_region_roles,
            spatial_policy_order = scaffold.spatial_policy_order,
            coverage_status =
                scaffold.coverage.coverage_complete ?
                :coverage_complete :
                :coverage_incomplete,
            coverage_complete = scaffold.coverage.coverage_complete,
            parent_qw_basis_object_available = true,
            parent_axis_bundle_object_available = true,
        )
    catch error
        error isa ArgumentError || rethrow()
        return (;
            object_kind = :cartesian_shell_stage_atom_growth_plan_payload,
            status = :blocked_atom_growth_unsupported_regions,
            construction_plan = nothing,
            scaffold = nothing,
            construction_plan_available = false,
            scaffold_available = false,
            shellification_plan_materialization_available = false,
            missing_parent_objects = (),
            blocker = :blocked_atom_growth_unsupported_regions,
            error_message = sprint(showerror, error),
            region_count = 0,
            unsupported_region_count = nothing,
            materialization_dependency_counts = nothing,
            construction_region_order = (),
            ordered_region_roles = (),
            spatial_policy_order = :atom_outward,
            coverage_status = :not_checked_atom_growth_plan_precondition,
            coverage_complete = nothing,
            parent_qw_basis_object_available = true,
            parent_axis_bundle_object_available = true,
        )
    end
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

function _pqs_source_box_route_driver_terminal_stage_status(scaffold)
    scaffold.materialization_status == :ready_supported_terminal_subset &&
        return :available_terminal_shellification_scaffold
    scaffold.materialization_status ==
    :deferred_pending_distorted_product_box_lowering &&
        return :deferred_terminal_shellification_distorted_product_box_lowering
    return scaffold.materialization_status
end

function _pqs_source_box_route_driver_shell_stage_terminal_shellification(parent)
    missing_parent_inputs = Symbol[]
    isnothing(parent.axis_counts) && push!(missing_parent_inputs, :parent_axis_counts)
    !hasproperty(parent, :standard_setup) &&
        push!(missing_parent_inputs, :standard_setup)
    if !isempty(missing_parent_inputs)
        return (;
            object_kind = :cartesian_shell_stage_terminal_shellification_payload,
            status = :blocked_terminal_shellification_missing_parent_inputs,
            plan = nothing,
            scaffold = nothing,
            plan_available = false,
            scaffold_available = false,
            shellification_plan_materialization_available = false,
            missing_parent_inputs = Tuple(missing_parent_inputs),
            blocker = :blocked_terminal_shellification_missing_parent_inputs,
            error_message = nothing,
            region_count = 0,
            materialization_dependency_counts = nothing,
            ordered_region_roles = (),
            spatial_policy_order = nothing,
            coverage_status = :not_checked_missing_parent_inputs,
            coverage_complete = nothing,
            materialization_status =
                :blocked_terminal_shellification_missing_parent_inputs,
            central_gap_region_count = 0,
            central_midpoint_slab_count = 0,
            central_distorted_product_box_count = 0,
            central_distorted_product_box_metadata = (),
            parent_axes_source = :unavailable,
            nuclear_positions_source = :unavailable,
        )
    end

    try
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
        plan = CartesianShellification.shellify(
            parent_axes,
            nuclear_positions,
            policy,
        )
        raw_plan = CartesianShellification.raw_plan(plan)
        scaffold = CartesianShellification.scaffold(
            plan;
            route_family = :white_lindsey_low_order,
        )
        status = _pqs_source_box_route_driver_terminal_stage_status(scaffold)
        materialization_available = false
        blocker =
            scaffold.materialization_status ==
            :deferred_pending_distorted_product_box_lowering ?
            :distorted_product_box_lowering_pending :
            nothing
        return (;
            object_kind = :cartesian_shell_stage_terminal_shellification_payload,
            status,
            policy,
            plan,
            raw_plan,
            scaffold,
            plan_available = true,
            scaffold_available = true,
            shellification_plan_materialization_available =
                materialization_available,
            missing_parent_inputs = (),
            blocker,
            error_message = nothing,
            region_count = scaffold.region_count,
            materialization_dependency_counts =
                scaffold.materialization_dependency_counts,
            ordered_region_roles = scaffold.ordered_region_roles,
            spatial_policy_order = scaffold.spatial_policy_order,
            coverage_status =
                scaffold.coverage.coverage_complete ?
                :coverage_complete :
                :coverage_incomplete,
            coverage_complete = scaffold.coverage.coverage_complete,
            materialization_status = status,
            central_gap_region_count = scaffold.central_gap_region_count,
            central_midpoint_slab_count = scaffold.central_midpoint_slab_count,
            central_distorted_product_box_count =
                scaffold.central_distorted_product_box_count,
            central_distorted_product_box_metadata =
                scaffold.central_distorted_product_box_metadata,
            parent_axes_source = :index_axes_from_parent_axis_counts,
            nuclear_positions_source =
                :index_positions_from_parent_counts_and_center_order,
        )
    catch error
        error isa ArgumentError || rethrow()
        return (;
            object_kind = :cartesian_shell_stage_terminal_shellification_payload,
            status = :blocked_terminal_shellification_geometry_precondition,
            plan = nothing,
            scaffold = nothing,
            plan_available = false,
            scaffold_available = false,
            shellification_plan_materialization_available = false,
            missing_parent_inputs = (),
            blocker = :blocked_terminal_shellification_geometry_precondition,
            error_message = sprint(showerror, error),
            region_count = 0,
            materialization_dependency_counts = nothing,
            ordered_region_roles = (),
            spatial_policy_order = nothing,
            coverage_status = :not_checked_terminal_geometry_precondition,
            coverage_complete = nothing,
            materialization_status =
                :blocked_terminal_shellification_geometry_precondition,
            central_gap_region_count = 0,
            central_midpoint_slab_count = 0,
            central_distorted_product_box_count = 0,
            central_distorted_product_box_metadata = (),
            parent_axes_source = :index_axes_from_parent_axis_counts,
            nuclear_positions_source =
                :index_positions_from_parent_counts_and_center_order,
        )
    end
end

function _pqs_source_box_route_driver_shell_stage_low_order_shellization(
    parent,
    recipe;
    low_order_shellization_policy = nothing,
    probe_route_configured_diatomic_atom_growth_materializer::Bool = false,
)
    policy =
        recipe.route_family == :white_lindsey_low_order ?
        _pqs_source_box_route_driver_low_order_shellization_policy(
            low_order_shellization_policy,
            probe_route_configured_diatomic_atom_growth_materializer,
        ) :
        (;
            low_order_shellization_policy_requested = low_order_shellization_policy,
            low_order_shellization_policy_resolved = :not_applicable,
            low_order_shellization_policy_source = :not_applicable,
            low_order_shellization_policy_status = :not_applicable,
            low_order_shellization_policy_blocker = nothing,
        )
    policy_resolved = policy.low_order_shellization_policy_resolved
    bond_aligned_diatomic =
        parent.system_classification == :bond_aligned_diatomic
    atom_growth_selected =
        recipe.route_family == :white_lindsey_low_order &&
        bond_aligned_diatomic &&
        policy.low_order_shellization_policy_status ==
        :available_low_order_shellization_policy &&
        policy_resolved == :atom_growth_complete_rectangular
    terminal_shellification_selected =
        (
            recipe.route_family == :white_lindsey_low_order &&
            parent.system_classification in (:one_center, :bond_aligned_diatomic) &&
            policy.low_order_shellization_policy_status ==
            :available_low_order_shellization_policy &&
            policy_resolved == :terminal_cartesian_shellification_geometry
        ) ||
        (
            recipe.route_family == :pqs_source_box &&
            parent.system_classification == :one_center
        )
    legacy_source_selected =
        recipe.route_family == :white_lindsey_low_order &&
        bond_aligned_diatomic &&
        policy.low_order_shellization_policy_status ==
        :available_low_order_shellization_policy &&
        policy_resolved == :legacy_diatomic_source
    shellization_source =
        atom_growth_selected ?
        :bond_aligned_diatomic_atom_growth_construction_plan :
        terminal_shellification_selected ?
        :terminal_cartesian_shellification_geometry :
        legacy_source_selected ?
        :route_configured_bond_aligned_diatomic_source :
        recipe.route_family == :white_lindsey_low_order ?
        (
            bond_aligned_diatomic ?
            :blocked_low_order_shellization_policy :
            :route_configured_low_order_non_diatomic_shellization
        ) :
        :not_applicable
    shellization_kind =
        atom_growth_selected ?
        :atom_growth_complete_rectangular :
        terminal_shellification_selected ?
        :terminal_cartesian_shellification_geometry :
        legacy_source_selected ?
        :legacy_diatomic_source :
        recipe.route_family == :white_lindsey_low_order ?
        :non_diatomic_low_order_shellization :
        :not_applicable
    materialization_required =
        atom_growth_selected || terminal_shellification_selected || legacy_source_selected
    materialization_status =
        atom_growth_selected ?
        :deferred_atom_growth_complete_rectangular_materialization :
        terminal_shellification_selected ?
        :deferred_terminal_cartesian_shellification_geometry :
        legacy_source_selected ?
        :deferred_legacy_diatomic_source_materialization :
        policy.low_order_shellization_policy_status
    atom_growth_plan_payload =
        atom_growth_selected ?
        _pqs_source_box_route_driver_shell_stage_atom_growth_plan(parent) :
        nothing
    terminal_shellification_payload =
        terminal_shellification_selected ?
        _pqs_source_box_route_driver_shell_stage_terminal_shellification(parent) :
        nothing
    atom_growth_plan_available =
        !isnothing(atom_growth_plan_payload) &&
        atom_growth_plan_payload.construction_plan_available
    atom_growth_scaffold_available =
        !isnothing(atom_growth_plan_payload) &&
        atom_growth_plan_payload.scaffold_available
    atom_growth_shellification_plan_available =
        !isnothing(atom_growth_plan_payload) &&
        atom_growth_plan_payload.status == :available_atom_growth_shellification_plan
    atom_growth_status =
        isnothing(atom_growth_plan_payload) ?
        nothing :
        atom_growth_plan_payload.status
    terminal_shellification_plan_available =
        !isnothing(terminal_shellification_payload) &&
        terminal_shellification_payload.plan_available
    terminal_shellification_scaffold_available =
        !isnothing(terminal_shellification_payload) &&
        terminal_shellification_payload.scaffold_available
    terminal_shellification_status =
        isnothing(terminal_shellification_payload) ?
        nothing :
        terminal_shellification_payload.status
    coverage_status =
        terminal_shellification_selected && !isnothing(terminal_shellification_payload) ?
        terminal_shellification_payload.coverage_status :
        isnothing(atom_growth_plan_payload) ?
        :not_checked_shell_stage_summary_only :
        atom_growth_plan_payload.coverage_status
    coverage_complete =
        terminal_shellification_selected && !isnothing(terminal_shellification_payload) ?
        terminal_shellification_payload.coverage_complete :
        isnothing(atom_growth_plan_payload) ?
        nothing :
        atom_growth_plan_payload.coverage_complete
    materialization_available =
        terminal_shellification_selected && !isnothing(terminal_shellification_payload) ?
        terminal_shellification_payload.shellification_plan_materialization_available :
        isnothing(atom_growth_plan_payload) ?
        false :
        atom_growth_plan_payload.shellification_plan_materialization_available
    materialization_status =
        terminal_shellification_selected && !isnothing(terminal_shellification_payload) ?
        terminal_shellification_payload.status :
        atom_growth_selected && !isnothing(atom_growth_plan_payload) ?
        atom_growth_plan_payload.status :
        materialization_status
    summary_only =
        terminal_shellification_selected ?
        !terminal_shellification_scaffold_available :
        atom_growth_selected ? !atom_growth_shellification_plan_available : true
    status =
        terminal_shellification_selected && !isnothing(terminal_shellification_payload) ?
        terminal_shellification_payload.status :
        atom_growth_selected && !isnothing(atom_growth_plan_payload) ?
        atom_growth_plan_payload.status :
        policy.low_order_shellization_policy_status ==
        :available_low_order_shellization_policy ?
        :available_shell_stage_low_order_shellization_summary :
        policy.low_order_shellization_policy_status

    return (;
        object_kind = :cartesian_shell_stage_low_order_shellization_summary,
        route_family = recipe.route_family,
        status,
        policy...,
        shellization_source,
        shellization_kind,
        atom_growth_selected,
        terminal_shellification_selected,
        legacy_source_selected,
        atom_growth_plan_summary_available = atom_growth_selected,
        atom_growth_plan_available,
        atom_growth_scaffold_available,
        atom_growth_shellification_plan_available,
        atom_growth_shellification_plan_status = atom_growth_status,
        atom_growth_plan_payload,
        atom_growth_construction_plan =
            isnothing(atom_growth_plan_payload) ?
            nothing :
            atom_growth_plan_payload.construction_plan,
        atom_growth_scaffold =
            isnothing(atom_growth_plan_payload) ?
            nothing :
            atom_growth_plan_payload.scaffold,
        atom_growth_missing_parent_objects =
            isnothing(atom_growth_plan_payload) ?
            () :
            atom_growth_plan_payload.missing_parent_objects,
        atom_growth_unsupported_region_count =
            isnothing(atom_growth_plan_payload) ?
            nothing :
            atom_growth_plan_payload.unsupported_region_count,
        atom_growth_region_count =
            isnothing(atom_growth_plan_payload) ?
            0 :
            atom_growth_plan_payload.region_count,
        atom_growth_spatial_policy_order =
            isnothing(atom_growth_plan_payload) ?
            nothing :
            atom_growth_plan_payload.spatial_policy_order,
        atom_growth_construction_region_order =
            isnothing(atom_growth_plan_payload) ?
            () :
            atom_growth_plan_payload.construction_region_order,
        atom_growth_materialization_dependency_counts =
            isnothing(atom_growth_plan_payload) ?
            nothing :
            atom_growth_plan_payload.materialization_dependency_counts,
        terminal_shellification_plan_available,
        terminal_shellification_scaffold_available,
        terminal_shellification_scaffold =
            isnothing(terminal_shellification_payload) ?
            nothing :
            terminal_shellification_payload.scaffold,
        terminal_shellification_plan =
            isnothing(terminal_shellification_payload) ?
            nothing :
            terminal_shellification_payload.plan,
        terminal_shellification_region_count =
            isnothing(terminal_shellification_payload) ?
            0 :
            terminal_shellification_payload.region_count,
        terminal_shellification_materialization_status =
            isnothing(terminal_shellification_payload) ?
            nothing :
            terminal_shellification_payload.materialization_status,
        terminal_shellification_spatial_policy_order =
            isnothing(terminal_shellification_payload) ?
            nothing :
            terminal_shellification_payload.spatial_policy_order,
        terminal_shellification_coverage_complete =
            isnothing(terminal_shellification_payload) ?
            nothing :
            terminal_shellification_payload.coverage_complete,
        terminal_shellification_central_gap_region_count =
            isnothing(terminal_shellification_payload) ?
            0 :
            terminal_shellification_payload.central_gap_region_count,
        terminal_shellification_central_midpoint_slab_count =
            isnothing(terminal_shellification_payload) ?
            0 :
            terminal_shellification_payload.central_midpoint_slab_count,
        terminal_shellification_central_distorted_product_box_count =
            isnothing(terminal_shellification_payload) ?
            0 :
            terminal_shellification_payload.central_distorted_product_box_count,
        terminal_shellification_central_distorted_product_box_metadata =
            isnothing(terminal_shellification_payload) ?
            () :
            terminal_shellification_payload.central_distorted_product_box_metadata,
        terminal_shellification_materialization_dependency_counts =
            isnothing(terminal_shellification_payload) ?
            nothing :
            terminal_shellification_payload.materialization_dependency_counts,
        terminal_shellification_payload,
        atom_growth_plan_authority = atom_growth_selected,
        terminal_shellification_authority = terminal_shellification_selected,
        active_source_authority = legacy_source_selected,
        legacy_source_authority = legacy_source_selected,
        coverage_status,
        coverage_complete,
        materialization_required,
        materialization_available,
        materialization_status,
        materialization_stage =
            atom_growth_shellification_plan_available ?
            :atom_growth_shellification_plan_stored_sequence_not_materialized :
            :deferred_to_route_materialization,
        materialized_sequence_available = false,
        materialized_operator_matrices_available = false,
        system_classification = parent.system_classification,
        system_classification_status = parent.system_classification_status,
        bond_axis = parent.bond_axis,
        private_development_only = true,
        full_plan_stored =
            atom_growth_plan_available || terminal_shellification_plan_available,
        scaffold_stored =
            atom_growth_scaffold_available ||
            terminal_shellification_scaffold_available,
        summary_only,
    )
end

function cartesian_shells(
    parent,
    spacing_inputs,
    recipe;
    low_order_shellization_policy = nothing,
    probe_route_configured_diatomic_atom_growth_materializer::Bool = false,
)
    route_skeleton =
        _pqs_source_box_route_driver_route_skeleton(
            parent.route_axis_counts, spacing_inputs, recipe)
    low_order_shellization =
        _pqs_source_box_route_driver_shell_stage_low_order_shellization(
            parent,
            recipe;
            low_order_shellization_policy,
            probe_route_configured_diatomic_atom_growth_materializer,
        )

    return (;
        object_kind = :cartesian_shells,
        status = route_skeleton.status,
        spacing_inputs,
        route_skeleton,
        low_order_shellization,
        low_order_shellization_policy_resolved =
            low_order_shellization.low_order_shellization_policy_resolved,
        low_order_shellization_policy_status =
            low_order_shellization.low_order_shellization_policy_status,
        shellization_source = low_order_shellization.shellization_source,
        shellization_kind = low_order_shellization.shellization_kind,
        atom_growth_plan_summary_available =
            low_order_shellization.atom_growth_plan_summary_available,
        atom_growth_plan_available =
            low_order_shellization.atom_growth_plan_available,
        atom_growth_scaffold_available =
            low_order_shellization.atom_growth_scaffold_available,
        atom_growth_shellification_plan_status =
            low_order_shellization.atom_growth_shellification_plan_status,
        atom_growth_region_count =
            low_order_shellization.atom_growth_region_count,
        atom_growth_unsupported_region_count =
            low_order_shellization.atom_growth_unsupported_region_count,
        terminal_shellification_selected =
            low_order_shellization.terminal_shellification_selected,
        terminal_shellification_plan_available =
            low_order_shellization.terminal_shellification_plan_available,
        terminal_shellification_scaffold_available =
            low_order_shellization.terminal_shellification_scaffold_available,
        terminal_shellification_scaffold =
            low_order_shellization.terminal_shellification_scaffold,
        terminal_shellification_region_count =
            low_order_shellization.terminal_shellification_region_count,
        terminal_shellification_materialization_status =
            low_order_shellization.terminal_shellification_materialization_status,
        terminal_shellification_spatial_policy_order =
            low_order_shellization.terminal_shellification_spatial_policy_order,
        terminal_shellification_coverage_complete =
            low_order_shellization.terminal_shellification_coverage_complete,
        terminal_shellification_central_gap_region_count =
            low_order_shellization.terminal_shellification_central_gap_region_count,
        terminal_shellification_central_midpoint_slab_count =
            low_order_shellization.terminal_shellification_central_midpoint_slab_count,
        terminal_shellification_central_distorted_product_box_count =
            low_order_shellization.terminal_shellification_central_distorted_product_box_count,
        coverage_complete = low_order_shellization.coverage_complete,
        materialization_available =
            low_order_shellization.materialization_available,
        full_plan_stored = low_order_shellization.full_plan_stored,
        scaffold_stored = low_order_shellization.scaffold_stored,
        summary_only = low_order_shellization.summary_only,
        atom_growth_plan_authority =
            low_order_shellization.atom_growth_plan_authority,
        terminal_shellification_authority =
            low_order_shellization.terminal_shellification_authority,
        active_source_authority = low_order_shellization.active_source_authority,
        route_shape = route_skeleton.route_shape,
        source_boxes = route_skeleton.source_boxes,
        shellization_stage =
            _pqs_source_box_route_driver_shellization_stage(
                low_order_shellization,
            ),
    )
end

function _pqs_source_box_route_driver_shellization_stage(low_order_shellization)
    low_order_shellization.terminal_shellification_selected &&
        return low_order_shellization.status
    low_order_shellization.atom_growth_selected && return (
        low_order_shellization.atom_growth_shellification_plan_available ?
        :available_atom_growth_shellification_plan :
        low_order_shellization.status
    )
    low_order_shellization.legacy_source_selected &&
        return :legacy_diatomic_source_summary
    return :represented_by_route_skeleton
end

function _pqs_source_box_route_driver_terminal_lowering_family(low_order_shellization)
    hasproperty(low_order_shellization, :route_family) || return nothing
    low_order_shellization.route_family == :white_lindsey_low_order &&
        return :white_lindsey_low_order
    low_order_shellization.route_family == :pqs_source_box && return :pqs
    return nothing
end

function _pqs_source_box_route_driver_empty_terminal_lowering_contract_kind_counts()
    return _cartesian_terminal_region_lowering_contract_kind_counts(())
end

function _pqs_source_box_route_driver_terminal_lowering_policy(
    route_lowering_family,
    low_order_shellization,
)
    route_lowering_family == :white_lindsey_low_order &&
        return CartesianTerminalLowering.WhiteLindseyLowering()
    if route_lowering_family == :pqs
        plan =
            hasproperty(low_order_shellization, :terminal_shellification_plan) ?
            low_order_shellization.terminal_shellification_plan :
            nothing
        raw_plan =
            plan isa CartesianShellification.ShellificationPlan ?
            CartesianShellification.raw_plan(plan) :
            nothing
        q = !isnothing(raw_plan) && hasproperty(raw_plan, :q) ? raw_plan.q : nothing
        isnothing(q) &&
            throw(ArgumentError("PQS terminal lowering requires shellification q"))
        return CartesianTerminalLowering.PQSLowering(q = q)
    end
    return nothing
end

function _pqs_source_box_route_driver_terminal_lowering_plan(
    low_order_shellization,
    route_lowering_family,
)
    plan =
        hasproperty(low_order_shellization, :terminal_shellification_plan) ?
        low_order_shellization.terminal_shellification_plan :
        nothing
    plan isa CartesianShellification.ShellificationPlan || return nothing

    policy =
        _pqs_source_box_route_driver_terminal_lowering_policy(
            route_lowering_family,
            low_order_shellization,
        )
    isnothing(policy) && return nothing
    return CartesianTerminalLowering.lower_terminal_regions(plan, policy)
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

function _pqs_source_box_route_driver_report_selected_terminal_lowering_fields(
    low_order_route_summary,
)
    return (;
        low_order_terminal_shellification_selected_lowering_contract_inventory_available =
            low_order_route_summary.terminal_shellification_selected_lowering_contract_inventory_available,
        low_order_terminal_shellification_selected_lowering_contract_inventory_status =
            low_order_route_summary.terminal_shellification_selected_lowering_contract_inventory_status,
        low_order_terminal_shellification_selected_lowering_contract_inventory =
            low_order_route_summary.terminal_shellification_selected_lowering_contract_inventory,
        low_order_terminal_shellification_selected_lowering_family =
            low_order_route_summary.terminal_shellification_selected_lowering_family,
        low_order_terminal_shellification_selected_contract_count =
            low_order_route_summary.terminal_shellification_selected_contract_count,
        low_order_terminal_shellification_selected_contract_kinds =
            low_order_route_summary.terminal_shellification_selected_contract_kinds,
        low_order_terminal_shellification_selected_contract_kind_counts =
            low_order_route_summary.terminal_shellification_selected_contract_kind_counts,
        low_order_terminal_shellification_selected_contract_counts_by_unit =
            low_order_route_summary.terminal_shellification_selected_contract_counts_by_unit,
        low_order_terminal_shellification_all_units_have_exactly_one_selected_contract =
            low_order_route_summary.terminal_shellification_all_units_have_exactly_one_selected_contract,
        low_order_terminal_shellification_unselected_contract_count =
            low_order_route_summary.terminal_shellification_unselected_contract_count,
        low_order_terminal_shellification_unselected_contract_kinds =
            low_order_route_summary.terminal_shellification_unselected_contract_kinds,
    )
end

function _pqs_source_box_route_driver_selected_terminal_crc_sidecar_summary(
    selected_terminal_lowering_contract_inventory,
)
    sidecar_inventory =
        _cartesian_route_core_selected_terminal_lowering_sidecar_inventory(
            selected_terminal_lowering_contract_inventory,
        )
    missing_entries = Tuple(
        entry for entry in sidecar_inventory.sidecar_entries
        if !entry.route_core_sidecar_available
    )
    return (;
        object_kind =
            :cartesian_unit_stage_selected_terminal_lowering_crc_sidecar_summary,
        status = sidecar_inventory.status,
        private_development_only = true,
        sidecar_inventory,
        selected_contract_count = sidecar_inventory.selected_contract_count,
        sidecar_available_count = sidecar_inventory.sidecar_available_count,
        sidecar_missing_count = sidecar_inventory.sidecar_missing_count,
        sidecar_inventory_complete =
            sidecar_inventory.route_core_sidecar_inventory_complete,
        missing_sidecar_reasons =
            Tuple(entry.missing_route_core_sidecar_reason for entry in missing_entries),
        missing_sidecar_kinds =
            Tuple(entry.lowering_contract_kind for entry in missing_entries),
        final_retained_unit_inventory_available =
            sidecar_inventory.final_retained_unit_inventory_available,
        pair_inventory_available = sidecar_inventory.pair_inventory_available,
        pair_inventory_status = sidecar_inventory.pair_inventory_status,
        operator_blocks_materialized = false,
        pair_operator_blocks_materialized =
            sidecar_inventory.pair_operator_blocks_materialized,
        hamiltonian_data_materialized =
            sidecar_inventory.hamiltonian_data_materialized,
        artifacts_materialized = sidecar_inventory.artifacts_materialized,
    )
end

function _pqs_source_box_route_driver_source_selected_crc_sidecar_summary(
    source,
    selected::Bool,
)
    !selected &&
        return _pqs_source_box_route_driver_selected_terminal_crc_sidecar_summary(nothing)
    if hasproperty(source, :terminal_shellification_selected_crc_sidecar_summary)
        return source.terminal_shellification_selected_crc_sidecar_summary
    end
    if hasproperty(source, :terminal_route_state) &&
       hasproperty(source.terminal_route_state, :selected_crc_sidecar_summary)
        return source.terminal_route_state.selected_crc_sidecar_summary
    end
    return _pqs_source_box_route_driver_selected_terminal_crc_sidecar_summary(nothing)
end

function _pqs_source_box_route_driver_terminal_shellification_alias_fields(
    source,
    selected::Bool;
    selected_status = :not_selected,
    include_crc_sidecar_summary::Bool = true,
)
    selected_fields =
        _pqs_source_box_route_driver_selected_terminal_lowering_fields(
            selected ?
            source.terminal_shellification_selected_lowering_contract_inventory :
            nothing,
            selected ?
            source.terminal_shellification_selected_lowering_contract_inventory_status :
            selected_status,
            selected ?
            source.terminal_shellification_selected_lowering_family :
            nothing,
        )
    crc_sidecar_fields =
        include_crc_sidecar_summary ?
        (;
            terminal_shellification_selected_crc_sidecar_summary =
                _pqs_source_box_route_driver_source_selected_crc_sidecar_summary(
                    source,
                    selected,
                ),
        ) :
        (;)

    return merge(
        (;
            terminal_shellification_scaffold_available =
                selected && source.terminal_shellification_scaffold_available,
            terminal_shellification_scaffold =
                selected ? source.terminal_shellification_scaffold : nothing,
            terminal_shellification_region_count =
                selected ? source.terminal_shellification_region_count : 0,
            terminal_shellification_unit_inventory_available =
                selected && source.terminal_shellification_unit_inventory_available,
            terminal_shellification_unit_inventory =
                selected ? source.terminal_shellification_unit_inventory : nothing,
            terminal_shellification_unit_count =
                selected ? source.terminal_shellification_unit_count : 0,
            terminal_shellification_unit_keys =
                selected ? source.terminal_shellification_unit_keys : (),
            terminal_shellification_unit_roles =
                selected ? source.terminal_shellification_unit_roles : (),
            terminal_shellification_unit_kinds =
                selected ? source.terminal_shellification_unit_kinds : (),
            terminal_shellification_unit_support_counts =
                selected ? source.terminal_shellification_unit_support_counts : (),
            terminal_shellification_lowering_contract_inventory_available =
                selected &&
                source.terminal_shellification_lowering_contract_inventory_available,
            terminal_shellification_lowering_contract_inventory_status =
                selected ?
                source.terminal_shellification_lowering_contract_inventory_status :
                :not_selected,
            terminal_shellification_lowering_contract_inventory =
                selected ?
                source.terminal_shellification_lowering_contract_inventory :
                nothing,
            terminal_shellification_lowering_contract_count =
                selected ?
                source.terminal_shellification_lowering_contract_count :
                0,
            terminal_shellification_lowering_contract_kinds =
                selected ?
                source.terminal_shellification_lowering_contract_kinds :
                (),
            terminal_shellification_lowering_contract_kind_counts =
                selected ?
                source.terminal_shellification_lowering_contract_kind_counts :
                _pqs_source_box_route_driver_empty_terminal_lowering_contract_kind_counts(),
        ),
        selected_fields,
        crc_sidecar_fields,
        (;
            terminal_shellification_contract_counts_by_unit =
                selected ?
                source.terminal_shellification_contract_counts_by_unit :
                (),
            terminal_shellification_lw_complete_shell_cpb_count =
                selected ?
                source.terminal_shellification_lw_complete_shell_cpb_count :
                0,
            terminal_shellification_lw_complete_shell_cpb_family_counts =
                selected ?
                source.terminal_shellification_lw_complete_shell_cpb_family_counts :
                (facet_cpb = 0, edge_cpb = 0, corner_cpb = 0),
            terminal_shellification_final_retained_unit_inventory_available =
                selected &&
                source.terminal_shellification_final_retained_unit_inventory_available,
            terminal_shellification_central_gap_region_count =
                selected ?
                source.terminal_shellification_central_gap_region_count :
                0,
            terminal_shellification_central_midpoint_slab_count =
                selected ?
                source.terminal_shellification_central_midpoint_slab_count :
                0,
            terminal_shellification_central_distorted_product_box_count =
                selected ?
                source.terminal_shellification_central_distorted_product_box_count :
                0,
            terminal_shellification_central_distorted_product_box_metadata =
                selected ?
                source.terminal_shellification_central_distorted_product_box_metadata :
                (),
        ),
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
    selected_crc_sidecar_summary =
        isnothing(selected_crc_sidecar_summary) ?
        _pqs_source_box_route_driver_selected_terminal_crc_sidecar_summary(
            nothing,
        ) :
        selected_crc_sidecar_summary
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
    low_order_shellization =
        hasproperty(shells, :low_order_shellization) ?
        shells.low_order_shellization :
        nothing
    if isnothing(low_order_shellization)
        terminal_route_state =
            _pqs_source_box_route_driver_terminal_route_state_unavailable(
                :not_available_missing_shell_stage_summary,
                :missing_shell_stage_low_order_summary,
            )
        return (;
            object_kind = :cartesian_unit_stage_low_order_summary,
            status = :not_available_missing_shell_stage_summary,
            terminal_route_state,
            terminal_route_summary = terminal_route_state.summary,
            low_order_shellization_policy_requested = nothing,
            low_order_shellization_policy_resolved = :not_available,
            low_order_shellization_policy_source = :not_available,
            low_order_shellization_policy_status =
                :not_available_missing_shell_stage_summary,
            low_order_shellization_policy_blocker =
                :missing_shell_stage_low_order_summary,
            shellization_source = :not_available,
            shellization_kind = :not_available,
            unit_route_kind = :not_available,
            atom_growth_units_selected = false,
            terminal_shellification_units_selected = false,
            legacy_source_units_selected = false,
            atom_growth_unit_summary_available = false,
            terminal_shellification_unit_summary_available = false,
            terminal_shellification_scaffold_available = false,
            terminal_shellification_scaffold = nothing,
            terminal_shellification_region_count = 0,
            terminal_shellification_unit_inventory_available = false,
            terminal_shellification_unit_inventory_status = :not_available,
            terminal_shellification_unit_inventory = nothing,
            terminal_shellification_unit_count = 0,
            terminal_shellification_unit_keys = (),
            terminal_shellification_unit_roles = (),
            terminal_shellification_unit_kinds = (),
            terminal_shellification_unit_support_counts = (),
            terminal_shellification_lowering_contract_inventory_available = false,
            terminal_shellification_lowering_contract_inventory_status =
                :not_available,
            terminal_shellification_lowering_contract_inventory = nothing,
            terminal_shellification_lowering_contract_count = 0,
            terminal_shellification_lowering_contract_kinds = (),
            terminal_shellification_lowering_contract_kind_counts =
                (
                    direct_core_identity_cpb_count = 0,
                    direct_slab_identity_cpb_count = 0,
                    direct_boundary_slab_identity_cpb_count = 0,
                    white_lindsey_boundary_strata_count = 0,
                    pqs_filled_source_cpb_count = 0,
                    distorted_product_box_comx_count = 0,
                ),
            terminal_shellification_lowering_plan_available = false,
            terminal_shellification_lowering_plan_status = :not_available,
            terminal_shellification_lowering_plan = nothing,
            terminal_shellification_lowering_summary = nothing,
            _pqs_source_box_route_driver_selected_terminal_lowering_fields(
                nothing,
                :not_available,
                nothing,
            )...,
            terminal_shellification_contract_counts_by_unit = (),
            terminal_shellification_lw_complete_shell_cpb_count = 0,
            terminal_shellification_lw_complete_shell_cpb_family_counts =
                (facet_cpb = 0, edge_cpb = 0, corner_cpb = 0),
            terminal_shellification_final_retained_unit_inventory_available = false,
            terminal_shellification_pair_inventory_available = false,
            terminal_shellification_central_gap_region_count = 0,
            terminal_shellification_central_midpoint_slab_count = 0,
            terminal_shellification_central_distorted_product_box_count = 0,
            terminal_shellification_central_distorted_product_box_metadata = (),
            materialized_units_available = false,
            materialization_status = :not_available,
            retained_unit_dimensions_known = false,
            retained_unit_ranges_known = false,
            retained_dimension_known = false,
            retained_dimension = nothing,
            plan_authority = false,
            active_source_authority = false,
            legacy_source_authority = false,
            unit_inventory_source = :not_available,
            unit_inventory_status = :not_available,
            atom_growth_unit_inventory_available = false,
            plan_unit_inventory_available = false,
            plan_unit_inventory = nothing,
            plan_unit_count = 0,
            plan_unit_roles = (),
            plan_unit_keys = (),
            plan_unit_support_counts = (),
            source_backed_region_count = 0,
            cpb_contract_stage = :not_available,
            shellification_regions_are_cpbs = false,
            owned_support_available = false,
            lowering_source_cpbs_available = false,
            source_cpb_count = 0,
            pqs_lowering_prototype_available = false,
            pqs_lowering_prototype = nothing,
            pqs_lowering_prototype_unit_key = nothing,
            route_core_sidecar_inventory_available = false,
            route_core_sidecar_inventory =
                _cartesian_route_core_sidecar_inventory(nothing),
            route_skeleton_unit_fields_preserved = false,
            route_skeleton_unit_inventory_source = :not_available,
            summary_only = true,
        )
    end

    atom_growth_units_selected = low_order_shellization.atom_growth_selected
    terminal_shellification_units_selected =
        low_order_shellization.terminal_shellification_selected
    legacy_source_units_selected = low_order_shellization.legacy_source_selected
    unit_route_kind =
        atom_growth_units_selected ?
        :atom_growth_complete_rectangular_low_order_units :
        terminal_shellification_units_selected ?
        :terminal_shellification_low_order_units :
        legacy_source_units_selected ?
        :legacy_diatomic_source_low_order_units :
        :not_selected
    terminal_shellification_scaffold_available =
        terminal_shellification_units_selected &&
        low_order_shellization.terminal_shellification_scaffold_available
    terminal_shellification_unit_summary_available =
        terminal_shellification_units_selected
    terminal_region_unit_inventory =
        terminal_shellification_units_selected &&
        terminal_shellification_scaffold_available ?
        _cartesian_terminal_shellification_region_unit_inventory(
            low_order_shellization.terminal_shellification_scaffold,
        ) :
        nothing
    terminal_region_unit_inventory_available =
        !isnothing(terminal_region_unit_inventory) &&
        terminal_region_unit_inventory.status ==
        :available_terminal_region_unit_inventory
    terminal_shellification_unit_inventory_status =
        terminal_shellification_units_selected ?
        (
            terminal_region_unit_inventory_available ?
            terminal_region_unit_inventory.status :
            :deferred_terminal_shellification_unit_inventory
        ) :
        :not_selected
    route_lowering_family =
        terminal_shellification_units_selected ?
        _pqs_source_box_route_driver_terminal_lowering_family(low_order_shellization) :
        nothing
    terminal_lowering_plan =
        terminal_shellification_units_selected && !isnothing(route_lowering_family) ?
        _pqs_source_box_route_driver_terminal_lowering_plan(
            low_order_shellization,
            route_lowering_family,
        ) :
        nothing
    terminal_lowering_plan_available =
        terminal_lowering_plan isa CartesianTerminalLowering.TerminalLoweringPlan
    terminal_lowering_summary =
        terminal_lowering_plan_available ?
        CartesianTerminalLowering.summary(terminal_lowering_plan) :
        nothing
    terminal_lowering_plan_status =
        terminal_shellification_units_selected ?
        (
            terminal_lowering_plan_available ?
            terminal_lowering_summary.status :
            isnothing(route_lowering_family) ?
            :not_available_terminal_lowering_family :
            :deferred_terminal_lowering_plan
        ) :
        :not_selected
    terminal_region_lowering_contract_inventory =
        terminal_region_unit_inventory_available && terminal_lowering_plan_available ?
        _pqs_source_box_route_driver_terminal_lowering_contract_inventory_from_plan(
            terminal_lowering_plan,
            terminal_region_unit_inventory,
        ) :
        nothing
    terminal_region_lowering_contract_inventory_available =
        !isnothing(terminal_region_lowering_contract_inventory) &&
        terminal_region_lowering_contract_inventory.status ==
        :available_terminal_region_lowering_contract_inventory
    terminal_shellification_lowering_contract_inventory_status =
        terminal_shellification_units_selected ?
        (
            terminal_region_lowering_contract_inventory_available ?
            terminal_region_lowering_contract_inventory.status :
            terminal_lowering_plan_available ?
            :deferred_terminal_lowering_plan_compatibility_inventory :
            :deferred_terminal_shellification_lowering_contract_inventory
        ) :
        :not_selected
    selected_terminal_lowering_contract_inventory =
        terminal_region_lowering_contract_inventory_available &&
        terminal_lowering_plan_available &&
        !isnothing(route_lowering_family) ?
        _pqs_source_box_route_driver_selected_terminal_lowering_contract_inventory_from_plan(
            terminal_lowering_plan,
            terminal_region_lowering_contract_inventory,
            route_lowering_family,
        ) :
        nothing
    selected_terminal_lowering_contract_inventory_available =
        !isnothing(selected_terminal_lowering_contract_inventory) &&
        selected_terminal_lowering_contract_inventory.status ==
        :available_selected_terminal_lowering_contract_inventory
    selected_terminal_lowering_contract_inventory_status =
        terminal_shellification_units_selected ?
        (
            selected_terminal_lowering_contract_inventory_available ?
            selected_terminal_lowering_contract_inventory.status :
            isnothing(route_lowering_family) ?
            :not_available_selected_terminal_lowering_family :
            terminal_region_lowering_contract_inventory_available ?
            :deferred_terminal_shellification_selected_lowering_contract_inventory :
            :deferred_terminal_shellification_lowering_contract_inventory
        ) :
        :not_selected
    selected_terminal_lowering_fields =
        _pqs_source_box_route_driver_selected_terminal_lowering_fields(
            selected_terminal_lowering_contract_inventory,
            selected_terminal_lowering_contract_inventory_status,
            route_lowering_family,
            terminal_lowering_plan,
        )
    selected_terminal_crc_sidecar_summary =
        _pqs_source_box_route_driver_selected_terminal_crc_sidecar_summary(
            selected_terminal_lowering_contract_inventory,
        )
    terminal_route_state =
        _pqs_source_box_route_driver_terminal_route_state(;
            status = terminal_shellification_units_selected ?
                     terminal_lowering_plan_status :
                     :not_selected,
            selected = terminal_shellification_units_selected,
            route_lowering_family,
            shellification_plan =
                terminal_shellification_units_selected ?
                low_order_shellization.terminal_shellification_plan :
                nothing,
            unit_inventory = terminal_region_unit_inventory,
            lowering_plan = terminal_lowering_plan,
            lowering_summary = terminal_lowering_summary,
            lowering_contract_inventory =
                terminal_region_lowering_contract_inventory,
            selected_contract_inventory =
                selected_terminal_lowering_contract_inventory,
            selected_crc_sidecar_summary = selected_terminal_crc_sidecar_summary,
            blocker =
                terminal_shellification_units_selected &&
                !terminal_lowering_plan_available ?
                terminal_lowering_plan_status :
                nothing,
        )
    atom_growth_plan_unit_inventory =
        atom_growth_units_selected ?
        _pqs_source_box_route_driver_atom_growth_plan_unit_inventory(
            low_order_shellization,
        ) :
        nothing
    atom_growth_plan_unit_inventory_available =
        !isnothing(atom_growth_plan_unit_inventory) &&
        atom_growth_plan_unit_inventory.status ==
        :available_atom_growth_plan_unit_inventory
    plan_unit_inventory =
        atom_growth_plan_unit_inventory_available ?
        atom_growth_plan_unit_inventory :
        terminal_region_unit_inventory_available ?
        terminal_region_unit_inventory :
        nothing
    plan_unit_inventory_available = !isnothing(plan_unit_inventory)
    unit_inventory_source =
        atom_growth_plan_unit_inventory_available ?
        atom_growth_plan_unit_inventory.unit_inventory_source :
        terminal_region_unit_inventory_available ?
        terminal_region_unit_inventory.inventory_source :
        atom_growth_units_selected ?
        :blocked_atom_growth_shellification_plan :
        terminal_shellification_units_selected ?
        :terminal_shellification_scaffold :
        legacy_source_units_selected ?
        :legacy_diatomic_source_summary :
        :route_skeleton_compatibility_fields
    unit_inventory_status =
        terminal_shellification_units_selected ?
        terminal_shellification_unit_inventory_status :
        isnothing(atom_growth_plan_unit_inventory) ?
            unit_inventory_source :
            atom_growth_plan_unit_inventory.status
    source_backed_region_count =
        isnothing(atom_growth_plan_unit_inventory) ?
        0 :
        atom_growth_plan_unit_inventory.source_backed_region_count
    cpb_contract_stage =
        atom_growth_plan_unit_inventory_available ?
        atom_growth_plan_unit_inventory.cpb_contract_stage :
        terminal_region_unit_inventory_available ?
        :terminal_region_unit_inventory_metadata :
        :not_available
    source_cpb_count =
        atom_growth_plan_unit_inventory_available ?
        atom_growth_plan_unit_inventory.source_cpb_count :
        0
    route_core_sidecar_inventory =
        atom_growth_plan_unit_inventory_available ?
        _cartesian_route_core_sidecar_inventory(atom_growth_plan_unit_inventory) :
        _cartesian_route_core_sidecar_inventory(nothing)
    status =
        atom_growth_plan_unit_inventory_available ||
        terminal_region_unit_inventory_available ?
        :available_unit_stage_low_order_summary :
        terminal_shellification_units_selected &&
        terminal_shellification_scaffold_available ?
        :deferred_terminal_shellification_unit_inventory :
        low_order_shellization.status ==
        :available_shell_stage_low_order_shellization_summary ?
        :available_unit_stage_low_order_summary :
        low_order_shellization.status

    return (;
        object_kind = :cartesian_unit_stage_low_order_summary,
        status,
        terminal_route_state,
        terminal_route_summary = terminal_route_state.summary,
        low_order_shellization_policy_requested =
            low_order_shellization.low_order_shellization_policy_requested,
        low_order_shellization_policy_resolved =
            low_order_shellization.low_order_shellization_policy_resolved,
        low_order_shellization_policy_source =
            low_order_shellization.low_order_shellization_policy_source,
        low_order_shellization_policy_status =
            low_order_shellization.low_order_shellization_policy_status,
        low_order_shellization_policy_blocker =
            low_order_shellization.low_order_shellization_policy_blocker,
        shellization_source = low_order_shellization.shellization_source,
        shellization_kind = low_order_shellization.shellization_kind,
        unit_route_kind,
        atom_growth_units_selected,
        terminal_shellification_units_selected,
        legacy_source_units_selected,
        atom_growth_unit_summary_available = atom_growth_units_selected,
        terminal_shellification_unit_summary_available,
        terminal_shellification_scaffold_available,
        terminal_shellification_scaffold =
            terminal_shellification_units_selected ?
            low_order_shellization.terminal_shellification_scaffold :
            nothing,
        terminal_shellification_region_count =
            terminal_shellification_units_selected ?
            low_order_shellization.terminal_shellification_region_count :
            0,
        terminal_shellification_unit_inventory_available =
            terminal_region_unit_inventory_available,
        terminal_shellification_unit_inventory_status,
        terminal_shellification_unit_inventory = terminal_region_unit_inventory,
        terminal_shellification_unit_count =
            terminal_region_unit_inventory_available ?
            terminal_region_unit_inventory.unit_count :
            0,
        terminal_shellification_unit_keys =
            terminal_region_unit_inventory_available ?
            terminal_region_unit_inventory.unit_keys :
            (),
        terminal_shellification_unit_roles =
            terminal_region_unit_inventory_available ?
            terminal_region_unit_inventory.unit_roles :
            (),
        terminal_shellification_unit_kinds =
            terminal_region_unit_inventory_available ?
            terminal_region_unit_inventory.unit_kinds :
            (),
        terminal_shellification_unit_support_counts =
            terminal_region_unit_inventory_available ?
            terminal_region_unit_inventory.support_counts :
            (),
        terminal_shellification_lowering_contract_inventory_available =
            terminal_region_lowering_contract_inventory_available,
        terminal_shellification_lowering_contract_inventory_status =
            terminal_shellification_lowering_contract_inventory_status,
        terminal_shellification_lowering_contract_inventory =
            terminal_region_lowering_contract_inventory,
        terminal_shellification_lowering_plan_available =
            terminal_lowering_plan_available,
        terminal_shellification_lowering_plan_status =
            terminal_lowering_plan_status,
        terminal_shellification_lowering_plan = terminal_lowering_plan,
        terminal_shellification_lowering_summary = terminal_lowering_summary,
        terminal_shellification_lowering_contract_count =
            terminal_lowering_plan_available ?
            terminal_lowering_summary.available_contract_count :
            0,
        terminal_shellification_lowering_contract_kinds =
            terminal_lowering_plan_available ?
            terminal_lowering_summary.available_contract_kinds :
            (),
        terminal_shellification_lowering_contract_kind_counts =
            terminal_lowering_plan_available ?
            _pqs_source_box_route_driver_terminal_lowering_kind_counts(
                CartesianTerminalLowering.available_contracts(terminal_lowering_plan),
            ) :
            (
                direct_core_identity_cpb_count = 0,
                direct_slab_identity_cpb_count = 0,
                direct_boundary_slab_identity_cpb_count = 0,
                white_lindsey_boundary_strata_count = 0,
                pqs_filled_source_cpb_count = 0,
                distorted_product_box_comx_count = 0,
            ),
        selected_terminal_lowering_fields...,
        terminal_shellification_contract_counts_by_unit =
            terminal_region_lowering_contract_inventory_available ?
            terminal_region_lowering_contract_inventory.contract_counts_by_unit :
            (),
        terminal_shellification_lw_complete_shell_cpb_count =
            terminal_region_lowering_contract_inventory_available ?
            terminal_region_lowering_contract_inventory.lw_complete_shell_cpb_count :
            0,
        terminal_shellification_lw_complete_shell_cpb_family_counts =
            terminal_region_lowering_contract_inventory_available ?
            terminal_region_lowering_contract_inventory.lw_complete_shell_cpb_family_counts :
            (facet_cpb = 0, edge_cpb = 0, corner_cpb = 0),
        terminal_shellification_final_retained_unit_inventory_available =
            terminal_region_unit_inventory_available &&
            terminal_region_unit_inventory.final_retained_unit_inventory_available,
        terminal_shellification_pair_inventory_available =
            terminal_region_unit_inventory_available &&
            terminal_region_unit_inventory.pair_inventory_available,
        terminal_shellification_central_gap_region_count =
            terminal_shellification_units_selected ?
            low_order_shellization.terminal_shellification_central_gap_region_count :
            0,
        terminal_shellification_central_midpoint_slab_count =
            terminal_shellification_units_selected ?
            low_order_shellization.terminal_shellification_central_midpoint_slab_count :
            0,
        terminal_shellification_central_distorted_product_box_count =
            terminal_shellification_units_selected ?
            low_order_shellization.terminal_shellification_central_distorted_product_box_count :
            0,
        terminal_shellification_central_distorted_product_box_metadata =
            terminal_shellification_units_selected ?
            low_order_shellization.terminal_shellification_central_distorted_product_box_metadata :
            (),
        materialized_units_available = false,
        materialization_status =
            atom_growth_units_selected ?
            :deferred_atom_growth_complete_rectangular_unit_materialization :
            terminal_shellification_units_selected ?
            :deferred_terminal_shellification_unit_materialization :
            legacy_source_units_selected ?
            :deferred_legacy_diatomic_source_unit_materialization :
            low_order_shellization.materialization_status,
        retained_unit_dimensions_known = false,
        retained_unit_ranges_known = false,
        retained_dimension_known = false,
        retained_dimension = nothing,
        plan_authority =
            terminal_shellification_units_selected ?
            low_order_shellization.terminal_shellification_authority :
            low_order_shellization.atom_growth_plan_authority,
        active_source_authority =
            terminal_shellification_units_selected ?
            false :
            low_order_shellization.active_source_authority,
        legacy_source_authority =
            terminal_shellification_units_selected ?
            false :
            low_order_shellization.legacy_source_authority,
        unit_inventory_source,
        unit_inventory_status,
        atom_growth_unit_inventory_available =
            atom_growth_plan_unit_inventory_available,
        plan_unit_inventory_available,
        plan_unit_inventory,
        plan_unit_count =
            plan_unit_inventory_available ? plan_unit_inventory.unit_count : 0,
        plan_unit_roles =
            plan_unit_inventory_available ? plan_unit_inventory.unit_roles : (),
        plan_unit_keys =
            plan_unit_inventory_available ? plan_unit_inventory.unit_keys : (),
        plan_unit_support_counts =
            plan_unit_inventory_available ?
            plan_unit_inventory.support_counts :
            (),
        source_backed_region_count,
        cpb_contract_stage,
        shellification_regions_are_cpbs =
            atom_growth_plan_unit_inventory_available ?
            atom_growth_plan_unit_inventory.shellification_regions_are_cpbs :
            false,
        owned_support_available =
            atom_growth_plan_unit_inventory_available &&
            atom_growth_plan_unit_inventory.owned_support_available,
        lowering_source_cpbs_available =
            atom_growth_plan_unit_inventory_available &&
            atom_growth_plan_unit_inventory.lowering_source_cpbs_available,
        source_cpb_count,
        pqs_lowering_prototype_available =
            atom_growth_plan_unit_inventory_available &&
            atom_growth_plan_unit_inventory.pqs_lowering_prototype_available,
        pqs_lowering_prototype =
            atom_growth_plan_unit_inventory_available ?
            atom_growth_plan_unit_inventory.pqs_lowering_prototype :
            nothing,
        pqs_lowering_prototype_unit_key =
            atom_growth_plan_unit_inventory_available ?
            atom_growth_plan_unit_inventory.pqs_lowering_prototype_unit_key :
            nothing,
        route_core_sidecar_inventory_available =
            route_core_sidecar_inventory.status in (
                :available_route_core_sidecar_inventory,
                :blocked_incomplete_route_core_sidecar_inventory,
            ),
        route_core_sidecar_inventory,
        route_skeleton_unit_fields_preserved = true,
        route_skeleton_unit_inventory_source =
            :route_skeleton_compatibility_fields,
        summary_only =
            terminal_shellification_units_selected || !plan_unit_inventory_available,
    )
end

function cartesian_units(parent, shells, route_inputs, recipe)
    raw_box =
        _pqs_source_box_route_driver_raw_box_probe(
            parent.standard_setup, shells.route_skeleton, parent.parent_axis,
            parent.route_axis_counts, route_inputs, recipe)
    low_order_units =
        _pqs_source_box_route_driver_unit_stage_low_order_summary(shells)

    return (;
        object_kind = :cartesian_units,
        status = shells.route_skeleton.status,
        route_inputs,
        route_skeleton = shells.route_skeleton,
        raw_box,
        low_order_units,
        terminal_route_state = low_order_units.terminal_route_state,
        terminal_route_summary = low_order_units.terminal_route_summary,
        low_order_unit_route_kind = low_order_units.unit_route_kind,
        atom_growth_unit_summary_available =
            low_order_units.atom_growth_unit_summary_available,
        atom_growth_units_selected = low_order_units.atom_growth_units_selected,
        atom_growth_unit_inventory_available =
            low_order_units.atom_growth_unit_inventory_available,
        plan_unit_inventory_available =
            low_order_units.plan_unit_inventory_available,
        plan_unit_inventory = low_order_units.plan_unit_inventory,
        unit_inventory_source = low_order_units.unit_inventory_source,
        unit_inventory_status = low_order_units.unit_inventory_status,
        cpb_contract_stage = low_order_units.cpb_contract_stage,
        shellification_regions_are_cpbs =
            low_order_units.shellification_regions_are_cpbs,
        owned_support_available = low_order_units.owned_support_available,
        lowering_source_cpbs_available =
            low_order_units.lowering_source_cpbs_available,
        source_cpb_count = low_order_units.source_cpb_count,
        pqs_lowering_prototype_available =
            low_order_units.pqs_lowering_prototype_available,
        pqs_lowering_prototype = low_order_units.pqs_lowering_prototype,
        pqs_lowering_prototype_unit_key =
            low_order_units.pqs_lowering_prototype_unit_key,
        route_core_sidecar_inventory_available =
            low_order_units.route_core_sidecar_inventory_available,
        route_core_sidecar_inventory = low_order_units.route_core_sidecar_inventory,
        materialized_units_available =
            low_order_units.materialized_units_available,
        retained_unit_dimensions_known =
            low_order_units.retained_unit_dimensions_known,
        retained_unit_ranges_known =
            low_order_units.retained_unit_ranges_known,
        retained_dimension_known = low_order_units.retained_dimension_known,
        summary_only = low_order_units.summary_only,
        active_source_authority = low_order_units.active_source_authority,
        source_boxes = shells.route_skeleton.source_boxes,
        source_dimensions = shells.route_skeleton.source_dimensions,
        retained_units = shells.route_skeleton.retained_units,
        retained_unit_order = shells.route_skeleton.retained_unit_order,
        unit_stage = :broken_into_source_units,
    )
end

function _pqs_source_box_route_driver_transform_stage_low_order_summary(units)
    low_order_units =
        hasproperty(units, :low_order_units) ?
        units.low_order_units :
        nothing
    if isnothing(low_order_units)
        terminal_route_state =
            _pqs_source_box_route_driver_terminal_route_state_unavailable(
                :not_available_missing_unit_stage_summary,
                :missing_unit_stage_low_order_summary,
            )
        return (;
            object_kind = :cartesian_transform_stage_low_order_summary,
            status = :not_available_missing_unit_stage_summary,
            terminal_route_state,
            terminal_route_summary = terminal_route_state.summary,
            low_order_shellization_policy_requested = nothing,
            low_order_shellization_policy_resolved = :not_available,
            low_order_shellization_policy_source = :not_available,
            low_order_shellization_policy_status =
                :not_available_missing_unit_stage_summary,
            low_order_shellization_policy_blocker =
                :missing_unit_stage_low_order_summary,
            shellization_source = :not_available,
            shellization_kind = :not_available,
            unit_route_kind = :not_available,
            transform_route_kind = :not_available,
            atom_growth_transforms_selected = false,
            terminal_shellification_transforms_selected = false,
            legacy_source_transforms_selected = false,
            terminal_shellification_transform_summary_available = false,
            terminal_shellification_scaffold_available = false,
            terminal_shellification_scaffold = nothing,
            terminal_shellification_region_count = 0,
            terminal_shellification_unit_inventory_available = false,
            terminal_shellification_unit_inventory = nothing,
            terminal_shellification_unit_count = 0,
            terminal_shellification_unit_keys = (),
            terminal_shellification_unit_roles = (),
            terminal_shellification_unit_kinds = (),
            terminal_shellification_unit_support_counts = (),
            terminal_shellification_lowering_contract_inventory_available = false,
            terminal_shellification_lowering_contract_inventory_status =
                :not_available,
            terminal_shellification_lowering_contract_inventory = nothing,
            terminal_shellification_lowering_contract_count = 0,
            terminal_shellification_lowering_contract_kinds = (),
            terminal_shellification_lowering_contract_kind_counts =
                (
                    direct_core_identity_cpb_count = 0,
                    direct_slab_identity_cpb_count = 0,
                    direct_boundary_slab_identity_cpb_count = 0,
                    white_lindsey_boundary_strata_count = 0,
                    pqs_filled_source_cpb_count = 0,
                    distorted_product_box_comx_count = 0,
                ),
            _pqs_source_box_route_driver_selected_terminal_lowering_fields(
                nothing,
                :not_available,
                nothing,
            )...,
            terminal_shellification_contract_counts_by_unit = (),
            terminal_shellification_lw_complete_shell_cpb_count = 0,
            terminal_shellification_lw_complete_shell_cpb_family_counts =
                (facet_cpb = 0, edge_cpb = 0, corner_cpb = 0),
            terminal_shellification_final_retained_unit_inventory_available = false,
            terminal_shellification_transform_contracts_available = false,
            terminal_shellification_transform_materialization_status =
                :not_available,
            terminal_shellification_central_gap_region_count = 0,
            terminal_shellification_central_midpoint_slab_count = 0,
            terminal_shellification_central_distorted_product_box_count = 0,
            terminal_shellification_central_distorted_product_box_metadata = (),
            coefficient_transforms_materialized = false,
            coefficient_maps_materialized = false,
            transform_materialization_status = :not_available,
            retained_unit_dimensions_known = false,
            retained_unit_ranges_known = false,
            retained_dimension_known = false,
            retained_dimension = nothing,
            plan_authority = false,
            active_source_authority = false,
            legacy_source_authority = false,
            transform_contract_source = :not_available,
            transform_contract_status = :not_available,
            atom_growth_transform_contracts_available = false,
            transform_contract_inventory_available = false,
            transform_contract_inventory = nothing,
            lw_complete_shell_cpb_enumeration_available = false,
            lw_complete_shell_region_count = 0,
            lw_complete_shell_cpb_count = 0,
            lw_complete_shell_cpb_family_counts =
                (facet_cpb = 0, edge_cpb = 0, corner_cpb = 0),
            lw_complete_shell_enumeration_policy = nothing,
            lw_complete_shell_coefficient_maps_materialized = false,
            lw_complete_shell_operator_blocks_materialized = false,
            lw_complete_shell_pair_operator_blocks_materialized = false,
            lw_complete_shell_hamiltonian_data_materialized = false,
            pqs_transform_prototype_available = false,
            pqs_transform_prototype = nothing,
            source_lowering_prototype_unit_key = nothing,
            transform_contract_count = 0,
            transform_contract_unit_keys = (),
            transform_contract_unit_roles = (),
            transform_contract_names = (),
            source_backed_contract_count = 0,
            cpb_contract_stage = :not_available,
            transform_contracts_derive_from_lowering = false,
            final_retained_units_are_pair_planning_inputs = false,
            transform_fields_preserved = false,
            route_skeleton_transform_inventory_source = :not_available,
            summary_only = true,
        )
    end

    atom_growth_transforms_selected = low_order_units.atom_growth_units_selected
    terminal_shellification_transforms_selected =
        low_order_units.terminal_shellification_units_selected
    legacy_source_transforms_selected = low_order_units.legacy_source_units_selected
    transform_route_kind =
        atom_growth_transforms_selected ?
        :atom_growth_complete_rectangular_low_order_transforms :
        terminal_shellification_transforms_selected ?
        :terminal_shellification_low_order_transforms :
        legacy_source_transforms_selected ?
        :legacy_diatomic_source_low_order_transforms :
        :not_selected
    terminal_shellification_scaffold_available =
        terminal_shellification_transforms_selected &&
        low_order_units.terminal_shellification_scaffold_available
    terminal_shellification_transform_summary_available =
        terminal_shellification_transforms_selected
    terminal_shellification_transform_materialization_status =
        terminal_shellification_transforms_selected ?
        :deferred_terminal_shellification_transform_contracts :
        :not_selected
    transform_contract_inventory =
        atom_growth_transforms_selected ?
        _pqs_source_box_route_driver_atom_growth_transform_contract_inventory(
            low_order_units.plan_unit_inventory,
        ) :
        nothing
    transform_contract_inventory_available =
        !isnothing(transform_contract_inventory) &&
        transform_contract_inventory.status ==
        :available_atom_growth_transform_contract_inventory
    transform_contract_source =
        transform_contract_inventory_available ?
        transform_contract_inventory.transform_contract_source :
        atom_growth_transforms_selected ?
        :blocked_atom_growth_plan_unit_inventory :
        terminal_shellification_transforms_selected ?
        :terminal_shellification_scaffold :
        legacy_source_transforms_selected ?
        :legacy_diatomic_source_summary :
        :route_skeleton_compatibility_fields
    transform_contract_status =
        terminal_shellification_transforms_selected ?
        terminal_shellification_transform_materialization_status :
        isnothing(transform_contract_inventory) ?
            transform_contract_source :
            transform_contract_inventory.status
    source_backed_contract_count =
        isnothing(transform_contract_inventory) ?
        0 :
        transform_contract_inventory.source_backed_contract_count
    transform_contracts_derive_from_lowering =
        transform_contract_inventory_available &&
        all(
            contract -> contract.final_unit_downstream_of_lowering,
            transform_contract_inventory.transform_contracts,
        )
    final_retained_units_are_pair_planning_inputs =
        transform_contract_inventory_available &&
        all(
            contract -> contract.final_retained_unit.pair_planning_input,
            transform_contract_inventory.transform_contracts,
        )
    status =
        transform_contract_inventory_available ?
        :available_transform_stage_low_order_summary :
        terminal_shellification_transforms_selected &&
        terminal_shellification_scaffold_available ?
        :deferred_terminal_shellification_transform_contracts :
        !isnothing(transform_contract_inventory) ?
        transform_contract_inventory.status :
        low_order_units.status == :available_unit_stage_low_order_summary ?
        :available_transform_stage_low_order_summary :
        low_order_units.status

    return (;
        object_kind = :cartesian_transform_stage_low_order_summary,
        status,
        terminal_route_state = low_order_units.terminal_route_state,
        terminal_route_summary = low_order_units.terminal_route_summary,
        low_order_shellization_policy_requested =
            low_order_units.low_order_shellization_policy_requested,
        low_order_shellization_policy_resolved =
            low_order_units.low_order_shellization_policy_resolved,
        low_order_shellization_policy_source =
            low_order_units.low_order_shellization_policy_source,
        low_order_shellization_policy_status =
            low_order_units.low_order_shellization_policy_status,
        low_order_shellization_policy_blocker =
            low_order_units.low_order_shellization_policy_blocker,
        shellization_source = low_order_units.shellization_source,
        shellization_kind = low_order_units.shellization_kind,
        unit_route_kind = low_order_units.unit_route_kind,
        transform_route_kind,
        atom_growth_transforms_selected,
        terminal_shellification_transforms_selected,
        legacy_source_transforms_selected,
        terminal_shellification_transform_summary_available,
        _pqs_source_box_route_driver_terminal_shellification_alias_fields(
            low_order_units,
            terminal_shellification_transforms_selected;
            include_crc_sidecar_summary = false,
        )...,
        terminal_shellification_transform_contracts_available = false,
        terminal_shellification_transform_materialization_status,
        coefficient_transforms_materialized = false,
        coefficient_maps_materialized = false,
        transform_materialization_status =
            atom_growth_transforms_selected ?
            :deferred_atom_growth_complete_rectangular_transform_materialization :
            terminal_shellification_transforms_selected ?
            terminal_shellification_transform_materialization_status :
            legacy_source_transforms_selected ?
            :deferred_legacy_diatomic_source_transform_materialization :
            low_order_units.materialization_status,
        retained_unit_dimensions_known = low_order_units.retained_unit_dimensions_known,
        retained_unit_ranges_known = low_order_units.retained_unit_ranges_known,
        retained_dimension_known = low_order_units.retained_dimension_known,
        retained_dimension = low_order_units.retained_dimension,
        plan_authority = low_order_units.plan_authority,
        active_source_authority = low_order_units.active_source_authority,
        legacy_source_authority = low_order_units.legacy_source_authority,
        transform_contract_source,
        transform_contract_status,
        atom_growth_transform_contracts_available =
            transform_contract_inventory_available,
        transform_contract_inventory_available,
        transform_contract_inventory,
        lw_complete_shell_cpb_enumeration_available =
            transform_contract_inventory_available &&
            transform_contract_inventory.lw_complete_shell_cpb_enumeration_available,
        lw_complete_shell_region_count =
            transform_contract_inventory_available ?
            transform_contract_inventory.lw_complete_shell_region_count :
            0,
        lw_complete_shell_cpb_count =
            transform_contract_inventory_available ?
            transform_contract_inventory.lw_complete_shell_cpb_count :
            0,
        lw_complete_shell_cpb_family_counts =
            transform_contract_inventory_available ?
            transform_contract_inventory.lw_complete_shell_cpb_family_counts :
            (facet_cpb = 0, edge_cpb = 0, corner_cpb = 0),
        lw_complete_shell_enumeration_policy =
            transform_contract_inventory_available ?
            transform_contract_inventory.lw_complete_shell_enumeration_policy :
            nothing,
        lw_complete_shell_coefficient_maps_materialized =
            transform_contract_inventory_available &&
            transform_contract_inventory.lw_complete_shell_coefficient_maps_materialized,
        lw_complete_shell_operator_blocks_materialized =
            transform_contract_inventory_available &&
            transform_contract_inventory.lw_complete_shell_operator_blocks_materialized,
        lw_complete_shell_pair_operator_blocks_materialized =
            transform_contract_inventory_available &&
            transform_contract_inventory.lw_complete_shell_pair_operator_blocks_materialized,
        lw_complete_shell_hamiltonian_data_materialized =
            transform_contract_inventory_available &&
            transform_contract_inventory.lw_complete_shell_hamiltonian_data_materialized,
        pqs_transform_prototype_available =
            transform_contract_inventory_available &&
            transform_contract_inventory.pqs_transform_prototype_available,
        pqs_transform_prototype =
            transform_contract_inventory_available ?
            transform_contract_inventory.pqs_transform_prototype :
            nothing,
        source_lowering_prototype_unit_key =
            transform_contract_inventory_available ?
            transform_contract_inventory.source_lowering_prototype_unit_key :
            nothing,
        transform_contract_count =
            transform_contract_inventory_available ?
            transform_contract_inventory.contract_count :
            0,
        transform_contract_unit_keys =
            transform_contract_inventory_available ?
            transform_contract_inventory.unit_keys :
            (),
        transform_contract_unit_roles =
            transform_contract_inventory_available ?
            transform_contract_inventory.unit_roles :
            (),
        transform_contract_names =
            transform_contract_inventory_available ?
            transform_contract_inventory.contract_names :
            (),
        source_backed_contract_count,
        cpb_contract_stage =
            transform_contract_inventory_available ?
            :construction_transform_contract :
            :not_available,
        transform_contracts_derive_from_lowering,
        final_retained_units_are_pair_planning_inputs,
        transform_fields_preserved = true,
        route_skeleton_transform_inventory_source =
            :route_skeleton_compatibility_fields,
        summary_only =
            terminal_shellification_transforms_selected ||
            !transform_contract_inventory_available,
    )
end

function cartesian_transforms(units, recipe)
    retained_units = units.retained_units
    low_order_transforms =
        _pqs_source_box_route_driver_transform_stage_low_order_summary(units)
    return (;
        object_kind = :cartesian_transforms,
        status = units.status,
        route_family = recipe.route_family,
        retained_units,
        low_order_transforms,
        terminal_route_state = low_order_transforms.terminal_route_state,
        terminal_route_summary = low_order_transforms.terminal_route_summary,
        low_order_transform_route_kind =
            low_order_transforms.transform_route_kind,
        atom_growth_transforms_selected =
            low_order_transforms.atom_growth_transforms_selected,
        coefficient_transforms_materialized =
            low_order_transforms.coefficient_transforms_materialized,
        coefficient_maps_materialized =
            low_order_transforms.coefficient_maps_materialized,
        atom_growth_transform_contracts_available =
            low_order_transforms.atom_growth_transform_contracts_available,
        transform_contract_inventory_available =
            low_order_transforms.transform_contract_inventory_available,
        transform_contract_inventory =
            low_order_transforms.transform_contract_inventory,
        lw_complete_shell_cpb_enumeration_available =
            low_order_transforms.lw_complete_shell_cpb_enumeration_available,
        lw_complete_shell_region_count =
            low_order_transforms.lw_complete_shell_region_count,
        lw_complete_shell_cpb_count =
            low_order_transforms.lw_complete_shell_cpb_count,
        lw_complete_shell_cpb_family_counts =
            low_order_transforms.lw_complete_shell_cpb_family_counts,
        lw_complete_shell_enumeration_policy =
            low_order_transforms.lw_complete_shell_enumeration_policy,
        lw_complete_shell_coefficient_maps_materialized =
            low_order_transforms.lw_complete_shell_coefficient_maps_materialized,
        lw_complete_shell_operator_blocks_materialized =
            low_order_transforms.lw_complete_shell_operator_blocks_materialized,
        lw_complete_shell_pair_operator_blocks_materialized =
            low_order_transforms.lw_complete_shell_pair_operator_blocks_materialized,
        lw_complete_shell_hamiltonian_data_materialized =
            low_order_transforms.lw_complete_shell_hamiltonian_data_materialized,
        pqs_transform_prototype_available =
            low_order_transforms.pqs_transform_prototype_available,
        pqs_transform_prototype =
            low_order_transforms.pqs_transform_prototype,
        source_lowering_prototype_unit_key =
            low_order_transforms.source_lowering_prototype_unit_key,
        transform_contract_source =
            low_order_transforms.transform_contract_source,
        transform_contract_status =
            low_order_transforms.transform_contract_status,
        cpb_contract_stage = low_order_transforms.cpb_contract_stage,
        transform_contracts_derive_from_lowering =
            low_order_transforms.transform_contracts_derive_from_lowering,
        final_retained_units_are_pair_planning_inputs =
            low_order_transforms.final_retained_units_are_pair_planning_inputs,
        summary_only = low_order_transforms.summary_only,
        active_source_authority = low_order_transforms.active_source_authority,
        retained_counts =
            _pqs_source_box_route_driver_named_tuple_from_units(
                retained_units, :retained_count),
        ranges =
            _pqs_source_box_route_driver_named_tuple_from_units(
                retained_units, :retained_range),
        retained_dimension =
            _pqs_source_box_route_driver_inventory_retained_dimension(retained_units),
        transform_stage = :unit_retained_transforms_described,
    )
end

function _pqs_source_box_route_driver_route_core_pair_family_counts(
    route_core_pair_inventory,
)
    isnothing(route_core_pair_inventory) && return ()
    pair_entries = CartesianRouteCore.pair_entries(route_core_pair_inventory)
    family_keys = Tuple(
        (
            CartesianRouteCore.lowering_recipe(pair.left),
            CartesianRouteCore.lowering_recipe(pair.right),
        ) for pair in pair_entries
    )
    unique_family_keys = Tuple(unique(family_keys))
    return Tuple(
        (;
            pair_family,
            pair_count = count(==(pair_family), family_keys),
            pair_family_source = :crc_final_unit_lowering_recipes,
        ) for pair_family in unique_family_keys
    )
end

function _pqs_source_box_route_driver_route_core_typed_pair_operator_plan_unavailable_metadata(
    status,
    blocker,
)
    return (;
        route_core_typed_pair_operator_plan_inventory_available = false,
        route_core_typed_pair_operator_plan_inventory_status = status,
        route_core_typed_pair_operator_plan_blocker = blocker,
        route_core_typed_pair_operator_plan_count = 0,
        route_core_typed_pair_operator_plan_blocked_count = 0,
        route_core_typed_pair_operator_plan_materialized = false,
        route_core_typed_pair_operator_source_path_counts = (),
        route_core_typed_pair_operator_final_block_path_counts = (),
        route_core_typed_pair_operator_materialization_status_counts = (),
        route_core_typed_pair_operator_blocker_counts = (),
        route_core_typed_pair_operator_plan_family_counts = (),
        route_core_typed_pair_operator_materialization_ready = false,
        route_core_typed_pair_operator_materialization_readiness_status =
            status,
        route_core_typed_pair_operator_materialization_readiness_blocker =
            blocker,
        route_core_typed_pair_operator_materialization_readiness_requirements =
            CartesianRouteCore.pair_operator_materialization_readiness_requirements(),
        route_core_typed_pair_operator_materialization_readiness_plan_count = 0,
        route_core_typed_pair_operator_materialization_readiness_blocked_count = 0,
        route_core_typed_pair_operator_materialization_readiness_materialized_count = 0,
    )
end

function _pqs_source_box_route_driver_route_core_typed_pair_operator_plan_metadata(
    route_core_sidecar_inventory,
)
    if isnothing(route_core_sidecar_inventory) ||
       !hasproperty(
           route_core_sidecar_inventory,
           :crc_pair_operator_plan_inventory_available,
       )
        return _pqs_source_box_route_driver_route_core_typed_pair_operator_plan_unavailable_metadata(
            :blocked_missing_route_core_sidecar_inventory,
            :missing_route_core_sidecar_inventory,
        )
    end

    if !route_core_sidecar_inventory.crc_pair_operator_plan_inventory_available
        return _pqs_source_box_route_driver_route_core_typed_pair_operator_plan_unavailable_metadata(
            route_core_sidecar_inventory.crc_pair_operator_plan_inventory_status,
            route_core_sidecar_inventory.crc_pair_operator_plan_blocker,
        )
    end

    plan_inventory =
        route_core_sidecar_inventory.crc_pair_operator_plan_inventory
    if isnothing(plan_inventory)
        return _pqs_source_box_route_driver_route_core_typed_pair_operator_plan_unavailable_metadata(
            :blocked_missing_route_core_typed_pair_operator_plan_inventory,
            :missing_route_core_typed_pair_operator_plan_inventory,
        )
    end

    return (;
        route_core_typed_pair_operator_plan_inventory_available = true,
        route_core_typed_pair_operator_plan_inventory_status =
            route_core_sidecar_inventory.crc_pair_operator_plan_inventory_status,
        route_core_typed_pair_operator_plan_blocker =
            route_core_sidecar_inventory.crc_pair_operator_plan_blocker,
        route_core_typed_pair_operator_plan_count =
            route_core_sidecar_inventory.crc_pair_operator_plan_count,
        route_core_typed_pair_operator_plan_blocked_count =
            route_core_sidecar_inventory.crc_pair_operator_plan_blocked_count,
        route_core_typed_pair_operator_plan_materialized =
            route_core_sidecar_inventory.crc_pair_operator_plan_materialized,
        route_core_typed_pair_operator_source_path_counts =
            CartesianRouteCore.pair_operator_source_path_counts(plan_inventory),
        route_core_typed_pair_operator_final_block_path_counts =
            CartesianRouteCore.pair_operator_final_block_path_counts(plan_inventory),
        route_core_typed_pair_operator_materialization_status_counts =
            CartesianRouteCore.pair_operator_materialization_status_counts(
                plan_inventory,
            ),
        route_core_typed_pair_operator_blocker_counts =
            CartesianRouteCore.pair_operator_blocker_counts(plan_inventory),
        route_core_typed_pair_operator_plan_family_counts =
            CartesianRouteCore.pair_operator_plan_family_counts(plan_inventory),
        route_core_typed_pair_operator_materialization_ready =
            route_core_sidecar_inventory.crc_pair_operator_materialization_ready,
        route_core_typed_pair_operator_materialization_readiness_status =
            route_core_sidecar_inventory.crc_pair_operator_materialization_readiness_status,
        route_core_typed_pair_operator_materialization_readiness_blocker =
            route_core_sidecar_inventory.crc_pair_operator_materialization_readiness_blocker,
        route_core_typed_pair_operator_materialization_readiness_requirements =
            route_core_sidecar_inventory.crc_pair_operator_materialization_readiness_requirements,
        route_core_typed_pair_operator_materialization_readiness_plan_count =
            route_core_sidecar_inventory.crc_pair_operator_materialization_readiness_plan_count,
        route_core_typed_pair_operator_materialization_readiness_blocked_count =
            route_core_sidecar_inventory.crc_pair_operator_materialization_readiness_blocked_count,
        route_core_typed_pair_operator_materialization_readiness_materialized_count =
            route_core_sidecar_inventory.crc_pair_operator_materialization_readiness_materialized_count,
    )
end

function _pqs_source_box_route_driver_route_core_readiness_requirements()
    return (
        :complete_crc_final_unit_inventory,
        :available_crc_pair_inventory,
        :positive_crc_pair_count,
        :crc_pair_order_matches_staged,
        :crc_pair_family_metadata_available,
    )
end

function _pqs_source_box_route_driver_route_core_pair_operator_readiness(;
    route_core_sidecar_inventory,
    route_core_pair_inventory_available::Bool,
    route_core_pair_count::Int,
    route_core_pair_order_matches_staged::Bool,
    route_core_pair_family_counts,
)
    requirements =
        _pqs_source_box_route_driver_route_core_readiness_requirements()
    if isnothing(route_core_sidecar_inventory) ||
       !hasproperty(route_core_sidecar_inventory, :status)
        return (;
            route_core_pair_operator_ready = false,
            route_core_pair_operator_readiness_status =
                :blocked_missing_route_core_sidecar_inventory,
            route_core_pair_operator_blocker =
                :missing_route_core_sidecar_inventory,
            route_core_pair_operator_readiness_requirements = requirements,
        )
    end
    if route_core_sidecar_inventory.status !=
       :available_route_core_sidecar_inventory ||
       route_core_sidecar_inventory.unsupported_unit_count != 0
        return (;
            route_core_pair_operator_ready = false,
            route_core_pair_operator_readiness_status =
                :blocked_incomplete_route_core_final_unit_inventory,
            route_core_pair_operator_blocker =
                :incomplete_route_core_final_unit_inventory,
            route_core_pair_operator_readiness_requirements = requirements,
        )
    end
    if !route_core_pair_inventory_available
        return (;
            route_core_pair_operator_ready = false,
            route_core_pair_operator_readiness_status =
                :blocked_missing_route_core_pair_inventory,
            route_core_pair_operator_blocker =
                :missing_route_core_pair_inventory,
            route_core_pair_operator_readiness_requirements = requirements,
        )
    end
    if route_core_pair_count <= 0
        return (;
            route_core_pair_operator_ready = false,
            route_core_pair_operator_readiness_status =
                :blocked_empty_route_core_pair_inventory,
            route_core_pair_operator_blocker =
                :empty_route_core_pair_inventory,
            route_core_pair_operator_readiness_requirements = requirements,
        )
    end
    if !route_core_pair_order_matches_staged
        return (;
            route_core_pair_operator_ready = false,
            route_core_pair_operator_readiness_status =
                :blocked_route_core_pair_order_mismatch,
            route_core_pair_operator_blocker = :route_core_pair_order_mismatch,
            route_core_pair_operator_readiness_requirements = requirements,
        )
    end
    if isempty(route_core_pair_family_counts)
        return (;
            route_core_pair_operator_ready = false,
            route_core_pair_operator_readiness_status =
                :blocked_missing_route_core_pair_family_metadata,
            route_core_pair_operator_blocker =
                :missing_route_core_pair_family_metadata,
            route_core_pair_operator_readiness_requirements = requirements,
        )
    end
    return (;
        route_core_pair_operator_ready = true,
        route_core_pair_operator_readiness_status =
            :ready_route_core_pair_operator_metadata,
        route_core_pair_operator_blocker = nothing,
        route_core_pair_operator_readiness_requirements = requirements,
    )
end

function _pqs_source_box_route_driver_route_core_pair_operator_preflight_metadata(
    route_core_pair_metadata,
)
    ready = route_core_pair_metadata.route_core_pair_operator_ready
    status =
        ready ?
        :ready_route_core_pair_operator_preflight :
        :blocked_route_core_pair_operator_preflight
    blocker =
        ready ? nothing : route_core_pair_metadata.route_core_pair_operator_blocker
    preflight = (;
        object_kind = :route_core_pair_operator_preflight,
        status,
        private_development_only = true,
        route_core_pair_operator_ready = ready,
        route_core_pair_operator_readiness_status =
            route_core_pair_metadata.route_core_pair_operator_readiness_status,
        route_core_pair_operator_blocker =
            route_core_pair_metadata.route_core_pair_operator_blocker,
        route_core_pair_operator_readiness_requirements =
            route_core_pair_metadata.route_core_pair_operator_readiness_requirements,
        route_core_final_unit_count =
            route_core_pair_metadata.route_core_final_unit_count,
        route_core_pair_count =
            route_core_pair_metadata.route_core_pair_count,
        route_core_pair_family_counts =
            route_core_pair_metadata.route_core_pair_family_counts,
        route_core_pair_order_matches_staged =
            route_core_pair_metadata.route_core_pair_order_matches_staged,
        route_core_pair_order_comparison_source =
            route_core_pair_metadata.route_core_pair_order_comparison_source,
        source_operator_blocks_materialized = false,
        operator_blocks_materialized = false,
        pair_operator_blocks_materialized = false,
        hamiltonian_matrices_materialized = false,
        artifact_export_materialized = false,
        product_doside_coefficient_maps_materialized = false,
        pqs_shell_realization_matrices_materialized = false,
    )
    return (;
        route_core_pair_operator_preflight_available = true,
        route_core_pair_operator_preflight_status = status,
        route_core_pair_operator_preflight = preflight,
        route_core_pair_operator_preflight_blocker = blocker,
    )
end

function _pqs_source_box_route_driver_route_core_operator_block_family_plan(
    preflight,
)
    return Tuple(
        (;
            pair_family = family.pair_family,
            pair_count = family.pair_count,
            operator_block_family =
                :metadata_only_route_core_pair_operator_family,
            operator_block_family_source =
                :crc_final_unit_lowering_recipe_pair_family,
        ) for family in preflight.route_core_pair_family_counts
    )
end

function _pqs_source_box_route_driver_route_core_pair_operator_plan_metadata(
    preflight_metadata,
)
    if !preflight_metadata.route_core_pair_operator_preflight_available
        return (;
            route_core_pair_operator_plan_available = false,
            route_core_pair_operator_plan_status = :not_available,
            route_core_pair_operator_plan = nothing,
            route_core_pair_operator_plan_blocker = :not_available,
        )
    end
    preflight = preflight_metadata.route_core_pair_operator_preflight
    ready = preflight.route_core_pair_operator_ready
    status =
        ready ?
        :ready_route_core_pair_operator_plan :
        :blocked_route_core_pair_operator_plan
    blocker =
        ready ? nothing : preflight_metadata.route_core_pair_operator_preflight_blocker
    plan = (;
        object_kind = :route_core_pair_operator_plan,
        status,
        private_development_only = true,
        preflight_status =
            preflight_metadata.route_core_pair_operator_preflight_status,
        preflight_blocker =
            preflight_metadata.route_core_pair_operator_preflight_blocker,
        route_core_pair_operator_ready =
            preflight.route_core_pair_operator_ready,
        route_core_pair_operator_readiness_status =
            preflight.route_core_pair_operator_readiness_status,
        route_core_final_unit_count =
            preflight.route_core_final_unit_count,
        route_core_pair_count =
            preflight.route_core_pair_count,
        route_core_pair_family_counts =
            preflight.route_core_pair_family_counts,
        route_core_pair_order_matches_staged =
            preflight.route_core_pair_order_matches_staged,
        route_core_pair_order_comparison_source =
            preflight.route_core_pair_order_comparison_source,
        operator_block_family_plan =
            _pqs_source_box_route_driver_route_core_operator_block_family_plan(
                preflight,
            ),
        operator_block_family_plan_source =
            :crc_pair_operator_preflight_family_counts,
        source_operator_blocks_materialized = false,
        pair_operator_blocks_materialized = false,
        operator_blocks_materialized = false,
        hamiltonian_matrices_materialized = false,
        artifact_export_materialized = false,
        product_doside_coefficient_maps_materialized = false,
        pqs_shell_realization_matrices_materialized = false,
    )
    return (;
        route_core_pair_operator_plan_available = true,
        route_core_pair_operator_plan_status = status,
        route_core_pair_operator_plan = plan,
        route_core_pair_operator_plan_blocker = blocker,
    )
end

function _pqs_source_box_route_driver_pair_keys_from_entries(pair_entries)
    return Tuple(pair.pair_key for pair in pair_entries)
end

function _pqs_source_box_route_driver_route_core_pair_comparison(
    route_core_sidecar_inventory,
    atom_growth_pair_inventory,
    route_skeleton,
)
    if !isnothing(atom_growth_pair_inventory) &&
       atom_growth_pair_inventory.status == :available_atom_growth_pair_inventory
        return (;
            pair_keys =
                _pqs_source_box_route_driver_pair_keys_from_entries(
                    atom_growth_pair_inventory.pair_entries,
                ),
            source = :atom_growth_pair_inventory,
        )
    end
    if !isnothing(route_core_sidecar_inventory) &&
       hasproperty(route_core_sidecar_inventory, :staged_pair_keys)
        return (;
            pair_keys = route_core_sidecar_inventory.staged_pair_keys,
            source = :route_core_sidecar_staged_pair_keys,
        )
    end
    if !isnothing(route_skeleton) && hasproperty(route_skeleton, :pair_entries)
        return (;
            pair_keys =
                _pqs_source_box_route_driver_pair_keys_from_entries(
                    route_skeleton.pair_entries,
                ),
            source = :route_skeleton_pair_entries,
        )
    end
    return (; pair_keys = (), source = :not_available)
end

function _pqs_source_box_route_driver_route_core_pair_stage_metadata(
    route_core_sidecar_inventory,
    atom_growth_pair_inventory,
    route_skeleton,
)
    if isnothing(route_core_sidecar_inventory) ||
       !hasproperty(route_core_sidecar_inventory, :crc_pair_inventory_available)
        return (;
            route_core_final_unit_count = 0,
            route_core_pair_inventory_available = false,
            route_core_pair_inventory_status =
                :blocked_missing_route_core_sidecar_inventory,
            route_core_pair_inventory = nothing,
            route_core_pair_count = 0,
            route_core_pair_keys = (),
            route_core_pair_order_matches_staged = false,
            route_core_pair_order_comparison_source = :not_available,
            route_core_pair_family_counts = (),
            route_core_pair_family_count_source = :not_available,
            route_core_summary_status =
                :blocked_missing_route_core_sidecar_inventory,
            route_core_pair_operator_ready = false,
            route_core_pair_operator_readiness_status =
                :blocked_missing_route_core_sidecar_inventory,
            route_core_pair_operator_blocker =
                :missing_route_core_sidecar_inventory,
            route_core_pair_operator_readiness_requirements =
                _pqs_source_box_route_driver_route_core_readiness_requirements(),
            _pqs_source_box_route_driver_route_core_typed_pair_operator_plan_unavailable_metadata(
                :blocked_missing_route_core_sidecar_inventory,
                :missing_route_core_sidecar_inventory,
            )...,
        )
    end

    route_core_typed_pair_operator_plan_metadata =
        _pqs_source_box_route_driver_route_core_typed_pair_operator_plan_metadata(
            route_core_sidecar_inventory,
        )
    comparison = _pqs_source_box_route_driver_route_core_pair_comparison(
        route_core_sidecar_inventory,
        atom_growth_pair_inventory,
        route_skeleton,
    )
    route_core_pair_inventory_available =
        route_core_sidecar_inventory.crc_pair_inventory_available
    route_core_pair_inventory =
        route_core_pair_inventory_available ?
        route_core_sidecar_inventory.crc_pair_inventory :
        nothing
    route_core_pair_keys =
        route_core_pair_inventory_available ?
        route_core_sidecar_inventory.crc_pair_keys :
        ()
    route_core_pair_order_matches_staged =
        route_core_pair_inventory_available &&
        !isempty(comparison.pair_keys) &&
        route_core_pair_keys == comparison.pair_keys
    if route_core_pair_inventory_available &&
       !isempty(comparison.pair_keys) &&
       !route_core_pair_order_matches_staged
        throw(ArgumentError("CRC pair keys do not match pair-stage comparison keys"))
    end
    route_core_pair_count =
        route_core_pair_inventory_available ?
        route_core_sidecar_inventory.crc_pair_count :
        0
    route_core_pair_family_counts =
        _pqs_source_box_route_driver_route_core_pair_family_counts(
            route_core_pair_inventory,
        )
    route_core_pair_operator_readiness =
        _pqs_source_box_route_driver_route_core_pair_operator_readiness(
            ;
            route_core_sidecar_inventory,
            route_core_pair_inventory_available,
            route_core_pair_count,
            route_core_pair_order_matches_staged,
            route_core_pair_family_counts,
        )

    return (;
        route_core_final_unit_count =
            route_core_pair_inventory_available ?
            route_core_sidecar_inventory.final_unit_count :
            0,
        route_core_pair_inventory_available,
        route_core_pair_inventory_status =
            route_core_sidecar_inventory.crc_pair_inventory_status,
        route_core_pair_inventory,
        route_core_pair_count,
        route_core_pair_keys,
        route_core_pair_order_matches_staged,
        route_core_pair_order_comparison_source = comparison.source,
        route_core_pair_family_counts,
        route_core_pair_family_count_source =
            route_core_pair_inventory_available ?
            :crc_final_unit_lowering_recipes :
            :not_available,
        route_core_summary_status =
            route_core_pair_inventory_available ?
            :available_route_core_unit_pair_summary :
            route_core_sidecar_inventory.crc_pair_inventory_status,
        route_core_pair_operator_readiness...,
        route_core_typed_pair_operator_plan_metadata...,
    )
end

function _pqs_source_box_route_driver_pair_stage_low_order_summary(
    units,
    transforms,
    route_skeleton,
)
    route_core_pair_unavailable = (;
        route_core_final_unit_count = 0,
        route_core_pair_inventory_available = false,
        route_core_pair_inventory_status = :not_available,
        route_core_pair_inventory = nothing,
        route_core_pair_count = 0,
        route_core_pair_keys = (),
        route_core_pair_order_matches_staged = false,
        route_core_pair_order_comparison_source = :not_available,
        route_core_pair_family_counts = (),
        route_core_pair_family_count_source = :not_available,
        route_core_summary_status = :not_available,
        route_core_pair_operator_ready = false,
        route_core_pair_operator_readiness_status = :not_available,
        route_core_pair_operator_blocker = :not_available,
        route_core_pair_operator_readiness_requirements =
            _pqs_source_box_route_driver_route_core_readiness_requirements(),
        route_core_pair_operator_preflight_available = false,
        route_core_pair_operator_preflight_status = :not_available,
        route_core_pair_operator_preflight = nothing,
        route_core_pair_operator_preflight_blocker = :not_available,
        route_core_pair_operator_plan_available = false,
        route_core_pair_operator_plan_status = :not_available,
        route_core_pair_operator_plan = nothing,
        route_core_pair_operator_plan_blocker = :not_available,
        _pqs_source_box_route_driver_route_core_typed_pair_operator_plan_unavailable_metadata(
            :not_available,
            :not_available,
        )...,
    )
    low_order_transforms =
        hasproperty(transforms, :low_order_transforms) ?
        transforms.low_order_transforms :
        nothing
    if isnothing(low_order_transforms)
        terminal_route_state =
            _pqs_source_box_route_driver_terminal_route_state_unavailable(
                :not_available_missing_transform_stage_summary,
                :missing_transform_stage_low_order_summary,
            )
        return (;
            object_kind = :cartesian_pair_stage_low_order_summary,
            status = :not_available_missing_transform_stage_summary,
            terminal_route_state,
            terminal_route_summary = terminal_route_state.summary,
            low_order_shellization_policy_requested = nothing,
            low_order_shellization_policy_resolved = :not_available,
            low_order_shellization_policy_source = :not_available,
            low_order_shellization_policy_status =
                :not_available_missing_transform_stage_summary,
            low_order_shellization_policy_blocker =
                :missing_transform_stage_low_order_summary,
            shellization_source = :not_available,
            shellization_kind = :not_available,
            unit_route_kind = :not_available,
            transform_route_kind = :not_available,
            pair_route_kind = :not_available,
            atom_growth_pairs_selected = false,
            terminal_shellification_pairs_selected = false,
            legacy_source_pairs_selected = false,
            terminal_shellification_pair_summary_available = false,
            terminal_shellification_scaffold_available = false,
            terminal_shellification_scaffold = nothing,
            terminal_shellification_region_count = 0,
            terminal_shellification_unit_inventory_available = false,
            terminal_shellification_unit_inventory = nothing,
            terminal_shellification_unit_count = 0,
            terminal_shellification_unit_keys = (),
            terminal_shellification_unit_roles = (),
            terminal_shellification_unit_kinds = (),
            terminal_shellification_unit_support_counts = (),
            terminal_shellification_lowering_contract_inventory_available = false,
            terminal_shellification_lowering_contract_inventory_status =
                :not_available,
            terminal_shellification_lowering_contract_inventory = nothing,
            terminal_shellification_lowering_contract_count = 0,
            terminal_shellification_lowering_contract_kinds = (),
            terminal_shellification_lowering_contract_kind_counts =
                (
                    direct_core_identity_cpb_count = 0,
                    direct_slab_identity_cpb_count = 0,
                    direct_boundary_slab_identity_cpb_count = 0,
                    white_lindsey_boundary_strata_count = 0,
                    pqs_filled_source_cpb_count = 0,
                    distorted_product_box_comx_count = 0,
                ),
            _pqs_source_box_route_driver_selected_terminal_lowering_fields(
                nothing,
                :not_available,
                nothing,
            )...,
            terminal_shellification_contract_counts_by_unit = (),
            terminal_shellification_lw_complete_shell_cpb_count = 0,
            terminal_shellification_lw_complete_shell_cpb_family_counts =
                (facet_cpb = 0, edge_cpb = 0, corner_cpb = 0),
            terminal_shellification_final_retained_unit_inventory_available = false,
            terminal_shellification_transform_contracts_available = false,
            terminal_shellification_pair_inventory_available = false,
            terminal_shellification_pair_inventory_status = :not_available,
            terminal_shellification_pair_materialization_status = :not_available,
            terminal_shellification_central_gap_region_count = 0,
            terminal_shellification_central_midpoint_slab_count = 0,
            terminal_shellification_central_distorted_product_box_count = 0,
            terminal_shellification_central_distorted_product_box_metadata = (),
            pair_operator_blocks_materialized = false,
            pair_inventory_known = false,
            pair_inventory_source = :not_available,
            pair_inventory = nothing,
            pair_entries = (),
            pair_count = 0,
            pair_family_counts = nothing,
            helper_by_pair_family = nothing,
            pair_operator_helper_by_family = nothing,
            pair_helper_status_by_family = nothing,
            operator_pairs_materialized = false,
            route_skeleton_pair_entry_count = 0,
            route_skeleton_pair_family_counts = nothing,
            route_skeleton_pair_entries = (),
            route_skeleton_helper_by_pair_family = nothing,
            independent_atom_growth_pair_inventory_available = false,
            plan_authority = false,
            active_source_authority = false,
            legacy_source_authority = false,
            pair_stage_fields_preserved = false,
            route_skeleton_pair_inventory_source = :not_available,
            route_core_pair_unavailable...,
            summary_only = true,
        )
    end

    atom_growth_pairs_selected =
        low_order_transforms.atom_growth_transforms_selected
    terminal_shellification_pairs_selected =
        low_order_transforms.terminal_shellification_transforms_selected
    legacy_source_pairs_selected =
        low_order_transforms.legacy_source_transforms_selected
    pair_route_kind =
        atom_growth_pairs_selected ?
        :atom_growth_complete_rectangular_low_order_pairs :
        terminal_shellification_pairs_selected ?
        :terminal_shellification_low_order_pairs :
        legacy_source_pairs_selected ?
        :legacy_diatomic_source_low_order_pairs :
        :not_selected
    terminal_shellification_scaffold_available =
        terminal_shellification_pairs_selected &&
        low_order_transforms.terminal_shellification_scaffold_available
    terminal_shellification_pair_summary_available =
        terminal_shellification_pairs_selected
    terminal_shellification_pair_inventory_status =
        terminal_shellification_pairs_selected ?
        :deferred_terminal_shellification_pair_inventory :
        :not_selected
    terminal_shellification_pair_materialization_status =
        terminal_shellification_pairs_selected ?
        :deferred_terminal_shellification_pair_materialization :
        :not_selected
    atom_growth_plan_unit_inventory =
        atom_growth_pairs_selected && hasproperty(units, :plan_unit_inventory) ?
        units.plan_unit_inventory :
        nothing
    atom_growth_pair_inventory =
        atom_growth_pairs_selected ?
        _pqs_source_box_route_driver_atom_growth_pair_inventory(
            atom_growth_plan_unit_inventory,
        ) :
        nothing
    route_core_pair_metadata =
        atom_growth_pairs_selected ?
        _pqs_source_box_route_driver_route_core_pair_stage_metadata(
            hasproperty(units, :route_core_sidecar_inventory) ?
            units.route_core_sidecar_inventory :
            nothing,
            atom_growth_pair_inventory,
            route_skeleton,
        ) :
        terminal_shellification_pairs_selected ?
        merge(
            route_core_pair_unavailable,
            (;
                route_core_pair_inventory_status =
                    :deferred_terminal_shellification_pair_inventory,
                route_core_summary_status =
                    :deferred_terminal_shellification_pair_inventory,
                route_core_pair_operator_ready = false,
                route_core_pair_operator_readiness_status =
                    :deferred_terminal_shellification_pair_inventory,
                route_core_pair_operator_blocker =
                    :deferred_terminal_shellification_pair_inventory,
                _pqs_source_box_route_driver_route_core_typed_pair_operator_plan_unavailable_metadata(
                    :deferred_terminal_shellification_typed_pair_operator_plan_inventory,
                    :deferred_terminal_shellification_pair_inventory,
                )...,
            ),
        ) :
        merge(
            route_core_pair_unavailable,
            (;
                route_core_pair_inventory_status =
                    :not_selected_legacy_source_pairs,
                route_core_summary_status =
                    :not_selected_legacy_source_pairs,
                route_core_pair_operator_ready = false,
                route_core_pair_operator_readiness_status =
                    :not_selected_legacy_source_pairs,
                route_core_pair_operator_blocker =
                    :not_selected_legacy_source_pairs,
                route_core_typed_pair_operator_plan_inventory_status =
                    :not_selected_legacy_source_pairs,
                route_core_typed_pair_operator_plan_blocker =
                    :not_selected_legacy_source_pairs,
            ),
        )
    route_core_pair_operator_preflight_metadata =
        _pqs_source_box_route_driver_route_core_pair_operator_preflight_metadata(
            route_core_pair_metadata,
        )
    route_core_pair_operator_plan_metadata =
        _pqs_source_box_route_driver_route_core_pair_operator_plan_metadata(
            route_core_pair_operator_preflight_metadata,
        )
    independent_atom_growth_pair_inventory_available =
        !isnothing(atom_growth_pair_inventory) &&
        atom_growth_pair_inventory.status == :available_atom_growth_pair_inventory
    pair_entries =
        independent_atom_growth_pair_inventory_available ?
        atom_growth_pair_inventory.pair_entries :
        atom_growth_pairs_selected ?
        () :
        terminal_shellification_pairs_selected ?
        () :
        route_skeleton.pair_entries
    pair_family_counts =
        independent_atom_growth_pair_inventory_available ?
        atom_growth_pair_inventory.pair_family_counts :
        atom_growth_pairs_selected ?
        (white_lindsey_low_order_atom_growth_unit_pair = 0,) :
        terminal_shellification_pairs_selected ?
        () :
        route_skeleton.pair_family_counts
    pair_inventory_source =
        independent_atom_growth_pair_inventory_available ?
        atom_growth_pair_inventory.pair_inventory_source :
        atom_growth_pairs_selected ?
        atom_growth_pair_inventory.pair_inventory_source :
        terminal_shellification_pairs_selected ?
        :terminal_shellification_scaffold :
        :route_skeleton_pair_entries_only
    pair_inventory_known =
        atom_growth_pairs_selected ?
        independent_atom_growth_pair_inventory_available :
        terminal_shellification_pairs_selected ?
        false :
        !isempty(pair_entries)
    atom_growth_pair_helper_status_by_family = (
        white_lindsey_low_order_atom_growth_unit_pair =
            :deferred_no_pair_operator_block_helper,
    )
    pair_operator_helper_by_family =
        atom_growth_pairs_selected ?
        atom_growth_pair_helper_status_by_family :
        terminal_shellification_pairs_selected ?
        () :
        route_skeleton.helper_by_pair_family
    pair_helper_status_by_family =
        atom_growth_pairs_selected ?
        atom_growth_pair_helper_status_by_family :
        terminal_shellification_pairs_selected ?
        () :
        nothing
    operator_pairs_materialized =
        independent_atom_growth_pair_inventory_available ?
        atom_growth_pair_inventory.operator_pairs_materialized :
        false
    terminal_route_state = low_order_transforms.terminal_route_state
    terminal_route_summary = low_order_transforms.terminal_route_summary

    return (;
        object_kind = :cartesian_pair_stage_low_order_summary,
        status =
            terminal_shellification_pairs_selected &&
            terminal_shellification_scaffold_available ?
            :deferred_terminal_shellification_pair_inventory :
            low_order_transforms.status ==
            :available_transform_stage_low_order_summary ?
            :available_pair_stage_low_order_summary :
            low_order_transforms.status,
        terminal_route_state,
        terminal_route_summary,
        low_order_shellization_policy_requested =
            low_order_transforms.low_order_shellization_policy_requested,
        low_order_shellization_policy_resolved =
            low_order_transforms.low_order_shellization_policy_resolved,
        low_order_shellization_policy_source =
            low_order_transforms.low_order_shellization_policy_source,
        low_order_shellization_policy_status =
            low_order_transforms.low_order_shellization_policy_status,
        low_order_shellization_policy_blocker =
            low_order_transforms.low_order_shellization_policy_blocker,
        shellization_source = low_order_transforms.shellization_source,
        shellization_kind = low_order_transforms.shellization_kind,
        unit_route_kind = low_order_transforms.unit_route_kind,
        transform_route_kind = low_order_transforms.transform_route_kind,
        pair_route_kind,
        atom_growth_pairs_selected,
        terminal_shellification_pairs_selected,
        legacy_source_pairs_selected,
        terminal_shellification_pair_summary_available,
        _pqs_source_box_route_driver_terminal_shellification_alias_fields(
            low_order_transforms,
            terminal_shellification_pairs_selected;
            include_crc_sidecar_summary = false,
        )...,
        terminal_shellification_transform_contracts_available =
            terminal_shellification_pairs_selected &&
            low_order_transforms.terminal_shellification_transform_contracts_available,
        terminal_shellification_pair_inventory_available = false,
        terminal_shellification_pair_inventory_status,
        terminal_shellification_pair_materialization_status,
        pair_operator_blocks_materialized = false,
        operator_pairs_materialized,
        pair_inventory_known,
        pair_inventory_source,
        pair_inventory = atom_growth_pair_inventory,
        pair_entries,
        pair_count = length(pair_entries),
        pair_family_counts,
        helper_by_pair_family = pair_operator_helper_by_family,
        pair_operator_helper_by_family,
        pair_helper_status_by_family,
        route_skeleton_pair_entry_count = length(route_skeleton.pair_entries),
        route_skeleton_pair_family_counts = route_skeleton.pair_family_counts,
        route_skeleton_pair_entries = route_skeleton.pair_entries,
        route_skeleton_helper_by_pair_family = route_skeleton.helper_by_pair_family,
        independent_atom_growth_pair_inventory_available,
        plan_authority = low_order_transforms.plan_authority,
        active_source_authority = low_order_transforms.active_source_authority,
        legacy_source_authority = low_order_transforms.legacy_source_authority,
        pair_stage_fields_preserved = true,
        route_skeleton_pair_inventory_source =
            :route_skeleton_compatibility_fields,
        route_core_pair_metadata...,
        route_core_pair_operator_preflight_metadata...,
        route_core_pair_operator_plan_metadata...,
        summary_only =
            terminal_shellification_pairs_selected ||
            !independent_atom_growth_pair_inventory_available,
    )
end

function cartesian_pair_terms(units, transforms, recipe)
    route_skeleton = units.route_skeleton
    low_order_pairs =
        _pqs_source_box_route_driver_pair_stage_low_order_summary(
            units,
            transforms,
            route_skeleton,
        )

    return (;
        object_kind = :cartesian_pair_terms,
        status = units.status,
        route_family = recipe.route_family,
        retained_dimension = transforms.retained_dimension,
        low_order_pairs,
        terminal_route_state = low_order_pairs.terminal_route_state,
        terminal_route_summary = low_order_pairs.terminal_route_summary,
        low_order_pair_route_kind = low_order_pairs.pair_route_kind,
        atom_growth_pairs_selected = low_order_pairs.atom_growth_pairs_selected,
        pair_operator_blocks_materialized =
            low_order_pairs.pair_operator_blocks_materialized,
        operator_pairs_materialized = low_order_pairs.operator_pairs_materialized,
        independent_atom_growth_pair_inventory_available =
            low_order_pairs.independent_atom_growth_pair_inventory_available,
        pair_inventory = low_order_pairs.pair_inventory,
        pair_inventory_source = low_order_pairs.pair_inventory_source,
        route_core_final_unit_count =
            low_order_pairs.route_core_final_unit_count,
        route_core_pair_inventory_available =
            low_order_pairs.route_core_pair_inventory_available,
        route_core_pair_inventory_status =
            low_order_pairs.route_core_pair_inventory_status,
        route_core_pair_count = low_order_pairs.route_core_pair_count,
        route_core_pair_order_matches_staged =
            low_order_pairs.route_core_pair_order_matches_staged,
        route_core_pair_order_comparison_source =
            low_order_pairs.route_core_pair_order_comparison_source,
        route_core_pair_family_counts =
            low_order_pairs.route_core_pair_family_counts,
        route_core_summary_status =
            low_order_pairs.route_core_summary_status,
        active_source_authority = low_order_pairs.active_source_authority,
        pair_entries = low_order_pairs.pair_entries,
        pair_family_counts = low_order_pairs.pair_family_counts,
        helper_by_pair_family = low_order_pairs.helper_by_pair_family,
        pair_operator_helper_by_family =
            low_order_pairs.pair_operator_helper_by_family,
        pair_helper_status_by_family = low_order_pairs.pair_helper_status_by_family,
        route_skeleton_pair_entries = route_skeleton.pair_entries,
        route_skeleton_pair_family_counts = route_skeleton.pair_family_counts,
        route_skeleton_helper_by_pair_family = route_skeleton.helper_by_pair_family,
        pair_stage = :pair_operator_terms_described,
    )
end

function _pqs_source_box_route_driver_assembly_stage_low_order_summary(pairs)
    low_order_pairs =
        hasproperty(pairs, :low_order_pairs) ? pairs.low_order_pairs : nothing
    if isnothing(low_order_pairs)
        terminal_route_state =
            _pqs_source_box_route_driver_terminal_route_state_unavailable(
                :not_available_missing_pair_stage_summary,
                :missing_pair_stage_low_order_summary,
            )
        return (;
            object_kind = :cartesian_assembly_stage_low_order_summary,
            status = :not_available_missing_pair_stage_summary,
            terminal_route_state,
            terminal_route_summary = terminal_route_state.summary,
            low_order_shellization_policy_requested = nothing,
            low_order_shellization_policy_resolved = :not_available,
            low_order_shellization_policy_source = :not_available,
            low_order_shellization_policy_status =
                :not_available_missing_pair_stage_summary,
            low_order_shellization_policy_blocker =
                :missing_pair_stage_low_order_summary,
            shellization_source = :not_available,
            shellization_kind = :not_available,
            unit_route_kind = :not_available,
            transform_route_kind = :not_available,
            pair_route_kind = :not_available,
            assembly_source = :not_available,
            assembly_route_kind = :not_available,
            assembly_kind = :not_available,
            atom_growth_assembly_selected = false,
            terminal_shellification_assembly_selected = false,
            legacy_source_assembly_selected = false,
            terminal_shellification_assembly_summary_available = false,
            terminal_shellification_scaffold_available = false,
            terminal_shellification_scaffold = nothing,
            terminal_shellification_region_count = 0,
            terminal_shellification_unit_inventory_available = false,
            terminal_shellification_unit_inventory = nothing,
            terminal_shellification_unit_count = 0,
            terminal_shellification_unit_keys = (),
            terminal_shellification_unit_roles = (),
            terminal_shellification_unit_kinds = (),
            terminal_shellification_unit_support_counts = (),
            terminal_shellification_lowering_contract_inventory_available = false,
            terminal_shellification_lowering_contract_inventory_status =
                :not_available,
            terminal_shellification_lowering_contract_inventory = nothing,
            terminal_shellification_lowering_contract_count = 0,
            terminal_shellification_lowering_contract_kinds = (),
            terminal_shellification_lowering_contract_kind_counts =
                (
                    direct_core_identity_cpb_count = 0,
                    direct_slab_identity_cpb_count = 0,
                    direct_boundary_slab_identity_cpb_count = 0,
                    white_lindsey_boundary_strata_count = 0,
                    pqs_filled_source_cpb_count = 0,
                    distorted_product_box_comx_count = 0,
                ),
            _pqs_source_box_route_driver_selected_terminal_lowering_fields(
                nothing,
                :not_available,
                nothing,
            )...,
            terminal_shellification_contract_counts_by_unit = (),
            terminal_shellification_lw_complete_shell_cpb_count = 0,
            terminal_shellification_lw_complete_shell_cpb_family_counts =
                (facet_cpb = 0, edge_cpb = 0, corner_cpb = 0),
            terminal_shellification_final_retained_unit_inventory_available = false,
            terminal_shellification_transform_contracts_available = false,
            terminal_shellification_pair_inventory_available = false,
            terminal_shellification_pair_inventory_status = :not_available,
            terminal_shellification_pair_materialization_status = :not_available,
            terminal_shellification_assembly_materialization_status =
                :not_available,
            terminal_shellification_central_gap_region_count = 0,
            terminal_shellification_central_midpoint_slab_count = 0,
            terminal_shellification_central_distorted_product_box_count = 0,
            terminal_shellification_central_distorted_product_box_metadata = (),
            hamiltonian_matrices_materialized = false,
            operator_matrices_materialized = false,
            pair_operator_blocks_materialized = false,
            pair_operator_blocks_available = false,
            pair_inventory_source = :not_available,
            pair_inventory_known = false,
            independent_atom_growth_pair_inventory_available = false,
            pair_count = 0,
            pair_family_counts = nothing,
            route_core_final_unit_count = 0,
            route_core_pair_inventory_available = false,
            route_core_pair_inventory_status = :not_available,
            route_core_pair_count = 0,
            route_core_pair_order_matches_staged = false,
            route_core_pair_order_comparison_source = :not_available,
            route_core_pair_family_counts = (),
            route_core_summary_status = :not_available,
            route_core_pair_operator_ready = false,
            route_core_pair_operator_readiness_status = :not_available,
            route_core_pair_operator_blocker = :not_available,
            route_core_pair_operator_readiness_requirements =
                _pqs_source_box_route_driver_route_core_readiness_requirements(),
            route_core_pair_operator_preflight_available = false,
            route_core_pair_operator_preflight_status = :not_available,
            route_core_pair_operator_preflight = nothing,
            route_core_pair_operator_preflight_blocker = :not_available,
            route_core_pair_operator_plan_available = false,
            route_core_pair_operator_plan_status = :not_available,
            route_core_pair_operator_plan = nothing,
            route_core_pair_operator_plan_blocker = :not_available,
            _pqs_source_box_route_driver_route_core_typed_pair_operator_plan_unavailable_metadata(
                :not_available,
                :not_available,
            )...,
            helper_by_pair_family = nothing,
            pair_operator_helper_by_family = nothing,
            pair_helper_status_by_family = nothing,
            assembly_can_proceed_from_current_staged_data = false,
            assembly_requires_materialization = true,
            assembly_materialization_status =
                :blocked_missing_pair_stage_summary,
            assembly_blocker = :missing_pair_stage_low_order_summary,
            plan_authority = false,
            active_source_authority = false,
            legacy_source_authority = false,
            assembly_stage_fields_preserved = false,
            summary_only = true,
        )
    end

    atom_growth_assembly_selected =
        low_order_pairs.atom_growth_pairs_selected
    terminal_shellification_assembly_selected =
        low_order_pairs.terminal_shellification_pairs_selected
    legacy_source_assembly_selected =
        low_order_pairs.legacy_source_pairs_selected
    assembly_route_kind =
        atom_growth_assembly_selected ?
        :atom_growth_complete_rectangular_low_order_assembly :
        terminal_shellification_assembly_selected ?
        :terminal_shellification_low_order_assembly :
        legacy_source_assembly_selected ?
        :legacy_diatomic_source_low_order_assembly :
        :not_selected
    assembly_source =
        atom_growth_assembly_selected ?
        :atom_growth_complete_rectangular_low_order_pair_terms :
        terminal_shellification_assembly_selected ?
        :terminal_shellification_pair_terms :
        legacy_source_assembly_selected ?
        :legacy_diatomic_source_pair_terms :
        :not_selected
    assembly_kind =
        atom_growth_assembly_selected ?
        :atom_growth_complete_rectangular_low_order :
        terminal_shellification_assembly_selected ?
        :terminal_shellification_low_order :
        legacy_source_assembly_selected ?
        :legacy_diatomic_source_low_order :
        :not_selected
    terminal_shellification_scaffold_available =
        terminal_shellification_assembly_selected &&
        low_order_pairs.terminal_shellification_scaffold_available
    terminal_shellification_assembly_summary_available =
        terminal_shellification_assembly_selected
    terminal_shellification_assembly_materialization_status =
        terminal_shellification_assembly_selected ?
        :deferred_terminal_shellification_assembly_materialization :
        :not_selected
    pair_operator_blocks_materialized =
        low_order_pairs.pair_operator_blocks_materialized
    assembly_can_proceed_from_current_staged_data =
        pair_operator_blocks_materialized
    assembly_requires_materialization =
        !assembly_can_proceed_from_current_staged_data
    assembly_materialization_status =
        pair_operator_blocks_materialized ?
        :ready_for_low_order_operator_matrix_assembly :
        atom_growth_assembly_selected ?
        :deferred_atom_growth_complete_rectangular_pair_block_materialization :
        terminal_shellification_assembly_selected ?
        terminal_shellification_assembly_materialization_status :
        legacy_source_assembly_selected ?
        :deferred_legacy_diatomic_source_pair_block_materialization :
        :not_selected
    assembly_blocker =
        assembly_requires_materialization ?
        (
            terminal_shellification_assembly_selected ?
            :terminal_shellification_pair_blocks_deferred :
            :pair_operator_blocks_deferred
        ) :
        nothing
    terminal_route_state = low_order_pairs.terminal_route_state
    terminal_route_summary = low_order_pairs.terminal_route_summary

    return (;
        object_kind = :cartesian_assembly_stage_low_order_summary,
        status =
            terminal_shellification_assembly_selected &&
            terminal_shellification_scaffold_available ?
            :deferred_terminal_shellification_assembly_materialization :
            low_order_pairs.status == :available_pair_stage_low_order_summary ?
            :available_assembly_stage_low_order_summary :
            low_order_pairs.status,
        terminal_route_state,
        terminal_route_summary,
        low_order_shellization_policy_requested =
            low_order_pairs.low_order_shellization_policy_requested,
        low_order_shellization_policy_resolved =
            low_order_pairs.low_order_shellization_policy_resolved,
        low_order_shellization_policy_source =
            low_order_pairs.low_order_shellization_policy_source,
        low_order_shellization_policy_status =
            low_order_pairs.low_order_shellization_policy_status,
        low_order_shellization_policy_blocker =
            low_order_pairs.low_order_shellization_policy_blocker,
        shellization_source = low_order_pairs.shellization_source,
        shellization_kind = low_order_pairs.shellization_kind,
        unit_route_kind = low_order_pairs.unit_route_kind,
        transform_route_kind = low_order_pairs.transform_route_kind,
        pair_route_kind = low_order_pairs.pair_route_kind,
        assembly_source,
        assembly_route_kind,
        assembly_kind,
        atom_growth_assembly_selected,
        terminal_shellification_assembly_selected,
        legacy_source_assembly_selected,
        terminal_shellification_assembly_summary_available,
        _pqs_source_box_route_driver_terminal_shellification_alias_fields(
            low_order_pairs,
            terminal_shellification_assembly_selected;
            include_crc_sidecar_summary = false,
        )...,
        terminal_shellification_transform_contracts_available =
            terminal_shellification_assembly_selected &&
            low_order_pairs.terminal_shellification_transform_contracts_available,
        terminal_shellification_pair_inventory_available =
            terminal_shellification_assembly_selected &&
            low_order_pairs.terminal_shellification_pair_inventory_available,
        terminal_shellification_pair_inventory_status =
            terminal_shellification_assembly_selected ?
            low_order_pairs.terminal_shellification_pair_inventory_status :
            :not_selected,
        terminal_shellification_pair_materialization_status =
            terminal_shellification_assembly_selected ?
            low_order_pairs.terminal_shellification_pair_materialization_status :
            :not_selected,
        terminal_shellification_assembly_materialization_status,
        hamiltonian_matrices_materialized = false,
        operator_matrices_materialized = false,
        pair_operator_blocks_materialized,
        pair_operator_blocks_available = pair_operator_blocks_materialized,
        pair_inventory_source = low_order_pairs.pair_inventory_source,
        pair_inventory_known = low_order_pairs.pair_inventory_known,
        independent_atom_growth_pair_inventory_available =
            low_order_pairs.independent_atom_growth_pair_inventory_available,
        pair_count = low_order_pairs.pair_count,
        pair_family_counts = low_order_pairs.pair_family_counts,
        route_core_final_unit_count =
            low_order_pairs.route_core_final_unit_count,
        route_core_pair_inventory_available =
            low_order_pairs.route_core_pair_inventory_available,
        route_core_pair_inventory_status =
            low_order_pairs.route_core_pair_inventory_status,
        route_core_pair_count = low_order_pairs.route_core_pair_count,
        route_core_pair_order_matches_staged =
            low_order_pairs.route_core_pair_order_matches_staged,
        route_core_pair_order_comparison_source =
            low_order_pairs.route_core_pair_order_comparison_source,
        route_core_pair_family_counts =
            low_order_pairs.route_core_pair_family_counts,
        route_core_summary_status = low_order_pairs.route_core_summary_status,
        route_core_pair_operator_ready =
            low_order_pairs.route_core_pair_operator_ready,
        route_core_pair_operator_readiness_status =
            low_order_pairs.route_core_pair_operator_readiness_status,
        route_core_pair_operator_blocker =
            low_order_pairs.route_core_pair_operator_blocker,
        route_core_pair_operator_readiness_requirements =
            low_order_pairs.route_core_pair_operator_readiness_requirements,
        route_core_pair_operator_preflight_available =
            low_order_pairs.route_core_pair_operator_preflight_available,
        route_core_pair_operator_preflight_status =
            low_order_pairs.route_core_pair_operator_preflight_status,
        route_core_pair_operator_preflight =
            low_order_pairs.route_core_pair_operator_preflight,
        route_core_pair_operator_preflight_blocker =
            low_order_pairs.route_core_pair_operator_preflight_blocker,
        route_core_pair_operator_plan_available =
            low_order_pairs.route_core_pair_operator_plan_available,
        route_core_pair_operator_plan_status =
            low_order_pairs.route_core_pair_operator_plan_status,
        route_core_pair_operator_plan =
            low_order_pairs.route_core_pair_operator_plan,
        route_core_pair_operator_plan_blocker =
            low_order_pairs.route_core_pair_operator_plan_blocker,
        route_core_typed_pair_operator_plan_inventory_available =
            low_order_pairs.route_core_typed_pair_operator_plan_inventory_available,
        route_core_typed_pair_operator_plan_inventory_status =
            low_order_pairs.route_core_typed_pair_operator_plan_inventory_status,
        route_core_typed_pair_operator_plan_blocker =
            low_order_pairs.route_core_typed_pair_operator_plan_blocker,
        route_core_typed_pair_operator_plan_count =
            low_order_pairs.route_core_typed_pair_operator_plan_count,
        route_core_typed_pair_operator_plan_blocked_count =
            low_order_pairs.route_core_typed_pair_operator_plan_blocked_count,
        route_core_typed_pair_operator_plan_materialized =
            low_order_pairs.route_core_typed_pair_operator_plan_materialized,
        route_core_typed_pair_operator_source_path_counts =
            low_order_pairs.route_core_typed_pair_operator_source_path_counts,
        route_core_typed_pair_operator_final_block_path_counts =
            low_order_pairs.route_core_typed_pair_operator_final_block_path_counts,
        route_core_typed_pair_operator_materialization_status_counts =
            low_order_pairs.route_core_typed_pair_operator_materialization_status_counts,
        route_core_typed_pair_operator_blocker_counts =
            low_order_pairs.route_core_typed_pair_operator_blocker_counts,
        route_core_typed_pair_operator_plan_family_counts =
            low_order_pairs.route_core_typed_pair_operator_plan_family_counts,
        route_core_typed_pair_operator_materialization_ready =
            low_order_pairs.route_core_typed_pair_operator_materialization_ready,
        route_core_typed_pair_operator_materialization_readiness_status =
            low_order_pairs.route_core_typed_pair_operator_materialization_readiness_status,
        route_core_typed_pair_operator_materialization_readiness_blocker =
            low_order_pairs.route_core_typed_pair_operator_materialization_readiness_blocker,
        route_core_typed_pair_operator_materialization_readiness_requirements =
            low_order_pairs.route_core_typed_pair_operator_materialization_readiness_requirements,
        route_core_typed_pair_operator_materialization_readiness_plan_count =
            low_order_pairs.route_core_typed_pair_operator_materialization_readiness_plan_count,
        route_core_typed_pair_operator_materialization_readiness_blocked_count =
            low_order_pairs.route_core_typed_pair_operator_materialization_readiness_blocked_count,
        route_core_typed_pair_operator_materialization_readiness_materialized_count =
            low_order_pairs.route_core_typed_pair_operator_materialization_readiness_materialized_count,
        helper_by_pair_family = low_order_pairs.helper_by_pair_family,
        pair_operator_helper_by_family =
            low_order_pairs.pair_operator_helper_by_family,
        pair_helper_status_by_family = low_order_pairs.pair_helper_status_by_family,
        assembly_can_proceed_from_current_staged_data,
        assembly_requires_materialization,
        assembly_materialization_status,
        assembly_blocker,
        plan_authority = low_order_pairs.plan_authority,
        active_source_authority = low_order_pairs.active_source_authority,
        legacy_source_authority = low_order_pairs.legacy_source_authority,
        assembly_stage_fields_preserved = true,
        summary_only = true,
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
    diatomic_physical_gausslet_rhf_input_contract =
        _pqs_source_box_route_driver_diatomic_physical_gausslet_rhf_input_contract(
            parent,
            route_skeleton,
            recipe,
            diatomic_physical_gausslet_source_plan_payload,
            diatomic_physical_gausslet_final_basis_payload,
            diatomic_physical_gausslet_h1_payload,
            diatomic_physical_gausslet_h1_j_payload,
        )
    diatomic_physical_gausslet_rhf_execution_payload =
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
        object_kind = :cartesian_assembly,
        status = route_skeleton.status,
        spacing_inputs = shells.spacing_inputs,
        route_inputs = units.route_inputs,
        route_skeleton,
        raw_box = units.raw_box,
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
        complete_core_shell_h1_j_diagnostic_summary =
            complete_core_shell_h1_j_diagnostic_payload.summary,
        complete_core_shell_h1_j_diagnostic_status =
            complete_core_shell_h1_j_diagnostic_payload.status,
        complete_core_shell_h1_j_diagnostic_blocker =
            complete_core_shell_h1_j_diagnostic_payload.blocker,
        terminal_route_state = low_order_assembly.terminal_route_state,
        terminal_route_summary = low_order_assembly.terminal_route_summary,
        low_order_assembly_source = low_order_assembly.assembly_source,
        low_order_assembly_route_kind =
            low_order_assembly.assembly_route_kind,
        atom_growth_assembly_selected =
            low_order_assembly.atom_growth_assembly_selected,
        hamiltonian_matrices_materialized =
            low_order_assembly.hamiltonian_matrices_materialized,
        operator_matrices_materialized =
            low_order_assembly.operator_matrices_materialized,
        pair_operator_blocks_materialized =
            low_order_assembly.pair_operator_blocks_materialized,
        low_order_pair_inventory_source =
            low_order_assembly.pair_inventory_source,
        low_order_pair_inventory_known =
            low_order_assembly.pair_inventory_known,
        low_order_independent_atom_growth_pair_inventory_available =
            low_order_assembly.independent_atom_growth_pair_inventory_available,
        low_order_pair_count = low_order_assembly.pair_count,
        low_order_pair_family_counts = low_order_assembly.pair_family_counts,
        low_order_route_core_final_unit_count =
            low_order_assembly.route_core_final_unit_count,
        low_order_route_core_pair_inventory_available =
            low_order_assembly.route_core_pair_inventory_available,
        low_order_route_core_pair_inventory_status =
            low_order_assembly.route_core_pair_inventory_status,
        low_order_route_core_pair_count =
            low_order_assembly.route_core_pair_count,
        low_order_route_core_pair_order_matches_staged =
            low_order_assembly.route_core_pair_order_matches_staged,
        low_order_route_core_pair_order_comparison_source =
            low_order_assembly.route_core_pair_order_comparison_source,
        low_order_route_core_pair_family_counts =
            low_order_assembly.route_core_pair_family_counts,
        low_order_route_core_summary_status =
            low_order_assembly.route_core_summary_status,
        low_order_pair_operator_helper_by_family =
            low_order_assembly.pair_operator_helper_by_family,
        low_order_pair_helper_status_by_family =
            low_order_assembly.pair_helper_status_by_family,
        assembly_requires_materialization =
            low_order_assembly.assembly_requires_materialization,
        active_source_authority = low_order_assembly.active_source_authority,
        assembly_stage = :assembled_report_inputs,
    )
end

function _pqs_source_box_route_driver_report_stage_pqs_prototype_summary(
    assembly,
)
    units = hasproperty(assembly, :units) ? assembly.units : nothing
    transforms =
        hasproperty(assembly, :transforms) ? assembly.transforms : nothing
    lowering_available =
        !isnothing(units) &&
        hasproperty(units, :pqs_lowering_prototype_available) &&
        units.pqs_lowering_prototype_available
    transform_available =
        !isnothing(transforms) &&
        hasproperty(transforms, :pqs_transform_prototype_available) &&
        transforms.pqs_transform_prototype_available
    lowering_prototype =
        lowering_available ? units.pqs_lowering_prototype : nothing
    transform_prototype =
        transform_available ? transforms.pqs_transform_prototype : nothing
    prototype =
        transform_available ?
        transform_prototype :
        lowering_available ?
        lowering_prototype :
        nothing
    source_cpb =
        isnothing(prototype) || !hasproperty(prototype, :source_cpb) ?
        nothing :
        prototype.source_cpb
    source_cpb_support_count =
        isnothing(prototype) ||
        !hasproperty(prototype, :source_cpb_support_count) ?
        nothing :
        prototype.source_cpb_support_count
    owned_support_count =
        isnothing(prototype) ||
        !hasproperty(prototype, :owned_support_count) ?
        nothing :
        prototype.owned_support_count
    intermediate_retained_space =
        isnothing(prototype) ||
        !hasproperty(prototype, :intermediate_retained_space) ?
        nothing :
        prototype.intermediate_retained_space

    return (;
        pqs_lowering_prototype_available = lowering_available,
        pqs_transform_prototype_available = transform_available,
        pqs_lowering_prototype = lowering_prototype,
        pqs_transform_prototype = transform_prototype,
        pqs_prototype_unit_key =
            isnothing(prototype) ? nothing : prototype.unit_key,
        pqs_prototype_stage =
            isnothing(prototype) ? :not_available : :metadata_only,
        pqs_prototype_source =
            transform_available ?
            :transform_stage_pqs_transform_prototype :
            lowering_available ?
            :unit_stage_pqs_lowering_prototype :
            :not_available,
        pqs_prototype_source_cpb_kind =
            isnothing(source_cpb) ? nothing : source_cpb.cpb_family,
        pqs_prototype_source_cpb_support_count =
            source_cpb_support_count,
        pqs_prototype_source_cpb_support_count_source =
            isnothing(prototype) ||
            !hasproperty(prototype, :source_cpb_support_count_source) ?
            nothing :
            prototype.source_cpb_support_count_source,
        pqs_prototype_owned_support_count = owned_support_count,
        pqs_prototype_owned_support_count_source =
            isnothing(prototype) ||
            !hasproperty(prototype, :owned_support_count_source) ?
            nothing :
            prototype.owned_support_count_source,
        pqs_prototype_source_count_distinct_from_owned_support_count =
            !isnothing(source_cpb_support_count) &&
            !isnothing(owned_support_count) &&
            source_cpb_support_count != owned_support_count,
        pqs_prototype_owned_support_is_cpb =
            isnothing(prototype) ||
            !hasproperty(prototype, :owned_support) ?
            false :
            prototype.owned_support.owned_support_is_cpb,
        pqs_prototype_intermediate_retained_space =
            isnothing(intermediate_retained_space) ?
            nothing :
            intermediate_retained_space.retained_rule,
        pqs_prototype_shell_realization =
            transform_available ?
            :shell_projection_lowdin_deferred :
            lowering_available ?
            :shell_projection_lowdin_deferred :
            nothing,
        pqs_prototype_coefficient_maps_materialized = false,
        pqs_prototype_coefficient_transform_materialized = false,
        pqs_prototype_numerical_transform_materialized = false,
        pqs_prototype_source_operator_blocks_materialized = false,
        pqs_prototype_operator_blocks_materialized = false,
        pqs_prototype_pair_operator_blocks_materialized = false,
        pqs_prototype_hamiltonian_data_materialized = false,
        pqs_prototype_artifacts_materialized = false,
    )
end

function _pqs_source_box_route_driver_report_stage_lw_complete_shell_summary(
    assembly,
)
    transforms =
        hasproperty(assembly, :transforms) ? assembly.transforms : nothing
    available =
        !isnothing(transforms) &&
        hasproperty(transforms, :lw_complete_shell_cpb_enumeration_available) &&
        transforms.lw_complete_shell_cpb_enumeration_available

    return (;
        lw_complete_shell_cpb_enumeration_available = available,
        lw_complete_shell_region_count =
            available ? transforms.lw_complete_shell_region_count : 0,
        lw_complete_shell_cpb_count =
            available ? transforms.lw_complete_shell_cpb_count : 0,
        lw_complete_shell_cpb_family_counts =
            available ?
            transforms.lw_complete_shell_cpb_family_counts :
            (facet_cpb = 0, edge_cpb = 0, corner_cpb = 0),
        lw_complete_shell_enumeration_policy =
            available ? transforms.lw_complete_shell_enumeration_policy : nothing,
        lw_complete_shell_coefficient_maps_materialized =
            available &&
            transforms.lw_complete_shell_coefficient_maps_materialized,
        lw_complete_shell_operator_blocks_materialized =
            available &&
            transforms.lw_complete_shell_operator_blocks_materialized,
        lw_complete_shell_pair_operator_blocks_materialized =
            available &&
            transforms.lw_complete_shell_pair_operator_blocks_materialized,
        lw_complete_shell_hamiltonian_data_materialized =
            available &&
            transforms.lw_complete_shell_hamiltonian_data_materialized,
    )
end

function _pqs_source_box_route_driver_complete_core_shell_h1_j_summary(payload)
    isnothing(payload) && return (;
        status = :not_available_missing_complete_core_shell_h1_j_payload,
        blocker = :missing_complete_core_shell_h1_j_payload,
        final_dimension = nothing,
        h1_energy = nothing,
        self_coulomb = nothing,
        density_gauge = nothing,
        missing_inputs = (:complete_core_shell_h1_j_payload,),
        signed_final_weight_division_used = false,
        raw_no_division_used = false,
        density_normalized_pair_terms_used_as_authority = false,
        driver_route_materialized = false,
        rhf_materialized = false,
        gto_materialized = false,
        exports_materialized = false,
        artifacts_materialized = false,
    )
    return payload.summary
end

function _pqs_source_box_route_driver_complete_core_shell_h1_j_report_fields(
    assembly,
)
    payload =
        hasproperty(assembly, :complete_core_shell_h1_j_diagnostic_payload) ?
        assembly.complete_core_shell_h1_j_diagnostic_payload :
        nothing
    summary =
        _pqs_source_box_route_driver_complete_core_shell_h1_j_summary(payload)
    final_basis =
        isnothing(payload) || !hasproperty(payload, :final_basis) ?
        nothing :
        payload.final_basis
    source_plan =
        isnothing(payload) || !hasproperty(payload, :source_plan) ?
        nothing :
        payload.source_plan
    shell_records =
        isnothing(source_plan) || !hasproperty(source_plan, :shell_records) ?
        () :
        source_plan.shell_records
    return (;
        complete_core_shell_h1_j_diagnostic_summary = summary,
        complete_core_shell_h1_j_diagnostic_status = summary.status,
        complete_core_shell_h1_j_diagnostic_blocker = summary.blocker,
        complete_core_shell_h1_j_final_dimension = summary.final_dimension,
        complete_core_shell_h1_j_h1_energy = summary.h1_energy,
        complete_core_shell_h1_j_self_coulomb = summary.self_coulomb,
        complete_core_shell_h1_j_density_gauge = summary.density_gauge,
        complete_core_shell_h1_j_driver_route_materialized =
            summary.driver_route_materialized,
        complete_core_shell_core_support_count =
            isnothing(final_basis) ? nothing : final_basis.core_support_count,
        complete_core_shell_shell_support_count =
            isnothing(final_basis) ? nothing : final_basis.shell_support_count,
        complete_core_shell_shell_layer_count = length(shell_records),
        complete_core_shell_retained_per_shell =
            Tuple(record.retained_count for record in shell_records),
        complete_core_shell_shell_final_retained_count =
            isnothing(final_basis) ? nothing : final_basis.shell_final_retained_count,
        complete_core_shell_final_overlap_identity_error =
            isnothing(final_basis) ? nothing : final_basis.final_overlap_identity_error,
        complete_core_shell_raw_pair_factor_convention =
            get(summary, :raw_pair_factor_convention, nothing),
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
    wl_total = get(recipe, :wl_rhf_total, nothing)
    if isnothing(payload)
        return (;
            private_rhf_requested = false,
            private_rhf_summary = nothing,
            private_rhf_status = :not_requested,
            private_rhf_blocker = nothing,
            private_rhf_total_energy = nothing,
            private_rhf_iteration_count = nothing,
            private_rhf_converged = false,
            private_rhf_residual = nothing,
            private_rhf_mixing_kind = nothing,
            private_rhf_wl_total = wl_total,
            private_rhf_delta = nothing,
        )
    end
    summary = payload.summary
    total_energy = get(summary, :final_total_energy, nothing)
    residual_diagnostics = get(summary, :residual_diagnostics, (;))
    return (;
        private_rhf_requested = true,
        private_rhf_summary = summary,
        private_rhf_status = get(summary, :status, payload.status),
        private_rhf_blocker = get(summary, :blocker, payload.blocker),
        private_rhf_total_energy = total_energy,
        private_rhf_iteration_count = get(summary, :iteration_count, nothing),
        private_rhf_converged = get(summary, :rhf_converged, false),
        private_rhf_residual =
            get(residual_diagnostics, :commutator_residual, nothing),
        private_rhf_mixing_kind = get(summary, :mixing_kind, nothing),
        private_rhf_wl_total = wl_total,
        private_rhf_delta =
            isnothing(total_energy) || isnothing(wl_total) ?
            nothing :
            total_energy - wl_total,
    )
end

function _pqs_source_box_route_driver_diatomic_complete_core_shell_report_fields(
    assembly,
)
    readiness =
        hasproperty(assembly, :diatomic_complete_core_shell_ham_readiness_payload) ?
        assembly.diatomic_complete_core_shell_ham_readiness_payload :
        nothing
    summary =
        isnothing(readiness) || !hasproperty(readiness, :summary) ?
        (;
            status = :not_available_missing_diatomic_complete_core_shell_readiness,
            blocker = :missing_diatomic_complete_core_shell_readiness,
            source_plan_status = :not_available,
            final_basis_status = :not_available,
            h1_status = :not_available,
            h1_payload_status = :not_available,
            ham_input_payload_status = :not_available,
            hamiltonian_handoff_payload_status = :not_available,
            hamiltonian_consumer_contract_payload_status = :not_available,
            final_basis_materialized = false,
            h1_materialized = false,
            h1_j_materialized = false,
            ham_input_materialized = false,
            hamiltonian_handoff_materialized = false,
            hamiltonian_consumer_contract_materialized = false,
            rhf_materialized = false,
            public_api = false,
            exports_materialized = false,
            artifacts_materialized = false,
            missing_objects = (:diatomic_complete_core_shell_readiness,),
        ) :
        readiness.summary
    final_basis_payload =
        hasproperty(assembly, :diatomic_complete_core_shell_final_basis_payload) ?
        assembly.diatomic_complete_core_shell_final_basis_payload :
        nothing
    h1_payload =
        hasproperty(assembly, :diatomic_complete_core_shell_h1_payload) ?
        assembly.diatomic_complete_core_shell_h1_payload :
        nothing
    h1_hamiltonian =
        isnothing(h1_payload) || !hasproperty(h1_payload, :final_hamiltonian) ?
        nothing :
        h1_payload.final_hamiltonian
    ham_input_payload =
        hasproperty(assembly, :diatomic_complete_core_shell_ham_input_payload) ?
        assembly.diatomic_complete_core_shell_ham_input_payload :
        nothing
    handoff_payload =
        hasproperty(assembly, :diatomic_complete_core_shell_hamiltonian_handoff_payload) ?
        assembly.diatomic_complete_core_shell_hamiltonian_handoff_payload :
        nothing
    diagnostic_summary = (;
        final_dimension =
            isnothing(final_basis_payload) ||
            !hasproperty(final_basis_payload, :summary) ?
            nothing :
            get(final_basis_payload.summary, :final_dimension, nothing),
        final_overlap_identity_error =
            isnothing(final_basis_payload) ||
            !hasproperty(final_basis_payload, :final_basis) ||
            isnothing(final_basis_payload.final_basis) ||
            !hasproperty(final_basis_payload.final_basis, :final_overlap_identity_error) ?
            nothing :
            final_basis_payload.final_basis.final_overlap_identity_error,
        h1_lowest_energy =
            isnothing(h1_payload) || !hasproperty(h1_payload, :summary) ?
            nothing :
            get(h1_payload.summary, :lowest_energy, nothing),
        h1_hamiltonian_matrix_finite =
            isnothing(h1_hamiltonian) ?
            nothing :
            get(h1_hamiltonian, :hamiltonian_matrix_finite, nothing),
        h1_hamiltonian_symmetry_error =
            isnothing(h1_hamiltonian) ?
            nothing :
            get(h1_hamiltonian, :hamiltonian_matrix_symmetry_error, nothing),
        density_gauge =
            isnothing(ham_input_payload) || !hasproperty(ham_input_payload, :summary) ?
            nothing :
            get(ham_input_payload.summary, :density_gauge, nothing),
        raw_pair_factor_convention =
            isnothing(ham_input_payload) || !hasproperty(ham_input_payload, :summary) ?
            nothing :
            get(ham_input_payload.summary, :raw_pair_factor_convention, nothing),
        nuclear_repulsion =
            isnothing(handoff_payload) || !hasproperty(handoff_payload, :summary) ?
            nothing :
            get(handoff_payload.summary, :nuclear_repulsion, nothing),
        electron_count =
            isnothing(handoff_payload) || !hasproperty(handoff_payload, :summary) ?
            nothing :
            get(handoff_payload.summary, :electron_count, nothing),
        spin_sector =
            isnothing(handoff_payload) || !hasproperty(handoff_payload, :summary) ?
            nothing :
            get(handoff_payload.summary, :spin_sector, nothing),
    )
    return (;
        diatomic_complete_core_shell_readiness_summary =
            merge(summary, diagnostic_summary),
    )
end

function _pqs_source_box_route_driver_physical_gausslet_target_report_fields(
    assembly,
)
    payload =
        hasproperty(assembly, :diatomic_physical_gausslet_target_payload) ?
        assembly.diatomic_physical_gausslet_target_payload :
        nothing
    supplement_preflight_payload =
        hasproperty(assembly, :diatomic_physical_gausslet_supplement_preflight_payload) ?
        assembly.diatomic_physical_gausslet_supplement_preflight_payload :
        nothing
    supplement_request_payload =
        hasproperty(assembly, :diatomic_physical_gausslet_supplement_request_payload) ?
        assembly.diatomic_physical_gausslet_supplement_request_payload :
        nothing
    supplement_representation_payload =
        hasproperty(assembly, :diatomic_physical_gausslet_supplement_representation_payload) ?
        assembly.diatomic_physical_gausslet_supplement_representation_payload :
        nothing
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
    rhf_input_contract =
        hasproperty(assembly, :diatomic_physical_gausslet_rhf_input_contract) ?
        assembly.diatomic_physical_gausslet_rhf_input_contract :
        nothing
    rhf_execution_payload =
        hasproperty(assembly, :diatomic_physical_gausslet_rhf_execution_payload) ?
        assembly.diatomic_physical_gausslet_rhf_execution_payload :
        nothing
    wl_reference_candidate =
        hasproperty(assembly, :h2_wl_gausslet_only_reference_candidate) ?
        assembly.h2_wl_gausslet_only_reference_candidate :
        nothing
    summary =
        isnothing(payload) ?
        (;
            status = :not_available_missing_physical_gausslet_target_payload,
            blocker = :missing_physical_gausslet_target_payload,
        ) :
        payload.summary
    if !isnothing(supplement_request_payload)
        request_summary = supplement_request_payload.summary
        summary = merge(
            summary,
            (;
                supplement_request_status = request_summary.status,
                supplement_request_blocker = request_summary.blocker,
                supplement_request_fixture_label = request_summary.fixture_label,
                supplement_request_basis_name = request_summary.basis_name,
                supplement_request_lmax = request_summary.lmax,
                supplement_request_uncontracted = request_summary.uncontracted,
                supplement_request_atom_symbols = request_summary.atom_symbols,
                supplement_request_nuclear_charges = request_summary.nuclear_charges,
                supplement_request_bond_axis = request_summary.bond_axis,
                supplement_request_bond_length = request_summary.bond_length,
                supplement_request_required_provider_blocks =
                    request_summary.required_provider_blocks,
                supplement_request_missing_fact_labels =
                    request_summary.missing_fact_labels,
                supplement_request_matrices_materialized =
                    request_summary.matrices_materialized,
            ),
        )
    end
    if !isnothing(supplement_representation_payload)
        representation_summary = supplement_representation_payload.summary
        summary = merge(
            summary,
            (;
                supplement_representation_status = representation_summary.status,
                supplement_representation_blocker = representation_summary.blocker,
                supplement_representation_object_kind =
                    representation_summary.object_kind,
                supplement_representation_basis_name =
                    representation_summary.basis_name,
                supplement_representation_lmax = representation_summary.lmax,
                supplement_representation_atom_symbols =
                    representation_summary.atom_symbols,
                supplement_representation_center_count =
                    representation_summary.center_count,
                supplement_representation_orbital_count =
                    representation_summary.orbital_count,
                supplement_representation_matrices_materialized =
                    representation_summary.matrices_materialized,
                supplement_representation_provider_blocks_materialized =
                    representation_summary.provider_blocks_materialized,
            ),
        )
    end
    if !isnothing(supplement_preflight_payload)
        preflight_summary = supplement_preflight_payload.summary
        summary = merge(
            summary,
            (;
                supplement_preflight_status = preflight_summary.status,
                supplement_preflight_blocker = preflight_summary.blocker,
                supplement_preflight_fixture_label = preflight_summary.fixture_label,
                supplement_preflight_retained_transform_kind =
                    preflight_summary.retained_transform_kind,
                supplement_preflight_gausslet_final_dimension =
                    preflight_summary.gausslet_final_dimension,
                supplement_preflight_required_fact_labels =
                    preflight_summary.required_fact_labels,
                supplement_preflight_available_fact_labels =
                    preflight_summary.available_fact_labels,
                supplement_preflight_missing_fact_labels =
                    preflight_summary.missing_fact_labels,
                supplement_preflight_matrices_materialized =
                    preflight_summary.matrices_materialized,
                supplement_preflight_supplemented_values_materialized =
                    preflight_summary.supplemented_values_materialized,
            ),
        )
    end
    if !isnothing(source_plan_payload)
        independent_source_plan =
            get(summary, :source_plan_role, nothing) ===
            :independent_pqs_source_box_construction
        summary = merge(
            summary,
            (;
                source_plan_status = source_plan_payload.status,
                source_plan_blocker = source_plan_payload.blocker,
                source_plan_authority_status =
                    source_plan_payload.summary.source_plan_authority_status,
                source_plan_descriptor_status =
                    get(
                        source_plan_payload.summary,
                        :source_plan_descriptor_status,
                        :not_available,
                    ),
                source_plan_descriptor_blocker =
                    get(
                        source_plan_payload.summary,
                        :source_plan_descriptor_blocker,
                        nothing,
                    ),
                source_plan_family =
                    get(source_plan_payload.summary, :source_plan_family, :not_available),
                shared_shell_realization_status =
                    get(
                        source_plan_payload.summary,
                        :shared_shell_realization_status,
                        :not_available,
                    ),
                shared_shell_realization_blocker =
                    get(
                        source_plan_payload.summary,
                        :shared_shell_realization_blocker,
                        nothing,
                    ),
                shared_shell_realization_counts =
                    get(source_plan_payload.summary, :shared_shell_realization_counts, ()),
                shared_shell_realization_identity_errors =
                    get(
                        source_plan_payload.summary,
                        :shared_shell_realization_identity_errors,
                        (),
                    ),
                source_coefficients_materialized =
                    get(
                        source_plan_payload.summary,
                        :source_coefficients_materialized,
                        false,
                    ),
                independent_source_plan_blocker =
                    independent_source_plan ?
                    source_plan_payload.blocker :
                    get(summary, :independent_source_plan_blocker, nothing),
            ),
        )
    end
    if !isnothing(final_basis_payload)
        summary = merge(
            summary,
            (;
                final_basis_status = final_basis_payload.final_basis_status,
                final_basis_blocker = final_basis_payload.blocker,
                final_basis_materialized =
                    final_basis_payload.summary.final_basis_materialized,
                final_dimension = final_basis_payload.summary.final_dimension,
                final_overlap_identity_error =
                    final_basis_payload.summary.final_overlap_identity_error,
                physics_endpoint_blocker =
                    final_basis_payload.summary.endpoint_blocker,
            ),
        )
    end
    if !isnothing(h1_payload)
        summary = merge(
            summary,
            (;
                h1_status = h1_payload.h1_status,
                h1_blocker = h1_payload.blocker,
                h1_materialized = h1_payload.summary.h1_materialized,
                h1_lowest_energy = h1_payload.summary.lowest_energy,
                h1_hamiltonian_matrix_finite =
                    h1_payload.summary.h1_hamiltonian_matrix_finite,
                h1_hamiltonian_symmetry_error =
                    h1_payload.summary.h1_hamiltonian_symmetry_error,
                support_kinetic_status = h1_payload.support_kinetic_status,
                support_electron_nuclear_status =
                    h1_payload.support_electron_nuclear_status,
                final_kinetic_status = h1_payload.final_kinetic_status,
                final_electron_nuclear_status =
                    h1_payload.final_electron_nuclear_status,
                physics_endpoint_blocker = h1_payload.summary.endpoint_blocker,
            ),
        )
    end
    if !isnothing(h1_j_payload)
        summary = merge(
            summary,
            (;
                h1_j_status = h1_j_payload.h1_j_status,
                h1_j_blocker = h1_j_payload.blocker,
                h1_j_materialized = h1_j_payload.summary.h1_j_materialized,
                density_interaction_materialized =
                    h1_j_payload.summary.density_interaction_materialized,
                density_interaction_status =
                    h1_j_payload.summary.density_interaction_status,
                density_gauge = h1_j_payload.summary.density_gauge,
                raw_pair_factor_convention =
                    h1_j_payload.summary.raw_pair_factor_convention,
                support_weight_count = h1_j_payload.summary.support_weight_count,
                support_weights_all_positive =
                    h1_j_payload.summary.support_weights_all_positive,
                support_raw_pair_shape =
                    h1_j_payload.summary.support_raw_pair_shape,
                support_raw_pair_finite =
                    h1_j_payload.summary.support_raw_pair_finite,
                pre_final_pair_matrix_shape =
                    h1_j_payload.summary.pre_final_pair_matrix_shape,
                pre_final_pair_matrix_finite =
                    h1_j_payload.summary.pre_final_pair_matrix_finite,
                pre_final_pair_matrix_symmetry_error =
                    h1_j_payload.summary.pre_final_pair_matrix_symmetry_error,
                h1_j_self_coulomb = h1_j_payload.summary.self_coulomb,
                physics_endpoint_blocker = h1_j_payload.summary.endpoint_blocker,
            ),
        )
    end
    if !isnothing(rhf_input_contract)
        contract_summary = rhf_input_contract.summary
        summary = merge(
            summary,
            (;
                private_rhf_input_contract_status = contract_summary.status,
                private_rhf_input_contract_blocker = contract_summary.blocker,
                private_rhf_input_contract_available =
                    contract_summary.input_contract_available,
                private_rhf_electron_count = contract_summary.electron_count,
                private_rhf_occupation_policy =
                    contract_summary.occupation_policy,
                private_rhf_occupation_nocc = contract_summary.occupation_nocc,
                private_rhf_h1_matrix_available =
                    contract_summary.h1_matrix_available,
                private_rhf_h1_matrix_finite =
                    contract_summary.h1_matrix_finite,
                private_rhf_h1_matrix_symmetry_error =
                    contract_summary.h1_matrix_symmetry_error,
                private_rhf_density_interaction_available =
                    contract_summary.density_interaction_available,
                private_rhf_final_to_pre_final_transform_available =
                    contract_summary.final_to_pre_final_transform_available,
                private_rhf_pre_final_pair_matrix_available =
                    contract_summary.pre_final_pair_matrix_available,
                private_rhf_pre_final_pair_matrix_finite =
                    contract_summary.pre_final_pair_matrix_finite,
                private_rhf_pre_final_pair_matrix_symmetry_error =
                    contract_summary.pre_final_pair_matrix_symmetry_error,
                private_rhf_materialized =
                    contract_summary.private_rhf_materialized,
                physics_endpoint_blocker = contract_summary.endpoint_blocker,
            ),
        )
    end
    if !isnothing(rhf_execution_payload)
        execution_summary = rhf_execution_payload.summary
        summary = merge(
            summary,
            (;
                private_rhf_execution_status = execution_summary.status,
                private_rhf_execution_blocker = execution_summary.blocker,
                private_rhf_executed = execution_summary.executed,
                private_rhf_materialized = execution_summary.materialized,
                private_rhf_converged = execution_summary.converged,
                private_rhf_total_energy = execution_summary.total_energy,
                private_rhf_one_body_energy = execution_summary.one_body_energy,
                private_rhf_two_body_energy = execution_summary.two_body_energy,
                private_rhf_iteration_count = execution_summary.iteration_count,
                private_rhf_density_trace = execution_summary.density_trace,
                private_rhf_idempotency_residual =
                    execution_summary.idempotency_residual,
                private_rhf_commutator_residual =
                    execution_summary.commutator_residual,
                private_rhf_energy_delta = execution_summary.energy_delta,
                private_rhf_final_density_one_step_consistency_status =
                    execution_summary.final_density_one_step_consistency_status,
                physics_endpoint_blocker = execution_summary.endpoint_blocker,
            ),
        )
    end
    if !isnothing(wl_reference_candidate)
        candidate_summary = wl_reference_candidate.summary
        endpoint_blocker =
            candidate_summary.status ===
            :available_wl_h2_gausslet_only_reference_candidate &&
            get(summary, :private_rhf_materialized, false) ?
            :missing_wl_h2_gausslet_only_reference_values :
            get(summary, :physics_endpoint_blocker, candidate_summary.blocker)
        summary = merge(
            summary,
            (;
                wl_reference_candidate_status = candidate_summary.status,
                wl_reference_candidate_blocker = candidate_summary.blocker,
                wl_reference_final_dimension = candidate_summary.final_dimension,
                wl_reference_retained_transform_kind =
                    candidate_summary.retained_transform_kind,
                wl_reference_supplement_policy =
                    candidate_summary.supplement_policy,
                wl_reference_label = candidate_summary.label,
                wl_reference_mismatches = candidate_summary.mismatches,
                old_supplemented_wl_qw_scalar_references_blocked =
                    candidate_summary.old_supplemented_scalar_references_blocked,
                physics_endpoint_blocker = endpoint_blocker,
            ),
        )
    end
    return (;
        physical_gausslet_target_summary = summary,
        physical_gausslet_target_status = summary.status,
        physical_gausslet_target_blocker = summary.blocker,
    )
end

function _pqs_source_box_route_driver_pair_operator_report_count_entries(
    entries,
)
    isnothing(entries) && return ()
    return Tuple(
        _pqs_source_box_route_driver_pair_operator_report_count_entry(entry)
        for entry in entries
    )
end

function _pqs_source_box_route_driver_pair_operator_report_count_entry(entry)
    names = propertynames(entry)
    (:pair_count in names || !(:count in names)) && return entry
    value_names = Tuple(name for name in names if name != :count)
    values = Tuple(getproperty(entry, name) for name in value_names)
    return merge(
        NamedTuple{value_names}(values),
        (; pair_count = getproperty(entry, :count)),
    )
end

function _pqs_source_box_route_driver_pair_operator_report_aliases(
    terminal_route_state,
    fallback_source,
)
    summary =
        !isnothing(terminal_route_state) &&
        hasproperty(terminal_route_state, :pair_operator_summary) ?
        terminal_route_state.pair_operator_summary :
        nothing
    if !isnothing(summary) &&
       hasproperty(summary, :route_core_pair_operator_plan_inventory_available)
        return (;
            route_core_typed_pair_operator_plan_inventory_available =
                summary.route_core_pair_operator_plan_inventory_available,
            route_core_typed_pair_operator_plan_inventory_status =
                fallback_source.route_core_typed_pair_operator_plan_inventory_status,
            route_core_typed_pair_operator_plan_blocker =
                fallback_source.route_core_typed_pair_operator_plan_blocker,
            route_core_typed_pair_operator_plan_count =
                summary.route_core_pair_operator_plan_count,
            route_core_typed_pair_operator_plan_blocked_count =
                summary.route_core_pair_operator_plan_blocked_count,
            route_core_typed_pair_operator_plan_materialized =
                summary.materialized,
            route_core_typed_pair_operator_source_path_counts =
                _pqs_source_box_route_driver_pair_operator_report_count_entries(
                    summary.source_operator_path_counts,
                ),
            route_core_typed_pair_operator_final_block_path_counts =
                _pqs_source_box_route_driver_pair_operator_report_count_entries(
                    summary.final_block_path_counts,
                ),
            route_core_typed_pair_operator_materialization_status_counts =
                _pqs_source_box_route_driver_pair_operator_report_count_entries(
                    summary.materialization_status_counts,
                ),
            route_core_typed_pair_operator_blocker_counts =
                _pqs_source_box_route_driver_pair_operator_report_count_entries(
                    summary.blocker_counts,
                ),
        )
    end

    return (;
        route_core_typed_pair_operator_plan_inventory_available =
            fallback_source.route_core_typed_pair_operator_plan_inventory_available,
        route_core_typed_pair_operator_plan_inventory_status =
            fallback_source.route_core_typed_pair_operator_plan_inventory_status,
        route_core_typed_pair_operator_plan_blocker =
            fallback_source.route_core_typed_pair_operator_plan_blocker,
        route_core_typed_pair_operator_plan_count =
            fallback_source.route_core_typed_pair_operator_plan_count,
        route_core_typed_pair_operator_plan_blocked_count =
            fallback_source.route_core_typed_pair_operator_plan_blocked_count,
        route_core_typed_pair_operator_plan_materialized =
            fallback_source.route_core_typed_pair_operator_plan_materialized,
        route_core_typed_pair_operator_source_path_counts =
            fallback_source.route_core_typed_pair_operator_source_path_counts,
        route_core_typed_pair_operator_final_block_path_counts =
            fallback_source.route_core_typed_pair_operator_final_block_path_counts,
        route_core_typed_pair_operator_materialization_status_counts =
            fallback_source.route_core_typed_pair_operator_materialization_status_counts,
        route_core_typed_pair_operator_blocker_counts =
            fallback_source.route_core_typed_pair_operator_blocker_counts,
    )
end

function _pqs_source_box_route_driver_report_stage_low_order_route_summary(
    assembly,
)
    low_order_assembly =
        hasproperty(assembly, :low_order_assembly) ?
        assembly.low_order_assembly :
        nothing
    pqs_prototype_summary =
        _pqs_source_box_route_driver_report_stage_pqs_prototype_summary(
            assembly,
        )
    lw_complete_shell_summary =
        _pqs_source_box_route_driver_report_stage_lw_complete_shell_summary(
            assembly,
        )
    if isnothing(low_order_assembly)
        return (;
            object_kind = :cartesian_report_stage_low_order_route_summary,
            status = :not_available_missing_assembly_stage_summary,
            low_order_shellization_policy_requested = nothing,
            low_order_shellization_policy_resolved = :not_available,
            low_order_shellization_policy_source = :not_available,
            low_order_shellization_policy_status =
                :not_available_missing_assembly_stage_summary,
            low_order_shellization_policy_blocker =
                :missing_assembly_stage_low_order_summary,
            shellization_source = :not_available,
            shellization_kind = :not_available,
            unit_route_kind = :not_available,
            transform_route_kind = :not_available,
            pair_route_kind = :not_available,
            assembly_source = :not_available,
            assembly_route_kind = :not_available,
            assembly_kind = :not_available,
            atom_growth_selected = false,
            terminal_shellification_selected = false,
            legacy_source_selected = false,
            terminal_shellification_summary_available = false,
            terminal_shellification_scaffold_available = false,
            terminal_shellification_scaffold = nothing,
            terminal_shellification_region_count = 0,
            terminal_shellification_unit_inventory_available = false,
            terminal_shellification_unit_inventory = nothing,
            terminal_shellification_unit_count = 0,
            terminal_shellification_unit_keys = (),
            terminal_shellification_unit_roles = (),
            terminal_shellification_unit_kinds = (),
            terminal_shellification_unit_support_counts = (),
            terminal_shellification_lowering_contract_inventory_available = false,
            terminal_shellification_lowering_contract_inventory_status =
                :not_available,
            terminal_shellification_lowering_contract_inventory = nothing,
            terminal_shellification_lowering_plan_available = false,
            terminal_shellification_lowering_plan_status = :not_available,
            terminal_shellification_lowering_plan = nothing,
            terminal_shellification_lowering_summary = nothing,
            terminal_shellification_lowering_contract_count = 0,
            terminal_shellification_lowering_contract_kinds = (),
            terminal_shellification_lowering_contract_kind_counts =
                (
                    direct_core_identity_cpb_count = 0,
                    direct_slab_identity_cpb_count = 0,
                    direct_boundary_slab_identity_cpb_count = 0,
                    white_lindsey_boundary_strata_count = 0,
                    pqs_filled_source_cpb_count = 0,
                    distorted_product_box_comx_count = 0,
                ),
            _pqs_source_box_route_driver_selected_terminal_lowering_fields(
                nothing,
                :not_available,
                nothing,
            )...,
            terminal_shellification_contract_counts_by_unit = (),
            terminal_shellification_lw_complete_shell_cpb_count = 0,
            terminal_shellification_lw_complete_shell_cpb_family_counts =
                (facet_cpb = 0, edge_cpb = 0, corner_cpb = 0),
            terminal_shellification_final_retained_unit_inventory_available = false,
            terminal_shellification_transform_contracts_available = false,
            terminal_shellification_pair_inventory_available = false,
            terminal_shellification_pair_inventory_status = :not_available,
            terminal_shellification_pair_materialization_status = :not_available,
            terminal_shellification_assembly_materialization_status =
                :not_available,
            terminal_shellification_central_gap_region_count = 0,
            terminal_shellification_central_midpoint_slab_count = 0,
            terminal_shellification_central_distorted_product_box_count = 0,
            terminal_shellification_central_distorted_product_box_metadata = (),
            materialization_required = true,
            materialization_status =
                :blocked_missing_assembly_stage_summary,
            materialization_blocker =
                :missing_assembly_stage_low_order_summary,
            hamiltonian_matrices_materialized = false,
            operator_matrices_materialized = false,
            pair_operator_blocks_materialized = false,
            pair_inventory_source = :not_available,
            pair_inventory_known = false,
            independent_atom_growth_pair_inventory_available = false,
            pair_count = 0,
            pair_family_counts = nothing,
            route_core_final_unit_count = 0,
            route_core_pair_inventory_available = false,
            route_core_pair_inventory_status = :not_available,
            route_core_pair_count = 0,
            route_core_pair_order_matches_staged = false,
            route_core_pair_order_comparison_source = :not_available,
            route_core_pair_family_counts = (),
            route_core_summary_status = :not_available,
            route_core_pair_operator_ready = false,
            route_core_pair_operator_readiness_status = :not_available,
            route_core_pair_operator_blocker = :not_available,
            route_core_pair_operator_readiness_requirements =
                _pqs_source_box_route_driver_route_core_readiness_requirements(),
            route_core_pair_operator_preflight_available = false,
            route_core_pair_operator_preflight_status = :not_available,
            route_core_pair_operator_preflight = nothing,
            route_core_pair_operator_preflight_blocker = :not_available,
            route_core_pair_operator_plan_available = false,
            route_core_pair_operator_plan_status = :not_available,
            route_core_pair_operator_plan = nothing,
            route_core_pair_operator_plan_blocker = :not_available,
            _pqs_source_box_route_driver_route_core_typed_pair_operator_plan_unavailable_metadata(
                :not_available,
                :not_available,
            )...,
            plan_authority = false,
            active_source_authority = false,
            legacy_source_authority = false,
            pqs_prototype_summary...,
            lw_complete_shell_summary...,
            report_stage_fields_preserved = false,
            summary_only = true,
        )
    end

    pair_operator_report_aliases =
        _pqs_source_box_route_driver_pair_operator_report_aliases(
            low_order_assembly.terminal_route_state,
            low_order_assembly,
        )

    return (;
        object_kind = :cartesian_report_stage_low_order_route_summary,
        status =
            low_order_assembly.status ==
            :available_assembly_stage_low_order_summary ?
            :available_report_stage_low_order_route_summary :
            low_order_assembly.status,
        low_order_shellization_policy_requested =
            low_order_assembly.low_order_shellization_policy_requested,
        low_order_shellization_policy_resolved =
            low_order_assembly.low_order_shellization_policy_resolved,
        low_order_shellization_policy_source =
            low_order_assembly.low_order_shellization_policy_source,
        low_order_shellization_policy_status =
            low_order_assembly.low_order_shellization_policy_status,
        low_order_shellization_policy_blocker =
            low_order_assembly.low_order_shellization_policy_blocker,
        shellization_source = low_order_assembly.shellization_source,
        shellization_kind = low_order_assembly.shellization_kind,
        unit_route_kind = low_order_assembly.unit_route_kind,
        transform_route_kind = low_order_assembly.transform_route_kind,
        pair_route_kind = low_order_assembly.pair_route_kind,
        assembly_source = low_order_assembly.assembly_source,
        assembly_route_kind = low_order_assembly.assembly_route_kind,
        assembly_kind = low_order_assembly.assembly_kind,
        atom_growth_selected = low_order_assembly.atom_growth_assembly_selected,
        terminal_shellification_selected =
            low_order_assembly.terminal_shellification_assembly_selected,
        legacy_source_selected =
            low_order_assembly.legacy_source_assembly_selected,
        terminal_shellification_summary_available =
            low_order_assembly.terminal_shellification_assembly_summary_available,
        terminal_shellification_lowering_plan_available =
            low_order_assembly.terminal_route_summary.lowering_plan_available,
        terminal_shellification_lowering_plan_status =
            low_order_assembly.terminal_route_summary.lowering_plan_status,
        terminal_shellification_lowering_plan =
            low_order_assembly.terminal_route_state.lowering_plan,
        terminal_shellification_lowering_summary =
            low_order_assembly.terminal_route_state.lowering_summary,
        _pqs_source_box_route_driver_terminal_shellification_alias_fields(
            low_order_assembly,
            low_order_assembly.terminal_shellification_assembly_selected;
            include_crc_sidecar_summary = false,
        )...,
        terminal_shellification_transform_contracts_available =
            low_order_assembly.terminal_shellification_transform_contracts_available,
        terminal_shellification_pair_inventory_available =
            low_order_assembly.terminal_shellification_pair_inventory_available,
        terminal_shellification_pair_inventory_status =
            low_order_assembly.terminal_shellification_pair_inventory_status,
        terminal_shellification_pair_materialization_status =
            low_order_assembly.terminal_shellification_pair_materialization_status,
        terminal_shellification_assembly_materialization_status =
            low_order_assembly.terminal_shellification_assembly_materialization_status,
        materialization_required =
            low_order_assembly.assembly_requires_materialization,
        materialization_status =
            low_order_assembly.assembly_materialization_status,
        materialization_blocker = low_order_assembly.assembly_blocker,
        hamiltonian_matrices_materialized =
            low_order_assembly.hamiltonian_matrices_materialized,
        operator_matrices_materialized =
            low_order_assembly.operator_matrices_materialized,
        pair_operator_blocks_materialized =
            low_order_assembly.pair_operator_blocks_materialized,
        pair_inventory_source = low_order_assembly.pair_inventory_source,
        pair_inventory_known = low_order_assembly.pair_inventory_known,
        independent_atom_growth_pair_inventory_available =
            low_order_assembly.independent_atom_growth_pair_inventory_available,
        pair_count = low_order_assembly.pair_count,
        pair_family_counts = low_order_assembly.pair_family_counts,
        route_core_final_unit_count =
            low_order_assembly.route_core_final_unit_count,
        route_core_pair_inventory_available =
            low_order_assembly.route_core_pair_inventory_available,
        route_core_pair_inventory_status =
            low_order_assembly.route_core_pair_inventory_status,
        route_core_pair_count = low_order_assembly.route_core_pair_count,
        route_core_pair_order_matches_staged =
            low_order_assembly.route_core_pair_order_matches_staged,
        route_core_pair_order_comparison_source =
            low_order_assembly.route_core_pair_order_comparison_source,
        route_core_pair_family_counts =
            low_order_assembly.route_core_pair_family_counts,
        route_core_summary_status = low_order_assembly.route_core_summary_status,
        route_core_pair_operator_ready =
            low_order_assembly.route_core_pair_operator_ready,
        route_core_pair_operator_readiness_status =
            low_order_assembly.route_core_pair_operator_readiness_status,
        route_core_pair_operator_blocker =
            low_order_assembly.route_core_pair_operator_blocker,
        route_core_pair_operator_readiness_requirements =
            low_order_assembly.route_core_pair_operator_readiness_requirements,
        route_core_pair_operator_preflight_available =
            low_order_assembly.route_core_pair_operator_preflight_available,
        route_core_pair_operator_preflight_status =
            low_order_assembly.route_core_pair_operator_preflight_status,
        route_core_pair_operator_preflight =
            low_order_assembly.route_core_pair_operator_preflight,
        route_core_pair_operator_preflight_blocker =
            low_order_assembly.route_core_pair_operator_preflight_blocker,
        route_core_pair_operator_plan_available =
            low_order_assembly.route_core_pair_operator_plan_available,
        route_core_pair_operator_plan_status =
            low_order_assembly.route_core_pair_operator_plan_status,
        route_core_pair_operator_plan =
            low_order_assembly.route_core_pair_operator_plan,
        route_core_pair_operator_plan_blocker =
            low_order_assembly.route_core_pair_operator_plan_blocker,
        pair_operator_report_aliases...,
        route_core_typed_pair_operator_plan_family_counts =
            low_order_assembly.route_core_typed_pair_operator_plan_family_counts,
        route_core_typed_pair_operator_materialization_ready =
            low_order_assembly.route_core_typed_pair_operator_materialization_ready,
        route_core_typed_pair_operator_materialization_readiness_status =
            low_order_assembly.route_core_typed_pair_operator_materialization_readiness_status,
        route_core_typed_pair_operator_materialization_readiness_blocker =
            low_order_assembly.route_core_typed_pair_operator_materialization_readiness_blocker,
        route_core_typed_pair_operator_materialization_readiness_requirements =
            low_order_assembly.route_core_typed_pair_operator_materialization_readiness_requirements,
        route_core_typed_pair_operator_materialization_readiness_plan_count =
            low_order_assembly.route_core_typed_pair_operator_materialization_readiness_plan_count,
        route_core_typed_pair_operator_materialization_readiness_blocked_count =
            low_order_assembly.route_core_typed_pair_operator_materialization_readiness_blocked_count,
        route_core_typed_pair_operator_materialization_readiness_materialized_count =
            low_order_assembly.route_core_typed_pair_operator_materialization_readiness_materialized_count,
        plan_authority = low_order_assembly.plan_authority,
        active_source_authority = low_order_assembly.active_source_authority,
        legacy_source_authority = low_order_assembly.legacy_source_authority,
        pqs_prototype_summary...,
        lw_complete_shell_summary...,
        report_stage_fields_preserved = true,
        summary_only = true,
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

function _pqs_source_box_route_driver_pqs_prototype_print_line(report)
    hasproperty(report, :low_order_pqs_transform_prototype_available) ||
        return nothing
    report.low_order_pqs_transform_prototype_available || return nothing
    source =
        report.low_order_pqs_prototype_source_cpb_kind == :filled_source_cpb ?
        "filled CPB" :
        string(report.low_order_pqs_prototype_source_cpb_kind)
    owned_support =
        report.low_order_pqs_prototype_owned_support_is_cpb ? "CPB" : "shell"
    retained_space =
        report.low_order_pqs_prototype_intermediate_retained_space ==
        :boundary_comx_product_mode_selection ?
        "boundary COMX products" :
        string(report.low_order_pqs_prototype_intermediate_retained_space)
    realization =
        report.low_order_pqs_prototype_shell_realization ==
        :shell_projection_lowdin_deferred ?
        "shell projection + Lowdin deferred" :
        string(report.low_order_pqs_prototype_shell_realization)
    materialized =
        report.low_order_pqs_prototype_source_operator_blocks_materialized ||
        report.low_order_pqs_prototype_operator_blocks_materialized ||
        report.low_order_pqs_prototype_pair_operator_blocks_materialized ||
        report.low_order_pqs_prototype_hamiltonian_data_materialized ||
        report.low_order_pqs_prototype_artifacts_materialized
    return string(
        "CPB/PQS prototype: metadata-only, unit=",
        report.low_order_pqs_prototype_unit_key,
        ", source=",
        source,
        ", owned support is ",
        owned_support,
        ", retained space=",
        retained_space,
        ", realization=",
        realization,
        ", operators/materialization=",
        materialized ? "yes" : "no",
    )
end

function _pqs_source_box_route_driver_print_pqs_prototype_summary(report)
    line = _pqs_source_box_route_driver_pqs_prototype_print_line(report)
    isnothing(line) || println(line)
    return nothing
end

function _pqs_source_box_route_driver_terminal_shellification_print_line(report)
    hasproperty(report, :low_order_terminal_shellification_selected) ||
        return nothing
    report.low_order_terminal_shellification_selected || return nothing
    materialized =
        report.low_order_operator_matrices_materialized ||
        report.low_order_pair_operator_blocks_materialized ||
        report.low_order_hamiltonian_matrices_materialized
    return string(
        "Terminal shellification: selected, regions ",
        report.low_order_terminal_shellification_region_count,
        ", central gaps ",
        report.low_order_terminal_shellification_central_gap_region_count,
        ", midpoint slabs ",
        report.low_order_terminal_shellification_central_midpoint_slab_count,
        ", central distorted product boxes ",
        report.low_order_terminal_shellification_central_distorted_product_box_count,
        ", pair inventory ",
        repr(report.low_order_terminal_shellification_pair_inventory_status),
        ", assembly/materialization ",
        repr(
            report.low_order_terminal_shellification_assembly_materialization_status,
        ),
        ", operators/materialization=",
        materialized ? "yes" : "no",
    )
end

function _pqs_source_box_route_driver_print_terminal_shellification_summary(report)
    line =
        _pqs_source_box_route_driver_terminal_shellification_print_line(report)
    isnothing(line) || println(line)
    return nothing
end

function _pqs_source_box_route_driver_crc_pair_family_label(pair_family)
    length(pair_family) == 2 || return string(pair_family)
    return string(pair_family[1], " / ", pair_family[2])
end

function _pqs_source_box_route_driver_crc_pair_family_print_line(report)
    hasproperty(report, :low_order_route_core_pair_family_counts) ||
        return nothing
    families = report.low_order_route_core_pair_family_counts
    isempty(families) && return nothing
    return string(
        "CRC pair families: ",
        join(
            (
                string(
                    _pqs_source_box_route_driver_crc_pair_family_label(
                        family.pair_family,
                    ),
                    "=",
                    family.pair_count,
                ) for family in families
            ),
            ", ",
        ),
    )
end

function _pqs_source_box_route_driver_crc_sidecar_print_line(report)
    hasproperty(report, :low_order_route_core_summary_status) || return nothing
    report.low_order_route_core_pair_inventory_available || return string(
        "CRC sidecars: ",
        report.low_order_route_core_summary_status,
    )
    order_match =
        report.low_order_route_core_pair_order_matches_staged ? "yes" : "no"
    return string(
        "CRC sidecars: final units ",
        report.low_order_route_core_final_unit_count,
        ", pairs ",
        report.low_order_route_core_pair_count,
        ", order match ",
        order_match,
    )
end

function _pqs_source_box_route_driver_print_crc_sidecar_summary(report)
    line = _pqs_source_box_route_driver_crc_sidecar_print_line(report)
    isnothing(line) || println(line)
    family_line = _pqs_source_box_route_driver_crc_pair_family_print_line(report)
    isnothing(family_line) || println(family_line)
    return nothing
end

function _pqs_source_box_route_driver_crc_operator_plan_print_line(report)
    hasproperty(report, :low_order_route_core_pair_operator_plan_status) ||
        return nothing
    status = report.low_order_route_core_pair_operator_plan_status
    if status == :ready_route_core_pair_operator_plan
        return string(
            "CRC pair-operator plan: ready metadata plan, final units ",
            report.low_order_route_core_pair_operator_plan.route_core_final_unit_count,
            ", pairs ",
            report.low_order_route_core_pair_operator_plan.route_core_pair_count,
            ", operator blocks materialized no",
        )
    end
    if status == :blocked_route_core_pair_operator_plan
        return string(
            "CRC pair-operator plan: blocked metadata plan, blocker ",
            report.low_order_route_core_pair_operator_plan_blocker,
            ", operator blocks materialized no",
        )
    end
    return string("CRC pair-operator plan: ", status)
end

function _pqs_source_box_route_driver_print_crc_operator_plan_summary(report)
    line = _pqs_source_box_route_driver_crc_operator_plan_print_line(report)
    isnothing(line) || println(line)
    return nothing
end

function _pqs_source_box_route_driver_crc_typed_pair_operator_plan_print_line(
    report,
)
    hasproperty(
        report,
        :low_order_route_core_typed_pair_operator_plan_inventory_available,
    ) || return nothing
    availability =
        report.low_order_route_core_typed_pair_operator_plan_inventory_available ?
        "available" :
        "unavailable"
    materialized =
        report.low_order_route_core_typed_pair_operator_plan_materialized ?
        "yes" :
        "no"
    return string(
        "CRC typed pair-operator inventory: ",
        availability,
        " (",
        repr(
            report.low_order_route_core_typed_pair_operator_plan_inventory_status,
        ),
        "), typed plans ",
        report.low_order_route_core_typed_pair_operator_plan_count,
        ", blocked ",
        report.low_order_route_core_typed_pair_operator_plan_blocked_count,
        ", materialized ",
        materialized,
        ", blocker ",
        repr(report.low_order_route_core_typed_pair_operator_plan_blocker),
    )
end

function _pqs_source_box_route_driver_print_crc_typed_pair_operator_plan_summary(
    report,
)
    line =
        _pqs_source_box_route_driver_crc_typed_pair_operator_plan_print_line(
            report,
        )
    isnothing(line) || println(line)
    return nothing
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
    @show report.low_order_shellization_policy_resolved
    @show report.low_order_shellization_policy_source
    @show report.atom_growth_low_order_route_selected
    @show report.low_order_active_source_authority
    @show report.low_order_materialization_required
    @show report.low_order_materialization_status
    @show report.low_order_pair_inventory_source
    @show report.low_order_pair_count
    _pqs_source_box_route_driver_print_terminal_shellification_summary(report)
    _pqs_source_box_route_driver_print_crc_sidecar_summary(report)
    _pqs_source_box_route_driver_print_crc_operator_plan_summary(report)
    _pqs_source_box_route_driver_print_crc_typed_pair_operator_plan_summary(
        report,
    )
    @show report.low_order_hamiltonian_matrices_materialized
    _pqs_source_box_route_driver_print_pqs_prototype_summary(report)
    @show materialization.basis_artifact_status materialization.basis_artifact_written
    @show materialization.status materialization.ham_bundle_export_status
    @show materialization.materialized_report_kind
    @show materialization.route_configured_system_classification
    @show materialization.route_configured_shellization_request_status
    @show materialization.route_configured_shellization_planning_family
    @show materialization.route_configured_midpoint_slab_status
    @show materialization.route_configured_primary_planned_helper
    @show materialization.route_configured_missing_input_count
    @show materialization.route_configured_input_readiness_status
    @show materialization.route_configured_available_fact_count
    @show materialization.route_configured_materializer_config_status
    @show materialization.route_configured_materializer_config_pending_input_count
    @show materialization.route_configured_materializer_backend_requested
    @show materialization.route_configured_materializer_backend_consumed
    @show materialization.route_configured_materializer_d_requested
    @show materialization.route_configured_materializer_d_consumed
    @show materialization.route_configured_materializer_nside_requested
    @show materialization.route_configured_materializer_nside_consumed
    @show materialization.low_order_shellization_policy_requested
    @show materialization.low_order_shellization_policy_resolved
    @show materialization.low_order_shellization_policy_source
    @show materialization.low_order_shellization_policy_status
    @show materialization.route_configured_one_center_materializer_probe_requested
    @show materialization.route_configured_one_center_materializer_probe_status
    @show materialization.route_configured_one_center_materializer_probe_blocker
    @show materialization.route_configured_diatomic_materializer_probe_requested
    @show materialization.route_configured_diatomic_materializer_probe_status
    @show materialization.route_configured_diatomic_materializer_probe_blocker
    @show materialization.route_configured_diatomic_atom_growth_materializer_probe_requested
    @show materialization.route_configured_diatomic_atom_growth_materializer_probe_status
    @show materialization.route_configured_diatomic_atom_growth_materializer_probe_consumed
    @show materialization.route_configured_diatomic_atom_growth_shellification_consumed
    @show materialization.route_configured_diatomic_atom_growth_basis_adapter_status
    @show materialization.route_configured_diatomic_atom_growth_final_integral_weights_status
    @show materialization.route_configured_diatomic_atom_growth_ham_adapter_status
    @show materialization.route_configured_legacy_diatomic_source_consumed
    @show materialization.shellization_source materialization.route_configured_shellization_consumed
    @show materialization.ham_artifact_status materialization.ham_artifact_written
    return nothing
end

function cartesian_print_details(report, materialization)
    _pqs_source_box_route_driver_print_details(report)
    if materialization.materialize_route_requested ||
       materialization.route_configured_one_center_materializer_probe_requested ||
       materialization.route_configured_diatomic_materializer_probe_requested ||
       materialization.route_configured_diatomic_atom_growth_materializer_probe_requested ||
       materialization.save_basis_artifact_requested ||
       materialization.save_ham_artifact_requested
        _pqs_source_box_route_driver_print_materialization(materialization)
    end
    return nothing
end

function cartesian_save(report, save_inputs, materialization)
    return _pqs_source_box_route_driver_save(
        report;
        save_inputs...,
        materialization,
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
    route_inputs = (;
        route_family, route_kind, route_shape, product_body_rule,
        pqs_retained_rule, product_retained_rule, terms, pair_factor_normalization,
        support_dense_direct_allowed, reference_only_authorities,
        white_lindsey_route_shape, white_lindsey_mapping_rule,
        white_lindsey_nesting_rule, white_lindsey_retained_rule,
        white_lindsey_operator_rule, white_lindsey_benchmark_role,
    )

    system = cartesian_system(system_inputs)
    recipe = cartesian_recipe(route_inputs)
    parent = cartesian_parent(system, spacing_inputs, probe_inputs, recipe)
    shells = cartesian_shells(parent, spacing_inputs, recipe)
    units = cartesian_units(parent, shells, probe_inputs, recipe)
    transforms = cartesian_transforms(units, recipe)
    pairs = cartesian_pair_terms(units, transforms, recipe)
    assembly = cartesian_assembly(parent, shells, units, transforms, pairs, recipe)
    return cartesian_report(system, parent, assembly, recipe)
end

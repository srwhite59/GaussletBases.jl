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

    family =
        hasproperty(parent_inputs, :parent_axis_probe_family) ?
        parent_inputs.parent_axis_probe_family :
        :G10
    axis_cache = Dict{Int, Any}()
    function _axis_for_count(count)
        axis_count = Int(count)
        get!(axis_cache, axis_count) do
            build_basis(MappedUniformBasisSpec(
                family;
                count = axis_count,
                mapping = IdentityMapping(),
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
                mapping = :IdentityMapping,
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

function _pqs_source_box_route_driver_one_center_materializer_probe(
    config;
    probe_route_configured_one_center_materializer::Bool = false,
    white_lindsey_expansion = nothing,
)
    if !probe_route_configured_one_center_materializer
        return (;
            object_kind = :route_configured_one_center_materializer_probe,
            requested = false,
            status = :not_requested,
            materialized = false,
            route_configured_shellization_consumed = false,
            blocker = nothing,
            materialization = nothing,
            error_message = nothing,
        )
    elseif config.system_classification != :one_center
        return (;
            object_kind = :route_configured_one_center_materializer_probe,
            requested = true,
            status = :blocked_not_one_center,
            materialized = false,
            route_configured_shellization_consumed = false,
            blocker = :route_config_not_one_center,
            materialization = nothing,
            error_message = nothing,
        )
    elseif config.route_family != :white_lindsey_low_order
        return (;
            object_kind = :route_configured_one_center_materializer_probe,
            requested = true,
            status = :blocked_not_white_lindsey_low_order,
            materialized = false,
            route_configured_shellization_consumed = false,
            blocker = :route_config_not_white_lindsey_low_order,
            materialization = nothing,
            error_message = nothing,
        )
    elseif !config.materializer_options_ready
        return (;
            object_kind = :route_configured_one_center_materializer_probe,
            requested = true,
            status = :blocked_missing_materializer_options,
            materialized = false,
            route_configured_shellization_consumed = false,
            blocker = config.materializer_option_blocker,
            materialization = nothing,
            error_message =
                "missing materializer options: $(config.missing_materializer_options)",
        )
    end

    expansion =
        isnothing(white_lindsey_expansion) ?
        coulomb_gaussian_expansion(doacc = false) :
        white_lindsey_expansion
    try
        materialization =
            _cartesian_shellization_route_materialize_one_center_low_order(
                config;
                expansion,
                gausslet_backend = config.materializer_backend_requested,
                d = config.materializer_d_requested,
                nside = config.materializer_nside_requested,
                reference_spacing = config.materializer_reference_spacing_requested,
                tail_spacing = config.materializer_tail_spacing_requested,
            )
        return (;
            object_kind = :route_configured_one_center_materializer_probe,
            requested = true,
            status = materialization.status,
            materialized = true,
            route_configured_shellization_consumed =
                materialization.route_configured_shellization_consumed,
            blocker = nothing,
            materialization,
            error_message = nothing,
        )
    catch error
        error isa ArgumentError || rethrow()
        return (;
            object_kind = :route_configured_one_center_materializer_probe,
            requested = true,
            status = :blocked_materializer_precondition,
            materialized = false,
            route_configured_shellization_consumed = false,
            blocker = :materializer_precondition_failed,
            materialization = nothing,
            error_message = sprint(showerror, error),
        )
    end
end

function _pqs_source_box_route_driver_diatomic_materializer_probe(
    config;
    probe_route_configured_diatomic_materializer::Bool = false,
    route_materializer_payload = nothing,
    white_lindsey_expansion = nothing,
    shared_shell_layer_policy = nothing,
    packet_kernel = nothing,
)
    parent_qw_basis_object =
        isnothing(route_materializer_payload) ?
        nothing :
        route_materializer_payload.parent_qw_basis_object
    parent_axis_bundle_object =
        isnothing(route_materializer_payload) ?
        nothing :
        route_materializer_payload.parent_axis_bundle_object
    axis_bundle_backend =
        isnothing(route_materializer_payload) ?
        nothing :
        route_materializer_payload.axis_bundle_backend
    parent_qw_basis_object_handoff_available = !isnothing(parent_qw_basis_object)
    parent_axis_bundle_object_handoff_available =
        !isnothing(parent_axis_bundle_object)
    axis_bundle_backend_handoff_available = !isnothing(axis_bundle_backend)

    if !probe_route_configured_diatomic_materializer
        return (;
            object_kind = :route_configured_diatomic_materializer_probe,
            requested = false,
            status = :not_requested,
            materialized = false,
            route_configured_shellization_consumed = false,
            blocker = nothing,
            missing_contract = (),
            parent_qw_basis_object_handoff_available,
            parent_axis_bundle_object_handoff_available,
            axis_bundle_backend_handoff_available,
            materialization = nothing,
            error_message = nothing,
        )
    elseif config.system_classification != :bond_aligned_diatomic
        return (;
            object_kind = :route_configured_diatomic_materializer_probe,
            requested = true,
            status = :blocked_not_bond_aligned_diatomic,
            materialized = false,
            route_configured_shellization_consumed = false,
            blocker = :route_config_not_bond_aligned_diatomic,
            missing_contract = (),
            parent_qw_basis_object_handoff_available,
            parent_axis_bundle_object_handoff_available,
            axis_bundle_backend_handoff_available,
            materialization = nothing,
            error_message = nothing,
        )
    elseif config.route_family != :white_lindsey_low_order
        return (;
            object_kind = :route_configured_diatomic_materializer_probe,
            requested = true,
            status = :blocked_not_white_lindsey_low_order,
            materialized = false,
            route_configured_shellization_consumed = false,
            blocker = :route_config_not_white_lindsey_low_order,
            missing_contract = (),
            parent_qw_basis_object_handoff_available,
            parent_axis_bundle_object_handoff_available,
            axis_bundle_backend_handoff_available,
            materialization = nothing,
            error_message = nothing,
        )
    end

    try
        materialization =
            _cartesian_shellization_route_materialize_bond_aligned_diatomic(
                config;
                parent_qw_basis_object,
                parent_axis_bundle_object,
                expansion = white_lindsey_expansion,
                shared_shell_layer_policy,
                packet_kernel,
                axis_bundle_backend,
            )
        return (;
            object_kind = :route_configured_diatomic_materializer_probe,
            requested = true,
            status = materialization.status,
            materialized = materialization.materialized,
            route_configured_shellization_consumed =
                materialization.route_configured_shellization_consumed,
            blocker = materialization.blocker,
            missing_contract = materialization.missing_contract,
            parent_qw_basis_object_handoff_available,
            parent_axis_bundle_object_handoff_available,
            axis_bundle_backend_handoff_available,
            materialization,
            error_message = nothing,
        )
    catch error
        error isa ArgumentError || rethrow()
        return (;
            object_kind = :route_configured_diatomic_materializer_probe,
            requested = true,
            status = :blocked_materializer_precondition,
            materialized = false,
            route_configured_shellization_consumed = false,
            blocker = :materializer_precondition_failed,
            missing_contract = (),
            parent_qw_basis_object_handoff_available,
            parent_axis_bundle_object_handoff_available,
            axis_bundle_backend_handoff_available,
            materialization = nothing,
            error_message = sprint(showerror, error),
        )
    end
end

function _pqs_source_box_route_driver_route_configured_one_center_report(
    materialization,
)
    fixture = materialization.fixture
    inventory = fixture.inventory
    route_units = _white_lindsey_low_order_materialized_seed_route_units(fixture)
    operator_inventory =
        _white_lindsey_low_order_materialized_seed_operator_inventory(fixture)
    operator_pairs_materialized =
        route_units.operator_pairs_materialized ||
        operator_inventory.operator_pairs_materialized
    shellization_summary = materialization.shellization_summary

    return (;
        object_kind = :white_lindsey_low_order_route_configured_one_center_report,
        route_family = :white_lindsey_low_order,
        status = :private_development_route_configured,
        private_development_only = true,
        materialization,
        fixture,
        inventory,
        route_units,
        operator_inventory,
        shellization_summary,
        shellization_summary_available = true,
        shellization_source = :route_configured_one_center_low_order,
        route_configured_shellization_consumed = true,
        materialized_shellization_stage = shellization_summary.shellization_stage,
        seed_materialization_status = :not_seed_route_configured_materialization,
        packet_kernel = fixture.packet_kernel,
        retained_dimension = route_units.retained_dimension,
        operator_pairs_materialized,
        electron_electron_materialized = operator_inventory.electron_electron_materialized,
        weight_semantics = :retained_basis_integral_weights,
    )
end

function _pqs_source_box_route_driver_materialization(
    report;
    materialize_route::Bool = false,
    probe_route_configured_one_center_materializer::Bool = false,
    save_basis_artifact::Bool = false,
    save_ham_artifact::Bool = false,
    basisfile::AbstractString = "cartesian_nesting_route_driver_basis_bundle.jld2",
    hamfile::AbstractString = "cartesian_nesting_route_driver_ham_bundle.jld2",
    materializer_backend = nothing,
    materializer_nside = nothing,
    white_lindsey_expansion = nothing,
    white_lindsey_Z = nothing,
)
    route_family = report.route_family
    route_configured_shellization_request =
        _cartesian_shellization_route_configured_request(
            report;
            materializer_backend,
            materializer_nside,
        )
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
    route_configured_shellization_helper_map =
        _cartesian_shellization_route_planning_helper_map(route_configured_shellization_plan)
    route_configured_shellization_helper_map_status =
        route_configured_shellization_helper_map.status
    route_configured_primary_planned_helper =
        route_configured_shellization_helper_map.primary_planned_helper
    route_configured_missing_input_count =
        route_configured_shellization_helper_map.missing_input_count
    route_configured_helper_map_blocker =
        route_configured_shellization_helper_map.blocker
    route_configured_input_readiness =
        _cartesian_shellization_route_materializer_input_readiness(
            route_configured_shellization_request,
            route_configured_shellization_plan,
            route_configured_shellization_helper_map,
        )
    route_configured_input_readiness_status = route_configured_input_readiness.status
    route_configured_available_fact_count =
        route_configured_input_readiness.available_fact_count
    route_configured_materializer_missing_input_count =
        route_configured_input_readiness.missing_input_count
    route_configured_input_readiness_blocker =
        route_configured_input_readiness.blocker
    route_configured_materializer_config =
        _cartesian_shellization_route_materializer_config(
            route_configured_shellization_request,
            route_configured_shellization_plan,
            route_configured_shellization_helper_map,
            route_configured_input_readiness,
        )
    route_configured_materializer_config_status =
        route_configured_materializer_config.status
    route_configured_materializer_config_planning_family =
        route_configured_materializer_config.planning_family
    route_configured_materializer_config_pending_input_count =
        route_configured_materializer_config.pending_input_count
    route_configured_one_center_materializer_requested =
        probe_route_configured_one_center_materializer ||
        (
            materialize_route &&
            route_family == :white_lindsey_low_order &&
            route_configured_system_classification == :one_center
        )
    route_configured_one_center_materializer_probe =
        _pqs_source_box_route_driver_one_center_materializer_probe(
            route_configured_materializer_config;
            probe_route_configured_one_center_materializer =
                route_configured_one_center_materializer_requested,
            white_lindsey_expansion,
        )
    route_configured_one_center_materializer_probe_requested =
        route_configured_one_center_materializer_probe.requested
    route_configured_one_center_materializer_probe_status =
        route_configured_one_center_materializer_probe.status
    route_configured_one_center_materializer_probe_materialized =
        route_configured_one_center_materializer_probe.materialized
    route_configured_one_center_materializer_probe_consumed =
        route_configured_one_center_materializer_probe.route_configured_shellization_consumed
    route_configured_one_center_materializer_probe_blocker =
        route_configured_one_center_materializer_probe.blocker
    route_configured_diatomic_materializer_requested =
        materialize_route &&
        route_family == :white_lindsey_low_order &&
        route_configured_system_classification == :bond_aligned_diatomic
    route_materializer_payload =
        hasproperty(report, :route_materializer_payload) ?
        report.route_materializer_payload :
        nothing
    route_configured_diatomic_shared_shell_layer_policy =
        route_configured_diatomic_materializer_requested &&
        route_configured_materializer_config.materializer_backend_requested ==
        :pgdg_localized_experimental ?
        :endcap_panel_owned :
        nothing
    route_configured_diatomic_packet_kernel =
        route_configured_diatomic_materializer_requested ?
        :factorized_direct :
        nothing
    route_configured_diatomic_policy_source =
        isnothing(route_configured_diatomic_shared_shell_layer_policy) ?
        nothing :
        :existing_endcap_panel_owned_pgdg_route
    route_configured_diatomic_materializer_probe =
        _pqs_source_box_route_driver_diatomic_materializer_probe(
            route_configured_materializer_config;
            probe_route_configured_diatomic_materializer =
                route_configured_diatomic_materializer_requested,
            route_materializer_payload,
            white_lindsey_expansion,
            shared_shell_layer_policy =
                route_configured_diatomic_shared_shell_layer_policy,
            packet_kernel = route_configured_diatomic_packet_kernel,
        )
    route_configured_diatomic_materializer_probe_requested =
        route_configured_diatomic_materializer_probe.requested
    route_configured_diatomic_materializer_probe_status =
        route_configured_diatomic_materializer_probe.status
    route_configured_diatomic_materializer_probe_materialized =
        route_configured_diatomic_materializer_probe.materialized
    route_configured_diatomic_materializer_probe_consumed =
        route_configured_diatomic_materializer_probe.route_configured_shellization_consumed
    route_configured_diatomic_materializer_probe_blocker =
        route_configured_diatomic_materializer_probe.blocker
    route_configured_diatomic_materializer_missing_contract =
        route_configured_diatomic_materializer_probe.missing_contract
    route_configured_diatomic_materializer_payload_available =
        !isnothing(route_materializer_payload)
    route_configured_diatomic_parent_qw_basis_object_handoff_available =
        route_configured_diatomic_materializer_probe.parent_qw_basis_object_handoff_available
    route_configured_diatomic_parent_axis_bundle_object_handoff_available =
        route_configured_diatomic_materializer_probe.parent_axis_bundle_object_handoff_available
    route_configured_diatomic_axis_bundle_backend_handoff_available =
        route_configured_diatomic_materializer_probe.axis_bundle_backend_handoff_available
    route_configured_diatomic_axis_bundle_backend_handoff =
        isnothing(route_materializer_payload) ?
        nothing :
        route_materializer_payload.axis_bundle_backend
    route_configured_diatomic_seed_fallback =
        route_configured_diatomic_materializer_probe_requested &&
        !route_configured_diatomic_materializer_probe_consumed
    route_configured_materializer_backend_requested =
        route_configured_materializer_config.materializer_backend_requested
    route_configured_materializer_backend_source =
        route_configured_materializer_config.materializer_backend_source
    route_configured_materializer_backend_status =
        route_configured_materializer_config.materializer_backend_status
    route_configured_materializer_d_requested =
        route_configured_materializer_config.materializer_d_requested
    route_configured_materializer_d_source =
        route_configured_materializer_config.materializer_d_source
    route_configured_materializer_d_status =
        route_configured_materializer_config.materializer_d_status
    route_configured_materializer_nside_requested =
        route_configured_materializer_config.materializer_nside_requested
    route_configured_materializer_nside_source =
        route_configured_materializer_config.materializer_nside_source
    route_configured_materializer_nside_status =
        route_configured_materializer_config.materializer_nside_status
    route_configured_materializer_reference_spacing_requested =
        route_configured_materializer_config.materializer_reference_spacing_requested
    route_configured_materializer_reference_spacing_source =
        route_configured_materializer_config.materializer_reference_spacing_source
    route_configured_materializer_reference_spacing_status =
        route_configured_materializer_config.materializer_reference_spacing_status
    route_configured_materializer_tail_spacing_requested =
        route_configured_materializer_config.materializer_tail_spacing_requested
    route_configured_materializer_tail_spacing_source =
        route_configured_materializer_config.materializer_tail_spacing_source
    route_configured_materializer_tail_spacing_status =
        route_configured_materializer_config.materializer_tail_spacing_status
    route_configured_materializer_options_ready =
        route_configured_materializer_config.materializer_options_ready
    route_configured_materializer_missing_options =
        route_configured_materializer_config.missing_materializer_options
    route_configured_materializer_option_blocker =
        route_configured_materializer_config.materializer_option_blocker
    route_configured_materializer_consumed_options =
        route_configured_one_center_materializer_probe_materialized ?
        route_configured_one_center_materializer_probe.materialization.materializer_options :
        route_configured_diatomic_materializer_probe_materialized ?
        route_configured_diatomic_materializer_probe.materialization.materializer_options :
        nothing
    route_configured_materializer_backend_consumed =
        isnothing(route_configured_materializer_consumed_options) ?
        nothing :
        route_configured_materializer_consumed_options.gausslet_backend
    route_configured_materializer_d_consumed =
        isnothing(route_configured_materializer_consumed_options) ?
        nothing :
        route_configured_materializer_consumed_options.d
    route_configured_materializer_nside_consumed =
        isnothing(route_configured_materializer_consumed_options) ?
        nothing :
        route_configured_materializer_consumed_options.nside
    route_configured_materializer_reference_spacing_consumed =
        isnothing(route_configured_materializer_consumed_options) ?
        nothing :
        route_configured_materializer_consumed_options.reference_spacing
    route_configured_materializer_tail_spacing_consumed =
        isnothing(route_configured_materializer_consumed_options) ?
        nothing :
        route_configured_materializer_consumed_options.tail_spacing
    route_configured_materializer_contract = (;
        route_configured_materializer_backend_requested,
        route_configured_materializer_backend_source,
        route_configured_materializer_backend_status,
        route_configured_materializer_backend_consumed,
        route_configured_materializer_d_requested,
        route_configured_materializer_d_source,
        route_configured_materializer_d_status,
        route_configured_materializer_d_consumed,
        route_configured_materializer_nside_requested,
        route_configured_materializer_nside_source,
        route_configured_materializer_nside_status,
        route_configured_materializer_nside_consumed,
        route_configured_materializer_reference_spacing_requested,
        route_configured_materializer_reference_spacing_source,
        route_configured_materializer_reference_spacing_status,
        route_configured_materializer_reference_spacing_consumed,
        route_configured_materializer_tail_spacing_requested,
        route_configured_materializer_tail_spacing_source,
        route_configured_materializer_tail_spacing_status,
        route_configured_materializer_tail_spacing_consumed,
        route_configured_materializer_options_ready,
        route_configured_materializer_missing_options,
        route_configured_materializer_option_blocker,
    )
    route_configured_diatomic_materializer_contract = (;
        route_configured_diatomic_materializer_probe,
        route_configured_diatomic_materializer_probe_requested,
        route_configured_diatomic_materializer_probe_status,
        route_configured_diatomic_materializer_probe_materialized,
        route_configured_diatomic_materializer_probe_consumed,
        route_configured_diatomic_materializer_probe_blocker,
        route_configured_diatomic_materializer_missing_contract,
        route_configured_diatomic_materializer_payload_available,
        route_configured_diatomic_parent_qw_basis_object_handoff_available,
        route_configured_diatomic_parent_axis_bundle_object_handoff_available,
        route_configured_diatomic_axis_bundle_backend_handoff_available,
        route_configured_diatomic_axis_bundle_backend_handoff,
        route_configured_diatomic_shared_shell_layer_policy,
        route_configured_diatomic_packet_kernel,
        route_configured_diatomic_policy_source,
        route_configured_diatomic_seed_fallback,
    )

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
            route_configured_shellization_helper_map,
            route_configured_shellization_helper_map_available = true,
            route_configured_shellization_helper_map_status,
            route_configured_primary_planned_helper,
            route_configured_missing_input_count,
            route_configured_helper_map_blocker,
            route_configured_input_readiness,
            route_configured_input_readiness_available = true,
            route_configured_input_readiness_status,
            route_configured_available_fact_count,
            route_configured_materializer_missing_input_count,
            route_configured_input_readiness_blocker,
            route_configured_materializer_config,
            route_configured_materializer_config_available = true,
            route_configured_materializer_config_status,
            route_configured_materializer_config_planning_family,
            route_configured_materializer_config_pending_input_count,
            route_configured_one_center_materializer_probe,
            route_configured_one_center_materializer_probe_requested,
            route_configured_one_center_materializer_probe_status,
            route_configured_one_center_materializer_probe_materialized,
            route_configured_one_center_materializer_probe_consumed,
            route_configured_one_center_materializer_probe_blocker,
            route_configured_diatomic_materializer_contract...,
            route_configured_materializer_contract...,
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
        route_configured_one_center_report_required =
            materialize_route &&
            route_configured_system_classification == :one_center
        if route_configured_one_center_report_required &&
           !route_configured_one_center_materializer_probe_materialized
            throw(
                ArgumentError(
                    "route-configured one-center materializer failed: " *
                    string(route_configured_one_center_materializer_probe_blocker) *
                    " " *
                    string(route_configured_one_center_materializer_probe.error_message),
                ),
            )
        end
        use_route_configured_one_center_report =
            route_configured_one_center_report_required &&
            route_configured_one_center_materializer_probe_materialized
        use_route_configured_diatomic_shellization =
            route_configured_diatomic_materializer_requested &&
            route_configured_diatomic_materializer_probe_consumed
        if use_route_configured_diatomic_shellization
            diatomic_materialization =
                route_configured_diatomic_materializer_probe.materialization
            shellization_summary = diatomic_materialization.shellization_summary
            shellization_summary_available = !isnothing(shellization_summary)
            basis_artifact_status =
                save_basis_artifact ?
                :not_written_route_configured_diatomic_basis_export_pending :
                :not_requested
            ham_artifact_status =
                save_ham_artifact ?
                :not_written_route_configured_diatomic_ham_export_pending :
                :not_requested
            ham_export_blocker =
                save_ham_artifact ?
                :pending_route_configured_diatomic_ham_export :
                nothing
            return (;
                object_kind = :cartesian_nesting_route_driver_materialization,
                route_family,
                private_development_only = true,
                materialize_route_requested = true,
                save_basis_artifact_requested = save_basis_artifact,
                save_ham_artifact_requested = save_ham_artifact,
                status = :materialized_route_configured_diatomic_shellization_available,
                materialized_report = nothing,
                materialized_report_kind = diatomic_materialization.object_kind,
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
                route_configured_shellization_helper_map,
                route_configured_shellization_helper_map_available = true,
                route_configured_shellization_helper_map_status,
                route_configured_primary_planned_helper,
                route_configured_missing_input_count,
                route_configured_helper_map_blocker,
                route_configured_input_readiness,
                route_configured_input_readiness_available = true,
                route_configured_input_readiness_status,
                route_configured_available_fact_count,
                route_configured_materializer_missing_input_count,
                route_configured_input_readiness_blocker,
                route_configured_materializer_config,
                route_configured_materializer_config_available = true,
                route_configured_materializer_config_status,
                route_configured_materializer_config_planning_family,
                route_configured_materializer_config_pending_input_count,
                route_configured_one_center_materializer_probe,
                route_configured_one_center_materializer_probe_requested,
                route_configured_one_center_materializer_probe_status,
                route_configured_one_center_materializer_probe_materialized,
                route_configured_one_center_materializer_probe_consumed,
                route_configured_one_center_materializer_probe_blocker,
                route_configured_diatomic_materializer_contract...,
                route_configured_materializer_contract...,
                shellization_summary,
                shellization_summary_available,
                shellization_source = :route_configured_bond_aligned_diatomic_source,
                route_configured_shellization_consumed = true,
                materialized_shellization_stage =
                    shellization_summary.shellization_stage,
                seed_materialization_status =
                    :not_seed_route_configured_diatomic_shellization,
                retained_dimension = diatomic_materialization.retained_dimension,
                final_integral_weights_status =
                    :pending_route_configured_diatomic_fixed_block_adapter,
                one_body_operator_status =
                    :pending_route_configured_diatomic_operator_inventory,
                basis_bundle_export_status =
                    :pending_route_configured_diatomic_basis_export,
                basis_artifact_status,
                basis_artifact_written = false,
                basisfile,
                basis_artifact_path = nothing,
                basis_export_blocker =
                    save_basis_artifact ?
                    :pending_route_configured_diatomic_basis_export :
                    nothing,
                ham_preflight_status =
                    save_ham_artifact ?
                    :blocked_route_configured_diatomic_ham_export_not_adopted :
                    :not_requested,
                ham_missing_builder =
                    save_ham_artifact ?
                    :pending_route_configured_diatomic_ham_builder :
                    nothing,
                ham_operator_payload_status =
                    :pending_route_configured_diatomic_operator_payload,
                ham_interaction_status =
                    :pending_route_configured_diatomic_density_density_builder,
                ham_bundle_export_status =
                    :pending_route_configured_diatomic_ham_export,
                ham_artifact_status,
                ham_artifact_written = false,
                hamfile,
                ham_export_blocker,
                ham_preflight = nothing,
                pqs_materialization_status = :not_applicable,
            )
        end
        materialized_report =
            use_route_configured_one_center_report ?
            _pqs_source_box_route_driver_route_configured_one_center_report(
                route_configured_one_center_materializer_probe.materialization,
            ) :
            _white_lindsey_low_order_materialized_seed_report()
        basis_export_status =
            use_route_configured_one_center_report ?
            :supported_route_configured_one_center_basis_only_fixed_block :
            :supported_basis_only_fixed_block
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
            save_basis_artifact ?
            (
                use_route_configured_one_center_report ?
                :written_route_configured_one_center_basis_only_bundle :
                :written_basis_only_bundle
            ) :
            :not_requested
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
                    route_configured_shellization_helper_map_status,
                    route_configured_primary_planned_helper,
                    route_configured_missing_input_count,
                    route_configured_helper_map_blocker,
                    route_configured_input_readiness_status,
                    route_configured_available_fact_count,
                    route_configured_materializer_missing_input_count,
                    route_configured_input_readiness_blocker,
                    route_configured_materializer_config_status,
                    route_configured_materializer_config_planning_family,
                    route_configured_materializer_config_pending_input_count,
                    route_configured_one_center_materializer_probe_requested,
                    route_configured_one_center_materializer_probe_status,
                    route_configured_one_center_materializer_probe_materialized,
                    route_configured_one_center_materializer_probe_consumed,
                    route_configured_one_center_materializer_probe_blocker,
                    route_configured_diatomic_materializer_contract...,
                    route_configured_materializer_contract...,
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
                    route_configured_shellization_helper_map_status,
                    route_configured_primary_planned_helper,
                    route_configured_missing_input_count,
                    route_configured_helper_map_blocker,
                    route_configured_input_readiness_status,
                    route_configured_available_fact_count,
                    route_configured_materializer_missing_input_count,
                    route_configured_input_readiness_blocker,
                    route_configured_materializer_config_status,
                    route_configured_materializer_config_planning_family,
                    route_configured_materializer_config_pending_input_count,
                    route_configured_one_center_materializer_probe_requested,
                    route_configured_one_center_materializer_probe_status,
                    route_configured_one_center_materializer_probe_materialized,
                    route_configured_one_center_materializer_probe_consumed,
                    route_configured_one_center_materializer_probe_blocker,
                    route_configured_diatomic_materializer_contract...,
                    route_configured_materializer_contract...,
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
            ham_artifact_status =
                use_route_configured_one_center_report ?
                :written_route_configured_one_center_ham_bundle :
                :written_white_lindsey_low_order_ham_bundle
        end
        return (;
            object_kind = :cartesian_nesting_route_driver_materialization,
            route_family,
            private_development_only = true,
            materialize_route_requested = true,
            save_basis_artifact_requested = save_basis_artifact,
            save_ham_artifact_requested = save_ham_artifact,
            status =
                use_route_configured_one_center_report ?
                :materialized_route_configured_one_center_report_available :
                :materialized_seed_report_available,
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
            route_configured_shellization_helper_map,
            route_configured_shellization_helper_map_available = true,
            route_configured_shellization_helper_map_status,
            route_configured_primary_planned_helper,
            route_configured_missing_input_count,
            route_configured_helper_map_blocker,
            route_configured_input_readiness,
            route_configured_input_readiness_available = true,
            route_configured_input_readiness_status,
            route_configured_available_fact_count,
            route_configured_materializer_missing_input_count,
            route_configured_input_readiness_blocker,
            route_configured_materializer_config,
            route_configured_materializer_config_available = true,
            route_configured_materializer_config_status,
            route_configured_materializer_config_planning_family,
            route_configured_materializer_config_pending_input_count,
            route_configured_one_center_materializer_probe,
            route_configured_one_center_materializer_probe_requested,
            route_configured_one_center_materializer_probe_status,
            route_configured_one_center_materializer_probe_materialized,
            route_configured_one_center_materializer_probe_consumed,
            route_configured_one_center_materializer_probe_blocker,
            route_configured_diatomic_materializer_contract...,
            route_configured_materializer_contract...,
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
        route_configured_shellization_helper_map,
        route_configured_shellization_helper_map_available = true,
        route_configured_shellization_helper_map_status,
        route_configured_primary_planned_helper,
        route_configured_missing_input_count,
        route_configured_helper_map_blocker,
        route_configured_input_readiness,
        route_configured_input_readiness_available = true,
        route_configured_input_readiness_status,
        route_configured_available_fact_count,
        route_configured_materializer_missing_input_count,
        route_configured_input_readiness_blocker,
        route_configured_materializer_config,
        route_configured_materializer_config_available = true,
        route_configured_materializer_config_status,
        route_configured_materializer_config_planning_family,
        route_configured_materializer_config_pending_input_count,
        route_configured_one_center_materializer_probe,
        route_configured_one_center_materializer_probe_requested,
        route_configured_one_center_materializer_probe_status,
        route_configured_one_center_materializer_probe_materialized,
        route_configured_one_center_materializer_probe_consumed,
        route_configured_one_center_materializer_probe_blocker,
        route_configured_diatomic_materializer_contract...,
        route_configured_materializer_contract...,
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

function cartesian_shells(parent, spacing_inputs, recipe)
    route_skeleton =
        _pqs_source_box_route_driver_route_skeleton(
            parent.route_axis_counts, spacing_inputs, recipe)

    return (;
        object_kind = :cartesian_shells,
        status = route_skeleton.status,
        spacing_inputs,
        route_skeleton,
        route_shape = route_skeleton.route_shape,
        source_boxes = route_skeleton.source_boxes,
        shellization_stage = :represented_by_route_skeleton,
    )
end

function cartesian_units(parent, shells, route_inputs, recipe)
    raw_box =
        _pqs_source_box_route_driver_raw_box_probe(
            parent.standard_setup, shells.route_skeleton, parent.parent_axis,
            parent.route_axis_counts, route_inputs, recipe)

    return (;
        object_kind = :cartesian_units,
        status = shells.route_skeleton.status,
        route_inputs,
        route_skeleton = shells.route_skeleton,
        raw_box,
        source_boxes = shells.route_skeleton.source_boxes,
        source_dimensions = shells.route_skeleton.source_dimensions,
        retained_units = shells.route_skeleton.retained_units,
        retained_unit_order = shells.route_skeleton.retained_unit_order,
        unit_stage = :broken_into_source_units,
    )
end

function cartesian_transforms(units, recipe)
    retained_units = units.retained_units
    return (;
        object_kind = :cartesian_transforms,
        status = units.status,
        route_family = recipe.route_family,
        retained_units,
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

function cartesian_pair_terms(units, transforms, recipe)
    route_skeleton = units.route_skeleton

    return (;
        object_kind = :cartesian_pair_terms,
        status = units.status,
        route_family = recipe.route_family,
        retained_dimension = transforms.retained_dimension,
        pair_entries = route_skeleton.pair_entries,
        pair_family_counts = route_skeleton.pair_family_counts,
        helper_by_pair_family = route_skeleton.helper_by_pair_family,
        pair_stage = :pair_operator_terms_described,
    )
end

function cartesian_assembly(parent, shells, units, transforms, pairs, recipe)
    route_skeleton = shells.route_skeleton
    route_facts = _pqs_source_box_route_driver_route_facts(route_skeleton)
    contract = _pqs_source_box_route_driver_contract_metadata(recipe)

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
        assembly_stage = :assembled_report_inputs,
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

    return _pqs_source_box_route_driver_report(
        standard_setup, parent, parent_axis, route_axis_counts, raw_box,
        system_metadata, recipe_metadata, parent_contract, parent_description,
        route_skeleton, route_facts, contract, diagnostics)
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
    @show materialization.basis_artifact_status materialization.basis_artifact_written
    @show materialization.status materialization.ham_bundle_export_status
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
    @show materialization.route_configured_one_center_materializer_probe_requested
    @show materialization.route_configured_one_center_materializer_probe_status
    @show materialization.route_configured_one_center_materializer_probe_blocker
    @show materialization.route_configured_diatomic_materializer_probe_requested
    @show materialization.route_configured_diatomic_materializer_probe_status
    @show materialization.route_configured_diatomic_materializer_probe_blocker
    @show materialization.shellization_source materialization.route_configured_shellization_consumed
    @show materialization.ham_artifact_status materialization.ham_artifact_written
    return nothing
end

function cartesian_print_details(report, materialization)
    _pqs_source_box_route_driver_print_details(report)
    if materialization.materialize_route_requested ||
       materialization.route_configured_one_center_materializer_probe_requested ||
       materialization.route_configured_diatomic_materializer_probe_requested ||
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

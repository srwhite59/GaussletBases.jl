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

function _pqs_source_box_route_driver_materializer_payload_property(
    payload,
    field::Symbol,
)
    return !isnothing(payload) && hasproperty(payload, field) ?
           getproperty(payload, field) :
           nothing
end

function _pqs_source_box_route_driver_diatomic_atom_growth_materializer_probe(
    config;
    probe_route_configured_diatomic_atom_growth_materializer::Bool = false,
    route_materializer_payload = nothing,
    white_lindsey_expansion = nothing,
    packet_kernel = nothing,
)
    parent_qw_basis_object =
        _pqs_source_box_route_driver_materializer_payload_property(
            route_materializer_payload,
            :parent_qw_basis_object,
        )
    parent_axis_bundle_object =
        _pqs_source_box_route_driver_materializer_payload_property(
            route_materializer_payload,
            :parent_axis_bundle_object,
        )
    axis_bundle_backend =
        _pqs_source_box_route_driver_materializer_payload_property(
            route_materializer_payload,
            :axis_bundle_backend,
        )
    parent_qw_basis_object_handoff_available = !isnothing(parent_qw_basis_object)
    parent_axis_bundle_object_handoff_available =
        !isnothing(parent_axis_bundle_object)
    axis_bundle_backend_handoff_available = !isnothing(axis_bundle_backend)

    blocked_probe(status, blocker; missing_contract = (), error_message = nothing) =
        (;
            object_kind =
                :route_configured_diatomic_atom_growth_materializer_probe,
            requested = probe_route_configured_diatomic_atom_growth_materializer,
            status,
            materialized = false,
            atom_growth_shellification_consumed = false,
            route_configured_diatomic_atom_growth_shellification_consumed = false,
            route_configured_shellization_consumed = false,
            route_configured_legacy_diatomic_source_consumed = false,
            blocker,
            missing_contract = Tuple(missing_contract),
            parent_qw_basis_object_handoff_available,
            parent_axis_bundle_object_handoff_available,
            axis_bundle_backend_handoff_available,
            axis_bundle_backend_handoff = axis_bundle_backend,
            materializer_options = nothing,
            construction_plan = nothing,
            scaffold = nothing,
            materialization = nothing,
            sequence_available = false,
            retained_dimension = nothing,
            support_count = nothing,
            coverage_audit = nothing,
            coverage_complete = false,
            atom_growth_construction_plan_authority = false,
            active_source_authority = false,
            active_source_oracle_comparison_run = false,
            route_behavior_changed = false,
            shellification_plan_path_used = false,
            calls_shellification_plan_materializer = false,
            fixed_block_available = false,
            fixed_block_status = :not_requested,
            basis_adapter = nothing,
            basis_adapter_summary = nothing,
            basis_adapter_status = :not_requested,
            basis_adapter_blocker = nothing,
            final_integral_weights_status = :not_requested,
            ham_adapter = nothing,
            ham_adapter_summary = nothing,
            ham_adapter_status = :not_requested,
            ham_adapter_blocker = nothing,
            error_message,
        )

    if !probe_route_configured_diatomic_atom_growth_materializer
        return blocked_probe(:not_requested, nothing)
    elseif config.system_classification != :bond_aligned_diatomic
        return blocked_probe(
            :blocked_not_bond_aligned_diatomic,
            :route_config_not_bond_aligned_diatomic,
        )
    elseif config.route_family != :white_lindsey_low_order
        return blocked_probe(
            :blocked_not_white_lindsey_low_order,
            :route_config_not_white_lindsey_low_order,
        )
    end

    missing_contract = Symbol[]
    !config.materializer_options_ready &&
        append!(missing_contract, config.missing_materializer_options)
    isnothing(parent_qw_basis_object) &&
        push!(missing_contract, :parent_qw_basis_object_handoff)
    isnothing(parent_axis_bundle_object) &&
        push!(missing_contract, :parent_axis_bundle_object_handoff)
    isnothing(axis_bundle_backend) &&
        push!(missing_contract, :axis_bundle_backend_provenance)
    if !isnothing(axis_bundle_backend) &&
       axis_bundle_backend != config.materializer_backend_requested
        push!(missing_contract, :axis_bundle_backend_mismatches_materializer_backend)
    end
    missing_contract = Tuple(unique(missing_contract))
    if !isempty(missing_contract)
        return blocked_probe(
            :blocked_missing_atom_growth_materializer_contract,
            :pending_route_configured_bond_aligned_diatomic_atom_growth_materializer_contract;
            missing_contract,
        )
    end

    expansion =
        isnothing(white_lindsey_expansion) ?
        coulomb_gaussian_expansion(doacc = false) :
        white_lindsey_expansion
    consumed_packet_kernel = isnothing(packet_kernel) ? :factorized_direct : packet_kernel
    materializer_options = (;
        nside = config.materializer_nside_requested,
        d = config.materializer_d_requested,
        reference_spacing = config.materializer_reference_spacing_requested,
        tail_spacing = config.materializer_tail_spacing_requested,
        gausslet_backend = config.materializer_backend_requested,
        axis_bundle_backend,
        packet_kernel = consumed_packet_kernel,
        term_coefficients_source = :coulomb_expansion_coefficients,
        plan_authority = :atom_growth,
    )

    try
        nside = config.materializer_nside_requested
        anatomy = _nested_bond_aligned_diatomic_atom_growth_anatomy(
            parent_qw_basis_object,
            parent_axis_bundle_object;
            bond_axis = config.bond_axis,
            protected_atom_side_count = nside,
        )
        construction_plan =
            _nested_bond_aligned_diatomic_atom_growth_construction_plan(anatomy)
        retention = _nested_resolve_complete_shell_retention(nside)
        protect_rows =
            _nested_diatomic_resolve_core_near_nucleus_protect_rows(:auto, nside)
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
                route_family = config.route_family,
            )
        materialization =
            _cartesian_materialize_atom_growth_complete_rectangular_shellification_low_order(
                scaffold,
                parent_qw_basis_object,
                parent_axis_bundle_object;
                term_coefficients = Float64.(expansion.coefficients),
                packet_kernel = consumed_packet_kernel,
            )
        sequence = materialization.sequence
        sequence_available = !isnothing(sequence)
        coverage_audit =
            sequence_available ?
            _nested_shell_sequence_contract_audit(
                sequence,
                Tuple(length.(construction_plan.anatomy.recipe.parent_box)),
            ) :
            nothing
        retained_dimension =
            sequence_available ? size(sequence.coefficient_matrix, 2) : nothing
        support_count =
            sequence_available ? length(sequence.support_indices) : nothing
        coverage_complete =
            !isnothing(coverage_audit) &&
            coverage_audit.full_parent_working_box &&
            coverage_audit.missing_row_count == 0 &&
            coverage_audit.ownership_unowned_row_count == 0 &&
            coverage_audit.ownership_multi_owned_row_count == 0
        materialized =
            materialization.status ==
            :materialized_supported_complete_rectangular_low_order
        basis_adapter =
            materialized ?
            _pqs_source_box_route_driver_diatomic_atom_growth_basis_adapter(
                materialization,
                construction_plan,
                scaffold,
                parent_qw_basis_object,
                config.materializer_backend_requested,
            ) :
            nothing
        basis_adapter_summary =
            isnothing(basis_adapter) ?
            nothing :
            _pqs_source_box_route_driver_route_configured_diatomic_basis_adapter_summary(
                basis_adapter,
            )
        basis_adapter_available =
            !isnothing(basis_adapter) &&
            basis_adapter.status ==
            :available_route_configured_diatomic_atom_growth_basis_adapter
        fixed_block_available =
            !isnothing(basis_adapter) && !isnothing(basis_adapter.fixed_block)
        ham_adapter =
            basis_adapter_available ?
            _pqs_source_box_route_driver_route_configured_diatomic_ham_adapter(
                basis_adapter,
                expansion;
                gausslet_backend = config.materializer_backend_requested,
                interaction_treatment = :ggt_nearest,
            ) :
            nothing
        ham_adapter_summary =
            isnothing(ham_adapter) ?
            nothing :
            _pqs_source_box_route_driver_route_configured_diatomic_ham_adapter_summary(
                ham_adapter,
            )
        ham_adapter_blocker =
            isnothing(ham_adapter) ?
            (
                basis_adapter_available ?
                :atom_growth_ham_operator_adapter_contract :
                isnothing(basis_adapter) ?
                :atom_growth_fixed_block_adapter_contract :
                basis_adapter.blocker
            ) :
            ham_adapter.blocker == :pending_route_configured_diatomic_ham_adapter_contract ?
            :atom_growth_ham_operator_adapter_contract :
            ham_adapter.blocker
        return (;
            object_kind =
                :route_configured_diatomic_atom_growth_materializer_probe,
            requested = true,
            status =
                materialized ?
                :materialized_route_configured_bond_aligned_diatomic_atom_growth_shellization :
                materialization.status,
            materialized,
            atom_growth_shellification_consumed = materialized,
            route_configured_diatomic_atom_growth_shellification_consumed =
                materialized,
            route_configured_shellization_consumed = materialized,
            route_configured_legacy_diatomic_source_consumed = false,
            blocker = materialized ? nothing : materialization.blocked_reason,
            missing_contract = (),
            parent_qw_basis_object_handoff_available,
            parent_axis_bundle_object_handoff_available,
            axis_bundle_backend_handoff_available,
            axis_bundle_backend_handoff = axis_bundle_backend,
            materializer_options,
            construction_plan,
            scaffold,
            materialization,
            sequence_available,
            retained_dimension,
            support_count,
            coverage_audit,
            coverage_complete,
            atom_growth_construction_plan_authority =
                scaffold.diagnostics.atom_growth_construction_plan_authority,
            active_source_authority =
                scaffold.diagnostics.active_source_authority ||
                materialization.active_source_authority,
            active_source_oracle_comparison_run =
                scaffold.diagnostics.active_source_oracle_comparison_run,
            route_behavior_changed =
                scaffold.diagnostics.materialization_behavior_changed ||
                materialization.route_behavior_changed,
            shellification_plan_path_used = true,
            calls_shellification_plan_materializer = true,
            fixed_block_available,
            fixed_block_status =
                fixed_block_available ?
                basis_adapter.fixed_block_status :
                isnothing(basis_adapter) ?
                :not_checked_missing_atom_growth_sequence :
                basis_adapter.fixed_block_status,
            basis_adapter,
            basis_adapter_summary,
            basis_adapter_status =
                isnothing(basis_adapter) ? :not_checked_missing_atom_growth_sequence :
                basis_adapter.status,
            basis_adapter_blocker =
                isnothing(basis_adapter) ? :atom_growth_fixed_block_adapter_contract :
                basis_adapter.blocker,
            final_integral_weights_status =
                isnothing(basis_adapter) ? :not_checked_missing_atom_growth_sequence :
                basis_adapter.final_integral_weights_status,
            ham_adapter,
            ham_adapter_summary,
            ham_adapter_status =
                isnothing(ham_adapter) ?
                :blocked_missing_route_configured_diatomic_atom_growth_basis_adapter :
                ham_adapter.status,
            ham_adapter_blocker,
            error_message = nothing,
        )
    catch error
        error isa ArgumentError || rethrow()
        return blocked_probe(
            :blocked_atom_growth_materializer_precondition,
            :atom_growth_materializer_precondition_failed;
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

function _pqs_source_box_route_driver_route_configured_diatomic_basis_adapter(
    materialization,
)
    if isnothing(materialization) || !materialization.materialized
        return (;
            object_kind = :route_configured_diatomic_basis_adapter,
            status = :blocked_missing_route_configured_diatomic_materialization,
            private_development_only = true,
            fixed_block = nothing,
            representation = nothing,
            final_integral_weights = nothing,
            retained_dimension = nothing,
            basis_metadata = nothing,
            grouping = nothing,
            label_status = :not_checked_no_materialization,
            grouping_status = :not_checked_no_materialization,
            final_integral_weights_status = :not_checked_no_materialization,
            missing_fields = (:route_configured_diatomic_materialization,),
            blocker = :missing_route_configured_diatomic_materialization,
        )
    elseif isnothing(materialization.source)
        return (;
            object_kind = :route_configured_diatomic_basis_adapter,
            status = :blocked_missing_route_configured_diatomic_source,
            private_development_only = true,
            fixed_block = nothing,
            representation = nothing,
            final_integral_weights = nothing,
            retained_dimension = materialization.retained_dimension,
            basis_metadata = nothing,
            grouping = nothing,
            label_status = :not_checked_no_source,
            grouping_status = :not_checked_no_source,
            final_integral_weights_status = :not_checked_no_source,
            missing_fields = (:route_configured_diatomic_source,),
            blocker = :missing_route_configured_diatomic_source,
        )
    end

    try
        source = materialization.source
        fixed_block = _nested_fixed_block(source)
        representation = basis_representation(fixed_block)
        final_integral_weights =
            _cartesian_bundle_integral_weights(fixed_block, representation)
        retained_dimension = representation.metadata.final_dimension
        labels = representation.metadata.basis_labels
        centers = representation.metadata.basis_centers
        weights_ready =
            length(final_integral_weights) == retained_dimension &&
            all(isfinite, final_integral_weights) &&
            all(>(0.0), final_integral_weights)
        labels_ready =
            length(labels) == retained_dimension &&
            all(!isempty, labels)
        centers_ready =
            size(centers) == (retained_dimension, 3) &&
            all(isfinite, centers)
        grouping = (;
            source_kind = :route_configured_bond_aligned_diatomic_source,
            shell_kind = representation.metadata.route_metadata.shell_kind,
            working_box_profile =
                representation.metadata.route_metadata.working_box_profile,
            support_count = representation.metadata.route_metadata.support_count,
            child_sequence_count = length(source.child_sequences),
            shared_shell_layer_count = length(source.shared_shell_layers),
            child_column_ranges = Tuple(source.child_column_ranges),
            midpoint_slab_column_range = source.midpoint_slab_column_range,
        )
        grouping_ready =
            grouping.child_sequence_count == length(grouping.child_column_ranges) &&
            grouping.support_count == length(fixed_block.support_indices)
        missing_fields = Symbol[]
        weights_ready || push!(
            missing_fields,
            :route_configured_diatomic_final_weight_contract,
        )
        labels_ready || push!(missing_fields, :route_configured_diatomic_basis_labels)
        centers_ready || push!(missing_fields, :route_configured_diatomic_basis_centers)
        grouping_ready || push!(missing_fields, :route_configured_diatomic_grouping)
        missing_fields = Tuple(missing_fields)

        return (;
            object_kind = :route_configured_diatomic_basis_adapter,
            status =
                isempty(missing_fields) ?
                :available_route_configured_diatomic_basis_adapter :
                :blocked_route_configured_diatomic_basis_adapter_contract,
            private_development_only = true,
            fixed_block,
            representation,
            final_integral_weights,
            retained_dimension,
            basis_metadata = (;
                basis_kind = representation.metadata.basis_kind,
                parent_kind = representation.metadata.parent_kind,
                axis_sharing = representation.metadata.axis_sharing,
                parent_axis_counts = representation.metadata.parent_axis_counts,
                parent_dimension = representation.metadata.parent_dimension,
                final_dimension = representation.metadata.final_dimension,
                label_count = length(labels),
                center_count = size(centers, 1),
            ),
            grouping,
            label_status =
                labels_ready ? :available_route_configured_diatomic_basis_labels :
                :blocked_route_configured_diatomic_basis_labels,
            grouping_status =
                grouping_ready ? :available_route_configured_diatomic_grouping :
                :blocked_route_configured_diatomic_grouping,
            final_integral_weights_status =
                weights_ready ?
                :available_retained_basis_integral_weights :
                :pending_route_configured_diatomic_final_weight_contract,
            missing_fields,
            blocker =
                isempty(missing_fields) ?
                nothing :
                :pending_route_configured_diatomic_basis_adapter_contract,
        )
    catch error
        error isa ArgumentError || rethrow()
        return (;
            object_kind = :route_configured_diatomic_basis_adapter,
            status = :blocked_route_configured_diatomic_basis_adapter_precondition,
            private_development_only = true,
            fixed_block = nothing,
            representation = nothing,
            final_integral_weights = nothing,
            retained_dimension = materialization.retained_dimension,
            basis_metadata = nothing,
            grouping = nothing,
            label_status = :not_checked_basis_adapter_precondition,
            grouping_status = :not_checked_basis_adapter_precondition,
            final_integral_weights_status =
                :not_checked_basis_adapter_precondition,
            missing_fields = (:route_configured_diatomic_basis_adapter_precondition,),
            blocker = :route_configured_diatomic_basis_adapter_precondition_failed,
            error_message = sprint(showerror, error),
        )
    end
end

function _pqs_source_box_route_driver_route_configured_diatomic_basis_adapter_summary(
    adapter,
)
    return (;
        object_kind = adapter.object_kind,
        status = adapter.status,
        private_development_only = adapter.private_development_only,
        retained_dimension = adapter.retained_dimension,
        basis_metadata = adapter.basis_metadata,
        grouping = adapter.grouping,
        label_status = adapter.label_status,
        grouping_status = adapter.grouping_status,
        final_integral_weights_status = adapter.final_integral_weights_status,
        final_integral_weight_count =
            isnothing(adapter.final_integral_weights) ?
            nothing :
            length(adapter.final_integral_weights),
        missing_fields = adapter.missing_fields,
        blocker = adapter.blocker,
    )
end

function _pqs_source_box_route_driver_diatomic_atom_growth_basis_adapter(
    materialization,
    construction_plan,
    scaffold,
    parent_qw_basis_object,
    gausslet_backend,
)
    if isnothing(materialization) || isnothing(materialization.sequence)
        return (;
            object_kind = :route_configured_diatomic_atom_growth_basis_adapter,
            status = :blocked_missing_route_configured_diatomic_atom_growth_sequence,
            private_development_only = true,
            fixed_block = nothing,
            fixed_block_status = :not_checked_missing_atom_growth_sequence,
            representation = nothing,
            final_integral_weights = nothing,
            retained_dimension = nothing,
            basis_metadata = nothing,
            grouping = nothing,
            label_status = :not_checked_missing_atom_growth_sequence,
            grouping_status = :not_checked_missing_atom_growth_sequence,
            final_integral_weights_status =
                :not_checked_missing_atom_growth_sequence,
            missing_fields = (:route_configured_diatomic_atom_growth_sequence,),
            blocker = :atom_growth_fixed_block_adapter_contract,
        )
    elseif isnothing(parent_qw_basis_object)
        return (;
            object_kind = :route_configured_diatomic_atom_growth_basis_adapter,
            status = :blocked_missing_atom_growth_parent_basis,
            private_development_only = true,
            fixed_block = nothing,
            fixed_block_status = :not_checked_missing_parent_basis,
            representation = nothing,
            final_integral_weights = nothing,
            retained_dimension = size(materialization.sequence.coefficient_matrix, 2),
            basis_metadata = nothing,
            grouping = nothing,
            label_status = :not_checked_missing_parent_basis,
            grouping_status = :not_checked_missing_parent_basis,
            final_integral_weights_status = :not_checked_missing_parent_basis,
            missing_fields = (:parent_qw_basis_object_handoff,),
            blocker = :atom_growth_fixed_block_adapter_contract,
        )
    end

    try
        sequence = materialization.sequence
        fixed_block = _nested_fixed_block(
            sequence,
            parent_qw_basis_object,
            gausslet_backend,
        )
        representation = basis_representation(fixed_block)
        final_integral_weights =
            _cartesian_bundle_integral_weights(fixed_block, representation)
        retained_dimension = representation.metadata.final_dimension
        labels = representation.metadata.basis_labels
        centers = representation.metadata.basis_centers
        route_metadata = representation.metadata.route_metadata
        weights_ready =
            length(final_integral_weights) == retained_dimension &&
            all(isfinite, final_integral_weights) &&
            all(>(0.0), final_integral_weights)
        labels_ready =
            length(labels) == retained_dimension &&
            all(!isempty, labels)
        centers_ready =
            size(centers) == (retained_dimension, 3) &&
            all(isfinite, centers)
        grouping = (;
            source_kind = :bond_aligned_diatomic_atom_growth_construction_plan,
            shell_kind = route_metadata.shell_kind,
            working_box_profile = route_metadata.working_box_profile,
            support_count = route_metadata.support_count,
            construction_region_order = Tuple(construction_plan.region_order),
            assembly_core_order = scaffold.assembly_core_order,
            assembly_shell_order = scaffold.assembly_shell_order,
            outer_mismatch_boundary_slab_set_count =
                length(scaffold.outer_mismatch_boundary_slab_sets),
            child_column_ranges = Tuple(materialization.assembly.child_column_ranges),
            contact_cap_column_range =
                materialization.assembly.contact_cap_column_range,
            shared_shell_column_ranges =
                Tuple(materialization.assembly.shared_shell_column_ranges),
        )
        grouping_ready =
            grouping.support_count == length(fixed_block.support_indices) &&
            grouping.support_count == length(sequence.support_indices)
        missing_fields = Symbol[]
        weights_ready ||
            push!(missing_fields, :atom_growth_final_integral_weight_contract)
        labels_ready ||
            push!(missing_fields, :atom_growth_basis_representation_contract)
        centers_ready ||
            push!(missing_fields, :atom_growth_basis_representation_contract)
        grouping_ready ||
            push!(missing_fields, :atom_growth_basis_representation_contract)
        missing_fields = Tuple(unique(missing_fields))
        blocker =
            isempty(missing_fields) ?
            nothing :
            in(:atom_growth_final_integral_weight_contract, missing_fields) ?
            :atom_growth_final_integral_weight_contract :
            :atom_growth_basis_representation_contract

        return (;
            object_kind = :route_configured_diatomic_atom_growth_basis_adapter,
            status =
                isempty(missing_fields) ?
                :available_route_configured_diatomic_atom_growth_basis_adapter :
                :blocked_route_configured_diatomic_atom_growth_basis_adapter_contract,
            private_development_only = true,
            fixed_block,
            fixed_block_status = :available_route_configured_diatomic_atom_growth_fixed_block,
            representation,
            final_integral_weights,
            retained_dimension,
            basis_metadata = (;
                basis_kind = representation.metadata.basis_kind,
                parent_kind = representation.metadata.parent_kind,
                axis_sharing = representation.metadata.axis_sharing,
                parent_axis_counts = representation.metadata.parent_axis_counts,
                parent_dimension = representation.metadata.parent_dimension,
                final_dimension = representation.metadata.final_dimension,
                label_count = length(labels),
                center_count = size(centers, 1),
            ),
            grouping,
            label_status =
                labels_ready ?
                :available_route_configured_diatomic_atom_growth_basis_labels :
                :blocked_atom_growth_basis_representation_contract,
            grouping_status =
                grouping_ready ?
                :available_route_configured_diatomic_atom_growth_grouping :
                :blocked_atom_growth_basis_representation_contract,
            final_integral_weights_status =
                weights_ready ?
                :available_retained_basis_integral_weights :
                :pending_atom_growth_final_integral_weight_contract,
            missing_fields,
            blocker,
        )
    catch error
        error isa ArgumentError || rethrow()
        return (;
            object_kind = :route_configured_diatomic_atom_growth_basis_adapter,
            status = :blocked_atom_growth_fixed_block_adapter_precondition,
            private_development_only = true,
            fixed_block = nothing,
            fixed_block_status = :blocked_atom_growth_fixed_block_adapter_precondition,
            representation = nothing,
            final_integral_weights = nothing,
            retained_dimension =
                !isnothing(materialization) && !isnothing(materialization.sequence) ?
                size(materialization.sequence.coefficient_matrix, 2) :
                nothing,
            basis_metadata = nothing,
            grouping = nothing,
            label_status = :not_checked_fixed_block_adapter_precondition,
            grouping_status = :not_checked_fixed_block_adapter_precondition,
            final_integral_weights_status =
                :not_checked_fixed_block_adapter_precondition,
            missing_fields = (:atom_growth_fixed_block_adapter_contract,),
            blocker = :atom_growth_fixed_block_adapter_contract,
            error_message = sprint(showerror, error),
        )
    end
end

function _pqs_source_box_route_driver_route_configured_diatomic_ham_adapter(
    basis_adapter,
    expansion;
    gausslet_backend,
    interaction_treatment::Symbol = :ggt_nearest,
    nuclear_term_storage::Symbol = :by_center,
)
    mwg_ida_treatments = (:mwg, :ida, :mwg_ida, :ida_mwg)
    available_basis_adapter_statuses = (
        :available_route_configured_diatomic_basis_adapter,
        :available_route_configured_diatomic_atom_growth_basis_adapter,
    )
    if !(basis_adapter.status in available_basis_adapter_statuses)
        return (;
            object_kind = :route_configured_diatomic_ham_adapter,
            status = :blocked_missing_route_configured_diatomic_basis_adapter,
            private_development_only = true,
            operators = nothing,
            retained_dimension = basis_adapter.retained_dimension,
            matrix_sizes = nothing,
            final_integral_weights_status = basis_adapter.final_integral_weights_status,
            nuclear_metadata_status = :not_checked_missing_basis_adapter,
            operator_payload_status = :not_checked_missing_basis_adapter,
            interaction_status = :not_checked_missing_basis_adapter,
            interaction_treatment_requested = interaction_treatment,
            interaction_treatment,
            missing_fields = (:route_configured_diatomic_basis_adapter,),
            blocker = :missing_route_configured_diatomic_basis_adapter,
        )
    elseif interaction_treatment in mwg_ida_treatments
        return (;
            object_kind = :route_configured_diatomic_ham_adapter,
            status =
                :blocked_route_configured_diatomic_ham_interaction_treatment,
            private_development_only = true,
            operators = nothing,
            retained_dimension = basis_adapter.retained_dimension,
            matrix_sizes = nothing,
            final_integral_weights_status = basis_adapter.final_integral_weights_status,
            nuclear_metadata_status = :not_checked_blocked_interaction_treatment,
            operator_payload_status = :blocked_route_configured_diatomic_operator_payload,
            interaction_status =
                :pending_route_configured_diatomic_mwg_operator_support,
            interaction_treatment_requested = interaction_treatment,
            interaction_treatment,
            gausslet_backend,
            nuclear_charges = nothing,
            nuclear_term_storage = nothing,
            nuclear_one_body_by_center_count = 0,
            missing_fields = (:pending_route_configured_diatomic_mwg_operator_support,),
            blocker = :pending_route_configured_diatomic_mwg_operator_support,
        )
    elseif isnothing(expansion)
        return (;
            object_kind = :route_configured_diatomic_ham_adapter,
            status = :blocked_missing_route_configured_diatomic_expansion,
            private_development_only = true,
            operators = nothing,
            retained_dimension = basis_adapter.retained_dimension,
            matrix_sizes = nothing,
            final_integral_weights_status = basis_adapter.final_integral_weights_status,
            nuclear_metadata_status = :not_checked_missing_expansion,
            operator_payload_status = :not_checked_missing_expansion,
            interaction_status = :not_checked_missing_expansion,
            interaction_treatment_requested = interaction_treatment,
            interaction_treatment,
            missing_fields = (:coulomb_expansion,),
            blocker = :missing_route_configured_diatomic_expansion,
        )
    end

    try
        operators = ordinary_cartesian_qiu_white_operators(
            basis_adapter.fixed_block;
            expansion,
            gausslet_backend,
            interaction_treatment,
            nuclear_term_storage,
        )
        retained_dimension = size(operators.overlap, 1)
        matrix_sizes = (;
            overlap = size(operators.overlap),
            one_body_hamiltonian = size(operators.one_body_hamiltonian),
            interaction_matrix = size(operators.interaction_matrix),
        )
        expected_matrix_size = (retained_dimension, retained_dimension)
        final_integral_weights =
            _cartesian_bundle_integral_weights(operators, basis_adapter.representation)
        matrix_size_ready =
            all(matrix_size -> matrix_size == expected_matrix_size, values(matrix_sizes))
        finite_ready =
            all(isfinite, operators.overlap) &&
            all(isfinite, operators.one_body_hamiltonian) &&
            all(isfinite, operators.interaction_matrix)
        weights_ready =
            length(final_integral_weights) == retained_dimension &&
            all(isfinite, final_integral_weights) &&
            all(>(0.0), final_integral_weights)
        parent_nuclear_charges = basis_adapter.fixed_block.parent_basis.nuclear_charges
        expected_nuclear_center_count = length(parent_nuclear_charges)
        nuclear_metadata_ready =
            !isnothing(operators.nuclear_charges) &&
            length(operators.nuclear_charges) == expected_nuclear_center_count &&
            operators.nuclear_term_storage == :by_center &&
            !isnothing(operators.nuclear_one_body_by_center) &&
            length(operators.nuclear_one_body_by_center) ==
            length(operators.nuclear_charges)
        missing_fields = Symbol[]
        matrix_size_ready || push!(missing_fields, :route_configured_diatomic_ham_matrix_sizes)
        finite_ready || push!(missing_fields, :route_configured_diatomic_ham_finite_matrices)
        weights_ready || push!(
            missing_fields,
            :route_configured_diatomic_ham_final_integral_weights,
        )
        nuclear_metadata_ready || push!(
            missing_fields,
            :route_configured_diatomic_ham_nuclear_metadata,
        )
        missing_fields = Tuple(missing_fields)

        return (;
            object_kind = :route_configured_diatomic_ham_adapter,
            status =
                isempty(missing_fields) ?
                :available_route_configured_diatomic_ham_adapter :
                :blocked_route_configured_diatomic_ham_adapter_contract,
            private_development_only = true,
            operators,
            retained_dimension,
            matrix_sizes,
            final_integral_weights_status =
                weights_ready ?
                :available_retained_basis_integral_weights :
                :pending_route_configured_diatomic_ham_final_integral_weights,
            nuclear_metadata_status =
                nuclear_metadata_ready ?
                :available_route_configured_diatomic_nuclear_metadata :
                :pending_route_configured_diatomic_nuclear_metadata,
            operator_payload_status =
                matrix_size_ready && finite_ready ?
                :available_route_configured_diatomic_operator_payload :
                :pending_route_configured_diatomic_operator_payload,
            interaction_status =
                matrix_size_ready && finite_ready ?
                :available_route_configured_diatomic_density_density_interaction_matrix :
                :pending_route_configured_diatomic_density_density_builder,
            interaction_treatment_requested = interaction_treatment,
            interaction_treatment = operators.interaction_treatment,
            gausslet_backend = operators.gausslet_backend,
            nuclear_charges = Tuple(operators.nuclear_charges),
            nuclear_term_storage = operators.nuclear_term_storage,
            nuclear_one_body_by_center_count =
                isnothing(operators.nuclear_one_body_by_center) ?
                0 :
                length(operators.nuclear_one_body_by_center),
            missing_fields,
            blocker =
                isempty(missing_fields) ?
                nothing :
                :pending_route_configured_diatomic_ham_adapter_contract,
        )
    catch error
        error isa ArgumentError || rethrow()
        return (;
            object_kind = :route_configured_diatomic_ham_adapter,
            status = :blocked_route_configured_diatomic_ham_adapter_precondition,
            private_development_only = true,
            operators = nothing,
            retained_dimension = basis_adapter.retained_dimension,
            matrix_sizes = nothing,
            final_integral_weights_status =
                :not_checked_ham_adapter_precondition,
            nuclear_metadata_status = :not_checked_ham_adapter_precondition,
            operator_payload_status = :not_checked_ham_adapter_precondition,
            interaction_status = :not_checked_ham_adapter_precondition,
            interaction_treatment_requested = interaction_treatment,
            interaction_treatment,
            missing_fields = (:route_configured_diatomic_ham_adapter_precondition,),
            blocker = :route_configured_diatomic_ham_adapter_precondition_failed,
            error_message = sprint(showerror, error),
        )
    end
end

function _pqs_source_box_route_driver_route_configured_diatomic_ham_adapter_summary(
    adapter,
)
    return (;
        object_kind = adapter.object_kind,
        status = adapter.status,
        private_development_only = adapter.private_development_only,
        retained_dimension = adapter.retained_dimension,
        matrix_sizes = adapter.matrix_sizes,
        final_integral_weights_status = adapter.final_integral_weights_status,
        nuclear_metadata_status = adapter.nuclear_metadata_status,
        operator_payload_status = adapter.operator_payload_status,
        interaction_status = adapter.interaction_status,
        interaction_treatment_requested =
            hasproperty(adapter, :interaction_treatment_requested) ?
            adapter.interaction_treatment_requested :
            nothing,
        interaction_treatment =
            hasproperty(adapter, :interaction_treatment) ?
            adapter.interaction_treatment :
            nothing,
        gausslet_backend =
            hasproperty(adapter, :gausslet_backend) ? adapter.gausslet_backend : nothing,
        nuclear_charges =
            hasproperty(adapter, :nuclear_charges) ? adapter.nuclear_charges : nothing,
        nuclear_term_storage =
            hasproperty(adapter, :nuclear_term_storage) ?
            adapter.nuclear_term_storage :
            nothing,
        nuclear_one_body_by_center_count =
            hasproperty(adapter, :nuclear_one_body_by_center_count) ?
            adapter.nuclear_one_body_by_center_count :
            nothing,
        missing_fields = adapter.missing_fields,
        blocker = adapter.blocker,
    )
end

function _pqs_source_box_route_driver_low_order_shellization_policy(
    requested_policy,
    probe_route_configured_diatomic_atom_growth_materializer::Bool,
)
    supported_policies = (
        :legacy_diatomic_source,
        :atom_growth_complete_rectangular,
        :terminal_cartesian_shellification_geometry,
    )
    explicit_policy_requested = !isnothing(requested_policy)
    resolved_policy =
        explicit_policy_requested ?
        requested_policy :
        probe_route_configured_diatomic_atom_growth_materializer ?
        :atom_growth_complete_rectangular :
        :legacy_diatomic_source
    policy_source =
        explicit_policy_requested ?
        :explicit_low_order_shellization_policy :
        probe_route_configured_diatomic_atom_growth_materializer ?
        :probe_route_configured_diatomic_atom_growth_materializer_alias :
        :default_legacy_diatomic_source
    supported = resolved_policy in supported_policies
    conflict =
        explicit_policy_requested &&
        probe_route_configured_diatomic_atom_growth_materializer &&
        resolved_policy != :atom_growth_complete_rectangular
    status =
        !supported ?
        :blocked_unsupported_low_order_shellization_policy :
        conflict ?
        :blocked_conflicting_low_order_shellization_policy :
        :available_low_order_shellization_policy
    blocker =
        status == :available_low_order_shellization_policy ?
        nothing :
        status

    return (;
        low_order_shellization_policy_requested = requested_policy,
        low_order_shellization_policy_resolved = resolved_policy,
        low_order_shellization_policy_source = policy_source,
        low_order_shellization_policy_status = status,
        low_order_shellization_policy_blocker = blocker,
    )
end

function _pqs_source_box_route_driver_materialization(
    report;
    materialize_route::Bool = false,
    probe_route_configured_one_center_materializer::Bool = false,
    probe_route_configured_diatomic_atom_growth_materializer::Bool = false,
    low_order_shellization_policy = nothing,
    save_basis_artifact::Bool = false,
    save_ham_artifact::Bool = false,
    basisfile::AbstractString = "cartesian_nesting_route_driver_basis_bundle.jld2",
    hamfile::AbstractString = "cartesian_nesting_route_driver_ham_bundle.jld2",
    materializer_backend = nothing,
    materializer_nside = nothing,
    route_configured_diatomic_ham_interaction_treatment::Symbol = :ggt_nearest,
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
    low_order_shellization_policy_contract =
        _pqs_source_box_route_driver_low_order_shellization_policy(
            low_order_shellization_policy,
            probe_route_configured_diatomic_atom_growth_materializer,
        )
    low_order_shellization_policy_requested =
        low_order_shellization_policy_contract.low_order_shellization_policy_requested
    low_order_shellization_policy_resolved =
        low_order_shellization_policy_contract.low_order_shellization_policy_resolved
    low_order_shellization_policy_source =
        low_order_shellization_policy_contract.low_order_shellization_policy_source
    low_order_shellization_policy_status =
        low_order_shellization_policy_contract.low_order_shellization_policy_status
    low_order_shellization_policy_blocker =
        low_order_shellization_policy_contract.low_order_shellization_policy_blocker
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
        route_configured_system_classification == :bond_aligned_diatomic &&
        low_order_shellization_policy_status == :available_low_order_shellization_policy &&
        low_order_shellization_policy_resolved == :legacy_diatomic_source
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
    route_configured_diatomic_atom_growth_materializer_requested =
        materialize_route &&
        route_family == :white_lindsey_low_order &&
        route_configured_system_classification == :bond_aligned_diatomic &&
        low_order_shellization_policy_status == :available_low_order_shellization_policy &&
        low_order_shellization_policy_resolved == :atom_growth_complete_rectangular
    route_configured_diatomic_atom_growth_materializer_probe_requested_input =
        low_order_shellization_policy_status == :available_low_order_shellization_policy &&
        (
            probe_route_configured_diatomic_atom_growth_materializer ||
            route_configured_diatomic_atom_growth_materializer_requested
        )
    route_configured_diatomic_atom_growth_materializer_probe =
        _pqs_source_box_route_driver_diatomic_atom_growth_materializer_probe(
            route_configured_materializer_config;
            probe_route_configured_diatomic_atom_growth_materializer =
                route_configured_diatomic_atom_growth_materializer_probe_requested_input,
            route_materializer_payload,
            white_lindsey_expansion,
            packet_kernel = route_configured_diatomic_packet_kernel,
        )
    route_configured_diatomic_atom_growth_materializer_probe_requested =
        route_configured_diatomic_atom_growth_materializer_probe.requested
    route_configured_diatomic_atom_growth_materializer_probe_status =
        route_configured_diatomic_atom_growth_materializer_probe.status
    route_configured_diatomic_atom_growth_materializer_probe_materialized =
        route_configured_diatomic_atom_growth_materializer_probe.materialized
    route_configured_diatomic_atom_growth_materializer_probe_consumed =
        route_configured_diatomic_atom_growth_materializer_probe.atom_growth_shellification_consumed
    route_configured_diatomic_atom_growth_shellification_consumed =
        route_configured_diatomic_atom_growth_materializer_probe.route_configured_diatomic_atom_growth_shellification_consumed
    route_configured_diatomic_atom_growth_materializer_probe_blocker =
        route_configured_diatomic_atom_growth_materializer_probe.blocker
    route_configured_diatomic_atom_growth_materializer_missing_contract =
        route_configured_diatomic_atom_growth_materializer_probe.missing_contract
    route_configured_diatomic_atom_growth_sequence_available =
        route_configured_diatomic_atom_growth_materializer_probe.sequence_available
    route_configured_diatomic_atom_growth_retained_dimension =
        route_configured_diatomic_atom_growth_materializer_probe.retained_dimension
    route_configured_diatomic_atom_growth_support_count =
        route_configured_diatomic_atom_growth_materializer_probe.support_count
    route_configured_diatomic_atom_growth_coverage_complete =
        route_configured_diatomic_atom_growth_materializer_probe.coverage_complete
    route_configured_diatomic_atom_growth_active_source_authority =
        route_configured_diatomic_atom_growth_materializer_probe.active_source_authority
    route_configured_diatomic_atom_growth_plan_authority =
        route_configured_diatomic_atom_growth_materializer_probe.atom_growth_construction_plan_authority
    route_configured_diatomic_atom_growth_calls_shellification_plan_materializer =
        route_configured_diatomic_atom_growth_materializer_probe.calls_shellification_plan_materializer
    route_configured_diatomic_atom_growth_fixed_block_available =
        route_configured_diatomic_atom_growth_materializer_probe.fixed_block_available
    route_configured_diatomic_atom_growth_fixed_block_status =
        route_configured_diatomic_atom_growth_materializer_probe.fixed_block_status
    route_configured_diatomic_atom_growth_basis_adapter_status =
        route_configured_diatomic_atom_growth_materializer_probe.basis_adapter_status
    route_configured_diatomic_atom_growth_basis_adapter_blocker =
        route_configured_diatomic_atom_growth_materializer_probe.basis_adapter_blocker
    route_configured_diatomic_atom_growth_final_integral_weights_status =
        route_configured_diatomic_atom_growth_materializer_probe.final_integral_weights_status
    route_configured_diatomic_atom_growth_ham_adapter_status =
        route_configured_diatomic_atom_growth_materializer_probe.ham_adapter_status
    route_configured_diatomic_atom_growth_ham_adapter_blocker =
        route_configured_diatomic_atom_growth_materializer_probe.ham_adapter_blocker
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
        route_configured_diatomic_atom_growth_materializer_probe_materialized ?
        route_configured_diatomic_atom_growth_materializer_probe.materializer_options :
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
        low_order_shellization_policy_requested,
        low_order_shellization_policy_resolved,
        low_order_shellization_policy_source,
        low_order_shellization_policy_status,
        low_order_shellization_policy_blocker,
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
    route_configured_diatomic_atom_growth_materializer_contract = (;
        route_configured_diatomic_atom_growth_materializer_probe_requested,
        route_configured_diatomic_atom_growth_materializer_probe_status,
        route_configured_diatomic_atom_growth_materializer_probe_materialized,
        route_configured_diatomic_atom_growth_materializer_probe_consumed,
        route_configured_diatomic_atom_growth_shellification_consumed,
        route_configured_diatomic_atom_growth_materializer_probe_blocker,
        route_configured_diatomic_atom_growth_materializer_missing_contract,
        route_configured_diatomic_atom_growth_sequence_available,
        route_configured_diatomic_atom_growth_retained_dimension,
        route_configured_diatomic_atom_growth_support_count,
        route_configured_diatomic_atom_growth_coverage_complete,
        route_configured_diatomic_atom_growth_active_source_authority,
        route_configured_diatomic_atom_growth_plan_authority,
        route_configured_diatomic_atom_growth_calls_shellification_plan_materializer,
        route_configured_diatomic_atom_growth_fixed_block_available,
        route_configured_diatomic_atom_growth_fixed_block_status,
        route_configured_diatomic_atom_growth_basis_adapter_status,
        route_configured_diatomic_atom_growth_basis_adapter_blocker,
        route_configured_diatomic_atom_growth_final_integral_weights_status,
        route_configured_diatomic_atom_growth_ham_adapter_status,
        route_configured_diatomic_atom_growth_ham_adapter_blocker,
    )

    if !materialize_route
        return (;
            object_kind = :cartesian_nesting_route_driver_materialization,
            route_family,
            private_development_only = true,
            materialize_route_requested = false,
            save_basis_artifact_requested = save_basis_artifact,
            save_ham_artifact_requested = save_ham_artifact,
            route_configured_diatomic_ham_interaction_treatment_requested = nothing,
            route_configured_diatomic_ham_interaction_treatment_consumed = nothing,
            route_configured_diatomic_ham_interaction_treatment_status =
                :not_applicable,
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
            route_configured_diatomic_atom_growth_materializer_probe,
            route_configured_diatomic_materializer_contract...,
            route_configured_diatomic_atom_growth_materializer_contract...,
            route_configured_materializer_contract...,
            shellization_summary = nothing,
            shellization_summary_available = false,
            shellization_source =
                route_family == :white_lindsey_low_order ?
                :white_lindsey_one_center_seed_not_materialized :
                nothing,
            route_configured_shellization_consumed = false,
            route_configured_legacy_diatomic_source_consumed = false,
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

    if low_order_shellization_policy_status != :available_low_order_shellization_policy
        return (;
            object_kind = :cartesian_nesting_route_driver_materialization,
            route_family,
            private_development_only = true,
            materialize_route_requested = true,
            save_basis_artifact_requested = save_basis_artifact,
            save_ham_artifact_requested = save_ham_artifact,
            route_configured_diatomic_ham_interaction_treatment_requested = nothing,
            route_configured_diatomic_ham_interaction_treatment_consumed = nothing,
            route_configured_diatomic_ham_interaction_treatment_status =
                :not_applicable,
            status = low_order_shellization_policy_status,
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
            route_configured_diatomic_atom_growth_materializer_probe,
            route_configured_diatomic_materializer_contract...,
            route_configured_diatomic_atom_growth_materializer_contract...,
            route_configured_materializer_contract...,
            shellization_summary = nothing,
            shellization_summary_available = false,
            shellization_source = :blocked_low_order_shellization_policy,
            route_configured_shellization_consumed = false,
            route_configured_legacy_diatomic_source_consumed = false,
            materialized_shellization_stage = :blocked_low_order_shellization_policy,
            seed_materialization_status = :not_applicable,
            retained_dimension = report.retained_dimension,
            final_integral_weights_status = :not_checked_low_order_policy_blocked,
            one_body_operator_status = :not_checked_low_order_policy_blocked,
            basis_bundle_export_status = :not_requested,
            basis_artifact_status =
                save_basis_artifact ?
                :not_written_invalid_low_order_shellization_policy :
                :not_requested,
            basis_artifact_written = false,
            basisfile,
            basis_artifact_path = nothing,
            basis_export_blocker =
                save_basis_artifact ? low_order_shellization_policy_blocker : nothing,
            ham_preflight_status = :not_checked_low_order_policy_blocked,
            ham_missing_builder =
                save_ham_artifact ? low_order_shellization_policy_blocker : nothing,
            ham_operator_payload_status = :not_checked_low_order_policy_blocked,
            ham_interaction_status = :not_checked_low_order_policy_blocked,
            ham_bundle_export_status = :not_requested,
            ham_artifact_status =
                save_ham_artifact ?
                :not_written_invalid_low_order_shellization_policy :
                :not_requested,
            ham_artifact_written = false,
            hamfile,
            ham_export_blocker =
                save_ham_artifact ? low_order_shellization_policy_blocker : nothing,
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
        route_configured_diatomic_atom_growth_materialization_requested =
            route_configured_diatomic_atom_growth_materializer_probe_requested
        if route_configured_diatomic_atom_growth_materialization_requested
            atom_growth_materialization_context = (;
                report,
                route_family,
                save_basis_artifact,
                save_ham_artifact,
                basisfile,
                hamfile,
                route_configured_diatomic_ham_interaction_treatment,
                route_configured_shellization_request,
                route_configured_shellization_request_status,
                route_configured_system_classification,
                route_configured_system_classification_status,
                route_configured_bond_axis,
                route_configured_shellization_plan,
                route_configured_shellization_plan_status,
                route_configured_shellization_planning_status,
                route_configured_shellization_planning_family,
                route_configured_midpoint_slab_status,
                route_configured_shellization_helper_map,
                route_configured_shellization_helper_map_status,
                route_configured_primary_planned_helper,
                route_configured_missing_input_count,
                route_configured_helper_map_blocker,
                route_configured_input_readiness,
                route_configured_input_readiness_status,
                route_configured_available_fact_count,
                route_configured_materializer_missing_input_count,
                route_configured_input_readiness_blocker,
                route_configured_materializer_config,
                route_configured_materializer_config_status,
                route_configured_materializer_config_planning_family,
                route_configured_materializer_config_pending_input_count,
                route_configured_one_center_materializer_probe,
                route_configured_one_center_materializer_probe_requested,
                route_configured_one_center_materializer_probe_status,
                route_configured_one_center_materializer_probe_materialized,
                route_configured_one_center_materializer_probe_consumed,
                route_configured_one_center_materializer_probe_blocker,
                route_configured_diatomic_atom_growth_materializer_probe,
                route_configured_diatomic_atom_growth_materializer_probe_consumed,
                route_configured_diatomic_atom_growth_materializer_probe_blocker,
                route_configured_diatomic_atom_growth_basis_adapter_blocker,
                route_configured_diatomic_atom_growth_final_integral_weights_status,
                route_configured_diatomic_atom_growth_ham_adapter_status,
                route_configured_diatomic_atom_growth_ham_adapter_blocker,
                route_configured_diatomic_materializer_contract,
                route_configured_diatomic_atom_growth_materializer_contract,
                route_configured_materializer_contract,
            )
            return _pqs_source_box_route_driver_diatomic_atom_growth_materialization(
                atom_growth_materialization_context,
            )
        end
        if use_route_configured_diatomic_shellization
            diatomic_materialization =
                route_configured_diatomic_materializer_probe.materialization
            diatomic_basis_adapter =
                _pqs_source_box_route_driver_route_configured_diatomic_basis_adapter(
                    diatomic_materialization,
                )
            diatomic_basis_adapter_summary =
                _pqs_source_box_route_driver_route_configured_diatomic_basis_adapter_summary(
                    diatomic_basis_adapter,
                )
            diatomic_basis_adapter_available =
                diatomic_basis_adapter.status ==
                :available_route_configured_diatomic_basis_adapter
            diatomic_ham_adapter =
                save_ham_artifact ?
                _pqs_source_box_route_driver_route_configured_diatomic_ham_adapter(
                    diatomic_basis_adapter,
                    white_lindsey_expansion;
                    gausslet_backend =
                        route_configured_materializer_backend_requested,
                    interaction_treatment =
                        route_configured_diatomic_ham_interaction_treatment,
                ) :
                nothing
            diatomic_ham_adapter_summary =
                isnothing(diatomic_ham_adapter) ?
                nothing :
                _pqs_source_box_route_driver_route_configured_diatomic_ham_adapter_summary(
                    diatomic_ham_adapter,
                )
            diatomic_ham_adapter_available =
                !isnothing(diatomic_ham_adapter) &&
                diatomic_ham_adapter.status ==
                :available_route_configured_diatomic_ham_adapter
            diatomic_ham_interaction_treatment_consumed =
                diatomic_ham_adapter_available ?
                diatomic_ham_adapter.interaction_treatment :
                nothing
            diatomic_ham_interaction_treatment_status =
                diatomic_ham_adapter_available ?
                :available_route_configured_diatomic_ham_interaction_treatment :
                save_ham_artifact && !isnothing(diatomic_ham_adapter) ?
                diatomic_ham_adapter.interaction_status :
                :not_requested
            shellization_summary = diatomic_materialization.shellization_summary
            shellization_summary_available = !isnothing(shellization_summary)
            basis_artifact_status =
                save_basis_artifact ?
                (
                    diatomic_basis_adapter_available ?
                    :written_route_configured_diatomic_basis_only_bundle :
                    :not_written_route_configured_diatomic_basis_adapter_blocked
                ) :
                :not_requested
            ham_artifact_status =
                save_ham_artifact ?
                (
                    diatomic_ham_adapter_available ?
                    :written_route_configured_diatomic_ham_bundle :
                    :not_written_route_configured_diatomic_ham_adapter_blocked
                ) :
                :not_requested
            diatomic_ham_adapter_blocker =
                isnothing(diatomic_ham_adapter) ? nothing : diatomic_ham_adapter.blocker
            diatomic_ham_bundle_export_status =
                diatomic_ham_adapter_available ?
                :available_route_configured_diatomic_ham_bundle_payload :
                save_ham_artifact ?
                something(
                    diatomic_ham_adapter_blocker,
                    :pending_route_configured_diatomic_ham_export,
                ) :
                :not_requested
            ham_export_blocker =
                diatomic_ham_adapter_available || !save_ham_artifact ?
                nothing :
                something(
                    diatomic_ham_adapter_blocker,
                    :pending_route_configured_diatomic_ham_export,
                )
            basis_companion_ham_artifact_status =
                save_ham_artifact ?
                (
                    diatomic_ham_adapter_available ?
                    :companion_route_configured_diatomic_ham_artifact_ready :
                    something(
                        diatomic_ham_adapter_blocker,
                        :pending_route_configured_diatomic_ham_export,
                    )
                ) :
                :not_requested
            basis_companion_ham_export_status =
                save_ham_artifact ? diatomic_ham_bundle_export_status : :not_requested
            basis_companion_ham_export_blocker =
                save_ham_artifact ? ham_export_blocker : nothing
            basis_artifact_written = false
            if save_basis_artifact && diatomic_basis_adapter_available
                write_cartesian_basis_bundle_jld2(
                    basisfile,
                    diatomic_basis_adapter.fixed_block;
                    include_ham = false,
                    meta = (;
                        route_family,
                        route_kind = report.recipe_metadata.route_kind,
                        benchmark_role = report.recipe_metadata.benchmark_role,
                        materialized_report_kind = diatomic_materialization.object_kind,
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
                        route_configured_diatomic_atom_growth_materializer_contract...,
                        route_configured_materializer_contract...,
                        route_configured_diatomic_basis_adapter_status =
                            diatomic_basis_adapter.status,
                        route_configured_diatomic_basis_adapter_retained_dimension =
                            diatomic_basis_adapter.retained_dimension,
                        route_configured_diatomic_basis_adapter_final_integral_weights_status =
                            diatomic_basis_adapter.final_integral_weights_status,
                        route_configured_diatomic_basis_adapter_label_status =
                            diatomic_basis_adapter.label_status,
                        route_configured_diatomic_basis_adapter_grouping_status =
                            diatomic_basis_adapter.grouping_status,
                        route_configured_diatomic_ham_interaction_treatment_requested =
                            route_configured_diatomic_ham_interaction_treatment,
                        route_configured_diatomic_ham_interaction_treatment_consumed =
                            diatomic_ham_interaction_treatment_consumed,
                        route_configured_diatomic_ham_interaction_treatment_status =
                            diatomic_ham_interaction_treatment_status,
                        shellization_summary_available,
                        shellization_source =
                            :route_configured_bond_aligned_diatomic_source,
                        route_configured_shellization_consumed = true,
                        route_configured_legacy_diatomic_source_consumed = true,
                        materialized_shellization_stage =
                            shellization_summary.shellization_stage,
                        seed_materialization_status =
                            :not_seed_route_configured_diatomic_shellization,
                        export_status = :basis_only,
                        basis_export_status =
                            :supported_route_configured_diatomic_basis_only_fixed_block,
                        ham_export_status =
                            :artifact_local_basis_only_no_ham_payload,
                        ham_export_blocker = nothing,
                        companion_ham_artifact_requested = save_ham_artifact,
                        companion_ham_artifact_status =
                            basis_companion_ham_artifact_status,
                        companion_ham_export_status =
                            basis_companion_ham_export_status,
                        companion_ham_export_blocker =
                            basis_companion_ham_export_blocker,
                        private_development_only = true,
                    ),
                )
                basis_artifact_written = true
            end
            ham_artifact_written = false
            if save_ham_artifact && diatomic_ham_adapter_available
                write_cartesian_basis_bundle_jld2(
                    hamfile,
                    diatomic_ham_adapter.operators;
                    include_ham = true,
                    meta = (;
                        route_family,
                        route_kind = report.recipe_metadata.route_kind,
                        benchmark_role = report.recipe_metadata.benchmark_role,
                        materialized_report_kind = diatomic_materialization.object_kind,
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
                        route_configured_diatomic_atom_growth_materializer_contract...,
                        route_configured_materializer_contract...,
                        route_configured_diatomic_basis_adapter_status =
                            diatomic_basis_adapter.status,
                        route_configured_diatomic_basis_adapter_retained_dimension =
                            diatomic_basis_adapter.retained_dimension,
                        route_configured_diatomic_basis_adapter_final_integral_weights_status =
                            diatomic_basis_adapter.final_integral_weights_status,
                        route_configured_diatomic_ham_adapter_status =
                            diatomic_ham_adapter.status,
                        route_configured_diatomic_ham_adapter_operator_payload_status =
                            diatomic_ham_adapter.operator_payload_status,
                        route_configured_diatomic_ham_adapter_interaction_status =
                            diatomic_ham_adapter.interaction_status,
                        route_configured_diatomic_ham_adapter_nuclear_metadata_status =
                            diatomic_ham_adapter.nuclear_metadata_status,
                        route_configured_diatomic_ham_interaction_treatment_requested =
                            route_configured_diatomic_ham_interaction_treatment,
                        route_configured_diatomic_ham_interaction_treatment_consumed =
                            diatomic_ham_interaction_treatment_consumed,
                        route_configured_diatomic_ham_interaction_treatment_status =
                            diatomic_ham_interaction_treatment_status,
                        shellization_summary_available,
                        shellization_source =
                            :route_configured_bond_aligned_diatomic_source,
                        route_configured_shellization_consumed = true,
                        route_configured_legacy_diatomic_source_consumed = true,
                        materialized_shellization_stage =
                            shellization_summary.shellization_stage,
                        seed_materialization_status =
                            :not_seed_route_configured_diatomic_shellization,
                        export_status = :basis_and_ham,
                        basis_export_status =
                            :supported_route_configured_diatomic_basis_only_fixed_block,
                        ham_preflight_status =
                            :available_route_configured_diatomic_ham_adapter,
                        ham_operator_payload_status =
                            :available_route_configured_diatomic_operator_payload,
                        ham_interaction_status =
                            :available_route_configured_diatomic_density_density_interaction_matrix,
                        ham_export_status =
                            :available_route_configured_diatomic_ham_bundle_payload,
                        ham_export_blocker = nothing,
                        private_development_only = true,
                    ),
                )
                ham_artifact_written = true
            end
            return (;
                object_kind = :cartesian_nesting_route_driver_materialization,
                route_family,
                private_development_only = true,
                materialize_route_requested = true,
                save_basis_artifact_requested = save_basis_artifact,
                save_ham_artifact_requested = save_ham_artifact,
                route_configured_diatomic_ham_interaction_treatment_requested =
                    route_configured_diatomic_ham_interaction_treatment,
                route_configured_diatomic_ham_interaction_treatment_consumed =
                    diatomic_ham_interaction_treatment_consumed,
                route_configured_diatomic_ham_interaction_treatment_status =
                    diatomic_ham_interaction_treatment_status,
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
                route_configured_diatomic_atom_growth_materializer_probe,
                route_configured_diatomic_materializer_contract...,
                route_configured_diatomic_atom_growth_materializer_contract...,
                route_configured_materializer_contract...,
                route_configured_diatomic_basis_adapter_summary =
                    diatomic_basis_adapter_summary,
                route_configured_diatomic_ham_adapter_summary =
                    diatomic_ham_adapter_summary,
                shellization_summary,
                shellization_summary_available,
                shellization_source = :route_configured_bond_aligned_diatomic_source,
                route_configured_shellization_consumed = true,
                route_configured_legacy_diatomic_source_consumed = true,
                materialized_shellization_stage =
                    shellization_summary.shellization_stage,
                seed_materialization_status =
                    :not_seed_route_configured_diatomic_shellization,
                retained_dimension = diatomic_materialization.retained_dimension,
                final_integral_weights_status =
                    diatomic_basis_adapter.final_integral_weights_status,
                one_body_operator_status =
                    :pending_route_configured_diatomic_operator_inventory,
                basis_bundle_export_status =
                    diatomic_basis_adapter_available ?
                    :supported_route_configured_diatomic_basis_only_fixed_block :
                    :pending_route_configured_diatomic_basis_export,
                basis_artifact_status,
                basis_artifact_written,
                basisfile,
                basis_artifact_path = basis_artifact_written ? basisfile : nothing,
                basis_export_blocker =
                    diatomic_basis_adapter_available ?
                    nothing :
                    :pending_route_configured_diatomic_basis_adapter_contract,
                ham_preflight_status =
                    diatomic_ham_adapter_available ?
                    :available_route_configured_diatomic_ham_adapter :
                    save_ham_artifact ?
                    :blocked_route_configured_diatomic_ham_adapter :
                    :not_requested,
                ham_missing_builder =
                    diatomic_ham_adapter_available ?
                    nothing :
                    save_ham_artifact ?
                    something(
                        diatomic_ham_adapter_blocker,
                        :pending_route_configured_diatomic_ham_adapter,
                    ) :
                    nothing,
                ham_operator_payload_status =
                    isnothing(diatomic_ham_adapter) ?
                    :not_requested :
                    diatomic_ham_adapter.operator_payload_status,
                ham_interaction_status =
                    isnothing(diatomic_ham_adapter) ?
                    :not_requested :
                    diatomic_ham_adapter.interaction_status,
                ham_bundle_export_status =
                    diatomic_ham_bundle_export_status,
                ham_artifact_status,
                ham_artifact_written,
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
                    route_configured_diatomic_atom_growth_materializer_contract...,
                    route_configured_materializer_contract...,
                    shellization_summary_available,
                    shellization_source,
                    route_configured_shellization_consumed,
                    route_configured_legacy_diatomic_source_consumed = false,
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
                    route_configured_diatomic_atom_growth_materializer_contract...,
                    route_configured_materializer_contract...,
                    shellization_summary_available,
                    shellization_source,
                    route_configured_shellization_consumed,
                    route_configured_legacy_diatomic_source_consumed = false,
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
            route_configured_diatomic_ham_interaction_treatment_requested = nothing,
            route_configured_diatomic_ham_interaction_treatment_consumed = nothing,
            route_configured_diatomic_ham_interaction_treatment_status =
                :not_applicable,
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
            route_configured_diatomic_atom_growth_materializer_probe,
            route_configured_diatomic_materializer_contract...,
            route_configured_diatomic_atom_growth_materializer_contract...,
            route_configured_materializer_contract...,
            shellization_summary,
            shellization_summary_available,
            shellization_source,
            route_configured_shellization_consumed,
            route_configured_legacy_diatomic_source_consumed = false,
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
        route_configured_diatomic_ham_interaction_treatment_requested = nothing,
        route_configured_diatomic_ham_interaction_treatment_consumed = nothing,
        route_configured_diatomic_ham_interaction_treatment_status = :not_applicable,
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
        route_configured_diatomic_atom_growth_materializer_probe,
        route_configured_diatomic_materializer_contract...,
        route_configured_diatomic_atom_growth_materializer_contract...,
        route_configured_materializer_contract...,
        shellization_summary = nothing,
        shellization_summary_available = false,
        shellization_source = :pending_source_box_route_shellization,
        route_configured_shellization_consumed = false,
        route_configured_legacy_diatomic_source_consumed = false,
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
        bond_axis =
            parent.system_classification == :bond_aligned_diatomic ?
            parent.bond_axis :
            :auto
        plan = _cartesian_terminal_shellification_geometry(
            parent_axes,
            nuclear_positions;
            core_side = parent.standard_setup.core_cube_side,
            q = parent.standard_setup.q,
            bond_axis,
        )
        scaffold = _cartesian_terminal_shellification_geometry_scaffold(
            plan;
            route_family = :white_lindsey_low_order,
        )
        status = _pqs_source_box_route_driver_terminal_stage_status(scaffold)
        materialization_available = false
        blocker =
            scaffold.materialization_status ==
            :deferred_pending_distorted_product_box_lowering ?
            :distorted_product_box_lowering_not_implemented :
            nothing
        return (;
            object_kind = :cartesian_shell_stage_terminal_shellification_payload,
            status,
            plan,
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
        recipe.route_family == :white_lindsey_low_order &&
        parent.system_classification in (:one_center, :bond_aligned_diatomic) &&
        policy.low_order_shellization_policy_status ==
        :available_low_order_shellization_policy &&
        policy_resolved == :terminal_cartesian_shellification_geometry
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

function _pqs_source_box_route_driver_unit_stage_low_order_summary(shells)
    low_order_shellization =
        hasproperty(shells, :low_order_shellization) ?
        shells.low_order_shellization :
        nothing
    if isnothing(low_order_shellization)
        return (;
            object_kind = :cartesian_unit_stage_low_order_summary,
            status = :not_available_missing_shell_stage_summary,
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
    terminal_region_lowering_contract_inventory =
        terminal_region_unit_inventory_available ?
        _cartesian_terminal_region_lowering_contract_inventory(
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
            :deferred_terminal_shellification_lowering_contract_inventory
        ) :
        :not_selected
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
        terminal_shellification_lowering_contract_count =
            terminal_region_lowering_contract_inventory_available ?
            terminal_region_lowering_contract_inventory.lowering_contract_count :
            0,
        terminal_shellification_lowering_contract_kinds =
            terminal_region_lowering_contract_inventory_available ?
            terminal_region_lowering_contract_inventory.lowering_contract_kinds :
            (),
        terminal_shellification_lowering_contract_kind_counts =
            terminal_region_lowering_contract_inventory_available ?
            terminal_region_lowering_contract_inventory.lowering_contract_kind_counts :
            (
                direct_core_identity_cpb_count = 0,
                direct_slab_identity_cpb_count = 0,
                direct_boundary_slab_identity_cpb_count = 0,
                white_lindsey_boundary_strata_count = 0,
                pqs_filled_source_cpb_count = 0,
                distorted_product_box_comx_count = 0,
            ),
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
        low_order_unit_route_kind = low_order_units.unit_route_kind,
        atom_growth_unit_summary_available =
            low_order_units.atom_growth_unit_summary_available,
        atom_growth_units_selected = low_order_units.atom_growth_units_selected,
        atom_growth_unit_inventory_available =
            low_order_units.atom_growth_unit_inventory_available,
        terminal_shellification_units_selected =
            low_order_units.terminal_shellification_units_selected,
        terminal_shellification_unit_summary_available =
            low_order_units.terminal_shellification_unit_summary_available,
        terminal_shellification_scaffold_available =
            low_order_units.terminal_shellification_scaffold_available,
        terminal_shellification_scaffold =
            low_order_units.terminal_shellification_scaffold,
        terminal_shellification_region_count =
            low_order_units.terminal_shellification_region_count,
        terminal_shellification_unit_inventory_available =
            low_order_units.terminal_shellification_unit_inventory_available,
        terminal_shellification_unit_inventory_status =
            low_order_units.terminal_shellification_unit_inventory_status,
        terminal_shellification_unit_inventory =
            low_order_units.terminal_shellification_unit_inventory,
        terminal_shellification_unit_count =
            low_order_units.terminal_shellification_unit_count,
        terminal_shellification_unit_keys =
            low_order_units.terminal_shellification_unit_keys,
        terminal_shellification_unit_roles =
            low_order_units.terminal_shellification_unit_roles,
        terminal_shellification_unit_kinds =
            low_order_units.terminal_shellification_unit_kinds,
        terminal_shellification_unit_support_counts =
            low_order_units.terminal_shellification_unit_support_counts,
        terminal_shellification_lowering_contract_inventory_available =
            low_order_units.terminal_shellification_lowering_contract_inventory_available,
        terminal_shellification_lowering_contract_inventory_status =
            low_order_units.terminal_shellification_lowering_contract_inventory_status,
        terminal_shellification_lowering_contract_inventory =
            low_order_units.terminal_shellification_lowering_contract_inventory,
        terminal_shellification_lowering_contract_count =
            low_order_units.terminal_shellification_lowering_contract_count,
        terminal_shellification_lowering_contract_kinds =
            low_order_units.terminal_shellification_lowering_contract_kinds,
        terminal_shellification_lowering_contract_kind_counts =
            low_order_units.terminal_shellification_lowering_contract_kind_counts,
        terminal_shellification_contract_counts_by_unit =
            low_order_units.terminal_shellification_contract_counts_by_unit,
        terminal_shellification_lw_complete_shell_cpb_count =
            low_order_units.terminal_shellification_lw_complete_shell_cpb_count,
        terminal_shellification_lw_complete_shell_cpb_family_counts =
            low_order_units.terminal_shellification_lw_complete_shell_cpb_family_counts,
        terminal_shellification_final_retained_unit_inventory_available =
            low_order_units.terminal_shellification_final_retained_unit_inventory_available,
        terminal_shellification_pair_inventory_available =
            low_order_units.terminal_shellification_pair_inventory_available,
        terminal_shellification_central_gap_region_count =
            low_order_units.terminal_shellification_central_gap_region_count,
        terminal_shellification_central_midpoint_slab_count =
            low_order_units.terminal_shellification_central_midpoint_slab_count,
        terminal_shellification_central_distorted_product_box_count =
            low_order_units.terminal_shellification_central_distorted_product_box_count,
        terminal_shellification_central_distorted_product_box_metadata =
            low_order_units.terminal_shellification_central_distorted_product_box_metadata,
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
        return (;
            object_kind = :cartesian_transform_stage_low_order_summary,
            status = :not_available_missing_unit_stage_summary,
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
            terminal_shellification_contract_counts_by_unit = (),
            terminal_shellification_lw_complete_shell_cpb_count = 0,
            terminal_shellification_lw_complete_shell_cpb_family_counts =
                (facet_cpb = 0, edge_cpb = 0, corner_cpb = 0),
            terminal_shellification_final_retained_unit_inventory_available = false,
            terminal_shellification_pair_inventory_available = false,
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
        terminal_shellification_scaffold_available,
        terminal_shellification_scaffold =
            terminal_shellification_transforms_selected ?
            low_order_units.terminal_shellification_scaffold :
            nothing,
        terminal_shellification_region_count =
            terminal_shellification_transforms_selected ?
            low_order_units.terminal_shellification_region_count :
            0,
        terminal_shellification_unit_inventory_available =
            terminal_shellification_transforms_selected &&
            low_order_units.terminal_shellification_unit_inventory_available,
        terminal_shellification_unit_inventory =
            terminal_shellification_transforms_selected ?
            low_order_units.terminal_shellification_unit_inventory :
            nothing,
        terminal_shellification_unit_count =
            terminal_shellification_transforms_selected ?
            low_order_units.terminal_shellification_unit_count :
            0,
        terminal_shellification_unit_keys =
            terminal_shellification_transforms_selected ?
            low_order_units.terminal_shellification_unit_keys :
            (),
        terminal_shellification_unit_roles =
            terminal_shellification_transforms_selected ?
            low_order_units.terminal_shellification_unit_roles :
            (),
        terminal_shellification_unit_kinds =
            terminal_shellification_transforms_selected ?
            low_order_units.terminal_shellification_unit_kinds :
            (),
        terminal_shellification_unit_support_counts =
            terminal_shellification_transforms_selected ?
            low_order_units.terminal_shellification_unit_support_counts :
            (),
        terminal_shellification_lowering_contract_inventory_available =
            terminal_shellification_transforms_selected &&
            low_order_units.terminal_shellification_lowering_contract_inventory_available,
        terminal_shellification_lowering_contract_inventory_status =
            terminal_shellification_transforms_selected ?
            low_order_units.terminal_shellification_lowering_contract_inventory_status :
            :not_selected,
        terminal_shellification_lowering_contract_inventory =
            terminal_shellification_transforms_selected ?
            low_order_units.terminal_shellification_lowering_contract_inventory :
            nothing,
        terminal_shellification_lowering_contract_count =
            terminal_shellification_transforms_selected ?
            low_order_units.terminal_shellification_lowering_contract_count :
            0,
        terminal_shellification_lowering_contract_kinds =
            terminal_shellification_transforms_selected ?
            low_order_units.terminal_shellification_lowering_contract_kinds :
            (),
        terminal_shellification_lowering_contract_kind_counts =
            terminal_shellification_transforms_selected ?
            low_order_units.terminal_shellification_lowering_contract_kind_counts :
            (
                direct_core_identity_cpb_count = 0,
                direct_slab_identity_cpb_count = 0,
                direct_boundary_slab_identity_cpb_count = 0,
                white_lindsey_boundary_strata_count = 0,
                pqs_filled_source_cpb_count = 0,
                distorted_product_box_comx_count = 0,
            ),
        terminal_shellification_contract_counts_by_unit =
            terminal_shellification_transforms_selected ?
            low_order_units.terminal_shellification_contract_counts_by_unit :
            (),
        terminal_shellification_lw_complete_shell_cpb_count =
            terminal_shellification_transforms_selected ?
            low_order_units.terminal_shellification_lw_complete_shell_cpb_count :
            0,
        terminal_shellification_lw_complete_shell_cpb_family_counts =
            terminal_shellification_transforms_selected ?
            low_order_units.terminal_shellification_lw_complete_shell_cpb_family_counts :
            (facet_cpb = 0, edge_cpb = 0, corner_cpb = 0),
        terminal_shellification_final_retained_unit_inventory_available =
            terminal_shellification_transforms_selected &&
            low_order_units.terminal_shellification_final_retained_unit_inventory_available,
        terminal_shellification_pair_inventory_available =
            terminal_shellification_transforms_selected &&
            low_order_units.terminal_shellification_pair_inventory_available,
        terminal_shellification_transform_contracts_available = false,
        terminal_shellification_transform_materialization_status,
        terminal_shellification_central_gap_region_count =
            terminal_shellification_transforms_selected ?
            low_order_units.terminal_shellification_central_gap_region_count :
            0,
        terminal_shellification_central_midpoint_slab_count =
            terminal_shellification_transforms_selected ?
            low_order_units.terminal_shellification_central_midpoint_slab_count :
            0,
        terminal_shellification_central_distorted_product_box_count =
            terminal_shellification_transforms_selected ?
            low_order_units.terminal_shellification_central_distorted_product_box_count :
            0,
        terminal_shellification_central_distorted_product_box_metadata =
            terminal_shellification_transforms_selected ?
            low_order_units.terminal_shellification_central_distorted_product_box_metadata :
            (),
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
        low_order_transform_route_kind =
            low_order_transforms.transform_route_kind,
        atom_growth_transforms_selected =
            low_order_transforms.atom_growth_transforms_selected,
        terminal_shellification_transforms_selected =
            low_order_transforms.terminal_shellification_transforms_selected,
        terminal_shellification_transform_summary_available =
            low_order_transforms.terminal_shellification_transform_summary_available,
        terminal_shellification_scaffold_available =
            low_order_transforms.terminal_shellification_scaffold_available,
        terminal_shellification_scaffold =
            low_order_transforms.terminal_shellification_scaffold,
        terminal_shellification_region_count =
            low_order_transforms.terminal_shellification_region_count,
        terminal_shellification_unit_inventory_available =
            low_order_transforms.terminal_shellification_unit_inventory_available,
        terminal_shellification_unit_inventory =
            low_order_transforms.terminal_shellification_unit_inventory,
        terminal_shellification_unit_count =
            low_order_transforms.terminal_shellification_unit_count,
        terminal_shellification_unit_keys =
            low_order_transforms.terminal_shellification_unit_keys,
        terminal_shellification_unit_roles =
            low_order_transforms.terminal_shellification_unit_roles,
        terminal_shellification_unit_kinds =
            low_order_transforms.terminal_shellification_unit_kinds,
        terminal_shellification_unit_support_counts =
            low_order_transforms.terminal_shellification_unit_support_counts,
        terminal_shellification_lowering_contract_inventory_available =
            low_order_transforms.terminal_shellification_lowering_contract_inventory_available,
        terminal_shellification_lowering_contract_inventory_status =
            low_order_transforms.terminal_shellification_lowering_contract_inventory_status,
        terminal_shellification_lowering_contract_inventory =
            low_order_transforms.terminal_shellification_lowering_contract_inventory,
        terminal_shellification_lowering_contract_count =
            low_order_transforms.terminal_shellification_lowering_contract_count,
        terminal_shellification_lowering_contract_kinds =
            low_order_transforms.terminal_shellification_lowering_contract_kinds,
        terminal_shellification_lowering_contract_kind_counts =
            low_order_transforms.terminal_shellification_lowering_contract_kind_counts,
        terminal_shellification_contract_counts_by_unit =
            low_order_transforms.terminal_shellification_contract_counts_by_unit,
        terminal_shellification_lw_complete_shell_cpb_count =
            low_order_transforms.terminal_shellification_lw_complete_shell_cpb_count,
        terminal_shellification_lw_complete_shell_cpb_family_counts =
            low_order_transforms.terminal_shellification_lw_complete_shell_cpb_family_counts,
        terminal_shellification_final_retained_unit_inventory_available =
            low_order_transforms.terminal_shellification_final_retained_unit_inventory_available,
        terminal_shellification_pair_inventory_available =
            low_order_transforms.terminal_shellification_pair_inventory_available,
        terminal_shellification_transform_contracts_available =
            low_order_transforms.terminal_shellification_transform_contracts_available,
        terminal_shellification_transform_materialization_status =
            low_order_transforms.terminal_shellification_transform_materialization_status,
        terminal_shellification_central_gap_region_count =
            low_order_transforms.terminal_shellification_central_gap_region_count,
        terminal_shellification_central_midpoint_slab_count =
            low_order_transforms.terminal_shellification_central_midpoint_slab_count,
        terminal_shellification_central_distorted_product_box_count =
            low_order_transforms.terminal_shellification_central_distorted_product_box_count,
        terminal_shellification_central_distorted_product_box_metadata =
            low_order_transforms.terminal_shellification_central_distorted_product_box_metadata,
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
        return (;
            object_kind = :cartesian_pair_stage_low_order_summary,
            status = :not_available_missing_transform_stage_summary,
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
        terminal_shellification_scaffold_available,
        terminal_shellification_scaffold =
            terminal_shellification_pairs_selected ?
            low_order_transforms.terminal_shellification_scaffold :
            nothing,
        terminal_shellification_region_count =
            terminal_shellification_pairs_selected ?
            low_order_transforms.terminal_shellification_region_count :
            0,
        terminal_shellification_unit_inventory_available =
            terminal_shellification_pairs_selected &&
            low_order_transforms.terminal_shellification_unit_inventory_available,
        terminal_shellification_unit_inventory =
            terminal_shellification_pairs_selected ?
            low_order_transforms.terminal_shellification_unit_inventory :
            nothing,
        terminal_shellification_unit_count =
            terminal_shellification_pairs_selected ?
            low_order_transforms.terminal_shellification_unit_count :
            0,
        terminal_shellification_unit_keys =
            terminal_shellification_pairs_selected ?
            low_order_transforms.terminal_shellification_unit_keys :
            (),
        terminal_shellification_unit_roles =
            terminal_shellification_pairs_selected ?
            low_order_transforms.terminal_shellification_unit_roles :
            (),
        terminal_shellification_unit_kinds =
            terminal_shellification_pairs_selected ?
            low_order_transforms.terminal_shellification_unit_kinds :
            (),
        terminal_shellification_unit_support_counts =
            terminal_shellification_pairs_selected ?
            low_order_transforms.terminal_shellification_unit_support_counts :
            (),
        terminal_shellification_final_retained_unit_inventory_available =
            terminal_shellification_pairs_selected &&
            low_order_transforms.terminal_shellification_final_retained_unit_inventory_available,
        terminal_shellification_transform_contracts_available =
            terminal_shellification_pairs_selected &&
            low_order_transforms.terminal_shellification_transform_contracts_available,
        terminal_shellification_pair_inventory_available = false,
        terminal_shellification_pair_inventory_status,
        terminal_shellification_pair_materialization_status,
        terminal_shellification_central_gap_region_count =
            terminal_shellification_pairs_selected ?
            low_order_transforms.terminal_shellification_central_gap_region_count :
            0,
        terminal_shellification_central_midpoint_slab_count =
            terminal_shellification_pairs_selected ?
            low_order_transforms.terminal_shellification_central_midpoint_slab_count :
            0,
        terminal_shellification_central_distorted_product_box_count =
            terminal_shellification_pairs_selected ?
            low_order_transforms.terminal_shellification_central_distorted_product_box_count :
            0,
        terminal_shellification_central_distorted_product_box_metadata =
            terminal_shellification_pairs_selected ?
            low_order_transforms.terminal_shellification_central_distorted_product_box_metadata :
            (),
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
        low_order_pair_route_kind = low_order_pairs.pair_route_kind,
        atom_growth_pairs_selected = low_order_pairs.atom_growth_pairs_selected,
        terminal_shellification_pairs_selected =
            low_order_pairs.terminal_shellification_pairs_selected,
        terminal_shellification_pair_summary_available =
            low_order_pairs.terminal_shellification_pair_summary_available,
        terminal_shellification_scaffold_available =
            low_order_pairs.terminal_shellification_scaffold_available,
        terminal_shellification_scaffold =
            low_order_pairs.terminal_shellification_scaffold,
        terminal_shellification_region_count =
            low_order_pairs.terminal_shellification_region_count,
        terminal_shellification_unit_inventory_available =
            low_order_pairs.terminal_shellification_unit_inventory_available,
        terminal_shellification_unit_inventory =
            low_order_pairs.terminal_shellification_unit_inventory,
        terminal_shellification_unit_count =
            low_order_pairs.terminal_shellification_unit_count,
        terminal_shellification_unit_keys =
            low_order_pairs.terminal_shellification_unit_keys,
        terminal_shellification_unit_roles =
            low_order_pairs.terminal_shellification_unit_roles,
        terminal_shellification_unit_kinds =
            low_order_pairs.terminal_shellification_unit_kinds,
        terminal_shellification_unit_support_counts =
            low_order_pairs.terminal_shellification_unit_support_counts,
        terminal_shellification_final_retained_unit_inventory_available =
            low_order_pairs.terminal_shellification_final_retained_unit_inventory_available,
        terminal_shellification_transform_contracts_available =
            low_order_pairs.terminal_shellification_transform_contracts_available,
        terminal_shellification_pair_inventory_available =
            low_order_pairs.terminal_shellification_pair_inventory_available,
        terminal_shellification_pair_inventory_status =
            low_order_pairs.terminal_shellification_pair_inventory_status,
        terminal_shellification_pair_materialization_status =
            low_order_pairs.terminal_shellification_pair_materialization_status,
        terminal_shellification_central_gap_region_count =
            low_order_pairs.terminal_shellification_central_gap_region_count,
        terminal_shellification_central_midpoint_slab_count =
            low_order_pairs.terminal_shellification_central_midpoint_slab_count,
        terminal_shellification_central_distorted_product_box_count =
            low_order_pairs.terminal_shellification_central_distorted_product_box_count,
        terminal_shellification_central_distorted_product_box_metadata =
            low_order_pairs.terminal_shellification_central_distorted_product_box_metadata,
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
        route_core_pair_inventory = low_order_pairs.route_core_pair_inventory,
        route_core_pair_count = low_order_pairs.route_core_pair_count,
        route_core_pair_keys = low_order_pairs.route_core_pair_keys,
        route_core_pair_order_matches_staged =
            low_order_pairs.route_core_pair_order_matches_staged,
        route_core_pair_order_comparison_source =
            low_order_pairs.route_core_pair_order_comparison_source,
        route_core_pair_family_counts =
            low_order_pairs.route_core_pair_family_counts,
        route_core_pair_family_count_source =
            low_order_pairs.route_core_pair_family_count_source,
        route_core_summary_status =
            low_order_pairs.route_core_summary_status,
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
        return (;
            object_kind = :cartesian_assembly_stage_low_order_summary,
            status = :not_available_missing_pair_stage_summary,
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

    return (;
        object_kind = :cartesian_assembly_stage_low_order_summary,
        status =
            terminal_shellification_assembly_selected &&
            terminal_shellification_scaffold_available ?
            :deferred_terminal_shellification_assembly_materialization :
            low_order_pairs.status == :available_pair_stage_low_order_summary ?
            :available_assembly_stage_low_order_summary :
            low_order_pairs.status,
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
        terminal_shellification_scaffold_available,
        terminal_shellification_scaffold =
            terminal_shellification_assembly_selected ?
            low_order_pairs.terminal_shellification_scaffold :
            nothing,
        terminal_shellification_region_count =
            terminal_shellification_assembly_selected ?
            low_order_pairs.terminal_shellification_region_count :
            0,
        terminal_shellification_unit_inventory_available =
            terminal_shellification_assembly_selected &&
            low_order_pairs.terminal_shellification_unit_inventory_available,
        terminal_shellification_unit_inventory =
            terminal_shellification_assembly_selected ?
            low_order_pairs.terminal_shellification_unit_inventory :
            nothing,
        terminal_shellification_unit_count =
            terminal_shellification_assembly_selected ?
            low_order_pairs.terminal_shellification_unit_count :
            0,
        terminal_shellification_unit_keys =
            terminal_shellification_assembly_selected ?
            low_order_pairs.terminal_shellification_unit_keys :
            (),
        terminal_shellification_unit_roles =
            terminal_shellification_assembly_selected ?
            low_order_pairs.terminal_shellification_unit_roles :
            (),
        terminal_shellification_unit_kinds =
            terminal_shellification_assembly_selected ?
            low_order_pairs.terminal_shellification_unit_kinds :
            (),
        terminal_shellification_unit_support_counts =
            terminal_shellification_assembly_selected ?
            low_order_pairs.terminal_shellification_unit_support_counts :
            (),
        terminal_shellification_final_retained_unit_inventory_available =
            terminal_shellification_assembly_selected &&
            low_order_pairs.terminal_shellification_final_retained_unit_inventory_available,
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
        terminal_shellification_central_gap_region_count =
            terminal_shellification_assembly_selected ?
            low_order_pairs.terminal_shellification_central_gap_region_count :
            0,
        terminal_shellification_central_midpoint_slab_count =
            terminal_shellification_assembly_selected ?
            low_order_pairs.terminal_shellification_central_midpoint_slab_count :
            0,
        terminal_shellification_central_distorted_product_box_count =
            terminal_shellification_assembly_selected ?
            low_order_pairs.terminal_shellification_central_distorted_product_box_count :
            0,
        terminal_shellification_central_distorted_product_box_metadata =
            terminal_shellification_assembly_selected ?
            low_order_pairs.terminal_shellification_central_distorted_product_box_metadata :
            (),
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

function cartesian_assembly(parent, shells, units, transforms, pairs, recipe)
    route_skeleton = shells.route_skeleton
    route_facts = _pqs_source_box_route_driver_route_facts(route_skeleton)
    contract = _pqs_source_box_route_driver_contract_metadata(recipe)
    low_order_assembly =
        _pqs_source_box_route_driver_assembly_stage_low_order_summary(pairs)

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
        low_order_assembly_source = low_order_assembly.assembly_source,
        low_order_assembly_route_kind =
            low_order_assembly.assembly_route_kind,
        atom_growth_assembly_selected =
            low_order_assembly.atom_growth_assembly_selected,
        terminal_shellification_assembly_selected =
            low_order_assembly.terminal_shellification_assembly_selected,
        terminal_shellification_assembly_summary_available =
            low_order_assembly.terminal_shellification_assembly_summary_available,
        terminal_shellification_scaffold_available =
            low_order_assembly.terminal_shellification_scaffold_available,
        terminal_shellification_scaffold =
            low_order_assembly.terminal_shellification_scaffold,
        terminal_shellification_region_count =
            low_order_assembly.terminal_shellification_region_count,
        terminal_shellification_unit_inventory_available =
            low_order_assembly.terminal_shellification_unit_inventory_available,
        terminal_shellification_unit_inventory =
            low_order_assembly.terminal_shellification_unit_inventory,
        terminal_shellification_unit_count =
            low_order_assembly.terminal_shellification_unit_count,
        terminal_shellification_unit_keys =
            low_order_assembly.terminal_shellification_unit_keys,
        terminal_shellification_unit_roles =
            low_order_assembly.terminal_shellification_unit_roles,
        terminal_shellification_unit_kinds =
            low_order_assembly.terminal_shellification_unit_kinds,
        terminal_shellification_unit_support_counts =
            low_order_assembly.terminal_shellification_unit_support_counts,
        terminal_shellification_final_retained_unit_inventory_available =
            low_order_assembly.terminal_shellification_final_retained_unit_inventory_available,
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
        terminal_shellification_central_gap_region_count =
            low_order_assembly.terminal_shellification_central_gap_region_count,
        terminal_shellification_central_midpoint_slab_count =
            low_order_assembly.terminal_shellification_central_midpoint_slab_count,
        terminal_shellification_central_distorted_product_box_count =
            low_order_assembly.terminal_shellification_central_distorted_product_box_count,
        terminal_shellification_central_distorted_product_box_metadata =
            low_order_assembly.terminal_shellification_central_distorted_product_box_metadata,
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
        low_order_route_core_pair_operator_ready =
            low_order_assembly.route_core_pair_operator_ready,
        low_order_route_core_pair_operator_readiness_status =
            low_order_assembly.route_core_pair_operator_readiness_status,
        low_order_route_core_pair_operator_blocker =
            low_order_assembly.route_core_pair_operator_blocker,
        low_order_route_core_pair_operator_readiness_requirements =
            low_order_assembly.route_core_pair_operator_readiness_requirements,
        low_order_route_core_pair_operator_preflight_available =
            low_order_assembly.route_core_pair_operator_preflight_available,
        low_order_route_core_pair_operator_preflight_status =
            low_order_assembly.route_core_pair_operator_preflight_status,
        low_order_route_core_pair_operator_preflight =
            low_order_assembly.route_core_pair_operator_preflight,
        low_order_route_core_pair_operator_preflight_blocker =
            low_order_assembly.route_core_pair_operator_preflight_blocker,
        low_order_route_core_pair_operator_plan_available =
            low_order_assembly.route_core_pair_operator_plan_available,
        low_order_route_core_pair_operator_plan_status =
            low_order_assembly.route_core_pair_operator_plan_status,
        low_order_route_core_pair_operator_plan =
            low_order_assembly.route_core_pair_operator_plan,
        low_order_route_core_pair_operator_plan_blocker =
            low_order_assembly.route_core_pair_operator_plan_blocker,
        low_order_route_core_typed_pair_operator_plan_inventory_available =
            low_order_assembly.route_core_typed_pair_operator_plan_inventory_available,
        low_order_route_core_typed_pair_operator_plan_inventory_status =
            low_order_assembly.route_core_typed_pair_operator_plan_inventory_status,
        low_order_route_core_typed_pair_operator_plan_blocker =
            low_order_assembly.route_core_typed_pair_operator_plan_blocker,
        low_order_route_core_typed_pair_operator_plan_count =
            low_order_assembly.route_core_typed_pair_operator_plan_count,
        low_order_route_core_typed_pair_operator_plan_blocked_count =
            low_order_assembly.route_core_typed_pair_operator_plan_blocked_count,
        low_order_route_core_typed_pair_operator_plan_materialized =
            low_order_assembly.route_core_typed_pair_operator_plan_materialized,
        low_order_route_core_typed_pair_operator_source_path_counts =
            low_order_assembly.route_core_typed_pair_operator_source_path_counts,
        low_order_route_core_typed_pair_operator_final_block_path_counts =
            low_order_assembly.route_core_typed_pair_operator_final_block_path_counts,
        low_order_route_core_typed_pair_operator_materialization_status_counts =
            low_order_assembly.route_core_typed_pair_operator_materialization_status_counts,
        low_order_route_core_typed_pair_operator_blocker_counts =
            low_order_assembly.route_core_typed_pair_operator_blocker_counts,
        low_order_route_core_typed_pair_operator_plan_family_counts =
            low_order_assembly.route_core_typed_pair_operator_plan_family_counts,
        low_order_route_core_typed_pair_operator_materialization_ready =
            low_order_assembly.route_core_typed_pair_operator_materialization_ready,
        low_order_route_core_typed_pair_operator_materialization_readiness_status =
            low_order_assembly.route_core_typed_pair_operator_materialization_readiness_status,
        low_order_route_core_typed_pair_operator_materialization_readiness_blocker =
            low_order_assembly.route_core_typed_pair_operator_materialization_readiness_blocker,
        low_order_route_core_typed_pair_operator_materialization_readiness_requirements =
            low_order_assembly.route_core_typed_pair_operator_materialization_readiness_requirements,
        low_order_route_core_typed_pair_operator_materialization_readiness_plan_count =
            low_order_assembly.route_core_typed_pair_operator_materialization_readiness_plan_count,
        low_order_route_core_typed_pair_operator_materialization_readiness_blocked_count =
            low_order_assembly.route_core_typed_pair_operator_materialization_readiness_blocked_count,
        low_order_route_core_typed_pair_operator_materialization_readiness_materialized_count =
            low_order_assembly.route_core_typed_pair_operator_materialization_readiness_materialized_count,
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
        terminal_shellification_scaffold_available =
            low_order_assembly.terminal_shellification_scaffold_available,
        terminal_shellification_scaffold =
            low_order_assembly.terminal_shellification_scaffold,
        terminal_shellification_region_count =
            low_order_assembly.terminal_shellification_region_count,
        terminal_shellification_unit_inventory_available =
            low_order_assembly.terminal_shellification_unit_inventory_available,
        terminal_shellification_unit_inventory =
            low_order_assembly.terminal_shellification_unit_inventory,
        terminal_shellification_unit_count =
            low_order_assembly.terminal_shellification_unit_count,
        terminal_shellification_unit_keys =
            low_order_assembly.terminal_shellification_unit_keys,
        terminal_shellification_unit_roles =
            low_order_assembly.terminal_shellification_unit_roles,
        terminal_shellification_unit_kinds =
            low_order_assembly.terminal_shellification_unit_kinds,
        terminal_shellification_unit_support_counts =
            low_order_assembly.terminal_shellification_unit_support_counts,
        terminal_shellification_final_retained_unit_inventory_available =
            low_order_assembly.terminal_shellification_final_retained_unit_inventory_available,
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
        terminal_shellification_central_gap_region_count =
            low_order_assembly.terminal_shellification_central_gap_region_count,
        terminal_shellification_central_midpoint_slab_count =
            low_order_assembly.terminal_shellification_central_midpoint_slab_count,
        terminal_shellification_central_distorted_product_box_count =
            low_order_assembly.terminal_shellification_central_distorted_product_box_count,
        terminal_shellification_central_distorted_product_box_metadata =
            low_order_assembly.terminal_shellification_central_distorted_product_box_metadata,
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
        route_core_typed_pair_operator_plan_inventory_available =
            low_order_assembly.route_core_typed_pair_operator_plan_inventory_available,
        route_core_typed_pair_operator_plan_inventory_status =
            low_order_assembly.route_core_typed_pair_operator_plan_inventory_status,
        route_core_typed_pair_operator_plan_blocker =
            low_order_assembly.route_core_typed_pair_operator_plan_blocker,
        route_core_typed_pair_operator_plan_count =
            low_order_assembly.route_core_typed_pair_operator_plan_count,
        route_core_typed_pair_operator_plan_blocked_count =
            low_order_assembly.route_core_typed_pair_operator_plan_blocked_count,
        route_core_typed_pair_operator_plan_materialized =
            low_order_assembly.route_core_typed_pair_operator_plan_materialized,
        route_core_typed_pair_operator_source_path_counts =
            low_order_assembly.route_core_typed_pair_operator_source_path_counts,
        route_core_typed_pair_operator_final_block_path_counts =
            low_order_assembly.route_core_typed_pair_operator_final_block_path_counts,
        route_core_typed_pair_operator_materialization_status_counts =
            low_order_assembly.route_core_typed_pair_operator_materialization_status_counts,
        route_core_typed_pair_operator_blocker_counts =
            low_order_assembly.route_core_typed_pair_operator_blocker_counts,
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

    return _pqs_source_box_route_driver_report(
        standard_setup, parent, parent_axis, route_axis_counts, raw_box,
        system_metadata, recipe_metadata, parent_contract, parent_description,
        route_skeleton, route_facts, contract, diagnostics,
        low_order_route_summary)
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

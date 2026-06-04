function _cartesian_shellization_layer_kind(layer::_AbstractCartesianNestedShellLayer3D)
    return Symbol(nameof(typeof(layer)))
end

function _cartesian_shellization_layer_column_ranges(
    sequence::_CartesianNestedShellSequence3D,
)
    return Tuple(sequence.layer_column_ranges)
end

function _cartesian_shellization_layer_kinds(sequence::_CartesianNestedShellSequence3D)
    return Tuple(_cartesian_shellization_layer_kind(layer) for layer in sequence.shell_layers)
end

function _cartesian_shellization_common_diagnostics(;
    source_kind::Symbol,
    shellization_role::Symbol,
)
    return (
        source_kind = source_kind,
        shellization_role = shellization_role,
        private_development_only = true,
        route_neutral_spatial_planning = true,
        lowering_applied_by_summary = false,
        white_lindsey_lowering_adopted_by_summary = false,
        pqs_lowering_adopted_by_summary = false,
        public_default_behavior_changed = false,
        hamiltonian_schema_changed = false,
        gto_supplement_semantics_changed = false,
        raw_or_diagnostic_weights_promoted = false,
    )
end

function _cartesian_shellization_sequence_summary(
    sequence::_CartesianNestedShellSequence3D;
    source_kind::Symbol,
    route_family::Symbol,
    shellization_role::Symbol,
    parent_box::NTuple{3,UnitRange{Int}} = sequence.working_box,
    split_status::Symbol = :not_applicable,
    bond_axis::Union{Nothing,Symbol} = nothing,
    midpoint_slab_present::Bool = false,
    child_sequence_count::Int = 0,
    shared_shell_layer_count::Int = length(sequence.shell_layers),
    child_column_ranges::Tuple = (),
    contact_or_merge_status::Symbol = :not_applicable,
)
    return (
        object_kind = :cartesian_shellization_route_summary,
        status = :private_development_summary,
        private_development_only = true,
        route_family = route_family,
        source_kind = source_kind,
        shellization_role = shellization_role,
        shellization_stage = :route_neutral_spatial_planning,
        lowering_stage = :not_lowered_by_shellization_summary,
        parent_box = parent_box,
        working_box = sequence.working_box,
        core_column_range = sequence.core_column_range,
        core_retained_count = length(sequence.core_column_range),
        shell_layer_count = length(sequence.shell_layers),
        shell_layer_kinds = _cartesian_shellization_layer_kinds(sequence),
        shell_layer_column_ranges = _cartesian_shellization_layer_column_ranges(sequence),
        retained_dimension = size(sequence.coefficient_matrix, 2),
        support_count = length(sequence.support_indices),
        split_status = split_status,
        bond_axis = bond_axis,
        midpoint_slab_present = midpoint_slab_present,
        child_sequence_count = child_sequence_count,
        shared_shell_layer_count = shared_shell_layer_count,
        child_column_ranges = child_column_ranges,
        contact_or_merge_status = contact_or_merge_status,
        diagnostics = _cartesian_shellization_common_diagnostics(
            source_kind = source_kind,
            shellization_role = shellization_role,
        ),
    )
end

function _cartesian_shellization_route_summary(
    sequence::_CartesianNestedShellSequence3D;
    route_family::Symbol = :one_center_full_parent,
    source_kind::Symbol = :one_center_shell_sequence,
    shellization_role::Symbol = :atom_local_full_parent_shellization,
)
    return _cartesian_shellization_sequence_summary(
        sequence;
        route_family = route_family,
        source_kind = source_kind,
        shellization_role = shellization_role,
    )
end

function _cartesian_shellization_route_summary(
    source::_CartesianNestedBondAlignedDiatomicSource3D;
    route_family::Symbol = :bond_aligned_diatomic,
    source_kind::Symbol = :bond_aligned_diatomic_source,
    shellization_role::Symbol = :bond_aligned_diatomic_shellization,
)
    geometry = source.geometry
    split_status = geometry.did_split ? :split : :no_split
    midpoint_slab_present = !isnothing(source.midpoint_slab_column_range)
    contact_or_merge_status =
        geometry.did_split ?
        (midpoint_slab_present ? :split_children_with_midpoint_slab : :split_children_merged) :
        :no_split_single_child

    return _cartesian_shellization_sequence_summary(
        source.sequence;
        route_family = route_family,
        source_kind = source_kind,
        shellization_role = shellization_role,
        parent_box = geometry.parent_box,
        split_status = split_status,
        bond_axis = geometry.bond_axis,
        midpoint_slab_present = midpoint_slab_present,
        child_sequence_count = length(source.child_sequences),
        shared_shell_layer_count = length(source.shared_shell_layers),
        child_column_ranges = Tuple(source.child_column_ranges),
        contact_or_merge_status = contact_or_merge_status,
    )
end

function _cartesian_shellization_route_location_tuple(location)
    length(location) == 3 || throw(
        ArgumentError("route shellization request expects three Cartesian coordinates"),
    )
    return (Float64(location[1]), Float64(location[2]), Float64(location[3]))
end

function _cartesian_shellization_route_diatomic_axis(atom_locations; atol::Float64 = 1.0e-12)
    first_location = _cartesian_shellization_route_location_tuple(atom_locations[1])
    second_location = _cartesian_shellization_route_location_tuple(atom_locations[2])
    deltas = abs.(
        (
            second_location[1] - first_location[1],
            second_location[2] - first_location[2],
            second_location[3] - first_location[3],
        ),
    )
    active_axes = Tuple(axis for (axis, delta) in zip((:x, :y, :z), deltas) if delta > atol)
    return length(active_axes) == 1 ? only(active_axes) : nothing
end

function _cartesian_shellization_route_system_classification(system_metadata)
    atom_count = length(system_metadata.atom_symbols)
    if atom_count == 1
        return (
            system_classification = :one_center,
            system_classification_status = :explicit_atom_count_one,
            bond_axis = nothing,
        )
    elseif atom_count == 2
        bond_axis = _cartesian_shellization_route_diatomic_axis(system_metadata.atom_locations)
        if isnothing(bond_axis)
            return (
                system_classification = :pending_system_classification,
                system_classification_status = :diatomic_not_axis_aligned_by_metadata,
                bond_axis = nothing,
            )
        end
        return (
            system_classification = :bond_aligned_diatomic,
            system_classification_status = :explicit_two_atom_single_axis_separation,
            bond_axis = bond_axis,
        )
    end
    return (
        system_classification = :unsupported_general_multi_atom,
        system_classification_status = :general_multi_atom_shellization_not_planned,
        bond_axis = nothing,
    )
end

function _cartesian_shellization_route_configured_request(report)
    system_metadata = report.system_metadata
    recipe_metadata = report.recipe_metadata
    classification = _cartesian_shellization_route_system_classification(system_metadata)
    current_materialization_source =
        report.route_family == :white_lindsey_low_order ?
        :white_lindsey_one_center_seed :
        :pending_source_box_route_materializer

    return (
        object_kind = :cartesian_shellization_route_configured_request,
        status = :metadata_only_pending_materializer,
        private_development_only = true,
        route_family = report.route_family,
        route_kind = recipe_metadata.route_kind,
        route_shape = recipe_metadata.route_shape,
        atom_count = length(system_metadata.atom_symbols),
        atom_symbols = Tuple(system_metadata.atom_symbols),
        atom_locations = Tuple(system_metadata.atom_locations),
        nuclear_charges = Tuple(system_metadata.nuclear_charges),
        parent_axis_counts = system_metadata.parent_axis_counts,
        parent_axis_counts_source = system_metadata.parent_axis_counts_source,
        parent_box = system_metadata.parent_box,
        requested_shellization_stage = :route_neutral_spatial_planning,
        expected_next_materializer_status =
            :pending_route_configured_shellization_materializer,
        current_materialization_source,
        route_configured_shellization_consumed = false,
        constructs_basis = false,
        constructs_shell_sequence = false,
        constructs_fixed_block = false,
        system_classification = classification.system_classification,
        system_classification_status = classification.system_classification_status,
        bond_axis = classification.bond_axis,
    )
end

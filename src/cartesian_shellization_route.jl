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

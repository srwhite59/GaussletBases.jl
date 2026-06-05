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

function _cartesian_shellification_box_point_count(box::NTuple{3,UnitRange{Int}})
    return prod(length.(box))
end

function _cartesian_shellification_region3d(;
    order::Int,
    role::Symbol,
    box::NTuple{3,UnitRange{Int}},
    lowering_family::Symbol,
    materialization_dependency::Symbol,
    provenance,
    retained_count::Int,
    source_point_count::Int = _cartesian_shellification_box_point_count(box),
    support_count::Int = source_point_count,
    column_range = nothing,
    next_inner_box = nothing,
)
    return (;
        object_kind = :cartesian_shellification_region3d,
        order,
        role,
        box,
        box_shape = Tuple(length.(box)),
        source_point_count,
        lowering_family,
        lowering_status = :planned_not_lowered,
        materialization_dependency,
        retained_count,
        support_count,
        column_range,
        next_inner_box,
        provenance,
    )
end

function _cartesian_shellification_materialization_dependency_counts(regions)
    plan_lowerable_dependencies = (
        :plan_lowerable_complete_shell,
        :plan_lowerable_direct_core,
        :plan_lowerable_direct_slab,
        :plan_lowerable_shared_complete_shell,
    )
    source_backed_dependencies = (
        :source_backed_shared_shell_layer,
        :source_backed_child_sequence,
        :source_box_direct_in_source_backed_adapter,
    )
    return (;
        object_kind = :cartesian_shellification_materialization_dependency_counts,
        allowed_values = (
            plan_lowerable_dependencies...,
            source_backed_dependencies...,
        ),
        ordered_materialization_dependencies =
            Tuple(region.materialization_dependency for region in regions),
        plan_lowerable_region_count = count(
            region -> region.materialization_dependency in plan_lowerable_dependencies,
            regions,
        ),
        source_backed_region_count = count(
            region -> region.materialization_dependency in source_backed_dependencies,
            regions,
        ),
        plan_lowerable_complete_shell_count = count(
            region -> region.materialization_dependency == :plan_lowerable_complete_shell,
            regions,
        ),
        plan_lowerable_direct_core_count = count(
            region -> region.materialization_dependency == :plan_lowerable_direct_core,
            regions,
        ),
        plan_lowerable_direct_slab_count = count(
            region -> region.materialization_dependency == :plan_lowerable_direct_slab,
            regions,
        ),
        plan_lowerable_shared_complete_shell_count = count(
            region ->
                region.materialization_dependency ==
                :plan_lowerable_shared_complete_shell,
            regions,
        ),
        source_backed_shared_shell_layer_count = count(
            region -> region.materialization_dependency == :source_backed_shared_shell_layer,
            regions,
        ),
        source_backed_child_sequence_count = count(
            region -> region.materialization_dependency == :source_backed_child_sequence,
            regions,
        ),
        source_box_direct_adapter_region_count = count(
            region ->
                region.materialization_dependency ==
                :source_box_direct_in_source_backed_adapter,
            regions,
        ),
    )
end

function _cartesian_shellification_layer_provenance(layer)
    hasproperty(layer, :provenance) || throw(
        ArgumentError("diatomic shellification plan requires shell-layer provenance"),
    )
    return getproperty(layer, :provenance)
end

function _cartesian_shellification_layer_source_box(layer)
    provenance = _cartesian_shellification_layer_provenance(layer)
    hasproperty(provenance, :source_box) && return provenance.source_box
    hasproperty(provenance, :current_box) && return provenance.current_box
    throw(
        ArgumentError("diatomic shellification plan requires shared shell-layer source/current box provenance"),
    )
end

function _cartesian_shellification_layer_next_inner_box(layer)
    provenance = _cartesian_shellification_layer_provenance(layer)
    hasproperty(provenance, :next_inner_box) && return provenance.next_inner_box
    hasproperty(provenance, :inner_box) && return provenance.inner_box
    throw(
        ArgumentError("diatomic shellification plan requires shared shell-layer next-inner/inner box provenance"),
    )
end

function _cartesian_shellification_layer_source_point_count(layer)
    provenance = _cartesian_shellification_layer_provenance(layer)
    hasproperty(provenance, :source_point_count) && return Int(provenance.source_point_count)
    if hasproperty(layer, :owned_units)
        owned_units = getproperty(layer, :owned_units)
        if hasproperty(owned_units, :audit) &&
           hasproperty(owned_units.audit, :expected_support_count)
            return Int(owned_units.audit.expected_support_count)
        end
    end
    hasproperty(layer, :support_indices) && return length(getproperty(layer, :support_indices))
    throw(
        ArgumentError("diatomic shellification plan requires shared shell-layer support/source count"),
    )
end

function _cartesian_shellification_layer_lowering_family(
    layer::_CartesianNestedCompleteShell3D,
)
    return :white_lindsey_complete_shell
end

function _cartesian_shellification_layer_lowering_family(
    layer::_CartesianNestedEndcapPanelShellLayer3D,
)
    return :white_lindsey_endcap_panel_owned_shell
end

function _cartesian_shellification_layer_lowering_family(
    layer::_CartesianNestedProjectedQShellLayer3D,
)
    return :projected_q_shell_boundary_comx_product_modes
end

function _cartesian_shellification_layer_lowering_family(
    layer::_AbstractCartesianNestedShellLayer3D,
)
    return _cartesian_shellization_layer_kind(layer)
end

function _cartesian_shellification_layer_owned_unit_summary(layer)
    hasproperty(layer, :owned_units) || return nothing
    owned_units = getproperty(layer, :owned_units)
    audit = hasproperty(owned_units, :audit) ? owned_units.audit : nothing
    return (;
        source = :shared_shell_layer_owned_units,
        unit_count = hasproperty(owned_units, :units) ? length(owned_units.units) : nothing,
        support_contract =
            hasproperty(owned_units, :support_contract) ? owned_units.support_contract : nothing,
        coefficient_contract =
            hasproperty(owned_units, :coefficient_contract) ?
            owned_units.coefficient_contract :
            nothing,
        expected_support_count =
            isnothing(audit) ? nothing : audit.expected_support_count,
        owned_support_count = isnothing(audit) ? nothing : audit.owned_support_count,
        duplicate_count = isnothing(audit) ? nothing : audit.duplicate_count,
        missing_count = isnothing(audit) ? nothing : audit.missing_count,
        outside_count = isnothing(audit) ? nothing : audit.outside_count,
        retained_count = isnothing(audit) ? nothing : audit.retained_count,
        coverage_ok = isnothing(audit) ? nothing : audit.coverage_ok,
    )
end

function _cartesian_shellification_shared_shell_region3d(
    layer::_AbstractCartesianNestedShellLayer3D;
    order::Int,
    shared_layer_index::Int,
    column_range::UnitRange{Int},
)
    box = _cartesian_shellification_layer_source_box(layer)
    next_inner_box = _cartesian_shellification_layer_next_inner_box(layer)
    support_count = length(getproperty(layer, :support_indices))
    retained_count = size(getproperty(layer, :coefficient_matrix), 2)
    return _cartesian_shellification_region3d(
        order = order,
        role = :shared_outer_shell,
        box = box,
        next_inner_box = next_inner_box,
        lowering_family = _cartesian_shellification_layer_lowering_family(layer),
        materialization_dependency = :source_backed_shared_shell_layer,
        provenance = (;
            source = :_CartesianNestedBondAlignedDiatomicSource3D_shared_shell_layer,
            shared_layer_index,
            layer_kind = _cartesian_shellization_layer_kind(layer),
            layer_provenance = _cartesian_shellification_layer_provenance(layer),
            owned_unit_summary = _cartesian_shellification_layer_owned_unit_summary(layer),
            sequence_column_range = column_range,
        ),
        retained_count = retained_count,
        source_point_count = _cartesian_shellification_layer_source_point_count(layer),
        support_count = support_count,
        column_range = column_range,
    )
end

function _cartesian_shellification_child_sequence_region3d(
    source::_CartesianNestedBondAlignedDiatomicSource3D,
    child_sequence::_CartesianNestedShellSequence3D;
    order::Int,
    child_index::Int,
    column_range::UnitRange{Int},
)
    geometry = source.geometry
    geometry_child_box =
        geometry.did_split ? geometry.child_boxes[child_index] : nothing
    if geometry.did_split && child_sequence.working_box != geometry_child_box
        throw(
            ArgumentError(
                "diatomic shellification plan requires child sequence working_box to match geometry child_box",
            ),
        )
    end
    support_count = length(child_sequence.support_indices)
    retained_count = size(child_sequence.coefficient_matrix, 2)
    return _cartesian_shellification_region3d(
        order = order,
        role = geometry.did_split ? :atom_local_subtree : :unsplit_child_subtree,
        box = child_sequence.working_box,
        lowering_family = :white_lindsey_child_shell_sequence,
        materialization_dependency = :source_backed_child_sequence,
        provenance = (;
            source = :_CartesianNestedBondAlignedDiatomicSource3D_child_sequence,
            child_index,
            split_status = geometry.did_split ? :split : :no_split,
            child_sequence_working_box = child_sequence.working_box,
            geometry_child_box,
            geometry_child_box_matches_sequence =
                !geometry.did_split || child_sequence.working_box == geometry_child_box,
            sequence_core_column_range = child_sequence.core_column_range,
            merged_sequence_column_range = column_range,
        ),
        retained_count = retained_count,
        source_point_count = support_count,
        support_count = support_count,
        column_range = column_range,
    )
end

function _cartesian_shellification_midpoint_slab_region3d(
    source::_CartesianNestedBondAlignedDiatomicSource3D;
    order::Int,
)
    column_range = source.midpoint_slab_column_range
    isnothing(column_range) && throw(
        ArgumentError("diatomic midpoint slab region requires midpoint slab column range"),
    )
    box = source.geometry.shared_midpoint_box
    isnothing(box) && throw(
        ArgumentError("diatomic midpoint slab region requires geometry shared_midpoint_box"),
    )
    source_point_count = _cartesian_shellification_box_point_count(box)
    return _cartesian_shellification_region3d(
        order = order,
        role = :midpoint_slab,
        box = box,
        lowering_family = :direct_midpoint_slab,
        materialization_dependency = :source_box_direct_in_source_backed_adapter,
        provenance = (;
            source = :_CartesianNestedBondAlignedDiatomicSource3D_midpoint_slab,
            box_source = :geometry_shared_midpoint_box,
            column_range_source = :midpoint_slab_column_range,
            split_index = source.geometry.split_index,
            bond_axis = source.geometry.bond_axis,
        ),
        retained_count = length(column_range),
        source_point_count = source_point_count,
        support_count = source_point_count,
        column_range = column_range,
    )
end

function _cartesian_shellification_plan_direct_midpoint_slab_region3d(
    bundles::_CartesianNestedAxisBundles3D,
    box::NTuple{3,UnitRange{Int}};
    order_index::Int,
    bond_axis::Symbol,
    split_index::Union{Nothing,Int} = nothing,
    column_range = nothing,
)
    bond_axis in (:x, :y, :z) || throw(
        ArgumentError("direct midpoint slab region requires bond_axis = :x, :y, or :z"),
    )
    dims = _nested_axis_lengths(bundles)
    support_indices = _nested_box_support_indices(box..., dims)
    source_point_count = _cartesian_shellification_box_point_count(box)
    return (;
        object_kind = :cartesian_atom_outward_direct_midpoint_slab_region3d,
        role = :midpoint_slab,
        order = order_index,
        order_index,
        box,
        outer_box = box,
        inner_exclusion_box = nothing,
        box_shape = Tuple(length.(box)),
        region_kind = :direct_contact_or_midpoint_slab,
        lowering_family = :direct_midpoint_slab,
        materialization_dependency = :plan_lowerable_direct_slab,
        source_point_count = source_point_count,
        support_indices = Int.(support_indices),
        support_count = length(support_indices),
        retained_count = source_point_count,
        column_range,
        provenance = (;
            source = :atom_outward_direct_midpoint_slab_region,
            box_source = :geometry_input,
            column_range_source = isnothing(column_range) ? nothing : :oracle_only_column_range,
            split_index,
            bond_axis,
        ),
        source_backed = false,
        missing_independent_lowering_reason = nothing,
        retirement_target = :already_plan_lowered_region,
    )
end

function _cartesian_shellification_plan_shared_complete_shell_region3d(
    bundles::_CartesianNestedAxisBundles3D,
    outer_box::NTuple{3,UnitRange{Int}},
    inner_exclusion_box::NTuple{3,UnitRange{Int}};
    order_index::Int,
    bond_axis::Symbol,
    nside::Int,
    retention_policy::CartesianNestedCompleteShellRetentionContract,
    shared_shell_angular_resolution_scale::Real,
    retained_count = nothing,
    column_range = nothing,
)
    bond_axis in (:x, :y, :z) || throw(
        ArgumentError("shared complete-shell region requires bond_axis = :x, :y, or :z"),
    )
    shared_shell_angular_resolution_scale > 0 || throw(
        ArgumentError("shared complete-shell region requires positive angular resolution scale"),
    )
    dims = _nested_axis_lengths(bundles)
    support_count =
        _cartesian_shellification_box_point_count(outer_box) -
        _cartesian_shellification_box_point_count(inner_exclusion_box)
    return (;
        object_kind = :cartesian_atom_outward_shared_complete_shell_region3d,
        role = :regular_shared_molecular_shell,
        order = order_index,
        order_index,
        outer_box,
        box = outer_box,
        inner_exclusion_box,
        next_inner_box = inner_exclusion_box,
        box_shape = Tuple(length.(outer_box)),
        inner_box_shape = Tuple(length.(inner_exclusion_box)),
        region_kind = :complete_rectangular_shell,
        lowering_family = :white_lindsey_adaptive_complete_shell,
        materialization_dependency = :plan_lowerable_shared_complete_shell,
        lowering_parameters = (;
            retention_policy,
            bond_axis,
            nside,
            shared_shell_angular_resolution_scale =
                Float64(shared_shell_angular_resolution_scale),
            enforce_symmetric_odd = false,
        ),
        coverage_metadata = (;
            object_kind = :cartesian_shared_complete_shell_region_coverage,
            parent_dims = dims,
            support_count,
        ),
        source_point_count = support_count,
        support_count,
        retained_count,
        column_range,
        provenance = (;
            source = :atom_outward_shared_complete_shell_region,
            box_source = :geometry_input,
            retained_count_source = isnothing(retained_count) ? nothing : :oracle_only,
            column_range_source = isnothing(column_range) ? nothing : :oracle_only,
            bond_axis,
        ),
        source_backed = false,
        missing_independent_lowering_reason = nothing,
        retirement_target = :already_plan_lowered_region,
    )
end

function _cartesian_shellification_region_support_count(regions, role::Symbol)
    return sum(region.support_count for region in regions if region.role == role; init = 0)
end

function _cartesian_shellification_diatomic_source_coverage(
    source::_CartesianNestedBondAlignedDiatomicSource3D,
    regions,
)
    audit = _nested_source_contract_audit(source)
    region_support_count = sum(region.support_count for region in regions; init = 0)
    coverage_complete =
        audit.missing_row_count == 0 &&
        audit.ownership_unowned_row_count == 0 &&
        audit.ownership_multi_owned_row_count == 0 &&
        region_support_count == audit.support_count
    return (;
        object_kind = :cartesian_shellification_plan_coverage3d,
        status = coverage_complete ? :coverage_complete : :coverage_incomplete,
        parent_point_count = _cartesian_shellification_box_point_count(source.geometry.parent_box),
        expected_support_count = audit.expected_support_count,
        support_count = audit.support_count,
        shared_shell_support_count =
            _cartesian_shellification_region_support_count(regions, :shared_outer_shell),
        child_subtree_support_count =
            _cartesian_shellification_region_support_count(regions, :atom_local_subtree) +
            _cartesian_shellification_region_support_count(regions, :unsplit_child_subtree),
        midpoint_slab_support_count =
            _cartesian_shellification_region_support_count(regions, :midpoint_slab),
        region_support_count,
        region_support_count_matches_sequence = region_support_count == audit.support_count,
        missing_row_count = audit.missing_row_count,
        ownership_group_count_min = audit.ownership_group_count_min,
        ownership_group_count_max = audit.ownership_group_count_max,
        ownership_unowned_row_count = audit.ownership_unowned_row_count,
        ownership_multi_owned_row_count = audit.ownership_multi_owned_row_count,
        coverage_complete,
    )
end

function _cartesian_shellification_one_center_low_order_coverage(
    parent_box::NTuple{3,UnitRange{Int}},
    regions,
)
    shell_source_point_count = sum(
        region.source_point_count for region in regions
        if region.role == :low_order_complete_shell;
        init = 0,
    )
    direct_core_source_point_count = sum(
        region.source_point_count for region in regions
        if region.role == :direct_core;
        init = 0,
    )
    total_source_point_count = shell_source_point_count + direct_core_source_point_count
    parent_point_count = _cartesian_shellification_box_point_count(parent_box)
    return (;
        object_kind = :cartesian_shellification_plan_coverage3d,
        status =
            total_source_point_count == parent_point_count ?
            :coverage_complete :
            :coverage_incomplete,
        parent_point_count,
        shell_source_point_count,
        direct_core_source_point_count,
        total_source_point_count,
        missing_point_count = parent_point_count - total_source_point_count,
        coverage_complete = total_source_point_count == parent_point_count,
    )
end

function _cartesian_shellification_plan_one_center_low_order(
    parent_side_count::Int;
    nside::Int,
    route_family::Symbol = :white_lindsey_low_order,
)
    parent_side_count >= nside || throw(
        ArgumentError("one-center shellification plan requires parent side count >= nside"),
    )
    parent_box = ntuple(_ -> 1:parent_side_count, 3)
    shell_layer_count, core_side_count =
        _one_center_atomic_shell_layer_count(parent_side_count, nside)
    retention = _one_center_atomic_complete_shell_retention(nside)

    regions = []
    for layer_index in 1:shell_layer_count
        provenance = _one_center_atomic_layer_provenance(
            parent_box,
            layer_index,
            retention.shell_increment,
        )
        push!(
            regions,
            _cartesian_shellification_region3d(
                order = layer_index,
                role = :low_order_complete_shell,
                box = provenance.source_box,
                next_inner_box = provenance.next_inner_box,
                lowering_family = :white_lindsey_complete_shell,
                materialization_dependency = :plan_lowerable_complete_shell,
                provenance = provenance,
                source_point_count = provenance.source_point_count,
                retained_count = retention.shell_increment,
            ),
        )
    end

    core_offset = shell_layer_count
    direct_core_box = ntuple(3) do axis
        (first(parent_box[axis]) + core_offset):(last(parent_box[axis]) - core_offset)
    end
    direct_core_point_count = _cartesian_shellification_box_point_count(direct_core_box)
    direct_core_point_count == core_side_count^3 || throw(
        ArgumentError("one-center shellification plan produced inconsistent direct core size"),
    )
    push!(
        regions,
        _cartesian_shellification_region3d(
            order = shell_layer_count + 1,
            role = :direct_core,
            box = direct_core_box,
            lowering_family = :direct_product_core,
            materialization_dependency = :plan_lowerable_direct_core,
            provenance = (;
                source = :one_center_atomic_full_parent_final_direct_core,
                shell_layer_count,
                core_side_count,
            ),
            retained_count = direct_core_point_count,
            source_point_count = direct_core_point_count,
        ),
    )
    regions = Tuple(regions)
    materialization_dependency_counts =
        _cartesian_shellification_materialization_dependency_counts(regions)
    coverage =
        _cartesian_shellification_one_center_low_order_coverage(parent_box, regions)

    return (;
        object_kind = :cartesian_shellification_plan3d,
        status = :planned_metadata_only,
        private_development_only = true,
        source_kind = :one_center_atomic_full_parent_shellification_plan,
        route_family,
        system_classification = :one_center,
        shellification_role = :atom_local_full_parent_shellification,
        shellification_stage = :route_neutral_spatial_planning,
        lowering_stage = :not_lowered_by_shellification_plan,
        parent_box,
        working_box = parent_box,
        full_parent_working_box = true,
        nside,
        core_side_count,
        direct_core_box,
        direct_core_point_count,
        shell_layer_count,
        shell_retention = retention,
        regions,
        region_count = length(regions),
        ordered_region_roles = Tuple(region.role for region in regions),
        ordered_region_boxes = Tuple(region.box for region in regions),
        ordered_materialization_dependencies =
            materialization_dependency_counts.ordered_materialization_dependencies,
        materialization_dependency_counts,
        shell_region_count = count(region -> region.role == :low_order_complete_shell, regions),
        direct_core_region_count = count(region -> region.role == :direct_core, regions),
        retained_dimension = sum(region.retained_count for region in regions),
        coverage,
        diagnostics = (;
            source = :current_one_center_atomic_full_parent_shell_sequence_contract,
            private_development_only = true,
            route_neutral_spatial_planning = true,
            lowering_applied_by_plan = false,
            materialization_behavior_changed = false,
            public_default_behavior_changed = false,
            shellification_rewrite = false,
            pqs_production_source_box_materialization_claimed = false,
            mwg_ida_semantics_changed = false,
        ),
    )
end

function _cartesian_shellification_plan_bond_aligned_diatomic_low_order(
    source::_CartesianNestedBondAlignedDiatomicSource3D;
    route_family::Symbol = :white_lindsey_low_order,
)
    length(source.sequence.shell_layers) == length(source.shared_shell_layers) || throw(
        ArgumentError("diatomic shellification plan requires merged sequence shell layers to match source shared shell layers"),
    )
    length(source.sequence.layer_column_ranges) == length(source.shared_shell_layers) || throw(
        ArgumentError("diatomic shellification plan requires one merged sequence column range per shared shell layer"),
    )
    geometry = source.geometry
    split_status = geometry.did_split ? :split : :no_split
    midpoint_slab_present = !isnothing(source.midpoint_slab_column_range)
    contact_or_merge_status =
        geometry.did_split ?
        (midpoint_slab_present ? :split_children_with_midpoint_slab : :split_children_merged) :
        :no_split_single_child

    regions = NamedTuple[]
    for (shared_layer_index, (layer, column_range)) in
        enumerate(zip(source.shared_shell_layers, source.sequence.layer_column_ranges))
        push!(
            regions,
            _cartesian_shellification_shared_shell_region3d(
                layer;
                order = length(regions) + 1,
                shared_layer_index,
                column_range,
            ),
        )
    end
    child_region = function (child_index::Int)
        child_sequence = source.child_sequences[child_index]
        column_range = source.child_column_ranges[child_index]
        push!(
            regions,
            _cartesian_shellification_child_sequence_region3d(
                source,
                child_sequence;
                order = length(regions) + 1,
                child_index,
                column_range,
            ),
        )
    end
    if geometry.did_split && midpoint_slab_present
        length(source.child_sequences) == 2 || throw(
            ArgumentError("diatomic midpoint slab plan requires exactly two split child sequences"),
        )
        child_region(1)
        push!(
            regions,
            _cartesian_shellification_midpoint_slab_region3d(
                source;
                order = length(regions) + 1,
            ),
        )
        child_region(2)
    else
        for child_index in eachindex(source.child_sequences)
            child_region(child_index)
        end
        if midpoint_slab_present
            push!(
                regions,
                _cartesian_shellification_midpoint_slab_region3d(
                    source;
                    order = length(regions) + 1,
                ),
            )
        end
    end
    regions = Tuple(regions)
    materialization_dependency_counts =
        _cartesian_shellification_materialization_dependency_counts(regions)
    coverage = _cartesian_shellification_diatomic_source_coverage(source, regions)

    return (;
        object_kind = :cartesian_shellification_plan3d,
        status = :planned_metadata_only,
        private_development_only = true,
        source_kind = :bond_aligned_diatomic_active_source_shellification_plan,
        route_family,
        system_classification = :bond_aligned_diatomic,
        shellification_role = :bond_aligned_diatomic_shellification,
        shellification_stage = :route_neutral_spatial_planning,
        lowering_stage = :not_lowered_by_shellification_plan,
        parent_box = geometry.parent_box,
        working_box = source.sequence.working_box,
        full_parent_working_box = source.sequence.working_box == geometry.parent_box,
        split_working_box = geometry.working_box,
        geometry_working_box = geometry.working_box,
        bond_axis = geometry.bond_axis,
        midpoint = geometry.midpoint,
        split_index = geometry.split_index,
        split_status,
        did_split = geometry.did_split,
        shared_midpoint_box = geometry.shared_midpoint_box,
        midpoint_slab_present,
        midpoint_slab_column_range = source.midpoint_slab_column_range,
        contact_or_merge_status,
        nside = source.nside,
        child_shell_retention_contract = source.child_shell_retention_contract,
        shared_shell_retention_contract = source.shared_shell_retention_contract,
        shared_shell_layer_count = length(source.shared_shell_layers),
        child_sequence_count = length(source.child_sequences),
        child_column_ranges = Tuple(source.child_column_ranges),
        merged_sequence_core_column_range = source.sequence.core_column_range,
        merged_sequence_shell_layer_column_ranges =
            Tuple(source.sequence.layer_column_ranges),
        merged_sequence_shell_layer_count = length(source.sequence.shell_layers),
        merged_sequence_retained_dimension = size(source.sequence.coefficient_matrix, 2),
        merged_sequence_support_count = length(source.sequence.support_indices),
        regions,
        region_count = length(regions),
        ordered_region_roles = Tuple(region.role for region in regions),
        ordered_region_boxes = Tuple(region.box for region in regions),
        ordered_materialization_dependencies =
            materialization_dependency_counts.ordered_materialization_dependencies,
        materialization_dependency_counts,
        shared_shell_region_count =
            count(region -> region.role == :shared_outer_shell, regions),
        child_subtree_region_count = count(
            region -> region.role in (:atom_local_subtree, :unsplit_child_subtree),
            regions,
        ),
        midpoint_slab_region_count = count(region -> region.role == :midpoint_slab, regions),
        retained_dimension = size(source.sequence.coefficient_matrix, 2),
        support_count = length(source.sequence.support_indices),
        coverage,
        diagnostics = (;
            source = :current_bond_aligned_diatomic_source_contract,
            private_development_only = true,
            route_neutral_spatial_planning = true,
            describes_existing_active_source = true,
            lowering_applied_by_plan = false,
            materialization_behavior_changed = false,
            public_default_behavior_changed = false,
            shellification_rewrite = false,
            pqs_production_source_box_materialization_claimed = false,
            mwg_ida_semantics_changed = false,
        ),
    )
end

function _cartesian_shellification_plan(
    plan_kind::Symbol,
    parent_side_count::Int;
    nside::Int,
    route_family::Symbol = :white_lindsey_low_order,
)
    plan_kind == :one_center_full_parent_low_order || throw(
        ArgumentError("unsupported one-center shellification plan kind: $(plan_kind)"),
    )
    return _cartesian_shellification_plan_one_center_low_order(
        parent_side_count;
        nside,
        route_family,
    )
end

function _cartesian_shellification_plan(
    plan_kind::Symbol,
    source::_CartesianNestedBondAlignedDiatomicSource3D;
    route_family::Symbol = :white_lindsey_low_order,
)
    plan_kind == :bond_aligned_diatomic_active_source_low_order || throw(
        ArgumentError("unsupported source-backed diatomic shellification plan kind: $(plan_kind)"),
    )
    return _cartesian_shellification_plan_bond_aligned_diatomic_low_order(
        source;
        route_family,
    )
end

function _cartesian_atom_growth_shellification_materialization_counts(regions)
    return (;
        object_kind = :cartesian_atom_growth_shellification_materialization_counts,
        plan_lowerable_region_count =
            count(region -> region.independently_lowerable, regions),
        source_backed_region_count = count(region -> region.source_backed, regions),
        unsupported_region_count =
            count(region -> !region.independently_lowerable, regions),
        atom_local_child_plan_count = count(
            region ->
                region.materialization_dependency ==
                :plan_lowerable_atom_local_child_shellification,
            regions,
        ),
        direct_contact_slab_count = count(
            region -> region.materialization_dependency == :plan_lowerable_direct_slab,
            regions,
        ),
        shared_complete_shell_count = count(
            region ->
                region.materialization_dependency ==
                :plan_lowerable_shared_complete_shell,
            regions,
        ),
        unsupported_outer_mismatch_count = count(
            region ->
                region.materialization_dependency ==
                :unsupported_outer_mismatch_adjustment,
            regions,
        ),
    )
end

function _cartesian_atom_growth_shellification_scaffold_region(
    construction_region::_BondAlignedDiatomicAtomGrowthConstructionRegion3D;
    lowering_piece = nothing,
    lowering_family::Symbol,
    materialization_dependency::Symbol,
    independently_lowerable::Bool,
    missing_independent_lowering_reason = nothing,
    retirement_target::Symbol,
)
    support_count = length(construction_region.support_indices)
    if !isnothing(lowering_piece)
        lowering_piece.support_count == support_count || throw(
            ArgumentError(
                "atom-growth shellification scaffold support count does not match lowering piece",
            ),
        )
    end
    return (;
        object_kind = :cartesian_atom_growth_shellification_region_scaffold3d,
        role = construction_region.role,
        order = construction_region.order_index,
        order_index = construction_region.order_index,
        box = construction_region.box,
        outer_box = construction_region.box,
        inner_exclusion_box = construction_region.inner_exclusion_box,
        box_shape = Tuple(length.(construction_region.box)),
        support_indices = Int.(construction_region.support_indices),
        support_count,
        source_point_count = support_count,
        construction_metadata = construction_region.metadata,
        lowering_family,
        materialization_dependency,
        lowering_piece,
        lowering_piece_object_kind =
            isnothing(lowering_piece) ? nothing : lowering_piece.object_kind,
        lowering_piece_role = isnothing(lowering_piece) ? nothing : lowering_piece.role,
        lowering_piece_support_count =
            isnothing(lowering_piece) ? nothing : lowering_piece.support_count,
        support_count_matches_lowering_piece =
            !isnothing(lowering_piece) && lowering_piece.support_count == support_count,
        independently_lowerable,
        source_backed = false,
        missing_independent_lowering_reason,
        retirement_target,
        provenance = (;
            source = :_BondAlignedDiatomicAtomGrowthConstructionPlan3D,
            construction_role = construction_region.role,
            construction_order_index = construction_region.order_index,
            construction_region_metadata = construction_region.metadata,
        ),
    )
end

function _cartesian_shellification_plan_atom_growth_complete_rectangular_low_order(
    construction_plan::_BondAlignedDiatomicAtomGrowthConstructionPlan3D,
    bundles::_CartesianNestedAxisBundles3D;
    nside::Int,
    child_retention_policy::CartesianNestedCompleteShellRetentionContract,
    shared_retention_policy::CartesianNestedCompleteShellRetentionContract,
    reference_fudge_factor::Real,
    core_near_nucleus_protect_rows::Int,
    shared_shell_angular_resolution_scale::Real,
    route_family::Symbol = :white_lindsey_low_order,
)
    parent_box = construction_plan.anatomy.recipe.parent_box
    parent_dims = Tuple(length.(parent_box))
    _nested_axis_lengths(bundles) == parent_dims || throw(
        ArgumentError("atom-growth shellification scaffold requires bundles matching construction-plan parent dimensions"),
    )
    bond_axis = construction_plan.anatomy.recipe.bond_axis
    bond_axis in (:x, :y, :z) || throw(
        ArgumentError("atom-growth shellification scaffold requires bond_axis = :x, :y, or :z"),
    )

    regions = NamedTuple[]
    shared_complete_shell_regions = NamedTuple[]
    unsupported_regions = NamedTuple[]
    left_child_plan = nothing
    right_child_plan = nothing
    contact_cap_region = nothing

    for construction_region in construction_plan.regions
        if construction_region.role in (:left_atom_box, :right_atom_box)
            atom_side = construction_region.role == :left_atom_box ? :left : :right
            child_plan = _cartesian_shellification_plan_atom_local_child_low_order(
                bundles,
                construction_region.box;
                order_index = construction_region.order_index,
                atom_side,
                bond_axis,
                nside,
                retention_policy = child_retention_policy,
                reference_fudge_factor,
                core_near_nucleus_protect_rows,
            )
            construction_region.role == :left_atom_box ?
                (left_child_plan = child_plan) :
                (right_child_plan = child_plan)
            push!(
                regions,
                _cartesian_atom_growth_shellification_scaffold_region(
                    construction_region;
                    lowering_piece = child_plan,
                    lowering_family = :white_lindsey_atom_local_child_shellification,
                    materialization_dependency =
                        :plan_lowerable_atom_local_child_shellification,
                    independently_lowerable = true,
                    missing_independent_lowering_reason = nothing,
                    retirement_target = :already_plan_lowered_region,
                ),
            )
        elseif construction_region.role == :contact_cap
            contact_cap_region =
                _cartesian_shellification_plan_direct_midpoint_slab_region3d(
                    bundles,
                    construction_region.box;
                    order_index = construction_region.order_index,
                    bond_axis,
                    split_index = nothing,
                    column_range = nothing,
                )
            push!(
                regions,
                _cartesian_atom_growth_shellification_scaffold_region(
                    construction_region;
                    lowering_piece = contact_cap_region,
                    lowering_family = :direct_contact_cap_slab,
                    materialization_dependency = :plan_lowerable_direct_slab,
                    independently_lowerable = true,
                    missing_independent_lowering_reason = nothing,
                    retirement_target = :already_plan_lowered_region,
                ),
            )
        elseif construction_region.role == :regular_shared_molecular_shell
            isnothing(construction_region.inner_exclusion_box) && throw(
                ArgumentError("regular shared molecular shell requires an inner exclusion box"),
            )
            shared_region =
                _cartesian_shellification_plan_shared_complete_shell_region3d(
                    bundles,
                    construction_region.box,
                    construction_region.inner_exclusion_box;
                    order_index = construction_region.order_index,
                    bond_axis,
                    nside,
                    retention_policy = shared_retention_policy,
                    shared_shell_angular_resolution_scale,
                )
            push!(shared_complete_shell_regions, shared_region)
            push!(
                regions,
                _cartesian_atom_growth_shellification_scaffold_region(
                    construction_region;
                    lowering_piece = shared_region,
                    lowering_family = :white_lindsey_adaptive_complete_shell,
                    materialization_dependency =
                        :plan_lowerable_shared_complete_shell,
                    independently_lowerable = true,
                    missing_independent_lowering_reason = nothing,
                    retirement_target = :already_plan_lowered_region,
                ),
            )
        elseif construction_region.role == :outer_mismatch_shared_molecular_shell
            unsupported = _cartesian_atom_growth_shellification_scaffold_region(
                construction_region;
                lowering_piece = nothing,
                lowering_family = :outer_mismatch_adjustment_region,
                materialization_dependency = :unsupported_outer_mismatch_adjustment,
                independently_lowerable = false,
                missing_independent_lowering_reason =
                    :outer_mismatch_adjustment_lowering_not_implemented_in_pass_005,
                retirement_target = :future_outer_mismatch_plan_lowered_region,
            )
            push!(regions, unsupported)
            push!(unsupported_regions, unsupported)
        else
            throw(
                ArgumentError(
                    "unsupported atom-growth construction region role $(construction_region.role)",
                ),
            )
        end
    end

    isnothing(left_child_plan) && throw(
        ArgumentError("atom-growth shellification scaffold requires a left atom box"),
    )
    isnothing(right_child_plan) && throw(
        ArgumentError("atom-growth shellification scaffold requires a right atom box"),
    )
    regions = Tuple(regions)
    shared_complete_shell_regions = Tuple(shared_complete_shell_regions)
    unsupported_regions = Tuple(unsupported_regions)
    materialization_dependency_counts =
        _cartesian_atom_growth_shellification_materialization_counts(regions)
    assembly_core_order = isnothing(contact_cap_region) ?
        (:left_atom_box, :right_atom_box) :
        (:left_atom_box, :contact_cap, :right_atom_box)
    materialization_status = isempty(unsupported_regions) ?
        :ready_supported_complete_rectangular_subset :
        :blocked_unsupported_outer_mismatch

    return (;
        object_kind = :cartesian_atom_growth_shellification_plan3d,
        status = :planned_metadata_only,
        materialization_status,
        private_development_only = true,
        source_kind = :bond_aligned_diatomic_atom_growth_construction_plan,
        route_family,
        system_classification = :bond_aligned_diatomic,
        shellification_role = :atom_growth_complete_rectangular_low_order_scaffold,
        shellification_stage = :route_neutral_spatial_planning,
        lowering_stage = :not_lowered_by_shellification_plan,
        parent_box,
        working_box = parent_box,
        full_parent_working_box = construction_plan.support_coverage.coverage_ok,
        bond_axis,
        nside,
        child_retention_policy,
        shared_retention_policy,
        reference_fudge_factor = Float64(reference_fudge_factor),
        core_near_nucleus_protect_rows,
        shared_shell_angular_resolution_scale =
            Float64(shared_shell_angular_resolution_scale),
        spatial_policy_order = :atom_outward,
        construction_region_order = Tuple(construction_plan.region_order),
        assembly_shell_order = :outside_in,
        assembly_core_order,
        construction_plan,
        regions,
        region_count = length(regions),
        ordered_region_roles = Tuple(region.role for region in regions),
        ordered_region_boxes = Tuple(region.box for region in regions),
        ordered_materialization_dependencies =
            Tuple(region.materialization_dependency for region in regions),
        materialization_dependency_counts,
        left_child_plan,
        right_child_plan,
        contact_cap_region,
        shared_complete_shell_regions,
        unsupported_regions,
        unsupported_region_count = length(unsupported_regions),
        coverage = (;
            object_kind = :cartesian_atom_growth_shellification_coverage3d,
            expected_support_count =
                construction_plan.support_coverage.expected_support_count,
            region_support_count =
                construction_plan.support_coverage.region_support_count,
            covered_support_count =
                construction_plan.support_coverage.covered_support_count,
            duplicate_count = construction_plan.support_coverage.duplicate_count,
            missing_count = construction_plan.support_coverage.missing_count,
            outside_count = construction_plan.support_coverage.outside_count,
            coverage_complete = construction_plan.support_coverage.coverage_ok,
        ),
        diagnostics = (;
            source = :_BondAlignedDiatomicAtomGrowthConstructionPlan3D,
            private_development_only = true,
            atom_growth_construction_plan_authority = true,
            active_source_authority = false,
            active_source_oracle_comparison_run = false,
            route_neutral_spatial_planning = true,
            lowering_applied_by_plan = false,
            materialization_behavior_changed = false,
            public_default_behavior_changed = false,
            shellification_rewrite = false,
            pqs_production_source_box_materialization_claimed = false,
            mwg_ida_semantics_changed = false,
        ),
    )
end

function _cartesian_shellification_plan(
    plan_kind::Symbol,
    construction_plan::_BondAlignedDiatomicAtomGrowthConstructionPlan3D,
    bundles::_CartesianNestedAxisBundles3D;
    kwargs...,
)
    plan_kind == :bond_aligned_diatomic_atom_growth_complete_rectangular_low_order ||
        throw(
            ArgumentError(
                "unsupported atom-growth shellification plan kind: $(plan_kind)",
            ),
        )
    return _cartesian_shellification_plan_atom_growth_complete_rectangular_low_order(
        construction_plan,
        bundles;
        kwargs...,
    )
end

function _cartesian_atom_local_child_shell_boxes(
    outer_box::NTuple{3,UnitRange{Int}},
    nside::Int,
)
    shell_boxes = NamedTuple[]
    current_box = outer_box
    while minimum(length.(current_box)) > nside
        _nested_can_shrink_box(current_box) || break
        inner_box = _nested_inner_box(current_box)
        push!(
            shell_boxes,
            (;
                outer_box = current_box,
                inner_exclusion_box = inner_box,
                source_point_count =
                    _cartesian_shellification_box_point_count(current_box) -
                    _cartesian_shellification_box_point_count(inner_box),
            ),
        )
        current_box = inner_box
    end
    return (;
        shell_boxes = Tuple(shell_boxes),
        core_box = current_box,
    )
end

function _cartesian_atom_local_child_shell_region3d(;
    order_index::Int,
    shell_box,
    retention_policy::CartesianNestedCompleteShellRetentionContract,
    nside::Int,
    reference_fudge_factor::Float64,
)
    return (;
        object_kind = :cartesian_atom_local_child_shell_region3d,
        role = :atom_local_complete_shell,
        order_index,
        outer_box = shell_box.outer_box,
        inner_exclusion_box = shell_box.inner_exclusion_box,
        region_kind = :complete_shell_layer,
        lowering_family = :white_lindsey_adaptive_complete_shell,
        lowering_parameters = (;
            retention_policy,
            nside,
            reference_fudge_factor,
            enforce_symmetric_odd = false,
        ),
        source_point_count = shell_box.source_point_count,
        source_backed = false,
        missing_independent_lowering_reason = nothing,
        retirement_target = :already_plan_lowered_region,
    )
end

function _cartesian_atom_local_child_core_region3d(;
    order_index::Int,
    core_box::NTuple{3,UnitRange{Int}},
    support_indices::AbstractVector{Int},
    bond_axis::Symbol,
    nside::Int,
    retention_policy::CartesianNestedCompleteShellRetentionContract,
    reference_fudge_factor::Float64,
    core_near_nucleus_protect_rows::Int,
)
    return (;
        object_kind = :cartesian_atom_local_child_core_region3d,
        role = :atom_local_child_core,
        order_index,
        outer_box = core_box,
        inner_exclusion_box = nothing,
        region_kind = :nonuniform_direct_core,
        lowering_family = :diatomic_nonuniform_bond_axis_core,
        core_policy = :diatomic_nonuniform_bond_axis_core,
        lowering_parameters = (;
            retention_policy,
            bond_axis,
            nside,
            minimum_parallel_retain =
                _nested_diatomic_axis_minimum_retain(retention_policy, bond_axis),
            reference_fudge_factor,
            core_near_nucleus_protect_rows,
        ),
        support_indices = Int.(support_indices),
        support_count = length(support_indices),
        source_backed = false,
        missing_independent_lowering_reason = nothing,
        retirement_target = :already_plan_lowered_region,
    )
end

function _cartesian_shellification_plan_atom_local_child_low_order(
    bundles::_CartesianNestedAxisBundles3D,
    outer_box::NTuple{3,UnitRange{Int}};
    order_index::Int,
    atom_side::Symbol,
    bond_axis::Symbol,
    nside::Int,
    retention_policy::CartesianNestedCompleteShellRetentionContract,
    reference_fudge_factor::Real,
    core_near_nucleus_protect_rows::Int,
    route_family::Symbol = :white_lindsey_low_order,
)
    atom_side in (:left, :right) || throw(
        ArgumentError("atom-local child shellification plan requires atom_side = :left or :right"),
    )
    bond_axis in (:x, :y, :z) || throw(
        ArgumentError("atom-local child shellification plan requires bond_axis = :x, :y, or :z"),
    )
    nside >= 1 || throw(
        ArgumentError("atom-local child shellification plan requires positive nside"),
    )
    reference_fudge_factor > 0 || throw(
        ArgumentError("atom-local child shellification plan requires positive reference_fudge_factor"),
    )

    dims = _nested_axis_lengths(bundles)
    support_indices = _nested_box_support_indices(outer_box..., dims)
    shell_plan = _cartesian_atom_local_child_shell_boxes(outer_box, nside)
    shell_regions = Tuple(
        _cartesian_atom_local_child_shell_region3d(
            order_index = index,
            shell_box = shell_box,
            retention_policy = retention_policy,
            nside = nside,
            reference_fudge_factor = Float64(reference_fudge_factor),
        ) for (index, shell_box) in enumerate(shell_plan.shell_boxes)
    )
    core_support_indices = _nested_box_support_indices(shell_plan.core_box..., dims)
    core_region = _cartesian_atom_local_child_core_region3d(
        order_index = length(shell_regions) + 1,
        core_box = shell_plan.core_box,
        support_indices = core_support_indices,
        bond_axis = bond_axis,
        nside = nside,
        retention_policy = retention_policy,
        reference_fudge_factor = Float64(reference_fudge_factor),
        core_near_nucleus_protect_rows = core_near_nucleus_protect_rows,
    )

    return (;
        object_kind = :cartesian_atom_local_child_shellification_plan3d,
        status = :planned_metadata_only,
        private_development_only = true,
        source_kind = :atom_outward_atom_local_child_shellification_plan,
        route_family,
        role = :atom_local_subtree,
        order_index,
        atom_side,
        bond_axis,
        nside,
        retention_policy,
        core_policy = :diatomic_nonuniform_bond_axis_core,
        reference_fudge_factor = Float64(reference_fudge_factor),
        core_near_nucleus_protect_rows,
        outer_box,
        working_box = outer_box,
        shell_regions,
        shell_layer_count = length(shell_regions),
        core_region,
        core_box = shell_plan.core_box,
        support_indices = Int.(support_indices),
        support_count = length(support_indices),
        coverage_metadata = (;
            object_kind = :cartesian_atom_local_child_shellification_coverage,
            expected_support_count = _cartesian_shellification_box_point_count(outer_box),
            support_count = length(support_indices),
            coverage_complete =
                length(support_indices) ==
                _cartesian_shellification_box_point_count(outer_box),
        ),
        source_backed = false,
        missing_independent_lowering_reason = nothing,
        retirement_target = :already_plan_lowered_region,
        diagnostics = (;
            private_development_only = true,
            atom_outward_policy = true,
            active_source_authority = false,
            active_source_oracle_only = true,
            independently_lowerable = true,
            shared_shell_policy_changed = false,
            route_behavior_changed = false,
            mwg_ida_semantics_changed = false,
        ),
    )
end

function _cartesian_materialize_atom_local_child_shellification_low_order(
    plan,
    basis,
    bundles::_CartesianNestedAxisBundles3D;
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
    packet_kernel::Symbol = :factorized_direct,
    build_packet::Bool = false,
)
    plan.object_kind == :cartesian_atom_local_child_shellification_plan3d || throw(
        ArgumentError("atom-local child materializer requires a cartesian_atom_local_child_shellification_plan3d"),
    )
    plan.source_backed && throw(
        ArgumentError("atom-local child materializer requires an independently lowerable plan"),
    )
    plan.core_policy == :diatomic_nonuniform_bond_axis_core || throw(
        ArgumentError("atom-local child materializer requires core_policy = :diatomic_nonuniform_bond_axis_core"),
    )
    isnothing(term_coefficients) && throw(
        ArgumentError("atom-local child materializer requires explicit term coefficients"),
    )

    current_box = plan.outer_box
    shell_layers = _CartesianNestedCompleteShell3D[]
    for region in plan.shell_regions
        region.source_backed && throw(
            ArgumentError("atom-local child shell region must be independently lowerable"),
        )
        region.outer_box == current_box || throw(
            ArgumentError("atom-local child shell region order is inconsistent with the plan"),
        )
        inner_box = region.inner_exclusion_box
        adaptive_retention = _nested_diatomic_adaptive_shell_retention(
            basis,
            bundles,
            current_box,
            inner_box,
            plan.retention_policy;
            nside = plan.nside,
            reference_fudge_factor = plan.reference_fudge_factor,
        )
        push!(
            shell_layers,
            _nested_complete_rectangular_shell(
                bundles,
                inner_box...;
                retain_xy = adaptive_retention.retain_xy,
                retain_xz = adaptive_retention.retain_xz,
                retain_yz = adaptive_retention.retain_yz,
                retain_x_edge = adaptive_retention.retain_x_edge,
                retain_y_edge = adaptive_retention.retain_y_edge,
                retain_z_edge = adaptive_retention.retain_z_edge,
                x_fixed = (first(current_box[1]), last(current_box[1])),
                y_fixed = (first(current_box[2]), last(current_box[2])),
                z_fixed = (first(current_box[3]), last(current_box[3])),
                enforce_symmetric_odd = false,
                term_coefficients = term_coefficients,
                packet_kernel = packet_kernel,
                verify_factorized_reconstruction = false,
            ),
        )
        current_box = inner_box
    end
    current_box == plan.core_box || throw(
        ArgumentError("atom-local child shell plan core box is inconsistent with shell regions"),
    )
    core_block = _nested_bond_aligned_diatomic_nonuniform_core_block(
        basis,
        bundles,
        plan.core_box;
        bond_axis = plan.bond_axis,
        nside = plan.nside,
        minimum_parallel_retain =
            _nested_diatomic_axis_minimum_retain(plan.retention_policy, plan.bond_axis),
        reference_fudge_factor = plan.reference_fudge_factor,
        core_near_nucleus_protect_rows = plan.core_near_nucleus_protect_rows,
    )
    return _nested_shell_sequence_from_core_block(
        bundles,
        core_block.support_indices,
        core_block.coefficient_matrix,
        shell_layers;
        term_coefficients = term_coefficients,
        packet_kernel = packet_kernel,
        build_packet = build_packet,
        verify_factorized_reconstruction = false,
    )
end

function _cartesian_materialize_shared_complete_shell_region(
    region,
    basis,
    bundles::_CartesianNestedAxisBundles3D;
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
    packet_kernel::Symbol = :factorized_direct,
)
    region.object_kind == :cartesian_atom_outward_shared_complete_shell_region3d ||
        throw(
            ArgumentError("shared complete-shell materializer requires a cartesian_atom_outward_shared_complete_shell_region3d"),
        )
    region.source_backed && throw(
        ArgumentError("shared complete-shell materializer requires an independently lowerable region"),
    )
    region.materialization_dependency == :plan_lowerable_shared_complete_shell || throw(
        ArgumentError("shared complete-shell materializer requires plan-lowerable dependency"),
    )
    region.lowering_family == :white_lindsey_adaptive_complete_shell || throw(
        ArgumentError("shared complete-shell materializer requires adaptive complete-shell lowering"),
    )
    isnothing(term_coefficients) && throw(
        ArgumentError("shared complete-shell materializer requires explicit term coefficients"),
    )

    adaptive_retention = _nested_diatomic_adaptive_shell_retention(
        basis,
        bundles,
        region.outer_box,
        region.inner_exclusion_box,
        region.lowering_parameters.retention_policy;
        nside = region.lowering_parameters.nside,
        shared_shell_angular_resolution_scale =
            region.lowering_parameters.shared_shell_angular_resolution_scale,
    )
    shell = _nested_complete_rectangular_shell(
        bundles,
        region.inner_exclusion_box...;
        retain_xy = adaptive_retention.retain_xy,
        retain_xz = adaptive_retention.retain_xz,
        retain_yz = adaptive_retention.retain_yz,
        retain_x_edge = adaptive_retention.retain_x_edge,
        retain_y_edge = adaptive_retention.retain_y_edge,
        retain_z_edge = adaptive_retention.retain_z_edge,
        x_fixed = (first(region.outer_box[1]), last(region.outer_box[1])),
        y_fixed = (first(region.outer_box[2]), last(region.outer_box[2])),
        z_fixed = (first(region.outer_box[3]), last(region.outer_box[3])),
        enforce_symmetric_odd = false,
        term_coefficients = term_coefficients,
        packet_kernel = packet_kernel,
        verify_factorized_reconstruction = false,
    )
    length(shell.support_indices) == region.support_count || throw(
        ArgumentError("shared complete-shell materializer produced support count inconsistent with the region"),
    )
    if !isnothing(region.retained_count)
        size(shell.coefficient_matrix, 2) == region.retained_count || throw(
            ArgumentError("shared complete-shell materializer produced retained count inconsistent with the region"),
        )
    end
    return shell
end

function _cartesian_materialize_split_complete_rectangular_shellification_low_order(
    left_child_plan,
    right_child_plan,
    shared_shell_regions,
    basis,
    bundles::_CartesianNestedAxisBundles3D;
    midpoint_region = nothing,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
    packet_kernel::Symbol = :factorized_direct,
)
    isnothing(term_coefficients) && throw(
        ArgumentError("split complete-rectangular materializer requires explicit term coefficients"),
    )
    left_child_plan.source_backed && throw(
        ArgumentError("left atom-local child plan must be independently lowerable"),
    )
    right_child_plan.source_backed && throw(
        ArgumentError("right atom-local child plan must be independently lowerable"),
    )

    left_sequence = _cartesian_materialize_atom_local_child_shellification_low_order(
        left_child_plan,
        basis,
        bundles;
        term_coefficients,
        packet_kernel,
        build_packet = false,
    )
    right_sequence = _cartesian_materialize_atom_local_child_shellification_low_order(
        right_child_plan,
        basis,
        bundles;
        term_coefficients,
        packet_kernel,
        build_packet = false,
    )
    midpoint_data =
        isnothing(midpoint_region) ?
        nothing :
        _cartesian_materialize_direct_box_region(midpoint_region, bundles)
    shared_layers = [
        _cartesian_materialize_shared_complete_shell_region(
            region,
            basis,
            bundles;
            term_coefficients,
            packet_kernel,
        ) for region in shared_shell_regions
    ]

    core_support_blocks = Vector{Vector{Int}}()
    core_coefficient_blocks = _CartesianCoefficientMap[]
    push!(core_support_blocks, left_sequence.support_indices)
    push!(core_coefficient_blocks, left_sequence.coefficient_matrix)
    if !isnothing(midpoint_data)
        push!(core_support_blocks, midpoint_data.support_indices)
        push!(core_coefficient_blocks, midpoint_data.coefficient_matrix)
    end
    push!(core_support_blocks, right_sequence.support_indices)
    push!(core_coefficient_blocks, right_sequence.coefficient_matrix)
    core_support = vcat(core_support_blocks...)
    core_coefficients = _nested_hcat_coefficient_maps(core_coefficient_blocks)
    sequence = _nested_shell_sequence_from_core_block(
        bundles,
        core_support,
        core_coefficients,
        shared_layers;
        term_coefficients,
        packet_kernel,
        verify_factorized_reconstruction = false,
    )

    column_start = first(sequence.core_column_range)
    left_columns = size(left_sequence.coefficient_matrix, 2)
    left_column_range = column_start:(column_start + left_columns - 1)
    column_start = last(left_column_range) + 1
    midpoint_column_range = nothing
    if !isnothing(midpoint_data)
        slab_columns = size(midpoint_data.coefficient_matrix, 2)
        midpoint_column_range = column_start:(column_start + slab_columns - 1)
        column_start = last(midpoint_column_range) + 1
    end
    right_columns = size(right_sequence.coefficient_matrix, 2)
    right_column_range = column_start:(column_start + right_columns - 1)

    return (;
        object_kind = :cartesian_atom_outward_split_complete_rectangular_materialization,
        status = :materialized_split_complete_rectangular_low_order,
        private_development_only = true,
        active_source_authority = false,
        route_behavior_changed = false,
        left_child_plan,
        right_child_plan,
        midpoint_region,
        shared_shell_regions = Tuple(shared_shell_regions),
        child_sequences = (left_sequence, right_sequence),
        midpoint_slab_data = midpoint_data,
        shared_shell_layers = Tuple(shared_layers),
        sequence,
        child_column_ranges = (left_column_range, right_column_range),
        midpoint_slab_column_range = midpoint_column_range,
    )
end

function _cartesian_materialize_atom_growth_complete_rectangular_shellification_low_order(
    plan,
    basis,
    bundles::_CartesianNestedAxisBundles3D;
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
    packet_kernel::Symbol = :factorized_direct,
)
    plan.object_kind == :cartesian_atom_growth_shellification_plan3d || throw(
        ArgumentError("atom-growth materializer requires a cartesian_atom_growth_shellification_plan3d"),
    )
    if !isempty(plan.unsupported_regions)
        return (;
            object_kind = :cartesian_atom_growth_shellification_materialization_result,
            status = :blocked_unsupported_regions,
            materialization_status = plan.materialization_status,
            blocked_reason =
                :outer_mismatch_adjustment_lowering_not_implemented_in_pass_005,
            unsupported_regions = plan.unsupported_regions,
            unsupported_region_count = length(plan.unsupported_regions),
            sequence = nothing,
            assembly = nothing,
            active_source_authority = false,
            route_behavior_changed = false,
        )
    end
    isnothing(term_coefficients) && throw(
        ArgumentError("atom-growth materializer requires explicit term coefficients"),
    )
    assembly = _cartesian_materialize_split_complete_rectangular_shellification_low_order(
        plan.left_child_plan,
        plan.right_child_plan,
        plan.shared_complete_shell_regions,
        basis,
        bundles;
        midpoint_region = plan.contact_cap_region,
        term_coefficients,
        packet_kernel,
    )
    return (;
        object_kind = :cartesian_atom_growth_shellification_materialization_result,
        status = :materialized_supported_complete_rectangular_low_order,
        materialization_status = plan.materialization_status,
        blocked_reason = nothing,
        unsupported_regions = (),
        unsupported_region_count = 0,
        sequence = assembly.sequence,
        assembly,
        active_source_authority = false,
        route_behavior_changed = false,
    )
end

function _cartesian_materialize_direct_box_region(
    region,
    bundles::_CartesianNestedAxisBundles3D,
)
    if region.role == :direct_core
        if region.materialization_dependency != :plan_lowerable_direct_core
            throw(
                ArgumentError("direct-core region must be plan-lowerable before direct-box materialization"),
            )
        end
        if region.lowering_family != :direct_product_core
            throw(ArgumentError("direct-core region must use :direct_product_core lowering"))
        end
    elseif region.role == :midpoint_slab
        if region.materialization_dependency ==
           :source_box_direct_in_source_backed_adapter
            nothing
        elseif region.materialization_dependency == :plan_lowerable_direct_slab
            if hasproperty(region, :source_backed) && region.source_backed
                throw(
                    ArgumentError("plan-lowerable midpoint-slab region must not be source-backed"),
                )
            end
        else
            throw(
                ArgumentError("midpoint-slab region must be source-backed direct or plan-lowerable direct"),
            )
        end
        if region.lowering_family != :direct_midpoint_slab
            throw(
                ArgumentError("midpoint-slab region must use :direct_midpoint_slab lowering"),
            )
        end
    else
        throw(
            ArgumentError("direct-box region materialization does not support role $(region.role)"),
        )
    end

    direct_box = _nested_direct_box_coefficients(bundles, region.box)
    length(direct_box.support_indices) == region.support_count || throw(
        ArgumentError("direct-box materialization produced support count inconsistent with the region"),
    )
    length(direct_box.support_indices) == region.retained_count || throw(
        ArgumentError("direct-box materialization produced retained count inconsistent with the region"),
    )
    size(direct_box.coefficient_matrix, 2) == region.retained_count || throw(
        ArgumentError("direct-box materialization produced coefficient columns inconsistent with the region"),
    )
    return direct_box
end

function _cartesian_materialize_shellification_low_order(
    plan,
    bundle::_MappedOrdinaryGausslet1DBundle;
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    packet_kernel::Symbol = :factorized_direct,
    verify_factorized_reconstruction::Bool = true,
)
    term_coefficients = _one_center_atomic_term_coefficients(bundle, expansion)
    bundles = _CartesianNestedAxisBundles3D(bundle, bundle, bundle)
    return _cartesian_materialize_shellification_low_order(
        plan,
        bundles;
        packet_kernel,
        term_coefficients,
        verify_factorized_reconstruction,
    )
end

function _cartesian_materialize_shellification_low_order(
    plan,
    bundles::_CartesianNestedAxisBundles3D;
    packet_kernel::Symbol = :factorized_direct,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
    verify_factorized_reconstruction::Bool = true,
)
    plan.object_kind == :cartesian_shellification_plan3d || throw(
        ArgumentError("low-order shellification materializer requires a cartesian_shellification_plan3d"),
    )
    plan.system_classification == :one_center || throw(
        ArgumentError("low-order shellification materializer currently supports only one-center plans"),
    )
    plan.source_kind == :one_center_atomic_full_parent_shellification_plan || throw(
        ArgumentError("low-order shellification materializer currently supports only the full-parent one-center plan"),
    )
    plan.route_family == :white_lindsey_low_order || throw(
        ArgumentError("low-order shellification materializer requires route_family = :white_lindsey_low_order"),
    )

    dims = _nested_axis_lengths(bundles)
    Tuple(length.(plan.parent_box)) == dims || throw(
        ArgumentError("low-order shellification materializer requires plan parent_box dimensions to match axis bundles"),
    )
    plan.full_parent_working_box || throw(
        ArgumentError("low-order one-center shellification materializer requires full-parent working box"),
    )
    regions = Tuple(plan.regions)
    !isempty(regions) || throw(
        ArgumentError("low-order shellification materializer requires at least one plan region"),
    )
    last(regions).role == :direct_core || throw(
        ArgumentError("low-order shellification materializer requires the final plan region to be :direct_core"),
    )
    all(region.role == :low_order_complete_shell for region in regions[1:(end - 1)]) ||
        throw(
            ArgumentError("low-order one-center shellification materializer supports only :low_order_complete_shell regions before the direct core"),
        )

    retention = plan.shell_retention
    shell_layers = _CartesianNestedCompleteShell3D[]
    for region in regions[1:(end - 1)]
        isnothing(region.next_inner_box) && throw(
            ArgumentError("low-order shell region materialization requires next_inner_box"),
        )
        shell = _nested_complete_rectangular_shell(
            bundles,
            region.next_inner_box...;
            retain_xy = retention.retain_xy,
            retain_xz = retention.retain_xz,
            retain_yz = retention.retain_yz,
            retain_x_edge = retention.retain_x_edge,
            retain_y_edge = retention.retain_y_edge,
            retain_z_edge = retention.retain_z_edge,
            x_fixed = (first(region.box[1]), last(region.box[1])),
            y_fixed = (first(region.box[2]), last(region.box[2])),
            z_fixed = (first(region.box[3]), last(region.box[3])),
            packet_kernel = packet_kernel,
            term_coefficients = term_coefficients,
            verify_factorized_reconstruction = verify_factorized_reconstruction,
        )
        shell.provenance == region.provenance || throw(
            ArgumentError("low-order shell region materialization produced provenance inconsistent with the plan"),
        )
        push!(shell_layers, shell)
    end

    direct_core_region = last(regions)
    direct_core = _cartesian_materialize_direct_box_region(direct_core_region, bundles)
    return _nested_shell_sequence_from_core_block(
        bundles,
        direct_core.support_indices,
        direct_core.coefficient_matrix,
        shell_layers;
        packet_kernel,
        term_coefficients,
        verify_factorized_reconstruction,
    )
end

function _cartesian_assert_diatomic_plan_describes_source(
    plan,
    source::_CartesianNestedBondAlignedDiatomicSource3D,
)
    plan.object_kind == :cartesian_shellification_plan3d || throw(
        ArgumentError("source-backed diatomic shellification materializer requires a cartesian_shellification_plan3d"),
    )
    plan.system_classification == :bond_aligned_diatomic || throw(
        ArgumentError("source-backed diatomic shellification materializer requires a bond-aligned diatomic plan"),
    )
    plan.source_kind == :bond_aligned_diatomic_active_source_shellification_plan ||
        throw(
            ArgumentError("source-backed diatomic shellification materializer requires an active-source diatomic plan"),
        )
    plan.parent_box == source.geometry.parent_box || throw(
        ArgumentError("source-backed diatomic shellification plan parent_box does not match the source"),
    )
    plan.split_working_box == source.geometry.working_box || throw(
        ArgumentError("source-backed diatomic shellification plan split_working_box does not match the source"),
    )
    plan.working_box == source.sequence.working_box || throw(
        ArgumentError("source-backed diatomic shellification plan working_box does not match the merged source sequence"),
    )
    plan.bond_axis == source.geometry.bond_axis || throw(
        ArgumentError("source-backed diatomic shellification plan bond_axis does not match the source"),
    )
    plan.did_split == source.geometry.did_split || throw(
        ArgumentError("source-backed diatomic shellification plan split status does not match the source"),
    )
    plan.shared_shell_layer_count == length(source.shared_shell_layers) || throw(
        ArgumentError("source-backed diatomic shellification plan shared-shell count does not match the source"),
    )
    plan.child_sequence_count == length(source.child_sequences) || throw(
        ArgumentError("source-backed diatomic shellification plan child-sequence count does not match the source"),
    )
    plan.child_column_ranges == Tuple(source.child_column_ranges) || throw(
        ArgumentError("source-backed diatomic shellification plan child column ranges do not match the source"),
    )
    plan.midpoint_slab_column_range == source.midpoint_slab_column_range || throw(
        ArgumentError("source-backed diatomic shellification plan midpoint slab column range does not match the source"),
    )
    return nothing
end

function _cartesian_materialize_source_backed_shellification_low_order(
    plan,
    source::_CartesianNestedBondAlignedDiatomicSource3D;
    packet_kernel::Symbol = :factorized_direct,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
    verify_factorized_reconstruction::Bool = false,
)
    _cartesian_assert_diatomic_plan_describes_source(plan, source)
    length(source.sequence.shell_layers) == length(source.shared_shell_layers) || throw(
        ArgumentError("source-backed diatomic materialization requires merged sequence shell layers to be the source shared-shell layers"),
    )
    length(source.sequence.layer_column_ranges) == length(source.shared_shell_layers) || throw(
        ArgumentError("source-backed diatomic materialization requires one merged sequence column range per shared-shell layer"),
    )

    shared_regions =
        Tuple(region for region in plan.regions if region.role == :shared_outer_shell)
    length(shared_regions) == length(source.shared_shell_layers) || throw(
        ArgumentError("source-backed diatomic plan shared-shell regions do not match source shared-shell layers"),
    )
    for (region, layer, column_range) in
        zip(shared_regions, source.shared_shell_layers, source.sequence.layer_column_ranges)
        region.box == _cartesian_shellification_layer_source_box(layer) || throw(
            ArgumentError("source-backed diatomic shared-shell region box does not match source layer provenance"),
        )
        region.next_inner_box == _cartesian_shellification_layer_next_inner_box(layer) ||
            throw(
                ArgumentError("source-backed diatomic shared-shell region inner box does not match source layer provenance"),
            )
        region.column_range == column_range || throw(
            ArgumentError("source-backed diatomic shared-shell region column range does not match source sequence"),
        )
    end

    core_support_blocks = Vector{Vector{Int}}()
    core_coefficient_blocks = _CartesianCoefficientMap[]
    for region in plan.regions
        if region.role in (:atom_local_subtree, :unsplit_child_subtree)
            child_index = region.provenance.child_index
            child_sequence = source.child_sequences[child_index]
            region.column_range == source.child_column_ranges[child_index] || throw(
                ArgumentError("source-backed diatomic child region column range does not match source child range"),
            )
            region.box == child_sequence.working_box || throw(
                ArgumentError("source-backed diatomic child region box does not match source child sequence"),
            )
            push!(core_support_blocks, child_sequence.support_indices)
            push!(core_coefficient_blocks, child_sequence.coefficient_matrix)
        elseif region.role == :midpoint_slab
            region.column_range == source.midpoint_slab_column_range || throw(
                ArgumentError("source-backed diatomic midpoint region column range does not match source midpoint range"),
            )
            region.box == source.geometry.shared_midpoint_box || throw(
                ArgumentError("source-backed diatomic midpoint region box does not match source midpoint box"),
            )
            slab_data =
                _cartesian_materialize_direct_box_region(region, source.axis_bundles)
            push!(core_support_blocks, slab_data.support_indices)
            push!(core_coefficient_blocks, slab_data.coefficient_matrix)
        elseif region.role != :shared_outer_shell
            throw(
                ArgumentError("source-backed diatomic materializer does not support region role $(region.role)"),
            )
        end
    end
    isempty(core_support_blocks) && throw(
        ArgumentError("source-backed diatomic materializer requires at least one child or midpoint core block"),
    )
    core_support = vcat(core_support_blocks...)
    core_coefficients = _nested_hcat_coefficient_maps(core_coefficient_blocks)
    return _nested_shell_sequence_from_core_block(
        source.axis_bundles,
        core_support,
        core_coefficients,
        source.shared_shell_layers;
        packet_kernel,
        term_coefficients,
        verify_factorized_reconstruction,
    )
end
function _cartesian_shellification_plan_private_summary(plan)
    if plan.system_classification == :one_center
        return (;
            object_kind = :cartesian_shellification_plan_private_summary,
            status = plan.status,
            private_development_only = true,
            source_kind = plan.source_kind,
            route_family = plan.route_family,
            system_classification = plan.system_classification,
            shellification_stage = plan.shellification_stage,
            lowering_stage = plan.lowering_stage,
            region_count = plan.region_count,
            ordered_region_roles = plan.ordered_region_roles,
            ordered_materialization_dependencies =
                plan.ordered_materialization_dependencies,
            materialization_dependency_counts = plan.materialization_dependency_counts,
            plan_lowerable_region_count =
                plan.materialization_dependency_counts.plan_lowerable_region_count,
            source_backed_region_count =
                plan.materialization_dependency_counts.source_backed_region_count,
            source_box_direct_adapter_region_count =
                plan.materialization_dependency_counts.source_box_direct_adapter_region_count,
            shell_region_count = plan.shell_region_count,
            direct_core_region_count = plan.direct_core_region_count,
            retained_dimension = plan.retained_dimension,
            coverage_status = plan.coverage.status,
            coverage_complete = plan.coverage.coverage_complete,
        )
    elseif plan.system_classification == :bond_aligned_diatomic
        return (;
            object_kind = :cartesian_shellification_plan_private_summary,
            status = plan.status,
            private_development_only = true,
            source_kind = plan.source_kind,
            route_family = plan.route_family,
            system_classification = plan.system_classification,
            shellification_stage = plan.shellification_stage,
            lowering_stage = plan.lowering_stage,
            split_status = plan.split_status,
            did_split = plan.did_split,
            midpoint_slab_present = plan.midpoint_slab_present,
            region_count = plan.region_count,
            ordered_region_roles = plan.ordered_region_roles,
            ordered_materialization_dependencies =
                plan.ordered_materialization_dependencies,
            materialization_dependency_counts = plan.materialization_dependency_counts,
            plan_lowerable_region_count =
                plan.materialization_dependency_counts.plan_lowerable_region_count,
            source_backed_region_count =
                plan.materialization_dependency_counts.source_backed_region_count,
            source_box_direct_adapter_region_count =
                plan.materialization_dependency_counts.source_box_direct_adapter_region_count,
            shared_shell_region_count = plan.shared_shell_region_count,
            child_subtree_region_count = plan.child_subtree_region_count,
            midpoint_slab_region_count = plan.midpoint_slab_region_count,
            shared_shell_layer_count = plan.shared_shell_layer_count,
            child_sequence_count = plan.child_sequence_count,
            retained_dimension = plan.retained_dimension,
            support_count = plan.support_count,
            coverage_status = plan.coverage.status,
            coverage_complete = plan.coverage.coverage_complete,
            region_support_count_matches_sequence =
                plan.coverage.region_support_count_matches_sequence,
        )
    end
    throw(
        ArgumentError(
            "private shellification plan summary does not support system_classification $(plan.system_classification)",
        ),
    )
end

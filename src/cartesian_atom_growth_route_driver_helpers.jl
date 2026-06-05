# Private route-report helpers for the atom-growth Cartesian diatomic driver path.
# This file is included after the generic route-driver helpers and White-Lindsey seed helpers.

function _pqs_source_box_route_driver_atom_growth_unit_key(prefix::Symbol, index::Int)
    return index == 1 ? prefix : Symbol(string(prefix), "_", index)
end

function _pqs_source_box_route_driver_atom_growth_unit_field(
    object,
    field::Symbol,
    default = nothing,
)
    return hasproperty(object, field) ? getproperty(object, field) : default
end

function _pqs_source_box_route_driver_atom_growth_lowering_parameters(region)
    piece = region.lowering_piece
    if !isnothing(piece) && hasproperty(piece, :lowering_parameters)
        return piece.lowering_parameters
    end
    return (;
        lowering_piece_object_kind = region.lowering_piece_object_kind,
        lowering_piece_role = region.lowering_piece_role,
        lowering_piece_support_count = region.lowering_piece_support_count,
    )
end

function _pqs_source_box_route_driver_cpb_codimension(
    box::NTuple{3,UnitRange{Int}},
)
    return count(==(1), length.(box))
end

function _pqs_source_box_route_driver_cpb_geometry_kind(codimension::Int)
    codimension == 0 && return :filled_volume_cpb
    codimension == 1 && return :facet_or_slab_cpb
    codimension == 2 && return :edge_cpb
    codimension == 3 && return :corner_cpb
    return :invalid_cpb_codimension
end

function _pqs_source_box_route_driver_cpb_descriptor(;
    box::NTuple{3,UnitRange{Int}},
    role::Symbol,
    cpb_family::Symbol,
    source::Symbol,
    support_count = nothing,
    metadata = nothing,
)
    codimension = _pqs_source_box_route_driver_cpb_codimension(box)
    return (;
        object_kind = :cartesian_coordinate_product_box3d,
        coordinate_product_box = true,
        intervals = box,
        box = box,
        dimensions = Tuple(length.(box)),
        codimension,
        geometry_kind = _pqs_source_box_route_driver_cpb_geometry_kind(codimension),
        role,
        cpb_family,
        source,
        support_count,
        metadata,
    )
end

function _pqs_source_box_route_driver_complete_shell_condition(region)
    outer_box =
        _pqs_source_box_route_driver_atom_growth_unit_field(
            region,
            :outer_box,
            region.box,
        )
    inner_box =
        _pqs_source_box_route_driver_atom_growth_unit_field(
            region,
            :inner_exclusion_box,
            nothing,
        )
    if isnothing(inner_box)
        return (;
            complete_shell_supported = false,
            complete_shell_status = :unsupported_complete_shell,
            unsupported_reasons = (:missing_inner_exclusion_box,),
            outer_box,
            inner_box,
        )
    end

    reasons = Symbol[]
    for axis_index in 1:3
        outer_axis = outer_box[axis_index]
        inner_axis = inner_box[axis_index]
        if isempty(inner_axis)
            push!(reasons, :empty_inner_interval)
        elseif first(inner_axis) != first(outer_axis) + 1 ||
               last(inner_axis) != last(outer_axis) - 1
            push!(reasons, :inner_interval_not_single_layer_strict_interior)
        end
    end
    unsupported_reasons = Tuple(unique(reasons))
    return (;
        complete_shell_supported = isempty(unsupported_reasons),
        complete_shell_status =
            isempty(unsupported_reasons) ?
            :supported_complete_shell :
            :unsupported_complete_shell,
        unsupported_reasons,
        outer_box,
        inner_box,
    )
end

function _pqs_source_box_route_driver_complete_shell_stratum_role(
    axis_indices,
    sides,
    stratum_kind::Symbol,
)
    axis_symbols = (:x, :y, :z)
    parts = String[]
    for (axis_index, side) in zip(axis_indices, sides)
        push!(parts, string(axis_symbols[axis_index]))
        push!(parts, string(side))
    end
    push!(parts, string(stratum_kind))
    return Symbol(join(parts, "_"))
end

function _pqs_source_box_route_driver_complete_shell_stratum_box(
    outer_box,
    inner_box,
    axis_indices,
    sides,
)
    ranges = Vector{UnitRange{Int}}(undef, 3)
    for axis_index in 1:3
        side_index = findfirst(==(axis_index), axis_indices)
        if isnothing(side_index)
            ranges[axis_index] = inner_box[axis_index]
        else
            side = sides[side_index]
            boundary_point =
                side == :low ?
                first(outer_box[axis_index]) :
                last(outer_box[axis_index])
            ranges[axis_index] = boundary_point:boundary_point
        end
    end
    return (ranges[1], ranges[2], ranges[3])
end

function _pqs_source_box_route_driver_complete_shell_stratum_cpb(
    outer_box,
    inner_box,
    axis_indices,
    sides,
    cpb_family::Symbol,
    stratum_kind::Symbol,
)
    box = _pqs_source_box_route_driver_complete_shell_stratum_box(
        outer_box,
        inner_box,
        axis_indices,
        sides,
    )
    return _pqs_source_box_route_driver_cpb_descriptor(
        box = box,
        role =
            _pqs_source_box_route_driver_complete_shell_stratum_role(
                axis_indices,
                sides,
                stratum_kind,
            ),
        cpb_family = cpb_family,
        source = :white_lindsey_complete_shell_boundary_strata,
        support_count = prod(length.(box)),
        metadata = (;
            enumeration_policy =
                :white_lindsey_complete_shell_boundary_strata,
            shellification_authority_scope = :owned_support_only,
            lowering_geometry_only = true,
            boundary_axis_indices = axis_indices,
            boundary_sides = sides,
            coefficient_maps_materialized = false,
            operator_blocks_materialized = false,
            pair_operator_blocks_materialized = false,
            hamiltonian_data_materialized = false,
        ),
    )
end

function _pqs_source_box_route_driver_complete_shell_source_cpbs(region)
    condition =
        _pqs_source_box_route_driver_complete_shell_condition(region)
    condition.complete_shell_supported || return ()

    cpbs = NamedTuple[]
    sides = (:low, :high)
    for axis_index in 1:3
        for side in sides
            push!(
                cpbs,
                _pqs_source_box_route_driver_complete_shell_stratum_cpb(
                    condition.outer_box,
                    condition.inner_box,
                    (axis_index,),
                    (side,),
                    :facet_cpb,
                    :facet,
                ),
            )
        end
    end
    for axis_indices in ((1, 2), (1, 3), (2, 3))
        for first_side in sides
            for second_side in sides
                push!(
                    cpbs,
                    _pqs_source_box_route_driver_complete_shell_stratum_cpb(
                        condition.outer_box,
                        condition.inner_box,
                        axis_indices,
                        (first_side, second_side),
                        :edge_cpb,
                        :edge,
                    ),
                )
            end
        end
    end
    for x_side in sides
        for y_side in sides
            for z_side in sides
                push!(
                    cpbs,
                    _pqs_source_box_route_driver_complete_shell_stratum_cpb(
                        condition.outer_box,
                        condition.inner_box,
                        (1, 2, 3),
                        (x_side, y_side, z_side),
                        :corner_cpb,
                        :corner,
                    ),
                )
            end
        end
    end
    return Tuple(cpbs)
end

function _pqs_source_box_route_driver_cpb_family_counts(source_cpbs)
    return (;
        facet_cpb = count(cpb -> cpb.cpb_family == :facet_cpb, source_cpbs),
        edge_cpb = count(cpb -> cpb.cpb_family == :edge_cpb, source_cpbs),
        corner_cpb = count(cpb -> cpb.cpb_family == :corner_cpb, source_cpbs),
    )
end

function _pqs_source_box_route_driver_lowering_source_enumeration_metadata(
    region,
    source_cpbs,
)
    if region.lowering_family == :white_lindsey_adaptive_complete_shell
        condition =
            _pqs_source_box_route_driver_complete_shell_condition(region)
        cpb_family_counts =
            _pqs_source_box_route_driver_cpb_family_counts(source_cpbs)
        source_cpb_enumeration_status =
            condition.complete_shell_supported && length(source_cpbs) == 26 ?
            :explicit_complete_shell_boundary_strata :
            :planned_cpb_families_not_enumerated
        return (;
            source_cpb_enumeration_status,
            source_cpb_enumeration_reason =
                source_cpb_enumeration_status ==
                :explicit_complete_shell_boundary_strata ?
                nothing :
                condition.unsupported_reasons,
            enumeration_policy =
                :white_lindsey_complete_shell_boundary_strata,
            complete_shell_condition_status =
                condition.complete_shell_status,
            complete_shell_unsupported_reasons =
                condition.unsupported_reasons,
            complete_shell_cpb_family_counts = cpb_family_counts,
        )
    end
    return (;
        source_cpb_enumeration_status =
            isempty(source_cpbs) ?
            :planned_cpb_families_not_enumerated :
            :explicit_source_cpbs,
        source_cpb_enumeration_reason = nothing,
        enumeration_policy = nothing,
        complete_shell_condition_status = :not_applicable,
        complete_shell_unsupported_reasons = (),
        complete_shell_cpb_family_counts =
            _pqs_source_box_route_driver_cpb_family_counts(source_cpbs),
    )
end

function _pqs_source_box_route_driver_owned_support(region)
    inner_exclusion_box =
        _pqs_source_box_route_driver_atom_growth_unit_field(
            region,
            :inner_exclusion_box,
            nothing,
        )
    outer_box =
        _pqs_source_box_route_driver_atom_growth_unit_field(
            region,
            :outer_box,
            region.box,
        )
    outer_cpb = _pqs_source_box_route_driver_cpb_descriptor(
        box = outer_box,
        role = :owned_support_outer_extent,
        cpb_family = :owned_support_extent_cpb,
        source = :atom_growth_shellification_region,
    )
    inner_exclusion_cpb =
        isnothing(inner_exclusion_box) ?
        nothing :
        _pqs_source_box_route_driver_cpb_descriptor(
            box = inner_exclusion_box,
            role = :owned_support_inner_exclusion,
            cpb_family = :owned_support_exclusion_cpb,
            source = :atom_growth_shellification_region,
        )
    owned_support_kind =
        isnothing(inner_exclusion_box) ?
        :coordinate_product_owned_support :
        :shell_difference_owned_support
    return (;
        object_kind = :cartesian_owned_support_region3d,
        role = region.role,
        support_kind = owned_support_kind,
        owned_support_authority = :shellification_region,
        shellification_authority_scope = :owned_support_only,
        shellification_region_is_cpb = false,
        shellification_region_is_lowering_source = false,
        owned_support_is_cpb = isnothing(inner_exclusion_box),
        owned_support_is_coordinate_product = isnothing(inner_exclusion_box),
        owned_support_is_shell_difference = !isnothing(inner_exclusion_box),
        box_difference_is_cpb = false,
        outer_cpb,
        inner_exclusion_cpb,
        support_count = region.support_count,
        support_count_source = :atom_growth_shellification_region,
    )
end

function _pqs_source_box_route_driver_lowering_source_cpbs(region)
    piece = region.lowering_piece
    if region.lowering_family == :white_lindsey_adaptive_complete_shell
        return _pqs_source_box_route_driver_complete_shell_source_cpbs(region)
    end
    if !isnothing(piece) &&
       region.lowering_piece_object_kind == :cartesian_outer_mismatch_boundary_slab_set3d
        slab_pieces =
            _pqs_source_box_route_driver_atom_growth_unit_field(
                piece,
                :slab_pieces,
                (),
            )
        return Tuple(
            _pqs_source_box_route_driver_cpb_descriptor(
                box = slab_piece.box,
                role = slab_piece.role,
                cpb_family = :direct_boundary_slab_cpb,
                source = :outer_mismatch_boundary_slab_set_lowering,
                support_count = slab_piece.support_count,
                metadata =
                    _pqs_source_box_route_driver_atom_growth_unit_field(
                        slab_piece,
                        :metadata,
                        nothing,
                    ),
            ) for slab_piece in slab_pieces
        )
    end
    if region.materialization_dependency == :plan_lowerable_direct_slab
        return (
            _pqs_source_box_route_driver_cpb_descriptor(
                box = region.box,
                role = region.role,
                cpb_family = :direct_slab_cpb,
                source = :direct_slab_lowering,
                support_count = region.support_count,
            ),
        )
    end
    return ()
end

function _pqs_source_box_route_driver_lowering_source_families(region)
    if region.lowering_family == :outer_mismatch_boundary_slab_set
        return (:direct_boundary_slab_cpb,)
    end
    if region.lowering_family == :white_lindsey_adaptive_complete_shell
        return (:facet_cpb, :edge_cpb, :corner_cpb)
    end
    if region.lowering_family == :white_lindsey_atom_local_child_shellification
        return (:facet_cpb, :edge_cpb, :corner_cpb, :direct_core_cpb)
    end
    if region.materialization_dependency == :plan_lowerable_direct_slab
        return (:direct_slab_cpb,)
    end
    return (:pending_cpb_lowering_family,)
end

function _pqs_source_box_route_driver_lowering_recipe(region, source_cpbs)
    source_families =
        _pqs_source_box_route_driver_lowering_source_families(region)
    enumeration_metadata =
        _pqs_source_box_route_driver_lowering_source_enumeration_metadata(
            region,
            source_cpbs,
        )
    return (;
        object_kind = :cartesian_cpb_lowering_recipe,
        lowering_stage = :coordinate_product_box_lowering,
        shellification_policy = :atom_growth_complete_rectangular,
        atom_growth_is_shellification_policy = true,
        lowering_family = region.lowering_family,
        materialization_dependency = region.materialization_dependency,
        source_cpb_families = source_families,
        source_cpb_count = length(source_cpbs),
        enumeration_metadata...,
        owned_support_authority = :shellification_region,
        shellification_authority_scope = :owned_support_only,
        shellification_region_is_lowering_source = false,
        lowering_source_authority = :lowering_recipe_cpbs,
        white_lindsey_lowering_sources_are_cpbs =
            all(!=(:pending_cpb_lowering_family), source_families),
        coefficient_maps_materialized = false,
        operator_blocks_materialized = false,
        pair_operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        shellification_region_authority = false,
    )
end

function _pqs_source_box_route_driver_intermediate_retained_space(unit_key)
    return (;
        object_kind = :cartesian_intermediate_retained_space_contract,
        unit_key,
        construction_stage = :intermediate_retained_space,
        status = :deferred_not_materialized,
        retained_dimension_known = false,
        retained_count = nothing,
        retained_range = nothing,
        direct_shellification_region_alias = false,
    )
end

function _pqs_source_box_route_driver_shell_realization_contract(unit_key)
    return (;
        object_kind = :cartesian_shell_realization_contract,
        unit_key,
        construction_stage = :shell_realization,
        status = :deferred_or_trivial_for_white_lindsey_lowering,
        shell_realization_materialized = false,
        shell_row_oracle_authority = false,
    )
end

function _pqs_source_box_route_driver_final_retained_unit_contract(
    unit_key,
    unit_role,
)
    return (;
        object_kind = :cartesian_final_retained_unit_contract,
        unit_key,
        unit_role,
        construction_stage = :final_retained_unit,
        status = :deferred_until_materialization,
        downstream_of_lowering = true,
        direct_shellification_region_alias = false,
        retained_count_known = false,
        retained_range_known = false,
        pair_planning_input = true,
    )
end

function _pqs_source_box_route_driver_pqs_lowering_prototype(unit)
    unit.owned_support.owned_support_is_shell_difference || return nothing
    source_box = unit.owned_support.outer_cpb.box
    source_cpb_support_count = prod(length.(source_box))
    owned_support_count = unit.support_count
    source_cpb = _pqs_source_box_route_driver_cpb_descriptor(
        box = source_box,
        role = :pqs_filled_source_cpb,
        cpb_family = :filled_source_cpb,
        source = :metadata_only_pqs_lowering_prototype,
        support_count = source_cpb_support_count,
    )
    return (;
        object_kind = :cartesian_pqs_lowering_metadata_prototype,
        status = :metadata_only_planned,
        private_development_only = true,
        unit_key = unit.unit_key,
        unit_role = unit.unit_role,
        shellification_policy = unit.shellification_policy,
        shellification_authority_scope = :owned_support_only,
        owned_support_authority = :shellification_region,
        owned_support = unit.owned_support,
        owned_support_is_cpb = unit.owned_support.owned_support_is_cpb,
        owned_support_is_shell_difference =
            unit.owned_support.owned_support_is_shell_difference,
        owned_support_count,
        owned_support_count_source = :shellification_region,
        source_cpb,
        source_cpb_support_count,
        source_cpb_support_count_source = :filled_coordinate_product_box,
        source_cpbs = (source_cpb,),
        source_cpb_count = 1,
        lowering_source_authority = :pqs_lowering_recipe_filled_source_cpb,
        shellification_region_is_lowering_source = false,
        lowering_recipe = :pqs_filled_source_cpb,
        lowering_recipe_contract = (;
            object_kind = :cartesian_cpb_lowering_recipe,
            lowering_stage = :coordinate_product_box_lowering,
            lowering_family = :projected_q_shell,
            lowering_recipe = :pqs_filled_source_cpb,
            source_cpb_families = (:filled_source_cpb,),
            retained_rule = :boundary_comx_product_mode_selection,
            shellification_region_authority = false,
            shellification_authority_scope = :owned_support_only,
            lowering_source_authority = :pqs_lowering_recipe_filled_source_cpb,
        ),
        intermediate_retained_space = (;
            object_kind =
                :pqs_boundary_comx_product_intermediate_retained_space,
            construction_stage = :intermediate_retained_space,
            status = :planned_deferred,
            retained_rule = :boundary_comx_product_mode_selection,
            selected_modes = :boundary_comx_product_modes,
            source_space_operator_blocks_planned = true,
            retained_dimension_known = false,
        ),
        shell_realization = (;
            object_kind = :pqs_shell_projection_lowdin_realization,
            construction_stage = :shell_realization,
            status = :planned_deferred,
            shell_projection_planned = true,
            lowdin_cleanup_planned = true,
            materialized = false,
        ),
        dense_parent_space_operators_are_algorithm = false,
        shell_row_operator_algorithm = false,
        coefficient_maps_materialized = false,
        operator_blocks_materialized = false,
        pair_operator_blocks_materialized = false,
        public_route_adoption = false,
        default_route_behavior_changed = false,
    )
end

function _pqs_source_box_route_driver_atom_growth_source_descriptor(region)
    piece = region.lowering_piece
    if !isnothing(piece) &&
       region.lowering_piece_object_kind == :cartesian_outer_mismatch_boundary_slab_set3d
        slab_pieces =
            _pqs_source_box_route_driver_atom_growth_unit_field(
                piece,
                :slab_pieces,
                (),
            )
        return (;
            object_kind = :atom_growth_plan_unit_slab_set_descriptor,
            descriptor_kind = :outer_mismatch_boundary_slab_set,
            box = region.box,
            box_shape = region.box_shape,
            support_count = region.support_count,
            slab_piece_count = length(slab_pieces),
            slab_piece_roles = Tuple(
                _pqs_source_box_route_driver_atom_growth_unit_field(
                    slab_piece,
                    :role,
                ) for slab_piece in slab_pieces
            ),
            slab_piece_support_counts = Tuple(
                _pqs_source_box_route_driver_atom_growth_unit_field(
                    slab_piece,
                    :support_count,
                ) for slab_piece in slab_pieces
            ),
            final_column_ranges_available = false,
        )
    end
    return (;
        object_kind = :atom_growth_plan_unit_box_descriptor,
        descriptor_kind = :source_box,
        box = region.box,
        box_shape = region.box_shape,
        support_count = region.support_count,
        final_column_ranges_available = false,
    )
end

function _pqs_source_box_route_driver_atom_growth_plan_unit_record(
    region;
    unit_key::Symbol,
)
    source_cpbs = _pqs_source_box_route_driver_lowering_source_cpbs(region)
    lowering_recipe =
        _pqs_source_box_route_driver_lowering_recipe(region, source_cpbs)
    owned_support = _pqs_source_box_route_driver_owned_support(region)
    intermediate_retained_space =
        _pqs_source_box_route_driver_intermediate_retained_space(unit_key)
    shell_realization =
        _pqs_source_box_route_driver_shell_realization_contract(unit_key)
    final_retained_unit =
        _pqs_source_box_route_driver_final_retained_unit_contract(
            unit_key,
            region.role,
        )
    return (;
        object_kind = :cartesian_atom_growth_plan_unit,
        unit_key,
        unit_role = region.role,
        region_order_index = region.order_index,
        shellification_policy = :atom_growth_complete_rectangular,
        shellification_region_object_kind = region.object_kind,
        shellification_region_is_cpb = false,
        owned_support,
        source_descriptor =
            _pqs_source_box_route_driver_atom_growth_source_descriptor(region),
        source_cpbs,
        source_cpb_count = length(source_cpbs),
        source_box = region.box,
        source_dimensions = region.box_shape,
        source_dimension = region.support_count,
        support_count = region.support_count,
        lowering_stage = :coordinate_product_box_lowering,
        lowering_recipe,
        lowering_family = region.lowering_family,
        lowering_parameters =
            _pqs_source_box_route_driver_atom_growth_lowering_parameters(region),
        lowering_piece_object_kind = region.lowering_piece_object_kind,
        materialization_dependency = region.materialization_dependency,
        intermediate_retained_space,
        shell_realization,
        final_retained_unit,
        source_backed = region.source_backed,
        independently_lowerable = region.independently_lowerable,
        retirement_target = region.retirement_target,
        retained_count = nothing,
        retained_range = nothing,
        retained_dimension = nothing,
        retained_count_known = false,
        retained_range_known = false,
        retained_dimension_known = false,
        materialized_units_available = false,
        cpb_contract = (;
            object_kind = :cartesian_cpb_stage_contract,
            shellification_layer = :owned_support,
            lowering_layer = :source_cpbs,
            construction_layer = :intermediate_to_final_retained_unit,
            pair_planning_layer = :final_retained_unit_pairs,
            shellification_region_is_cpb = false,
            final_unit_downstream_of_lowering = true,
        ),
        provenance_label = :bond_aligned_diatomic_atom_growth_shellification_plan,
    )
end

function _pqs_source_box_route_driver_atom_growth_plan_unit_inventory(
    low_order_shellization,
)
    if !low_order_shellization.atom_growth_scaffold_available ||
       isnothing(low_order_shellization.atom_growth_scaffold)
        return (;
            object_kind = :cartesian_atom_growth_plan_unit_inventory,
            status = low_order_shellization.status,
            private_development_only = true,
            unit_inventory_source = :blocked_atom_growth_shellification_plan,
            plan_units = (),
            unit_count = 0,
            unit_keys = (),
            unit_roles = (),
            support_counts = (),
            materialization_dependencies = (),
            source_backed_region_count = 0,
            source_backed_unit_count = 0,
            plan_lowerable_unit_count = 0,
            retained_unit_dimensions_known = false,
            retained_unit_ranges_known = false,
            retained_dimension_known = false,
            retained_dimension = nothing,
            materialized_units_available = false,
            route_skeleton_authority = false,
            blocker =
                low_order_shellization.atom_growth_shellification_plan_status,
        )
    end

    scaffold = low_order_shellization.atom_growth_scaffold
    role_counts = Dict{Symbol,Int}()
    plan_units = NamedTuple[]
    for region in scaffold.regions
        role_index = get(role_counts, region.role, 0) + 1
        role_counts[region.role] = role_index
        push!(
            plan_units,
            _pqs_source_box_route_driver_atom_growth_plan_unit_record(
                region;
                unit_key =
                    _pqs_source_box_route_driver_atom_growth_unit_key(
                        region.role,
                        role_index,
                    ),
            ),
        )
    end

    plan_units = Tuple(plan_units)
    unit_keys = Tuple(unit.unit_key for unit in plan_units)
    support_counts =
        NamedTuple{unit_keys}(Tuple(unit.support_count for unit in plan_units))
    source_backed_unit_count = count(unit -> unit.source_backed, plan_units)
    source_cpb_count = sum(unit.source_cpb_count for unit in plan_units; init = 0)
    pqs_lowering_prototype_unit =
        findfirst(
            unit -> unit.unit_role == :regular_shared_molecular_shell,
            plan_units,
        )
    pqs_lowering_prototype =
        isnothing(pqs_lowering_prototype_unit) ?
        nothing :
        _pqs_source_box_route_driver_pqs_lowering_prototype(
            plan_units[pqs_lowering_prototype_unit],
        )
    return (;
        object_kind = :cartesian_atom_growth_plan_unit_inventory,
        status = :available_atom_growth_plan_unit_inventory,
        private_development_only = true,
        unit_inventory_source = :atom_growth_shellification_plan,
        cpb_contract_stage = :shellification_to_cpb_lowering_to_construction,
        shellification_regions_are_cpbs = false,
        owned_support_available = true,
        lowering_source_cpbs_available = true,
        source_cpb_count,
        pqs_lowering_prototype_available = !isnothing(pqs_lowering_prototype),
        pqs_lowering_prototype,
        pqs_lowering_prototype_unit_key =
            isnothing(pqs_lowering_prototype) ?
            nothing :
            pqs_lowering_prototype.unit_key,
        plan_units,
        unit_count = length(plan_units),
        unit_keys,
        unit_roles = Tuple(unit.unit_role for unit in plan_units),
        support_counts,
        materialization_dependencies =
            Tuple(unit.materialization_dependency for unit in plan_units),
        source_backed_region_count =
            scaffold.materialization_dependency_counts.source_backed_region_count,
        source_backed_unit_count,
        plan_lowerable_unit_count =
            scaffold.materialization_dependency_counts.plan_lowerable_region_count,
        retained_unit_dimensions_known = false,
        retained_unit_ranges_known = false,
        retained_dimension_known = false,
        retained_dimension = nothing,
        materialized_units_available = false,
        route_skeleton_authority = false,
        blocker = nothing,
    )
end

function _pqs_source_box_route_driver_atom_growth_transform_contract_name(unit)
    unit.lowering_family == :white_lindsey_atom_local_child_shellification &&
        return :atom_local_child_shellification_sequence
    unit.materialization_dependency == :plan_lowerable_direct_slab &&
        return :direct_identity_selector
    unit.lowering_family == :white_lindsey_adaptive_complete_shell &&
        return :adaptive_complete_shell_layer
    unit.lowering_family == :outer_mismatch_boundary_slab_set &&
        return :outer_mismatch_boundary_slab_set
    return nothing
end

function _pqs_source_box_route_driver_atom_growth_transform_contract_record(
    unit,
)
    transform_contract =
        _pqs_source_box_route_driver_atom_growth_transform_contract_name(unit)
    isnothing(transform_contract) && throw(
        ArgumentError(
            "atom-growth plan unit $(unit.unit_key) has no transform contract for lowering family $(unit.lowering_family)",
        ),
    )
    return (;
        object_kind = :cartesian_atom_growth_transform_contract,
        unit_key = unit.unit_key,
        unit_role = unit.unit_role,
        transform_contract,
        contract_source = :atom_growth_plan_unit_inventory,
        cpb_contract_stage = :construction_transform_contract,
        owned_support = unit.owned_support,
        source_cpbs = unit.source_cpbs,
        source_cpb_count = unit.source_cpb_count,
        lowering_recipe = unit.lowering_recipe,
        intermediate_retained_space = unit.intermediate_retained_space,
        shell_realization = unit.shell_realization,
        final_retained_unit = unit.final_retained_unit,
        final_unit_downstream_of_lowering =
            unit.final_retained_unit.downstream_of_lowering,
        lowering_family = unit.lowering_family,
        materialization_dependency = unit.materialization_dependency,
        source_backed = unit.source_backed,
        independently_lowerable = unit.independently_lowerable,
        coefficient_transform_materialized = false,
        coefficient_map_materialized = false,
        retained_count_known = unit.retained_count_known,
        retained_range_known = unit.retained_range_known,
        retained_dimension_known = unit.retained_dimension_known,
    )
end

function _pqs_source_box_route_driver_pqs_transform_prototype(
    lowering_prototype,
    transform_contract,
)
    return (;
        object_kind = :cartesian_pqs_transform_metadata_prototype,
        status = :metadata_only_planned,
        private_development_only = true,
        transform_stage = :construction_transform_contract,
        transform_plan = (
            :source_retained_modes,
            :shell_projection,
            :lowdin_cleanup,
            :final_retained_unit,
        ),
        source_lowering_prototype = lowering_prototype,
        source_lowering_prototype_unit_key = lowering_prototype.unit_key,
        unit_key = lowering_prototype.unit_key,
        unit_role = lowering_prototype.unit_role,
        source_cpb = lowering_prototype.source_cpb,
        source_cpb_support_count =
            lowering_prototype.source_cpb_support_count,
        source_cpb_support_count_source =
            lowering_prototype.source_cpb_support_count_source,
        owned_support = lowering_prototype.owned_support,
        owned_support_count = lowering_prototype.owned_support_count,
        owned_support_count_source =
            lowering_prototype.owned_support_count_source,
        intermediate_retained_space =
            lowering_prototype.intermediate_retained_space,
        shell_realization = lowering_prototype.shell_realization,
        final_retained_unit = transform_contract.final_retained_unit,
        coefficient_transform_materialized = false,
        coefficient_maps_materialized = false,
        numerical_transform_materialized = false,
        source_operator_blocks_materialized = false,
        operator_blocks_materialized = false,
        pair_operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        dense_parent_space_operators_are_algorithm = false,
    )
end

function _pqs_source_box_route_driver_atom_growth_transform_contract_inventory(
    plan_unit_inventory,
)
    if isnothing(plan_unit_inventory) ||
       plan_unit_inventory.status != :available_atom_growth_plan_unit_inventory
        return (;
            object_kind = :cartesian_atom_growth_transform_contract_inventory,
            status = :blocked_missing_atom_growth_plan_unit_inventory,
            private_development_only = true,
            transform_contract_source = :blocked_missing_plan_unit_inventory,
            transform_contracts = (),
            contract_count = 0,
            unit_keys = (),
            unit_roles = (),
            contract_names = (),
            source_backed_contract_count = 0,
            coefficient_transforms_materialized = false,
            coefficient_maps_materialized = false,
            pqs_transform_prototype_available = false,
            pqs_transform_prototype = nothing,
            source_lowering_prototype_unit_key = nothing,
            retained_unit_dimensions_known = false,
            retained_unit_ranges_known = false,
            retained_dimension_known = false,
            retained_dimension = nothing,
            route_skeleton_authority = false,
            blocker = :missing_atom_growth_plan_unit_inventory,
        )
    end

    unsupported_unit_keys = Tuple(
        unit.unit_key for unit in plan_unit_inventory.plan_units
        if isnothing(
            _pqs_source_box_route_driver_atom_growth_transform_contract_name(unit),
        )
    )
    if !isempty(unsupported_unit_keys)
        return (;
            object_kind = :cartesian_atom_growth_transform_contract_inventory,
            status = :blocked_unknown_atom_growth_transform_contract,
            private_development_only = true,
            transform_contract_source = :atom_growth_plan_unit_inventory,
            transform_contracts = (),
            contract_count = 0,
            unit_keys = (),
            unit_roles = (),
            contract_names = (),
            source_backed_contract_count = 0,
            coefficient_transforms_materialized = false,
            coefficient_maps_materialized = false,
            pqs_transform_prototype_available = false,
            pqs_transform_prototype = nothing,
            source_lowering_prototype_unit_key = nothing,
            retained_unit_dimensions_known = false,
            retained_unit_ranges_known = false,
            retained_dimension_known = false,
            retained_dimension = nothing,
            route_skeleton_authority = false,
            blocker = (;
                reason = :unknown_atom_growth_transform_contract,
                unit_keys = unsupported_unit_keys,
            ),
        )
    end

    transform_contracts = Tuple(
        _pqs_source_box_route_driver_atom_growth_transform_contract_record(unit)
        for unit in plan_unit_inventory.plan_units
    )
    pqs_lowering_prototype =
        hasproperty(plan_unit_inventory, :pqs_lowering_prototype_available) &&
        plan_unit_inventory.pqs_lowering_prototype_available ?
        plan_unit_inventory.pqs_lowering_prototype :
        nothing
    pqs_transform_contract_index =
        isnothing(pqs_lowering_prototype) ?
        nothing :
        findfirst(
            contract -> contract.unit_key == pqs_lowering_prototype.unit_key,
            transform_contracts,
        )
    pqs_transform_prototype =
        isnothing(pqs_transform_contract_index) ?
        nothing :
        _pqs_source_box_route_driver_pqs_transform_prototype(
            pqs_lowering_prototype,
            transform_contracts[pqs_transform_contract_index],
        )
    return (;
        object_kind = :cartesian_atom_growth_transform_contract_inventory,
        status = :available_atom_growth_transform_contract_inventory,
        private_development_only = true,
        transform_contract_source = :atom_growth_plan_unit_inventory,
        transform_contracts,
        contract_count = length(transform_contracts),
        unit_keys = Tuple(contract.unit_key for contract in transform_contracts),
        unit_roles = Tuple(contract.unit_role for contract in transform_contracts),
        contract_names =
            Tuple(contract.transform_contract for contract in transform_contracts),
        source_backed_contract_count =
            count(contract -> contract.source_backed, transform_contracts),
        coefficient_transforms_materialized = false,
        coefficient_maps_materialized = false,
        pqs_transform_prototype_available = !isnothing(pqs_transform_prototype),
        pqs_transform_prototype,
        source_lowering_prototype_unit_key =
            isnothing(pqs_transform_prototype) ?
            nothing :
            pqs_transform_prototype.source_lowering_prototype_unit_key,
        retained_unit_dimensions_known = false,
        retained_unit_ranges_known = false,
        retained_dimension_known = false,
        retained_dimension = nothing,
        route_skeleton_authority = false,
        blocker = nothing,
    )
end

function _pqs_source_box_route_driver_atom_growth_pair_record(
    left_unit,
    right_unit;
    left_unit_index::Int,
    right_unit_index::Int,
)
    return (;
        object_kind = :cartesian_atom_growth_plan_pair,
        pair_key = (left_unit.unit_key, right_unit.unit_key),
        left_unit_key = left_unit.unit_key,
        right_unit_key = right_unit.unit_key,
        left_unit_role = left_unit.unit_role,
        right_unit_role = right_unit.unit_role,
        left_unit_index,
        right_unit_index,
        pair_planning_source = :final_retained_units,
        left_final_retained_unit = left_unit.final_retained_unit,
        right_final_retained_unit = right_unit.final_retained_unit,
        left_owned_support = left_unit.owned_support,
        right_owned_support = right_unit.owned_support,
        pair_family = :white_lindsey_low_order_atom_growth_unit_pair,
        pair_contract = :planned_low_order_unit_pair_operator_block,
        pair_inventory_source = :atom_growth_unit_inventory,
        operator_pair_block_materialized = false,
        operator_block_materialized = false,
        retained_block_dimensions_known = false,
        operator_block_dimensions_known = false,
        retained_block_dimensions = nothing,
        operator_block_dimensions = nothing,
    )
end

function _pqs_source_box_route_driver_atom_growth_pair_inventory(
    plan_unit_inventory,
)
    if isnothing(plan_unit_inventory) ||
       plan_unit_inventory.status != :available_atom_growth_plan_unit_inventory
        return (;
            object_kind = :cartesian_atom_growth_plan_pair_inventory,
            status = :blocked_missing_atom_growth_plan_unit_inventory,
            private_development_only = true,
            pair_inventory_source = :blocked_missing_plan_unit_inventory,
            unit_count = 0,
            pair_count = 0,
            pair_entries = (),
            pair_family_counts =
                (white_lindsey_low_order_atom_growth_unit_pair = 0,),
            upper_triangular_unit_pairs = false,
            operator_pairs_materialized = false,
            pair_operator_blocks_materialized = false,
            operator_blocks_materialized = false,
            retained_block_dimensions_known = false,
            operator_block_dimensions_known = false,
            blocker = :missing_atom_growth_plan_unit_inventory,
        )
    end

    plan_units = plan_unit_inventory.plan_units
    pair_entries = NamedTuple[]
    for left_index in eachindex(plan_units)
        for right_index in left_index:length(plan_units)
            push!(
                pair_entries,
                _pqs_source_box_route_driver_atom_growth_pair_record(
                    plan_units[left_index],
                    plan_units[right_index];
                    left_unit_index = left_index,
                    right_unit_index = right_index,
                ),
            )
        end
    end

    pair_entries = Tuple(pair_entries)
    pair_count = length(pair_entries)
    return (;
        object_kind = :cartesian_atom_growth_plan_pair_inventory,
        status = :available_atom_growth_pair_inventory,
        private_development_only = true,
        pair_inventory_source = :atom_growth_unit_inventory,
        unit_count = length(plan_units),
        pair_count,
        pair_entries,
        pair_family_counts =
            (white_lindsey_low_order_atom_growth_unit_pair = pair_count,),
        upper_triangular_unit_pairs = true,
        operator_pairs_materialized = false,
        pair_operator_blocks_materialized = false,
        operator_blocks_materialized = false,
        retained_block_dimensions_known = false,
        operator_block_dimensions_known = false,
        blocker = nothing,
    )
end

function _pqs_source_box_route_driver_atom_growth_unit_record(;
    unit_key,
    unit_role,
    retained_unit_kind,
    source_family,
    source_box,
    retained_rule_kind,
    retained_rule_derivation,
    retained_range,
    provenance_label,
    source_dimensions = isnothing(source_box) ? nothing : Tuple(length.(source_box)),
    source_dimension =
        isnothing(source_dimensions) ?
        (isnothing(retained_range) ? nothing : length(retained_range)) :
        prod(source_dimensions),
)
    return _pqs_source_box_route_driver_unit_record(
        unit_key = unit_key,
        unit_role = unit_role,
        retained_unit_kind = retained_unit_kind,
        source_family = source_family,
        source_box = source_box,
        source_dimensions = source_dimensions,
        source_dimension = source_dimension,
        retained_rule_kind = retained_rule_kind,
        retained_rule_derivation = retained_rule_derivation,
        retained_range = retained_range,
        retained_count = isnothing(retained_range) ? nothing : length(retained_range),
        provenance_label = provenance_label,
        weight_semantics = :retained_basis_integral_weights,
    )
end

function _pqs_source_box_route_driver_diatomic_atom_growth_route_units(probe)
    if !probe.materialized || isnothing(probe.materialization) ||
       isnothing(probe.materialization.assembly)
        return (;
            object_kind = :white_lindsey_low_order_diatomic_atom_growth_route_units,
            route_family = :white_lindsey_low_order,
            status = :blocked_missing_atom_growth_assembly,
            private_development_only = true,
            retained_units = (),
            unit_inventory = nothing,
            standard_unit_inventory = nothing,
            retained_dimension = probe.retained_dimension,
            pair_entries = (),
            pair_family_counts = (white_lindsey_low_order_atom_growth = 0,),
            operator_pairs_materialized = false,
            weight_semantics = :retained_basis_integral_weights,
            blocker = :missing_atom_growth_assembly,
        )
    end

    scaffold = probe.scaffold
    assembly = probe.materialization.assembly
    retained_units = NamedTuple[]

    for (index, slab_set) in enumerate(scaffold.outer_mismatch_boundary_slab_sets)
        retained_range = assembly.outer_mismatch_column_ranges[index]
        push!(
            retained_units,
            _pqs_source_box_route_driver_atom_growth_unit_record(
                unit_key = _pqs_source_box_route_driver_atom_growth_unit_key(
                    :outer_mismatch_shared_molecular_shell,
                    index,
                ),
                unit_role = :outer_mismatch_shared_molecular_shell,
                retained_unit_kind =
                    :atom_growth_outer_mismatch_boundary_slab_set,
                source_family =
                    :white_lindsey_low_order_atom_growth_outer_mismatch,
                source_box = nothing,
                source_dimensions = nothing,
                source_dimension = slab_set.support_count,
                retained_rule_kind = :direct_boundary_slab_parent_sites,
                retained_rule_derivation =
                    :atom_growth_outer_mismatch_boundary_slab_set,
                retained_range = retained_range,
                provenance_label =
                    :bond_aligned_diatomic_atom_growth_outer_mismatch,
            ),
        )
    end

    push!(
        retained_units,
        _pqs_source_box_route_driver_atom_growth_unit_record(
            unit_key = :left_atom_box,
            unit_role = :left_atom_box,
            retained_unit_kind = :atom_growth_atom_local_child_box,
            source_family = :white_lindsey_low_order_atom_growth_atom_box,
            source_box = scaffold.left_child_plan.outer_box,
            retained_rule_kind = :atom_local_child_shellification_plan,
            retained_rule_derivation =
                :atom_growth_complete_rectangular_left_child_plan,
            retained_range = assembly.child_column_ranges[1],
            provenance_label = :bond_aligned_diatomic_atom_growth_left_atom_box,
        ),
    )

    if !isnothing(scaffold.contact_cap_region)
        push!(
            retained_units,
            _pqs_source_box_route_driver_atom_growth_unit_record(
                unit_key = :contact_cap,
                unit_role = :contact_cap,
                retained_unit_kind = :atom_growth_direct_contact_cap,
                source_family = :white_lindsey_low_order_atom_growth_contact_cap,
                source_box = scaffold.contact_cap_region.box,
                retained_rule_kind = :direct_contact_cap_parent_sites,
                retained_rule_derivation =
                    :atom_growth_complete_rectangular_contact_cap,
                retained_range = assembly.contact_cap_column_range,
                provenance_label =
                    :bond_aligned_diatomic_atom_growth_contact_cap,
            ),
        )
    end

    push!(
        retained_units,
        _pqs_source_box_route_driver_atom_growth_unit_record(
            unit_key = :right_atom_box,
            unit_role = :right_atom_box,
            retained_unit_kind = :atom_growth_atom_local_child_box,
            source_family = :white_lindsey_low_order_atom_growth_atom_box,
            source_box = scaffold.right_child_plan.outer_box,
            retained_rule_kind = :atom_local_child_shellification_plan,
            retained_rule_derivation =
                :atom_growth_complete_rectangular_right_child_plan,
            retained_range = assembly.child_column_ranges[2],
            provenance_label = :bond_aligned_diatomic_atom_growth_right_atom_box,
        ),
    )

    for (index, region) in enumerate(scaffold.shared_complete_shell_regions)
        retained_range = assembly.shared_shell_column_ranges[index]
        push!(
            retained_units,
            _pqs_source_box_route_driver_atom_growth_unit_record(
                unit_key = _pqs_source_box_route_driver_atom_growth_unit_key(
                    :regular_shared_molecular_shell,
                    index,
                ),
                unit_role = :regular_shared_molecular_shell,
                retained_unit_kind = :atom_growth_shared_complete_rectangular_shell,
                source_family =
                    :white_lindsey_low_order_atom_growth_shared_shell,
                source_box = region.outer_box,
                source_dimensions = nothing,
                source_dimension = region.support_count,
                retained_rule_kind = :shared_complete_rectangular_shell_plan,
                retained_rule_derivation =
                    :atom_growth_complete_rectangular_shared_shell_region,
                retained_range = retained_range,
                provenance_label =
                    :bond_aligned_diatomic_atom_growth_shared_shell,
            ),
        )
    end

    retained_units = Tuple(retained_units)
    unit_inventory = _pqs_source_box_route_driver_unit_inventory(retained_units)
    pair_entries = ()
    pair_family_counts = (white_lindsey_low_order_atom_growth = 0,)
    route_facts = (;
        source_dimensions = unit_inventory.source_dimensions,
        retained_units,
        retained_counts = unit_inventory.retained_counts,
        ranges = unit_inventory.ranges,
        retained_dimension = unit_inventory.retained_dimension,
        pair_entries,
        pair_family_counts,
    )

    return (;
        object_kind = :white_lindsey_low_order_diatomic_atom_growth_route_units,
        route_family = :white_lindsey_low_order,
        status = :available_atom_growth_retained_unit_inventory,
        private_development_only = true,
        retained_units,
        unit_inventory,
        standard_unit_inventory =
            _pqs_source_box_route_driver_standard_unit_inventory_summary(route_facts),
        retained_dimension = unit_inventory.retained_dimension,
        pair_entries,
        pair_family_counts,
        operator_pairs_materialized = false,
        pair_inventory_status =
            :assembled_sequence_payload_not_pair_decomposed,
        weight_semantics = :retained_basis_integral_weights,
        blocker = nothing,
    )
end

function _pqs_source_box_route_driver_diatomic_atom_growth_transform_inventory(
    probe,
)
    sequence = probe.materialization.sequence
    weights =
        isnothing(probe.basis_adapter) ?
        nothing :
        probe.basis_adapter.final_integral_weights
    retained_dimension = size(sequence.coefficient_matrix, 2)
    final_integral_weights_ready =
        !isnothing(weights) &&
        length(weights) == retained_dimension &&
        all(isfinite, weights) &&
        all(>(0.0), weights)

    return (;
        object_kind =
            :white_lindsey_low_order_diatomic_atom_growth_transform_inventory,
        route_family = :white_lindsey_low_order,
        status =
            final_integral_weights_ready ?
            :available_atom_growth_transform_inventory :
            :blocked_atom_growth_transform_inventory_contract,
        private_development_only = true,
        transform_source = :atom_growth_shell_sequence_coefficient_matrix,
        coefficient_matrix_size = size(sequence.coefficient_matrix),
        coefficient_matrix_finite = all(isfinite, sequence.coefficient_matrix),
        retained_dimension,
        support_count = length(sequence.support_indices),
        final_integral_weight_count = isnothing(weights) ? 0 : length(weights),
        final_integral_weights_status = probe.final_integral_weights_status,
        final_integral_weights_ready,
        weight_semantics = :retained_basis_integral_weights,
        blocker =
            final_integral_weights_ready ?
            nothing :
            :atom_growth_final_integral_weight_contract,
    )
end

function _pqs_source_box_route_driver_diatomic_atom_growth_operator_inventory(
    probe,
)
    fixed_block = probe.basis_adapter.fixed_block
    fixed_block_matrices =
        _white_lindsey_low_order_materialized_seed_operator_matrices(fixed_block)
    fixed_block_matrix_sizes =
        _white_lindsey_low_order_operator_matrix_sizes(fixed_block_matrices)
    fixed_block_finite_ready =
        _white_lindsey_low_order_operator_finite_ready(fixed_block_matrices)
    ham_adapter_available =
        !isnothing(probe.ham_adapter) &&
        probe.ham_adapter.status == :available_route_configured_diatomic_ham_adapter
    ham_matrices =
        ham_adapter_available ?
        (;
            overlap = probe.ham_adapter.operators.overlap,
            one_body_hamiltonian =
                probe.ham_adapter.operators.one_body_hamiltonian,
            interaction_matrix = probe.ham_adapter.operators.interaction_matrix,
        ) :
        nothing
    ham_matrix_sizes =
        ham_adapter_available ?
        _white_lindsey_low_order_operator_matrix_sizes(ham_matrices) :
        nothing
    ham_finite_ready =
        ham_adapter_available ?
        _white_lindsey_low_order_operator_finite_ready(ham_matrices) :
        nothing
    ham_all_finite =
        ham_adapter_available ? all(values(ham_finite_ready)) : false

    return (;
        object_kind =
            :white_lindsey_low_order_diatomic_atom_growth_operator_inventory,
        route_family = :white_lindsey_low_order,
        status =
            ham_adapter_available && ham_all_finite ?
            :available_atom_growth_operator_inventory :
            :blocked_atom_growth_operator_inventory_contract,
        private_development_only = true,
        operator_source = :atom_growth_fixed_block_and_ham_adapter,
        fixed_block_matrix_sizes,
        fixed_block_finite_ready,
        fixed_block_all_finite = all(values(fixed_block_finite_ready)),
        ham_adapter_status = probe.ham_adapter_status,
        ham_matrix_sizes,
        ham_finite_ready,
        ham_all_finite,
        final_integral_weights_status = probe.final_integral_weights_status,
        operator_pairs_materialized = false,
        pair_inventory_status =
            :assembled_sequence_payload_not_pair_decomposed,
        electron_electron_materialized = ham_adapter_available,
        overlap_materialized = ham_adapter_available,
        one_body_hamiltonian_materialized = ham_adapter_available,
        density_density_interaction_materialized = ham_adapter_available,
        blocker =
            ham_adapter_available && ham_all_finite ?
            nothing :
            :atom_growth_ham_operator_adapter_contract,
    )
end

function _pqs_source_box_route_driver_route_configured_diatomic_atom_growth_report(
    probe;
    basis_artifact_status,
    ham_artifact_status,
    basis_bundle_export_status,
    ham_bundle_export_status,
)
    if !probe.materialized || isnothing(probe.materialization)
        return (;
            object_kind =
                :white_lindsey_low_order_route_configured_diatomic_atom_growth_report,
            route_family = :white_lindsey_low_order,
            status = :blocked_missing_atom_growth_materialization,
            private_development_only = true,
            shellization_source =
                :bond_aligned_diatomic_atom_growth_construction_plan,
            shellization_authority =
                :bond_aligned_diatomic_atom_growth_construction_plan,
            active_source_authority = false,
            route_default_behavior_changed = false,
            sequence_available = false,
            retained_dimension = probe.retained_dimension,
            support_count = probe.support_count,
            route_units = nothing,
            transform_inventory = nothing,
            operator_inventory = nothing,
            final_integral_weights_status = probe.final_integral_weights_status,
            basis_adapter_status = probe.basis_adapter_status,
            ham_adapter_status = probe.ham_adapter_status,
            basis_artifact_status,
            ham_artifact_status,
            basis_bundle_export_status,
            ham_bundle_export_status,
            blocker = :missing_atom_growth_materialization,
        )
    elseif isnothing(probe.basis_adapter) ||
           isnothing(probe.basis_adapter.fixed_block)
        return (;
            object_kind =
                :white_lindsey_low_order_route_configured_diatomic_atom_growth_report,
            route_family = :white_lindsey_low_order,
            status = :blocked_missing_atom_growth_basis_adapter,
            private_development_only = true,
            shellization_source =
                :bond_aligned_diatomic_atom_growth_construction_plan,
            shellization_authority =
                :bond_aligned_diatomic_atom_growth_construction_plan,
            active_source_authority = false,
            route_default_behavior_changed = false,
            sequence_available = probe.sequence_available,
            retained_dimension = probe.retained_dimension,
            support_count = probe.support_count,
            route_units = nothing,
            transform_inventory = nothing,
            operator_inventory = nothing,
            final_integral_weights_status = probe.final_integral_weights_status,
            basis_adapter_status = probe.basis_adapter_status,
            ham_adapter_status = probe.ham_adapter_status,
            basis_artifact_status,
            ham_artifact_status,
            basis_bundle_export_status,
            ham_bundle_export_status,
            blocker = :atom_growth_basis_representation_contract,
        )
    end

    route_units =
        _pqs_source_box_route_driver_diatomic_atom_growth_route_units(probe)
    transform_inventory =
        _pqs_source_box_route_driver_diatomic_atom_growth_transform_inventory(probe)
    operator_inventory =
        _pqs_source_box_route_driver_diatomic_atom_growth_operator_inventory(probe)
    retained_dimension = probe.retained_dimension
    inventory_ready =
        route_units.status == :available_atom_growth_retained_unit_inventory &&
        transform_inventory.status == :available_atom_growth_transform_inventory &&
        operator_inventory.status == :available_atom_growth_operator_inventory

    return (;
        object_kind =
            :white_lindsey_low_order_route_configured_diatomic_atom_growth_report,
        route_family = :white_lindsey_low_order,
        status =
            inventory_ready ?
            :private_development_route_configured_atom_growth :
            :blocked_atom_growth_route_report_contract,
        private_development_only = true,
        materialization_status = probe.materialization.status,
        shellization_source =
            :bond_aligned_diatomic_atom_growth_construction_plan,
        shellization_authority =
            :bond_aligned_diatomic_atom_growth_construction_plan,
        active_source_authority = false,
        route_default_behavior_changed = false,
        sequence_available = probe.sequence_available,
        retained_dimension,
        support_count = probe.support_count,
        route_units,
        transform_inventory,
        operator_inventory,
        basis_adapter_summary = probe.basis_adapter_summary,
        ham_adapter_summary = probe.ham_adapter_summary,
        final_integral_weights_status = probe.final_integral_weights_status,
        final_integral_weights_ready =
            transform_inventory.final_integral_weights_ready,
        basis_adapter_status = probe.basis_adapter_status,
        ham_adapter_status = probe.ham_adapter_status,
        basis_artifact_status,
        ham_artifact_status,
        basis_bundle_export_status,
        ham_bundle_export_status,
        operator_pairs_materialized = route_units.operator_pairs_materialized,
        electron_electron_materialized = operator_inventory.electron_electron_materialized,
        weight_semantics = :retained_basis_integral_weights,
        blocker = inventory_ready ? nothing : :atom_growth_route_report_contract,
    )
end

function _pqs_source_box_route_driver_diatomic_atom_growth_materialization(
    context,
)
    (;
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
    ) = context

    atom_growth_probe =
        route_configured_diatomic_atom_growth_materializer_probe
    atom_growth_basis_adapter = atom_growth_probe.basis_adapter
    atom_growth_ham_adapter = atom_growth_probe.ham_adapter
    atom_growth_basis_adapter_available =
        !isnothing(atom_growth_basis_adapter) &&
        atom_growth_basis_adapter.status ==
        :available_route_configured_diatomic_atom_growth_basis_adapter
    atom_growth_ham_adapter_available =
        !isnothing(atom_growth_ham_adapter) &&
        atom_growth_ham_adapter.status ==
        :available_route_configured_diatomic_ham_adapter
    atom_growth_materialized =
        route_configured_diatomic_atom_growth_materializer_probe_consumed
    atom_growth_export_blocker =
        !atom_growth_materialized ?
        something(
            route_configured_diatomic_atom_growth_materializer_probe_blocker,
            :atom_growth_materializer_not_consumed,
        ) :
        !atom_growth_basis_adapter_available ?
        something(
            route_configured_diatomic_atom_growth_basis_adapter_blocker,
            :atom_growth_basis_representation_contract,
        ) :
        nothing
    atom_growth_ham_export_blocker =
        !save_ham_artifact ?
        nothing :
        !atom_growth_basis_adapter_available ?
        atom_growth_export_blocker :
        !atom_growth_ham_adapter_available ?
        something(
            route_configured_diatomic_atom_growth_ham_adapter_blocker,
            :atom_growth_ham_operator_adapter_contract,
        ) :
        nothing
    atom_growth_artifact_export_requested =
        save_basis_artifact || save_ham_artifact
    atom_growth_retained_dimension =
        atom_growth_basis_adapter_available ?
        atom_growth_basis_adapter.retained_dimension :
        atom_growth_probe.retained_dimension
    atom_growth_basis_artifact_status =
        save_basis_artifact ?
        (
            atom_growth_basis_adapter_available ?
            :written_route_configured_diatomic_atom_growth_basis_only_bundle :
            :not_written_route_configured_diatomic_atom_growth_basis_adapter_blocked
        ) :
        :not_requested
    atom_growth_ham_artifact_status =
        save_ham_artifact ?
        (
            atom_growth_ham_adapter_available ?
            :written_route_configured_diatomic_atom_growth_ham_bundle :
            :not_written_route_configured_diatomic_atom_growth_ham_adapter_blocked
        ) :
        :not_requested
    atom_growth_basis_bundle_export_status =
        atom_growth_artifact_export_requested ?
        (
            atom_growth_basis_adapter_available ?
            :supported_route_configured_diatomic_atom_growth_basis_only_fixed_block :
            :pending_route_configured_diatomic_atom_growth_basis_export
        ) :
        :not_requested
    atom_growth_ham_bundle_export_status =
        atom_growth_artifact_export_requested && atom_growth_ham_adapter_available ?
        :available_route_configured_diatomic_atom_growth_ham_bundle_payload :
        save_ham_artifact ?
        something(
            atom_growth_ham_export_blocker,
            :pending_route_configured_diatomic_atom_growth_ham_export,
        ) :
        :not_requested
    atom_growth_ham_preflight_status =
        atom_growth_artifact_export_requested && atom_growth_ham_adapter_available ?
        :available_route_configured_diatomic_atom_growth_ham_adapter :
        save_ham_artifact ?
        route_configured_diatomic_atom_growth_ham_adapter_status :
        :not_requested
    atom_growth_ham_interaction_treatment_consumed =
        atom_growth_artifact_export_requested && atom_growth_ham_adapter_available ?
        atom_growth_ham_adapter.interaction_treatment :
        nothing
    atom_growth_ham_interaction_treatment_status =
        atom_growth_artifact_export_requested && atom_growth_ham_adapter_available ?
        :available_route_configured_diatomic_ham_interaction_treatment :
        save_ham_artifact ?
        route_configured_diatomic_atom_growth_ham_adapter_status :
        :not_requested
    atom_growth_materialized_report =
        atom_growth_materialized ?
        _pqs_source_box_route_driver_route_configured_diatomic_atom_growth_report(
            atom_growth_probe;
            basis_artifact_status = atom_growth_basis_artifact_status,
            ham_artifact_status = atom_growth_ham_artifact_status,
            basis_bundle_export_status =
                atom_growth_basis_bundle_export_status,
            ham_bundle_export_status =
                atom_growth_ham_bundle_export_status,
        ) :
        nothing
    atom_growth_materialized_report_kind =
        isnothing(atom_growth_materialized_report) ?
        nothing :
        atom_growth_materialized_report.object_kind
    atom_growth_artifact_meta = (;
        route_family,
        route_kind = report.recipe_metadata.route_kind,
        benchmark_role = report.recipe_metadata.benchmark_role,
        materialized_report_kind =
            something(
                atom_growth_materialized_report_kind,
                :not_materialized_atom_growth_probe,
            ),
        shellification_materialization_kind =
            atom_growth_materialized ?
            atom_growth_probe.materialization.object_kind :
            :not_materialized_atom_growth_probe,
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
        shellization_summary_available = false,
        shellization_source =
            :bond_aligned_diatomic_atom_growth_construction_plan,
        shellization_authority =
            :bond_aligned_diatomic_atom_growth_construction_plan,
        active_source_authority = false,
        route_configured_shellization_consumed = atom_growth_materialized,
        route_configured_diatomic_atom_growth_probe_consumed =
            route_configured_diatomic_atom_growth_materializer_probe_consumed,
        route_configured_legacy_diatomic_source_consumed = false,
        route_default_behavior_changed = false,
        materialized_shellization_stage =
            :atom_growth_complete_rectangular_low_order,
        seed_materialization_status =
            :not_seed_route_configured_diatomic_atom_growth_shellization,
        private_development_only = true,
    )
    atom_growth_basis_artifact_written = false
    if save_basis_artifact && atom_growth_basis_adapter_available
        write_cartesian_basis_bundle_jld2(
            basisfile,
            atom_growth_basis_adapter.fixed_block;
            include_ham = false,
            meta = (;
                atom_growth_artifact_meta...,
                export_status = :basis_only,
                basis_export_status =
                    atom_growth_basis_bundle_export_status,
                ham_export_status =
                    :artifact_local_basis_only_no_ham_payload,
                ham_export_blocker = nothing,
                companion_ham_artifact_requested = save_ham_artifact,
                companion_ham_artifact_status =
                    atom_growth_ham_artifact_status,
                companion_ham_export_status =
                    atom_growth_ham_bundle_export_status,
                companion_ham_export_blocker =
                    atom_growth_ham_export_blocker,
            ),
        )
        atom_growth_basis_artifact_written = true
    end
    atom_growth_ham_artifact_written = false
    if save_ham_artifact && atom_growth_ham_adapter_available
        write_cartesian_basis_bundle_jld2(
            hamfile,
            atom_growth_ham_adapter.operators;
            include_ham = true,
            meta = (;
                atom_growth_artifact_meta...,
                export_status = :basis_and_ham,
                basis_export_status =
                    atom_growth_basis_bundle_export_status,
                ham_preflight_status = atom_growth_ham_preflight_status,
                ham_operator_payload_status =
                    atom_growth_ham_adapter.operator_payload_status,
                ham_interaction_status =
                    atom_growth_ham_adapter.interaction_status,
                ham_export_status =
                    atom_growth_ham_bundle_export_status,
                ham_export_blocker = nothing,
            ),
        )
        atom_growth_ham_artifact_written = true
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
            atom_growth_ham_interaction_treatment_consumed,
        route_configured_diatomic_ham_interaction_treatment_status =
            atom_growth_ham_interaction_treatment_status,
        status =
            atom_growth_export_blocker === nothing &&
            (!save_ham_artifact || atom_growth_ham_export_blocker === nothing) ?
            (
                atom_growth_artifact_export_requested ?
                :materialized_route_configured_diatomic_atom_growth_artifacts_available :
                :materialized_route_configured_diatomic_atom_growth_report_available
            ) :
            (
                atom_growth_artifact_export_requested ?
                :blocked_route_configured_diatomic_atom_growth_artifact_export :
                :blocked_route_configured_diatomic_atom_growth_materialization_report
            ),
        materialized_report = atom_growth_materialized_report,
        materialized_report_kind = atom_growth_materialized_report_kind,
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
            atom_growth_probe.basis_adapter_summary,
        route_configured_diatomic_ham_adapter_summary =
            atom_growth_probe.ham_adapter_summary,
        shellization_summary = nothing,
        shellization_summary_available = false,
        shellization_source =
            :bond_aligned_diatomic_atom_growth_construction_plan,
        route_configured_shellization_consumed = atom_growth_materialized,
        route_configured_legacy_diatomic_source_consumed = false,
        materialized_shellization_stage =
            :atom_growth_complete_rectangular_low_order,
        seed_materialization_status =
            :not_seed_route_configured_diatomic_atom_growth_shellization,
        retained_dimension = atom_growth_retained_dimension,
        final_integral_weights_status =
            route_configured_diatomic_atom_growth_final_integral_weights_status,
        one_body_operator_status =
            atom_growth_ham_adapter_available ?
            :available_route_configured_diatomic_operator_payload :
            :not_requested,
        basis_bundle_export_status =
            atom_growth_basis_bundle_export_status,
        basis_artifact_status = atom_growth_basis_artifact_status,
        basis_artifact_written = atom_growth_basis_artifact_written,
        basisfile,
        basis_artifact_path =
            atom_growth_basis_artifact_written ? basisfile : nothing,
        basis_export_blocker =
            atom_growth_basis_adapter_available ?
            nothing :
            atom_growth_export_blocker,
        ham_preflight_status = atom_growth_ham_preflight_status,
        ham_missing_builder = atom_growth_ham_export_blocker,
        ham_operator_payload_status =
            atom_growth_ham_adapter_available ?
            atom_growth_ham_adapter.operator_payload_status :
            route_configured_diatomic_atom_growth_ham_adapter_status,
        ham_interaction_status =
            atom_growth_ham_adapter_available ?
            atom_growth_ham_adapter.interaction_status :
            route_configured_diatomic_atom_growth_ham_adapter_status,
        ham_bundle_export_status =
            atom_growth_ham_bundle_export_status,
        ham_artifact_status = atom_growth_ham_artifact_status,
        ham_artifact_written = atom_growth_ham_artifact_written,
        hamfile,
        ham_export_blocker = atom_growth_ham_export_blocker,
        ham_preflight = nothing,
        pqs_materialization_status = :not_applicable,
    )
end

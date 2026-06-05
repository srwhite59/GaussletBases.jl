# Convert existing private staged route records into CartesianRouteCore
# sidecars. Behavior-neutral adapter layer only.

function _cartesian_route_core_staged_field(record, field::Symbol, default = nothing)
    return hasproperty(record, field) ? getproperty(record, field) : default
end

function _cartesian_route_core_named_metadata(value)
    isnothing(value) && return (;)
    value isa NamedTuple && return value
    return (; staged_metadata = value)
end

function _cartesian_route_core_merge_metadata(record; extra = (;))
    return merge(
        _cartesian_route_core_named_metadata(
            _cartesian_route_core_staged_field(record, :metadata, nothing),
        ),
        NamedTuple(extra),
    )
end

function _cartesian_route_core_optional_positive_int(value, label::AbstractString)
    isnothing(value) && return nothing
    value isa Integer && value > 0 ||
        throw(ArgumentError("$label must be a positive integer or nothing"))
    return Int(value)
end

function _cartesian_route_core_is_complete_shell(
    outer_box::CartesianRouteCore.CoordinateProductBox,
    inner_box::CartesianRouteCore.CoordinateProductBox,
)
    CartesianRouteCore.codimension(outer_box) == 0 || return false
    CartesianRouteCore.codimension(inner_box) == 0 || return false
    for axis_index in 1:3
        outer = CartesianRouteCore.intervals(outer_box)[axis_index]
        inner = CartesianRouteCore.intervals(inner_box)[axis_index]
        length(outer) >= 3 || return false
        first(inner) == first(outer) + 1 || return false
        last(inner) == last(outer) - 1 || return false
    end
    return true
end

function _cartesian_route_core_cpb_from_staged(cpb_record)
    cpb_record isa CartesianRouteCore.CoordinateProductBox && return cpb_record
    intervals =
        _cartesian_route_core_staged_field(
            cpb_record,
            :intervals,
            _cartesian_route_core_staged_field(cpb_record, :box, nothing),
        )
    isnothing(intervals) &&
        throw(ArgumentError("staged CPB record is missing intervals/box"))
    role =
        _cartesian_route_core_staged_field(
            cpb_record,
            :role,
            :staged_coordinate_product_box,
        )
    metadata = _cartesian_route_core_merge_metadata(
        cpb_record;
        extra = (;
            staged_object_kind =
                _cartesian_route_core_staged_field(cpb_record, :object_kind, nothing),
            cpb_family =
                _cartesian_route_core_staged_field(cpb_record, :cpb_family, nothing),
            source =
                _cartesian_route_core_staged_field(cpb_record, :source, nothing),
            staged_support_count =
                _cartesian_route_core_staged_field(cpb_record, :support_count, nothing),
        ),
    )
    box = CartesianRouteCore.cpb(intervals; role, metadata)
    staged_support_count =
        _cartesian_route_core_staged_field(cpb_record, :support_count, nothing)
    if !isnothing(staged_support_count) &&
       CartesianRouteCore.support_count(box) != staged_support_count
        throw(ArgumentError("staged CPB support count does not match intervals"))
    end
    return box
end

function _cartesian_route_core_filled_cpb_from_staged(cpb_record)
    box = _cartesian_route_core_cpb_from_staged(cpb_record)
    CartesianRouteCore.codimension(box) == 0 ||
        throw(ArgumentError("staged filled source CPB must have codimension 0"))
    return box
end

function _cartesian_route_core_owned_support_from_staged(owned_support_record)
    outer_record =
        _cartesian_route_core_staged_field(owned_support_record, :outer_cpb, nothing)
    isnothing(outer_record) &&
        throw(ArgumentError("staged owned support is missing outer_cpb"))
    inner_record =
        _cartesian_route_core_staged_field(
            owned_support_record,
            :inner_exclusion_cpb,
            nothing,
        )
    outer_box =
        isnothing(inner_record) ?
        _cartesian_route_core_cpb_from_staged(outer_record) :
        _cartesian_route_core_filled_cpb_from_staged(outer_record)
    support_kind =
        _cartesian_route_core_staged_field(
            owned_support_record,
            :support_kind,
            :staged_owned_support,
        )
    metadata = _cartesian_route_core_merge_metadata(
        owned_support_record;
        extra = (;
            staged_object_kind =
                _cartesian_route_core_staged_field(
                    owned_support_record,
                    :object_kind,
                    nothing,
                ),
            owned_support_authority =
                _cartesian_route_core_staged_field(
                    owned_support_record,
                    :owned_support_authority,
                    nothing,
                ),
            shellification_authority_scope =
                _cartesian_route_core_staged_field(
                    owned_support_record,
                    :shellification_authority_scope,
                    nothing,
                ),
            staged_support_count =
                _cartesian_route_core_staged_field(
                    owned_support_record,
                    :support_count,
                    nothing,
                ),
        ),
    )

    support =
        isnothing(inner_record) ?
        CartesianRouteCore.owned_cpb(outer_box; support_kind, metadata) :
        begin
            inner_box = _cartesian_route_core_filled_cpb_from_staged(inner_record)
            _cartesian_route_core_is_complete_shell(outer_box, inner_box) ?
            CartesianRouteCore.complete_shell_support(
                outer_box,
                inner_box;
                support_kind,
                metadata,
            ) :
            CartesianRouteCore.OwnedSupport(
                support_kind,
                (),
                outer_box,
                inner_box,
                metadata,
            )
        end

    staged_support_count =
        _cartesian_route_core_staged_field(owned_support_record, :support_count, nothing)
    if !isnothing(staged_support_count) &&
       CartesianRouteCore.support_count(support) != staged_support_count
        throw(ArgumentError("staged owned support count does not match CRC support"))
    end
    return support
end

function _cartesian_route_core_shellification_region_from_staged(record)
    owned_support_record =
        _cartesian_route_core_staged_field(record, :owned_support, nothing)
    isnothing(owned_support_record) &&
        throw(ArgumentError("staged record is missing owned_support"))
    unit_role =
        _cartesian_route_core_staged_field(
            record,
            :unit_role,
            _cartesian_route_core_staged_field(record, :role, :staged_unit),
        )
    unit_key = _cartesian_route_core_staged_field(record, :unit_key, nothing)
    metadata = _cartesian_route_core_merge_metadata(
        record;
        extra = (;
            staged_object_kind =
                _cartesian_route_core_staged_field(record, :object_kind, nothing),
            unit_key,
            shellification_policy =
                _cartesian_route_core_staged_field(
                    record,
                    :shellification_policy,
                    nothing,
                ),
            source = :staged_route_record,
        ),
    )
    return CartesianRouteCore.shellification_region(
        unit_role,
        _cartesian_route_core_owned_support_from_staged(owned_support_record);
        metadata,
    )
end

function _cartesian_route_core_sidecar_object(;
    sidecar_source::Symbol,
    shellification_region,
    lowering_source,
    intermediate_retained_space,
    shell_realization,
    final_retained_unit,
)
    return (;
        object_kind = :cartesian_route_core_sidecar,
        status = :metadata_only_adapter_sidecar,
        private_development_only = true,
        sidecar_source,
        shellification_region,
        lowering_source,
        intermediate_retained_space,
        shell_realization,
        final_retained_unit,
        numerical_behavior_changed = false,
        materialization_behavior_changed = false,
    )
end

function _cartesian_route_core_lw_complete_shell_sidecar(unit)
    _cartesian_route_core_staged_field(unit, :lowering_family, nothing) ===
    :white_lindsey_adaptive_complete_shell ||
        throw(ArgumentError("LW CRC sidecar expects a complete-shell staged unit"))
    source_cpbs = Tuple(
        _cartesian_route_core_cpb_from_staged(cpb)
        for cpb in _cartesian_route_core_staged_field(unit, :source_cpbs, ())
    )
    length(source_cpbs) == 26 ||
        throw(ArgumentError("LW complete-shell CRC sidecar expects 26 source CPBs"))

    family_counts = (;
        facet_cpb = count(
            cpb -> cpb.metadata.cpb_family === :facet_cpb,
            source_cpbs,
        ),
        edge_cpb = count(cpb -> cpb.metadata.cpb_family === :edge_cpb, source_cpbs),
        corner_cpb = count(
            cpb -> cpb.metadata.cpb_family === :corner_cpb,
            source_cpbs,
        ),
    )
    family_counts == (facet_cpb = 6, edge_cpb = 12, corner_cpb = 8) ||
        throw(ArgumentError("LW complete-shell CRC sidecar expects 6/12/8 CPB strata"))

    shellification_region =
        _cartesian_route_core_shellification_region_from_staged(unit)
    stratum_count = sum(CartesianRouteCore.support_count, source_cpbs; init = 0)
    stratum_count == CartesianRouteCore.support_count(shellification_region) ||
        throw(ArgumentError("LW source CPB strata do not match owned shell support"))

    lowering_source = CartesianRouteCore.white_lindsey_boundary_strata_lowering(
        shellification_region,
        source_cpbs;
        metadata = (;
            source = :atom_growth_plan_unit,
            unit_key = unit.unit_key,
            staged_lowering_family = unit.lowering_family,
            source_cpb_family_counts = family_counts,
            coefficient_maps_materialized = false,
            operator_blocks_materialized = false,
            pair_operator_blocks_materialized = false,
            hamiltonian_data_materialized = false,
        ),
    )
    dimension =
        _cartesian_route_core_optional_positive_int(
            _cartesian_route_core_staged_field(unit, :retained_dimension, nothing),
            "retained_dimension",
        )
    intermediate_retained_space = CartesianRouteCore.intermediate_retained_space(
        lowering_source;
        retained_rule = :white_lindsey_boundary_stratum_product,
        dimension,
        materialized = false,
        metadata = (;
            source = :atom_growth_plan_unit,
            unit_key = unit.unit_key,
            status = :deferred_not_materialized,
        ),
    )
    shell_realization = CartesianRouteCore.trivial_shell_realization(
        intermediate_retained_space,
        shellification_region;
        status = :deferred_or_trivial_for_white_lindsey_lowering,
        final_dimension = dimension,
        metadata = (;
            source = :atom_growth_plan_unit,
            unit_key = unit.unit_key,
            shell_row_oracle_authority = false,
        ),
    )
    final_retained_unit = CartesianRouteCore.final_retained_unit(
        unit.unit_key,
        unit.unit_role,
        lowering_source,
        intermediate_retained_space,
        shell_realization;
        column_range = _cartesian_route_core_staged_field(unit, :retained_range, nothing),
        dimension,
        metadata = (;
            source = :atom_growth_plan_unit,
            status = :deferred_until_materialization,
            pair_planning_input = true,
        ),
    )

    return _cartesian_route_core_sidecar_object(
        ;
        sidecar_source = :atom_growth_lw_complete_shell_plan_unit,
        shellification_region,
        lowering_source,
        intermediate_retained_space,
        shell_realization,
        final_retained_unit,
    )
end

function _cartesian_route_core_ranges_intersect(
    left::UnitRange{Int},
    right::UnitRange{Int},
)
    return max(first(left), first(right)) <= min(last(left), last(right))
end

function _cartesian_route_core_cpbs_overlap(left, right)
    return all(
        axis_index -> _cartesian_route_core_ranges_intersect(
            CartesianRouteCore.intervals(left)[axis_index],
            CartesianRouteCore.intervals(right)[axis_index],
        ),
        1:3,
    )
end

function _cartesian_route_core_cpbs_are_disjoint(cpbs)
    for left_index in eachindex(cpbs)
        for right_index in (left_index + 1):length(cpbs)
            _cartesian_route_core_cpbs_overlap(
                cpbs[left_index],
                cpbs[right_index],
            ) && return false
        end
    end
    return true
end

function _cartesian_route_core_boundary_slab_set_sidecar(unit)
    _cartesian_route_core_staged_field(unit, :lowering_family, nothing) ===
    :outer_mismatch_boundary_slab_set ||
        throw(ArgumentError("boundary slab-set CRC sidecar expects outer mismatch unit"))
    source_cpbs = Tuple(
        _cartesian_route_core_cpb_from_staged(cpb)
        for cpb in _cartesian_route_core_staged_field(unit, :source_cpbs, ())
    )
    length(source_cpbs) >= 2 ||
        throw(ArgumentError("boundary slab-set CRC sidecar expects slab-piece source CPBs"))
    _cartesian_route_core_cpbs_are_disjoint(source_cpbs) ||
        throw(ArgumentError("boundary slab-set source CPBs must be disjoint"))

    shellification_region =
        _cartesian_route_core_shellification_region_from_staged(unit)
    owned_support = CartesianRouteCore.owned_support(shellification_region)
    isempty(owned_support.cpbs) ||
        throw(ArgumentError("outer mismatch support must not be a single owned CPB"))
    !isnothing(owned_support.outer_box) && !isnothing(owned_support.inner_exclusion_box) ||
        throw(ArgumentError("outer mismatch support must be a shell-difference support"))

    source_support_count =
        sum(CartesianRouteCore.support_count, source_cpbs; init = 0)
    source_support_count == CartesianRouteCore.support_count(shellification_region) ||
        throw(ArgumentError("boundary slab-set source CPBs do not match owned support"))

    lowering_source = CartesianRouteCore.lowering_source(
        :direct_boundary_slab_set,
        shellification_region,
        source_cpbs;
        metadata = (;
            source = :atom_growth_plan_unit,
            unit_key = unit.unit_key,
            staged_lowering_family = unit.lowering_family,
            source_cpb_count = length(source_cpbs),
            source_support_count,
            coefficient_maps_materialized = false,
            operator_blocks_materialized = false,
            pair_operator_blocks_materialized = false,
            hamiltonian_data_materialized = false,
        ),
    )
    intermediate_retained_space = CartesianRouteCore.intermediate_retained_space(
        lowering_source;
        retained_rule = :direct_slab_set_identity_modes,
        dimension = source_support_count,
        materialized = false,
        metadata = (;
            source = :atom_growth_plan_unit,
            unit_key = unit.unit_key,
            status = :deferred_not_materialized,
        ),
    )
    shell_realization = CartesianRouteCore.trivial_shell_realization(
        intermediate_retained_space,
        shellification_region;
        status = :direct_or_trivial,
        final_dimension = source_support_count,
        metadata = (;
            source = :atom_growth_plan_unit,
            unit_key = unit.unit_key,
            shell_row_oracle_authority = false,
        ),
    )
    final_retained_unit = CartesianRouteCore.final_retained_unit(
        unit.unit_key,
        unit.unit_role,
        lowering_source,
        intermediate_retained_space,
        shell_realization;
        dimension = source_support_count,
        metadata = (;
            source = :atom_growth_plan_unit,
            status = :deferred_until_materialization,
            pair_planning_input = true,
        ),
    )

    return _cartesian_route_core_sidecar_object(
        ;
        sidecar_source = :atom_growth_outer_mismatch_boundary_slab_set_plan_unit,
        shellification_region,
        lowering_source,
        intermediate_retained_space,
        shell_realization,
        final_retained_unit,
    )
end

function _cartesian_route_core_atom_side_from_unit(unit)
    unit_key = _cartesian_route_core_staged_field(unit, :unit_key, nothing)
    unit_key === :left_atom_box && return :left
    unit_key === :right_atom_box && return :right
    throw(ArgumentError("atom-local child CRC sidecar requires left/right atom box unit_key"))
end

function _cartesian_route_core_atom_local_child_source_cpb(unit)
    source_box = _cartesian_route_core_staged_field(unit, :source_box, nothing)
    isnothing(source_box) &&
        throw(ArgumentError("atom-local child CRC sidecar is missing source_box"))
    descriptor =
        _cartesian_route_core_staged_field(unit, :source_descriptor, nothing)
    isnothing(descriptor) &&
        throw(ArgumentError("atom-local child CRC sidecar is missing source_descriptor"))
    _cartesian_route_core_staged_field(descriptor, :descriptor_kind, nothing) ===
    :source_box ||
        throw(ArgumentError("atom-local child CRC sidecar expects a source-box descriptor"))

    unit_key = _cartesian_route_core_staged_field(unit, :unit_key, nothing)
    descriptor_support_count =
        _cartesian_route_core_staged_field(descriptor, :support_count, nothing)
    metadata = _cartesian_route_core_merge_metadata(
        descriptor;
        extra = (;
            staged_object_kind =
                _cartesian_route_core_staged_field(descriptor, :object_kind, nothing),
            descriptor_kind =
                _cartesian_route_core_staged_field(descriptor, :descriptor_kind, nothing),
            source = :atom_growth_atom_local_child_source_box,
            unit_key,
            staged_support_count = descriptor_support_count,
            child_source_cpbs_enumerated = false,
        ),
    )
    source_cpb = CartesianRouteCore.cpb(
        source_box;
        role = Symbol(String(unit_key), "_source_cpb"),
        metadata,
    )
    if !isnothing(descriptor_support_count) &&
       CartesianRouteCore.support_count(source_cpb) != descriptor_support_count
        throw(ArgumentError("atom-local child source-box descriptor support count does not match source_box"))
    end
    return source_cpb
end

function _cartesian_route_core_atom_local_child_shellification_sidecar(unit)
    _cartesian_route_core_staged_field(unit, :lowering_family, nothing) ===
    :white_lindsey_atom_local_child_shellification ||
        throw(ArgumentError("atom-local child CRC sidecar expects atom-local child shellification"))
    _cartesian_route_core_staged_field(unit, :materialization_dependency, nothing) ===
    :plan_lowerable_atom_local_child_shellification ||
        throw(ArgumentError("atom-local child CRC sidecar expects the atom-local child materialization dependency"))

    atom_side = _cartesian_route_core_atom_side_from_unit(unit)
    shellification_region =
        _cartesian_route_core_shellification_region_from_staged(unit)
    owned_support = CartesianRouteCore.owned_support(shellification_region)
    length(owned_support.cpbs) == 1 ||
        throw(ArgumentError("atom-local child CRC sidecar expects coordinate-product owned support"))
    isnothing(owned_support.inner_exclusion_box) ||
        throw(ArgumentError("atom-local child CRC sidecar does not expect a shell-difference owned support"))

    source_cpb = _cartesian_route_core_atom_local_child_source_cpb(unit)
    owned_cpb = only(owned_support.cpbs)
    CartesianRouteCore.intervals(source_cpb) == CartesianRouteCore.intervals(owned_cpb) ||
        throw(ArgumentError("atom-local child source CPB must match owned support CPB"))
    source_support_count = CartesianRouteCore.support_count(source_cpb)
    source_support_count == CartesianRouteCore.support_count(shellification_region) ||
        throw(ArgumentError("atom-local child source CPB support count does not match owned support"))

    lowering_parameters =
        _cartesian_route_core_staged_field(unit, :lowering_parameters, (;))
    lowering_recipe =
        _cartesian_route_core_staged_field(unit, :lowering_recipe, (;))
    lowering_source = CartesianRouteCore.lowering_source(
        :white_lindsey_atom_local_child_shellification,
        shellification_region,
        (source_cpb,);
        metadata = (;
            source = :atom_growth_plan_unit,
            unit_key = unit.unit_key,
            atom_side,
            staged_lowering_family = unit.lowering_family,
            staged_materialization_dependency = unit.materialization_dependency,
            route_level_source_cpb_count = 1,
            route_level_source_support_count = source_support_count,
            staged_child_source_cpb_count =
                _cartesian_route_core_staged_field(unit, :source_cpb_count, nothing),
            child_source_cpbs_enumerated = false,
            child_shellification_plan_object_kind =
                _cartesian_route_core_staged_field(
                    lowering_parameters,
                    :lowering_piece_object_kind,
                    nothing,
                ),
            child_shellification_plan_role =
                _cartesian_route_core_staged_field(
                    lowering_parameters,
                    :lowering_piece_role,
                    nothing,
                ),
            child_shellification_plan_support_count =
                _cartesian_route_core_staged_field(
                    lowering_parameters,
                    :lowering_piece_support_count,
                    nothing,
                ),
            planned_child_source_cpb_families =
                _cartesian_route_core_staged_field(
                    lowering_recipe,
                    :source_cpb_families,
                    (),
                ),
            staged_source_cpb_enumeration_status =
                _cartesian_route_core_staged_field(
                    lowering_recipe,
                    :source_cpb_enumeration_status,
                    nothing,
                ),
            coefficient_maps_materialized = false,
            source_operator_blocks_materialized = false,
            operator_blocks_materialized = false,
            pair_operator_blocks_materialized = false,
            hamiltonian_data_materialized = false,
        ),
    )
    intermediate_retained_space = CartesianRouteCore.intermediate_retained_space(
        lowering_source;
        retained_rule = :atom_local_child_shellification_sequence,
        source_mode_dims = CartesianRouteCore.shape(source_cpb),
        materialized = false,
        metadata = (;
            source = :atom_growth_plan_unit,
            unit_key = unit.unit_key,
            atom_side,
            status = :deferred_not_materialized,
            child_shellification_sequence_planned = true,
            child_source_cpbs_enumerated = false,
            child_core_policy_available = false,
            coefficient_maps_materialized = false,
        ),
    )
    shell_realization = CartesianRouteCore.trivial_shell_realization(
        intermediate_retained_space,
        shellification_region;
        status = :route_level_atom_local_child_shellification_planned,
        final_dimension = nothing,
        metadata = (;
            source = :atom_growth_plan_unit,
            unit_key = unit.unit_key,
            atom_side,
            shell_row_oracle_authority = false,
            child_coefficient_maps_materialized = false,
        ),
    )
    final_retained_unit = CartesianRouteCore.final_retained_unit(
        unit.unit_key,
        unit.unit_role,
        lowering_source,
        intermediate_retained_space,
        shell_realization;
        metadata = (;
            source = :atom_growth_plan_unit,
            status = :deferred_until_materialization,
            pair_planning_input = true,
            route_level_final_unit_preserves_staged_unit_key = true,
            coefficient_maps_materialized = false,
        ),
    )

    return _cartesian_route_core_sidecar_object(
        ;
        sidecar_source = :atom_growth_atom_local_child_shellification_plan_unit,
        shellification_region,
        lowering_source,
        intermediate_retained_space,
        shell_realization,
        final_retained_unit,
    )
end

function _cartesian_route_core_direct_cpb_sidecar(unit)
    _cartesian_route_core_staged_field(unit, :materialization_dependency, nothing) ===
    :plan_lowerable_direct_slab ||
        throw(ArgumentError("direct CRC sidecar expects a direct slab staged unit"))
    _cartesian_route_core_staged_field(
        _cartesian_route_core_staged_field(unit, :owned_support, (;)),
        :owned_support_is_cpb,
        false,
    ) ||
        throw(ArgumentError("direct CRC sidecar expects coordinate-product owned support"))

    shellification_region =
        _cartesian_route_core_shellification_region_from_staged(unit)
    source_records = _cartesian_route_core_staged_field(unit, :source_cpbs, ())
    length(source_records) == 1 ||
        throw(ArgumentError("direct CRC sidecar expects one source CPB"))
    source_cpb = _cartesian_route_core_cpb_from_staged(only(source_records))
    owned_cpb = only(CartesianRouteCore.owned_support(shellification_region).cpbs)
    CartesianRouteCore.intervals(source_cpb) == CartesianRouteCore.intervals(owned_cpb) ||
        throw(ArgumentError("direct CRC source CPB does not match owned support CPB"))

    lowering_source = CartesianRouteCore.lowering_source(
        :direct_identity_cpb,
        shellification_region,
        (source_cpb,);
        metadata = (;
            source = :atom_growth_plan_unit,
            unit_key = unit.unit_key,
            staged_lowering_family = unit.lowering_family,
            coefficient_maps_materialized = false,
            operator_blocks_materialized = false,
            pair_operator_blocks_materialized = false,
            hamiltonian_data_materialized = false,
        ),
    )
    intermediate_retained_space = CartesianRouteCore.intermediate_retained_space(
        lowering_source;
        retained_rule = :identity_source_modes,
        source_mode_dims = CartesianRouteCore.shape(source_cpb),
        materialized = false,
        metadata = (;
            source = :atom_growth_plan_unit,
            unit_key = unit.unit_key,
            status = :deferred_not_materialized,
        ),
    )
    shell_realization = CartesianRouteCore.trivial_shell_realization(
        intermediate_retained_space,
        shellification_region;
        status = :direct_or_trivial,
        metadata = (;
            source = :atom_growth_plan_unit,
            unit_key = unit.unit_key,
            shell_row_oracle_authority = false,
        ),
    )
    final_retained_unit = CartesianRouteCore.final_retained_unit(
        unit.unit_key,
        unit.unit_role,
        lowering_source,
        intermediate_retained_space,
        shell_realization;
        metadata = (;
            source = :atom_growth_plan_unit,
            status = :deferred_until_materialization,
            pair_planning_input = true,
        ),
    )

    return _cartesian_route_core_sidecar_object(
        ;
        sidecar_source = :atom_growth_direct_cpb_plan_unit,
        shellification_region,
        lowering_source,
        intermediate_retained_space,
        shell_realization,
        final_retained_unit,
    )
end

function _cartesian_route_core_pqs_prototype_sidecar(prototype)
    _cartesian_route_core_staged_field(prototype, :object_kind, nothing) ===
    :cartesian_pqs_lowering_metadata_prototype ||
        throw(ArgumentError("PQS CRC sidecar expects a PQS lowering prototype"))
    shellification_region =
        _cartesian_route_core_shellification_region_from_staged(prototype)
    source_cpb =
        _cartesian_route_core_filled_cpb_from_staged(prototype.source_cpb)
    lowering_source = CartesianRouteCore.pqs_filled_source_lowering(
        shellification_region,
        source_cpb;
        metadata = (;
            source = :pqs_lowering_metadata_prototype,
            unit_key = prototype.unit_key,
            staged_lowering_recipe = prototype.lowering_recipe,
            coefficient_maps_materialized = false,
            operator_blocks_materialized = false,
            pair_operator_blocks_materialized = false,
            hamiltonian_data_materialized = false,
        ),
    )
    intermediate_retained_space = CartesianRouteCore.intermediate_retained_space(
        lowering_source;
        retained_rule = :pqs_boundary_comx_product_modes,
        source_mode_dims = CartesianRouteCore.shape(source_cpb),
        materialized = false,
        metadata = (;
            source = :pqs_lowering_metadata_prototype,
            unit_key = prototype.unit_key,
            selected_modes = :boundary_comx_product_modes,
            status =
                prototype.intermediate_retained_space.status,
            source_space_operator_blocks_planned = true,
        ),
    )
    shell_realization = CartesianRouteCore.pqs_shell_realization(
        intermediate_retained_space,
        shellification_region;
        status = prototype.shell_realization.status,
        metadata = (;
            source = :pqs_lowering_metadata_prototype,
            unit_key = prototype.unit_key,
            shell_projection_planned = true,
            lowdin_cleanup_planned = true,
            materialized = false,
        ),
    )
    final_retained_unit = CartesianRouteCore.final_retained_unit(
        prototype.unit_key,
        prototype.unit_role,
        lowering_source,
        intermediate_retained_space,
        shell_realization;
        metadata = (;
            source = :pqs_lowering_metadata_prototype,
            status = :metadata_only_planned,
            pair_planning_input = true,
        ),
    )

    return _cartesian_route_core_sidecar_object(
        ;
        sidecar_source = :pqs_lowering_metadata_prototype,
        shellification_region,
        lowering_source,
        intermediate_retained_space,
        shell_realization,
        final_retained_unit,
    )
end

function _cartesian_route_core_sidecar(record)
    object_kind = _cartesian_route_core_staged_field(record, :object_kind, nothing)
    if object_kind === :cartesian_atom_growth_plan_unit
        reason = _cartesian_route_core_missing_plan_unit_sidecar_reason(record)
        isnothing(reason) ||
            throw(ArgumentError("unsupported staged plan unit for CRC sidecar: $reason"))
        if _cartesian_route_core_staged_field(
            record,
            :lowering_family,
            nothing,
        ) === :white_lindsey_adaptive_complete_shell
            return _cartesian_route_core_lw_complete_shell_sidecar(record)
        end
        if _cartesian_route_core_staged_field(
            record,
            :lowering_family,
            nothing,
        ) === :outer_mismatch_boundary_slab_set
            return _cartesian_route_core_boundary_slab_set_sidecar(record)
        end
        if _cartesian_route_core_staged_field(
            record,
            :lowering_family,
            nothing,
        ) === :white_lindsey_atom_local_child_shellification
            return _cartesian_route_core_atom_local_child_shellification_sidecar(record)
        end
        return _cartesian_route_core_direct_cpb_sidecar(record)
    end
    if object_kind === :cartesian_pqs_lowering_metadata_prototype
        return _cartesian_route_core_pqs_prototype_sidecar(record)
    end
    throw(ArgumentError("unsupported staged record for CRC sidecar: $object_kind"))
end

function _cartesian_route_core_missing_plan_unit_sidecar_reason(unit)
    _cartesian_route_core_staged_field(unit, :object_kind, nothing) ===
    :cartesian_atom_growth_plan_unit ||
        return :not_an_atom_growth_plan_unit

    lowering_family =
        _cartesian_route_core_staged_field(unit, :lowering_family, nothing)
    materialization_dependency =
        _cartesian_route_core_staged_field(unit, :materialization_dependency, nothing)
    source_cpb_count =
        _cartesian_route_core_staged_field(unit, :source_cpb_count, 0)
    owned_support =
        _cartesian_route_core_staged_field(unit, :owned_support, (;))
    owned_support_is_cpb =
        _cartesian_route_core_staged_field(
            owned_support,
            :owned_support_is_cpb,
            false,
        )

    lowering_family === :white_lindsey_adaptive_complete_shell && return nothing
    if lowering_family === :outer_mismatch_boundary_slab_set &&
       materialization_dependency ===
       :plan_lowerable_outer_mismatch_boundary_slab_set &&
       source_cpb_count > 0
        return nothing
    end
    if materialization_dependency === :plan_lowerable_direct_slab &&
       owned_support_is_cpb &&
       source_cpb_count == 1
        return nothing
    end
    if lowering_family === :white_lindsey_atom_local_child_shellification &&
       materialization_dependency ===
       :plan_lowerable_atom_local_child_shellification &&
       owned_support_is_cpb &&
       !isnothing(_cartesian_route_core_staged_field(unit, :source_box, nothing)) &&
       !isnothing(
           _cartesian_route_core_staged_field(unit, :source_descriptor, nothing),
       )
        return nothing
    end
    lowering_family === :outer_mismatch_boundary_slab_set &&
        return :outer_mismatch_boundary_slab_set_not_yet_mapped_to_crc_final_unit
    lowering_family === :white_lindsey_atom_local_child_shellification &&
        return :atom_local_child_shellification_missing_source_box_descriptor
    materialization_dependency === :plan_lowerable_direct_slab &&
        return :direct_slab_without_single_owned_source_cpb
    return :unsupported_lowering_family
end

function _cartesian_route_core_sidecar_inventory_entry(unit)
    reason = _cartesian_route_core_missing_plan_unit_sidecar_reason(unit)
    if isnothing(reason)
        return (;
            object_kind = :cartesian_route_core_sidecar_inventory_entry,
            unit_key = unit.unit_key,
            unit_role = unit.unit_role,
            lowering_family = unit.lowering_family,
            retirement_target = unit.retirement_target,
            route_core_sidecar_available = true,
            route_core_sidecar = _cartesian_route_core_sidecar(unit),
            missing_route_core_sidecar_reason = nothing,
        )
    end
    return (;
        object_kind = :cartesian_route_core_sidecar_inventory_entry,
        unit_key = unit.unit_key,
        unit_role = unit.unit_role,
        lowering_family = unit.lowering_family,
        retirement_target = unit.retirement_target,
        route_core_sidecar_available = false,
        route_core_sidecar = nothing,
        missing_route_core_sidecar_reason = reason,
    )
end

function _cartesian_route_core_expected_pair_keys(unit_keys)
    keys_tuple = Tuple(unit_keys)
    return Tuple(
        (keys_tuple[left_index], keys_tuple[right_index])
        for left_index in eachindex(keys_tuple)
        for right_index in left_index:length(keys_tuple)
    )
end

function _cartesian_route_core_pair_operator_plan_inventory_metadata(
    crc_pair_inventory,
)
    readiness_requirements =
        CartesianRouteCore.pair_operator_materialization_readiness_requirements()
    if isnothing(crc_pair_inventory)
        return (;
            crc_pair_operator_plan_inventory_available = false,
            crc_pair_operator_plan_inventory_status =
                :blocked_missing_route_core_unit_pair_inventory,
            crc_pair_operator_plan_inventory = nothing,
            crc_pair_operator_plan_count = 0,
            crc_pair_operator_plan_blocker =
                :missing_route_core_unit_pair_inventory,
            crc_pair_operator_plan_blocked_count = 0,
            crc_pair_operator_plan_materialized = false,
            crc_pair_operator_materialization_ready = false,
            crc_pair_operator_materialization_readiness_status =
                :blocked_missing_route_core_unit_pair_inventory,
            crc_pair_operator_materialization_readiness_blocker =
                :missing_route_core_unit_pair_inventory,
            crc_pair_operator_materialization_readiness_requirements =
                readiness_requirements,
            crc_pair_operator_materialization_readiness_plan_count = 0,
            crc_pair_operator_materialization_readiness_blocked_count = 0,
            crc_pair_operator_materialization_readiness_materialized_count = 0,
        )
    end

    plan_inventory = CartesianRouteCore.pair_operator_plan_inventory(
        crc_pair_inventory;
        metadata = (;
            source = :cartesian_route_core_sidecar_inventory,
            status = :metadata_only_adapter_sidecar,
        ),
    )
    plans = CartesianRouteCore.pair_operator_plans(plan_inventory)
    blocked_count =
        count(plan -> !isnothing(CartesianRouteCore.blocker(plan)), plans)
    plan_materialized = any(plan -> plan.materialized, plans)
    readiness =
        CartesianRouteCore.pair_operator_materialization_readiness(plan_inventory)

    return (;
        crc_pair_operator_plan_inventory_available = true,
        crc_pair_operator_plan_inventory_status =
            blocked_count == 0 ?
            :available_route_core_pair_operator_plan_inventory :
            :blocked_route_core_pair_operator_plan_inventory,
        crc_pair_operator_plan_inventory = plan_inventory,
        crc_pair_operator_plan_count =
            CartesianRouteCore.pair_operator_plan_count(plan_inventory),
        crc_pair_operator_plan_blocker =
            CartesianRouteCore.blocker(plan_inventory),
        crc_pair_operator_plan_blocked_count = blocked_count,
        crc_pair_operator_plan_materialized = plan_materialized,
        crc_pair_operator_materialization_ready = readiness.ready,
        crc_pair_operator_materialization_readiness_status =
            readiness.status,
        crc_pair_operator_materialization_readiness_blocker =
            readiness.blocker,
        crc_pair_operator_materialization_readiness_requirements =
            readiness.requirements,
        crc_pair_operator_materialization_readiness_plan_count =
            readiness.plan_count,
        crc_pair_operator_materialization_readiness_blocked_count =
            readiness.blocked_count,
        crc_pair_operator_materialization_readiness_materialized_count =
            readiness.materialized_count,
    )
end

function _cartesian_route_core_blocked_sidecar_inventory(status, blocker)
    crc_pair_operator_plan_metadata =
        _cartesian_route_core_pair_operator_plan_inventory_metadata(nothing)
    return (;
        object_kind = :cartesian_route_core_sidecar_inventory,
        status,
        private_development_only = true,
        sidecar_inventory_source = :blocked_missing_atom_growth_plan_unit_inventory,
        unit_count = 0,
        supported_unit_count = 0,
        unsupported_unit_count = 0,
        entries = (),
        supported_entries = (),
        unsupported_entries = (),
        unsupported_unit_keys = (),
        unsupported_unit_roles = (),
        unsupported_lowering_families = (),
        missing_route_core_sidecar_reasons = (),
        final_units = (),
        final_unit_count = 0,
        crc_pair_inventory_available = false,
        crc_pair_inventory_status = :blocked_missing_atom_growth_plan_unit_inventory,
        crc_pair_inventory = nothing,
        crc_pair_count = 0,
        crc_pair_keys = (),
        staged_pair_keys = (),
        pair_inventory_order_matches_staged = false,
        crc_pair_operator_plan_metadata...,
        pqs_prototype_sidecar_available = false,
        pqs_prototype_sidecar = nothing,
        blocker,
        numerical_behavior_changed = false,
        materialization_behavior_changed = false,
    )
end

function _cartesian_route_core_sidecar_inventory(plan_inventory)
    if isnothing(plan_inventory) ||
       _cartesian_route_core_staged_field(plan_inventory, :status, nothing) !==
       :available_atom_growth_plan_unit_inventory
        return _cartesian_route_core_blocked_sidecar_inventory(
            :blocked_missing_atom_growth_plan_unit_inventory,
            :missing_atom_growth_plan_unit_inventory,
        )
    end

    entries = Tuple(
        _cartesian_route_core_sidecar_inventory_entry(unit)
        for unit in plan_inventory.plan_units
    )
    supported_entries =
        Tuple(entry for entry in entries if entry.route_core_sidecar_available)
    unsupported_entries =
        Tuple(entry for entry in entries if !entry.route_core_sidecar_available)
    final_units = Tuple(
        entry.route_core_sidecar.final_retained_unit
        for entry in supported_entries
    )
    staged_pair_keys =
        _cartesian_route_core_expected_pair_keys(plan_inventory.unit_keys)

    crc_pair_inventory =
        isempty(unsupported_entries) ?
        CartesianRouteCore.unit_pair_inventory(
            final_units;
            metadata = (;
                source = :cartesian_route_core_sidecar_inventory,
                staged_pair_order = :symmetric_upper_triangle,
            ),
        ) :
        nothing
    crc_pair_keys =
        isnothing(crc_pair_inventory) ?
        () :
        CartesianRouteCore.pair_keys(crc_pair_inventory)
    crc_pair_operator_plan_metadata =
        _cartesian_route_core_pair_operator_plan_inventory_metadata(
            crc_pair_inventory,
        )

    pqs_prototype_sidecar =
        _cartesian_route_core_staged_field(
            plan_inventory,
            :pqs_lowering_prototype_available,
            false,
        ) ?
        _cartesian_route_core_sidecar(plan_inventory.pqs_lowering_prototype) :
        nothing

    return (;
        object_kind = :cartesian_route_core_sidecar_inventory,
        status =
            isempty(unsupported_entries) ?
            :available_route_core_sidecar_inventory :
            :blocked_incomplete_route_core_sidecar_inventory,
        private_development_only = true,
        sidecar_inventory_source = :atom_growth_plan_unit_inventory,
        unit_count = length(entries),
        supported_unit_count = length(supported_entries),
        unsupported_unit_count = length(unsupported_entries),
        entries,
        supported_entries,
        unsupported_entries,
        unsupported_unit_keys = Tuple(entry.unit_key for entry in unsupported_entries),
        unsupported_unit_roles = Tuple(entry.unit_role for entry in unsupported_entries),
        unsupported_lowering_families =
            Tuple(entry.lowering_family for entry in unsupported_entries),
        missing_route_core_sidecar_reasons =
            Tuple(entry.missing_route_core_sidecar_reason for entry in unsupported_entries),
        final_units,
        final_unit_count = length(final_units),
        crc_pair_inventory_available = !isnothing(crc_pair_inventory),
        crc_pair_inventory_status =
            isnothing(crc_pair_inventory) ?
            :blocked_missing_unit_sidecars :
            :available_route_core_unit_pair_inventory,
        crc_pair_inventory,
        crc_pair_count =
            isnothing(crc_pair_inventory) ?
            0 :
            length(CartesianRouteCore.pair_entries(crc_pair_inventory)),
        crc_pair_keys,
        staged_pair_keys,
        pair_inventory_order_matches_staged =
            !isnothing(crc_pair_inventory) && crc_pair_keys == staged_pair_keys,
        crc_pair_operator_plan_metadata...,
        pqs_prototype_sidecar_available = !isnothing(pqs_prototype_sidecar),
        pqs_prototype_sidecar,
        blocker =
            isempty(unsupported_entries) ?
            nothing :
            :unsupported_unit_sidecars,
        numerical_behavior_changed = false,
        materialization_behavior_changed = false,
    )
end

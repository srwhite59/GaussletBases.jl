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
    outer_box = _cartesian_route_core_filled_cpb_from_staged(outer_record)
    inner_record =
        _cartesian_route_core_staged_field(
            owned_support_record,
            :inner_exclusion_cpb,
            nothing,
        )
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
        CartesianRouteCore.complete_shell_support(
            outer_box,
            _cartesian_route_core_filled_cpb_from_staged(inner_record);
            support_kind,
            metadata,
        )

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
        return _cartesian_route_core_lw_complete_shell_sidecar(record)
    end
    if object_kind === :cartesian_pqs_lowering_metadata_prototype
        return _cartesian_route_core_pqs_prototype_sidecar(record)
    end
    throw(ArgumentError("unsupported staged record for CRC sidecar: $object_kind"))
end

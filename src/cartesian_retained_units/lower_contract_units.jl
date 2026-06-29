# Selected lowering-contract to retained-unit metadata rules.

const _DIRECT_LOWERING_KINDS = (
    :direct_core_identity_cpb,
)

"""
    retained_unit_plan(lowering_plan; policy = MetadataOnlyRetainedUnits())

Convert selected terminal-lowering contracts into metadata-only retained-unit
records. Only `CartesianTerminalLowering.selected_contracts(lowering_plan)` are
used; unselected available alternatives are not carried downstream as retained
units.
"""
function retained_unit_plan(
    lowering_plan::CartesianTerminalLowering.TerminalLoweringPlan;
    policy::RetainedUnitPolicy = MetadataOnlyRetainedUnits(),
    metadata = (;),
)
    selected = CartesianTerminalLowering.selected_contracts(lowering_plan)
    planned_units = RetainedUnitRecord[]
    for contract in selected
        _append_retained_units_for_contract!(planned_units, contract, policy)
    end
    plan_summary = _retained_unit_plan_summary(policy, lowering_plan, planned_units)
    return RetainedUnitPlan(
        policy,
        lowering_plan,
        planned_units,
        plan_summary,
        NamedTuple(metadata),
    )
end

function _append_retained_units_for_contract!(
    planned_units::Vector{RetainedUnitRecord},
    contract::CartesianTerminalLowering.TerminalLoweringContract,
    policy::MetadataOnlyRetainedUnits,
)
    kind = CartesianTerminalLowering.lowering_kind(contract)
    if kind in _DIRECT_LOWERING_KINDS
        push!(planned_units, _direct_retained_unit(contract, length(planned_units) + 1))
        return planned_units
    end
    if kind === :white_lindsey_boundary_strata
        for (index, source_cpb) in enumerate(contract.source_cpbs)
            push!(
                planned_units,
                _white_lindsey_boundary_stratum_unit(
                    contract,
                    source_cpb,
                    index,
                    length(planned_units) + 1,
                ),
            )
        end
        return planned_units
    end
    if kind === :pqs_filled_source_cpb
        push!(planned_units, _pqs_shell_retained_unit(contract, length(planned_units) + 1))
        return planned_units
    end
    if kind === :distorted_product_box_comx
        push!(
            planned_units,
            _distorted_product_retained_unit(contract, length(planned_units) + 1),
        )
        return planned_units
    end
    if kind === :compact_thin_slab_product_cpb
        push!(planned_units, _compact_thin_slab_retained_unit(contract, length(planned_units) + 1))
        return planned_units
    end
    throw(ArgumentError("unsupported retained-unit lowering kind $kind"))
end

function _direct_retained_unit(
    contract::CartesianTerminalLowering.TerminalLoweringContract,
    unit_index::Int,
)
    return _make_retained_unit(
        contract,
        _unit_key(contract, :retained_unit),
        unit_index,
        :direct_cpb_retained_unit,
        contract.source_cpbs,
        nothing;
        owned_support = contract.owned_support,
        dimension_status = :not_materialized,
        metadata = (;
            identity_like = true,
            selected_contract_granularity = :one_retained_unit,
        ),
    )
end

function _white_lindsey_boundary_stratum_unit(
    contract::CartesianTerminalLowering.TerminalLoweringContract,
    source_cpb,
    source_cpb_index::Int,
    unit_index::Int,
)
    return _make_retained_unit(
        contract,
        _unit_key(contract, :stratum_unit, source_cpb_index),
        unit_index,
        :white_lindsey_boundary_stratum_retained_unit,
        (source_cpb,),
        source_cpb_index;
        owned_support = CartesianRouteCore.owned_cpb(
            source_cpb;
            support_kind = :white_lindsey_boundary_stratum_owned_support,
            metadata = (;
                parent_contract_key = contract.contract_key,
                source_cpb_index,
                terminal_region_key = contract.terminal_region_key,
            ),
        ),
        retained_rule = :white_lindsey_boundary_stratum_product,
        realization_rule = :direct_or_trivial_embedding,
        dimension_status = :not_materialized,
        metadata = (;
            parent_contract_key = contract.contract_key,
            source_cpb_index,
            stratum_kind = _stratum_kind(source_cpb),
            child_of_complete_shell = true,
        ),
    )
end

function _pqs_shell_retained_unit(
    contract::CartesianTerminalLowering.TerminalLoweringContract,
    unit_index::Int,
)
    return _make_retained_unit(
        contract,
        _unit_key(contract, :retained_unit),
        unit_index,
        :pqs_shell_retained_unit,
        contract.source_cpbs,
        nothing;
        owned_support = contract.owned_support,
        retained_rule = :pqs_boundary_comx_product_modes,
        realization_rule = :shell_projection_lowdin,
        dimension_status = :planned_not_materialized,
        metadata = (;
            shell_projection_lowdin_planned = true,
            face_edge_corner_expansion_used = false,
        ),
    )
end

function _distorted_product_retained_unit(
    contract::CartesianTerminalLowering.TerminalLoweringContract,
    unit_index::Int,
)
    retained_metadata = (;
        q = _metadata_value(contract.metadata, :q),
        L = _metadata_value(contract.metadata, :L),
        source_mode_shape = _metadata_value(contract.metadata, :source_mode_shape),
        aspect_ratio = _metadata_value(contract.metadata, :aspect_ratio),
        identity_like = false,
    )
    return _make_retained_unit(
        contract,
        _unit_key(contract, :retained_unit),
        unit_index,
        :distorted_product_box_retained_unit,
        contract.source_cpbs,
        nothing;
        owned_support = contract.owned_support,
        retained_rule = :distorted_product_comx_all_axes,
        realization_rule = contract.realization_rule,
        dimension_status = :planned_not_materialized,
        metadata = retained_metadata,
    )
end

function _compact_thin_slab_retained_unit(
    contract::CartesianTerminalLowering.TerminalLoweringContract,
    unit_index::Int,
)
    return _make_retained_unit(
        contract,
        _unit_key(contract, :retained_unit),
        unit_index,
        :compact_thin_slab_retained_unit,
        contract.source_cpbs,
        nothing;
        owned_support = contract.owned_support,
        retained_rule = :compact_thin_slab_face_product,
        realization_rule = :compact_thin_slab_face_stack,
        dimension_status = :planned_not_materialized,
    )
end

function _make_retained_unit(
    contract::CartesianTerminalLowering.TerminalLoweringContract,
    unit_key::Symbol,
    unit_index::Int,
    unit_kind::Symbol,
    source_cpbs,
    source_cpb_index;
    owned_support,
    retained_rule::Symbol = contract.retained_rule,
    realization_rule = contract.realization_rule,
    dimension_status::Symbol,
    dimension = nothing,
    column_range_status::Symbol = :not_materialized,
    column_range = nothing,
    metadata = (;),
)
    cpb_tuple = Tuple(source_cpbs)
    sidecar = _route_core_final_unit_sidecar(
        unit_key,
        unit_kind,
        contract,
        owned_support,
        cpb_tuple,
        retained_rule,
        realization_rule,
        metadata,
    )
    record_metadata = _merge_metadata(
        contract.metadata,
        metadata,
        (;
            route_core_sidecar_status = sidecar.status,
            route_core_sidecar_blocker = sidecar.blocker,
        ),
    )

    return RetainedUnitRecord(
        unit_key,
        unit_index,
        unit_kind,
        contract.contract_key,
        contract.terminal_region_key,
        contract.terminal_region_role,
        contract.terminal_region_kind,
        contract.lowering_kind,
        retained_rule,
        realization_rule,
        owned_support,
        cpb_tuple,
        source_cpb_index,
        dimension_status,
        dimension,
        column_range_status,
        column_range,
        sidecar.unit,
        false,
        record_metadata,
    )
end

function _route_core_final_unit_sidecar(
    unit_key::Symbol,
    unit_kind::Symbol,
    contract::CartesianTerminalLowering.TerminalLoweringContract,
    owned_support,
    source_cpbs::Tuple,
    retained_rule::Symbol,
    realization_rule,
    metadata,
)
    try
        owned_region = CartesianRouteCore.shellification_region(
            unit_kind,
            owned_support;
            metadata = (;
                terminal_region_key = contract.terminal_region_key,
                terminal_region_role = contract.terminal_region_role,
                terminal_region_kind = contract.terminal_region_kind,
                source_contract_key = contract.contract_key,
            ),
        )

        lowering = _route_core_lowering_source(
            contract,
            owned_region,
            source_cpbs,
            metadata,
        )
        intermediate = CartesianRouteCore.intermediate_retained_space(
            lowering;
            retained_rule,
            dimension = nothing,
            source_mode_dims = nothing,
            materialized = false,
            metadata = (;
                source_contract_key = contract.contract_key,
                dimension_status = :not_materialized,
            ),
        )
        realization = _route_core_shell_realization(
            intermediate,
            owned_region,
            realization_rule,
            metadata,
        )
        final_unit = CartesianRouteCore.final_retained_unit(
            unit_key,
            unit_kind,
            lowering,
            intermediate,
            realization;
            column_range = nothing,
            dimension = nothing,
            metadata = (;
                source_contract_key = contract.contract_key,
                materialized = false,
            ),
        )
        return (; unit = final_unit, status = :available, blocker = nothing)
    catch err
        return (; unit = nothing, status = :blocked, blocker = sprint(showerror, err))
    end
end

function _route_core_lowering_source(
    contract::CartesianTerminalLowering.TerminalLoweringContract,
    owned_region::CartesianRouteCore.ShellificationRegion,
    source_cpbs::Tuple,
    metadata,
)
    if contract.lowering_kind === :pqs_filled_source_cpb && length(source_cpbs) == 1
        return CartesianRouteCore.pqs_filled_source_lowering(
            owned_region,
            only(source_cpbs);
            metadata = _merge_metadata(
                metadata,
                (; source_contract_key = contract.contract_key),
            ),
        )
    end
    return CartesianRouteCore.lowering_source(
        contract.lowering_kind,
        owned_region,
        source_cpbs;
        metadata = _merge_metadata(
            metadata,
            (; source_contract_key = contract.contract_key),
        ),
    )
end

function _route_core_shell_realization(
    intermediate::CartesianRouteCore.IntermediateRetainedSpace,
    owned_region::CartesianRouteCore.ShellificationRegion,
    realization_rule,
    metadata,
)
    if realization_rule === :shell_projection_lowdin
        return CartesianRouteCore.pqs_shell_realization(
            intermediate,
            owned_region;
            status = :planned_shell_projection_lowdin,
            final_dimension = nothing,
            metadata,
        )
    end
    if realization_rule === :direct_or_trivial_embedding
        return CartesianRouteCore.trivial_shell_realization(
            intermediate,
            owned_region;
            status = :planned_not_materialized,
            final_dimension = nothing,
            metadata,
        )
    end
    return CartesianRouteCore.shell_realization(
        intermediate,
        owned_region;
        realization_kind = isnothing(realization_rule) ? :planned_realization : realization_rule,
        status = :planned_not_materialized,
        final_dimension = nothing,
        metadata,
    )
end

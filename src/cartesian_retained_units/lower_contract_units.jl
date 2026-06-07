# Selected lowering-contract to retained-unit metadata rules.

const _DIRECT_LOWERING_KINDS = (
    :direct_core_identity_cpb,
    :direct_slab_identity_cpb,
    :direct_boundary_slab_identity_cpb,
)

"""
    retained_unit_plan(lowering_plan; policy = MetadataOnlyRetainedUnits())

Convert selected terminal-lowering contracts into metadata-only retained-unit
records. Only `CartesianTerminalLowering.selected_contracts(lowering_plan)` are
used; unselected available alternatives are not carried downstream as retained
units.
"""
function retained_unit_plan(
    lowering_plan::CTL.TerminalLoweringPlan;
    policy::RetainedUnitPolicy = MetadataOnlyRetainedUnits(),
    metadata = (;),
)
    selected = CTL.selected_contracts(lowering_plan)
    unit_chunks = Tuple(
        _retained_units_for_contract(contract, policy)
        for contract in selected
    )
    planned_units = Tuple(unit for chunk in unit_chunks for unit in chunk)
    reindexed_units = _reindex_units(planned_units)
    plan_summary = _retained_unit_plan_summary(policy, lowering_plan, reindexed_units)
    return RetainedUnitPlan(
        policy,
        lowering_plan,
        reindexed_units,
        plan_summary,
        NamedTuple(metadata),
    )
end

function _retained_units_for_contract(
    contract::CTL.TerminalLoweringContract,
    ::MetadataOnlyRetainedUnits,
)
    kind = CTL.lowering_kind(contract)
    kind in _DIRECT_LOWERING_KINDS &&
        return (_direct_retained_unit(contract),)
    kind === :white_lindsey_boundary_strata &&
        return _white_lindsey_boundary_stratum_units(contract)
    kind === :pqs_filled_source_cpb &&
        return (_pqs_shell_retained_unit(contract),)
    kind === :distorted_product_box_comx &&
        return (_distorted_product_retained_unit(contract),)
    throw(ArgumentError("unsupported retained-unit lowering kind $kind"))
end

function _reindex_units(planned_units)
    return Tuple(_replace_unit_index(unit, index) for (index, unit) in enumerate(planned_units))
end

function _replace_unit_index(unit::RetainedUnitRecord, unit_index::Int)
    return RetainedUnitRecord(
        unit.unit_key,
        unit_index,
        unit.unit_kind,
        unit.source_contract_key,
        unit.terminal_region_key,
        unit.terminal_region_role,
        unit.terminal_region_kind,
        unit.lowering_kind,
        unit.retained_rule,
        unit.realization_rule,
        unit.owned_support,
        unit.source_cpbs,
        unit.source_cpb_index,
        unit.dimension_status,
        unit.dimension,
        unit.column_range_status,
        unit.column_range,
        unit.route_core_final_unit,
        unit.materialized,
        unit.metadata,
    )
end

function _direct_retained_unit(contract::CTL.TerminalLoweringContract)
    return _make_retained_unit(
        contract,
        _unit_key(contract, :retained_unit),
        0,
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

function _white_lindsey_boundary_stratum_units(contract::CTL.TerminalLoweringContract)
    return Tuple(
        _make_retained_unit(
            contract,
            _unit_key(contract, :stratum_unit, index),
            0,
            :white_lindsey_boundary_stratum_retained_unit,
            (source_cpb,),
            index;
            owned_support = CRC.owned_cpb(
                source_cpb;
                support_kind = :white_lindsey_boundary_stratum_owned_support,
                metadata = (;
                    parent_contract_key = contract.contract_key,
                    source_cpb_index = index,
                    terminal_region_key = contract.terminal_region_key,
                ),
            ),
            retained_rule = :white_lindsey_boundary_stratum_product,
            realization_rule = :direct_or_trivial_embedding,
            dimension_status = :not_materialized,
            metadata = (;
                parent_contract_key = contract.contract_key,
                source_cpb_index = index,
                stratum_kind = _stratum_kind(source_cpb),
                child_of_complete_shell = true,
            ),
        )
        for (index, source_cpb) in enumerate(contract.source_cpbs)
    )
end

function _pqs_shell_retained_unit(contract::CTL.TerminalLoweringContract)
    return _make_retained_unit(
        contract,
        _unit_key(contract, :retained_unit),
        0,
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

function _distorted_product_retained_unit(contract::CTL.TerminalLoweringContract)
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
        0,
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

function _make_retained_unit(
    contract::CTL.TerminalLoweringContract,
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
    contract::CTL.TerminalLoweringContract,
    owned_support,
    source_cpbs::Tuple,
    retained_rule::Symbol,
    realization_rule,
    metadata,
)
    try
        owned_region = CRC.shellification_region(
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
        intermediate = CRC.intermediate_retained_space(
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
        final_unit = CRC.final_retained_unit(
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
    contract::CTL.TerminalLoweringContract,
    owned_region::CRC.ShellificationRegion,
    source_cpbs::Tuple,
    metadata,
)
    if contract.lowering_kind === :pqs_filled_source_cpb && length(source_cpbs) == 1
        return CRC.pqs_filled_source_lowering(
            owned_region,
            only(source_cpbs);
            metadata = _merge_metadata(
                metadata,
                (; source_contract_key = contract.contract_key),
            ),
        )
    end
    return CRC.lowering_source(
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
    intermediate::CRC.IntermediateRetainedSpace,
    owned_region::CRC.ShellificationRegion,
    realization_rule,
    metadata,
)
    if realization_rule === :shell_projection_lowdin
        return CRC.pqs_shell_realization(
            intermediate,
            owned_region;
            status = :planned_shell_projection_lowdin,
            final_dimension = nothing,
            metadata,
        )
    end
    if realization_rule === :direct_or_trivial_embedding
        return CRC.trivial_shell_realization(
            intermediate,
            owned_region;
            status = :planned_not_materialized,
            final_dimension = nothing,
            metadata,
        )
    end
    return CRC.shell_realization(
        intermediate,
        owned_region;
        realization_kind = isnothing(realization_rule) ? :planned_realization : realization_rule,
        status = :planned_not_materialized,
        final_dimension = nothing,
        metadata,
    )
end

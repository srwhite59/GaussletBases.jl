# Retained-unit to transform-contract metadata rules.

"""
    retained_unit_transform_contract_plan(retained_unit_plan; policy)

Convert retained-unit records into one metadata-only transform contract per
retained unit. This is a planning contract only: no transform matrices,
coefficient maps, Lowdin objects, dimensions, or column ranges are built here.
"""
function retained_unit_transform_contract_plan(
    retained_unit_plan::CartesianRetainedUnits.RetainedUnitPlan;
    policy::RetainedUnitTransformContractPolicy =
        MetadataOnlyRetainedUnitTransformContracts(),
    metadata = (;),
)
    units = CartesianRetainedUnits.units(retained_unit_plan)
    contracts = Tuple(
        _retained_unit_transform_contract(unit, policy)
        for unit in units
    )
    plan_summary =
        _retained_unit_transform_contract_plan_summary(
            policy,
            retained_unit_plan,
            contracts,
        )
    return RetainedUnitTransformContractPlan(
        policy,
        retained_unit_plan,
        contracts,
        plan_summary,
        NamedTuple(metadata),
    )
end

function _retained_unit_transform_contract(
    unit::CartesianRetainedUnits.RetainedUnitRecord,
    ::MetadataOnlyRetainedUnitTransformContracts,
)
    transform_path, realization_path, blocker =
        _retained_unit_transform_paths(unit.unit_kind)
    contract_metadata = _merge_metadata(
        unit.metadata,
        (;
            source = :cartesian_retained_unit_transform_contracts,
            source_contract_key = unit.source_contract_key,
            terminal_region_key = unit.terminal_region_key,
            terminal_region_role = unit.terminal_region_role,
            terminal_region_kind = unit.terminal_region_kind,
            source_cpb_index = unit.source_cpb_index,
            retained_unit_materialized = unit.materialized,
        ),
        _raw_product_source_contract_metadata(unit),
    )

    return RetainedUnitTransformContract(
        unit.unit_key,
        unit.unit_index,
        unit.unit_kind,
        unit.lowering_kind,
        unit.retained_rule,
        unit.realization_rule,
        unit.source_cpbs,
        transform_path,
        realization_path,
        unit.dimension_status,
        unit.column_range_status,
        false,
        blocker,
        contract_metadata,
    )
end

function _metadata_value(metadata::NamedTuple, key::Symbol, default = nothing)
    return haskey(metadata, key) ? getfield(metadata, key) : default
end

function _valid_source_mode_dims(value)
    value isa Tuple || return nothing
    length(value) == 3 || return nothing
    all(dim -> dim isa Integer && dim > 0, value) || return nothing
    return Tuple(Int(dim) for dim in value)::NTuple{3,Int}
end

function _pqs_source_mode_dims(unit::CartesianRetainedUnits.RetainedUnitRecord)
    source_mode_shape = _valid_source_mode_dims(
        _metadata_value(unit.metadata, :source_mode_shape),
    )
    !isnothing(source_mode_shape) && return source_mode_shape

    q = _metadata_value(unit.metadata, :q)
    q isa Integer && q > 0 && return (Int(q), Int(q), Int(q))
    return nothing
end

function _pqs_source_cpb(unit::CartesianRetainedUnits.RetainedUnitRecord)
    length(unit.source_cpbs) == 1 || return nothing
    source_cpb = only(unit.source_cpbs)
    source_cpb isa CartesianRawProductSources.CPB.CoordinateProductBox || return nothing
    CartesianRawProductSources.CPB.codimension(source_cpb) == 0 || return nothing
    return source_cpb
end

function _raw_product_source_unavailable_metadata(status::Symbol, blocker::Symbol)
    return (;
        raw_product_source_plan = nothing,
        raw_product_source_summary =
            CartesianRawProductSources.unavailable_summary(status, blocker),
        raw_product_source_plan_status = status,
        raw_product_source_retained_rule = nothing,
        raw_product_source_retained_rule_summary = nothing,
    )
end

function _raw_product_source_contract_metadata(unit::CartesianRetainedUnits.RetainedUnitRecord)
    unit.unit_kind === :pqs_shell_retained_unit || return (;)

    dims = _pqs_source_mode_dims(unit)
    isnothing(dims) && return _raw_product_source_unavailable_metadata(
        :blocked_missing_source_mode_dims,
        :missing_pqs_source_mode_dims,
    )

    source_cpb = _pqs_source_cpb(unit)
    isnothing(source_cpb) && return _raw_product_source_unavailable_metadata(
        :blocked_missing_pqs_source_cpb,
        :missing_single_filled_pqs_source_cpb,
    )

    raw_plan = CartesianRawProductSources.raw_product_box_plan(
        source_cpb;
        source_key = unit.unit_key,
        source_mode_dims = dims,
        metadata = (;
            source = :cartesian_retained_unit_transform_contracts,
            unit_key = unit.unit_key,
            source_contract_key = unit.source_contract_key,
        ),
    )
    raw_summary = CartesianRawProductSources.summary(raw_plan)
    retained_rule =
        CartesianRawProductSources.pqs_boundary_product_mode_retained_rule(raw_plan)
    retained_rule_summary = CartesianRawProductSources.summary(retained_rule)
    return (;
        raw_product_source_plan = raw_plan,
        raw_product_source_summary = raw_summary,
        raw_product_source_plan_status = raw_summary.status,
        raw_product_source_retained_rule = retained_rule,
        raw_product_source_retained_rule_summary = retained_rule_summary,
    )
end

function _retained_unit_transform_paths(unit_kind::Symbol)
    unit_kind === :direct_cpb_retained_unit && return (
        :direct_identity_transform_contract,
        :identity_or_trivial_embedding,
        nothing,
    )
    unit_kind === :white_lindsey_boundary_stratum_retained_unit && return (
        :white_lindsey_boundary_stratum_product_contract,
        :identity_or_trivial_embedding,
        nothing,
    )
    unit_kind === :pqs_shell_retained_unit && return (
        :pqs_source_modes_boundary_selection_shell_realization_contract,
        :shell_projection_lowdin_planned,
        nothing,
    )
    unit_kind === :distorted_product_box_retained_unit && return (
        :distorted_product_comx_contract,
        :distorted_product_realization_planned,
        nothing,
    )
    return (
        :pending_retained_unit_transform_contract,
        :pending_retained_unit_realization_contract,
        :unclassified_retained_unit_transform_contract,
    )
end

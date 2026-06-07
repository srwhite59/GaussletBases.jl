# Retained-unit to transform-contract metadata rules.

"""
    retained_unit_transform_contract_plan(retained_unit_plan; policy)

Convert retained-unit records into one metadata-only transform contract per
retained unit. This is a planning contract only: no transform matrices,
coefficient maps, Lowdin objects, dimensions, or column ranges are built here.
"""
function retained_unit_transform_contract_plan(
    retained_unit_plan::CRU.RetainedUnitPlan;
    policy::RetainedUnitTransformContractPolicy =
        MetadataOnlyRetainedUnitTransformContracts(),
    metadata = (;),
)
    units = CRU.units(retained_unit_plan)
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
    unit::CRU.RetainedUnitRecord,
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

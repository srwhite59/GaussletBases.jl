"""
    _cartesian_terminal_shellification_geometry(parent_axes, nuclear_positions; kwargs...)

Compatibility wrapper for `CartesianShellification.raw_terminal_geometry`.
The implementation lives in `src/cartesian_shellification/`.
"""
function _cartesian_terminal_shellification_geometry(
    parent_axes::NTuple{3,<:AbstractVector},
    nuclear_positions;
    core_side::Int = 5,
    q::Int = core_side,
    bond_axis::Symbol = :auto,
    audit_coverage::Bool = true,
)
    return CartesianShellification.raw_terminal_geometry(
        parent_axes,
        nuclear_positions;
        core_side,
        q,
        bond_axis,
        audit_coverage,
    )
end

function _cartesian_terminal_shellification_geometry_private_summary(plan)
    return CartesianShellification.private_summary(plan)
end

function _cartesian_terminal_shellification_geometry_scaffold(
    plan;
    route_family::Symbol = :white_lindsey_low_order,
)
    return CartesianShellification.scaffold(plan; route_family)
end

function _cartesian_terminal_region_unit_key(region)
    return Symbol(
        "terminal_region_",
        string(region.order_index),
        "_",
        String(region.role),
    )
end

function _cartesian_terminal_region_unit_role(region)
    return Symbol(String(region.role), "_unit")
end

function _cartesian_terminal_region_unit_required_metadata(region, field::Symbol)
    hasproperty(region.metadata, field) ||
        throw(
            ArgumentError(
                "terminal region $(region.role) is missing required metadata field $field",
            ),
        )
    return getproperty(region.metadata, field)
end

function _cartesian_terminal_region_complete_shell_lw_metadata(region)
    isnothing(region.inner_exclusion_box) &&
        throw(ArgumentError("complete shell terminal region is missing inner exclusion box"))
    return (;
        status = :planned_not_materialized,
        lowering_family_planned = :white_lindsey_boundary_stratum_cpbs,
        enumeration_policy = :complete_shell_boundary_strata_metadata,
        facet_count = 6,
        edge_count = 12,
        corner_count = 8,
        total_source_cpb_count = 26,
        source_cpbs_materialized = false,
        coefficient_maps_materialized = false,
        transform_contracts_materialized = false,
        operator_blocks_materialized = false,
        pair_operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
    )
end

function _cartesian_terminal_region_complete_shell_pqs_metadata(region)
    isnothing(region.inner_exclusion_box) &&
        throw(ArgumentError("complete shell terminal region is missing inner exclusion box"))
    return (;
        status = :planned_not_materialized,
        lowering_family_planned = :pqs_filled_source_cpb,
        source_cpb_plan_status = :planned_not_materialized,
        source_cpb_plan_box = region.outer_box,
        source_cpb_plan_kind = :filled_source_cpb,
        source_cpb_count = 1,
        retained_rule = :boundary_comx_product_mode_selection,
        intermediate_retained_space_status = :not_materialized,
        shell_realization_status = :projection_lowdin_planned_not_materialized,
        face_edge_corner_decomposition_required = false,
        coefficient_maps_materialized = false,
        transform_contracts_materialized = false,
        operator_blocks_materialized = false,
        pair_operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
    )
end

function _cartesian_terminal_region_direct_source_metadata(
    region,
    lowering_family_planned::Symbol,
)
    return (;
        status = :planned_not_materialized,
        lowering_family_planned,
        source_cpb_plan_status = :planned_not_materialized,
        source_cpb_plan_box = region.outer_box,
        source_cpb_plan_kind = :direct_coordinate_product_source,
        source_cpb_plan_equals_owned_support = true,
        source_cpb_count = 1,
        retained_rule = :direct_source_modes,
        coefficient_maps_materialized = false,
        transform_contracts_materialized = false,
        operator_blocks_materialized = false,
        pair_operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
    )
end

function _cartesian_terminal_region_unit_mapping(region)
    if region.region_kind in (:direct_core, :direct_atom_contact_core)
        unit_kind =
            region.region_kind == :direct_atom_contact_core ?
            :direct_atom_contact_core_unit :
            :direct_core_unit
        return (;
            unit_kind,
            lowering_family_planned = :direct_core_identity_cpb,
            identity_lowering_planned = true,
            owned_support_is_cpb = false,
            owned_support_status = :owned_support_equals_region_box,
            source_cpb_plan =
                _cartesian_terminal_region_direct_source_metadata(
                    region,
                    :direct_core_identity_cpb,
                ),
            lw_complete_shell_lowering = nothing,
            pqs_complete_shell_lowering = nothing,
        )
    end
    if region.region_kind == :complete_shell
        unit_kind =
            region.role == :shared_molecular_shell ?
            :shared_molecular_complete_shell_unit :
            :atom_local_complete_shell_unit
        return (;
            unit_kind,
            lowering_family_planned = :complete_shell_lowering_recipe_choice_pending,
            identity_lowering_planned = false,
            owned_support_is_cpb = false,
            owned_support_status = :owned_shell_support_outer_minus_inner_exclusion,
            source_cpb_plan = nothing,
            lw_complete_shell_lowering =
                _cartesian_terminal_region_complete_shell_lw_metadata(region),
            pqs_complete_shell_lowering =
                _cartesian_terminal_region_complete_shell_pqs_metadata(region),
        )
    end
    if region.region_kind == :direct_midpoint_slab
        return (;
            unit_kind = :direct_midpoint_slab_unit,
            lowering_family_planned = :direct_slab_identity_cpb,
            identity_lowering_planned = true,
            owned_support_is_cpb = false,
            owned_support_status = :owned_support_equals_region_box,
            source_cpb_plan =
                _cartesian_terminal_region_direct_source_metadata(
                    region,
                    :direct_slab_identity_cpb,
                ),
            lw_complete_shell_lowering = nothing,
            pqs_complete_shell_lowering = nothing,
        )
    end
    if region.region_kind == :central_distorted_product_box
        source_mode_shape =
            _cartesian_terminal_region_unit_required_metadata(
                region,
                :source_mode_shape,
            )
        q = _cartesian_terminal_region_unit_required_metadata(region, :q)
        L = _cartesian_terminal_region_unit_required_metadata(region, :L)
        aspect_ratio =
            _cartesian_terminal_region_unit_required_metadata(region, :aspect_ratio)
        return (;
            unit_kind = :central_distorted_product_box_unit,
            lowering_family_planned = :distorted_product_box_comx,
            identity_lowering_planned = false,
            owned_support_is_cpb = false,
            owned_support_status = :owned_support_equals_region_box,
            source_cpb_plan = merge(
                _cartesian_terminal_region_direct_source_metadata(
                    region,
                    :distorted_product_box_comx,
                ),
                (;
                    source_cpb_plan_kind = :central_distorted_product_box_source,
                    retained_rule = :distorted_product_comx_all_axes,
                    source_mode_shape,
                    q,
                    L,
                    aspect_ratio,
                ),
            ),
            lw_complete_shell_lowering = nothing,
            pqs_complete_shell_lowering = nothing,
        )
    end
    if region.region_kind == :outer_mismatch_slab
        return (;
            unit_kind = :outer_mismatch_boundary_slab_unit,
            lowering_family_planned = :direct_boundary_slab_identity_cpb,
            identity_lowering_planned = true,
            owned_support_is_cpb = false,
            owned_support_status = :owned_support_equals_region_box,
            source_cpb_plan =
                _cartesian_terminal_region_direct_source_metadata(
                    region,
                    :direct_boundary_slab_identity_cpb,
                ),
            lw_complete_shell_lowering = nothing,
            pqs_complete_shell_lowering = nothing,
        )
    end
    throw(
        ArgumentError(
            "unsupported terminal shellification region kind $(region.region_kind)",
        ),
    )
end

function _cartesian_terminal_region_unit_record(region)
    mapping = _cartesian_terminal_region_unit_mapping(region)
    source_mode_shape =
        region.region_kind == :central_distorted_product_box ?
        mapping.source_cpb_plan.source_mode_shape :
        nothing
    q =
        region.region_kind == :central_distorted_product_box ?
        mapping.source_cpb_plan.q :
        nothing
    L =
        region.region_kind == :central_distorted_product_box ?
        mapping.source_cpb_plan.L :
        nothing
    aspect_ratio =
        region.region_kind == :central_distorted_product_box ?
        mapping.source_cpb_plan.aspect_ratio :
        nothing

    return (;
        object_kind = :cartesian_terminal_region_unit_record,
        unit_index = region.order_index,
        unit_key = _cartesian_terminal_region_unit_key(region),
        unit_role = _cartesian_terminal_region_unit_role(region),
        unit_kind = mapping.unit_kind,
        terminal_region_order_index = region.order_index,
        terminal_region_role = region.role,
        terminal_region_kind = region.region_kind,
        outer_box = region.outer_box,
        box = region.outer_box,
        inner_exclusion_box = region.inner_exclusion_box,
        support_count = region.support_count,
        owned_support_status = mapping.owned_support_status,
        owned_support_is_cpb = mapping.owned_support_is_cpb,
        shellification_region_is_cpb = false,
        shellification_region_is_lowering_source = false,
        lowering_family_planned = mapping.lowering_family_planned,
        lowering_recipe_status = :planned_not_materialized,
        source_cpb_plan_status = :planned_not_materialized,
        source_cpb_plan = mapping.source_cpb_plan,
        lw_complete_shell_lowering = mapping.lw_complete_shell_lowering,
        pqs_complete_shell_lowering = mapping.pqs_complete_shell_lowering,
        pqs_filled_source_cpb_plan =
            isnothing(mapping.pqs_complete_shell_lowering) ?
            nothing :
            mapping.pqs_complete_shell_lowering,
        source_mode_shape,
        q,
        L,
        aspect_ratio,
        identity_lowering_planned = mapping.identity_lowering_planned,
        retained_space_status = :not_materialized,
        coefficient_maps_materialized = false,
        transform_contracts_materialized = false,
        operator_blocks_materialized = false,
        pair_operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        final_retained_unit_status =
            :not_available_terminal_region_metadata_only,
        metadata = region.metadata,
        provenance = (;
            source = :_cartesian_terminal_shellification_region_unit_inventory,
            scaffold_region_order_index = region.order_index,
            terminal_region_role = region.role,
            terminal_region_kind = region.region_kind,
        ),
    )
end

function _cartesian_terminal_shellification_region_unit_inventory(scaffold)
    scaffold.object_kind == :cartesian_terminal_shellification_scaffold3d ||
        throw(
            ArgumentError(
                "terminal-region unit inventory requires a cartesian_terminal_shellification_scaffold3d",
            ),
        )
    records = Tuple(
        _cartesian_terminal_region_unit_record(region)
        for region in scaffold.regions
    )
    complete_shell_records =
        Tuple(record for record in records if record.terminal_region_kind == :complete_shell)
    lw_complete_shell_cpb_count =
        sum(
            record -> record.lw_complete_shell_lowering.total_source_cpb_count,
            complete_shell_records;
            init = 0,
        )

    return (;
        object_kind = :cartesian_terminal_region_unit_inventory,
        status = :available_terminal_region_unit_inventory,
        inventory_source = :terminal_shellification_scaffold,
        private_development_only = true,
        terminal_region_count = scaffold.region_count,
        unit_count = length(records),
        terminal_region_units = records,
        unit_records = records,
        unit_keys = Tuple(record.unit_key for record in records),
        unit_roles = Tuple(record.unit_role for record in records),
        unit_kinds = Tuple(record.unit_kind for record in records),
        terminal_region_roles =
            Tuple(record.terminal_region_role for record in records),
        terminal_region_kinds =
            Tuple(record.terminal_region_kind for record in records),
        support_counts = Tuple(record.support_count for record in records),
        all_units_from_terminal_regions = length(records) == scaffold.region_count,
        complete_shell_unit_count = length(complete_shell_records),
        lw_complete_shell_cpb_enumeration_available =
            !isempty(complete_shell_records),
        lw_complete_shell_region_count = length(complete_shell_records),
        lw_complete_shell_cpb_count,
        lw_complete_shell_cpb_family_counts =
            (
                facet_cpb = 6 * length(complete_shell_records),
                edge_cpb = 12 * length(complete_shell_records),
                corner_cpb = 8 * length(complete_shell_records),
            ),
        final_retained_unit_inventory_available = false,
        pair_inventory_available = false,
        pair_inventory_status = :not_available_terminal_region_metadata_only,
        coefficient_maps_materialized = false,
        transform_contracts_materialized = false,
        retained_spaces_materialized = false,
        operator_blocks_materialized = false,
        pair_operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
        diagnostics = (;
            source = :_cartesian_terminal_shellification_region_unit_inventory,
            private_development_only = true,
            terminal_region_metadata_only = true,
            all_units_from_terminal_regions =
                length(records) == scaffold.region_count,
            no_aggregate_atom_box_units =
                all(
                    !(record.unit_key in (:left_atom_box, :right_atom_box))
                    for record in records
                ),
            no_atom_growth_route_units =
                all(
                    record.lowering_family_planned !=
                    :white_lindsey_atom_local_child_shellification
                    for record in records
                ),
            shellification_regions_are_cpbs = false,
            final_retained_unit_inventory_available = false,
            pair_inventory_available = false,
            materialization_behavior_changed = false,
            public_default_behavior_changed = false,
        ),
    )
end

function _cartesian_terminal_region_lowering_contract_base(
    unit,
    contract_index::Int,
    lowering_contract_kind::Symbol,
)
    return (;
        object_kind = :cartesian_terminal_region_lowering_contract,
        status = :planned_not_materialized,
        contract_index,
        contract_key =
            Symbol(String(unit.unit_key), "_", String(lowering_contract_kind)),
        unit_index = unit.unit_index,
        unit_key = unit.unit_key,
        unit_role = unit.unit_role,
        unit_kind = unit.unit_kind,
        terminal_region_order_index = unit.terminal_region_order_index,
        terminal_region_role = unit.terminal_region_role,
        terminal_region_kind = unit.terminal_region_kind,
        lowering_contract_kind,
        outer_box = unit.outer_box,
        box = unit.outer_box,
        inner_exclusion_box = unit.inner_exclusion_box,
        support_count = unit.support_count,
        owned_support_status = unit.owned_support_status,
        owned_support_is_cpb = unit.owned_support_is_cpb,
        shellification_region_is_cpb = unit.shellification_region_is_cpb,
        shellification_region_is_lowering_source =
            unit.shellification_region_is_lowering_source,
        coefficient_maps_materialized = false,
        transform_contracts_materialized = false,
        retained_spaces_materialized = false,
        operator_blocks_materialized = false,
        pair_operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
        final_retained_unit_status = :not_materialized,
        final_retained_unit_records_materialized = false,
        metadata = unit.metadata,
        provenance = (;
            source = :_cartesian_terminal_region_lowering_contract_inventory,
            unit_inventory_source =
                :_cartesian_terminal_shellification_region_unit_inventory,
            unit_key = unit.unit_key,
            terminal_region_role = unit.terminal_region_role,
            terminal_region_kind = unit.terminal_region_kind,
        ),
    )
end

function _cartesian_terminal_region_direct_lowering_contract(
    unit,
    contract_index::Int,
)
    plan = unit.source_cpb_plan
    return merge(
        _cartesian_terminal_region_lowering_contract_base(
            unit,
            contract_index,
            unit.lowering_family_planned,
        ),
        (;
            source_cpb_plan_status = plan.source_cpb_plan_status,
            source_cpb_plan_box = plan.source_cpb_plan_box,
            source_cpb_plan_kind = plan.source_cpb_plan_kind,
            source_cpb_plan_equals_owned_support =
                plan.source_cpb_plan_equals_owned_support,
            source_cpb_count = plan.source_cpb_count,
            source_cpb_family_counts =
                (facet_cpb = 0, edge_cpb = 0, corner_cpb = 0),
            identity_like_source_contract = true,
            retained_rule = plan.retained_rule,
            intermediate_retained_space_status = :not_materialized,
            shell_realization_status = :not_materialized,
            face_edge_corner_decomposition_required = false,
            final_unit_count_planned = 1,
            final_unit_granularity = :single_direct_source_cpb,
        ),
    )
end

function _cartesian_terminal_region_lw_complete_shell_lowering_contract(
    unit,
    contract_index::Int,
)
    lowering = unit.lw_complete_shell_lowering
    return merge(
        _cartesian_terminal_region_lowering_contract_base(
            unit,
            contract_index,
            :white_lindsey_boundary_strata,
        ),
        (;
            source_cpb_plan_status = :planned_not_materialized,
            source_cpb_plan_box = nothing,
            source_cpb_plan_kind = :complete_shell_boundary_strata,
            source_cpb_plan_equals_owned_support = false,
            source_cpb_count = lowering.total_source_cpb_count,
            source_cpb_family_counts =
                (
                    facet_cpb = lowering.facet_count,
                    edge_cpb = lowering.edge_count,
                    corner_cpb = lowering.corner_count,
                ),
            source_cpbs_materialized = false,
            identity_like_source_contract = false,
            retained_rule = :white_lindsey_boundary_strata_children,
            intermediate_retained_space_status = :not_materialized,
            shell_realization_status = :not_materialized,
            face_edge_corner_decomposition_required = true,
            final_unit_count_planned = lowering.total_source_cpb_count,
            final_unit_granularity =
                :white_lindsey_boundary_strata_children_planned,
        ),
    )
end

function _cartesian_terminal_region_pqs_complete_shell_lowering_contract(
    unit,
    contract_index::Int,
)
    lowering = unit.pqs_complete_shell_lowering
    return merge(
        _cartesian_terminal_region_lowering_contract_base(
            unit,
            contract_index,
            :pqs_filled_source_cpb,
        ),
        (;
            source_cpb_plan_status = lowering.source_cpb_plan_status,
            source_cpb_plan_box = lowering.source_cpb_plan_box,
            source_cpb_plan_kind = lowering.source_cpb_plan_kind,
            source_cpb_plan_equals_owned_support = false,
            source_cpb_count = lowering.source_cpb_count,
            source_cpb_family_counts =
                (facet_cpb = 0, edge_cpb = 0, corner_cpb = 0),
            identity_like_source_contract = false,
            retained_rule = lowering.retained_rule,
            intermediate_retained_space_status = :planned_not_materialized,
            shell_realization_status = lowering.shell_realization_status,
            face_edge_corner_decomposition_required =
                lowering.face_edge_corner_decomposition_required,
            final_unit_count_planned = 1,
            final_unit_granularity = :single_pqs_shell_realized_unit_planned,
        ),
    )
end

function _cartesian_terminal_region_distorted_product_box_lowering_contract(
    unit,
    contract_index::Int,
)
    plan = unit.source_cpb_plan
    return merge(
        _cartesian_terminal_region_lowering_contract_base(
            unit,
            contract_index,
            :distorted_product_box_comx,
        ),
        (;
            source_cpb_plan_status = plan.source_cpb_plan_status,
            source_cpb_plan_box = plan.source_cpb_plan_box,
            source_cpb_plan_kind = plan.source_cpb_plan_kind,
            source_cpb_plan_equals_owned_support =
                plan.source_cpb_plan_equals_owned_support,
            source_cpb_count = plan.source_cpb_count,
            source_cpb_family_counts =
                (facet_cpb = 0, edge_cpb = 0, corner_cpb = 0),
            identity_like_source_contract = false,
            retained_rule = plan.retained_rule,
            source_mode_shape = plan.source_mode_shape,
            q = plan.q,
            L = plan.L,
            aspect_ratio = plan.aspect_ratio,
            intermediate_retained_space_status = :planned_not_materialized,
            shell_realization_status = :not_materialized,
            face_edge_corner_decomposition_required = false,
            final_unit_count_planned = 1,
            final_unit_granularity =
                :single_distorted_product_box_unit_planned,
        ),
    )
end

function _cartesian_terminal_region_lowering_contracts_for_unit(
    unit,
    first_contract_index::Int,
)
    records = NamedTuple[]

    if unit.terminal_region_kind in (
        :direct_core,
        :direct_midpoint_slab,
        :outer_mismatch_slab,
    )
        push!(
            records,
            _cartesian_terminal_region_direct_lowering_contract(
                unit,
                first_contract_index,
            ),
        )
    elseif unit.terminal_region_kind == :complete_shell
        push!(
            records,
            _cartesian_terminal_region_lw_complete_shell_lowering_contract(
                unit,
                first_contract_index,
            ),
        )
        push!(
            records,
            _cartesian_terminal_region_pqs_complete_shell_lowering_contract(
                unit,
                first_contract_index + 1,
            ),
        )
    elseif unit.terminal_region_kind == :central_distorted_product_box
        push!(
            records,
            _cartesian_terminal_region_distorted_product_box_lowering_contract(
                unit,
                first_contract_index,
            ),
        )
    else
        throw(
            ArgumentError(
                "unsupported terminal-region unit kind $(unit.terminal_region_kind)",
            ),
        )
    end

    return Tuple(records)
end

function _cartesian_terminal_region_lowering_contract_kind_counts(contracts)
    return (;
        direct_core_identity_cpb_count =
            count(
                contract ->
                    contract.lowering_contract_kind == :direct_core_identity_cpb,
                contracts,
            ),
        direct_slab_identity_cpb_count =
            count(
                contract ->
                    contract.lowering_contract_kind == :direct_slab_identity_cpb,
                contracts,
            ),
        direct_boundary_slab_identity_cpb_count =
            count(
                contract ->
                    contract.lowering_contract_kind ==
                    :direct_boundary_slab_identity_cpb,
                contracts,
            ),
        white_lindsey_boundary_strata_count =
            count(
                contract ->
                    contract.lowering_contract_kind ==
                    :white_lindsey_boundary_strata,
                contracts,
            ),
        pqs_filled_source_cpb_count =
            count(
                contract ->
                    contract.lowering_contract_kind == :pqs_filled_source_cpb,
                contracts,
            ),
        distorted_product_box_comx_count =
            count(
                contract ->
                    contract.lowering_contract_kind == :distorted_product_box_comx,
                contracts,
            ),
    )
end

function _cartesian_terminal_region_lowering_contract_inventory(unit_inventory)
    unit_inventory.object_kind == :cartesian_terminal_region_unit_inventory ||
        throw(
            ArgumentError(
                "terminal-region lowering contracts require a cartesian_terminal_region_unit_inventory",
            ),
        )

    units = unit_inventory.terminal_region_units
    contracts = NamedTuple[]
    for unit in units
        append!(
            contracts,
            _cartesian_terminal_region_lowering_contracts_for_unit(
                unit,
                length(contracts) + 1,
            ),
        )
    end
    lowering_contracts = Tuple(contracts)
    contract_counts_by_unit = Tuple(
        (
            unit_key = unit.unit_key,
            lowering_contract_count =
                count(contract -> contract.unit_key == unit.unit_key, lowering_contracts),
        ) for unit in units
    )
    all(entry -> entry.lowering_contract_count >= 1, contract_counts_by_unit) ||
        throw(ArgumentError("every terminal-region unit must have a lowering contract"))

    complete_shell_contracts = Tuple(
        contract for contract in lowering_contracts
        if contract.terminal_region_kind == :complete_shell
    )
    lw_contracts = Tuple(
        contract for contract in lowering_contracts
        if contract.lowering_contract_kind == :white_lindsey_boundary_strata
    )

    return (;
        object_kind = :cartesian_terminal_region_lowering_contract_inventory,
        status = :available_terminal_region_lowering_contract_inventory,
        inventory_source = :terminal_region_unit_inventory,
        source_object_kind = unit_inventory.object_kind,
        private_development_only = true,
        terminal_region_unit_count = unit_inventory.unit_count,
        lowering_contract_count = length(lowering_contracts),
        lowering_contracts,
        contract_records = lowering_contracts,
        unit_keys = unit_inventory.unit_keys,
        unit_roles = unit_inventory.unit_roles,
        unit_kinds = unit_inventory.unit_kinds,
        terminal_region_roles = unit_inventory.terminal_region_roles,
        terminal_region_kinds = unit_inventory.terminal_region_kinds,
        support_counts = unit_inventory.support_counts,
        contract_counts_by_unit,
        lowering_contract_kinds =
            Tuple(contract.lowering_contract_kind for contract in lowering_contracts),
        lowering_contract_kind_counts =
            _cartesian_terminal_region_lowering_contract_kind_counts(
                lowering_contracts,
            ),
        complete_shell_unit_count = unit_inventory.complete_shell_unit_count,
        complete_shell_lowering_contract_count = length(complete_shell_contracts),
        lw_complete_shell_lowering_contract_count = length(lw_contracts),
        lw_complete_shell_cpb_count =
            sum(contract -> contract.source_cpb_count, lw_contracts; init = 0),
        lw_complete_shell_cpb_family_counts =
            (
                facet_cpb =
                    sum(
                        contract -> contract.source_cpb_family_counts.facet_cpb,
                        lw_contracts;
                        init = 0,
                    ),
                edge_cpb =
                    sum(
                        contract -> contract.source_cpb_family_counts.edge_cpb,
                        lw_contracts;
                        init = 0,
                    ),
                corner_cpb =
                    sum(
                        contract -> contract.source_cpb_family_counts.corner_cpb,
                        lw_contracts;
                        init = 0,
                    ),
            ),
        final_retained_unit_inventory_available = false,
        pair_inventory_available = false,
        pair_inventory_status = :not_available_lowering_contract_metadata_only,
        coefficient_maps_materialized = false,
        transform_contracts_materialized = false,
        retained_spaces_materialized = false,
        operator_blocks_materialized = false,
        pair_operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
        diagnostics = (;
            source = :_cartesian_terminal_region_lowering_contract_inventory,
            private_development_only = true,
            terminal_region_metadata_only = true,
            lowering_contracts_metadata_only = true,
            all_units_have_lowering_contracts = true,
            final_retained_unit_inventory_available = false,
            pair_inventory_available = false,
            coefficient_maps_materialized = false,
            transform_contracts_materialized = false,
            retained_spaces_materialized = false,
            operator_blocks_materialized = false,
            pair_operator_blocks_materialized = false,
            hamiltonian_data_materialized = false,
            artifacts_materialized = false,
            materialization_behavior_changed = false,
            public_default_behavior_changed = false,
        ),
    )
end

function _cartesian_selected_terminal_lowering_contract_kind(
    terminal_region_kind::Symbol,
    route_lowering_family::Symbol,
)
    terminal_region_kind in (:direct_core, :direct_atom_contact_core) &&
        return :direct_core_identity_cpb
    terminal_region_kind == :direct_midpoint_slab &&
        return :direct_slab_identity_cpb
    terminal_region_kind == :outer_mismatch_slab &&
        return :direct_boundary_slab_identity_cpb
    terminal_region_kind == :central_distorted_product_box &&
        return :distorted_product_box_comx
    if terminal_region_kind == :complete_shell
        route_lowering_family == :white_lindsey_low_order &&
            return :white_lindsey_boundary_strata
        route_lowering_family == :pqs && return :pqs_filled_source_cpb
        throw(
            ArgumentError(
                "complete-shell selected lowering requires route_lowering_family :white_lindsey_low_order or :pqs",
            ),
        )
    end
    throw(
        ArgumentError(
            "unsupported terminal-region kind $terminal_region_kind for selected lowering",
        ),
    )
end

function _cartesian_selected_terminal_unit_keys(contracts)
    unit_keys = Symbol[]
    for contract in contracts
        contract.object_kind == :cartesian_terminal_region_lowering_contract ||
            throw(
                ArgumentError(
                    "selected terminal lowering requires terminal-region lowering contract records",
                ),
            )
        contract.unit_key in unit_keys || push!(unit_keys, contract.unit_key)
    end
    return Tuple(unit_keys)
end

function _cartesian_selected_terminal_lowering_contract_for_unit(
    unit_contracts,
    route_lowering_family::Symbol,
)
    isempty(unit_contracts) &&
        throw(ArgumentError("each terminal-region unit must have contract candidates"))
    terminal_region_kind = first(unit_contracts).terminal_region_kind
    all(contract -> contract.terminal_region_kind == terminal_region_kind, unit_contracts) ||
        throw(ArgumentError("contract candidates for one unit disagree on terminal-region kind"))
    selected_kind =
        _cartesian_selected_terminal_lowering_contract_kind(
            terminal_region_kind,
            route_lowering_family,
        )
    selected = Tuple(
        contract for contract in unit_contracts
        if contract.lowering_contract_kind == selected_kind
    )
    length(selected) == 1 ||
        throw(
            ArgumentError(
                "terminal-region unit $(first(unit_contracts).unit_key) must have exactly one selected $selected_kind contract",
            ),
        )
    return only(selected)
end

function _cartesian_selected_terminal_lowering_contract_inventory(
    lowering_contract_inventory,
    route_lowering_family::Symbol,
)
    route_lowering_family in (:white_lindsey_low_order, :pqs) ||
        throw(
            ArgumentError(
                "route_lowering_family must be :white_lindsey_low_order or :pqs",
            ),
        )
    lowering_contract_inventory.object_kind ==
    :cartesian_terminal_region_lowering_contract_inventory ||
        throw(
            ArgumentError(
                "selected terminal lowering contracts require a cartesian_terminal_region_lowering_contract_inventory",
            ),
        )

    lowering_contracts = lowering_contract_inventory.lowering_contracts
    unit_keys = _cartesian_selected_terminal_unit_keys(lowering_contracts)
    terminal_region_unit_count = length(unit_keys)
    selected_contracts = Tuple(
        _cartesian_selected_terminal_lowering_contract_for_unit(
            Tuple(
                contract for contract in lowering_contracts
                if contract.unit_key == unit_key
            ),
            route_lowering_family,
        ) for unit_key in unit_keys
    )
    selected_contract_keys =
        Tuple(contract.contract_key for contract in selected_contracts)
    unselected_contracts = Tuple(
        contract for contract in lowering_contracts
        if !(contract.contract_key in selected_contract_keys)
    )
    selected_contract_counts_by_unit = Tuple(
        (;
            unit_key,
            selected_contract_count =
                count(contract -> contract.unit_key == unit_key, selected_contracts),
        ) for unit_key in unit_keys
    )
    all_units_have_exactly_one_selected_contract =
        all(entry -> entry.selected_contract_count == 1, selected_contract_counts_by_unit)

    return (;
        object_kind = :cartesian_selected_terminal_lowering_contract_inventory,
        status = :available_selected_terminal_lowering_contract_inventory,
        inventory_source = :terminal_region_lowering_contract_inventory,
        source_object_kind = lowering_contract_inventory.object_kind,
        private_development_only = true,
        route_lowering_family,
        terminal_region_unit_count,
        selected_contract_count = length(selected_contracts),
        selected_contracts,
        selected_contract_records = selected_contracts,
        selected_contract_kinds =
            Tuple(contract.lowering_contract_kind for contract in selected_contracts),
        selected_contract_kind_counts =
            _cartesian_terminal_region_lowering_contract_kind_counts(
                selected_contracts,
            ),
        selected_contract_counts_by_unit,
        all_units_have_exactly_one_selected_contract,
        unselected_contract_count = length(unselected_contracts),
        unselected_contracts,
        unselected_contract_kinds =
            Tuple(contract.lowering_contract_kind for contract in unselected_contracts),
        final_retained_unit_inventory_available = false,
        pair_inventory_available = false,
        pair_inventory_status = :not_available_selected_lowering_metadata_only,
        coefficient_maps_materialized = false,
        transform_contracts_materialized = false,
        retained_spaces_materialized = false,
        operator_blocks_materialized = false,
        pair_operator_blocks_materialized = false,
        hamiltonian_data_materialized = false,
        artifacts_materialized = false,
        diagnostics = (;
            source =
                :_cartesian_selected_terminal_lowering_contract_inventory,
            private_development_only = true,
            selected_lowering_metadata_only = true,
            all_units_have_exactly_one_selected_contract,
            final_retained_unit_inventory_available = false,
            pair_inventory_available = false,
            coefficient_maps_materialized = false,
            transform_contracts_materialized = false,
            retained_spaces_materialized = false,
            operator_blocks_materialized = false,
            pair_operator_blocks_materialized = false,
            hamiltonian_data_materialized = false,
            artifacts_materialized = false,
            materialization_behavior_changed = false,
            public_default_behavior_changed = false,
        ),
    )
end

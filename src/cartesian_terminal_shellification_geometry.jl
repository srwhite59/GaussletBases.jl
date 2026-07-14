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

function _cartesian_terminal_region_thin_slab_metadata(region)
    return merge(
        _cartesian_terminal_region_direct_source_metadata(region, :compact_thin_slab_product_cpb),
        (;
            source_cpb_plan_kind = :compact_thin_slab_source,
            retained_rule = :compact_thin_slab_face_product,
        ),
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
            lowering_family_planned = :compact_thin_slab_product_cpb,
            identity_lowering_planned = false,
            owned_support_is_cpb = false,
            owned_support_status = :owned_support_equals_region_box,
            source_cpb_plan = _cartesian_terminal_region_thin_slab_metadata(region),
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
            lowering_family_planned = :compact_thin_slab_product_cpb,
            identity_lowering_planned = false,
            owned_support_is_cpb = false,
            owned_support_status = :owned_support_equals_region_box,
            source_cpb_plan = _cartesian_terminal_region_thin_slab_metadata(region),
            lw_complete_shell_lowering = nothing,
            pqs_complete_shell_lowering = nothing,
        )
    end
    if region.region_kind == :angular_z_extension_slab
        return (;
            unit_kind = :angular_z_extension_slab_unit,
            lowering_family_planned = :compact_thin_slab_product_cpb,
            identity_lowering_planned = false,
            owned_support_is_cpb = false,
            owned_support_status = :owned_support_equals_region_box,
            source_cpb_plan = _cartesian_terminal_region_thin_slab_metadata(region),
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
    records = [
        _cartesian_terminal_region_unit_record(region)
        for region in scaffold.regions
    ]
    complete_shell_records =
        [record for record in records if record.terminal_region_kind == :complete_shell]
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
        unit_keys = [record.unit_key for record in records],
        unit_roles = [record.unit_role for record in records],
        unit_kinds = [record.unit_kind for record in records],
        terminal_region_roles =
            [record.terminal_region_role for record in records],
        terminal_region_kinds =
            [record.terminal_region_kind for record in records],
        support_counts = [record.support_count for record in records],
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

function _cartesian_terminal_region_lowering_contract_kind_counts(contracts)
    return (;
        direct_core_identity_cpb_count =
            count(
                contract ->
                    contract.lowering_contract_kind == :direct_core_identity_cpb,
                contracts,
            ),
        compact_thin_slab_product_cpb_count =
            count(
                contract ->
                    contract.lowering_contract_kind ==
                    :compact_thin_slab_product_cpb,
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

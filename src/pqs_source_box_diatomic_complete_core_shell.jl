function _pqs_source_box_route_driver_box_support_indices(
    box::NTuple{3,UnitRange{Int}},
    inner_box,
    parent_dims::NTuple{3,Int},
)
    outer_indices =
        _nested_box_support_indices(box[1], box[2], box[3], parent_dims)
    isnothing(inner_box) && return outer_indices
    inner = Set(
        _nested_box_support_indices(
            inner_box[1],
            inner_box[2],
            inner_box[3],
            parent_dims,
        ),
    )
    return Int[index for index in outer_indices if !(index in inner)]
end

function _pqs_source_box_route_driver_terminal_support_record(
    contract,
    parent_dims::NTuple{3,Int},
)
    support_indices = _pqs_source_box_route_driver_box_support_indices(
        contract.outer_box,
        contract.inner_exclusion_box,
        parent_dims,
    )
    length(support_indices) == contract.support_count ||
        throw(
            ArgumentError(
                "terminal support count mismatch for $(contract.unit_key)",
            ),
        )
    return (;
        unit_key = contract.unit_key,
        unit_role = contract.unit_role,
        terminal_region_key = contract.terminal_region_key,
        terminal_region_role = contract.terminal_region_role,
        terminal_region_kind = contract.terminal_region_kind,
        lowering_contract_kind = contract.lowering_contract_kind,
        outer_box = contract.outer_box,
        inner_exclusion_box = contract.inner_exclusion_box,
        support_indices,
        support_states = NTuple{3,Int}[
            _cartesian_unflat_index(index, parent_dims) for index in support_indices
        ],
        support_count = length(support_indices),
        source_cpb = hasproperty(contract, :source_cpb) ? contract.source_cpb : nothing,
        source_cpb_plan_box =
            hasproperty(contract, :source_cpb_plan_box) ?
            contract.source_cpb_plan_box :
            nothing,
        source_mode_shape =
            hasproperty(contract, :source_mode_shape) ?
            contract.source_mode_shape :
            nothing,
        retained_rule =
            hasproperty(contract, :retained_rule) ?
            contract.retained_rule :
            nothing,
        bond_axis = get(contract.metadata, :bond_axis, nothing),
    )
end

function _pqs_source_box_route_driver_terminal_support_coverage(records, parent_dims)
    expected_count = prod(parent_dims)
    seen = Set{Int}()
    duplicate_count = 0
    outside_count = 0
    for record in records
        for index in record.support_indices
            if index < 1 || index > expected_count
                outside_count += 1
                continue
            end
            if index in seen
                duplicate_count += 1
            else
                push!(seen, index)
            end
        end
    end
    missing_count = expected_count - length(seen)
    return (;
        coverage_complete =
            missing_count == 0 && duplicate_count == 0 && outside_count == 0,
        duplicate_count,
        missing_count,
        outside_count,
    )
end

function _pqs_source_box_route_driver_terminal_topology_support_region_plan(
    parent,
    low_order_assembly,
)
    blocked(blocker) = (;
        status = :blocked_terminal_topology_support_region_plan,
        blocker,
        authority = :terminal_lowering_contract_inventory,
        support_order = (),
        support_counts = (;),
        counts_generated = false,
        counts_source = :terminal_topology_support_region_plan_blocked,
        coverage_complete = false,
        duplicate_count = nothing,
        missing_count = nothing,
        outside_count = nothing,
    )
    isnothing(low_order_assembly) &&
        return blocked(:missing_low_order_assembly_terminal_topology)
    parent_axis_counts =
        hasproperty(parent, :axis_counts) ? parent.axis_counts : nothing
    parent_dims =
        _pqs_source_box_route_driver_axis_counts_tuple(parent_axis_counts)
    isnothing(parent_dims) && return blocked(:missing_terminal_parent_axis_counts)
    lowering_inventory =
        hasproperty(low_order_assembly, :lowering_contract_inventory) ?
        low_order_assembly.lowering_contract_inventory :
        nothing
    isnothing(lowering_inventory) &&
        return blocked(:missing_terminal_lowering_contract_inventory)
    supported_lowering_kinds = (
        :direct_core_identity_cpb,
        :direct_slab_identity_cpb,
        :direct_boundary_slab_identity_cpb,
        :pqs_filled_source_cpb,
        :distorted_product_box_comx,
    )
    unsupported_lowering_kinds = Tuple(
        unique(
            contract.lowering_contract_kind for
            contract in lowering_inventory.lowering_contracts
            if !(contract.lowering_contract_kind in supported_lowering_kinds)
        ),
    )
    isempty(unsupported_lowering_kinds) ||
        return merge(
            blocked(:unsupported_terminal_lowering_kind),
            (; unsupported_lowering_kinds),
        )
    contracts = Tuple(lowering_inventory.lowering_contracts)
    isempty(contracts) && return blocked(:missing_terminal_support_contracts)
    records = Tuple(
        _pqs_source_box_route_driver_terminal_support_record(
            contract,
            parent_dims,
        ) for contract in contracts
    )
    support_order = Tuple(record.unit_key for record in records)
    length(unique(support_order)) == length(support_order) ||
        return blocked(:duplicate_terminal_support_unit_keys)
    support_counts =
        NamedTuple{support_order}(Tuple(record.support_count for record in records))
    coverage =
        _pqs_source_box_route_driver_terminal_support_coverage(records, parent_dims)
    coverage.coverage_complete ||
        return merge(
            blocked(:terminal_support_coverage_incomplete),
            (;
                duplicate_count = coverage.duplicate_count,
                missing_count = coverage.missing_count,
                outside_count = coverage.outside_count,
            ),
        )
    bond_axes = unique(
        Any[record.bond_axis for record in records if !isnothing(record.bond_axis)],
    )
    length(bond_axes) <= 1 || return blocked(:inconsistent_terminal_bond_axes)
    bond_axis = isempty(bond_axes) ? nothing : only(bond_axes)
    return (;
        status = :available_terminal_topology_support_region_plan,
        blocker = nothing,
        authority = :terminal_lowering_contract_inventory,
        support_order,
        support_counts,
        terminal_support_order = support_order,
        terminal_support_counts = support_counts,
        terminal_support_records = records,
        terminal_region_roles =
            Tuple(record.terminal_region_role for record in records),
        terminal_region_kinds =
            Tuple(record.terminal_region_kind for record in records),
        lowering_contract_kinds =
            Tuple(record.lowering_contract_kind for record in records),
        counts_generated = true,
        counts_source = :terminal_lowering_contract_inventory,
        coverage_complete = coverage.coverage_complete,
        duplicate_count = coverage.duplicate_count,
        missing_count = coverage.missing_count,
        outside_count = coverage.outside_count,
        parent_dims,
        bond_axis,
    )
end

function _pqs_source_box_route_driver_terminal_retained_rule_record(
    support_record,
    retained_unit,
    transform_contract,
    order_index::Int,
)
    kind = support_record.lowering_contract_kind
    blocked(blocker) = (;
        order_index,
        support_record,
        retained_unit_key = isnothing(retained_unit) ? nothing : retained_unit.unit_key,
        transform_contract_unit_key =
            isnothing(transform_contract) ? nothing : transform_contract.unit_key,
        role = support_record.terminal_region_role,
        lowering_contract_kind = kind,
        support_count = support_record.support_count,
        retained_count = nothing,
        transform_kind = :not_available,
        realization_status = :blocked_terminal_retained_rule,
        blocker,
    )
    if kind in (
        :direct_core_identity_cpb,
        :direct_slab_identity_cpb,
        :direct_boundary_slab_identity_cpb,
    )
        return (;
            order_index,
            support_record,
            retained_unit_key = retained_unit.unit_key,
            transform_contract_unit_key = transform_contract.unit_key,
            role = support_record.terminal_region_role,
            lowering_contract_kind = kind,
            support_count = support_record.support_count,
            retained_count = support_record.support_count,
            transform_kind = transform_contract.transform_path,
            realization_status = :identity_retained_sector_available,
            blocker = nothing,
        )
    elseif kind === :pqs_filled_source_cpb
        retained_rule =
            get(transform_contract.metadata, :raw_product_source_retained_rule, nothing)
        isnothing(retained_rule) &&
            return blocked(:missing_raw_product_source_retained_rule)
        return (;
            order_index,
            support_record,
            retained_unit_key = retained_unit.unit_key,
            transform_contract_unit_key = transform_contract.unit_key,
            role = support_record.terminal_region_role,
            lowering_contract_kind = kind,
            support_count = support_record.support_count,
            retained_count = Int(retained_rule.retained_count),
            transform_kind = transform_contract.transform_path,
            realization_status = :retained_rule_available_not_realized,
            blocker = nothing,
        )
    elseif kind === :distorted_product_box_comx
        return blocked(:distorted_product_realization_missing)
    end
    return blocked(:unsupported_terminal_lowering_kind)
end

function _pqs_source_box_route_driver_terminal_retained_records_by_region(records)
    by_region = Dict{Symbol,Vector{Any}}()
    for record in records
        units = get!(by_region, record.terminal_region_key, Any[])
        push!(units, record)
    end
    return by_region
end

function _pqs_source_box_route_driver_terminal_transform_contracts_by_unit(
    contracts,
)
    by_unit = Dict{Symbol,Any}()
    for contract in contracts
        haskey(by_unit, contract.unit_key) &&
            throw(ArgumentError("duplicate retained-unit transform contract key"))
        by_unit[contract.unit_key] = contract
    end
    return by_unit
end

function _pqs_source_box_route_driver_terminal_retained_rule_dimension_budget(records)
    budget = NamedTuple[]
    cumulative = 0
    for record in records
        retained_count = record.retained_count
        if !isnothing(retained_count)
            cumulative += retained_count
        end
        push!(
            budget,
            (;
                order_index = record.order_index,
                role = record.role,
                lowering_contract_kind = record.lowering_contract_kind,
                support_count = record.support_count,
                retained_count,
                transform_kind = record.transform_kind,
                cumulative_retained_dimension = cumulative,
                blocker = record.blocker,
            ),
        )
    end
    return Tuple(budget), cumulative
end

function _pqs_source_box_route_driver_terminal_retained_rule_plan(
    parent,
    low_order_stage,
    retained_unit_plan,
    retained_unit_transform_contract_plan,
)
    blocked(blocker) = (;
        status = :blocked_terminal_retained_rule_plan,
        blocker,
        authority = :cartesian_retained_unit_transform_contract_plan,
        records = (),
        dimension_budget = (),
        retained_order = (),
        retained_counts = (;),
        total_retained_dimension = nothing,
    )
    support_plan =
        _pqs_source_box_route_driver_terminal_topology_support_region_plan(
            parent,
            low_order_stage,
        )
    support_plan.status === :available_terminal_topology_support_region_plan ||
        return blocked(support_plan.blocker)
    retained_unit_plan isa CartesianRetainedUnits.RetainedUnitPlan ||
        return blocked(:missing_retained_unit_plan)
    retained_unit_transform_contract_plan isa
    CartesianRetainedUnitTransformContracts.RetainedUnitTransformContractPlan ||
        return blocked(:missing_retained_unit_transform_contract_plan)
    support_records = support_plan.terminal_support_records
    retained_units_by_region =
        _pqs_source_box_route_driver_terminal_retained_records_by_region(
            CartesianRetainedUnits.units(retained_unit_plan),
        )
    transform_contracts_by_unit =
        _pqs_source_box_route_driver_terminal_transform_contracts_by_unit(
            CartesianRetainedUnitTransformContracts.transform_contracts(
                retained_unit_transform_contract_plan,
            ),
        )
    records = Tuple(
        begin
            retained_units =
                get(retained_units_by_region, support_record.terminal_region_key, Any[])
            if length(retained_units) != 1
                (;
                    order_index,
                    support_record,
                    retained_unit_key = nothing,
                    transform_contract_unit_key = nothing,
                    role = support_record.terminal_region_role,
                    lowering_contract_kind = support_record.lowering_contract_kind,
                    support_count = support_record.support_count,
                    retained_count = nothing,
                    transform_kind = :not_available,
                    realization_status = :blocked_terminal_retained_rule,
                    blocker = :terminal_retained_unit_join_mismatch,
                )
            else
                retained_unit = only(retained_units)
                transform_contract =
                    get(transform_contracts_by_unit, retained_unit.unit_key, nothing)
                if isnothing(transform_contract)
                    (;
                        order_index,
                        support_record,
                        retained_unit_key = retained_unit.unit_key,
                        transform_contract_unit_key = nothing,
                        role = support_record.terminal_region_role,
                        lowering_contract_kind = support_record.lowering_contract_kind,
                        support_count = support_record.support_count,
                        retained_count = nothing,
                        transform_kind = :not_available,
                        realization_status = :blocked_terminal_retained_rule,
                        blocker = :missing_retained_unit_transform_contract,
                    )
                else
                    _pqs_source_box_route_driver_terminal_retained_rule_record(
                        support_record,
                        retained_unit,
                        transform_contract,
                        order_index,
                    )
                end
            end
        end for (order_index, support_record) in pairs(support_records)
    )
    blocked_records = Tuple(record for record in records if !isnothing(record.blocker))
    budget, total_retained_dimension =
        _pqs_source_box_route_driver_terminal_retained_rule_dimension_budget(records)
    retained_order = Tuple(record.support_record.unit_key for record in records)
    retained_counts = NamedTuple{retained_order}(
        Tuple(record.retained_count for record in records),
    )
    isempty(blocked_records) ||
        return (;
            status = :blocked_terminal_retained_rule_plan,
            blocker = first(blocked_records).blocker,
            authority = :cartesian_retained_unit_transform_contract_plan,
            records,
            dimension_budget = budget,
            retained_order,
            retained_counts,
            total_retained_dimension,
            support_plan,
        )
    return (;
        status = :available_terminal_retained_rule_plan,
        blocker = nothing,
        authority = :cartesian_retained_unit_transform_contract_plan,
        records,
        dimension_budget = budget,
        retained_order,
        retained_counts,
        total_retained_dimension,
        support_plan,
    )
end

function _pqs_source_box_route_driver_axis_counts_tuple(axis_counts)
    isnothing(axis_counts) && return nothing
    if axis_counts isa NamedTuple
        return (Int(axis_counts.x), Int(axis_counts.y), Int(axis_counts.z))
    end
    axis_counts_tuple = Tuple(axis_counts)
    length(axis_counts_tuple) == 3 || return nothing
    return ntuple(axis -> Int(axis_counts_tuple[axis]), 3)
end

function _pqs_source_box_route_driver_diatomic_center_summary(parent)
    center_table =
        hasproperty(parent, :center_table) ? Tuple(parent.center_table) : ()
    return (;
        center_count = length(center_table),
        center_keys = Tuple(
            _pqs_source_box_route_driver_descriptor_property(
                center,
                :center_key,
                index,
            ) for (index, center) in pairs(center_table)
        ),
        nuclear_charges = Tuple(
            _pqs_source_box_route_driver_descriptor_property(
                center,
                :nuclear_charge,
            ) for center in center_table
        ),
        locations = Tuple(
            _pqs_source_box_route_driver_descriptor_property(center, :location)
            for center in center_table
        ),
    )
end
struct _PQSDiatomicPhysicalGaussletCoreShellTargetPayload
    status::Symbol
    blocker
    route_family::Symbol
    route_kind::Symbol
    parent_axis_counts
    support_units::Tuple
    retained_units::Tuple
    support_counts
    retained_counts
    retained_order::Tuple
    expected_final_dimension
    retained_atom_core_interiors::Bool
    supplement_policy::Symbol
    summary
    metadata
end

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

struct _PQSDiatomicPhysicalGaussletSupplementRequestPayload
    status::Symbol
    blocker
    route_family::Symbol
    route_kind::Symbol
    fixture_label::Symbol
    supplement_policy::Symbol
    atom_symbols::Tuple
    nuclear_charges::Tuple
    atom_locations::Tuple
    bond_axis
    bond_length
    basis_name
    lmax
    uncontracted
    residual_keep_policy::Symbol
    residual_drop_tolerance
    required_provider_blocks::Tuple
    summary
    metadata
end

struct _PQSDiatomicPhysicalGaussletSupplementRepresentationPayload
    status::Symbol
    blocker
    route_family::Symbol
    route_kind::Symbol
    fixture_label::Symbol
    supplement_policy::Symbol
    object_kind::Symbol
    supplement
    representation
    basis_name
    lmax
    uncontracted
    atom_symbols::Tuple
    center_count
    orbital_count
    summary
    metadata
end

struct _PQSDiatomicPhysicalGaussletSupplementPreflightPayload
    status::Symbol
    blocker
    route_family::Symbol
    route_kind::Symbol
    fixture_label::Symbol
    support_counts
    retained_counts
    retained_order::Tuple
    retained_transform_kind::Symbol
    gausslet_final_dimension
    supplement_policy::Symbol
    summary
    metadata
end

struct _PQSDiatomicPhysicalGaussletCoreShellSourcePlanPayload
    status::Symbol
    blocker
    route_family::Symbol
    source_plan
    summary
    metadata
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

function _pqs_source_box_route_driver_ordered_count_tuple(counts, order)
    isnothing(counts) && return nothing
    isempty(order) && return ()
    counts isa NamedTuple || return nothing
    all(hasproperty(counts, key) for key in order) || return nothing
    return Tuple(Int(getproperty(counts, key)) for key in order)
end

function _pqs_source_box_route_driver_location_tuple(location)
    isnothing(location) && return nothing
    location_tuple = Tuple(location)
    length(location_tuple) == 3 || return nothing
    return ntuple(axis -> Float64(location_tuple[axis]), 3)
end

function _pqs_source_box_route_driver_h2_r4_center_match(parent)
    summary = _pqs_source_box_route_driver_diatomic_center_summary(parent)
    summary.center_count == 2 || return false
    all(charge -> isapprox(Float64(charge), 1.0; atol = 1.0e-12),
        summary.nuclear_charges) || return false
    locations =
        _pqs_source_box_route_driver_location_tuple.(summary.locations)
    all(!isnothing, locations) || return false
    ordered = sort(collect(locations); by = location -> location[3])
    expected = ((0.0, 0.0, -2.0), (0.0, 0.0, 2.0))
    return all(
        isapprox(ordered[index][axis], expected[index][axis]; atol = 1.0e-12)
        for index in 1:2 for axis in 1:3
    )
end

function _pqs_source_box_route_driver_diatomic_physical_gausslet_target_payload(
    parent,
    route_skeleton,
    recipe,
    low_order_assembly = nothing,
)
    route_family =
        hasproperty(route_skeleton, :route_family) ?
        route_skeleton.route_family :
        recipe.route_family
    route_kind =
        hasproperty(route_skeleton, :route_kind) ?
        route_skeleton.route_kind :
        recipe.route_kind
    system_classification =
        hasproperty(parent, :system_classification) ?
        parent.system_classification :
        nothing
    parent_axis_counts =
        hasproperty(route_skeleton, :parent_axis_counts) ?
        route_skeleton.parent_axis_counts :
        hasproperty(parent, :axis_counts) ?
        parent.axis_counts :
        nothing
    parent_axis_count_tuple =
        _pqs_source_box_route_driver_axis_counts_tuple(parent_axis_counts)
    inventory =
        hasproperty(route_skeleton, :physical_target_inventory) ?
        route_skeleton.physical_target_inventory :
        nothing
    independent_target = route_kind ===
                         :bond_aligned_diatomic_independent_pqs_source_box_core_shell
    expected_route_kinds =
        (:bond_aligned_diatomic_independent_pqs_source_box_core_shell,)

    if route_family !== :pqs_source_box
        status = :not_applicable_physical_gausslet_target_non_pqs_route
        blocker = nothing
    elseif !(route_kind in expected_route_kinds)
        status = :not_applicable_physical_gausslet_target_route_kind
        blocker = nothing
    elseif system_classification !== :bond_aligned_diatomic
        status = :blocked_physical_gausslet_target_inventory
        blocker = :not_bond_aligned_diatomic
    elseif isnothing(inventory)
        status = :blocked_physical_gausslet_target_inventory
        blocker = :missing_physical_gausslet_target_inventory
    else
        status = :available_physical_gausslet_core_shell_target_inventory
        blocker = nothing
    end

    support_units =
        isnothing(inventory) ? () : Tuple(inventory.support_units)
    retained_units =
        isnothing(inventory) ? () : Tuple(inventory.retained_units)
    generated_support_plan =
        independent_target ?
        _pqs_source_box_route_driver_terminal_topology_support_region_plan(
            parent,
            low_order_assembly,
        ) :
        nothing
    generated_support_available =
        !isnothing(generated_support_plan) &&
        generated_support_plan.status ===
        :available_terminal_topology_support_region_plan
    retained_rule_plan = nothing
    terminal_retained_rule_plan =
        independent_target &&
        !isnothing(low_order_assembly) &&
        hasproperty(low_order_assembly, :terminal_retained_rule_plan) ?
        low_order_assembly.terminal_retained_rule_plan :
        nothing
    terminal_retained_available =
        !isnothing(terminal_retained_rule_plan) &&
        terminal_retained_rule_plan.status === :available_terminal_retained_rule_plan
    terminal_retained_blocker =
        isnothing(terminal_retained_rule_plan) ?
        :missing_terminal_retained_rule_plan :
        terminal_retained_rule_plan.blocker
    if independent_target && terminal_retained_available
        status = :blocked_physical_gausslet_target_inventory
        blocker = :missing_terminal_source_plan_realization
    elseif independent_target && generated_support_available
        status = :blocked_physical_gausslet_target_inventory
        blocker = terminal_retained_blocker
    end
    support_units =
        generated_support_available ? generated_support_plan.support_order :
        support_units
    support_counts =
        generated_support_available ? generated_support_plan.support_counts :
        isnothing(inventory) ? (;) : inventory.support_counts
    retained_units =
        terminal_retained_available ? terminal_retained_rule_plan.retained_order :
        retained_units
    retained_counts =
        terminal_retained_available ? terminal_retained_rule_plan.retained_counts :
        isnothing(inventory) ? (;) : inventory.retained_counts
    retained_order =
        terminal_retained_available ? terminal_retained_rule_plan.retained_order :
        isnothing(inventory) ? () : Tuple(inventory.retained_order)
    expected_final_dimension =
        terminal_retained_available ?
        terminal_retained_rule_plan.total_retained_dimension :
        isnothing(inventory) ? nothing : inventory.expected_final_dimension
    retained_atom_core_interiors =
        (!isnothing(inventory) && inventory.retained_atom_core_interiors)
    source_plan_blocker =
        terminal_retained_available ?
        :missing_terminal_source_plan_realization :
        generated_support_available ?
        terminal_retained_blocker :
        nothing
    retained_transform_authority =
        terminal_retained_available ?
        :terminal_retained_rule_preflight :
        generated_support_available ?
        terminal_retained_blocker :
        isnothing(inventory) ?
        :not_available :
        :pqs_source_box_construction
    primary_blocker =
        generated_support_available ?
        source_plan_blocker :
        blocker
    secondary_blocker =
        nothing
    support_plan =
        !isnothing(generated_support_plan) ? generated_support_plan :
        isnothing(inventory) ? nothing : get(inventory, :support_plan, nothing)
    supplement_policy =
        independent_target ?
        something(get(recipe, :supplement_policy, nothing), :none) :
        isnothing(inventory) ?
        :not_available :
        something(get(recipe, :supplement_policy, nothing), inventory.supplement_policy)
    summary = (;
        status,
        blocker,
        route_family,
        route_kind,
        parent_axis_counts,
        support_units,
        support_counts,
        retained_units,
        retained_counts,
        retained_order,
        expected_final_dimension,
        retained_atom_core_interiors,
        supplement_policy,
        target_inventory_available =
            status === :available_physical_gausslet_core_shell_target_inventory,
        source_plan_materialized = false,
        retained_transform_authority,
        primary_blocker,
        secondary_blocker,
        independent_source_plan_blocker = source_plan_blocker,
        support_plan,
        terminal_retained_rule_plan,
        retained_rule_plan = nothing,
    )
    metadata = (;
        source =
            :pqs_source_box_route_driver_diatomic_physical_gausslet_target_payload,
        route_kind,
        provenance =
            isnothing(inventory) ? :missing_inventory : inventory.provenance,
        target_inventory_hard_coded_from_reviewed_contract = !independent_target,
        reviewed_contract_pass = independent_target ? nothing : 200,
        old_wl_qw_fixed_block_size = independent_target ? nothing : (1215, 463),
    )

    return _PQSDiatomicPhysicalGaussletCoreShellTargetPayload(
        status,
        blocker,
        route_family,
        route_kind,
        parent_axis_counts,
        support_units,
        retained_units,
        support_counts,
        retained_counts,
        retained_order,
        expected_final_dimension,
        retained_atom_core_interiors,
        supplement_policy,
        summary,
        metadata,
    )
end

function _pqs_source_box_route_driver_diatomic_bond_length(atom_locations)
    length(atom_locations) == 2 || return nothing
    a = atom_locations[1]
    b = atom_locations[2]
    length(a) == 3 && length(b) == 3 || return nothing
    return sqrt(sum((Float64(a[axis]) - Float64(b[axis]))^2 for axis in 1:3))
end

function _pqs_source_box_route_driver_supplement_constructor_basis_name(basis_name)
    isnothing(basis_name) && return basis_name
    parts = split(String(basis_name), '/')
    return String(parts[end])
end

function _pqs_source_box_route_driver_diatomic_fixture_label(target_payload)
    isnothing(target_payload) && return :not_available
    route_kind =
        hasproperty(target_payload, :route_kind) ?
        target_payload.route_kind :
        :not_available
    if route_kind === :bond_aligned_diatomic_independent_pqs_source_box_core_shell
        return :bond_aligned_diatomic_terminal_topology
    end
    return :h2_r4_physical_gausslet_q5
end

function _pqs_source_box_route_driver_missing_legacy_basis_error(error)
    error isa ArgumentError || return false
    return occursin("could not find legacy basis block", sprint(showerror, error))
end

function _pqs_source_box_route_driver_diatomic_physical_gausslet_supplement_request_payload(
    parent,
    target_payload,
    recipe = (;),
)
    target_available =
        !isnothing(target_payload) &&
        target_payload.status === :available_physical_gausslet_core_shell_target_inventory
    center_table =
        hasproperty(parent, :center_table) ? Tuple(parent.center_table) : ()
    atom_symbols =
        Tuple(
            _pqs_source_box_route_driver_descriptor_property(
                center,
                :atom_symbol,
                :unknown,
            ) for center in center_table
        )
    nuclear_charges =
        Tuple(
            _pqs_source_box_route_driver_descriptor_property(
                center,
                :nuclear_charge,
                nothing,
            ) for center in center_table
        )
    atom_locations =
        Tuple(
            _pqs_source_box_route_driver_location_tuple(
                _pqs_source_box_route_driver_descriptor_property(
                    center,
                    :location,
                    nothing,
                ),
            ) for center in center_table
        )
    bond_axis = hasproperty(parent, :bond_axis) ? parent.bond_axis : nothing
    bond_length = _pqs_source_box_route_driver_diatomic_bond_length(atom_locations)
    fixture_label =
        _pqs_source_box_route_driver_diatomic_fixture_label(target_payload)
    supplement_policy =
        target_available ? target_payload.supplement_policy : :not_available
    route_family = target_available ? target_payload.route_family : :not_available
    route_kind = target_available ? target_payload.route_kind : :not_available
    basis_name =
        supplement_policy === :mwg_residual_gto ?
        something(get(recipe, :supplement_basis, nothing), "H/cc-pVTZ") :
        nothing
    lmax =
        supplement_policy === :mwg_residual_gto ?
        something(get(recipe, :supplement_lmax, nothing), 1) :
        nothing
    uncontracted =
        supplement_policy === :mwg_residual_gto ?
        something(get(recipe, :supplement_uncontracted, nothing), false) :
        nothing
    residual_keep_policy =
        supplement_policy === :mwg_residual_gto ?
        :route_private_mwg_residual_gto_preflight_only :
        :not_applicable
    residual_drop_tolerance =
        supplement_policy === :mwg_residual_gto ? 1.0e-10 : nothing
    required_provider_blocks =
        supplement_policy === :mwg_residual_gto ?
        (
            :mixed_gausslet_gto_blocks,
            :gto_gto_blocks,
            :combined_raw_moment_matrices,
            :residual_mwg_representation,
            :combined_density_density_readiness,
        ) :
        ()

    if !target_available
        status = :blocked_pqs_physical_gausslet_supplement_request
        blocker = :missing_physical_gausslet_target_inventory
    elseif supplement_policy === :none
        status = :not_requested
        blocker = nothing
    elseif supplement_policy === :mwg_residual_gto
        status = :available_pqs_physical_gausslet_supplement_request
        blocker = nothing
    else
        status = :blocked_pqs_physical_gausslet_supplement_request
        blocker = :unsupported_physical_gausslet_supplement_policy
    end

    summary = (;
        status,
        blocker,
        route_family,
        route_kind,
        fixture_label,
        supplement_policy,
        atom_symbols,
        nuclear_charges,
        atom_locations,
        bond_axis,
        bond_length,
        basis_name,
        lmax,
        uncontracted,
        residual_keep_policy,
        residual_drop_tolerance,
        required_provider_blocks,
    )
    metadata = (;
        source =
            :pqs_source_box_route_driver_diatomic_physical_gausslet_supplement_request_payload,
        boundary = :request_metadata_only,
        route_private = true,
        gto_mwg_materialization = false,
    )
    return _PQSDiatomicPhysicalGaussletSupplementRequestPayload(
        status,
        blocker,
        route_family,
        route_kind,
        fixture_label,
        supplement_policy,
        atom_symbols,
        nuclear_charges,
        atom_locations,
        bond_axis,
        bond_length,
        basis_name,
        lmax,
        uncontracted,
        residual_keep_policy,
        residual_drop_tolerance,
        required_provider_blocks,
        summary,
        metadata,
    )
end

function _pqs_source_box_route_driver_diatomic_supplement_atom_string(atom)
    atom isa AbstractString && return String(atom)
    atom isa Symbol && return String(atom)
    return String(atom)
end

function _pqs_source_box_route_driver_diatomic_physical_gausslet_supplement_representation_payload(
    request_payload,
)
    if isnothing(request_payload)
        status = :not_available
        blocker = :missing_physical_gausslet_supplement_request_payload
        route_family = :not_available
        route_kind = :not_available
        fixture_label = :not_available
        supplement_policy = :not_available
        object_kind = :not_available
        supplement = nothing
        representation = nothing
        basis_name = nothing
        lmax = nothing
        uncontracted = nothing
        atom_symbols = ()
    else
        route_family = request_payload.route_family
        route_kind = request_payload.route_kind
        fixture_label = request_payload.fixture_label
        supplement_policy = request_payload.supplement_policy
        basis_name = request_payload.basis_name
        lmax = request_payload.lmax
        uncontracted = request_payload.uncontracted
        atom_symbols = request_payload.atom_symbols
        if supplement_policy === :none
            status = :not_requested
            blocker = nothing
            object_kind = :not_requested
            supplement = nothing
            representation = nothing
        elseif request_payload.status !==
               :available_pqs_physical_gausslet_supplement_request
            status = :blocked_pqs_physical_gausslet_gto_supplement_representation
            blocker = request_payload.blocker
            object_kind = :not_available
            supplement = nothing
            representation = nothing
        else
            locations = request_payload.atom_locations
            locations_valid =
                length(locations) == 2 &&
                length(atom_symbols) == 2 &&
                all(
                    location -> !isnothing(location) && length(location) == 3,
                    locations,
                )
            if !locations_valid
                status = :blocked_pqs_physical_gausslet_gto_supplement_representation
                blocker = :missing_h2_supplement_nuclei
                object_kind = :not_available
                supplement = nothing
                representation = nothing
            else
                nuclei = NTuple{3,Float64}[
                    ntuple(axis -> Float64(location[axis]), 3)
                    for location in locations
                ]
                constructor_basis =
                    _pqs_source_box_route_driver_supplement_constructor_basis_name(
                        basis_name,
                    )
                atom_name =
                    _pqs_source_box_route_driver_diatomic_supplement_atom_string.(
                        atom_symbols,
                    )
                try
                    supplement =
                        atom_name[1] == atom_name[2] ?
                        legacy_bond_aligned_diatomic_gaussian_supplement(
                            atom_name[1],
                            constructor_basis,
                            nuclei;
                            lmax,
                            uncontracted,
                        ) :
                        legacy_bond_aligned_heteronuclear_gaussian_supplement(
                            atom_name[1],
                            constructor_basis,
                            atom_name[2],
                            constructor_basis,
                            nuclei;
                            lmax,
                            uncontracted,
                        )
                    representation = basis_representation(supplement)
                    status =
                        :available_pqs_physical_gausslet_gto_supplement_representation
                    blocker = nothing
                    object_kind =
                        :cartesian_gaussian_shell_supplement_representation
                catch err
                    if _pqs_source_box_route_driver_missing_legacy_basis_error(err)
                        status =
                            :blocked_pqs_physical_gausslet_gto_supplement_representation
                        blocker = :missing_gto_supplement_basis
                        object_kind = :not_available
                        supplement = nothing
                        representation = nothing
                    else
                        rethrow()
                    end
                end
            end
        end
    end

    center_count = isnothing(representation) ? 0 : length(representation.metadata.nuclei)
    orbital_count = isnothing(representation) ? 0 : length(representation.orbitals)
    summary = (;
        status,
        blocker,
        route_family,
        route_kind,
        fixture_label,
        supplement_policy,
        object_kind,
        basis_name,
        lmax,
        uncontracted,
        atom_symbols,
        center_count,
        orbital_count,
    )
    metadata = (;
        source =
            :pqs_source_box_route_driver_diatomic_physical_gausslet_supplement_representation_payload,
        boundary = :route_private_matrix_free_representation,
    )
    return _PQSDiatomicPhysicalGaussletSupplementRepresentationPayload(
        status,
        blocker,
        route_family,
        route_kind,
        fixture_label,
        supplement_policy,
        object_kind,
        supplement,
        representation,
        basis_name,
        lmax,
        uncontracted,
        atom_symbols,
        center_count,
        orbital_count,
        summary,
        metadata,
    )
end

function _pqs_source_box_route_driver_diatomic_physical_gausslet_supplement_preflight_payload(
    target_payload,
    request_payload = nothing,
    representation_payload = nothing,
)
    target_available =
        !isnothing(target_payload) &&
        target_payload.status === :available_physical_gausslet_core_shell_target_inventory
    fixture_label =
        _pqs_source_box_route_driver_diatomic_fixture_label(target_payload)
    retained_transform_kind = :pqs
    if !target_available
        status = :blocked_pqs_physical_gausslet_supplement_preflight
        blocker = :missing_physical_gausslet_target_inventory
    elseif target_payload.supplement_policy === :none
        status = :not_requested
        blocker = nothing
    elseif target_payload.supplement_policy === :mwg_residual_gto
        status = :blocked_pqs_physical_gausslet_mwg_residual_gto_preflight
        request_blocker =
            !isnothing(request_payload) &&
            request_payload.status !== :available_pqs_physical_gausslet_supplement_request ?
            request_payload.blocker :
            nothing
        representation_blocker =
            isnothing(request_blocker) &&
            (
                isnothing(representation_payload) ||
                representation_payload.status !==
                :available_pqs_physical_gausslet_gto_supplement_representation
            ) ?
            :missing_gto_supplement_representation :
            nothing
        blocker =
            !isnothing(request_blocker) ?
            request_blocker :
            !isnothing(representation_blocker) ?
            representation_blocker :
            :missing_provider_gto_supplement_blocks
    else
        status = :blocked_pqs_physical_gausslet_supplement_preflight
        blocker = :unsupported_physical_gausslet_supplement_policy
    end
    support_counts =
        target_available ?
        _pqs_source_box_route_driver_ordered_count_tuple(
            target_payload.support_counts,
            target_payload.support_units,
        ) :
        nothing
    retained_counts =
        target_available ?
        _pqs_source_box_route_driver_ordered_count_tuple(
            target_payload.retained_counts,
            target_payload.retained_order,
        ) :
        nothing
    retained_order = target_available ? target_payload.retained_order : ()
    gausslet_final_dimension =
        target_available ? target_payload.expected_final_dimension : nothing
    supplement_policy =
        target_available ? target_payload.supplement_policy : :not_available
    summary = (;
        status,
        blocker,
        route_family =
            target_available ? target_payload.route_family : :not_available,
        route_kind =
            target_available ? target_payload.route_kind : :not_available,
        fixture_label,
        support_counts,
        retained_counts,
        retained_order,
        retained_transform_kind,
        gausslet_final_dimension,
        supplement_policy,
    )
    metadata = (;
        source =
            :pqs_source_box_route_driver_diatomic_physical_gausslet_supplement_preflight_payload,
        boundary = :metadata_readiness_only,
        gto_mwg_materialization = false,
    )
    return _PQSDiatomicPhysicalGaussletSupplementPreflightPayload(
        status,
        blocker,
        summary.route_family,
        summary.route_kind,
        fixture_label,
        support_counts,
        retained_counts,
        retained_order,
        retained_transform_kind,
        gausslet_final_dimension,
        supplement_policy,
        summary,
        metadata,
    )
end

function _pqs_source_box_route_driver_diatomic_physical_gausslet_source_plan_payload(
    parent,
    target_payload,
    low_order_assembly = nothing,
)
    target_available =
        !isnothing(target_payload) &&
        target_payload.status === :available_physical_gausslet_core_shell_target_inventory
    independent_target =
        !isnothing(target_payload) &&
        target_payload.route_kind ===
        :bond_aligned_diatomic_independent_pqs_source_box_core_shell
    source_plan = nothing
    independent_source_plan_available = false
    independent_source_plan_blocker =
        !isnothing(target_payload) ?
        get(
            target_payload.summary,
            :independent_source_plan_blocker,
            :missing_independent_pqs_physical_source_plan_materializer,
        ) :
        :missing_independent_pqs_physical_source_plan_materializer
    status =
        isnothing(source_plan) ?
        :blocked_pqs_diatomic_physical_gausslet_core_shell_source_plan :
        source_plan.status
    blocker =
        !isnothing(source_plan) ?
        nothing :
        independent_target ?
        independent_source_plan_blocker :
        target_available ?
        :missing_atom_contact_core_support_rows :
        :missing_physical_gausslet_target_inventory
    summary = (;
        object_kind = :pqs_diatomic_physical_gausslet_core_shell_source_plan,
        status,
        blocker,
        target_status = isnothing(target_payload) ? :not_available : target_payload.status,
        parent_axis_counts = isnothing(target_payload) ? nothing : target_payload.parent_axis_counts,
        support_order = isnothing(target_payload) ? () : target_payload.support_units,
        retained_order = isnothing(target_payload) ? () : target_payload.retained_order,
        support_counts = isnothing(target_payload) ? (;) : target_payload.support_counts,
        retained_counts = isnothing(target_payload) ? (;) : target_payload.retained_counts,
        expected_final_dimension =
            isnothing(target_payload) ? nothing : target_payload.expected_final_dimension,
        retained_atom_core_interiors =
            !isnothing(target_payload) && target_payload.retained_atom_core_interiors,
        supplement_policy =
            isnothing(target_payload) ? :not_available : target_payload.supplement_policy,
        source_plan_materialized = !isnothing(source_plan),
        retained_transform_authority =
            isnothing(target_payload) ?
            :not_available :
            get(target_payload.summary, :retained_transform_authority, :not_available),
        secondary_blocker =
            isnothing(target_payload) ? nothing : get(target_payload.summary, :secondary_blocker, nothing),
        independent_source_plan_blocker =
            isnothing(target_payload) ?
            nothing :
            get(target_payload.summary, :independent_source_plan_blocker, nothing),
    )
    metadata = (;
        source =
            :pqs_source_box_route_driver_diatomic_physical_gausslet_source_plan_payload,
        route_owned = true,
        diagnostic_221_source_plan_reused = false,
        placeholders_synthesized = false,
    )

    return _PQSDiatomicPhysicalGaussletCoreShellSourcePlanPayload(
        status,
        blocker,
        isnothing(target_payload) ? :not_available : target_payload.route_family,
        source_plan,
        summary,
        metadata,
    )
end

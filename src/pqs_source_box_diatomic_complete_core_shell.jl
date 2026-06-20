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

struct _PQSIndependentH2PQSSupplementSupportPartitionPayload
    status::Symbol
    blocker
    route_family::Symbol
    route_kind::Symbol
    support_order::Tuple
    retained_order::Tuple
    support_counts
    retained_counts
    retained_ranges
    unit_partitions::Tuple
    total_tile_count::Int
    total_support_count::Int
    coverage_complete::Bool
    duplicate_parent_row_count::Int
    missing_parent_row_count::Int
    outside_parent_row_count::Int
end

function _pqs_source_box_route_driver_support_tile(
    unit_key::Symbol,
    source_region_role::Symbol,
    tile_index::Int,
    box::NTuple{3,UnitRange{Int}},
    parent_dims::NTuple{3,Int},
)
    tile_key = Symbol(unit_key, "_", source_region_role, "_tile_", tile_index)
    support_indices = _nested_box_support_indices(box[1], box[2], box[3], parent_dims)
    support_states = NTuple{3,Int}[
        _cartesian_unflat_index(index, parent_dims) for index in support_indices
    ]
    cpb = CartesianCPB.cpb(
        box...;
        role = tile_key,
        metadata = (;
            unit_key,
            source_region_role,
            source = :independent_h2_pqs_supplement_support_partition,
        ),
    )
    return (;
        tile_key,
        unit_key,
        source_region_role,
        support_tile_role = :rectangular_parent_row_support_tile,
        cpb,
        cpb_role = tile_key,
        intervals = box,
        support_indices,
        support_states,
        support_count = length(support_indices),
        parent_row_min = minimum(support_indices),
        parent_row_max = maximum(support_indices),
        provider_tile_ready = true,
    )
end

function _pqs_source_box_route_driver_nonempty_box_push!(
    boxes,
    x::UnitRange{Int},
    y::UnitRange{Int},
    z::UnitRange{Int},
)
    (isempty(x) || isempty(y) || isempty(z)) && return boxes
    push!(boxes, (x, y, z))
    return boxes
end

function _pqs_source_box_route_driver_rectangular_difference_boxes(
    outer_box::NTuple{3,UnitRange{Int}},
    inner_box::NTuple{3,UnitRange{Int}},
)
    boxes = NTuple{3,UnitRange{Int}}[]
    _pqs_source_box_route_driver_nonempty_box_push!(
        boxes,
        first(outer_box[1]):(first(inner_box[1]) - 1),
        outer_box[2],
        outer_box[3],
    )
    _pqs_source_box_route_driver_nonempty_box_push!(
        boxes,
        (last(inner_box[1]) + 1):last(outer_box[1]),
        outer_box[2],
        outer_box[3],
    )
    _pqs_source_box_route_driver_nonempty_box_push!(
        boxes,
        inner_box[1],
        first(outer_box[2]):(first(inner_box[2]) - 1),
        outer_box[3],
    )
    _pqs_source_box_route_driver_nonempty_box_push!(
        boxes,
        inner_box[1],
        (last(inner_box[2]) + 1):last(outer_box[2]),
        outer_box[3],
    )
    _pqs_source_box_route_driver_nonempty_box_push!(
        boxes,
        inner_box[1],
        inner_box[2],
        first(outer_box[3]):(first(inner_box[3]) - 1),
    )
    _pqs_source_box_route_driver_nonempty_box_push!(
        boxes,
        inner_box[1],
        inner_box[2],
        (last(inner_box[3]) + 1):last(outer_box[3]),
    )
    return Tuple(boxes)
end

function _pqs_source_box_route_driver_shared_shell_support_tiles(
    role::Symbol,
    outer_box::NTuple{3,UnitRange{Int}},
    inner_box::NTuple{3,UnitRange{Int}},
    parent_dims::NTuple{3,Int},
)
    outer_cpb = CartesianCPB.filled_cpb(
        outer_box...;
        role = Symbol(role, "_outer_support_box"),
    )
    inner_cpb = CartesianCPB.filled_cpb(
        inner_box...;
        role = Symbol(role, "_inner_exclusion_box"),
    )
    boxes = try
        strata = CartesianCPB.complete_shell_boundary_strata(outer_cpb, inner_cpb)
        Tuple(CartesianCPB.intervals(cpb) for cpb in strata.all_strata)
    catch
        _pqs_source_box_route_driver_rectangular_difference_boxes(outer_box, inner_box)
    end
    return Tuple(
        _pqs_source_box_route_driver_support_tile(
            role,
            :shared_molecular_shell,
            index,
            box,
            parent_dims,
        ) for (index, box) in pairs(boxes)
    )
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
    seen = Set{Int}()
    duplicate_count = 0
    for record in records
        for index in record.support_indices
            if index in seen
                duplicate_count += 1
            else
                push!(seen, index)
            end
        end
    end
    expected_count = prod(parent_dims)
    return (;
        coverage_complete = length(seen) == expected_count && duplicate_count == 0,
        duplicate_count,
        missing_count = expected_count - length(seen),
        outside_count = count(index -> index < 1 || index > expected_count, seen),
    )
end

function _pqs_source_box_route_driver_h2_terminal_support_plan(records, parent_dims)
    contact_records = Tuple(
        record for record in records
        if record.terminal_region_role == :atom_contact_core &&
           record.lowering_contract_kind == :direct_core_identity_cpb
    )
    shared_records = sort(
        collect(
            record for record in records
            if record.terminal_region_role == :shared_molecular_shell &&
               record.lowering_contract_kind == :pqs_filled_source_cpb
        );
        by = record -> prod(length.(record.outer_box)),
        rev = true,
    )
    length(contact_records) == 1 && length(shared_records) == 2 &&
        length(records) == 3 ||
        return nothing
    support_order = (:atom_contact_core, :shared_shell_1, :shared_shell_2)
    support_counts = (;
        atom_contact_core = contact_records[1].support_count,
        shared_shell_1 = shared_records[1].support_count,
        shared_shell_2 = shared_records[2].support_count,
    )
    atom_record = contact_records[1]
    atom_contact_core_descriptor = (;
        role = :atom_contact_core,
        source_region_roles = (atom_record.terminal_region_role,),
        support_tiles = (
            _pqs_source_box_route_driver_support_tile(
                :atom_contact_core,
                atom_record.terminal_region_role,
                1,
                atom_record.outer_box,
                parent_dims,
            ),
        ),
        tile_count = 1,
        tile_support_counts = (atom_record.support_count,),
        support_indices = atom_record.support_indices,
        support_states = atom_record.support_states,
        support_count = atom_record.support_count,
        coefficient_representation = :sparse_parent_row_direct_selector,
    )
    shared_descriptor(role, record) = begin
        support_tiles =
            _pqs_source_box_route_driver_shared_shell_support_tiles(
                role,
                record.outer_box,
                record.inner_exclusion_box,
                parent_dims,
            )
        (;
            role,
            terminal_region_key = record.terminal_region_key,
            current_box = record.outer_box,
            inner_box = record.inner_exclusion_box,
            source_cpb = record.source_cpb,
            support_tiles,
            tile_count = length(support_tiles),
            tile_support_counts = Tuple(tile.support_count for tile in support_tiles),
            support_indices = record.support_indices,
            support_states = record.support_states,
            support_count = record.support_count,
        )
    end
    return (;
        support_order,
        support_counts,
        atom_contact_core_descriptor,
        shared_shell_descriptors = (;
            shared_shell_1 = shared_descriptor(:shared_shell_1, shared_records[1]),
            shared_shell_2 = shared_descriptor(:shared_shell_2, shared_records[2]),
        ),
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
    contracts = Tuple(
        contract for contract in lowering_inventory.lowering_contracts
        if contract.lowering_contract_kind in (
            :direct_core_identity_cpb,
            :direct_slab_identity_cpb,
            :direct_boundary_slab_identity_cpb,
            :pqs_filled_source_cpb,
        )
    )
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
    bond_axes = unique(
        Any[record.bond_axis for record in records if !isnothing(record.bond_axis)],
    )
    length(bond_axes) <= 1 || return blocked(:inconsistent_terminal_bond_axes)
    bond_axis = isempty(bond_axes) ? nothing : only(bond_axes)
    h2_plan =
        coverage.coverage_complete ?
        _pqs_source_box_route_driver_h2_terminal_support_plan(records, parent_dims) :
        nothing
    h2_available = !isnothing(h2_plan)
    status =
        h2_available ?
        :available_independent_pqs_support_region_plan :
        :available_terminal_topology_support_region_plan
    return merge(
        (;
            status,
            blocker = nothing,
            authority = :terminal_lowering_contract_inventory,
            support_order = h2_available ? h2_plan.support_order : support_order,
            support_counts = h2_available ? h2_plan.support_counts : support_counts,
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
        ),
        h2_available ?
        (;
            atom_contact_core_descriptor = h2_plan.atom_contact_core_descriptor,
            shared_shell_descriptors = h2_plan.shared_shell_descriptors,
        ) :
        (;),
    )
end

function _pqs_source_box_route_driver_independent_h2_retained_rule_plan(
    support_plan,
)
    support_order = (:atom_contact_core, :shared_shell_1, :shared_shell_2)
    isnothing(support_plan) && return nothing
    support_plan.status === :available_independent_pqs_support_region_plan ||
        return nothing
    support_plan.support_order == support_order ||
        throw(ArgumentError("independent H2 retained rule support order mismatch"))
    support_counts = support_plan.support_counts
    Tuple(values(support_counts)) == (275, 578, 362) ||
        throw(ArgumentError("independent H2 retained rule support count mismatch"))

    shared_shell_rule =
        CartesianRawProductSources.pqs_boundary_product_mode_retained_rule(
            (5, 5, 5);
            source_key = :independent_h2_shared_shell_q5_source,
            metadata = (;
                route_unit_role = :shared_molecular_shell,
                source = :independent_h2_retained_rule_readiness,
            ),
        )
    retained_counts = (;
        atom_contact_core = support_counts.atom_contact_core,
        shared_shell_1 = shared_shell_rule.retained_count,
        shared_shell_2 = shared_shell_rule.retained_count,
    )
    Tuple(values(retained_counts)) == (275, 98, 98) ||
        throw(ArgumentError("independent H2 retained rule count mismatch"))
    return (;
        support_order,
        support_counts,
        retained_order = support_order,
        retained_counts,
        expected_final_dimension = sum(values(retained_counts)),
        per_unit_transform_kind = (;
            atom_contact_core = :identity_source_modes,
            shared_shell_1 = shared_shell_rule.transform_kind,
            shared_shell_2 = shared_shell_rule.transform_kind,
        ),
        shared_shell_source_mode_dims = shared_shell_rule.source_mode_dims,
        shared_shell_retained_rule_kind = shared_shell_rule.retained_rule_kind,
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

struct _PQSH2WLGaussletOnlyReferenceCandidatePayload
    status::Symbol
    blocker
    route_family::Symbol
    route_kind::Symbol
    label::String
    retained_transform_kind::Symbol
    supplement_policy::Symbol
    final_dimension
    conditions
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

struct _PQSDiatomicPhysicalGaussletFinalBasisPayload
    status::Symbol
    blocker
    route_family::Symbol
    source_plan
    source_plan_status::Symbol
    final_basis
    final_basis_status::Symbol
    summary
    metadata
end

struct _PQSDiatomicPhysicalGaussletH1Payload
    status::Symbol
    blocker
    route_family::Symbol
    source_plan
    source_plan_status::Symbol
    final_basis
    final_basis_status::Symbol
    support_kinetic
    support_kinetic_status::Symbol
    support_electron_nuclear_by_center
    support_electron_nuclear_status::Symbol
    final_kinetic
    final_kinetic_status::Symbol
    final_electron_nuclear_by_center
    final_electron_nuclear_status::Symbol
    final_hamiltonian
    final_hamiltonian_status::Symbol
    h1
    h1_status::Symbol
    summary
    metadata
end

struct _PQSDiatomicPhysicalGaussletH1JPayload
    status::Symbol
    blocker
    route_family::Symbol
    source_plan
    source_plan_status::Symbol
    final_basis
    final_basis_status::Symbol
    h1_payload
    h1_payload_status::Symbol
    density_provenance
    density_provenance_status::Symbol
    support_weights
    support_weights_status::Symbol
    raw_pair_factor_terms
    raw_pair_factor_status::Symbol
    support_pair_raw_numerator
    support_pair_raw_numerator_status::Symbol
    density_interaction
    density_interaction_status::Symbol
    h1_j_diagnostic
    h1_j_status::Symbol
    summary
    metadata
end

struct _PQSDiatomicPhysicalGaussletRHFInputContractPayload
    status::Symbol
    blocker
    route_family::Symbol
    source_plan
    source_plan_status::Symbol
    final_basis
    final_basis_status::Symbol
    h1_payload
    h1_payload_status::Symbol
    h1_j_payload
    h1_j_payload_status::Symbol
    electron_count
    occupation
    summary
    metadata
end

struct _PQSDiatomicPhysicalGaussletRHFExecutionPayload
    status::Symbol
    blocker
    route_family::Symbol
    input_contract
    input_contract_status::Symbol
    h1_payload
    h1_payload_status::Symbol
    h1_j_payload
    h1_j_payload_status::Symbol
    summary
    metadata
end

struct _PQSDiatomicPhysicalGaussletSourcePlanCandidatePayload
    status::Symbol
    blocker
    candidate_source::Symbol
    candidate
    counts_match::Bool
    authority_status::Symbol
    summary
    metadata
end

struct _PQSDiatomicPhysicalGaussletCoreShellSourcePlan
    object_kind::Symbol
    status::Symbol
    parent_basis
    axis_bundles
    atom_contact_core_support_indices::Vector{Int}
    atom_contact_core_support_states
    shared_shell_support_indices::Tuple
    shared_shell_support_states::Tuple
    core_coefficient_matrix
    shared_shell_coefficient_matrices::Tuple
    support_order::Tuple
    retained_order::Tuple
    retained_ranges
    final_dimension::Int
    convention_labels
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
    expected_route_kinds = (
        :bond_aligned_diatomic_physical_gausslet_core_shell_pqs,
        :bond_aligned_diatomic_independent_pqs_source_box_core_shell,
    )

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
        generated_support_plan.status in (
            :available_independent_pqs_support_region_plan,
            :available_terminal_topology_support_region_plan,
        )
    retained_rule_plan =
        independent_target ?
        _pqs_source_box_route_driver_independent_h2_retained_rule_plan(
            generated_support_plan,
        ) :
        nothing
    if independent_target && generated_support_available
        status = :available_physical_gausslet_core_shell_target_inventory
        blocker = nothing
    end
    support_units =
        generated_support_available ? generated_support_plan.support_order :
        support_units
    support_counts =
        generated_support_available ? generated_support_plan.support_counts :
        isnothing(inventory) ? (;) : inventory.support_counts
    retained_units =
        !isnothing(retained_rule_plan) ? retained_rule_plan.retained_order :
        retained_units
    retained_counts =
        !isnothing(retained_rule_plan) ? retained_rule_plan.retained_counts :
        isnothing(inventory) ? (;) : inventory.retained_counts
    retained_order =
        !isnothing(retained_rule_plan) ? retained_rule_plan.retained_order :
        isnothing(inventory) ? () : Tuple(inventory.retained_order)
    expected_final_dimension =
        !isnothing(retained_rule_plan) ? retained_rule_plan.expected_final_dimension :
        isnothing(inventory) ? nothing : inventory.expected_final_dimension
    retained_atom_core_interiors =
        !isnothing(retained_rule_plan) ||
        (!isnothing(inventory) && inventory.retained_atom_core_interiors)
    source_backed_fixed_source_oracle_used =
        false
    source_plan_blocker =
        !isnothing(retained_rule_plan) ?
        nothing :
        generated_support_available ?
        :missing_independent_pqs_retained_rule_plan :
        nothing
    retained_transform_authority =
        !isnothing(retained_rule_plan) ?
        :pqs_source_box_construction :
        generated_support_available ?
        :missing_terminal_retained_rule_plan :
        isnothing(inventory) ?
        :not_available :
        :pqs_source_box_construction
    primary_blocker =
        !isnothing(retained_rule_plan) ?
        source_plan_blocker :
        generated_support_available ?
        source_plan_blocker :
        blocker
    secondary_blocker =
        !isnothing(retained_rule_plan) ?
        nothing :
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
        source_backed_fixed_source_oracle_used,
        retained_transform_authority,
        primary_blocker,
        secondary_blocker,
        independent_source_plan_blocker = source_plan_blocker,
        support_plan,
        retained_rule_plan,
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
    fixture_label = :h2_r4_physical_gausslet_q5
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
        fixture_label = :h2_r4_physical_gausslet_q5
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
                    if err isa ArgumentError
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
    fixture_label = :h2_r4_physical_gausslet_q5
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

function _pqs_source_box_route_driver_h2_wl_gausslet_only_reference_candidate(
    parent,
    route_skeleton,
    recipe,
    target_payload,
    final_basis_payload = nothing,
)
    route_family =
        hasproperty(route_skeleton, :route_family) ?
        route_skeleton.route_family :
        recipe.route_family
    route_kind =
        hasproperty(route_skeleton, :route_kind) ?
        route_skeleton.route_kind :
        recipe.route_kind
    final_basis_summary =
        isnothing(final_basis_payload) || !hasproperty(final_basis_payload, :summary) ?
        (;) :
        final_basis_payload.summary
    actual_final_dimension = get(final_basis_summary, :final_dimension, nothing)
    reference_label = "WL/QW H2 R=4 gausslet-only 463"
    retained_transform_kind = :white_lindsey_old_qw_gausslet_retained_transform
    expected_order = (:atom_contact_core, :shared_shell_1, :shared_shell_2)
    support_counts =
        _pqs_source_box_route_driver_ordered_count_tuple(
            target_payload.support_counts,
            expected_order,
        )
    retained_counts =
        _pqs_source_box_route_driver_ordered_count_tuple(
            target_payload.retained_counts,
            expected_order,
        )
    final_dimension =
        isnothing(actual_final_dimension) ?
        target_payload.expected_final_dimension :
        actual_final_dimension
    conditions = (;
        route_kind =
            route_family === :pqs_source_box &&
            route_kind === :bond_aligned_diatomic_physical_gausslet_core_shell_pqs,
        geometry = _pqs_source_box_route_driver_h2_r4_center_match(parent),
        parent_axis_counts =
            _pqs_source_box_route_driver_axis_counts_tuple(
                target_payload.parent_axis_counts,
            ) == (9, 9, 15),
        support_counts = support_counts == (275, 578, 362),
        retained_counts = retained_counts == (251, 98, 114),
        retained_order = target_payload.retained_order == expected_order,
        supplement_policy = target_payload.supplement_policy === :none,
        final_dimension =
            target_payload.expected_final_dimension == 463 &&
            (isnothing(actual_final_dimension) || actual_final_dimension == 463),
        retained_transform_kind =
            retained_transform_kind ===
            :white_lindsey_old_qw_gausslet_retained_transform,
        reference_label = reference_label == "WL/QW H2 R=4 gausslet-only 463",
    )
    mismatches = Tuple(key for key in keys(conditions) if !getproperty(conditions, key))
    status =
        isempty(mismatches) ?
        :available_wl_h2_gausslet_only_reference_candidate :
        :blocked_wl_h2_gausslet_only_reference_candidate
    blocker =
        isempty(mismatches) ?
        nothing :
        :wl_h2_gausslet_only_reference_candidate_mismatch
    summary = (;
        status,
        blocker,
        route_family,
        route_kind,
        geometry_label = :h2_r4,
        parent_axis_counts =
            _pqs_source_box_route_driver_axis_counts_tuple(
                target_payload.parent_axis_counts,
            ),
        support_counts,
        retained_counts,
        retained_order = target_payload.retained_order,
        supplement_policy = target_payload.supplement_policy,
        final_dimension,
        retained_transform_kind,
        label = reference_label,
        old_supplemented_scalar_references_blocked = true,
        conditions,
        mismatches,
    )
    metadata = (;
        source =
            :pqs_source_box_route_driver_h2_wl_gausslet_only_reference_candidate,
        old_wl_qw_fixed_block_size = (1215, 463),
        reviewed_contract_pass = 218,
        matrix_materialized = false,
        scalar_reference_materialized = false,
    )
    return _PQSH2WLGaussletOnlyReferenceCandidatePayload(
        status,
        blocker,
        route_family,
        route_kind,
        reference_label,
        retained_transform_kind,
        target_payload.supplement_policy,
        final_dimension,
        conditions,
        summary,
        metadata,
    )
end

function _pqs_source_box_route_driver_diatomic_physical_gausslet_source_plan_candidate_payload(
    parent,
    route_skeleton,
    _recipe,
    target_payload,
)
    expected_order = (:atom_contact_core, :shared_shell_1, :shared_shell_2)
    expected_support_counts = (275, 578, 362)
    expected_retained_counts = (251, 98, 114)
    expected_final_dimension = 463
    target_available =
        !isnothing(target_payload) &&
        target_payload.status === :available_physical_gausslet_core_shell_target_inventory
    parent_basis =
        hasproperty(parent, :parent_qw_basis_object) ?
        parent.parent_qw_basis_object :
        nothing
    if !target_available || isnothing(parent_basis)
        blocker =
            target_available ?
            :missing_physical_gausslet_source_plan_materializer :
            :missing_physical_gausslet_target_inventory
        summary = (;
            candidate_status = :not_available,
            candidate_source = :not_available,
            candidate_counts_match = false,
            source_plan_authority_status = :not_available,
            blocker,
        )
        return _PQSDiatomicPhysicalGaussletSourcePlanCandidatePayload(
            :not_available_physical_gausslet_source_plan_candidate,
            blocker,
            :not_available,
            nothing,
            false,
            :not_available,
            summary,
            (; source = :pqs_diatomic_physical_gausslet_source_plan_candidate),
        )
    end

    source = bond_aligned_diatomic_nested_fixed_source(parent_basis; nside = 5)
    atom_support_count = sum(
        child -> length(child.support_indices),
        source.child_sequences;
        init = 0,
    )
    if !isnothing(source.geometry.shared_midpoint_box)
        dims = Tuple(length.(source.geometry.parent_box))
        atom_support_count += length(
            _nested_box_support_indices(source.geometry.shared_midpoint_box..., dims),
        )
    end
    support_counts = (
        atom_support_count,
        (length(layer.support_indices) for layer in source.shared_shell_layers)...,
    )
    retained_counts = (
        length(source.sequence.core_column_range),
        (length(range) for range in source.sequence.layer_column_ranges)...,
    )
    final_dimension = size(source.sequence.coefficient_matrix, 2)
    counts_match =
        support_counts == expected_support_counts &&
        retained_counts == expected_retained_counts &&
        final_dimension == expected_final_dimension
    retained_order = expected_order
    order_match =
        target_payload.support_units == expected_order &&
        target_payload.retained_order == retained_order
    no_supplement = target_payload.supplement_policy === :none
    status =
        counts_match && order_match && no_supplement ?
        :available_physical_gausslet_source_plan_candidate :
        :blocked_physical_gausslet_source_plan_candidate
    blocker =
        status === :available_physical_gausslet_source_plan_candidate ?
        :source_plan_candidate_not_route_authority :
        :physical_gausslet_source_plan_count_mismatch
    authority_status = :candidate_not_route_authority
    summary = (;
        candidate_status = status,
        candidate_source = :source_backed_fixed_source_oracle,
        candidate_counts_match = counts_match,
        source_plan_authority_status = authority_status,
        blocker,
        support_order = expected_order,
        retained_order,
        support_counts,
        retained_counts,
        final_dimension,
        no_supplement,
    )
    return _PQSDiatomicPhysicalGaussletSourcePlanCandidatePayload(
        status,
        blocker,
        :source_backed_fixed_source_oracle,
        source,
        counts_match,
        authority_status,
        summary,
        (;
            source =
                :pqs_diatomic_physical_gausslet_source_plan_candidate,
            route_owned_authority = false,
            diagnostic_221_source_plan_reused = false,
        ),
    )
end

function _pqs_source_box_route_driver_physical_gausslet_source_plan_from_candidate(
    candidate_payload,
)
    source = candidate_payload.candidate
    shared_shell_support_indices =
        Tuple(Vector{Int}(layer.support_indices) for layer in source.shared_shell_layers)
    shared_shell_support_states =
        Tuple(layer.support_states for layer in source.shared_shell_layers)
    retained_ranges = (
        atom_contact_core = source.sequence.core_column_range,
        shared_shell_1 = source.sequence.layer_column_ranges[1],
        shared_shell_2 = source.sequence.layer_column_ranges[2],
    )
    core_coefficient_matrix =
        source.sequence.coefficient_matrix[:, retained_ranges.atom_contact_core]
    shared_shell_coefficient_matrices = (
        source.sequence.coefficient_matrix[:, retained_ranges.shared_shell_1],
        source.sequence.coefficient_matrix[:, retained_ranges.shared_shell_2],
    )
    support_order = (:atom_contact_core, :shared_shell_1, :shared_shell_2)
    retained_order = support_order
    support_counts = (
        length(source.sequence.core_indices),
        length(shared_shell_support_indices[1]),
        length(shared_shell_support_indices[2]),
    )
    retained_counts = (
        length(retained_ranges.atom_contact_core),
        length(retained_ranges.shared_shell_1),
        length(retained_ranges.shared_shell_2),
    )
    final_dimension = size(source.sequence.coefficient_matrix, 2)
    convention_labels = (;
        source_plan_family = :fake_pqs_source_backed_fixed_source_adapter,
        source_backed_adapter = true,
        source_backed_candidate_source = candidate_payload.candidate_source,
        retained_transform_kind = :wl_qw_source_backed_retained_transform,
        independent_pqs_transform = false,
        fake_pqs_source = :source_backed_fixed_source_oracle,
        fake_pqs_warning =
            :retained_transform_imported_from_wl_qw_fixed_source_oracle,
        route_owned_authority = true,
        supplement_policy = :none,
        h2_221_diagnostic_source_plan_reused = false,
        final_basis_materialized = false,
        h1_materialized = false,
        h1_j_materialized = false,
        rhf_materialized = false,
        exports_materialized = false,
    )
    summary = (;
        object_kind = :pqs_diatomic_physical_gausslet_core_shell_source_plan,
        status = :available_pqs_diatomic_physical_gausslet_core_shell_source_plan,
        support_order,
        retained_order,
        support_counts,
        retained_counts,
        final_dimension,
        retained_ranges,
        source_plan_authority_status = :fake_pqs_private_source_backed_adapter_authority,
        source_backed_candidate_source = candidate_payload.candidate_source,
        source_backed_adapter = true,
        retained_transform_kind = :wl_qw_source_backed_retained_transform,
        independent_pqs_transform = false,
        route_owned_authority = true,
        supplement_policy = :none,
    )
    return _PQSDiatomicPhysicalGaussletCoreShellSourcePlan(
        :pqs_diatomic_physical_gausslet_core_shell_source_plan,
        :available_pqs_diatomic_physical_gausslet_core_shell_source_plan,
        source.basis,
        source.axis_bundles,
        Vector{Int}(source.sequence.core_indices),
        source.sequence.core_states,
        shared_shell_support_indices,
        shared_shell_support_states,
        core_coefficient_matrix,
        shared_shell_coefficient_matrices,
        support_order,
        retained_order,
        retained_ranges,
        final_dimension,
        convention_labels,
        summary,
        (;
            source =
                :pqs_source_box_route_driver_physical_gausslet_source_plan_from_candidate,
            source_backed_adapter = true,
        ),
    )
end

function _pqs_source_box_route_driver_independent_h2_physical_source_plan_descriptor(
    target_payload,
)
    descriptor_blocker = :missing_independent_pqs_source_plan_numerical_materialization
    blocked(blocker) = (;
        object_kind = :pqs_independent_h2_physical_source_plan_descriptor,
        status = :blocked_independent_pqs_physical_source_plan_descriptor,
        blocker,
        source_plan_family = :independent_pqs_physical_source_box_descriptor,
        source_plan_authority_status =
            :blocked_independent_pqs_route_owned_source_plan_descriptor,
        source_coefficients_materialized = false,
        final_basis_materialized = false,
        source_backed_fixed_source_oracle_used = false,
        fake_pqs_enabled = false,
        missing_objects = (blocker,),
    )
    isnothing(target_payload) &&
        return blocked(:missing_independent_pqs_target_payload)
    target_payload.route_kind ===
        :bond_aligned_diatomic_independent_pqs_source_box_core_shell ||
        return blocked(:not_independent_pqs_route_kind)
    target_payload.status === :available_physical_gausslet_core_shell_target_inventory ||
        return blocked(target_payload.blocker)
    support_plan = get(target_payload.summary, :support_plan, nothing)
    retained_rule_plan = get(target_payload.summary, :retained_rule_plan, nothing)
    isnothing(support_plan) &&
        return blocked(:missing_independent_pqs_support_region_plan)
    isnothing(retained_rule_plan) &&
        return blocked(:missing_independent_pqs_retained_rule_plan)
    support_counts = target_payload.support_counts
    retained_counts = target_payload.retained_counts
    Tuple(values(support_counts)) == (275, 578, 362) ||
        return blocked(:independent_pqs_source_plan_support_count_mismatch)
    Tuple(values(retained_counts)) == (275, 98, 98) ||
        return blocked(:independent_pqs_source_plan_retained_count_mismatch)

    shared_source(role) = (;
        source_family = :pqs_filled_source_cpb,
        source_mode_dims = retained_rule_plan.shared_shell_source_mode_dims,
        source_mode_count = prod(retained_rule_plan.shared_shell_source_mode_dims),
        support_count = getproperty(support_counts, role),
        transform_kind = getproperty(retained_rule_plan.per_unit_transform_kind, role),
        coefficient_matrix_materialized = false,
    )
    retained_rule(role, rule, materialized) = (;
        retained_rule = rule,
        retained_count = getproperty(retained_counts, role),
        retained_column_selector_materialized = materialized,
    )
    descriptor_summary = (;
        source_plan_family = :independent_pqs_physical_source_box_descriptor,
        source_plan_authority_status =
            :independent_pqs_route_owned_source_plan_descriptor,
        support_counts,
        retained_counts,
        final_dimension = target_payload.expected_final_dimension,
        source_coefficients_materialized = false,
        final_basis_materialized = false,
    )
    return (;
        object_kind = :pqs_independent_h2_physical_source_plan_descriptor,
        status = :available_independent_pqs_physical_source_plan_descriptor,
        blocker = nothing,
        source_plan_family = :independent_pqs_physical_source_box_descriptor,
        source_plan_authority_status =
            :independent_pqs_route_owned_source_plan_descriptor,
        support_order = target_payload.support_units,
        retained_order = target_payload.retained_order,
        support_counts,
        retained_counts,
        final_dimension = target_payload.expected_final_dimension,
        expected_final_dimension = target_payload.expected_final_dimension,
        unit_source_descriptors = (;
            atom_contact_core = (;
                source_family = :direct_terminal_source_modes,
                source_region_roles =
                    support_plan.atom_contact_core_descriptor.source_region_roles,
                support_count = support_counts.atom_contact_core,
                source_mode_count = support_counts.atom_contact_core,
                transform_kind = :identity_source_modes,
                coefficient_matrix_materialized = false,
            ),
            shared_shell_1 = shared_source(:shared_shell_1),
            shared_shell_2 = shared_source(:shared_shell_2),
        ),
        unit_retained_rule_descriptors = (;
            atom_contact_core =
                retained_rule(:atom_contact_core, :direct_source_modes, false),
            shared_shell_1 = retained_rule(
                :shared_shell_1,
                retained_rule_plan.shared_shell_retained_rule_kind,
                true,
            ),
            shared_shell_2 = retained_rule(
                :shared_shell_2,
                retained_rule_plan.shared_shell_retained_rule_kind,
                true,
            ),
        ),
        source_coefficients_materialized = false,
        final_basis_materialized = false,
        source_backed_fixed_source_oracle_used = false,
        fake_pqs_enabled = false,
        next_blocker = descriptor_blocker,
        support_plan,
        retained_rule_plan,
        missing_objects = (descriptor_blocker,),
        summary = descriptor_summary,
        metadata = (;
            source =
                :pqs_source_box_route_driver_independent_h2_physical_source_plan_descriptor,
            descriptor_only = true,
        ),
    )
end

function _pqs_source_box_route_driver_independent_h2_shared_shell_realization_payload(
    parent,
    target_payload,
    descriptor,
)
    blocked(blocker) = (;
        object_kind = :pqs_independent_h2_shared_shell_realization_payload,
        status = :blocked_independent_pqs_shared_shell_realization_payload,
        blocker,
        shared_shell_realization_counts = (),
        shared_shell_realization_identity_errors = (),
        shared_shell_realization_materialized = false,
        source_backed_fixed_source_oracle_used = false,
        fake_pqs_enabled = false,
        shared_shells = (;),
        summary = (;
            status = :blocked_independent_pqs_shared_shell_realization_payload,
            blocker,
            shared_shell_realization_counts = (),
            shared_shell_realization_identity_errors = (),
            shared_shell_realization_materialized = false,
        ),
    )
    isnothing(descriptor) &&
        return blocked(:missing_independent_pqs_source_plan_descriptor)
    descriptor.status === :available_independent_pqs_physical_source_plan_descriptor ||
        return blocked(descriptor.blocker)
    support_plan = get(target_payload.summary, :support_plan, nothing)
    retained_rule_plan = get(target_payload.summary, :retained_rule_plan, nothing)
    isnothing(support_plan) &&
        return blocked(:missing_independent_pqs_support_region_plan)
    isnothing(retained_rule_plan) &&
        return blocked(:missing_independent_pqs_retained_rule_plan)
    hasproperty(support_plan, :shared_shell_descriptors) ||
        return blocked(:missing_independent_pqs_shared_shell_support_descriptors)
    bundles =
        hasproperty(parent, :parent_axis_bundle_object) ?
        parent.parent_axis_bundle_object :
        nothing
    isnothing(bundles) && return blocked(:missing_parent_axis_bundle_object)
    metrics = _pqs_multilayer_axis_metrics(bundles)
    source_mode_dims = retained_rule_plan.shared_shell_source_mode_dims
    shells = map(
        role -> _pqs_source_box_route_driver_independent_h2_shared_shell_realization(
            role,
            getproperty(support_plan.shared_shell_descriptors, role),
            bundles,
            metrics,
            source_mode_dims,
            support_plan.bond_axis,
        ),
        (:shared_shell_1, :shared_shell_2),
    )
    blocker = findfirst(shell -> shell.status !==
                                 :available_independent_pqs_shared_shell_realization,
                        shells)
    isnothing(blocker) || return blocked(shells[blocker].blocker)
    counts = Tuple(shell.retained_count for shell in shells)
    identity_errors = Tuple(shell.realized_overlap_identity_error for shell in shells)
    summary = (;
        status = :available_independent_pqs_shared_shell_realization_payload,
        blocker = nothing,
        shared_shell_realization_counts = counts,
        shared_shell_realization_identity_errors = identity_errors,
        coefficient_shapes = Tuple(shell.coefficient_shape for shell in shells),
        shell_projection_materialized = true,
        lowdin_cleanup_materialized = true,
        shared_shell_realization_materialized = true,
        source_backed_fixed_source_oracle_used = false,
        fake_pqs_enabled = false,
    )
    return (;
        object_kind = :pqs_independent_h2_shared_shell_realization_payload,
        status = summary.status,
        blocker = nothing,
        shared_shell_realization_counts = counts,
        shared_shell_realization_identity_errors = identity_errors,
        shared_shell_realization_materialized = true,
        source_backed_fixed_source_oracle_used = false,
        fake_pqs_enabled = false,
        shared_shells = (; shared_shell_1 = shells[1], shared_shell_2 = shells[2]),
        summary,
    )
end

function _pqs_source_box_route_driver_independent_h2_shared_shell_realization(
    role::Symbol,
    shell_descriptor,
    bundles,
    metrics,
    source_mode_dims::NTuple{3,Int},
    bond_axis::Symbol,
)
    blocked(blocker) = (;
        role,
        status = :blocked_independent_pqs_shared_shell_realization,
        blocker,
        retained_count = nothing,
        coefficient_shape = nothing,
        realized_overlap_identity_error = nothing,
    )
    try
        source_q = source_mode_dims[1]
        source_key = Symbol("independent_h2_", role, "_q", source_q, "_source")
        raw_plan = CartesianRawProductSources.raw_product_box_plan(
            shell_descriptor.source_cpb;
            source_mode_dims,
            source_key,
            metadata = (; route_unit_role = role),
        )
        retained_rule =
            CartesianRawProductSources.pqs_boundary_product_mode_retained_rule(
                raw_plan;
                metadata = (; route_unit_role = role),
            )
        retained_count = Int(retained_rule.retained_count)
        retained_count > 0 ||
            return blocked(:independent_pqs_shared_shell_retained_count_mismatch)
        layer = _nested_projected_q_shell_layer(
            bundles,
            shell_descriptor.current_box,
            shell_descriptor.inner_box;
            bond_axis,
            q = source_mode_dims[1],
            L = source_mode_dims[1],
            raw_source_dims = source_mode_dims,
            selected_q = source_mode_dims[1],
        )
        projected_descriptor = _nested_projected_q_shell_staged_unit_descriptor(layer)
        projected_descriptor.support_count == shell_descriptor.support_count ||
            return blocked(:independent_pqs_shared_shell_support_count_mismatch)
        shell_plan =
            CartesianContractedParentMetrics._pqs_shell_realization_plan(
                projected_descriptor,
                metrics,
            )
        coefficients = shell_plan.shell_projection_matrix * shell_plan.lowdin_cleanup
        size(coefficients) == (shell_descriptor.support_count, retained_count) ||
            return blocked(:independent_pqs_shared_shell_coefficient_shape_mismatch)
        return (;
            role,
            status = :available_independent_pqs_shared_shell_realization,
            blocker = nothing,
            raw_source_plan = raw_plan,
            retained_rule,
            support_indices = shell_descriptor.support_indices,
            support_states = shell_descriptor.support_states,
            shell_projection = shell_plan.shell_projection_matrix,
            lowdin_cleanup = shell_plan.lowdin_cleanup,
            shell_final_coefficients = coefficients,
            retained_count,
            coefficient_shape = size(coefficients),
            realized_overlap_identity_error = shell_plan.isometry_error,
            shell_projection_materialized = true,
            lowdin_cleanup_materialized = true,
            source_backed_fixed_source_oracle_used = false,
            fake_pqs_enabled = false,
        )
    catch err
        return blocked(Symbol(:independent_pqs_shared_shell_realization_error))
    end
end

function _pqs_source_box_route_driver_parent_row_sparse_coefficients(
    support_indices::AbstractVector{Int},
    local_coefficients::AbstractMatrix{<:Real},
    parent_count::Int,
)
    row_indices = Int[]
    col_indices = Int[]
    values = Float64[]
    for column in axes(local_coefficients, 2)
        for (local_row, parent_row) in pairs(support_indices)
            value = local_coefficients[local_row, column]
            iszero(value) && continue
            push!(row_indices, Int(parent_row))
            push!(col_indices, Int(column))
            push!(values, Float64(value))
        end
    end
    return _nested_sparse_coefficient_map(
        row_indices,
        col_indices,
        values,
        parent_count,
        size(local_coefficients, 2),
    )
end

function _pqs_source_box_route_driver_independent_h2_complete_core_shell_source_plan(
    parent,
    target_payload,
    descriptor,
    shared_shell_realization,
)
    blocked(blocker) = (;
        status = :blocked_pqs_diatomic_physical_gausslet_core_shell_source_plan,
        blocker,
        source_plan = nothing,
    )
    isnothing(target_payload) && return blocked(:missing_independent_pqs_target_payload)
    descriptor.status === :available_independent_pqs_physical_source_plan_descriptor ||
        return blocked(descriptor.blocker)
    shared_shell_realization.status ===
        :available_independent_pqs_shared_shell_realization_payload ||
        return blocked(shared_shell_realization.blocker)
    support_plan = get(target_payload.summary, :support_plan, nothing)
    isnothing(support_plan) &&
        return blocked(:missing_independent_pqs_support_region_plan)
    hasproperty(support_plan, :atom_contact_core_descriptor) ||
        return blocked(:missing_independent_pqs_atom_contact_core_descriptor)
    bundles =
        hasproperty(parent, :parent_axis_bundle_object) ?
        parent.parent_axis_bundle_object :
        nothing
    isnothing(bundles) && return blocked(:missing_parent_axis_bundle_object)

    support_order = (:atom_contact_core, :shared_shell_1, :shared_shell_2)
    retained_order = support_order
    support_counts =
        Tuple(getproperty(target_payload.support_counts, role) for role in support_order)
    retained_counts =
        Tuple(getproperty(target_payload.retained_counts, role) for role in retained_order)
    target_payload.support_units == support_order ||
        return blocked(:independent_pqs_source_plan_support_order_mismatch)
    target_payload.retained_order == retained_order ||
        return blocked(:independent_pqs_source_plan_retained_order_mismatch)
    Tuple(values(target_payload.support_counts)) == support_counts ||
        return blocked(:independent_pqs_source_plan_support_count_mismatch)
    Tuple(values(target_payload.retained_counts)) == retained_counts ||
        return blocked(:independent_pqs_source_plan_retained_count_mismatch)

    core = support_plan.atom_contact_core_descriptor
    core.support_count == retained_counts[1] ||
        return blocked(:independent_pqs_atom_contact_core_count_mismatch)
    shell_1 = shared_shell_realization.shared_shells.shared_shell_1
    shell_2 = shared_shell_realization.shared_shells.shared_shell_2
    parent_dims = _pqs_source_box_route_driver_axis_counts_tuple(
        target_payload.parent_axis_counts,
    )
    isnothing(parent_dims) &&
        return blocked(:missing_independent_pqs_parent_axis_counts)
    parent_count = prod(parent_dims)
    core_indices = Vector{Int}(core.support_indices)
    core_coefficient_matrix = _nested_sparse_coefficient_map(
        core_indices,
        collect(1:length(core_indices)),
        ones(Float64, length(core_indices)),
        parent_count,
        length(core_indices),
    )
    shared_shell_coefficient_matrices = (
        _pqs_source_box_route_driver_parent_row_sparse_coefficients(
            shell_1.support_indices,
            shell_1.shell_final_coefficients,
            parent_count,
        ),
        _pqs_source_box_route_driver_parent_row_sparse_coefficients(
            shell_2.support_indices,
            shell_2.shell_final_coefficients,
            parent_count,
        ),
    )
    retained_ranges = (;
        atom_contact_core = 1:retained_counts[1],
        shared_shell_1 = (retained_counts[1] + 1):sum(retained_counts[1:2]),
        shared_shell_2 = (sum(retained_counts[1:2]) + 1):sum(retained_counts),
    )
    convention_labels = (;
        source_plan_family = :independent_pqs_physical_source_box_core_shell,
        source_backed_adapter = false,
        source_backed_fixed_source_oracle_used = false,
        retained_transform_kind = :pqs_source_box_retained_transform,
        independent_pqs_transform = true,
        fake_pqs_enabled = false,
        retained_transform_authority = :pqs_source_box_construction,
        source_plan_authority_status = :independent_pqs_route_owned_source_plan,
        route_owned_authority = true,
        core_coefficient_representation = :sparse_parent_row_direct_selector,
        shared_shell_coefficient_representation = :sparse_parent_row_realization,
        final_basis_materialized = false,
        h1_materialized = false,
        h1_j_materialized = false,
        rhf_materialized = false,
        exports_materialized = false,
    )
    summary = (;
        object_kind = :pqs_diatomic_physical_gausslet_core_shell_source_plan,
        status = :available_pqs_diatomic_physical_gausslet_core_shell_source_plan,
        blocker = nothing,
        support_order,
        retained_order,
        support_counts,
        retained_counts,
        final_dimension = sum(retained_counts),
        retained_ranges,
        source_plan_family = :independent_pqs_physical_source_box_core_shell,
        source_plan_authority_status = :independent_pqs_route_owned_source_plan,
        source_backed_fixed_source_oracle_used = false,
        fake_pqs_enabled = false,
        retained_transform_authority = :pqs_source_box_construction,
        source_coefficients_materialized = true,
        final_basis_materialized = false,
        shared_shell_realization_counts =
            shared_shell_realization.shared_shell_realization_counts,
        shared_shell_realization_identity_errors =
            shared_shell_realization.shared_shell_realization_identity_errors,
    )
    parent_basis =
        hasproperty(parent, :parent_qw_basis_object) ? parent.parent_qw_basis_object :
        hasproperty(parent, :parent_basis_object) ? parent.parent_basis_object :
        nothing
    source_plan = _PQSDiatomicPhysicalGaussletCoreShellSourcePlan(
        :pqs_diatomic_physical_gausslet_core_shell_source_plan,
        :available_pqs_diatomic_physical_gausslet_core_shell_source_plan,
        parent_basis,
        bundles,
        core_indices,
        core.support_states,
        (Vector{Int}(shell_1.support_indices), Vector{Int}(shell_2.support_indices)),
        (shell_1.support_states, shell_2.support_states),
        core_coefficient_matrix,
        shared_shell_coefficient_matrices,
        support_order,
        retained_order,
        retained_ranges,
        sum(retained_counts),
        convention_labels,
        summary,
        (;
            source =
                :pqs_source_box_route_driver_independent_h2_complete_core_shell_source_plan,
            route_owned = true,
            source_backed_adapter = false,
        ),
    )
    return (;
        status = source_plan.status,
        blocker = nothing,
        source_plan,
    )
end

function _pqs_source_box_route_driver_diatomic_physical_gausslet_source_plan_payload(
    parent,
    target_payload,
    candidate_payload = nothing,
)
    target_available =
        !isnothing(target_payload) &&
        target_payload.status === :available_physical_gausslet_core_shell_target_inventory
    independent_target =
        !isnothing(target_payload) &&
        target_payload.route_kind ===
        :bond_aligned_diatomic_independent_pqs_source_box_core_shell
    source_plan =
        !isnothing(candidate_payload) &&
        candidate_payload.status === :available_physical_gausslet_source_plan_candidate ?
        _pqs_source_box_route_driver_physical_gausslet_source_plan_from_candidate(
            candidate_payload,
        ) :
        nothing
    descriptor =
        independent_target ?
        _pqs_source_box_route_driver_independent_h2_physical_source_plan_descriptor(
            target_payload,
        ) :
        nothing
    descriptor_available =
        !isnothing(descriptor) &&
        descriptor.status === :available_independent_pqs_physical_source_plan_descriptor
    shared_shell_realization =
        descriptor_available ?
        _pqs_source_box_route_driver_independent_h2_shared_shell_realization_payload(
            parent,
            target_payload,
            descriptor,
        ) :
        nothing
    shared_shell_realization_available =
        !isnothing(shared_shell_realization) &&
        shared_shell_realization.status ===
        :available_independent_pqs_shared_shell_realization_payload
    independent_source_plan =
        shared_shell_realization_available ?
        _pqs_source_box_route_driver_independent_h2_complete_core_shell_source_plan(
            parent,
            target_payload,
            descriptor,
            shared_shell_realization,
        ) :
        nothing
    independent_source_plan_available =
        !isnothing(independent_source_plan) &&
        independent_source_plan.status ===
        :available_pqs_diatomic_physical_gausslet_core_shell_source_plan
    source_plan =
        independent_source_plan_available ?
        independent_source_plan.source_plan :
        descriptor_available ? nothing : source_plan
    independent_source_plan_blocker =
        !isnothing(target_payload) ?
        get(
            target_payload.summary,
            :independent_source_plan_blocker,
            :missing_independent_pqs_physical_source_plan_materializer,
        ) :
        :missing_independent_pqs_physical_source_plan_materializer
    status =
        independent_source_plan_available ?
        source_plan.status :
        descriptor_available ?
        :blocked_pqs_diatomic_physical_gausslet_core_shell_source_plan :
        isnothing(source_plan) ?
        :blocked_pqs_diatomic_physical_gausslet_core_shell_source_plan :
        source_plan.status
    blocker =
        independent_source_plan_available ?
        nothing :
        descriptor_available ?
        shared_shell_realization_available ?
        independent_source_plan.blocker :
        shared_shell_realization.blocker :
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
        source_plan_descriptor_status =
            isnothing(descriptor) ? :not_available : descriptor.status,
        source_plan_descriptor_blocker =
            isnothing(descriptor) ? nothing : descriptor.blocker,
        shared_shell_realization_status =
            isnothing(shared_shell_realization) ?
            :not_available :
            shared_shell_realization.status,
        shared_shell_realization_blocker =
            isnothing(shared_shell_realization) ?
            nothing :
            shared_shell_realization.blocker,
        shared_shell_realization_counts =
            isnothing(shared_shell_realization) ?
            () :
            shared_shell_realization.shared_shell_realization_counts,
        shared_shell_realization_identity_errors =
            isnothing(shared_shell_realization) ?
            () :
            shared_shell_realization.shared_shell_realization_identity_errors,
        shared_shell_realization_materialized = shared_shell_realization_available,
        source_plan_descriptor_available = descriptor_available,
        source_plan_family =
            independent_source_plan_available ?
            source_plan.summary.source_plan_family :
            isnothing(source_plan) ?
            :not_available :
            hasproperty(source_plan, :source_plan_family) ?
            source_plan.source_plan_family :
            :not_available,
        source_plan_authority_status =
            independent_source_plan_available ?
            source_plan.summary.source_plan_authority_status :
            descriptor_available ?
            descriptor.source_plan_authority_status :
            isnothing(source_plan) ?
            independent_target ?
            :blocked_pqs_source_box_construction_authority :
            isnothing(candidate_payload) ? :not_available : candidate_payload.authority_status :
            :fake_pqs_private_source_backed_adapter_authority,
        source_backed_fixed_source_oracle_used =
            !isnothing(target_payload) &&
            get(target_payload.summary, :source_backed_fixed_source_oracle_used, false),
        retained_transform_authority =
            independent_source_plan_available ?
            source_plan.summary.retained_transform_authority :
            isnothing(target_payload) ?
            :not_available :
            get(target_payload.summary, :retained_transform_authority, :not_available),
        source_coefficients_materialized =
            independent_source_plan_available,
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

function _pqs_source_box_route_driver_support_partition_rows(
    source_plan,
    unit_key::Symbol,
)
    unit_key === :atom_contact_core &&
        return Vector{Int}(source_plan.atom_contact_core_support_indices)
    unit_key === :shared_shell_1 &&
        return Vector{Int}(source_plan.shared_shell_support_indices[1])
    unit_key === :shared_shell_2 &&
        return Vector{Int}(source_plan.shared_shell_support_indices[2])
    return Int[]
end

function _pqs_source_box_route_driver_support_partition_local_rows(
    support_indices::AbstractVector{Int},
    tile_indices::AbstractVector{Int},
)
    local_by_parent = Dict{Int,Int}(
        Int(parent_row) => Int(local_row) for (local_row, parent_row) in pairs(support_indices)
    )
    return Tuple(local_by_parent[Int(parent_row)] for parent_row in tile_indices)
end

function _pqs_source_box_route_driver_contiguous_unit_range(rows::Tuple)
    isempty(rows) && return nothing
    first_row = first(rows)
    expected = first_row:(first_row + length(rows) - 1)
    return Tuple(expected) == rows ? expected : nothing
end

function _pqs_source_box_route_driver_support_partition_unit(
    unit_key::Symbol,
    support_indices::AbstractVector{Int},
    retained_range,
    support_tiles::Tuple,
)
    observed = Int[]
    tile_fingerprints = map(support_tiles) do tile
        tile_rows = Vector{Int}(tile.support_indices)
        append!(observed, tile_rows)
        unit_local_rows =
            _pqs_source_box_route_driver_support_partition_local_rows(
                support_indices,
                tile_rows,
            )
        unit_local_row_range =
            _pqs_source_box_route_driver_contiguous_unit_range(unit_local_rows)
        (;
            tile_key = tile.tile_key,
            unit_key,
            source_region_role = tile.source_region_role,
            support_tile_role = tile.support_tile_role,
            cpb_role = tile.cpb_role,
            intervals = tile.intervals,
            support_count = tile.support_count,
            parent_row_count = length(tile_rows),
            parent_row_min = tile.parent_row_min,
            parent_row_max = tile.parent_row_max,
            parent_row_ids = Tuple(tile_rows),
            unit_local_row_range,
            unit_local_rows,
            retained_range,
            provider_tile_ready = tile.provider_tile_ready,
        )
    end
    expected_set = Set(Int.(support_indices))
    observed_set = Set(observed)
    duplicate_count = length(observed) - length(observed_set)
    missing_count = length(setdiff(expected_set, observed_set))
    outside_count = length(setdiff(observed_set, expected_set))
    coverage_complete =
        duplicate_count == 0 && missing_count == 0 && outside_count == 0 &&
        length(observed_set) == length(expected_set)
    return (;
        unit_key,
        support_count = length(support_indices),
        retained_range,
        tile_count = length(support_tiles),
        tile_support_counts = Tuple(tile.support_count for tile in support_tiles),
        parent_row_count = length(observed_set),
        parent_row_min = isempty(support_indices) ? nothing : minimum(support_indices),
        parent_row_max = isempty(support_indices) ? nothing : maximum(support_indices),
        duplicate_parent_row_count = duplicate_count,
        missing_parent_row_count = missing_count,
        outside_parent_row_count = outside_count,
        coverage_complete,
        tiles = support_tiles,
        tile_fingerprints = Tuple(tile_fingerprints),
    )
end

function _pqs_source_box_route_driver_independent_h2_pqs_supplement_support_partition_payload(
    target_payload,
    source_plan_payload,
)
    support_order = (:atom_contact_core, :shared_shell_1, :shared_shell_2)
    blocked(blocker) = _PQSIndependentH2PQSSupplementSupportPartitionPayload(
        :blocked_independent_h2_pqs_supplement_support_partition,
        blocker,
        isnothing(target_payload) ? :not_available : target_payload.route_family,
        isnothing(target_payload) ? :not_available : target_payload.route_kind,
        support_order,
        (),
        nothing,
        nothing,
        nothing,
        (),
        0,
        0,
        false,
        0,
        0,
        0,
    )
    isnothing(target_payload) && return blocked(:missing_independent_h2_target_payload)
    target_payload.status === :available_physical_gausslet_core_shell_target_inventory ||
        return blocked(target_payload.blocker)
    support_plan = get(target_payload.summary, :support_plan, nothing)
    isnothing(support_plan) && return blocked(:missing_independent_pqs_support_region_plan)
    source_plan =
        !isnothing(source_plan_payload) && hasproperty(source_plan_payload, :source_plan) ?
        source_plan_payload.source_plan :
        nothing
    isnothing(source_plan) && return blocked(:missing_independent_pqs_source_plan)
    source_plan.status ===
        :available_pqs_diatomic_physical_gausslet_core_shell_source_plan ||
        return blocked(source_plan.status)
    support_plan.support_order == support_order ||
        return blocked(:independent_pqs_support_partition_order_mismatch)
    source_plan.support_order == support_order ||
        return blocked(:independent_pqs_source_plan_support_order_mismatch)

    atom_descriptor = support_plan.atom_contact_core_descriptor
    shared_descriptors = support_plan.shared_shell_descriptors
    retained_ranges = source_plan.retained_ranges
    tile_sources = (;
        atom_contact_core = atom_descriptor.support_tiles,
        shared_shell_1 = shared_descriptors.shared_shell_1.support_tiles,
        shared_shell_2 = shared_descriptors.shared_shell_2.support_tiles,
    )
    unit_partitions = Tuple(
        _pqs_source_box_route_driver_support_partition_unit(
            unit_key,
            _pqs_source_box_route_driver_support_partition_rows(source_plan, unit_key),
            getproperty(retained_ranges, unit_key),
            getproperty(tile_sources, unit_key),
        ) for unit_key in support_order
    )
    duplicate_count =
        sum(partition.duplicate_parent_row_count for partition in unit_partitions)
    missing_count =
        sum(partition.missing_parent_row_count for partition in unit_partitions)
    outside_count =
        sum(partition.outside_parent_row_count for partition in unit_partitions)
    coverage_complete =
        all(partition.coverage_complete for partition in unit_partitions) &&
        duplicate_count == 0 && missing_count == 0 && outside_count == 0
    total_tile_count = sum(partition.tile_count for partition in unit_partitions)
    total_support_count = sum(partition.support_count for partition in unit_partitions)
    status =
        coverage_complete ?
        :available_independent_h2_pqs_supplement_support_partition :
        :blocked_independent_h2_pqs_supplement_support_partition
    blocker =
        coverage_complete ? nothing : :independent_h2_pqs_support_partition_coverage_mismatch
    support_counts =
        _pqs_source_box_route_driver_ordered_count_tuple(
            target_payload.support_counts,
            support_order,
        )
    retained_counts =
        _pqs_source_box_route_driver_ordered_count_tuple(
            target_payload.retained_counts,
            source_plan.retained_order,
        )
    return _PQSIndependentH2PQSSupplementSupportPartitionPayload(
        status,
        blocker,
        target_payload.route_family,
        target_payload.route_kind,
        support_order,
        source_plan.retained_order,
        support_counts,
        retained_counts,
        retained_ranges,
        unit_partitions,
        total_tile_count,
        total_support_count,
        coverage_complete,
        duplicate_count,
        missing_count,
        outside_count,
    )
end

function _pqs_source_box_route_driver_physical_gausslet_axis_metrics(axis_bundles)
    isnothing(axis_bundles) && return nothing
    pgdg_x = _nested_axis_pgdg(axis_bundles, :x)
    pgdg_y = _nested_axis_pgdg(axis_bundles, :y)
    pgdg_z = _nested_axis_pgdg(axis_bundles, :z)
    return (;
        x = (; overlap = pgdg_x.overlap, kinetic = pgdg_x.kinetic,
            source = :nested_pgdg_axis),
        y = (; overlap = pgdg_y.overlap, kinetic = pgdg_y.kinetic,
            source = :nested_pgdg_axis),
        z = (; overlap = pgdg_z.overlap, kinetic = pgdg_z.kinetic,
            source = :nested_pgdg_axis),
    )
end

function _pqs_source_box_route_driver_physical_gausslet_multilayer_shell_records(
    source_plan,
)
    records = NamedTuple[]
    for shell_index in eachindex(source_plan.shared_shell_support_indices)
        support_indices =
            Vector{Int}(source_plan.shared_shell_support_indices[shell_index])
        coefficient_matrix =
            source_plan.shared_shell_coefficient_matrices[shell_index]
        retained_count = size(coefficient_matrix, 2)
        local_coefficients =
            _pqs_source_box_route_driver_physical_gausslet_local_coefficients(
                coefficient_matrix,
                support_indices,
                retained_count,
                :physical_gausslet_shared_shell_coefficient_shape_mismatch,
            )
        isnothing(local_coefficients.blocker) ||
            throw(ArgumentError(string(local_coefficients.blocker)))
        push!(
            records,
            (;
                layer_index = shell_index,
                shell_support_indices = support_indices,
                shell_support_states =
                    source_plan.shared_shell_support_states[shell_index],
                shell_final_coefficients = local_coefficients.coefficients,
                support_count = length(support_indices),
                retained_count,
                provenance = :pqs_diatomic_physical_gausslet_source_plan,
            ),
        )
    end
    return Tuple(records)
end

function _pqs_source_box_route_driver_physical_gausslet_multilayer_plan(
    source_plan,
)
    source_plan.status ===
        :available_pqs_diatomic_physical_gausslet_core_shell_source_plan ||
        throw(ArgumentError("physical gausslet source plan must be available"))
    metrics =
        _pqs_source_box_route_driver_physical_gausslet_axis_metrics(
            source_plan.axis_bundles,
        )
    isnothing(metrics) &&
        throw(ArgumentError("physical gausslet source plan is missing axis metrics"))

    shell_records =
        _pqs_source_box_route_driver_physical_gausslet_multilayer_shell_records(
            source_plan,
        )
    shell_support_indices =
        reduce(vcat, (record.shell_support_indices for record in shell_records); init = Int[])
    shell_support_states =
        reduce(
            vcat,
            (record.shell_support_states for record in shell_records);
            init = Tuple{Int,Int,Int}[],
        )
    shell_final_coefficients =
        _pqs_multilayer_block_concatenate_shell_coefficients(shell_records)
    support_counts =
        (length(source_plan.atom_contact_core_support_indices),
            map(record -> record.support_count, shell_records)...)
    retained_counts =
        (length(source_plan.retained_ranges.atom_contact_core),
            map(record -> record.retained_count, shell_records)...)
    summary = (;
        status = :available_pqs_multilayer_shell_source_plan,
        blocker = nothing,
        source_plan_family =
            get(source_plan.convention_labels, :source_plan_family, :unknown),
        layer_count = length(shell_records),
        core_support_count = support_counts[1],
        shell_support_count = length(shell_support_indices),
        shell_final_retained_count = size(shell_final_coefficients, 2),
        support_counts,
        retained_counts,
        collapsed_shell_sector = true,
        final_basis_helper = :pqs_complete_core_shell_final_basis,
    )
    return (;
        object_kind = :pqs_multilayer_shell_source_plan,
        status = :available_pqs_multilayer_shell_source_plan,
        blocker = nothing,
        source_kind =
            :pqs_diatomic_physical_gausslet_core_shell_source_plan_adapter,
        bundles = source_plan.axis_bundles,
        metrics,
        core_box = nothing,
        outer_box = nothing,
        bond_axis = :z,
        layer_count = length(shell_records),
        shell_records,
        core_support_indices = source_plan.atom_contact_core_support_indices,
        core_support_states = source_plan.atom_contact_core_support_states,
        shell_support_indices,
        shell_support_states,
        shell_final_coefficients,
        summary,
        metadata = merge(
            NamedTuple(source_plan.metadata),
            (;
                source =
                    :pqs_source_box_route_driver_physical_gausslet_multilayer_plan,
                input_source_plan = source_plan.object_kind,
                route_owned_authority = true,
                collapsed_shell_sector = true,
            ),
        ),
    )
end

function _pqs_source_box_route_driver_physical_gausslet_common_final_basis(
    final_basis,
    source_plan,
)
    final_basis.status === :available_pqs_physical_gausslet_final_basis ||
        throw(ArgumentError("physical gausslet final basis must be available"))
    shell_support_indices =
        reduce(vcat, source_plan.shared_shell_support_indices; init = Int[])
    shell_retained_count =
        length(source_plan.retained_ranges.shared_shell_1) +
        length(source_plan.retained_ranges.shared_shell_2)
    return merge(
        NamedTuple(final_basis),
        (;
            object_kind = :pqs_complete_core_shell_final_basis,
            status = :available_pqs_complete_core_shell_final_basis,
            blocker = nothing,
            core_support_indices = source_plan.atom_contact_core_support_indices,
            shell_support_indices,
            support_row_order = :core_then_shell,
            core_support_count =
                length(source_plan.atom_contact_core_support_indices),
            shell_support_count = length(shell_support_indices),
            shell_final_retained_count = shell_retained_count,
            source_physical_final_basis_status = final_basis.status,
            metadata = merge(
                NamedTuple(final_basis.metadata),
                (;
                    source =
                        :pqs_source_box_route_driver_physical_gausslet_common_final_basis,
                    common_operator_input_adapter = true,
                    input_final_basis = :pqs_physical_gausslet_final_basis,
                ),
            ),
        ),
    )
end

function _pqs_source_box_route_driver_physical_gausslet_local_coefficients(
    coefficient_matrix,
    support_indices,
    retained_count::Int,
    blocker::Symbol,
)
    coefficients = Matrix{Float64}(coefficient_matrix)
    support = Vector{Int}(support_indices)
    shape = (length(support), retained_count)
    size(coefficients, 2) == retained_count || return (; blocker, coefficients = nothing)
    isempty(support) && return (; blocker, coefficients = nothing)
    maximum(support) <= size(coefficients, 1) || return (; blocker, coefficients = nothing)
    local_coefficients = coefficients[support, :]
    size(local_coefficients) == shape || return (; blocker, coefficients = nothing)
    all(isfinite, local_coefficients) || return (;
        blocker = :physical_gausslet_coefficient_nonfinite,
        coefficients = nothing,
    )
    return (; blocker = nothing, coefficients = local_coefficients)
end

function _pqs_source_box_route_driver_physical_gausslet_overlap_diagnostics(
    overlap::AbstractMatrix{<:Real};
    identity_atol::Real,
    rank_atol::Real,
)
    matrix = Matrix{Float64}(overlap)
    size(matrix, 1) == size(matrix, 2) || return (;
        blocker = :physical_gausslet_final_overlap_not_square,
        summary = (;),
        final_overlap = nothing,
        final_identity_error = nothing,
    )
    symmetry_error = norm(matrix - transpose(matrix), Inf)
    identity_error = norm(
        matrix - Matrix{Float64}(I, size(matrix, 1), size(matrix, 2)),
        Inf,
    )
    summary = (;
        shape = size(matrix),
        symmetry_error,
        identity_error,
        expected_dimension = size(matrix, 1),
        identity_atol = Float64(identity_atol),
    )
    blocker =
        identity_error <= Float64(identity_atol) ?
        nothing :
        :physical_gausslet_final_overlap_not_identity
    return (; blocker, summary, final_overlap = matrix, final_identity_error = identity_error)
end

function _pqs_source_box_route_driver_physical_gausslet_final_basis(
    source_plan;
    identity_atol::Real = 1.0e-8,
    rank_atol::Real = 1.0e-10,
)
    support_order = (:atom_contact_core, :shared_shell_1, :shared_shell_2)
    source_plan.support_order == support_order ||
        return (; status = :blocked_pqs_physical_gausslet_final_basis,
            blocker = :physical_gausslet_support_order_mismatch)
    source_plan.retained_order == support_order ||
        return (; status = :blocked_pqs_physical_gausslet_final_basis,
            blocker = :physical_gausslet_retained_order_mismatch)

    support_indices = (
        source_plan.atom_contact_core_support_indices,
        source_plan.shared_shell_support_indices[1],
        source_plan.shared_shell_support_indices[2],
    )
    support_states = (
        source_plan.atom_contact_core_support_states,
        source_plan.shared_shell_support_states[1],
        source_plan.shared_shell_support_states[2],
    )
    coefficient_matrices = (
        source_plan.core_coefficient_matrix,
        source_plan.shared_shell_coefficient_matrices[1],
        source_plan.shared_shell_coefficient_matrices[2],
    )
    retained_counts = (
        length(source_plan.retained_ranges.atom_contact_core),
        length(source_plan.retained_ranges.shared_shell_1),
        length(source_plan.retained_ranges.shared_shell_2),
    )
    support_counts = map(length, support_indices)
    support_counts == (275, 578, 362) ||
        return (; status = :blocked_pqs_physical_gausslet_final_basis,
            blocker = :physical_gausslet_support_count_mismatch)
    expected_retained_counts =
        get(source_plan.summary, :retained_counts, retained_counts)
    expected_retained_counts =
        expected_retained_counts isa Tuple ?
        Tuple(Int(count) for count in expected_retained_counts) :
        expected_retained_counts
    retained_counts == expected_retained_counts ||
        return (; status = :blocked_pqs_physical_gausslet_final_basis,
            blocker = :physical_gausslet_retained_count_mismatch)

    local_blocks = Matrix{Float64}[]
    for (index, matrix) in pairs(coefficient_matrices)
        blocker =
            index == 1 ?
            :physical_gausslet_core_coefficient_shape_mismatch :
            :physical_gausslet_shared_shell_coefficient_shape_mismatch
        local_coefficients =
            _pqs_source_box_route_driver_physical_gausslet_local_coefficients(
                matrix,
                support_indices[index],
                retained_counts[index],
                blocker,
            )
        isnothing(local_coefficients.blocker) || return (;
            status = :blocked_pqs_physical_gausslet_final_basis,
            blocker = local_coefficients.blocker,
        )
        push!(local_blocks, local_coefficients.coefficients)
    end

    metrics =
        _pqs_source_box_route_driver_physical_gausslet_axis_metrics(
            source_plan.axis_bundles,
        )
    isnothing(metrics) && return (;
        status = :blocked_pqs_physical_gausslet_final_basis,
        blocker = :missing_physical_gausslet_support_overlap,
    )

    combined_support_states = vcat(
        Vector{NTuple{3,Int}}(support_states[1]),
        Vector{NTuple{3,Int}}(support_states[2]),
        Vector{NTuple{3,Int}}(support_states[3]),
    )
    support_overlap = _pqs_multilayer_support_product_matrix(
        combined_support_states,
        combined_support_states,
        metrics.x.overlap,
        metrics.y.overlap,
        metrics.z.overlap,
    )
    total_support_count = sum(support_counts)
    final_dimension = sum(retained_counts)
    final_coefficients = zeros(Float64, total_support_count, final_dimension)
    support_row_ranges = (
        atom_contact_core = 1:support_counts[1],
        shared_shell_1 = (support_counts[1] + 1):(support_counts[1] + support_counts[2]),
        shared_shell_2 =
            (support_counts[1] + support_counts[2] + 1):total_support_count,
    )
    retained_ranges = (
        atom_contact_core = 1:retained_counts[1],
        shared_shell_1 = (retained_counts[1] + 1):(retained_counts[1] + retained_counts[2]),
        shared_shell_2 =
            (retained_counts[1] + retained_counts[2] + 1):final_dimension,
    )
    final_coefficients[support_row_ranges.atom_contact_core, retained_ranges.atom_contact_core] .=
        local_blocks[1]
    final_coefficients[support_row_ranges.shared_shell_1, retained_ranges.shared_shell_1] .=
        local_blocks[2]
    final_coefficients[support_row_ranges.shared_shell_2, retained_ranges.shared_shell_2] .=
        local_blocks[3]
    final_overlap =
        transpose(final_coefficients) * support_overlap * final_coefficients
    final_identity_error = norm(
        final_overlap -
        Matrix{Float64}(I, final_dimension, final_dimension),
        Inf,
    )
    diagnostics =
        _pqs_source_box_route_driver_physical_gausslet_overlap_diagnostics(
            final_overlap;
            identity_atol,
            rank_atol,
        )
    isnothing(diagnostics.blocker) || return (;
        object_kind = :pqs_physical_gausslet_final_basis,
        status = :blocked_pqs_physical_gausslet_final_basis,
        blocker = diagnostics.blocker,
        final_retained_count = final_dimension,
        final_overlap_identity_error = diagnostics.final_identity_error,
        support_counts,
        retained_counts,
        overlap_diagnostics = diagnostics.summary,
    )
    return (;
        object_kind = :pqs_physical_gausslet_final_basis,
        status = :available_pqs_physical_gausslet_final_basis,
        blocker = nothing,
        source_plan_object_kind = source_plan.object_kind,
        support_order,
        retained_order = support_order,
        support_counts,
        retained_counts,
        support_row_ranges,
        retained_ranges,
        final_retained_count = final_dimension,
        final_dimension,
        final_coefficients,
        final_overlap_identity_error = diagnostics.final_identity_error,
        final_overlap_is_identity = true,
        overlap_diagnostics = diagnostics.summary,
        transform_source_plan_provenance = source_plan.convention_labels,
        source_backed_fixed_source_oracle_used =
            get(source_plan.convention_labels, :source_backed_fixed_source_oracle_used, false),
        fake_pqs_enabled =
            get(source_plan.convention_labels, :fake_pqs_enabled, false),
        retained_transform_authority =
            get(
                source_plan.convention_labels,
                :retained_transform_authority,
                :not_available,
            ),
        final_basis_materialized = true,
        h1_materialized = false,
        h1_j_materialized = false,
        density_interaction_materialized = false,
        rhf_materialized = false,
        exports_materialized = false,
        metadata = (;
            source = :pqs_source_box_route_driver_physical_gausslet_final_basis,
            support_decomposition = :shared_physical_gausslet_core_shell,
            source_backed_adapter =
                get(source_plan.convention_labels, :source_backed_adapter, false),
            route_owned_authority = true,
            diagnostic_materialization = true,
            self_overlap_stored_for_downstream = false,
        ),
    )
end

function _pqs_source_box_route_driver_diatomic_physical_gausslet_final_basis_payload(
    route_skeleton,
    recipe,
    source_plan_payload = nothing,
)
    route_family =
        hasproperty(route_skeleton, :route_family) ?
        route_skeleton.route_family :
        recipe.route_family
    source_plan =
        !isnothing(source_plan_payload) && hasproperty(source_plan_payload, :source_plan) ?
        source_plan_payload.source_plan :
        nothing
    source_plan_status =
        isnothing(source_plan) ? :not_available : source_plan.status
    final_basis_requested =
        get(recipe, :run_final_basis, false) ||
        get(recipe, :run_h1, false) ||
        get(recipe, :run_h1_j, false) ||
        get(get(recipe, :private_rhf_inputs, (;)), :run_private_rhf, false)
    final_basis = nothing
    final_basis_status = :not_materialized_pqs_physical_gausslet_final_basis

    if route_family !== :pqs_source_box
        status = :not_applicable_pqs_physical_gausslet_final_basis_non_pqs_route
        blocker = nothing
    elseif source_plan_status !==
           :available_pqs_diatomic_physical_gausslet_core_shell_source_plan
        status = :blocked_pqs_physical_gausslet_final_basis_payload
        blocker = :missing_pqs_diatomic_physical_gausslet_source_plan
    elseif !final_basis_requested
        status = :not_requested_pqs_physical_gausslet_final_basis_payload
        blocker = :physical_gausslet_final_basis_request_not_enabled
    else
        final_basis =
            _pqs_source_box_route_driver_physical_gausslet_final_basis(
                source_plan,
            )
        final_basis_status = final_basis.status
        if final_basis_status === :available_pqs_physical_gausslet_final_basis
            status = :available_pqs_physical_gausslet_final_basis_payload
            blocker = nothing
        else
            status = :blocked_pqs_physical_gausslet_final_basis_payload
            blocker =
                isnothing(final_basis.blocker) ?
                :physical_gausslet_final_basis_blocked :
                final_basis.blocker
        end
    end

    support_counts =
        isnothing(source_plan) ? nothing : source_plan.summary.support_counts
    retained_counts =
        isnothing(source_plan) ? nothing : source_plan.summary.retained_counts
    final_dimension =
        !isnothing(final_basis) && hasproperty(final_basis, :final_dimension) ?
        final_basis.final_dimension :
        isnothing(source_plan) ? nothing : source_plan.final_dimension
    final_overlap_identity_error =
        !isnothing(final_basis) &&
        hasproperty(final_basis, :final_overlap_identity_error) ?
        final_basis.final_overlap_identity_error :
        nothing
    overlap_diagnostics =
        !isnothing(final_basis) &&
        hasproperty(final_basis, :overlap_diagnostics) ?
        final_basis.overlap_diagnostics :
        (;)
    summary = (;
        status,
        blocker,
        route_family,
        source_plan_status,
        final_basis_status,
        support_order =
            isnothing(source_plan) ? () : source_plan.support_order,
        retained_order =
            isnothing(source_plan) ? () : source_plan.retained_order,
        support_counts,
        retained_counts,
        final_dimension,
        final_overlap_identity_error,
        final_overlap_rank = get(overlap_diagnostics, :rank, nothing),
        final_overlap_full_rank = get(overlap_diagnostics, :full_rank, nothing),
        final_overlap_eigenvalue_min =
            get(overlap_diagnostics, :eigenvalue_min, nothing),
        final_overlap_eigenvalue_max =
            get(overlap_diagnostics, :eigenvalue_max, nothing),
        source_backed_fixed_source_oracle_used =
            !isnothing(final_basis) &&
            hasproperty(final_basis, :source_backed_fixed_source_oracle_used) ?
            final_basis.source_backed_fixed_source_oracle_used :
            false,
        fake_pqs_enabled =
            !isnothing(final_basis) &&
            hasproperty(final_basis, :fake_pqs_enabled) ?
            final_basis.fake_pqs_enabled :
            false,
        retained_transform_authority =
            !isnothing(final_basis) &&
            hasproperty(final_basis, :retained_transform_authority) ?
            final_basis.retained_transform_authority :
            isnothing(source_plan) ?
            :not_available :
            get(source_plan.summary, :retained_transform_authority, :not_available),
        final_basis_materialized =
            final_basis_status === :available_pqs_physical_gausslet_final_basis,
        endpoint_blocker =
            final_basis_status === :available_pqs_physical_gausslet_final_basis ?
            :missing_physical_gausslet_h1_builder :
            blocker,
        h1_materialized = false,
        h1_j_materialized = false,
        density_interaction_materialized = false,
        rhf_materialized = false,
        exports_materialized = false,
    )
    metadata = (;
        source =
            :pqs_source_box_route_driver_diatomic_physical_gausslet_final_basis_payload,
        route_kind = recipe.route_kind,
        route_family,
        source_plan_status,
        final_basis_status,
        support_decomposition = :shared_physical_gausslet_core_shell,
        h2_221_diagnostic_source_plan_reused = false,
    )
    return _PQSDiatomicPhysicalGaussletFinalBasisPayload(
        status,
        blocker,
        route_family,
        source_plan,
        source_plan_status,
        final_basis,
        final_basis_status,
        summary,
        metadata,
    )
end

function _pqs_source_box_route_driver_physical_gausslet_support_states(source_plan)
    return vcat(
        Vector{NTuple{3,Int}}(source_plan.atom_contact_core_support_states),
        Vector{NTuple{3,Int}}(source_plan.shared_shell_support_states[1]),
        Vector{NTuple{3,Int}}(source_plan.shared_shell_support_states[2]),
    )
end

function _pqs_source_box_route_driver_diatomic_physical_gausslet_h1_payload(
    parent,
    route_skeleton,
    recipe,
    source_plan_payload = nothing,
    final_basis_payload = nothing,
)
    route_family =
        hasproperty(route_skeleton, :route_family) ?
        route_skeleton.route_family :
        recipe.route_family
    source_plan =
        !isnothing(source_plan_payload) && hasproperty(source_plan_payload, :source_plan) ?
        source_plan_payload.source_plan :
        nothing
    final_basis =
        !isnothing(final_basis_payload) && hasproperty(final_basis_payload, :final_basis) ?
        final_basis_payload.final_basis :
        nothing
    source_plan_status = isnothing(source_plan) ? :not_available : source_plan.status
    final_basis_status = isnothing(final_basis) ? :not_available : final_basis.status
    h1_requested =
        get(recipe, :run_h1, false) ||
        get(recipe, :run_h1_j, false) ||
        get(get(recipe, :private_rhf_inputs, (;)), :run_private_rhf, false)
    support_kinetic = support_nuclear = final_kinetic = final_nuclear = nothing
    final_hamiltonian = h1 = nothing
    support_kinetic_status = :not_materialized_pqs_physical_gausslet_support_kinetic
    support_nuclear_status = :not_materialized_pqs_physical_gausslet_support_electron_nuclear
    final_kinetic_status = :not_materialized_pqs_physical_gausslet_final_kinetic
    final_nuclear_status = :not_materialized_pqs_physical_gausslet_final_electron_nuclear
    final_hamiltonian_status = :not_materialized_pqs_physical_gausslet_h1_hamiltonian
    h1_status = :not_materialized_pqs_physical_gausslet_h1_solve

    if route_family !== :pqs_source_box
        status = :not_applicable_pqs_physical_gausslet_h1_non_pqs_route
        blocker = nothing
    elseif source_plan_status !==
           :available_pqs_diatomic_physical_gausslet_core_shell_source_plan
        status = :blocked_pqs_physical_gausslet_h1_payload
        blocker = :missing_pqs_diatomic_physical_gausslet_source_plan
    elseif final_basis_status !== :available_pqs_physical_gausslet_final_basis
        status = :blocked_pqs_physical_gausslet_h1_payload
        blocker = :missing_pqs_physical_gausslet_final_basis
    elseif !h1_requested
        status = :not_requested_pqs_physical_gausslet_h1_payload
        blocker = :physical_gausslet_h1_request_not_enabled
    else
        center_records, missing_centers =
            _pqs_source_box_route_driver_complete_core_shell_center_records(parent)
        axis_layers, missing_axis_layers =
            _pqs_source_box_route_driver_complete_core_shell_axis_layers(
                source_plan.axis_bundles,
            )
        if !isempty((missing_centers..., missing_axis_layers...))
            status = :blocked_pqs_physical_gausslet_h1_payload
            blocker = :missing_physical_gausslet_h1_inputs
        else
            try
                coulomb_expansion = coulomb_gaussian_expansion(doacc = false)
                multilayer_plan =
                    _pqs_source_box_route_driver_physical_gausslet_multilayer_plan(
                        source_plan,
                    )
                common_final_basis =
                    _pqs_source_box_route_driver_physical_gausslet_common_final_basis(
                        final_basis,
                        source_plan,
                    )
                common_h1 = pqs_multilayer_complete_core_shell_h1_payload(
                    multilayer_plan;
                    final_basis = common_final_basis,
                    coulomb_expansion,
                    center_records,
                    axis_layers,
                    metadata = (;
                        source =
                            :pqs_source_box_route_driver_diatomic_physical_gausslet_h1_payload,
                        h1_operator_authority =
                            :pqs_multilayer_complete_core_shell_h1_payload,
                        input_source_plan = source_plan.object_kind,
                        input_final_basis = final_basis.object_kind,
                    ),
                )
                support_kinetic = common_h1.support_kinetic
                support_kinetic_status =
                    :materialized_pqs_multilayer_support_kinetic_matrix
                support_nuclear = common_h1.support_nuclear_by_center
                support_nuclear_status = support_nuclear.status
                final_kinetic = common_h1.final_kinetic
                final_kinetic_status = final_kinetic.status
                final_nuclear = common_h1.final_nuclear_by_center
                final_nuclear_status =
                    :materialized_pqs_physical_gausslet_final_electron_nuclear_by_center
                final_hamiltonian = common_h1.final_hamiltonian
                final_hamiltonian_status = final_hamiltonian.status
                h1 = common_h1.h1
                h1_status = h1.status
                if !(isfinite(h1.lowest_energy) && h1.lowest_energy < 0)
                    status = :blocked_pqs_physical_gausslet_h1_payload
                    blocker = :physical_gausslet_h1_lowest_energy_nonnegative
                else
                    status = :available_pqs_physical_gausslet_h1_payload
                    blocker = nothing
                end
            catch error
                error isa ArgumentError || error isa DimensionMismatch || rethrow()
                status = :blocked_pqs_physical_gausslet_h1_payload
                blocker = :physical_gausslet_h1_builder_error
                h1_status = :blocked_pqs_physical_gausslet_h1_solve
            end
        end
    end

    h1_materialized = status === :available_pqs_physical_gausslet_h1_payload
    summary = (;
        status,
        blocker,
        route_family,
        source_plan_status,
        final_basis_status,
        support_kinetic_status,
        support_electron_nuclear_status = support_nuclear_status,
        final_kinetic_status,
        final_electron_nuclear_status = final_nuclear_status,
        final_hamiltonian_status,
        h1_status,
        final_dimension =
            !isnothing(h1) ? h1.final_dimension :
            !isnothing(final_basis) ? final_basis.final_dimension : nothing,
        lowest_energy = isnothing(h1) ? nothing : h1.lowest_energy,
        h1_hamiltonian_matrix_finite =
            isnothing(final_hamiltonian) ? nothing :
            final_hamiltonian.hamiltonian_matrix_finite,
        h1_hamiltonian_symmetry_error =
            isnothing(final_hamiltonian) ? nothing :
            final_hamiltonian.hamiltonian_matrix_symmetry_error,
        center_count = isnothing(support_nuclear) ? nothing : support_nuclear.center_count,
        h1_materialized,
        endpoint_blocker =
            h1_materialized ? :missing_physical_gausslet_h1_j_builder : blocker,
        h1_j_materialized = false,
        density_interaction_materialized = false,
        rhf_materialized = false,
        exports_materialized = false,
    )
    metadata = (;
        source =
            :pqs_source_box_route_driver_diatomic_physical_gausslet_h1_payload,
        route_kind = recipe.route_kind,
        physical_final_basis = true,
        h1_operator_authority = :pqs_multilayer_complete_core_shell_h1_payload,
        h2_221_diagnostic_source_plan_reused = false,
    )
    return _PQSDiatomicPhysicalGaussletH1Payload(
        status,
        blocker,
        route_family,
        source_plan,
        source_plan_status,
        final_basis,
        final_basis_status,
        support_kinetic,
        support_kinetic_status,
        support_nuclear,
        support_nuclear_status,
        final_kinetic,
        final_kinetic_status,
        final_nuclear,
        final_nuclear_status,
        final_hamiltonian,
        final_hamiltonian_status,
        h1,
        h1_status,
        summary,
        metadata,
    )
end

function _pqs_source_box_route_driver_physical_gausslet_density_provenance(
    source_plan,
    coulomb_expansion,
)
    expected_term_count = length(coulomb_expansion.coefficients)
    provenance =
        CartesianContractedParentMetrics._pqs_source_box_ida_factor_provenance(
            source_plan.axis_bundles;
            expected_term_count,
        )
    hasproperty(provenance, :axis_weights) ||
        throw(ArgumentError("physical gausslet density provenance missing axis weights"))
    hasproperty(provenance, :raw_axis_pair_factor_terms) ||
        throw(ArgumentError("physical gausslet density provenance missing raw pair factor terms"))
    return (;
        object_kind = :pqs_physical_gausslet_density_provenance,
        status = :available_pqs_physical_gausslet_density_provenance,
        blocker = nothing,
        provenance,
        axis_weights = provenance.axis_weights,
        raw_pair_factor_terms = provenance.raw_axis_pair_factor_terms,
        term_count = provenance.term_count,
        factor_dimensions = provenance.factor_dimensions,
        metadata = (;
            source =
                :pqs_source_box_route_driver_physical_gausslet_density_provenance,
            provenance_source = :pqs_source_box_ida_factor_provenance,
            expected_term_count,
            retained_diagnostic_weights_are_ida_weights = false,
        ),
    )
end

function _pqs_source_box_route_driver_physical_gausslet_support_weights(
    source_plan;
    axis_weights,
)
    states = _pqs_source_box_route_driver_physical_gausslet_support_states(source_plan)
    weights_x, weights_y, weights_z =
        _pqs_multilayer_common_or_axis_tuple(axis_weights, "axis_weights")
    weights = Vector{Float64}(undef, length(states))
    @inbounds for (index, (ix, iy, iz)) in pairs(states)
        weights[index] =
            Float64(weights_x[ix]) * Float64(weights_y[iy]) * Float64(weights_z[iz])
    end
    return (;
        object_kind = :pqs_physical_gausslet_support_weights,
        status = :materialized_pqs_physical_gausslet_support_weights,
        blocker = nothing,
        support_weights = weights,
        support_weight_count = length(weights),
        support_weights_all_finite = all(isfinite, weights),
        support_weights_all_positive = all(>(0.0), weights),
        support_weight_min = minimum(weights),
        support_weight_max = maximum(weights),
        support_weight_sum = sum(weights),
        support_order = source_plan.support_order,
        metadata = (;
            source =
                :pqs_source_box_route_driver_physical_gausslet_support_weights,
            support_order = source_plan.support_order,
            density_gauge = :localized_ida,
        ),
    )
end

function _pqs_source_box_route_driver_physical_gausslet_support_pair_raw_numerator_matrix(
    source_plan;
    raw_pair_factor_terms,
    coulomb_expansion,
)
    states = _pqs_source_box_route_driver_physical_gausslet_support_states(source_plan)
    coefficients = Float64.(coulomb_expansion.coefficients)
    all(>(0.0), coefficients) ||
        throw(ArgumentError("physical gausslet raw pair numerator requires positive Coulomb expansion coefficients"))
    axis_terms =
        _pqs_multilayer_validate_raw_pair_factor_terms(
            _pqs_multilayer_raw_pair_factor_terms(raw_pair_factor_terms),
            length(coefficients),
        )
    support_pair_raw_numerator = zeros(Float64, length(states), length(states))
    for term_index in eachindex(coefficients)
        support_pair_raw_numerator .+=
            coefficients[term_index] *
            _pqs_multilayer_support_product_matrix(
                states,
                states,
                @view(axis_terms[1][term_index, :, :]),
                @view(axis_terms[2][term_index, :, :]),
                @view(axis_terms[3][term_index, :, :]),
            )
    end
    symmetry_error =
        norm(support_pair_raw_numerator - transpose(support_pair_raw_numerator), Inf)
    return (;
        object_kind = :pqs_physical_gausslet_support_raw_pair_numerator_matrix,
        status =
            :materialized_pqs_physical_gausslet_support_raw_pair_numerator_matrix,
        blocker = nothing,
        support_pair_raw_numerator,
        support_pair_raw_numerator_shape = size(support_pair_raw_numerator),
        support_pair_raw_numerator_finite = all(isfinite, support_pair_raw_numerator),
        support_pair_raw_numerator_symmetry_error = symmetry_error,
        support_order = source_plan.support_order,
        raw_pair_factor_convention = :raw_numerator,
        metadata = (;
            source =
                :pqs_source_box_route_driver_physical_gausslet_support_pair_raw_numerator_matrix,
            raw_pair_factor_convention = :raw_numerator,
        ),
    )
end

function _pqs_source_box_route_driver_physical_gausslet_ida_density_interaction(
    final_basis,
    support_pair_raw_numerator,
    support_weights;
    near_zero_atol::Real = 1.0e-12,
    symmetry_atol::Real = 1.0e-8,
    metadata = (;),
)
    get(final_basis, :object_kind, nothing) === :pqs_physical_gausslet_final_basis ||
        throw(ArgumentError("physical gausslet density interaction requires physical final basis"))
    get(final_basis, :status, nothing) === :available_pqs_physical_gausslet_final_basis ||
        throw(ArgumentError("physical gausslet density interaction requires available final basis"))
    get(final_basis, :final_basis_materialized, false) ||
        throw(ArgumentError("physical gausslet final basis is not materialized"))
    support_count = sum(final_basis.support_counts)
    raw_pair = Matrix{Float64}(support_pair_raw_numerator)
    size(raw_pair) == (support_count, support_count) ||
        throw(DimensionMismatch("physical gausslet raw pair numerator shape mismatch"))
    all(isfinite, raw_pair) ||
        throw(ArgumentError("physical gausslet raw pair numerator contains non-finite entries"))
    raw_pair_symmetry_error = norm(raw_pair - transpose(raw_pair), Inf)
    raw_pair_symmetry_error <= Float64(symmetry_atol) ||
        throw(ArgumentError("physical gausslet raw pair numerator must be symmetric"))

    weights = Float64[Float64(weight) for weight in support_weights]
    length(weights) == support_count ||
        throw(DimensionMismatch("physical gausslet support weight length mismatch"))
    all(isfinite, weights) ||
        throw(ArgumentError("physical gausslet support weights contain non-finite entries"))

    final_coefficients = Matrix{Float64}(final_basis.final_coefficients)
    size(final_coefficients, 1) == support_count ||
        throw(DimensionMismatch("physical gausslet final coefficient row mismatch"))
    size(final_coefficients, 2) == final_basis.final_retained_count ||
        throw(DimensionMismatch("physical gausslet final coefficient column mismatch"))

    ida_weights = vec(transpose(final_coefficients) * weights)
    near_zero_threshold = Float64(near_zero_atol)
    near_zero_count = count(weight -> abs(weight) <= near_zero_threshold, ida_weights)
    negative_count = count(<(0.0), ida_weights)
    positive_count = count(>(0.0), ida_weights)
    all_finite = all(isfinite, ida_weights)
    positive_weight_gauge =
        all_finite && near_zero_count == 0 && negative_count == 0 &&
        positive_count == length(ida_weights)

    if !positive_weight_gauge
        return (;
            object_kind = :pqs_physical_gausslet_ida_density_interaction,
            status =
                :blocked_pqs_physical_gausslet_ida_density_interaction,
            blocker = :physical_gausslet_ida_density_weights_not_positive,
            final_basis_object_kind = final_basis.object_kind,
            final_basis_status = final_basis.status,
            support_order = final_basis.support_order,
            density_gauge = :localized_ida,
            support_weight_count = length(weights),
            ida_weight_count = length(ida_weights),
            ida_weights,
            ida_weight_min = minimum(ida_weights),
            ida_weight_max = maximum(ida_weights),
            ida_weight_sum = sum(ida_weights),
            ida_weight_positive_count = positive_count,
            ida_weight_negative_count = negative_count,
            ida_weight_near_zero_count = near_zero_count,
            ida_weights_all_finite = all_finite,
            ida_weights_all_positive = false,
            raw_pair_numerator_shape = size(raw_pair),
            raw_pair_numerator_symmetry_error = raw_pair_symmetry_error,
            electron_electron_ida = nothing,
            electron_electron_ida_shape = nothing,
            electron_electron_ida_finite = false,
            electron_electron_ida_symmetry_error = nothing,
            metadata = merge(NamedTuple(metadata), (;
                source =
                    :pqs_source_box_route_driver_physical_gausslet_ida_density_interaction,
                raw_pair_factor_convention = :raw_numerator,
            )),
        )
    end

    weighted_coefficients =
        final_coefficients .* reshape(1.0 ./ ida_weights, 1, :)
    electron_electron_ida =
        transpose(weighted_coefficients) * raw_pair * weighted_coefficients
    ida_symmetry_error =
        norm(electron_electron_ida - transpose(electron_electron_ida), Inf)
    return (;
        object_kind = :pqs_physical_gausslet_ida_density_interaction,
        status = :materialized_pqs_physical_gausslet_ida_density_interaction,
        blocker = nothing,
        final_basis_object_kind = final_basis.object_kind,
        final_basis_status = final_basis.status,
        support_order = final_basis.support_order,
        density_gauge = :localized_ida,
        support_weight_count = length(weights),
        ida_weight_count = length(ida_weights),
        ida_weights,
        ida_weight_min = minimum(ida_weights),
        ida_weight_max = maximum(ida_weights),
        ida_weight_sum = sum(ida_weights),
        ida_weight_positive_count = positive_count,
        ida_weight_negative_count = negative_count,
        ida_weight_near_zero_count = near_zero_count,
        ida_weights_all_finite = all_finite,
        ida_weights_all_positive = true,
        raw_pair_numerator_shape = size(raw_pair),
        raw_pair_numerator_symmetry_error = raw_pair_symmetry_error,
        electron_electron_ida,
        electron_electron_ida_shape = size(electron_electron_ida),
        electron_electron_ida_finite = all(isfinite, electron_electron_ida),
        electron_electron_ida_symmetry_error = ida_symmetry_error,
        ida_weight_division_applied = true,
        signed_final_weight_division_used = false,
        raw_no_division_used = false,
        fixed_block_pair_data_authority_used = false,
        density_density_materialized = true,
        metadata = merge(NamedTuple(metadata), (;
            source =
                :pqs_source_box_route_driver_physical_gausslet_ida_density_interaction,
            raw_pair_factor_convention = :raw_numerator,
            weight_application_stage = :localized_ida_density_interaction_boundary,
            final_orbital_consumption_rule = :localized_ida_coefficients,
        )),
    )
end

function _pqs_source_box_route_driver_physical_gausslet_h1_j_diagnostic(
    density_interaction,
    final_orbital_coefficients,
)
    get(density_interaction, :status, nothing) ===
        :materialized_pqs_physical_gausslet_ida_density_interaction ||
        throw(ArgumentError("physical gausslet H1-J diagnostic requires materialized density interaction"))
    coefficients = Float64[Float64(value) for value in final_orbital_coefficients]
    length(coefficients) == density_interaction.ida_weight_count ||
        throw(DimensionMismatch("physical gausslet H1 orbital coefficient length mismatch"))
    orbital = coefficients
    pair_matrix = Matrix{Float64}(density_interaction.electron_electron_ida)
    density = orbital * transpose(orbital)
    rho = 0.5 .* (density .+ transpose(density))
    v = 0.5 .* (pair_matrix .+ transpose(pair_matrix))
    occupations = vec(diag(rho))
    direct = 2.0 * dot(occupations, v * occupations)
    exchange = dot(vec(rho), vec(v .* rho))
    self_coulomb = direct - exchange
    return (;
        object_kind = :pqs_physical_gausslet_h1_j_diagnostic,
        status = :materialized_pqs_physical_gausslet_h1_j_diagnostic,
        blocker = nothing,
        density_gauge = density_interaction.density_gauge,
        raw_pair_factor_convention = :raw_numerator,
        self_coulomb,
        h1_j_self_coulomb = self_coulomb,
        final_orbital_coefficient_count = length(coefficients),
        ida_orbital_coefficient_count = length(orbital),
        h1_self_coulomb_contraction =
            :restricted_one_orbital_direct_minus_exchange,
    )
end

function _pqs_source_box_route_driver_diatomic_physical_gausslet_h1_j_payload(
    route_skeleton,
    recipe,
    source_plan_payload = nothing,
    final_basis_payload = nothing,
    h1_payload = nothing,
)
    route_family =
        hasproperty(route_skeleton, :route_family) ?
        route_skeleton.route_family :
        recipe.route_family
    source_plan =
        !isnothing(source_plan_payload) && hasproperty(source_plan_payload, :source_plan) ?
        source_plan_payload.source_plan :
        nothing
    final_basis =
        !isnothing(final_basis_payload) && hasproperty(final_basis_payload, :final_basis) ?
        final_basis_payload.final_basis :
        nothing
    source_plan_status = isnothing(source_plan) ? :not_available : source_plan.status
    final_basis_status = isnothing(final_basis) ? :not_available : final_basis.status
    h1_payload_status = isnothing(h1_payload) ? :not_available : h1_payload.status
    h1_j_requested =
        get(recipe, :run_h1_j, false) ||
        get(get(recipe, :private_rhf_inputs, (;)), :run_private_rhf, false)
    density_provenance = support_weights = raw_pair_factor_terms = nothing
    support_pair_raw_numerator = density_interaction = h1_j_diagnostic = nothing
    density_provenance_status =
        :not_materialized_pqs_physical_gausslet_density_provenance
    support_weights_status =
        :not_materialized_pqs_physical_gausslet_support_weights
    raw_pair_factor_status =
        :not_materialized_pqs_physical_gausslet_raw_pair_factor_terms
    support_pair_raw_numerator_status =
        :not_materialized_pqs_physical_gausslet_support_raw_pair_numerator_matrix
    density_interaction_status =
        :not_materialized_pqs_physical_gausslet_ida_density_interaction
    h1_j_status = :not_materialized_pqs_physical_gausslet_h1_j_payload

    if route_family !== :pqs_source_box
        status = :not_applicable_pqs_physical_gausslet_h1_j_non_pqs_route
        blocker = nothing
    elseif source_plan_status !==
           :available_pqs_diatomic_physical_gausslet_core_shell_source_plan
        status = :blocked_pqs_physical_gausslet_h1_j_payload
        blocker = :missing_pqs_diatomic_physical_gausslet_source_plan
    elseif final_basis_status !== :available_pqs_physical_gausslet_final_basis
        status = :blocked_pqs_physical_gausslet_h1_j_payload
        blocker = :missing_pqs_physical_gausslet_final_basis
    elseif h1_payload_status !== :available_pqs_physical_gausslet_h1_payload
        status = :blocked_pqs_physical_gausslet_h1_j_payload
        blocker = :missing_pqs_physical_gausslet_h1_payload
    elseif !h1_j_requested
        status = :not_requested_pqs_physical_gausslet_h1_j_payload
        blocker = :physical_gausslet_h1_j_request_not_enabled
    else
        coulomb_expansion = coulomb_gaussian_expansion(doacc = false)
        try
            density_provenance =
                _pqs_source_box_route_driver_physical_gausslet_density_provenance(
                    source_plan,
                    coulomb_expansion,
                )
            density_provenance_status = density_provenance.status
            raw_pair_factor_terms = density_provenance.raw_pair_factor_terms
            raw_pair_factor_status =
                :available_pqs_physical_gausslet_raw_pair_factor_terms
            support_weights =
                _pqs_source_box_route_driver_physical_gausslet_support_weights(
                    source_plan;
                    axis_weights = density_provenance.axis_weights,
                )
            support_weights_status = support_weights.status
            support_weights.support_weights_all_positive ||
                throw(ArgumentError("physical gausslet support weights must be positive"))
            support_pair_raw_numerator =
                _pqs_source_box_route_driver_physical_gausslet_support_pair_raw_numerator_matrix(
                    source_plan;
                    raw_pair_factor_terms,
                    coulomb_expansion,
                )
            support_pair_raw_numerator_status =
                support_pair_raw_numerator.status
            density_interaction =
                _pqs_source_box_route_driver_physical_gausslet_ida_density_interaction(
                    final_basis,
                    support_pair_raw_numerator.support_pair_raw_numerator,
                    support_weights.support_weights;
                    metadata = (;
                        source =
                            :pqs_source_box_route_driver_diatomic_physical_gausslet_h1_j_payload,
                        raw_pair_factor_convention = :raw_numerator,
                    ),
                )
            density_interaction_status = density_interaction.status
            if density_interaction_status !==
               :materialized_pqs_physical_gausslet_ida_density_interaction
                status = :blocked_pqs_physical_gausslet_h1_j_payload
                blocker = :physical_gausslet_ida_density_interaction_blocked
            else
                h1_j_diagnostic =
                    _pqs_source_box_route_driver_physical_gausslet_h1_j_diagnostic(
                        density_interaction,
                        h1_payload.h1.lowest_orbital_coefficients,
                    )
                h1_j_status = :materialized_pqs_physical_gausslet_h1_j_payload
                status = h1_j_status
                blocker = nothing
            end
        catch error
            error isa ArgumentError || error isa DimensionMismatch || rethrow()
            status = :blocked_pqs_physical_gausslet_h1_j_payload
            message = sprint(showerror, error)
            if occursin("axis weights", message)
                blocker = :missing_physical_gausslet_axis_weights
            elseif occursin("raw pair", message)
                blocker = :missing_physical_gausslet_raw_pair_factor_terms
            elseif occursin("support weights", message)
                blocker = :physical_gausslet_support_weight_nonpositive
            else
                blocker = :physical_gausslet_ida_density_interaction_blocked
            end
        end
    end

    h1_j_materialized =
        status === :materialized_pqs_physical_gausslet_h1_j_payload
    density_interaction_materialized =
        density_interaction_status ===
        :materialized_pqs_physical_gausslet_ida_density_interaction
    summary = (;
        status,
        blocker,
        route_family,
        source_plan_status,
        final_basis_status,
        h1_payload_status,
        density_provenance_status,
        support_weights_status,
        raw_pair_factor_status,
        support_pair_raw_numerator_status,
        density_interaction_status,
        h1_j_status,
        density_gauge =
            isnothing(density_interaction) ? nothing : density_interaction.density_gauge,
        raw_pair_factor_convention =
            isnothing(support_pair_raw_numerator) ? nothing :
            support_pair_raw_numerator.raw_pair_factor_convention,
        support_weight_count =
            isnothing(support_weights) ? nothing : support_weights.support_weight_count,
        support_weights_all_positive =
            isnothing(support_weights) ? nothing :
            support_weights.support_weights_all_positive,
        support_raw_pair_shape =
            isnothing(support_pair_raw_numerator) ? nothing :
            support_pair_raw_numerator.support_pair_raw_numerator_shape,
        support_raw_pair_finite =
            isnothing(support_pair_raw_numerator) ? nothing :
            support_pair_raw_numerator.support_pair_raw_numerator_finite,
        electron_electron_ida_shape =
            isnothing(density_interaction) ? nothing :
            density_interaction.electron_electron_ida_shape,
        electron_electron_ida_finite =
            isnothing(density_interaction) ? nothing :
            density_interaction.electron_electron_ida_finite,
        electron_electron_ida_symmetry_error =
            isnothing(density_interaction) ? nothing :
            density_interaction.electron_electron_ida_symmetry_error,
        self_coulomb =
            isnothing(h1_j_diagnostic) ? nothing : h1_j_diagnostic.self_coulomb,
        h1_j_materialized,
        density_interaction_materialized,
        endpoint_blocker =
            h1_j_materialized ?
            :missing_physical_gausslet_rhf_or_solver_contract :
            blocker,
        rhf_materialized = false,
        exports_materialized = false,
    )
    metadata = (;
        source =
            :pqs_source_box_route_driver_diatomic_physical_gausslet_h1_j_payload,
        route_kind = recipe.route_kind,
        physical_final_basis = true,
        density_gauge = summary.density_gauge,
        raw_pair_factor_convention = summary.raw_pair_factor_convention,
        h2_221_diagnostic_source_plan_reused = false,
    )
    return _PQSDiatomicPhysicalGaussletH1JPayload(
        status,
        blocker,
        route_family,
        source_plan,
        source_plan_status,
        final_basis,
        final_basis_status,
        h1_payload,
        h1_payload_status,
        density_provenance,
        density_provenance_status,
        support_weights,
        support_weights_status,
        raw_pair_factor_terms,
        raw_pair_factor_status,
        support_pair_raw_numerator,
        support_pair_raw_numerator_status,
        density_interaction,
        density_interaction_status,
        h1_j_diagnostic,
        h1_j_status,
        summary,
        metadata,
    )
end

function _pqs_source_box_route_driver_diatomic_physical_gausslet_rhf_input_contract(
    parent,
    route_skeleton,
    recipe,
    source_plan_payload = nothing,
    final_basis_payload = nothing,
    h1_payload = nothing,
    h1_j_payload = nothing,
)
    route_family =
        hasproperty(route_skeleton, :route_family) ?
        route_skeleton.route_family :
        recipe.route_family
    route_kind =
        hasproperty(route_skeleton, :route_kind) ?
        route_skeleton.route_kind :
        recipe.route_kind
    source_plan =
        !isnothing(source_plan_payload) && hasproperty(source_plan_payload, :source_plan) ?
        source_plan_payload.source_plan :
        nothing
    final_basis =
        !isnothing(final_basis_payload) && hasproperty(final_basis_payload, :final_basis) ?
        final_basis_payload.final_basis :
        nothing
    source_plan_status = isnothing(source_plan) ? :not_available : source_plan.status
    final_basis_status = isnothing(final_basis) ? :not_available : final_basis.status
    h1_payload_status = isnothing(h1_payload) ? :not_available : h1_payload.status
    h1_j_payload_status =
        isnothing(h1_j_payload) ? :not_available : h1_j_payload.status
    private_rhf_inputs = get(recipe, :private_rhf_inputs, (;))
    fixture_role = get(private_rhf_inputs, :private_rhf_fixture_role, :route_smoke)
    electron_count =
        get(private_rhf_inputs, :private_rhf_electron_count, nothing)
    occupation =
        electron_count === 2 ?
        (; policy = :closed_shell_rhf, nocc = 1, occupancy = 2.0) :
        nothing
    h1_matrix =
        !isnothing(h1_payload) &&
        hasproperty(h1_payload, :final_hamiltonian) &&
        !isnothing(h1_payload.final_hamiltonian) ?
        h1_payload.final_hamiltonian.hamiltonian_matrix :
        nothing
    density_interaction =
        !isnothing(h1_j_payload) && hasproperty(h1_j_payload, :density_interaction) ?
        h1_j_payload.density_interaction :
        nothing
    electron_electron_ida =
        !isnothing(density_interaction) &&
        hasproperty(density_interaction, :electron_electron_ida) ?
        density_interaction.electron_electron_ida :
        nothing
    if route_family !== :pqs_source_box
        status = :not_applicable_pqs_physical_gausslet_rhf_input_contract
        blocker = nothing
    elseif route_kind ∉ (
        :bond_aligned_diatomic_physical_gausslet_core_shell_pqs,
        :bond_aligned_diatomic_independent_pqs_source_box_core_shell,
    )
        status = :blocked_pqs_physical_gausslet_rhf_input_contract
        blocker = :unsupported_physical_gausslet_route_kind
    elseif source_plan_status !==
           :available_pqs_diatomic_physical_gausslet_core_shell_source_plan
        status = :blocked_pqs_physical_gausslet_rhf_input_contract
        blocker = :missing_pqs_diatomic_physical_gausslet_source_plan
    elseif final_basis_status !== :available_pqs_physical_gausslet_final_basis
        status = :blocked_pqs_physical_gausslet_rhf_input_contract
        blocker = :missing_pqs_physical_gausslet_final_basis
    elseif h1_payload_status !== :available_pqs_physical_gausslet_h1_payload ||
           isnothing(h1_matrix)
        status = :blocked_pqs_physical_gausslet_rhf_input_contract
        blocker = :missing_physical_gausslet_h1_payload
    elseif h1_j_payload_status !== :materialized_pqs_physical_gausslet_h1_j_payload ||
           isnothing(density_interaction)
        status = :blocked_pqs_physical_gausslet_rhf_input_contract
        blocker = :missing_physical_gausslet_density_interaction
    elseif isnothing(electron_electron_ida)
        status = :blocked_pqs_physical_gausslet_rhf_input_contract
        blocker = :missing_physical_gausslet_electron_electron_ida
    elseif isnothing(electron_count)
        status = :blocked_pqs_physical_gausslet_rhf_input_contract
        blocker = :missing_physical_gausslet_electron_count
    elseif electron_count !== 2 || isnothing(occupation)
        status = :blocked_pqs_physical_gausslet_rhf_input_contract
        blocker = :unsupported_physical_gausslet_electron_count
    else
        h1 = Matrix{Float64}(h1_matrix)
        pair = Matrix{Float64}(electron_electron_ida)
        h1_symmetry_error = norm(h1 - transpose(h1), Inf)
        pair_symmetry_error = norm(pair - transpose(pair), Inf)
        if size(h1) != (final_basis.final_dimension, final_basis.final_dimension) ||
           !all(isfinite, h1) ||
           h1_symmetry_error > 1.0e-8
            status = :blocked_pqs_physical_gausslet_rhf_input_contract
            blocker = :missing_physical_gausslet_h1_payload
        elseif get(density_interaction, :density_gauge, nothing) !== :localized_ida
            status = :blocked_pqs_physical_gausslet_rhf_input_contract
            blocker = :physical_gausslet_rhf_input_contract_unreviewed
        elseif get(density_interaction.metadata, :raw_pair_factor_convention, nothing) !==
               :raw_numerator
            status = :blocked_pqs_physical_gausslet_rhf_input_contract
            blocker = :physical_gausslet_rhf_input_contract_unreviewed
        elseif size(pair, 1) != size(pair, 2) ||
               size(pair, 1) != size(h1, 1) ||
               !all(isfinite, pair) ||
               pair_symmetry_error > 1.0e-8
            status = :blocked_pqs_physical_gausslet_rhf_input_contract
            blocker = :missing_physical_gausslet_electron_electron_ida
        else
            status = :available_pqs_physical_gausslet_rhf_input_contract
            blocker = nothing
        end
    end

    contract_available =
        status === :available_pqs_physical_gausslet_rhf_input_contract
    h1_matrix_finite =
        isnothing(h1_matrix) ? nothing : all(isfinite, Matrix{Float64}(h1_matrix))
    h1_matrix_symmetry_error =
        isnothing(h1_matrix) ? nothing :
        norm(Matrix{Float64}(h1_matrix) - transpose(Matrix{Float64}(h1_matrix)), Inf)
    electron_electron_ida_finite =
        isnothing(electron_electron_ida) ? nothing :
        all(isfinite, Matrix{Float64}(electron_electron_ida))
    electron_electron_ida_symmetry_error =
        isnothing(electron_electron_ida) ? nothing :
        norm(
            Matrix{Float64}(electron_electron_ida) -
            transpose(Matrix{Float64}(electron_electron_ida)),
            Inf,
        )
    summary = (;
        status,
        blocker,
        route_family,
        route_kind,
        fixture_role,
        source_plan_status,
        final_basis_status,
        h1_payload_status,
        h1_j_payload_status,
        input_contract_available = contract_available,
        electron_count,
        occupation_policy = isnothing(occupation) ? nothing : occupation.policy,
        occupation_nocc = isnothing(occupation) ? nothing : occupation.nocc,
        final_dimension =
            isnothing(final_basis) ? nothing : final_basis.final_dimension,
        h1_matrix_available = !isnothing(h1_matrix),
        h1_matrix_finite,
        h1_matrix_symmetry_error,
        density_interaction_available = !isnothing(density_interaction),
        density_gauge =
            isnothing(density_interaction) ? nothing : density_interaction.density_gauge,
        raw_pair_factor_convention =
            isnothing(density_interaction) ? nothing :
            get(density_interaction.metadata, :raw_pair_factor_convention, nothing),
        electron_electron_ida_available = !isnothing(electron_electron_ida),
        electron_electron_ida_finite,
        electron_electron_ida_symmetry_error,
        private_diagnostic_only = true,
        private_rhf_materialized = false,
        endpoint_blocker =
            contract_available ?
            :missing_physical_gausslet_private_rhf_execution :
            blocker,
    )
    metadata = (;
        source =
            :pqs_source_box_route_driver_diatomic_physical_gausslet_rhf_input_contract,
        existing_complete_core_shell_rhf_helpers_adapted = false,
        scf_executed = false,
        private_diagnostic_only = true,
    )
    return _PQSDiatomicPhysicalGaussletRHFInputContractPayload(
        status,
        blocker,
        route_family,
        source_plan,
        source_plan_status,
        final_basis,
        final_basis_status,
        h1_payload,
        h1_payload_status,
        h1_j_payload,
        h1_j_payload_status,
        electron_count,
        occupation,
        summary,
        metadata,
    )
end

function _pqs_source_box_route_driver_diatomic_physical_gausslet_rhf_execution_payload(
    input_contract = nothing,
    h1_payload = nothing,
    h1_j_payload = nothing,
)
    input_contract_status =
        isnothing(input_contract) ? :not_available : input_contract.status
    h1_payload_status = isnothing(h1_payload) ? :not_available : h1_payload.status
    h1_j_payload_status =
        isnothing(h1_j_payload) ? :not_available : h1_j_payload.status
    contract_summary =
        isnothing(input_contract) || !hasproperty(input_contract, :summary) ?
        (;) :
        input_contract.summary
    route_family =
        isnothing(input_contract) ? :unknown : input_contract.route_family
    electron_count = get(contract_summary, :electron_count, nothing)
    occupation_nocc = get(contract_summary, :occupation_nocc, nothing)
    scf_payload = nothing
    if input_contract_status !==
       :available_pqs_physical_gausslet_rhf_input_contract
        status = :blocked_pqs_physical_gausslet_private_rhf_execution
        blocker =
            isnothing(input_contract) ?
            :missing_physical_gausslet_rhf_input_contract :
            input_contract.blocker
        executed = false
        materialized = false
        converged = false
        total_energy = nothing
        one_body_energy = nothing
        two_body_energy = nothing
        iteration_count = 0
        density_trace = nothing
        idempotency_residual = nothing
        commutator_residual = nothing
        energy_delta = nothing
        consistency_status =
            :not_evaluated_missing_physical_gausslet_rhf_input_contract
        endpoint_blocker = blocker
    else
        scf_payload = _pqs_multilayer_complete_core_shell_rhf_scf_payload(
            ;
            input_contract,
            h1_payload,
            h1_j_payload,
            mixing_kind = :fock_diis,
            max_iterations = 25,
            density_atol = 1.0e-8,
            energy_atol = 1.0e-10,
            residual_atol = 1.0e-8,
            trace_atol = 1.0e-8,
            idempotency_atol = 1.0e-8,
            metadata = (;
                source =
                    :pqs_source_box_route_driver_diatomic_physical_gausslet_rhf_execution_payload,
                private_diagnostic_only = true,
            ),
        )
        scf_summary = scf_payload.summary
        residual_diagnostics = get(scf_summary, :residual_diagnostics, (;))
        final_one_step = scf_payload.final_one_step_payload
        converged =
            scf_payload.status ===
            :materialized_pqs_multilayer_complete_core_shell_rhf_scf_payload
        materialized = converged
        executed = true
        status =
            converged ?
            :materialized_pqs_physical_gausslet_private_rhf_execution :
            :blocked_pqs_physical_gausslet_private_rhf_execution
        blocker =
            converged ? nothing :
            scf_payload.blocker === :scf_not_converged ?
            :physical_gausslet_private_rhf_not_converged :
            scf_payload.blocker
        total_energy =
            isnothing(final_one_step) ? nothing : final_one_step.total_energy
        one_body_energy =
            isnothing(final_one_step) ? nothing : final_one_step.one_body_energy
        two_body_energy =
            isnothing(final_one_step) ? nothing : final_one_step.two_body_energy
        iteration_count = get(scf_summary, :iteration_count, 0)
        density_trace = get(residual_diagnostics, :density_trace, nothing)
        idempotency_residual =
            get(residual_diagnostics, :closed_shell_idempotency_error, nothing)
        commutator_residual =
            get(residual_diagnostics, :commutator_residual, nothing)
        energy_delta = get(scf_summary, :final_energy_change, nothing)
        consistency_status =
            get(scf_summary, :final_one_step_recomputed, false) ?
            :reviewed_recomputed :
            :not_evaluated_missing_recomputed_final_one_step
        endpoint_blocker =
            converged ? :missing_h2_gausslet_only_reference_comparison : blocker
    end
    summary = (;
        status,
        blocker,
        input_contract_status,
        h1_payload_status,
        h1_j_payload_status,
        input_contract_available =
            input_contract_status ===
            :available_pqs_physical_gausslet_rhf_input_contract,
        executed,
        materialized,
        converged,
        electron_count,
        electron_count_source = :explicit_private_rhf_electron_count,
        occupation_policy = get(contract_summary, :occupation_policy, nothing),
        occupation_nocc,
        total_energy,
        one_body_energy,
        two_body_energy,
        iteration_count,
        density_trace,
        idempotency_residual,
        commutator_residual,
        energy_delta,
        final_density_one_step_consistency_status = consistency_status,
        endpoint_blocker,
        private_diagnostic_only = true,
    )
    metadata = (;
        source =
            :pqs_source_box_route_driver_diatomic_physical_gausslet_rhf_execution_payload,
        existing_complete_core_shell_rhf_helpers_adapted = true,
        execution_attempted = summary.executed,
        scf_payload_status = isnothing(scf_payload) ? nothing : scf_payload.status,
        scf_payload_blocker =
            isnothing(scf_payload) ? nothing : scf_payload.blocker,
        private_diagnostic_only = true,
    )
    return _PQSDiatomicPhysicalGaussletRHFExecutionPayload(
        status,
        blocker,
        route_family,
        input_contract,
        input_contract_status,
        h1_payload,
        h1_payload_status,
        h1_j_payload,
        h1_j_payload_status,
        summary,
        metadata,
    )
end

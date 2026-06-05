# Private route-report helpers for the atom-growth Cartesian diatomic driver path.
# This file is included after the generic route-driver helpers and White-Lindsey seed helpers.

function _pqs_source_box_route_driver_atom_growth_unit_key(prefix::Symbol, index::Int)
    return index == 1 ? prefix : Symbol(string(prefix), "_", index)
end

function _pqs_source_box_route_driver_atom_growth_unit_field(
    object,
    field::Symbol,
    default = nothing,
)
    return hasproperty(object, field) ? getproperty(object, field) : default
end

function _pqs_source_box_route_driver_atom_growth_lowering_parameters(region)
    piece = region.lowering_piece
    if !isnothing(piece) && hasproperty(piece, :lowering_parameters)
        return piece.lowering_parameters
    end
    return (;
        lowering_piece_object_kind = region.lowering_piece_object_kind,
        lowering_piece_role = region.lowering_piece_role,
        lowering_piece_support_count = region.lowering_piece_support_count,
    )
end

function _pqs_source_box_route_driver_atom_growth_source_descriptor(region)
    piece = region.lowering_piece
    if !isnothing(piece) &&
       region.lowering_piece_object_kind == :cartesian_outer_mismatch_boundary_slab_set3d
        slab_pieces =
            _pqs_source_box_route_driver_atom_growth_unit_field(
                piece,
                :slab_pieces,
                (),
            )
        return (;
            object_kind = :atom_growth_plan_unit_slab_set_descriptor,
            descriptor_kind = :outer_mismatch_boundary_slab_set,
            box = region.box,
            box_shape = region.box_shape,
            support_count = region.support_count,
            slab_piece_count = length(slab_pieces),
            slab_piece_roles = Tuple(
                _pqs_source_box_route_driver_atom_growth_unit_field(
                    slab_piece,
                    :role,
                ) for slab_piece in slab_pieces
            ),
            slab_piece_support_counts = Tuple(
                _pqs_source_box_route_driver_atom_growth_unit_field(
                    slab_piece,
                    :support_count,
                ) for slab_piece in slab_pieces
            ),
            final_column_ranges_available = false,
        )
    end
    return (;
        object_kind = :atom_growth_plan_unit_box_descriptor,
        descriptor_kind = :source_box,
        box = region.box,
        box_shape = region.box_shape,
        support_count = region.support_count,
        final_column_ranges_available = false,
    )
end

function _pqs_source_box_route_driver_atom_growth_plan_unit_record(
    region;
    unit_key::Symbol,
)
    return (;
        object_kind = :cartesian_atom_growth_plan_unit,
        unit_key,
        unit_role = region.role,
        region_order_index = region.order_index,
        source_descriptor =
            _pqs_source_box_route_driver_atom_growth_source_descriptor(region),
        source_box = region.box,
        source_dimensions = region.box_shape,
        source_dimension = region.support_count,
        support_count = region.support_count,
        lowering_family = region.lowering_family,
        lowering_parameters =
            _pqs_source_box_route_driver_atom_growth_lowering_parameters(region),
        lowering_piece_object_kind = region.lowering_piece_object_kind,
        materialization_dependency = region.materialization_dependency,
        source_backed = region.source_backed,
        independently_lowerable = region.independently_lowerable,
        retirement_target = region.retirement_target,
        retained_count = nothing,
        retained_range = nothing,
        retained_dimension = nothing,
        retained_count_known = false,
        retained_range_known = false,
        retained_dimension_known = false,
        materialized_units_available = false,
        provenance_label = :bond_aligned_diatomic_atom_growth_shellification_plan,
    )
end

function _pqs_source_box_route_driver_atom_growth_plan_unit_inventory(
    low_order_shellization,
)
    if !low_order_shellization.atom_growth_scaffold_available ||
       isnothing(low_order_shellization.atom_growth_scaffold)
        return (;
            object_kind = :cartesian_atom_growth_plan_unit_inventory,
            status = low_order_shellization.status,
            private_development_only = true,
            unit_inventory_source = :blocked_atom_growth_shellification_plan,
            plan_units = (),
            unit_count = 0,
            unit_keys = (),
            unit_roles = (),
            support_counts = (),
            materialization_dependencies = (),
            source_backed_region_count = 0,
            source_backed_unit_count = 0,
            plan_lowerable_unit_count = 0,
            retained_unit_dimensions_known = false,
            retained_unit_ranges_known = false,
            retained_dimension_known = false,
            retained_dimension = nothing,
            materialized_units_available = false,
            route_skeleton_authority = false,
            blocker =
                low_order_shellization.atom_growth_shellification_plan_status,
        )
    end

    scaffold = low_order_shellization.atom_growth_scaffold
    role_counts = Dict{Symbol,Int}()
    plan_units = NamedTuple[]
    for region in scaffold.regions
        role_index = get(role_counts, region.role, 0) + 1
        role_counts[region.role] = role_index
        push!(
            plan_units,
            _pqs_source_box_route_driver_atom_growth_plan_unit_record(
                region;
                unit_key =
                    _pqs_source_box_route_driver_atom_growth_unit_key(
                        region.role,
                        role_index,
                    ),
            ),
        )
    end

    plan_units = Tuple(plan_units)
    unit_keys = Tuple(unit.unit_key for unit in plan_units)
    support_counts =
        NamedTuple{unit_keys}(Tuple(unit.support_count for unit in plan_units))
    source_backed_unit_count = count(unit -> unit.source_backed, plan_units)
    return (;
        object_kind = :cartesian_atom_growth_plan_unit_inventory,
        status = :available_atom_growth_plan_unit_inventory,
        private_development_only = true,
        unit_inventory_source = :atom_growth_shellification_plan,
        plan_units,
        unit_count = length(plan_units),
        unit_keys,
        unit_roles = Tuple(unit.unit_role for unit in plan_units),
        support_counts,
        materialization_dependencies =
            Tuple(unit.materialization_dependency for unit in plan_units),
        source_backed_region_count =
            scaffold.materialization_dependency_counts.source_backed_region_count,
        source_backed_unit_count,
        plan_lowerable_unit_count =
            scaffold.materialization_dependency_counts.plan_lowerable_region_count,
        retained_unit_dimensions_known = false,
        retained_unit_ranges_known = false,
        retained_dimension_known = false,
        retained_dimension = nothing,
        materialized_units_available = false,
        route_skeleton_authority = false,
        blocker = nothing,
    )
end

function _pqs_source_box_route_driver_atom_growth_transform_contract_name(unit)
    unit.lowering_family == :white_lindsey_atom_local_child_shellification &&
        return :atom_local_child_shellification_sequence
    unit.materialization_dependency == :plan_lowerable_direct_slab &&
        return :direct_identity_selector
    unit.lowering_family == :white_lindsey_adaptive_complete_shell &&
        return :adaptive_complete_shell_layer
    unit.lowering_family == :outer_mismatch_boundary_slab_set &&
        return :outer_mismatch_boundary_slab_set
    return nothing
end

function _pqs_source_box_route_driver_atom_growth_transform_contract_record(
    unit,
)
    transform_contract =
        _pqs_source_box_route_driver_atom_growth_transform_contract_name(unit)
    isnothing(transform_contract) && throw(
        ArgumentError(
            "atom-growth plan unit $(unit.unit_key) has no transform contract for lowering family $(unit.lowering_family)",
        ),
    )
    return (;
        object_kind = :cartesian_atom_growth_transform_contract,
        unit_key = unit.unit_key,
        unit_role = unit.unit_role,
        transform_contract,
        contract_source = :atom_growth_plan_unit_inventory,
        lowering_family = unit.lowering_family,
        materialization_dependency = unit.materialization_dependency,
        source_backed = unit.source_backed,
        independently_lowerable = unit.independently_lowerable,
        coefficient_transform_materialized = false,
        coefficient_map_materialized = false,
        retained_count_known = unit.retained_count_known,
        retained_range_known = unit.retained_range_known,
        retained_dimension_known = unit.retained_dimension_known,
    )
end

function _pqs_source_box_route_driver_atom_growth_transform_contract_inventory(
    plan_unit_inventory,
)
    if isnothing(plan_unit_inventory) ||
       plan_unit_inventory.status != :available_atom_growth_plan_unit_inventory
        return (;
            object_kind = :cartesian_atom_growth_transform_contract_inventory,
            status = :blocked_missing_atom_growth_plan_unit_inventory,
            private_development_only = true,
            transform_contract_source = :blocked_missing_plan_unit_inventory,
            transform_contracts = (),
            contract_count = 0,
            unit_keys = (),
            unit_roles = (),
            contract_names = (),
            source_backed_contract_count = 0,
            coefficient_transforms_materialized = false,
            coefficient_maps_materialized = false,
            retained_unit_dimensions_known = false,
            retained_unit_ranges_known = false,
            retained_dimension_known = false,
            retained_dimension = nothing,
            route_skeleton_authority = false,
            blocker = :missing_atom_growth_plan_unit_inventory,
        )
    end

    unsupported_unit_keys = Tuple(
        unit.unit_key for unit in plan_unit_inventory.plan_units
        if isnothing(
            _pqs_source_box_route_driver_atom_growth_transform_contract_name(unit),
        )
    )
    if !isempty(unsupported_unit_keys)
        return (;
            object_kind = :cartesian_atom_growth_transform_contract_inventory,
            status = :blocked_unknown_atom_growth_transform_contract,
            private_development_only = true,
            transform_contract_source = :atom_growth_plan_unit_inventory,
            transform_contracts = (),
            contract_count = 0,
            unit_keys = (),
            unit_roles = (),
            contract_names = (),
            source_backed_contract_count = 0,
            coefficient_transforms_materialized = false,
            coefficient_maps_materialized = false,
            retained_unit_dimensions_known = false,
            retained_unit_ranges_known = false,
            retained_dimension_known = false,
            retained_dimension = nothing,
            route_skeleton_authority = false,
            blocker = (;
                reason = :unknown_atom_growth_transform_contract,
                unit_keys = unsupported_unit_keys,
            ),
        )
    end

    transform_contracts = Tuple(
        _pqs_source_box_route_driver_atom_growth_transform_contract_record(unit)
        for unit in plan_unit_inventory.plan_units
    )
    return (;
        object_kind = :cartesian_atom_growth_transform_contract_inventory,
        status = :available_atom_growth_transform_contract_inventory,
        private_development_only = true,
        transform_contract_source = :atom_growth_plan_unit_inventory,
        transform_contracts,
        contract_count = length(transform_contracts),
        unit_keys = Tuple(contract.unit_key for contract in transform_contracts),
        unit_roles = Tuple(contract.unit_role for contract in transform_contracts),
        contract_names =
            Tuple(contract.transform_contract for contract in transform_contracts),
        source_backed_contract_count =
            count(contract -> contract.source_backed, transform_contracts),
        coefficient_transforms_materialized = false,
        coefficient_maps_materialized = false,
        retained_unit_dimensions_known = false,
        retained_unit_ranges_known = false,
        retained_dimension_known = false,
        retained_dimension = nothing,
        route_skeleton_authority = false,
        blocker = nothing,
    )
end

function _pqs_source_box_route_driver_atom_growth_unit_record(;
    unit_key,
    unit_role,
    retained_unit_kind,
    source_family,
    source_box,
    retained_rule_kind,
    retained_rule_derivation,
    retained_range,
    provenance_label,
    source_dimensions = isnothing(source_box) ? nothing : Tuple(length.(source_box)),
    source_dimension =
        isnothing(source_dimensions) ?
        (isnothing(retained_range) ? nothing : length(retained_range)) :
        prod(source_dimensions),
)
    return _pqs_source_box_route_driver_unit_record(
        unit_key = unit_key,
        unit_role = unit_role,
        retained_unit_kind = retained_unit_kind,
        source_family = source_family,
        source_box = source_box,
        source_dimensions = source_dimensions,
        source_dimension = source_dimension,
        retained_rule_kind = retained_rule_kind,
        retained_rule_derivation = retained_rule_derivation,
        retained_range = retained_range,
        retained_count = isnothing(retained_range) ? nothing : length(retained_range),
        provenance_label = provenance_label,
        weight_semantics = :retained_basis_integral_weights,
    )
end

function _pqs_source_box_route_driver_diatomic_atom_growth_route_units(probe)
    if !probe.materialized || isnothing(probe.materialization) ||
       isnothing(probe.materialization.assembly)
        return (;
            object_kind = :white_lindsey_low_order_diatomic_atom_growth_route_units,
            route_family = :white_lindsey_low_order,
            status = :blocked_missing_atom_growth_assembly,
            private_development_only = true,
            retained_units = (),
            unit_inventory = nothing,
            standard_unit_inventory = nothing,
            retained_dimension = probe.retained_dimension,
            pair_entries = (),
            pair_family_counts = (white_lindsey_low_order_atom_growth = 0,),
            operator_pairs_materialized = false,
            weight_semantics = :retained_basis_integral_weights,
            blocker = :missing_atom_growth_assembly,
        )
    end

    scaffold = probe.scaffold
    assembly = probe.materialization.assembly
    retained_units = NamedTuple[]

    for (index, slab_set) in enumerate(scaffold.outer_mismatch_boundary_slab_sets)
        retained_range = assembly.outer_mismatch_column_ranges[index]
        push!(
            retained_units,
            _pqs_source_box_route_driver_atom_growth_unit_record(
                unit_key = _pqs_source_box_route_driver_atom_growth_unit_key(
                    :outer_mismatch_shared_molecular_shell,
                    index,
                ),
                unit_role = :outer_mismatch_shared_molecular_shell,
                retained_unit_kind =
                    :atom_growth_outer_mismatch_boundary_slab_set,
                source_family =
                    :white_lindsey_low_order_atom_growth_outer_mismatch,
                source_box = nothing,
                source_dimensions = nothing,
                source_dimension = slab_set.support_count,
                retained_rule_kind = :direct_boundary_slab_parent_sites,
                retained_rule_derivation =
                    :atom_growth_outer_mismatch_boundary_slab_set,
                retained_range = retained_range,
                provenance_label =
                    :bond_aligned_diatomic_atom_growth_outer_mismatch,
            ),
        )
    end

    push!(
        retained_units,
        _pqs_source_box_route_driver_atom_growth_unit_record(
            unit_key = :left_atom_box,
            unit_role = :left_atom_box,
            retained_unit_kind = :atom_growth_atom_local_child_box,
            source_family = :white_lindsey_low_order_atom_growth_atom_box,
            source_box = scaffold.left_child_plan.outer_box,
            retained_rule_kind = :atom_local_child_shellification_plan,
            retained_rule_derivation =
                :atom_growth_complete_rectangular_left_child_plan,
            retained_range = assembly.child_column_ranges[1],
            provenance_label = :bond_aligned_diatomic_atom_growth_left_atom_box,
        ),
    )

    if !isnothing(scaffold.contact_cap_region)
        push!(
            retained_units,
            _pqs_source_box_route_driver_atom_growth_unit_record(
                unit_key = :contact_cap,
                unit_role = :contact_cap,
                retained_unit_kind = :atom_growth_direct_contact_cap,
                source_family = :white_lindsey_low_order_atom_growth_contact_cap,
                source_box = scaffold.contact_cap_region.box,
                retained_rule_kind = :direct_contact_cap_parent_sites,
                retained_rule_derivation =
                    :atom_growth_complete_rectangular_contact_cap,
                retained_range = assembly.contact_cap_column_range,
                provenance_label =
                    :bond_aligned_diatomic_atom_growth_contact_cap,
            ),
        )
    end

    push!(
        retained_units,
        _pqs_source_box_route_driver_atom_growth_unit_record(
            unit_key = :right_atom_box,
            unit_role = :right_atom_box,
            retained_unit_kind = :atom_growth_atom_local_child_box,
            source_family = :white_lindsey_low_order_atom_growth_atom_box,
            source_box = scaffold.right_child_plan.outer_box,
            retained_rule_kind = :atom_local_child_shellification_plan,
            retained_rule_derivation =
                :atom_growth_complete_rectangular_right_child_plan,
            retained_range = assembly.child_column_ranges[2],
            provenance_label = :bond_aligned_diatomic_atom_growth_right_atom_box,
        ),
    )

    for (index, region) in enumerate(scaffold.shared_complete_shell_regions)
        retained_range = assembly.shared_shell_column_ranges[index]
        push!(
            retained_units,
            _pqs_source_box_route_driver_atom_growth_unit_record(
                unit_key = _pqs_source_box_route_driver_atom_growth_unit_key(
                    :regular_shared_molecular_shell,
                    index,
                ),
                unit_role = :regular_shared_molecular_shell,
                retained_unit_kind = :atom_growth_shared_complete_rectangular_shell,
                source_family =
                    :white_lindsey_low_order_atom_growth_shared_shell,
                source_box = region.outer_box,
                source_dimensions = nothing,
                source_dimension = region.support_count,
                retained_rule_kind = :shared_complete_rectangular_shell_plan,
                retained_rule_derivation =
                    :atom_growth_complete_rectangular_shared_shell_region,
                retained_range = retained_range,
                provenance_label =
                    :bond_aligned_diatomic_atom_growth_shared_shell,
            ),
        )
    end

    retained_units = Tuple(retained_units)
    unit_inventory = _pqs_source_box_route_driver_unit_inventory(retained_units)
    pair_entries = ()
    pair_family_counts = (white_lindsey_low_order_atom_growth = 0,)
    route_facts = (;
        source_dimensions = unit_inventory.source_dimensions,
        retained_units,
        retained_counts = unit_inventory.retained_counts,
        ranges = unit_inventory.ranges,
        retained_dimension = unit_inventory.retained_dimension,
        pair_entries,
        pair_family_counts,
    )

    return (;
        object_kind = :white_lindsey_low_order_diatomic_atom_growth_route_units,
        route_family = :white_lindsey_low_order,
        status = :available_atom_growth_retained_unit_inventory,
        private_development_only = true,
        retained_units,
        unit_inventory,
        standard_unit_inventory =
            _pqs_source_box_route_driver_standard_unit_inventory_summary(route_facts),
        retained_dimension = unit_inventory.retained_dimension,
        pair_entries,
        pair_family_counts,
        operator_pairs_materialized = false,
        pair_inventory_status =
            :assembled_sequence_payload_not_pair_decomposed,
        weight_semantics = :retained_basis_integral_weights,
        blocker = nothing,
    )
end

function _pqs_source_box_route_driver_diatomic_atom_growth_transform_inventory(
    probe,
)
    sequence = probe.materialization.sequence
    weights =
        isnothing(probe.basis_adapter) ?
        nothing :
        probe.basis_adapter.final_integral_weights
    retained_dimension = size(sequence.coefficient_matrix, 2)
    final_integral_weights_ready =
        !isnothing(weights) &&
        length(weights) == retained_dimension &&
        all(isfinite, weights) &&
        all(>(0.0), weights)

    return (;
        object_kind =
            :white_lindsey_low_order_diatomic_atom_growth_transform_inventory,
        route_family = :white_lindsey_low_order,
        status =
            final_integral_weights_ready ?
            :available_atom_growth_transform_inventory :
            :blocked_atom_growth_transform_inventory_contract,
        private_development_only = true,
        transform_source = :atom_growth_shell_sequence_coefficient_matrix,
        coefficient_matrix_size = size(sequence.coefficient_matrix),
        coefficient_matrix_finite = all(isfinite, sequence.coefficient_matrix),
        retained_dimension,
        support_count = length(sequence.support_indices),
        final_integral_weight_count = isnothing(weights) ? 0 : length(weights),
        final_integral_weights_status = probe.final_integral_weights_status,
        final_integral_weights_ready,
        weight_semantics = :retained_basis_integral_weights,
        blocker =
            final_integral_weights_ready ?
            nothing :
            :atom_growth_final_integral_weight_contract,
    )
end

function _pqs_source_box_route_driver_diatomic_atom_growth_operator_inventory(
    probe,
)
    fixed_block = probe.basis_adapter.fixed_block
    fixed_block_matrices =
        _white_lindsey_low_order_materialized_seed_operator_matrices(fixed_block)
    fixed_block_matrix_sizes =
        _white_lindsey_low_order_operator_matrix_sizes(fixed_block_matrices)
    fixed_block_finite_ready =
        _white_lindsey_low_order_operator_finite_ready(fixed_block_matrices)
    ham_adapter_available =
        !isnothing(probe.ham_adapter) &&
        probe.ham_adapter.status == :available_route_configured_diatomic_ham_adapter
    ham_matrices =
        ham_adapter_available ?
        (;
            overlap = probe.ham_adapter.operators.overlap,
            one_body_hamiltonian =
                probe.ham_adapter.operators.one_body_hamiltonian,
            interaction_matrix = probe.ham_adapter.operators.interaction_matrix,
        ) :
        nothing
    ham_matrix_sizes =
        ham_adapter_available ?
        _white_lindsey_low_order_operator_matrix_sizes(ham_matrices) :
        nothing
    ham_finite_ready =
        ham_adapter_available ?
        _white_lindsey_low_order_operator_finite_ready(ham_matrices) :
        nothing
    ham_all_finite =
        ham_adapter_available ? all(values(ham_finite_ready)) : false

    return (;
        object_kind =
            :white_lindsey_low_order_diatomic_atom_growth_operator_inventory,
        route_family = :white_lindsey_low_order,
        status =
            ham_adapter_available && ham_all_finite ?
            :available_atom_growth_operator_inventory :
            :blocked_atom_growth_operator_inventory_contract,
        private_development_only = true,
        operator_source = :atom_growth_fixed_block_and_ham_adapter,
        fixed_block_matrix_sizes,
        fixed_block_finite_ready,
        fixed_block_all_finite = all(values(fixed_block_finite_ready)),
        ham_adapter_status = probe.ham_adapter_status,
        ham_matrix_sizes,
        ham_finite_ready,
        ham_all_finite,
        final_integral_weights_status = probe.final_integral_weights_status,
        operator_pairs_materialized = false,
        pair_inventory_status =
            :assembled_sequence_payload_not_pair_decomposed,
        electron_electron_materialized = ham_adapter_available,
        overlap_materialized = ham_adapter_available,
        one_body_hamiltonian_materialized = ham_adapter_available,
        density_density_interaction_materialized = ham_adapter_available,
        blocker =
            ham_adapter_available && ham_all_finite ?
            nothing :
            :atom_growth_ham_operator_adapter_contract,
    )
end

function _pqs_source_box_route_driver_route_configured_diatomic_atom_growth_report(
    probe;
    basis_artifact_status,
    ham_artifact_status,
    basis_bundle_export_status,
    ham_bundle_export_status,
)
    if !probe.materialized || isnothing(probe.materialization)
        return (;
            object_kind =
                :white_lindsey_low_order_route_configured_diatomic_atom_growth_report,
            route_family = :white_lindsey_low_order,
            status = :blocked_missing_atom_growth_materialization,
            private_development_only = true,
            shellization_source =
                :bond_aligned_diatomic_atom_growth_construction_plan,
            shellization_authority =
                :bond_aligned_diatomic_atom_growth_construction_plan,
            active_source_authority = false,
            route_default_behavior_changed = false,
            sequence_available = false,
            retained_dimension = probe.retained_dimension,
            support_count = probe.support_count,
            route_units = nothing,
            transform_inventory = nothing,
            operator_inventory = nothing,
            final_integral_weights_status = probe.final_integral_weights_status,
            basis_adapter_status = probe.basis_adapter_status,
            ham_adapter_status = probe.ham_adapter_status,
            basis_artifact_status,
            ham_artifact_status,
            basis_bundle_export_status,
            ham_bundle_export_status,
            blocker = :missing_atom_growth_materialization,
        )
    elseif isnothing(probe.basis_adapter) ||
           isnothing(probe.basis_adapter.fixed_block)
        return (;
            object_kind =
                :white_lindsey_low_order_route_configured_diatomic_atom_growth_report,
            route_family = :white_lindsey_low_order,
            status = :blocked_missing_atom_growth_basis_adapter,
            private_development_only = true,
            shellization_source =
                :bond_aligned_diatomic_atom_growth_construction_plan,
            shellization_authority =
                :bond_aligned_diatomic_atom_growth_construction_plan,
            active_source_authority = false,
            route_default_behavior_changed = false,
            sequence_available = probe.sequence_available,
            retained_dimension = probe.retained_dimension,
            support_count = probe.support_count,
            route_units = nothing,
            transform_inventory = nothing,
            operator_inventory = nothing,
            final_integral_weights_status = probe.final_integral_weights_status,
            basis_adapter_status = probe.basis_adapter_status,
            ham_adapter_status = probe.ham_adapter_status,
            basis_artifact_status,
            ham_artifact_status,
            basis_bundle_export_status,
            ham_bundle_export_status,
            blocker = :atom_growth_basis_representation_contract,
        )
    end

    route_units =
        _pqs_source_box_route_driver_diatomic_atom_growth_route_units(probe)
    transform_inventory =
        _pqs_source_box_route_driver_diatomic_atom_growth_transform_inventory(probe)
    operator_inventory =
        _pqs_source_box_route_driver_diatomic_atom_growth_operator_inventory(probe)
    retained_dimension = probe.retained_dimension
    inventory_ready =
        route_units.status == :available_atom_growth_retained_unit_inventory &&
        transform_inventory.status == :available_atom_growth_transform_inventory &&
        operator_inventory.status == :available_atom_growth_operator_inventory

    return (;
        object_kind =
            :white_lindsey_low_order_route_configured_diatomic_atom_growth_report,
        route_family = :white_lindsey_low_order,
        status =
            inventory_ready ?
            :private_development_route_configured_atom_growth :
            :blocked_atom_growth_route_report_contract,
        private_development_only = true,
        materialization_status = probe.materialization.status,
        shellization_source =
            :bond_aligned_diatomic_atom_growth_construction_plan,
        shellization_authority =
            :bond_aligned_diatomic_atom_growth_construction_plan,
        active_source_authority = false,
        route_default_behavior_changed = false,
        sequence_available = probe.sequence_available,
        retained_dimension,
        support_count = probe.support_count,
        route_units,
        transform_inventory,
        operator_inventory,
        basis_adapter_summary = probe.basis_adapter_summary,
        ham_adapter_summary = probe.ham_adapter_summary,
        final_integral_weights_status = probe.final_integral_weights_status,
        final_integral_weights_ready =
            transform_inventory.final_integral_weights_ready,
        basis_adapter_status = probe.basis_adapter_status,
        ham_adapter_status = probe.ham_adapter_status,
        basis_artifact_status,
        ham_artifact_status,
        basis_bundle_export_status,
        ham_bundle_export_status,
        operator_pairs_materialized = route_units.operator_pairs_materialized,
        electron_electron_materialized = operator_inventory.electron_electron_materialized,
        weight_semantics = :retained_basis_integral_weights,
        blocker = inventory_ready ? nothing : :atom_growth_route_report_contract,
    )
end

function _pqs_source_box_route_driver_diatomic_atom_growth_materialization(
    context,
)
    (;
        report,
        route_family,
        save_basis_artifact,
        save_ham_artifact,
        basisfile,
        hamfile,
        route_configured_diatomic_ham_interaction_treatment,
        route_configured_shellization_request,
        route_configured_shellization_request_status,
        route_configured_system_classification,
        route_configured_system_classification_status,
        route_configured_bond_axis,
        route_configured_shellization_plan,
        route_configured_shellization_plan_status,
        route_configured_shellization_planning_status,
        route_configured_shellization_planning_family,
        route_configured_midpoint_slab_status,
        route_configured_shellization_helper_map,
        route_configured_shellization_helper_map_status,
        route_configured_primary_planned_helper,
        route_configured_missing_input_count,
        route_configured_helper_map_blocker,
        route_configured_input_readiness,
        route_configured_input_readiness_status,
        route_configured_available_fact_count,
        route_configured_materializer_missing_input_count,
        route_configured_input_readiness_blocker,
        route_configured_materializer_config,
        route_configured_materializer_config_status,
        route_configured_materializer_config_planning_family,
        route_configured_materializer_config_pending_input_count,
        route_configured_one_center_materializer_probe,
        route_configured_one_center_materializer_probe_requested,
        route_configured_one_center_materializer_probe_status,
        route_configured_one_center_materializer_probe_materialized,
        route_configured_one_center_materializer_probe_consumed,
        route_configured_one_center_materializer_probe_blocker,
        route_configured_diatomic_atom_growth_materializer_probe,
        route_configured_diatomic_atom_growth_materializer_probe_consumed,
        route_configured_diatomic_atom_growth_materializer_probe_blocker,
        route_configured_diatomic_atom_growth_basis_adapter_blocker,
        route_configured_diatomic_atom_growth_final_integral_weights_status,
        route_configured_diatomic_atom_growth_ham_adapter_status,
        route_configured_diatomic_atom_growth_ham_adapter_blocker,
        route_configured_diatomic_materializer_contract,
        route_configured_diatomic_atom_growth_materializer_contract,
        route_configured_materializer_contract,
    ) = context

    atom_growth_probe =
        route_configured_diatomic_atom_growth_materializer_probe
    atom_growth_basis_adapter = atom_growth_probe.basis_adapter
    atom_growth_ham_adapter = atom_growth_probe.ham_adapter
    atom_growth_basis_adapter_available =
        !isnothing(atom_growth_basis_adapter) &&
        atom_growth_basis_adapter.status ==
        :available_route_configured_diatomic_atom_growth_basis_adapter
    atom_growth_ham_adapter_available =
        !isnothing(atom_growth_ham_adapter) &&
        atom_growth_ham_adapter.status ==
        :available_route_configured_diatomic_ham_adapter
    atom_growth_materialized =
        route_configured_diatomic_atom_growth_materializer_probe_consumed
    atom_growth_export_blocker =
        !atom_growth_materialized ?
        something(
            route_configured_diatomic_atom_growth_materializer_probe_blocker,
            :atom_growth_materializer_not_consumed,
        ) :
        !atom_growth_basis_adapter_available ?
        something(
            route_configured_diatomic_atom_growth_basis_adapter_blocker,
            :atom_growth_basis_representation_contract,
        ) :
        nothing
    atom_growth_ham_export_blocker =
        !save_ham_artifact ?
        nothing :
        !atom_growth_basis_adapter_available ?
        atom_growth_export_blocker :
        !atom_growth_ham_adapter_available ?
        something(
            route_configured_diatomic_atom_growth_ham_adapter_blocker,
            :atom_growth_ham_operator_adapter_contract,
        ) :
        nothing
    atom_growth_artifact_export_requested =
        save_basis_artifact || save_ham_artifact
    atom_growth_retained_dimension =
        atom_growth_basis_adapter_available ?
        atom_growth_basis_adapter.retained_dimension :
        atom_growth_probe.retained_dimension
    atom_growth_basis_artifact_status =
        save_basis_artifact ?
        (
            atom_growth_basis_adapter_available ?
            :written_route_configured_diatomic_atom_growth_basis_only_bundle :
            :not_written_route_configured_diatomic_atom_growth_basis_adapter_blocked
        ) :
        :not_requested
    atom_growth_ham_artifact_status =
        save_ham_artifact ?
        (
            atom_growth_ham_adapter_available ?
            :written_route_configured_diatomic_atom_growth_ham_bundle :
            :not_written_route_configured_diatomic_atom_growth_ham_adapter_blocked
        ) :
        :not_requested
    atom_growth_basis_bundle_export_status =
        atom_growth_artifact_export_requested ?
        (
            atom_growth_basis_adapter_available ?
            :supported_route_configured_diatomic_atom_growth_basis_only_fixed_block :
            :pending_route_configured_diatomic_atom_growth_basis_export
        ) :
        :not_requested
    atom_growth_ham_bundle_export_status =
        atom_growth_artifact_export_requested && atom_growth_ham_adapter_available ?
        :available_route_configured_diatomic_atom_growth_ham_bundle_payload :
        save_ham_artifact ?
        something(
            atom_growth_ham_export_blocker,
            :pending_route_configured_diatomic_atom_growth_ham_export,
        ) :
        :not_requested
    atom_growth_ham_preflight_status =
        atom_growth_artifact_export_requested && atom_growth_ham_adapter_available ?
        :available_route_configured_diatomic_atom_growth_ham_adapter :
        save_ham_artifact ?
        route_configured_diatomic_atom_growth_ham_adapter_status :
        :not_requested
    atom_growth_ham_interaction_treatment_consumed =
        atom_growth_artifact_export_requested && atom_growth_ham_adapter_available ?
        atom_growth_ham_adapter.interaction_treatment :
        nothing
    atom_growth_ham_interaction_treatment_status =
        atom_growth_artifact_export_requested && atom_growth_ham_adapter_available ?
        :available_route_configured_diatomic_ham_interaction_treatment :
        save_ham_artifact ?
        route_configured_diatomic_atom_growth_ham_adapter_status :
        :not_requested
    atom_growth_materialized_report =
        atom_growth_materialized ?
        _pqs_source_box_route_driver_route_configured_diatomic_atom_growth_report(
            atom_growth_probe;
            basis_artifact_status = atom_growth_basis_artifact_status,
            ham_artifact_status = atom_growth_ham_artifact_status,
            basis_bundle_export_status =
                atom_growth_basis_bundle_export_status,
            ham_bundle_export_status =
                atom_growth_ham_bundle_export_status,
        ) :
        nothing
    atom_growth_materialized_report_kind =
        isnothing(atom_growth_materialized_report) ?
        nothing :
        atom_growth_materialized_report.object_kind
    atom_growth_artifact_meta = (;
        route_family,
        route_kind = report.recipe_metadata.route_kind,
        benchmark_role = report.recipe_metadata.benchmark_role,
        materialized_report_kind =
            something(
                atom_growth_materialized_report_kind,
                :not_materialized_atom_growth_probe,
            ),
        shellification_materialization_kind =
            atom_growth_materialized ?
            atom_growth_probe.materialization.object_kind :
            :not_materialized_atom_growth_probe,
        route_configured_shellization_request_status,
        route_configured_system_classification,
        route_configured_system_classification_status,
        route_configured_bond_axis,
        route_configured_shellization_plan_status,
        route_configured_shellization_planning_status,
        route_configured_shellization_planning_family,
        route_configured_midpoint_slab_status,
        route_configured_shellization_helper_map_status,
        route_configured_primary_planned_helper,
        route_configured_missing_input_count,
        route_configured_helper_map_blocker,
        route_configured_input_readiness_status,
        route_configured_available_fact_count,
        route_configured_materializer_missing_input_count,
        route_configured_input_readiness_blocker,
        route_configured_materializer_config_status,
        route_configured_materializer_config_planning_family,
        route_configured_materializer_config_pending_input_count,
        route_configured_one_center_materializer_probe_requested,
        route_configured_one_center_materializer_probe_status,
        route_configured_one_center_materializer_probe_materialized,
        route_configured_one_center_materializer_probe_consumed,
        route_configured_one_center_materializer_probe_blocker,
        route_configured_diatomic_materializer_contract...,
        route_configured_diatomic_atom_growth_materializer_contract...,
        route_configured_materializer_contract...,
        shellization_summary_available = false,
        shellization_source =
            :bond_aligned_diatomic_atom_growth_construction_plan,
        shellization_authority =
            :bond_aligned_diatomic_atom_growth_construction_plan,
        active_source_authority = false,
        route_configured_shellization_consumed = atom_growth_materialized,
        route_configured_diatomic_atom_growth_probe_consumed =
            route_configured_diatomic_atom_growth_materializer_probe_consumed,
        route_configured_legacy_diatomic_source_consumed = false,
        route_default_behavior_changed = false,
        materialized_shellization_stage =
            :atom_growth_complete_rectangular_low_order,
        seed_materialization_status =
            :not_seed_route_configured_diatomic_atom_growth_shellization,
        private_development_only = true,
    )
    atom_growth_basis_artifact_written = false
    if save_basis_artifact && atom_growth_basis_adapter_available
        write_cartesian_basis_bundle_jld2(
            basisfile,
            atom_growth_basis_adapter.fixed_block;
            include_ham = false,
            meta = (;
                atom_growth_artifact_meta...,
                export_status = :basis_only,
                basis_export_status =
                    atom_growth_basis_bundle_export_status,
                ham_export_status =
                    :artifact_local_basis_only_no_ham_payload,
                ham_export_blocker = nothing,
                companion_ham_artifact_requested = save_ham_artifact,
                companion_ham_artifact_status =
                    atom_growth_ham_artifact_status,
                companion_ham_export_status =
                    atom_growth_ham_bundle_export_status,
                companion_ham_export_blocker =
                    atom_growth_ham_export_blocker,
            ),
        )
        atom_growth_basis_artifact_written = true
    end
    atom_growth_ham_artifact_written = false
    if save_ham_artifact && atom_growth_ham_adapter_available
        write_cartesian_basis_bundle_jld2(
            hamfile,
            atom_growth_ham_adapter.operators;
            include_ham = true,
            meta = (;
                atom_growth_artifact_meta...,
                export_status = :basis_and_ham,
                basis_export_status =
                    atom_growth_basis_bundle_export_status,
                ham_preflight_status = atom_growth_ham_preflight_status,
                ham_operator_payload_status =
                    atom_growth_ham_adapter.operator_payload_status,
                ham_interaction_status =
                    atom_growth_ham_adapter.interaction_status,
                ham_export_status =
                    atom_growth_ham_bundle_export_status,
                ham_export_blocker = nothing,
            ),
        )
        atom_growth_ham_artifact_written = true
    end

    return (;
        object_kind = :cartesian_nesting_route_driver_materialization,
        route_family,
        private_development_only = true,
        materialize_route_requested = true,
        save_basis_artifact_requested = save_basis_artifact,
        save_ham_artifact_requested = save_ham_artifact,
        route_configured_diatomic_ham_interaction_treatment_requested =
            route_configured_diatomic_ham_interaction_treatment,
        route_configured_diatomic_ham_interaction_treatment_consumed =
            atom_growth_ham_interaction_treatment_consumed,
        route_configured_diatomic_ham_interaction_treatment_status =
            atom_growth_ham_interaction_treatment_status,
        status =
            atom_growth_export_blocker === nothing &&
            (!save_ham_artifact || atom_growth_ham_export_blocker === nothing) ?
            (
                atom_growth_artifact_export_requested ?
                :materialized_route_configured_diatomic_atom_growth_artifacts_available :
                :materialized_route_configured_diatomic_atom_growth_report_available
            ) :
            (
                atom_growth_artifact_export_requested ?
                :blocked_route_configured_diatomic_atom_growth_artifact_export :
                :blocked_route_configured_diatomic_atom_growth_materialization_report
            ),
        materialized_report = atom_growth_materialized_report,
        materialized_report_kind = atom_growth_materialized_report_kind,
        route_configured_shellization_request,
        route_configured_shellization_request_available = true,
        route_configured_shellization_request_status,
        route_configured_system_classification,
        route_configured_system_classification_status,
        route_configured_bond_axis,
        route_configured_shellization_plan,
        route_configured_shellization_plan_available = true,
        route_configured_shellization_plan_status,
        route_configured_shellization_planning_status,
        route_configured_shellization_planning_family,
        route_configured_midpoint_slab_status,
        route_configured_shellization_helper_map,
        route_configured_shellization_helper_map_available = true,
        route_configured_shellization_helper_map_status,
        route_configured_primary_planned_helper,
        route_configured_missing_input_count,
        route_configured_helper_map_blocker,
        route_configured_input_readiness,
        route_configured_input_readiness_available = true,
        route_configured_input_readiness_status,
        route_configured_available_fact_count,
        route_configured_materializer_missing_input_count,
        route_configured_input_readiness_blocker,
        route_configured_materializer_config,
        route_configured_materializer_config_available = true,
        route_configured_materializer_config_status,
        route_configured_materializer_config_planning_family,
        route_configured_materializer_config_pending_input_count,
        route_configured_one_center_materializer_probe,
        route_configured_one_center_materializer_probe_requested,
        route_configured_one_center_materializer_probe_status,
        route_configured_one_center_materializer_probe_materialized,
        route_configured_one_center_materializer_probe_consumed,
        route_configured_one_center_materializer_probe_blocker,
        route_configured_diatomic_atom_growth_materializer_probe,
        route_configured_diatomic_materializer_contract...,
        route_configured_diatomic_atom_growth_materializer_contract...,
        route_configured_materializer_contract...,
        route_configured_diatomic_basis_adapter_summary =
            atom_growth_probe.basis_adapter_summary,
        route_configured_diatomic_ham_adapter_summary =
            atom_growth_probe.ham_adapter_summary,
        shellization_summary = nothing,
        shellization_summary_available = false,
        shellization_source =
            :bond_aligned_diatomic_atom_growth_construction_plan,
        route_configured_shellization_consumed = atom_growth_materialized,
        route_configured_legacy_diatomic_source_consumed = false,
        materialized_shellization_stage =
            :atom_growth_complete_rectangular_low_order,
        seed_materialization_status =
            :not_seed_route_configured_diatomic_atom_growth_shellization,
        retained_dimension = atom_growth_retained_dimension,
        final_integral_weights_status =
            route_configured_diatomic_atom_growth_final_integral_weights_status,
        one_body_operator_status =
            atom_growth_ham_adapter_available ?
            :available_route_configured_diatomic_operator_payload :
            :not_requested,
        basis_bundle_export_status =
            atom_growth_basis_bundle_export_status,
        basis_artifact_status = atom_growth_basis_artifact_status,
        basis_artifact_written = atom_growth_basis_artifact_written,
        basisfile,
        basis_artifact_path =
            atom_growth_basis_artifact_written ? basisfile : nothing,
        basis_export_blocker =
            atom_growth_basis_adapter_available ?
            nothing :
            atom_growth_export_blocker,
        ham_preflight_status = atom_growth_ham_preflight_status,
        ham_missing_builder = atom_growth_ham_export_blocker,
        ham_operator_payload_status =
            atom_growth_ham_adapter_available ?
            atom_growth_ham_adapter.operator_payload_status :
            route_configured_diatomic_atom_growth_ham_adapter_status,
        ham_interaction_status =
            atom_growth_ham_adapter_available ?
            atom_growth_ham_adapter.interaction_status :
            route_configured_diatomic_atom_growth_ham_adapter_status,
        ham_bundle_export_status =
            atom_growth_ham_bundle_export_status,
        ham_artifact_status = atom_growth_ham_artifact_status,
        ham_artifact_written = atom_growth_ham_artifact_written,
        hamfile,
        ham_export_blocker = atom_growth_ham_export_blocker,
        ham_preflight = nothing,
        pqs_materialization_status = :not_applicable,
    )
end

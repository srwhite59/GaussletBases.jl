function _cartesian_shellization_sequence_summary(
    sequence::_CartesianNestedShellSequence3D;
    source_kind::Symbol,
    route_family::Symbol,
    shellization_role::Symbol,
    parent_box::NTuple{3,UnitRange{Int}} = sequence.working_box,
    split_status::Symbol = :not_applicable,
    bond_axis::Union{Nothing,Symbol} = nothing,
    midpoint_slab_present::Bool = false,
    child_sequence_count::Int = 0,
    shared_shell_layer_count::Int = length(sequence.shell_layers),
    child_column_ranges::Tuple = (),
    contact_or_merge_status::Symbol = :not_applicable,
)
    return (
        object_kind = :cartesian_shellization_route_summary,
        status = :private_development_summary,
        private_development_only = true,
        route_family = route_family,
        source_kind = source_kind,
        shellization_role = shellization_role,
        shellization_stage = :route_neutral_spatial_planning,
        lowering_stage = :not_lowered_by_shellization_summary,
        parent_box = parent_box,
        working_box = sequence.working_box,
        core_column_range = sequence.core_column_range,
        core_retained_count = length(sequence.core_column_range),
        shell_layer_count = length(sequence.shell_layers),
        shell_layer_kinds = _cartesian_shellization_layer_kinds(sequence),
        shell_layer_column_ranges = _cartesian_shellization_layer_column_ranges(sequence),
        retained_dimension = size(sequence.coefficient_matrix, 2),
        support_count = length(sequence.support_indices),
        split_status = split_status,
        bond_axis = bond_axis,
        midpoint_slab_present = midpoint_slab_present,
        child_sequence_count = child_sequence_count,
        shared_shell_layer_count = shared_shell_layer_count,
        child_column_ranges = child_column_ranges,
        contact_or_merge_status = contact_or_merge_status,
        diagnostics = _cartesian_shellization_common_diagnostics(
            source_kind = source_kind,
            shellization_role = shellization_role,
        ),
    )
end

function _cartesian_shellization_route_summary(
    sequence::_CartesianNestedShellSequence3D;
    route_family::Symbol = :one_center_full_parent,
    source_kind::Symbol = :one_center_shell_sequence,
    shellization_role::Symbol = :atom_local_full_parent_shellization,
)
    return _cartesian_shellization_sequence_summary(
        sequence;
        route_family = route_family,
        source_kind = source_kind,
        shellization_role = shellization_role,
    )
end

function _cartesian_shellization_route_summary(
    source::_CartesianNestedBondAlignedDiatomicSource3D;
    route_family::Symbol = :bond_aligned_diatomic,
    source_kind::Symbol = :bond_aligned_diatomic_source,
    shellization_role::Symbol = :bond_aligned_diatomic_shellization,
)
    geometry = source.geometry
    split_status = geometry.did_split ? :split : :no_split
    midpoint_slab_present = !isnothing(source.midpoint_slab_column_range)
    contact_or_merge_status =
        geometry.did_split ?
        (midpoint_slab_present ? :split_children_with_midpoint_slab : :split_children_merged) :
        :no_split_single_child

    return _cartesian_shellization_sequence_summary(
        source.sequence;
        route_family = route_family,
        source_kind = source_kind,
        shellization_role = shellization_role,
        parent_box = geometry.parent_box,
        split_status = split_status,
        bond_axis = geometry.bond_axis,
        midpoint_slab_present = midpoint_slab_present,
        child_sequence_count = length(source.child_sequences),
        shared_shell_layer_count = length(source.shared_shell_layers),
        child_column_ranges = Tuple(source.child_column_ranges),
        contact_or_merge_status = contact_or_merge_status,
    )
end

function _cartesian_shellization_route_location_tuple(location)
    length(location) == 3 || throw(
        ArgumentError("route shellization request expects three Cartesian coordinates"),
    )
    return (Float64(location[1]), Float64(location[2]), Float64(location[3]))
end

function _cartesian_shellization_route_diatomic_axis(atom_locations; atol::Float64 = 1.0e-12)
    first_location = _cartesian_shellization_route_location_tuple(atom_locations[1])
    second_location = _cartesian_shellization_route_location_tuple(atom_locations[2])
    deltas = abs.(
        (
            second_location[1] - first_location[1],
            second_location[2] - first_location[2],
            second_location[3] - first_location[3],
        ),
    )
    active_axes = Tuple(axis for (axis, delta) in zip((:x, :y, :z), deltas) if delta > atol)
    return length(active_axes) == 1 ? only(active_axes) : nothing
end

function _cartesian_shellization_route_system_classification(system_metadata)
    atom_count = length(system_metadata.atom_symbols)
    if atom_count == 1
        return (
            system_classification = :one_center,
            system_classification_status = :explicit_atom_count_one,
            bond_axis = nothing,
        )
    elseif atom_count == 2
        bond_axis = _cartesian_shellization_route_diatomic_axis(system_metadata.atom_locations)
        if isnothing(bond_axis)
            return (
                system_classification = :pending_system_classification,
                system_classification_status = :diatomic_not_axis_aligned_by_metadata,
                bond_axis = nothing,
            )
        end
        return (
            system_classification = :bond_aligned_diatomic,
            system_classification_status = :explicit_two_atom_single_axis_separation,
            bond_axis = bond_axis,
        )
    end
    return (
        system_classification = :unsupported_general_multi_atom,
        system_classification_status = :general_multi_atom_shellization_not_planned,
        bond_axis = nothing,
    )
end

function _cartesian_shellization_route_parent_contract(report)
    return hasproperty(report, :parent_contract) ? report.parent_contract : nothing
end

function _cartesian_shellization_route_report_classification(report, system_metadata)
    parent_contract = _cartesian_shellization_route_parent_contract(report)
    if !isnothing(parent_contract)
        return (
            system_classification = parent_contract.system_classification,
            system_classification_status =
                parent_contract.system_classification_status,
            bond_axis = parent_contract.bond_axis,
            chain_axis = parent_contract.chain_axis,
            classification_source = :parent_contract,
            parent_contract_available = true,
            parent_contract_status = parent_contract.status,
            parent_basis_materialization_status =
                parent_contract.parent_basis_materialization_status,
        )
    end

    classification = _cartesian_shellization_route_system_classification(system_metadata)
    return (
        system_classification = classification.system_classification,
        system_classification_status = classification.system_classification_status,
        bond_axis = classification.bond_axis,
        chain_axis = nothing,
        classification_source = :system_metadata_fallback,
        parent_contract_available = false,
        parent_contract_status = :not_available_legacy_report,
        parent_basis_materialization_status = :not_available_legacy_report,
    )
end

function _cartesian_shellization_route_get(values, field::Symbol, default = nothing)
    return hasproperty(values, field) ? getproperty(values, field) : default
end

function _cartesian_shellization_route_materializer_backend_request(
    system_metadata,
    recipe_metadata,
    materializer_backend,
)
    if !isnothing(materializer_backend)
        return (
            requested = materializer_backend,
            source = :materialization_input,
            status = :explicit_materializer_backend,
        )
    end

    map_backend = _cartesian_shellization_route_get(system_metadata, :map_backend)
    if !isnothing(map_backend)
        return (
            requested = map_backend,
            source = :system_map_backend,
            status = :defaulted_from_system_map_backend,
        )
    end

    parent_axis_probe_backend =
        _cartesian_shellization_route_get(recipe_metadata, :parent_axis_probe_backend)
    if !isnothing(parent_axis_probe_backend)
        return (
            requested = parent_axis_probe_backend,
            source = :parent_axis_probe_backend,
            status = :defaulted_from_parent_axis_probe_backend,
        )
    end

    raw_product_box_probe_backend =
        _cartesian_shellization_route_get(recipe_metadata, :raw_product_box_probe_backend)
    if !isnothing(raw_product_box_probe_backend)
        return (
            requested = raw_product_box_probe_backend,
            source = :raw_product_box_probe_backend,
            status = :defaulted_from_raw_product_box_probe_backend,
        )
    end

    return (
        requested = nothing,
        source = :unavailable,
        status = :missing_materializer_backend,
    )
end

function _cartesian_shellization_route_materializer_d_request(recipe_metadata)
    core_spacing = _cartesian_shellization_route_get(recipe_metadata, :core_spacing)
    if !isnothing(core_spacing)
        return (
            requested = Float64(core_spacing),
            source = :recipe_core_spacing,
            status = :defaulted_from_route_core_spacing,
        )
    end

    return (
        requested = nothing,
        source = :unavailable,
        status = :missing_route_core_spacing,
    )
end

function _cartesian_shellization_route_materializer_nside_request(
    recipe_metadata,
    materializer_nside,
)
    if !isnothing(materializer_nside)
        return (
            requested = Int(materializer_nside),
            source = :materialization_input,
            status = :explicit_materializer_nside,
        )
    end

    recipe_n_s = _cartesian_shellization_route_get(recipe_metadata, :n_s)
    if !isnothing(recipe_n_s)
        return (
            requested = Int(recipe_n_s),
            source = :recipe_n_s,
            status = :defaulted_from_route_n_s,
        )
    end

    return (
        requested = nothing,
        source = :unavailable,
        status = :missing_materializer_nside,
    )
end

function _cartesian_shellization_route_materializer_spacing_request(
    recipe_metadata,
    field::Symbol,
    source::Symbol,
    missing_status::Symbol,
)
    value = _cartesian_shellization_route_get(recipe_metadata, field)
    if !isnothing(value)
        return (
            requested = Float64(value),
            source = source,
            status = Symbol(:defaulted_from_, source),
        )
    end

    return (
        requested = nothing,
        source = :unavailable,
        status = missing_status,
    )
end

function _cartesian_shellization_route_configured_request(
    report;
    materializer_backend = nothing,
    materializer_nside = nothing,
)
    system_metadata = report.system_metadata
    recipe_metadata = report.recipe_metadata
    parent_contract = _cartesian_shellization_route_parent_contract(report)
    classification =
        _cartesian_shellization_route_report_classification(report, system_metadata)
    atom_count =
        isnothing(parent_contract) ?
        length(system_metadata.atom_symbols) :
        parent_contract.atom_count
    atom_symbols =
        isnothing(parent_contract) ?
        Tuple(system_metadata.atom_symbols) :
        parent_contract.atom_symbols
    atom_locations =
        isnothing(parent_contract) ?
        Tuple(system_metadata.atom_locations) :
        parent_contract.atom_locations
    nuclear_charges =
        isnothing(parent_contract) ?
        Tuple(system_metadata.nuclear_charges) :
        parent_contract.nuclear_charges
    parent_axis_counts =
        isnothing(parent_contract) ?
        system_metadata.parent_axis_counts :
        parent_contract.parent_axis_counts
    parent_axis_counts_source =
        isnothing(parent_contract) ?
        system_metadata.parent_axis_counts_source :
        parent_contract.parent_axis_counts_source
    parent_box =
        isnothing(parent_contract) ?
        system_metadata.parent_box :
        parent_contract.parent_box
    current_materialization_source =
        report.route_family == :white_lindsey_low_order ?
        :white_lindsey_one_center_seed :
        :pending_source_box_route_materializer
    backend_request = _cartesian_shellization_route_materializer_backend_request(
        system_metadata,
        recipe_metadata,
        materializer_backend,
    )
    d_request = _cartesian_shellization_route_materializer_d_request(recipe_metadata)
    nside_request =
        _cartesian_shellization_route_materializer_nside_request(
            recipe_metadata,
            materializer_nside,
        )
    reference_spacing_request =
        _cartesian_shellization_route_materializer_spacing_request(
            recipe_metadata,
            :reference_spacing,
            :recipe_reference_spacing,
            :missing_reference_spacing,
        )
    tail_spacing_request =
        _cartesian_shellization_route_materializer_spacing_request(
            recipe_metadata,
            :tail_spacing,
            :recipe_tail_spacing,
            :missing_tail_spacing,
        )
    missing_materializer_options = Tuple(
        option for (option, value) in (
            :gausslet_backend => backend_request.requested,
            :d => d_request.requested,
            :nside => nside_request.requested,
            :reference_spacing => reference_spacing_request.requested,
            :tail_spacing => tail_spacing_request.requested,
        ) if isnothing(value)
    )
    materializer_options_ready = isempty(missing_materializer_options)

    return (
        object_kind = :cartesian_shellization_route_configured_request,
        status = :metadata_only_pending_materializer,
        private_development_only = true,
        route_family = report.route_family,
        route_kind = recipe_metadata.route_kind,
        route_shape = recipe_metadata.route_shape,
        atom_count,
        atom_symbols,
        atom_locations,
        nuclear_charges,
        parent_axis_counts,
        parent_axis_counts_source,
        parent_box,
        requested_shellization_stage = :route_neutral_spatial_planning,
        expected_next_materializer_status =
            :pending_route_configured_shellization_materializer,
        current_materialization_source,
        materializer_backend_requested = backend_request.requested,
        materializer_backend_source = backend_request.source,
        materializer_backend_status = backend_request.status,
        materializer_d_requested = d_request.requested,
        materializer_d_source = d_request.source,
        materializer_d_status = d_request.status,
        materializer_nside_requested = nside_request.requested,
        materializer_nside_source = nside_request.source,
        materializer_nside_status = nside_request.status,
        materializer_reference_spacing_requested = reference_spacing_request.requested,
        materializer_reference_spacing_source = reference_spacing_request.source,
        materializer_reference_spacing_status = reference_spacing_request.status,
        materializer_tail_spacing_requested = tail_spacing_request.requested,
        materializer_tail_spacing_source = tail_spacing_request.source,
        materializer_tail_spacing_status = tail_spacing_request.status,
        materializer_options_ready,
        missing_materializer_options,
        materializer_option_blocker =
            materializer_options_ready ? nothing : :missing_route_materializer_options,
        route_configured_shellization_consumed = false,
        constructs_basis = false,
        constructs_shell_sequence = false,
        constructs_fixed_block = false,
        system_classification = classification.system_classification,
        system_classification_status = classification.system_classification_status,
        bond_axis = classification.bond_axis,
        chain_axis = classification.chain_axis,
        classification_source = classification.classification_source,
        parent_contract_available = classification.parent_contract_available,
        parent_contract_status = classification.parent_contract_status,
        parent_basis_materialization_status =
            classification.parent_basis_materialization_status,
    )
end

function _cartesian_shellization_route_planning_stage_statuses(request)
    if request.system_classification == :one_center
        return (
            planning_family = :one_center_atomic_shellization,
            planning_status = :planned_metadata_only,
            atom_local_core_status = :planned_single_atom_uncontracted_core,
            atom_local_shell_status = :planned_single_atom_local_shells,
            contact_merge_status = :not_applicable_one_center,
            midpoint_slab_status = :not_applicable_one_center,
            outer_rectangular_shell_status = :planned_after_atom_local_region,
            boundary_edge_adjustment_status = :planned_final_adjustment,
        )
    elseif request.system_classification == :bond_aligned_diatomic
        return (
            planning_family = :bond_aligned_diatomic_shellization,
            planning_status = :planned_metadata_only,
            atom_local_core_status = :planned_two_atom_uncontracted_cores,
            atom_local_shell_status = :planned_two_atom_local_shells,
            contact_merge_status = :planned_contact_or_merge_decision,
            midpoint_slab_status = :conditional_diatomic_midpoint_slab,
            outer_rectangular_shell_status = :planned_after_contact_merge_region,
            boundary_edge_adjustment_status = :planned_final_adjustment,
        )
    elseif request.system_classification == :pending_system_classification
        return (
            planning_family = :pending_system_classification,
            planning_status = :blocked_pending_system_classification,
            atom_local_core_status = :pending_system_classification,
            atom_local_shell_status = :pending_system_classification,
            contact_merge_status = :pending_system_classification,
            midpoint_slab_status = :pending_system_classification,
            outer_rectangular_shell_status = :pending_system_classification,
            boundary_edge_adjustment_status = :pending_system_classification,
        )
    end
    return (
        planning_family = :unsupported_general_multi_atom_shellization,
        planning_status = :unsupported_general_multi_atom,
        atom_local_core_status = :unsupported_general_multi_atom,
        atom_local_shell_status = :unsupported_general_multi_atom,
        contact_merge_status = :unsupported_general_multi_atom,
        midpoint_slab_status = :unsupported_general_multi_atom,
        outer_rectangular_shell_status = :unsupported_general_multi_atom,
        boundary_edge_adjustment_status = :unsupported_general_multi_atom,
    )
end

function _cartesian_shellization_route_planning_stub(request)
    stage_statuses = _cartesian_shellization_route_planning_stage_statuses(request)
    spatial_stage_order = (
        :atom_local_uncontracted_cores,
        :atom_local_shells,
        :contact_merge,
        :optional_midpoint_slab,
        :outer_rectangular_shell_boxes,
        :final_boundary_edge_adjustment,
    )

    return (
        object_kind = :cartesian_shellization_route_planning_stub,
        status = :metadata_only_pending_materializer,
        planning_status = stage_statuses.planning_status,
        private_development_only = true,
        request_object_kind = request.object_kind,
        route_family = request.route_family,
        route_kind = request.route_kind,
        system_classification = request.system_classification,
        system_classification_status = request.system_classification_status,
        bond_axis = request.bond_axis,
        planning_family = stage_statuses.planning_family,
        spatial_stage_order,
        atom_local_core_status = stage_statuses.atom_local_core_status,
        atom_local_shell_status = stage_statuses.atom_local_shell_status,
        contact_merge_status = stage_statuses.contact_merge_status,
        midpoint_slab_status = stage_statuses.midpoint_slab_status,
        outer_rectangular_shell_status = stage_statuses.outer_rectangular_shell_status,
        boundary_edge_adjustment_status = stage_statuses.boundary_edge_adjustment_status,
        expected_next_materializer_status = request.expected_next_materializer_status,
        constructs_basis = false,
        constructs_shell_sequence = false,
        constructs_fixed_block = false,
        route_configured_shellization_consumed = false,
    )
end

function _cartesian_shellization_route_planning_helper_map(plan)
    if plan.planning_family == :one_center_atomic_shellization
        missing_inputs = (
            :route_configured_one_center_basis,
            :coulomb_expansion,
            :refinement_levels,
            :packet_kernel,
            :fixed_block_materialization_handoff,
            :basis_or_ham_export_handoff,
        )
        return (
            object_kind = :cartesian_shellization_route_planning_helper_map,
            status = :metadata_only_pending_materializer_inputs,
            private_development_only = true,
            planning_family = plan.planning_family,
            primary_planned_helper = :_cartesian_materialize_shellification_low_order,
            helper_chain = (
                :_cartesian_shellification_plan_one_center_low_order,
                :_cartesian_materialize_shellification_low_order,
                :_nested_complete_rectangular_shell,
                :_nested_shell_sequence_from_core_block,
                :_nested_fixed_block,
            ),
            missing_inputs,
            missing_input_count = length(missing_inputs),
            blocker = :pending_route_configured_one_center_materializer_inputs,
            route_configured_shellization_consumed = false,
            calls_mapped_helpers = false,
        )
    elseif plan.planning_family == :bond_aligned_diatomic_shellization
        missing_inputs = (
            :bond_aligned_basis_or_source_object,
            :axis_bundles,
            :coulomb_expansion_coefficients,
            :split_options,
            :shared_shell_layer_policy,
            :packet_kernel,
            :fixed_block_adapter_handoff,
            :basis_or_ham_export_handoff,
        )
        return (
            object_kind = :cartesian_shellization_route_planning_helper_map,
            status = :metadata_only_pending_materializer_inputs,
            private_development_only = true,
            planning_family = plan.planning_family,
            primary_planned_helper = :_nested_bond_aligned_diatomic_source,
            helper_chain = (
                :_nested_bond_aligned_diatomic_source,
                :_nested_bond_aligned_diatomic_split_geometry,
                :_nested_bond_aligned_diatomic_sequence_for_box,
                :_nested_shell_sequence_from_core_block,
                :_nested_fixed_block,
            ),
            missing_inputs,
            missing_input_count = length(missing_inputs),
            blocker = :pending_route_configured_bond_aligned_diatomic_materializer_inputs,
            route_configured_shellization_consumed = false,
            calls_mapped_helpers = false,
        )
    elseif plan.planning_family == :pending_system_classification
        return (
            object_kind = :cartesian_shellization_route_planning_helper_map,
            status = :blocked_pending_system_classification,
            private_development_only = true,
            planning_family = plan.planning_family,
            primary_planned_helper = nothing,
            helper_chain = (),
            missing_inputs = (),
            missing_input_count = 0,
            blocker = :pending_system_classification,
            route_configured_shellization_consumed = false,
            calls_mapped_helpers = false,
        )
    end
    return (
        object_kind = :cartesian_shellization_route_planning_helper_map,
        status = :blocked_unsupported_general_multi_atom,
        private_development_only = true,
        planning_family = plan.planning_family,
        primary_planned_helper = nothing,
        helper_chain = (),
        missing_inputs = (),
        missing_input_count = 0,
        blocker = :unsupported_general_multi_atom_shellization,
        route_configured_shellization_consumed = false,
        calls_mapped_helpers = false,
    )
end

function _cartesian_shellization_route_materializer_input_readiness(
    request,
    plan,
    helper_map,
)
    available_facts = (
        :atom_symbols,
        :atom_locations,
        :nuclear_charges,
        :route_family,
        :route_kind,
        :parent_axis_counts,
        :parent_box,
        :requested_shellization_stage,
        :planning_family,
        :planned_helper_chain,
    )
    missing_inputs = helper_map.missing_inputs
    materializer_ready = isempty(missing_inputs) && isnothing(helper_map.blocker)
    status =
        materializer_ready ?
        :ready_for_route_configured_materializer :
        helper_map.status == :metadata_only_pending_materializer_inputs ?
        :blocked_missing_materializer_inputs :
        helper_map.status

    return (
        object_kind = :cartesian_shellization_route_materializer_input_readiness,
        status,
        private_development_only = true,
        route_family = request.route_family,
        route_kind = request.route_kind,
        planning_family = plan.planning_family,
        available_facts,
        available_fact_count = length(available_facts),
        missing_inputs,
        missing_input_count = length(missing_inputs),
        blocker = helper_map.blocker,
        materializer_backend_requested = request.materializer_backend_requested,
        materializer_backend_source = request.materializer_backend_source,
        materializer_backend_status = request.materializer_backend_status,
        materializer_d_requested = request.materializer_d_requested,
        materializer_d_source = request.materializer_d_source,
        materializer_d_status = request.materializer_d_status,
        materializer_nside_requested = request.materializer_nside_requested,
        materializer_nside_source = request.materializer_nside_source,
        materializer_nside_status = request.materializer_nside_status,
        materializer_reference_spacing_requested =
            request.materializer_reference_spacing_requested,
        materializer_reference_spacing_source =
            request.materializer_reference_spacing_source,
        materializer_reference_spacing_status =
            request.materializer_reference_spacing_status,
        materializer_tail_spacing_requested = request.materializer_tail_spacing_requested,
        materializer_tail_spacing_source = request.materializer_tail_spacing_source,
        materializer_tail_spacing_status = request.materializer_tail_spacing_status,
        materializer_options_ready = request.materializer_options_ready,
        missing_materializer_options = request.missing_materializer_options,
        materializer_option_blocker = request.materializer_option_blocker,
        driver_defaults_not_materializer_contract = (
            :refinement_options,
            :packet_kernel,
            :fixed_block_or_export_handoff,
        ),
        materializer_ready,
        route_configured_shellization_consumed = false,
        constructs_basis = false,
        constructs_axis_bundles = false,
        constructs_shell_sequence = false,
        constructs_fixed_block = false,
    )
end

function _cartesian_shellization_route_materializer_config(
    request,
    plan,
    helper_map,
    readiness,
)
    return (
        object_kind = :cartesian_shellization_route_materializer_config,
        status = readiness.status,
        private_development_only = true,
        route_family = request.route_family,
        route_kind = request.route_kind,
        system_classification = request.system_classification,
        system_classification_status = request.system_classification_status,
        bond_axis = request.bond_axis,
        atom_symbols = request.atom_symbols,
        atom_locations = request.atom_locations,
        nuclear_charges = request.nuclear_charges,
        parent_axis_counts = request.parent_axis_counts,
        parent_axis_counts_source = request.parent_axis_counts_source,
        parent_box = request.parent_box,
        materializer_backend_requested = readiness.materializer_backend_requested,
        materializer_backend_source = readiness.materializer_backend_source,
        materializer_backend_status = readiness.materializer_backend_status,
        materializer_d_requested = readiness.materializer_d_requested,
        materializer_d_source = readiness.materializer_d_source,
        materializer_d_status = readiness.materializer_d_status,
        materializer_nside_requested = readiness.materializer_nside_requested,
        materializer_nside_source = readiness.materializer_nside_source,
        materializer_nside_status = readiness.materializer_nside_status,
        materializer_reference_spacing_requested =
            readiness.materializer_reference_spacing_requested,
        materializer_reference_spacing_source =
            readiness.materializer_reference_spacing_source,
        materializer_reference_spacing_status =
            readiness.materializer_reference_spacing_status,
        materializer_tail_spacing_requested =
            readiness.materializer_tail_spacing_requested,
        materializer_tail_spacing_source = readiness.materializer_tail_spacing_source,
        materializer_tail_spacing_status = readiness.materializer_tail_spacing_status,
        materializer_options_ready = readiness.materializer_options_ready,
        missing_materializer_options = readiness.missing_materializer_options,
        materializer_option_blocker = readiness.materializer_option_blocker,
        planning_family = plan.planning_family,
        primary_planned_helper = helper_map.primary_planned_helper,
        helper_chain = helper_map.helper_chain,
        pending_inputs = readiness.missing_inputs,
        pending_input_count = readiness.missing_input_count,
        blocker = readiness.blocker,
        driver_defaults_not_materializer_contract =
            readiness.driver_defaults_not_materializer_contract,
        materializer_ready = readiness.materializer_ready,
        route_configured_shellization_consumed = false,
        constructs_basis = false,
        constructs_axis_bundles = false,
        constructs_shell_sequence = false,
        constructs_fixed_block = false,
    )
end

function _cartesian_shellization_route_axis_count_values(parent_axis_counts)
    if hasproperty(parent_axis_counts, :x)
        return (Int(parent_axis_counts.x), Int(parent_axis_counts.y), Int(parent_axis_counts.z))
    end
    length(parent_axis_counts) == 3 || throw(
        ArgumentError("one-center shellization materializer requires three parent axis counts"),
    )
    return (Int(parent_axis_counts[1]), Int(parent_axis_counts[2]), Int(parent_axis_counts[3]))
end

function _cartesian_shellization_route_origin_centered(location; atol::Float64 = 1.0e-12)
    coords = _cartesian_shellization_route_location_tuple(location)
    return all(coord -> isapprox(coord, 0.0; atol = atol, rtol = 0.0), coords)
end
function _cartesian_shellization_route_materialize_one_center_low_order(
    config;
    nside::Int = 5,
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    gausslet_backend::Symbol = :numerical_reference,
    refinement_levels::Integer = 0,
    packet_kernel::Symbol = :factorized_direct,
    d::Real = 0.2,
    tail_spacing::Real = 10.0,
    basis_family::Symbol = :G10,
    reference_spacing::Real = 1.0,
)
    config.system_classification == :one_center || throw(
        ArgumentError("one-center shellization materializer requires config.system_classification = :one_center"),
    )
    config.route_family == :white_lindsey_low_order || throw(
        ArgumentError("one-center shellization materializer is private for :white_lindsey_low_order configs"),
    )
    length(config.nuclear_charges) == 1 || throw(
        ArgumentError("one-center shellization materializer requires exactly one nuclear charge"),
    )
    length(config.atom_locations) == 1 || throw(
        ArgumentError("one-center shellization materializer requires exactly one atom location"),
    )
    _cartesian_shellization_route_origin_centered(only(config.atom_locations)) || throw(
        ArgumentError("one-center shellization materializer currently requires an origin-centered atom"),
    )
    axis_counts = _cartesian_shellization_route_axis_count_values(config.parent_axis_counts)
    axis_counts[1] == axis_counts[2] == axis_counts[3] || throw(
        ArgumentError("one-center shellization materializer currently requires cubic parent axis counts"),
    )

    parent_side_count = axis_counts[1]
    Z = Float64(only(config.nuclear_charges))
    basis = build_basis(
        MappedUniformBasisSpec(
            basis_family;
            count = parent_side_count,
            mapping = white_lindsey_atomic_mapping(;
                Z = Z,
                d = d,
                tail_spacing = tail_spacing,
            ),
            reference_spacing = reference_spacing,
        ),
    )
    bundle = _mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = expansion.exponents,
        center = 0.0,
        backend = gausslet_backend,
        refinement_levels = refinement_levels,
    )
    shellification_plan = _cartesian_shellification_plan_one_center_low_order(
        parent_side_count;
        nside,
        route_family = config.route_family,
    )
    shellification_plan_summary =
        _cartesian_shellification_plan_private_summary(shellification_plan)
    sequence = _cartesian_materialize_shellification_low_order(
        shellification_plan,
        bundle;
        expansion,
        packet_kernel,
    )
    fixed_block = _nested_fixed_block(sequence, basis, gausslet_backend)
    structure = one_center_atomic_nested_structure_diagnostics(
        sequence;
        parent_side_count = parent_side_count,
        nside = nside,
    )
    shellization_summary = _cartesian_shellization_route_summary(
        sequence;
        route_family = config.route_family,
        source_kind = :route_configured_one_center_low_order,
        shellization_role = :route_configured_one_center_full_parent_shellization,
    )
    inventory = _white_lindsey_low_order_materialized_seed_inventory(
        sequence,
        fixed_block,
        structure,
        packet_kernel = packet_kernel,
    )
    fixture = (;
        parent_side_count,
        nside,
        packet_kernel,
        basis,
        shellification_plan_summary,
        sequence,
        fixed_block,
        structure,
        shellization_summary,
        inventory,
        source_kind = :route_configured_one_center_low_order,
    )

    return (
        object_kind = :cartesian_shellization_route_one_center_materialization,
        status = :materialized_route_configured_one_center_low_order,
        private_development_only = true,
        route_family = config.route_family,
        route_kind = config.route_kind,
        planning_family = config.planning_family,
        consumed_config_fields = (
            :route_family,
            :route_kind,
            :system_classification,
            :atom_locations,
            :nuclear_charges,
            :parent_axis_counts,
            :materializer_backend_requested,
            :materializer_d_requested,
            :materializer_nside_requested,
            :materializer_reference_spacing_requested,
            :materializer_tail_spacing_requested,
            :planning_family,
            :primary_planned_helper,
            :helper_chain,
        ),
        materializer_options = (
            parent_side_count = parent_side_count,
            nside = nside,
            Z = Z,
            d = Float64(d),
            tail_spacing = Float64(tail_spacing),
            basis_family = basis_family,
            reference_spacing = Float64(reference_spacing),
            gausslet_backend = gausslet_backend,
            refinement_levels = Int(refinement_levels),
            packet_kernel = packet_kernel,
        ),
        fixture,
        sequence,
        fixed_block,
        shellization_summary,
        shellification_plan_summary,
        retained_dimension = size(sequence.coefficient_matrix, 2),
        route_configured_shellization_consumed = true,
        shellification_plan_path_used = true,
        calls_shellification_plan_materializer = true,
        calls_build_one_center_atomic_full_parent_shell_sequence = false,
        calls_white_lindsey_seed_fixture = false,
        calls_lower_level_one_center_helpers_directly = true,
        public_default_behavior_changed = false,
    )
end

function _cartesian_shellization_route_materialize_bond_aligned_diatomic(
    config;
    parent_qw_basis_object = nothing,
    parent_axis_bundle_object = nothing,
    expansion::Union{Nothing,CoulombGaussianExpansion} = nothing,
    term_coefficients::Union{Nothing,AbstractVector{<:Real}} = nothing,
    shared_shell_layer_policy = nothing,
    packet_kernel = nothing,
    axis_bundle_backend = nothing,
)
    config.system_classification == :bond_aligned_diatomic || throw(
        ArgumentError("diatomic shellization materializer requires config.system_classification = :bond_aligned_diatomic"),
    )
    config.route_family == :white_lindsey_low_order || throw(
        ArgumentError("diatomic shellization materializer is private for :white_lindsey_low_order configs"),
    )

    missing_contract = Symbol[]
    !config.materializer_options_ready &&
        append!(missing_contract, config.missing_materializer_options)
    isnothing(parent_qw_basis_object) &&
        push!(missing_contract, :parent_qw_basis_object_handoff)
    isnothing(parent_axis_bundle_object) &&
        push!(missing_contract, :parent_axis_bundle_object_handoff)
    isnothing(axis_bundle_backend) &&
        push!(missing_contract, :axis_bundle_backend_provenance)
    isnothing(shared_shell_layer_policy) &&
        push!(missing_contract, :shared_shell_layer_policy)
    isnothing(packet_kernel) && push!(missing_contract, :packet_kernel)
    isnothing(expansion) && isnothing(term_coefficients) &&
        push!(missing_contract, :coulomb_expansion_or_term_coefficients)
    if config.materializer_backend_requested == :pgdg_localized_experimental &&
       shared_shell_layer_policy != :endcap_panel_owned
        push!(missing_contract, :pgdg_requires_endcap_panel_owned_shared_shell_policy)
    end
    if !isnothing(axis_bundle_backend) &&
       axis_bundle_backend != config.materializer_backend_requested
        push!(missing_contract, :axis_bundle_backend_mismatches_materializer_backend)
    end
    missing_contract = Tuple(unique(missing_contract))

    materializer_options = (
        nside = config.materializer_nside_requested,
        d = config.materializer_d_requested,
        reference_spacing = config.materializer_reference_spacing_requested,
        tail_spacing = config.materializer_tail_spacing_requested,
        gausslet_backend = config.materializer_backend_requested,
        axis_bundle_backend = axis_bundle_backend,
        shared_shell_layer_policy = shared_shell_layer_policy,
        packet_kernel = packet_kernel,
        term_coefficients_source =
            isnothing(term_coefficients) ?
            (isnothing(expansion) ? nothing : :coulomb_expansion_coefficients) :
            :caller_supplied_term_coefficients,
    )

    if !isempty(missing_contract)
        return (
            object_kind = :cartesian_shellization_route_bond_aligned_diatomic_materialization,
            status = :blocked_missing_diatomic_materializer_contract,
            private_development_only = true,
            route_family = config.route_family,
            route_kind = config.route_kind,
            planning_family = config.planning_family,
            materialized = false,
            blocker = :pending_route_configured_bond_aligned_diatomic_materializer_contract,
            missing_contract,
            materializer_options,
            source = nothing,
            shellization_summary = nothing,
            shellification_plan_summary = nothing,
            retained_dimension = nothing,
            route_configured_shellization_consumed = false,
            shellification_plan_path_available = false,
            shellification_plan_path_used = false,
            calls_shellification_plan_materializer = false,
            public_default_behavior_changed = false,
        )
    end

    coefficients =
        isnothing(term_coefficients) ?
        Float64[Float64(value) for value in expansion.coefficients] :
        Float64[Float64(value) for value in term_coefficients]
    source = _nested_bond_aligned_diatomic_source(
        parent_qw_basis_object,
        parent_axis_bundle_object;
        bond_axis = config.bond_axis,
        nside = config.materializer_nside_requested,
        term_coefficients = coefficients,
        shared_shell_layer_policy = shared_shell_layer_policy,
        packet_kernel = packet_kernel,
    )
    shellization_summary = _cartesian_shellization_route_summary(
        source;
        route_family = config.route_family,
        source_kind = :route_configured_bond_aligned_diatomic_source,
        shellization_role = :route_configured_bond_aligned_diatomic_shellization,
    )
    shellification_plan =
        _cartesian_shellification_plan_bond_aligned_diatomic_low_order(
            source;
            route_family = config.route_family,
        )
    shellification_plan_summary =
        _cartesian_shellification_plan_private_summary(shellification_plan)

    return (
        object_kind = :cartesian_shellization_route_bond_aligned_diatomic_materialization,
        status = :materialized_route_configured_bond_aligned_diatomic_shellization,
        private_development_only = true,
        route_family = config.route_family,
        route_kind = config.route_kind,
        planning_family = config.planning_family,
        materialized = true,
        blocker = nothing,
        missing_contract = (),
        materializer_options,
        source,
        shellization_summary,
        shellification_plan_summary,
        retained_dimension = shellization_summary.retained_dimension,
        route_configured_shellization_consumed = true,
        shellification_plan_path_available = true,
        shellification_plan_path_used = false,
        calls_shellification_plan_materializer = false,
        public_default_behavior_changed = false,
    )
end

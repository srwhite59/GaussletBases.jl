function _pqs_source_box_route_driver_white_lindsey_preflight_fixed_block(seed_or_fixed_block)
    if hasproperty(seed_or_fixed_block, :fixture) &&
       hasproperty(seed_or_fixed_block.fixture, :fixed_block)
        return seed_or_fixed_block.fixture.fixed_block
    elseif hasproperty(seed_or_fixed_block, :fixed_block)
        return seed_or_fixed_block.fixed_block
    end
    return seed_or_fixed_block
end

function _pqs_source_box_route_driver_white_lindsey_ham_preflight(
    seed_or_fixed_block;
    ham_bundle_adapter = nothing,
)
    fixed_block =
        isnothing(ham_bundle_adapter) ?
        _pqs_source_box_route_driver_white_lindsey_preflight_fixed_block(
            seed_or_fixed_block,
        ) :
        ham_bundle_adapter.fixed_block
    ordinary_qw_fixed_block_applicable =
        applicable(ordinary_cartesian_qiu_white_operators, fixed_block)
    nested_cartesian_fixed_block_applicable =
        applicable(nested_cartesian_operators, fixed_block)
    ida_builder_name_defined = isdefined(@__MODULE__, :ordinary_cartesian_ida_operators)
    ordinary_cartesian_ida_fixed_block_applicable =
        ida_builder_name_defined &&
        applicable(getfield(@__MODULE__, :ordinary_cartesian_ida_operators), fixed_block)
    bundle_object = isnothing(ham_bundle_adapter) ? fixed_block : ham_bundle_adapter
    basis_bundle_payload = cartesian_basis_bundle_payload(bundle_object; include_ham = true)
    basis_bundle_ham_payload_available = !isnothing(basis_bundle_payload.ham)
    pure_operator_payload_available =
        ordinary_qw_fixed_block_applicable ||
        nested_cartesian_fixed_block_applicable ||
        ordinary_cartesian_ida_fixed_block_applicable ||
        basis_bundle_ham_payload_available
    missing_builder =
        pure_operator_payload_available ?
        nothing :
        :missing_pure_low_order_fixed_block_density_density_interaction_builder
    status =
        basis_bundle_ham_payload_available && !isnothing(ham_bundle_adapter) ?
        :available_private_low_order_ham_bundle_adapter :
        pure_operator_payload_available ?
        :available_pure_low_order_operator_payload :
        :blocked_missing_pure_low_order_interaction_builder
    ham_operator_payload_status =
        pure_operator_payload_available ?
        :available_low_order_operator_payload :
        :pending_low_order_operator_payload
    ham_interaction_status =
        basis_bundle_ham_payload_available ?
        :available_low_order_density_density_interaction_matrix :
        :pending_low_order_density_density_interaction_matrix
    ham_bundle_export_status =
        basis_bundle_ham_payload_available ?
        :available_low_order_ham_bundle_payload :
        :pending_low_order_density_density_interaction_matrix

    return (;
        object_kind = :white_lindsey_low_order_ham_preflight,
        route_family = :white_lindsey_low_order,
        fixed_block_type_label = string(typeof(fixed_block)),
        parent_basis_type_label =
            hasproperty(fixed_block, :parent_basis) ?
            string(typeof(fixed_block.parent_basis)) :
            "unavailable",
        ordinary_qw_fixed_block_applicable,
        nested_cartesian_fixed_block_applicable,
        ordinary_cartesian_ida_builder_name_defined = ida_builder_name_defined,
        ordinary_cartesian_ida_fixed_block_applicable,
        basis_bundle_include_ham_checked = true,
        basis_bundle_ham_payload_available,
        basis_bundle_ham_payload_status =
            basis_bundle_ham_payload_available ?
            (
                isnothing(ham_bundle_adapter) ?
                :available :
                :available_private_writer_adapter
            ) :
            :absent_for_fixed_block,
        pure_operator_payload_available,
        status,
        required_builder_contract =
            :white_lindsey_low_order_fixed_block_density_density_builder,
        ham_operator_payload_status,
        ham_interaction_status,
        ham_bundle_export_status,
        missing_builder,
        supplement_required_paths_policy = :diagnostic_only_not_benchmark_route,
        full_ham_export_ready = basis_bundle_ham_payload_available,
        private_writer_adapter_used = !isnothing(ham_bundle_adapter),
        private_payload_candidate_status =
            isnothing(ham_bundle_adapter) ? nothing : ham_bundle_adapter.candidate.status,
    )
end

function _pqs_source_box_route_driver_one_center_materializer_probe(
    config;
    probe_route_configured_one_center_materializer::Bool = false,
    white_lindsey_expansion = nothing,
)
    if !probe_route_configured_one_center_materializer
        return (;
            object_kind = :route_configured_one_center_materializer_probe,
            requested = false,
            status = :not_requested,
            materialized = false,
            route_configured_shellization_consumed = false,
            blocker = nothing,
            materialization = nothing,
            error_message = nothing,
        )
    elseif config.system_classification != :one_center
        return (;
            object_kind = :route_configured_one_center_materializer_probe,
            requested = true,
            status = :blocked_not_one_center,
            materialized = false,
            route_configured_shellization_consumed = false,
            blocker = :route_config_not_one_center,
            materialization = nothing,
            error_message = nothing,
        )
    elseif config.route_family != :white_lindsey_low_order
        return (;
            object_kind = :route_configured_one_center_materializer_probe,
            requested = true,
            status = :blocked_not_white_lindsey_low_order,
            materialized = false,
            route_configured_shellization_consumed = false,
            blocker = :route_config_not_white_lindsey_low_order,
            materialization = nothing,
            error_message = nothing,
        )
    elseif !config.materializer_options_ready
        return (;
            object_kind = :route_configured_one_center_materializer_probe,
            requested = true,
            status = :blocked_missing_materializer_options,
            materialized = false,
            route_configured_shellization_consumed = false,
            blocker = config.materializer_option_blocker,
            materialization = nothing,
            error_message =
                "missing materializer options: $(config.missing_materializer_options)",
        )
    end

    expansion =
        isnothing(white_lindsey_expansion) ?
        coulomb_gaussian_expansion(doacc = false) :
        white_lindsey_expansion
    try
        materialization =
            _cartesian_shellization_route_materialize_one_center_low_order(
                config;
                expansion,
                gausslet_backend = config.materializer_backend_requested,
                d = config.materializer_d_requested,
                nside = config.materializer_nside_requested,
                reference_spacing = config.materializer_reference_spacing_requested,
                tail_spacing = config.materializer_tail_spacing_requested,
            )
        return (;
            object_kind = :route_configured_one_center_materializer_probe,
            requested = true,
            status = materialization.status,
            materialized = true,
            route_configured_shellization_consumed =
                materialization.route_configured_shellization_consumed,
            blocker = nothing,
            materialization,
            error_message = nothing,
        )
    catch error
        error isa ArgumentError || rethrow()
        return (;
            object_kind = :route_configured_one_center_materializer_probe,
            requested = true,
            status = :blocked_materializer_precondition,
            materialized = false,
            route_configured_shellization_consumed = false,
            blocker = :materializer_precondition_failed,
            materialization = nothing,
            error_message = sprint(showerror, error),
        )
    end
end

function _pqs_source_box_route_driver_diatomic_materializer_probe(
    config;
    probe_route_configured_diatomic_materializer::Bool = false,
    route_materializer_payload = nothing,
    white_lindsey_expansion = nothing,
    shared_shell_layer_policy = nothing,
    packet_kernel = nothing,
)
    parent_qw_basis_object =
        isnothing(route_materializer_payload) ?
        nothing :
        route_materializer_payload.parent_qw_basis_object
    parent_axis_bundle_object =
        isnothing(route_materializer_payload) ?
        nothing :
        route_materializer_payload.parent_axis_bundle_object
    axis_bundle_backend =
        isnothing(route_materializer_payload) ?
        nothing :
        route_materializer_payload.axis_bundle_backend
    parent_qw_basis_object_handoff_available = !isnothing(parent_qw_basis_object)
    parent_axis_bundle_object_handoff_available =
        !isnothing(parent_axis_bundle_object)
    axis_bundle_backend_handoff_available = !isnothing(axis_bundle_backend)

    if !probe_route_configured_diatomic_materializer
        return (;
            object_kind = :route_configured_diatomic_materializer_probe,
            requested = false,
            status = :not_requested,
            materialized = false,
            route_configured_shellization_consumed = false,
            blocker = nothing,
            missing_contract = (),
            parent_qw_basis_object_handoff_available,
            parent_axis_bundle_object_handoff_available,
            axis_bundle_backend_handoff_available,
            materialization = nothing,
            error_message = nothing,
        )
    elseif config.system_classification != :bond_aligned_diatomic
        return (;
            object_kind = :route_configured_diatomic_materializer_probe,
            requested = true,
            status = :blocked_not_bond_aligned_diatomic,
            materialized = false,
            route_configured_shellization_consumed = false,
            blocker = :route_config_not_bond_aligned_diatomic,
            missing_contract = (),
            parent_qw_basis_object_handoff_available,
            parent_axis_bundle_object_handoff_available,
            axis_bundle_backend_handoff_available,
            materialization = nothing,
            error_message = nothing,
        )
    elseif config.route_family != :white_lindsey_low_order
        return (;
            object_kind = :route_configured_diatomic_materializer_probe,
            requested = true,
            status = :blocked_not_white_lindsey_low_order,
            materialized = false,
            route_configured_shellization_consumed = false,
            blocker = :route_config_not_white_lindsey_low_order,
            missing_contract = (),
            parent_qw_basis_object_handoff_available,
            parent_axis_bundle_object_handoff_available,
            axis_bundle_backend_handoff_available,
            materialization = nothing,
            error_message = nothing,
        )
    end

    try
        materialization =
            _cartesian_shellization_route_materialize_bond_aligned_diatomic(
                config;
                parent_qw_basis_object,
                parent_axis_bundle_object,
                expansion = white_lindsey_expansion,
                shared_shell_layer_policy,
                packet_kernel,
                axis_bundle_backend,
            )
        return (;
            object_kind = :route_configured_diatomic_materializer_probe,
            requested = true,
            status = materialization.status,
            materialized = materialization.materialized,
            route_configured_shellization_consumed =
                materialization.route_configured_shellization_consumed,
            blocker = materialization.blocker,
            missing_contract = materialization.missing_contract,
            parent_qw_basis_object_handoff_available,
            parent_axis_bundle_object_handoff_available,
            axis_bundle_backend_handoff_available,
            materialization,
            error_message = nothing,
        )
    catch error
        error isa ArgumentError || rethrow()
        return (;
            object_kind = :route_configured_diatomic_materializer_probe,
            requested = true,
            status = :blocked_materializer_precondition,
            materialized = false,
            route_configured_shellization_consumed = false,
            blocker = :materializer_precondition_failed,
            missing_contract = (),
            parent_qw_basis_object_handoff_available,
            parent_axis_bundle_object_handoff_available,
            axis_bundle_backend_handoff_available,
            materialization = nothing,
            error_message = sprint(showerror, error),
        )
    end
end

function _pqs_source_box_route_driver_materializer_payload_property(
    payload,
    field::Symbol,
)
    return !isnothing(payload) && hasproperty(payload, field) ?
           getproperty(payload, field) :
           nothing
end

function _pqs_source_box_route_driver_diatomic_atom_growth_materializer_probe(
    config;
    probe_route_configured_diatomic_atom_growth_materializer::Bool = false,
    route_materializer_payload = nothing,
    white_lindsey_expansion = nothing,
    packet_kernel = nothing,
)
    parent_qw_basis_object =
        _pqs_source_box_route_driver_materializer_payload_property(
            route_materializer_payload,
            :parent_qw_basis_object,
        )
    parent_axis_bundle_object =
        _pqs_source_box_route_driver_materializer_payload_property(
            route_materializer_payload,
            :parent_axis_bundle_object,
        )
    axis_bundle_backend =
        _pqs_source_box_route_driver_materializer_payload_property(
            route_materializer_payload,
            :axis_bundle_backend,
        )
    parent_qw_basis_object_handoff_available = !isnothing(parent_qw_basis_object)
    parent_axis_bundle_object_handoff_available =
        !isnothing(parent_axis_bundle_object)
    axis_bundle_backend_handoff_available = !isnothing(axis_bundle_backend)

    blocked_probe(status, blocker; missing_contract = (), error_message = nothing) =
        (;
            object_kind =
                :route_configured_diatomic_atom_growth_materializer_probe,
            requested = probe_route_configured_diatomic_atom_growth_materializer,
            status,
            materialized = false,
            atom_growth_shellification_consumed = false,
            route_configured_diatomic_atom_growth_shellification_consumed = false,
            route_configured_shellization_consumed = false,
            route_configured_legacy_diatomic_source_consumed = false,
            blocker,
            missing_contract = Tuple(missing_contract),
            parent_qw_basis_object_handoff_available,
            parent_axis_bundle_object_handoff_available,
            axis_bundle_backend_handoff_available,
            axis_bundle_backend_handoff = axis_bundle_backend,
            materializer_options = nothing,
            construction_plan = nothing,
            scaffold = nothing,
            materialization = nothing,
            sequence_available = false,
            retained_dimension = nothing,
            support_count = nothing,
            coverage_audit = nothing,
            coverage_complete = false,
            atom_growth_construction_plan_authority = false,
            active_source_authority = false,
            active_source_oracle_comparison_run = false,
            route_behavior_changed = false,
            shellification_plan_path_used = false,
            calls_shellification_plan_materializer = false,
            fixed_block_available = false,
            fixed_block_status = :not_requested,
            basis_adapter = nothing,
            basis_adapter_summary = nothing,
            basis_adapter_status = :not_requested,
            basis_adapter_blocker = nothing,
            final_integral_weights_status = :not_requested,
            ham_adapter = nothing,
            ham_adapter_summary = nothing,
            ham_adapter_status = :not_requested,
            ham_adapter_blocker = nothing,
            error_message,
        )

    if !probe_route_configured_diatomic_atom_growth_materializer
        return blocked_probe(:not_requested, nothing)
    elseif config.system_classification != :bond_aligned_diatomic
        return blocked_probe(
            :blocked_not_bond_aligned_diatomic,
            :route_config_not_bond_aligned_diatomic,
        )
    elseif config.route_family != :white_lindsey_low_order
        return blocked_probe(
            :blocked_not_white_lindsey_low_order,
            :route_config_not_white_lindsey_low_order,
        )
    end

    missing_contract = Symbol[]
    !config.materializer_options_ready &&
        append!(missing_contract, config.missing_materializer_options)
    isnothing(parent_qw_basis_object) &&
        push!(missing_contract, :parent_qw_basis_object_handoff)
    isnothing(parent_axis_bundle_object) &&
        push!(missing_contract, :parent_axis_bundle_object_handoff)
    isnothing(axis_bundle_backend) &&
        push!(missing_contract, :axis_bundle_backend_provenance)
    if !isnothing(axis_bundle_backend) &&
       axis_bundle_backend != config.materializer_backend_requested
        push!(missing_contract, :axis_bundle_backend_mismatches_materializer_backend)
    end
    missing_contract = Tuple(unique(missing_contract))
    if !isempty(missing_contract)
        return blocked_probe(
            :blocked_missing_atom_growth_materializer_contract,
            :pending_route_configured_bond_aligned_diatomic_atom_growth_materializer_contract;
            missing_contract,
        )
    end

    expansion =
        isnothing(white_lindsey_expansion) ?
        coulomb_gaussian_expansion(doacc = false) :
        white_lindsey_expansion
    consumed_packet_kernel = isnothing(packet_kernel) ? :factorized_direct : packet_kernel
    materializer_options = (;
        nside = config.materializer_nside_requested,
        d = config.materializer_d_requested,
        reference_spacing = config.materializer_reference_spacing_requested,
        tail_spacing = config.materializer_tail_spacing_requested,
        gausslet_backend = config.materializer_backend_requested,
        axis_bundle_backend,
        packet_kernel = consumed_packet_kernel,
        term_coefficients_source = :coulomb_expansion_coefficients,
        plan_authority = :atom_growth,
    )

    try
        nside = config.materializer_nside_requested
        anatomy = _nested_bond_aligned_diatomic_atom_growth_anatomy(
            parent_qw_basis_object,
            parent_axis_bundle_object;
            bond_axis = config.bond_axis,
            protected_atom_side_count = nside,
        )
        construction_plan =
            _nested_bond_aligned_diatomic_atom_growth_construction_plan(anatomy)
        retention = _nested_resolve_complete_shell_retention(nside)
        protect_rows =
            _nested_diatomic_resolve_core_near_nucleus_protect_rows(:auto, nside)
        scaffold =
            _cartesian_shellification_plan_atom_growth_complete_rectangular_low_order(
                construction_plan,
                parent_axis_bundle_object;
                nside,
                child_retention_policy = retention,
                shared_retention_policy = retention,
                reference_fudge_factor = 1.2,
                core_near_nucleus_protect_rows = protect_rows,
                shared_shell_angular_resolution_scale = 1.4,
                route_family = config.route_family,
            )
        materialization =
            _cartesian_materialize_atom_growth_complete_rectangular_shellification_low_order(
                scaffold,
                parent_qw_basis_object,
                parent_axis_bundle_object;
                term_coefficients = Float64.(expansion.coefficients),
                packet_kernel = consumed_packet_kernel,
            )
        sequence = materialization.sequence
        sequence_available = !isnothing(sequence)
        coverage_audit =
            sequence_available ?
            _nested_shell_sequence_contract_audit(
                sequence,
                Tuple(length.(construction_plan.anatomy.recipe.parent_box)),
            ) :
            nothing
        retained_dimension =
            sequence_available ? size(sequence.coefficient_matrix, 2) : nothing
        support_count =
            sequence_available ? length(sequence.support_indices) : nothing
        coverage_complete =
            !isnothing(coverage_audit) &&
            coverage_audit.full_parent_working_box &&
            coverage_audit.missing_row_count == 0 &&
            coverage_audit.ownership_unowned_row_count == 0 &&
            coverage_audit.ownership_multi_owned_row_count == 0
        materialized =
            materialization.status ==
            :materialized_supported_complete_rectangular_low_order
        basis_adapter =
            materialized ?
            _pqs_source_box_route_driver_diatomic_atom_growth_basis_adapter(
                materialization,
                construction_plan,
                scaffold,
                parent_qw_basis_object,
                config.materializer_backend_requested,
            ) :
            nothing
        basis_adapter_summary =
            isnothing(basis_adapter) ?
            nothing :
            _pqs_source_box_route_driver_route_configured_diatomic_basis_adapter_summary(
                basis_adapter,
            )
        basis_adapter_available =
            !isnothing(basis_adapter) &&
            basis_adapter.status ==
            :available_route_configured_diatomic_atom_growth_basis_adapter
        fixed_block_available =
            !isnothing(basis_adapter) && !isnothing(basis_adapter.fixed_block)
        ham_adapter =
            basis_adapter_available ?
            _pqs_source_box_route_driver_route_configured_diatomic_ham_adapter(
                basis_adapter,
                expansion;
                gausslet_backend = config.materializer_backend_requested,
                interaction_treatment = :ggt_nearest,
            ) :
            nothing
        ham_adapter_summary =
            isnothing(ham_adapter) ?
            nothing :
            _pqs_source_box_route_driver_route_configured_diatomic_ham_adapter_summary(
                ham_adapter,
            )
        ham_adapter_blocker =
            isnothing(ham_adapter) ?
            (
                basis_adapter_available ?
                :atom_growth_ham_operator_adapter_contract :
                isnothing(basis_adapter) ?
                :atom_growth_fixed_block_adapter_contract :
                basis_adapter.blocker
            ) :
            ham_adapter.blocker == :pending_route_configured_diatomic_ham_adapter_contract ?
            :atom_growth_ham_operator_adapter_contract :
            ham_adapter.blocker
        return (;
            object_kind =
                :route_configured_diatomic_atom_growth_materializer_probe,
            requested = true,
            status =
                materialized ?
                :materialized_route_configured_bond_aligned_diatomic_atom_growth_shellization :
                materialization.status,
            materialized,
            atom_growth_shellification_consumed = materialized,
            route_configured_diatomic_atom_growth_shellification_consumed =
                materialized,
            route_configured_shellization_consumed = materialized,
            route_configured_legacy_diatomic_source_consumed = false,
            blocker = materialized ? nothing : materialization.blocked_reason,
            missing_contract = (),
            parent_qw_basis_object_handoff_available,
            parent_axis_bundle_object_handoff_available,
            axis_bundle_backend_handoff_available,
            axis_bundle_backend_handoff = axis_bundle_backend,
            materializer_options,
            construction_plan,
            scaffold,
            materialization,
            sequence_available,
            retained_dimension,
            support_count,
            coverage_audit,
            coverage_complete,
            atom_growth_construction_plan_authority =
                scaffold.diagnostics.atom_growth_construction_plan_authority,
            active_source_authority =
                scaffold.diagnostics.active_source_authority ||
                materialization.active_source_authority,
            active_source_oracle_comparison_run =
                scaffold.diagnostics.active_source_oracle_comparison_run,
            route_behavior_changed =
                scaffold.diagnostics.materialization_behavior_changed ||
                materialization.route_behavior_changed,
            shellification_plan_path_used = true,
            calls_shellification_plan_materializer = true,
            fixed_block_available,
            fixed_block_status =
                fixed_block_available ?
                basis_adapter.fixed_block_status :
                isnothing(basis_adapter) ?
                :not_checked_missing_atom_growth_sequence :
                basis_adapter.fixed_block_status,
            basis_adapter,
            basis_adapter_summary,
            basis_adapter_status =
                isnothing(basis_adapter) ? :not_checked_missing_atom_growth_sequence :
                basis_adapter.status,
            basis_adapter_blocker =
                isnothing(basis_adapter) ? :atom_growth_fixed_block_adapter_contract :
                basis_adapter.blocker,
            final_integral_weights_status =
                isnothing(basis_adapter) ? :not_checked_missing_atom_growth_sequence :
                basis_adapter.final_integral_weights_status,
            ham_adapter,
            ham_adapter_summary,
            ham_adapter_status =
                isnothing(ham_adapter) ?
                :blocked_missing_route_configured_diatomic_atom_growth_basis_adapter :
                ham_adapter.status,
            ham_adapter_blocker,
            error_message = nothing,
        )
    catch error
        error isa ArgumentError || rethrow()
        return blocked_probe(
            :blocked_atom_growth_materializer_precondition,
            :atom_growth_materializer_precondition_failed;
            error_message = sprint(showerror, error),
        )
    end
end

function _pqs_source_box_route_driver_route_configured_one_center_report(
    materialization,
)
    fixture = materialization.fixture
    inventory = fixture.inventory
    route_units = _white_lindsey_low_order_materialized_seed_route_units(fixture)
    operator_inventory =
        _white_lindsey_low_order_materialized_seed_operator_inventory(fixture)
    operator_pairs_materialized =
        route_units.operator_pairs_materialized ||
        operator_inventory.operator_pairs_materialized
    shellization_summary = materialization.shellization_summary

    return (;
        object_kind = :white_lindsey_low_order_route_configured_one_center_report,
        route_family = :white_lindsey_low_order,
        status = :private_development_route_configured,
        private_development_only = true,
        materialization,
        fixture,
        inventory,
        route_units,
        operator_inventory,
        shellization_summary,
        shellization_summary_available = true,
        shellization_source = :route_configured_one_center_low_order,
        route_configured_shellization_consumed = true,
        materialized_shellization_stage = shellization_summary.shellization_stage,
        seed_materialization_status = :not_seed_route_configured_materialization,
        packet_kernel = fixture.packet_kernel,
        retained_dimension = route_units.retained_dimension,
        operator_pairs_materialized,
        electron_electron_materialized = operator_inventory.electron_electron_materialized,
        weight_semantics = :retained_basis_integral_weights,
    )
end

function _pqs_source_box_route_driver_route_configured_diatomic_basis_adapter(
    materialization,
)
    if isnothing(materialization) || !materialization.materialized
        return (;
            object_kind = :route_configured_diatomic_basis_adapter,
            status = :blocked_missing_route_configured_diatomic_materialization,
            private_development_only = true,
            fixed_block = nothing,
            representation = nothing,
            final_integral_weights = nothing,
            retained_dimension = nothing,
            basis_metadata = nothing,
            grouping = nothing,
            label_status = :not_checked_no_materialization,
            grouping_status = :not_checked_no_materialization,
            final_integral_weights_status = :not_checked_no_materialization,
            missing_fields = (:route_configured_diatomic_materialization,),
            blocker = :missing_route_configured_diatomic_materialization,
        )
    elseif isnothing(materialization.source)
        return (;
            object_kind = :route_configured_diatomic_basis_adapter,
            status = :blocked_missing_route_configured_diatomic_source,
            private_development_only = true,
            fixed_block = nothing,
            representation = nothing,
            final_integral_weights = nothing,
            retained_dimension = materialization.retained_dimension,
            basis_metadata = nothing,
            grouping = nothing,
            label_status = :not_checked_no_source,
            grouping_status = :not_checked_no_source,
            final_integral_weights_status = :not_checked_no_source,
            missing_fields = (:route_configured_diatomic_source,),
            blocker = :missing_route_configured_diatomic_source,
        )
    end

    try
        source = materialization.source
        fixed_block = _nested_fixed_block(source)
        representation = basis_representation(fixed_block)
        final_integral_weights =
            _cartesian_bundle_integral_weights(fixed_block, representation)
        retained_dimension = representation.metadata.final_dimension
        labels = representation.metadata.basis_labels
        centers = representation.metadata.basis_centers
        weights_ready =
            length(final_integral_weights) == retained_dimension &&
            all(isfinite, final_integral_weights) &&
            all(>(0.0), final_integral_weights)
        labels_ready =
            length(labels) == retained_dimension &&
            all(!isempty, labels)
        centers_ready =
            size(centers) == (retained_dimension, 3) &&
            all(isfinite, centers)
        grouping = (;
            source_kind = :route_configured_bond_aligned_diatomic_source,
            shell_kind = representation.metadata.route_metadata.shell_kind,
            working_box_profile =
                representation.metadata.route_metadata.working_box_profile,
            support_count = representation.metadata.route_metadata.support_count,
            child_sequence_count = length(source.child_sequences),
            shared_shell_layer_count = length(source.shared_shell_layers),
            child_column_ranges = Tuple(source.child_column_ranges),
            midpoint_slab_column_range = source.midpoint_slab_column_range,
        )
        grouping_ready =
            grouping.child_sequence_count == length(grouping.child_column_ranges) &&
            grouping.support_count == length(fixed_block.support_indices)
        missing_fields = Symbol[]
        weights_ready || push!(
            missing_fields,
            :route_configured_diatomic_final_weight_contract,
        )
        labels_ready || push!(missing_fields, :route_configured_diatomic_basis_labels)
        centers_ready || push!(missing_fields, :route_configured_diatomic_basis_centers)
        grouping_ready || push!(missing_fields, :route_configured_diatomic_grouping)
        missing_fields = Tuple(missing_fields)

        return (;
            object_kind = :route_configured_diatomic_basis_adapter,
            status =
                isempty(missing_fields) ?
                :available_route_configured_diatomic_basis_adapter :
                :blocked_route_configured_diatomic_basis_adapter_contract,
            private_development_only = true,
            fixed_block,
            representation,
            final_integral_weights,
            retained_dimension,
            basis_metadata = (;
                basis_kind = representation.metadata.basis_kind,
                parent_kind = representation.metadata.parent_kind,
                axis_sharing = representation.metadata.axis_sharing,
                parent_axis_counts = representation.metadata.parent_axis_counts,
                parent_dimension = representation.metadata.parent_dimension,
                final_dimension = representation.metadata.final_dimension,
                label_count = length(labels),
                center_count = size(centers, 1),
            ),
            grouping,
            label_status =
                labels_ready ? :available_route_configured_diatomic_basis_labels :
                :blocked_route_configured_diatomic_basis_labels,
            grouping_status =
                grouping_ready ? :available_route_configured_diatomic_grouping :
                :blocked_route_configured_diatomic_grouping,
            final_integral_weights_status =
                weights_ready ?
                :available_retained_basis_integral_weights :
                :pending_route_configured_diatomic_final_weight_contract,
            missing_fields,
            blocker =
                isempty(missing_fields) ?
                nothing :
                :pending_route_configured_diatomic_basis_adapter_contract,
        )
    catch error
        error isa ArgumentError || rethrow()
        return (;
            object_kind = :route_configured_diatomic_basis_adapter,
            status = :blocked_route_configured_diatomic_basis_adapter_precondition,
            private_development_only = true,
            fixed_block = nothing,
            representation = nothing,
            final_integral_weights = nothing,
            retained_dimension = materialization.retained_dimension,
            basis_metadata = nothing,
            grouping = nothing,
            label_status = :not_checked_basis_adapter_precondition,
            grouping_status = :not_checked_basis_adapter_precondition,
            final_integral_weights_status =
                :not_checked_basis_adapter_precondition,
            missing_fields = (:route_configured_diatomic_basis_adapter_precondition,),
            blocker = :route_configured_diatomic_basis_adapter_precondition_failed,
            error_message = sprint(showerror, error),
        )
    end
end

function _pqs_source_box_route_driver_route_configured_diatomic_basis_adapter_summary(
    adapter,
)
    return (;
        object_kind = adapter.object_kind,
        status = adapter.status,
        private_development_only = adapter.private_development_only,
        retained_dimension = adapter.retained_dimension,
        basis_metadata = adapter.basis_metadata,
        grouping = adapter.grouping,
        label_status = adapter.label_status,
        grouping_status = adapter.grouping_status,
        final_integral_weights_status = adapter.final_integral_weights_status,
        final_integral_weight_count =
            isnothing(adapter.final_integral_weights) ?
            nothing :
            length(adapter.final_integral_weights),
        missing_fields = adapter.missing_fields,
        blocker = adapter.blocker,
    )
end

function _pqs_source_box_route_driver_diatomic_atom_growth_basis_adapter(
    materialization,
    construction_plan,
    scaffold,
    parent_qw_basis_object,
    gausslet_backend,
)
    if isnothing(materialization) || isnothing(materialization.sequence)
        return (;
            object_kind = :route_configured_diatomic_atom_growth_basis_adapter,
            status = :blocked_missing_route_configured_diatomic_atom_growth_sequence,
            private_development_only = true,
            fixed_block = nothing,
            fixed_block_status = :not_checked_missing_atom_growth_sequence,
            representation = nothing,
            final_integral_weights = nothing,
            retained_dimension = nothing,
            basis_metadata = nothing,
            grouping = nothing,
            label_status = :not_checked_missing_atom_growth_sequence,
            grouping_status = :not_checked_missing_atom_growth_sequence,
            final_integral_weights_status =
                :not_checked_missing_atom_growth_sequence,
            missing_fields = (:route_configured_diatomic_atom_growth_sequence,),
            blocker = :atom_growth_fixed_block_adapter_contract,
        )
    elseif isnothing(parent_qw_basis_object)
        return (;
            object_kind = :route_configured_diatomic_atom_growth_basis_adapter,
            status = :blocked_missing_atom_growth_parent_basis,
            private_development_only = true,
            fixed_block = nothing,
            fixed_block_status = :not_checked_missing_parent_basis,
            representation = nothing,
            final_integral_weights = nothing,
            retained_dimension = size(materialization.sequence.coefficient_matrix, 2),
            basis_metadata = nothing,
            grouping = nothing,
            label_status = :not_checked_missing_parent_basis,
            grouping_status = :not_checked_missing_parent_basis,
            final_integral_weights_status = :not_checked_missing_parent_basis,
            missing_fields = (:parent_qw_basis_object_handoff,),
            blocker = :atom_growth_fixed_block_adapter_contract,
        )
    end

    try
        sequence = materialization.sequence
        fixed_block = _nested_fixed_block(
            sequence,
            parent_qw_basis_object,
            gausslet_backend,
        )
        representation = basis_representation(fixed_block)
        final_integral_weights =
            _cartesian_bundle_integral_weights(fixed_block, representation)
        retained_dimension = representation.metadata.final_dimension
        labels = representation.metadata.basis_labels
        centers = representation.metadata.basis_centers
        route_metadata = representation.metadata.route_metadata
        weights_ready =
            length(final_integral_weights) == retained_dimension &&
            all(isfinite, final_integral_weights) &&
            all(>(0.0), final_integral_weights)
        labels_ready =
            length(labels) == retained_dimension &&
            all(!isempty, labels)
        centers_ready =
            size(centers) == (retained_dimension, 3) &&
            all(isfinite, centers)
        grouping = (;
            source_kind = :bond_aligned_diatomic_atom_growth_construction_plan,
            shell_kind = route_metadata.shell_kind,
            working_box_profile = route_metadata.working_box_profile,
            support_count = route_metadata.support_count,
            construction_region_order = Tuple(construction_plan.region_order),
            assembly_core_order = scaffold.assembly_core_order,
            assembly_shell_order = scaffold.assembly_shell_order,
            outer_mismatch_boundary_slab_set_count =
                length(scaffold.outer_mismatch_boundary_slab_sets),
            child_column_ranges = Tuple(materialization.assembly.child_column_ranges),
            contact_cap_column_range =
                materialization.assembly.contact_cap_column_range,
            shared_shell_column_ranges =
                Tuple(materialization.assembly.shared_shell_column_ranges),
        )
        grouping_ready =
            grouping.support_count == length(fixed_block.support_indices) &&
            grouping.support_count == length(sequence.support_indices)
        missing_fields = Symbol[]
        weights_ready ||
            push!(missing_fields, :atom_growth_final_integral_weight_contract)
        labels_ready ||
            push!(missing_fields, :atom_growth_basis_representation_contract)
        centers_ready ||
            push!(missing_fields, :atom_growth_basis_representation_contract)
        grouping_ready ||
            push!(missing_fields, :atom_growth_basis_representation_contract)
        missing_fields = Tuple(unique(missing_fields))
        blocker =
            isempty(missing_fields) ?
            nothing :
            in(:atom_growth_final_integral_weight_contract, missing_fields) ?
            :atom_growth_final_integral_weight_contract :
            :atom_growth_basis_representation_contract

        return (;
            object_kind = :route_configured_diatomic_atom_growth_basis_adapter,
            status =
                isempty(missing_fields) ?
                :available_route_configured_diatomic_atom_growth_basis_adapter :
                :blocked_route_configured_diatomic_atom_growth_basis_adapter_contract,
            private_development_only = true,
            fixed_block,
            fixed_block_status = :available_route_configured_diatomic_atom_growth_fixed_block,
            representation,
            final_integral_weights,
            retained_dimension,
            basis_metadata = (;
                basis_kind = representation.metadata.basis_kind,
                parent_kind = representation.metadata.parent_kind,
                axis_sharing = representation.metadata.axis_sharing,
                parent_axis_counts = representation.metadata.parent_axis_counts,
                parent_dimension = representation.metadata.parent_dimension,
                final_dimension = representation.metadata.final_dimension,
                label_count = length(labels),
                center_count = size(centers, 1),
            ),
            grouping,
            label_status =
                labels_ready ?
                :available_route_configured_diatomic_atom_growth_basis_labels :
                :blocked_atom_growth_basis_representation_contract,
            grouping_status =
                grouping_ready ?
                :available_route_configured_diatomic_atom_growth_grouping :
                :blocked_atom_growth_basis_representation_contract,
            final_integral_weights_status =
                weights_ready ?
                :available_retained_basis_integral_weights :
                :pending_atom_growth_final_integral_weight_contract,
            missing_fields,
            blocker,
        )
    catch error
        error isa ArgumentError || rethrow()
        return (;
            object_kind = :route_configured_diatomic_atom_growth_basis_adapter,
            status = :blocked_atom_growth_fixed_block_adapter_precondition,
            private_development_only = true,
            fixed_block = nothing,
            fixed_block_status = :blocked_atom_growth_fixed_block_adapter_precondition,
            representation = nothing,
            final_integral_weights = nothing,
            retained_dimension =
                !isnothing(materialization) && !isnothing(materialization.sequence) ?
                size(materialization.sequence.coefficient_matrix, 2) :
                nothing,
            basis_metadata = nothing,
            grouping = nothing,
            label_status = :not_checked_fixed_block_adapter_precondition,
            grouping_status = :not_checked_fixed_block_adapter_precondition,
            final_integral_weights_status =
                :not_checked_fixed_block_adapter_precondition,
            missing_fields = (:atom_growth_fixed_block_adapter_contract,),
            blocker = :atom_growth_fixed_block_adapter_contract,
            error_message = sprint(showerror, error),
        )
    end
end

function _pqs_source_box_route_driver_route_configured_diatomic_ham_adapter(
    basis_adapter,
    expansion;
    gausslet_backend,
    interaction_treatment::Symbol = :ggt_nearest,
    nuclear_term_storage::Symbol = :by_center,
)
    mwg_ida_treatments = (:mwg, :ida, :mwg_ida, :ida_mwg)
    available_basis_adapter_statuses = (
        :available_route_configured_diatomic_basis_adapter,
        :available_route_configured_diatomic_atom_growth_basis_adapter,
    )
    if !(basis_adapter.status in available_basis_adapter_statuses)
        return (;
            object_kind = :route_configured_diatomic_ham_adapter,
            status = :blocked_missing_route_configured_diatomic_basis_adapter,
            private_development_only = true,
            operators = nothing,
            retained_dimension = basis_adapter.retained_dimension,
            matrix_sizes = nothing,
            final_integral_weights_status = basis_adapter.final_integral_weights_status,
            nuclear_metadata_status = :not_checked_missing_basis_adapter,
            operator_payload_status = :not_checked_missing_basis_adapter,
            interaction_status = :not_checked_missing_basis_adapter,
            interaction_treatment_requested = interaction_treatment,
            interaction_treatment,
            missing_fields = (:route_configured_diatomic_basis_adapter,),
            blocker = :missing_route_configured_diatomic_basis_adapter,
        )
    elseif interaction_treatment in mwg_ida_treatments
        return (;
            object_kind = :route_configured_diatomic_ham_adapter,
            status =
                :blocked_route_configured_diatomic_ham_interaction_treatment,
            private_development_only = true,
            operators = nothing,
            retained_dimension = basis_adapter.retained_dimension,
            matrix_sizes = nothing,
            final_integral_weights_status = basis_adapter.final_integral_weights_status,
            nuclear_metadata_status = :not_checked_blocked_interaction_treatment,
            operator_payload_status = :blocked_route_configured_diatomic_operator_payload,
            interaction_status =
                :pending_route_configured_diatomic_mwg_operator_support,
            interaction_treatment_requested = interaction_treatment,
            interaction_treatment,
            gausslet_backend,
            nuclear_charges = nothing,
            nuclear_term_storage = nothing,
            nuclear_one_body_by_center_count = 0,
            missing_fields = (:pending_route_configured_diatomic_mwg_operator_support,),
            blocker = :pending_route_configured_diatomic_mwg_operator_support,
        )
    elseif isnothing(expansion)
        return (;
            object_kind = :route_configured_diatomic_ham_adapter,
            status = :blocked_missing_route_configured_diatomic_expansion,
            private_development_only = true,
            operators = nothing,
            retained_dimension = basis_adapter.retained_dimension,
            matrix_sizes = nothing,
            final_integral_weights_status = basis_adapter.final_integral_weights_status,
            nuclear_metadata_status = :not_checked_missing_expansion,
            operator_payload_status = :not_checked_missing_expansion,
            interaction_status = :not_checked_missing_expansion,
            interaction_treatment_requested = interaction_treatment,
            interaction_treatment,
            missing_fields = (:coulomb_expansion,),
            blocker = :missing_route_configured_diatomic_expansion,
        )
    end

    try
        operators = ordinary_cartesian_qiu_white_operators(
            basis_adapter.fixed_block;
            expansion,
            gausslet_backend,
            interaction_treatment,
            nuclear_term_storage,
        )
        retained_dimension = size(operators.overlap, 1)
        matrix_sizes = (;
            overlap = size(operators.overlap),
            one_body_hamiltonian = size(operators.one_body_hamiltonian),
            interaction_matrix = size(operators.interaction_matrix),
        )
        expected_matrix_size = (retained_dimension, retained_dimension)
        final_integral_weights =
            _cartesian_bundle_integral_weights(operators, basis_adapter.representation)
        matrix_size_ready =
            all(matrix_size -> matrix_size == expected_matrix_size, values(matrix_sizes))
        finite_ready =
            all(isfinite, operators.overlap) &&
            all(isfinite, operators.one_body_hamiltonian) &&
            all(isfinite, operators.interaction_matrix)
        weights_ready =
            length(final_integral_weights) == retained_dimension &&
            all(isfinite, final_integral_weights) &&
            all(>(0.0), final_integral_weights)
        parent_nuclear_charges = basis_adapter.fixed_block.parent_basis.nuclear_charges
        expected_nuclear_center_count = length(parent_nuclear_charges)
        nuclear_metadata_ready =
            !isnothing(operators.nuclear_charges) &&
            length(operators.nuclear_charges) == expected_nuclear_center_count &&
            operators.nuclear_term_storage == :by_center &&
            !isnothing(operators.nuclear_one_body_by_center) &&
            length(operators.nuclear_one_body_by_center) ==
            length(operators.nuclear_charges)
        missing_fields = Symbol[]
        matrix_size_ready || push!(missing_fields, :route_configured_diatomic_ham_matrix_sizes)
        finite_ready || push!(missing_fields, :route_configured_diatomic_ham_finite_matrices)
        weights_ready || push!(
            missing_fields,
            :route_configured_diatomic_ham_final_integral_weights,
        )
        nuclear_metadata_ready || push!(
            missing_fields,
            :route_configured_diatomic_ham_nuclear_metadata,
        )
        missing_fields = Tuple(missing_fields)

        return (;
            object_kind = :route_configured_diatomic_ham_adapter,
            status =
                isempty(missing_fields) ?
                :available_route_configured_diatomic_ham_adapter :
                :blocked_route_configured_diatomic_ham_adapter_contract,
            private_development_only = true,
            operators,
            retained_dimension,
            matrix_sizes,
            final_integral_weights_status =
                weights_ready ?
                :available_retained_basis_integral_weights :
                :pending_route_configured_diatomic_ham_final_integral_weights,
            nuclear_metadata_status =
                nuclear_metadata_ready ?
                :available_route_configured_diatomic_nuclear_metadata :
                :pending_route_configured_diatomic_nuclear_metadata,
            operator_payload_status =
                matrix_size_ready && finite_ready ?
                :available_route_configured_diatomic_operator_payload :
                :pending_route_configured_diatomic_operator_payload,
            interaction_status =
                matrix_size_ready && finite_ready ?
                :available_route_configured_diatomic_density_density_interaction_matrix :
                :pending_route_configured_diatomic_density_density_builder,
            interaction_treatment_requested = interaction_treatment,
            interaction_treatment = operators.interaction_treatment,
            gausslet_backend = operators.gausslet_backend,
            nuclear_charges = Tuple(operators.nuclear_charges),
            nuclear_term_storage = operators.nuclear_term_storage,
            nuclear_one_body_by_center_count =
                isnothing(operators.nuclear_one_body_by_center) ?
                0 :
                length(operators.nuclear_one_body_by_center),
            missing_fields,
            blocker =
                isempty(missing_fields) ?
                nothing :
                :pending_route_configured_diatomic_ham_adapter_contract,
        )
    catch error
        error isa ArgumentError || rethrow()
        return (;
            object_kind = :route_configured_diatomic_ham_adapter,
            status = :blocked_route_configured_diatomic_ham_adapter_precondition,
            private_development_only = true,
            operators = nothing,
            retained_dimension = basis_adapter.retained_dimension,
            matrix_sizes = nothing,
            final_integral_weights_status =
                :not_checked_ham_adapter_precondition,
            nuclear_metadata_status = :not_checked_ham_adapter_precondition,
            operator_payload_status = :not_checked_ham_adapter_precondition,
            interaction_status = :not_checked_ham_adapter_precondition,
            interaction_treatment_requested = interaction_treatment,
            interaction_treatment,
            missing_fields = (:route_configured_diatomic_ham_adapter_precondition,),
            blocker = :route_configured_diatomic_ham_adapter_precondition_failed,
            error_message = sprint(showerror, error),
        )
    end
end

function _pqs_source_box_route_driver_route_configured_diatomic_ham_adapter_summary(
    adapter,
)
    return (;
        object_kind = adapter.object_kind,
        status = adapter.status,
        private_development_only = adapter.private_development_only,
        retained_dimension = adapter.retained_dimension,
        matrix_sizes = adapter.matrix_sizes,
        final_integral_weights_status = adapter.final_integral_weights_status,
        nuclear_metadata_status = adapter.nuclear_metadata_status,
        operator_payload_status = adapter.operator_payload_status,
        interaction_status = adapter.interaction_status,
        interaction_treatment_requested =
            hasproperty(adapter, :interaction_treatment_requested) ?
            adapter.interaction_treatment_requested :
            nothing,
        interaction_treatment =
            hasproperty(adapter, :interaction_treatment) ?
            adapter.interaction_treatment :
            nothing,
        gausslet_backend =
            hasproperty(adapter, :gausslet_backend) ? adapter.gausslet_backend : nothing,
        nuclear_charges =
            hasproperty(adapter, :nuclear_charges) ? adapter.nuclear_charges : nothing,
        nuclear_term_storage =
            hasproperty(adapter, :nuclear_term_storage) ?
            adapter.nuclear_term_storage :
            nothing,
        nuclear_one_body_by_center_count =
            hasproperty(adapter, :nuclear_one_body_by_center_count) ?
            adapter.nuclear_one_body_by_center_count :
            nothing,
        missing_fields = adapter.missing_fields,
        blocker = adapter.blocker,
    )
end

function _pqs_source_box_route_driver_low_order_shellization_policy(
    requested_policy,
    probe_route_configured_diatomic_atom_growth_materializer::Bool,
)
    supported_policies = (
        :legacy_diatomic_source,
        :atom_growth_complete_rectangular,
        :terminal_cartesian_shellification_geometry,
    )
    explicit_policy_requested = !isnothing(requested_policy)
    resolved_policy =
        explicit_policy_requested ?
        requested_policy :
        probe_route_configured_diatomic_atom_growth_materializer ?
        :atom_growth_complete_rectangular :
        :legacy_diatomic_source
    policy_source =
        explicit_policy_requested ?
        :explicit_low_order_shellization_policy :
        probe_route_configured_diatomic_atom_growth_materializer ?
        :probe_route_configured_diatomic_atom_growth_materializer_alias :
        :default_legacy_diatomic_source
    supported = resolved_policy in supported_policies
    conflict =
        explicit_policy_requested &&
        probe_route_configured_diatomic_atom_growth_materializer &&
        resolved_policy != :atom_growth_complete_rectangular
    status =
        !supported ?
        :blocked_unsupported_low_order_shellization_policy :
        conflict ?
        :blocked_conflicting_low_order_shellization_policy :
        :available_low_order_shellization_policy
    blocker =
        status == :available_low_order_shellization_policy ?
        nothing :
        status

    return (;
        low_order_shellization_policy_requested = requested_policy,
        low_order_shellization_policy_resolved = resolved_policy,
        low_order_shellization_policy_source = policy_source,
        low_order_shellization_policy_status = status,
        low_order_shellization_policy_blocker = blocker,
    )
end

function _pqs_source_box_route_driver_materialization(
    report;
    materialize_route::Bool = false,
    probe_route_configured_one_center_materializer::Bool = false,
    probe_route_configured_diatomic_atom_growth_materializer::Bool = false,
    low_order_shellization_policy = nothing,
    save_basis_artifact::Bool = false,
    save_ham_artifact::Bool = false,
    basisfile::AbstractString = "cartesian_nesting_route_driver_basis_bundle.jld2",
    hamfile::AbstractString = "cartesian_nesting_route_driver_ham_bundle.jld2",
    materializer_backend = nothing,
    materializer_nside = nothing,
    route_configured_diatomic_ham_interaction_treatment::Symbol = :ggt_nearest,
    white_lindsey_expansion = nothing,
    white_lindsey_Z = nothing,
)
    route_family = report.route_family
    route_configured_shellization_request =
        _cartesian_shellization_route_configured_request(
            report;
            materializer_backend,
            materializer_nside,
        )
    route_configured_shellization_request_status =
        route_configured_shellization_request.status
    route_configured_system_classification =
        route_configured_shellization_request.system_classification
    route_configured_system_classification_status =
        route_configured_shellization_request.system_classification_status
    route_configured_bond_axis = route_configured_shellization_request.bond_axis
    route_configured_shellization_plan =
        _cartesian_shellization_route_planning_stub(route_configured_shellization_request)
    route_configured_shellization_plan_status =
        route_configured_shellization_plan.status
    route_configured_shellization_planning_status =
        route_configured_shellization_plan.planning_status
    route_configured_shellization_planning_family =
        route_configured_shellization_plan.planning_family
    route_configured_midpoint_slab_status =
        route_configured_shellization_plan.midpoint_slab_status
    route_configured_shellization_helper_map =
        _cartesian_shellization_route_planning_helper_map(route_configured_shellization_plan)
    route_configured_shellization_helper_map_status =
        route_configured_shellization_helper_map.status
    route_configured_primary_planned_helper =
        route_configured_shellization_helper_map.primary_planned_helper
    route_configured_missing_input_count =
        route_configured_shellization_helper_map.missing_input_count
    route_configured_helper_map_blocker =
        route_configured_shellization_helper_map.blocker
    route_configured_input_readiness =
        _cartesian_shellization_route_materializer_input_readiness(
            route_configured_shellization_request,
            route_configured_shellization_plan,
            route_configured_shellization_helper_map,
        )
    route_configured_input_readiness_status = route_configured_input_readiness.status
    route_configured_available_fact_count =
        route_configured_input_readiness.available_fact_count
    route_configured_materializer_missing_input_count =
        route_configured_input_readiness.missing_input_count
    route_configured_input_readiness_blocker =
        route_configured_input_readiness.blocker
    route_configured_materializer_config =
        _cartesian_shellization_route_materializer_config(
            route_configured_shellization_request,
            route_configured_shellization_plan,
            route_configured_shellization_helper_map,
            route_configured_input_readiness,
        )
    route_configured_materializer_config_status =
        route_configured_materializer_config.status
    route_configured_materializer_config_planning_family =
        route_configured_materializer_config.planning_family
    route_configured_materializer_config_pending_input_count =
        route_configured_materializer_config.pending_input_count
    low_order_shellization_policy_contract =
        _pqs_source_box_route_driver_low_order_shellization_policy(
            low_order_shellization_policy,
            probe_route_configured_diatomic_atom_growth_materializer,
        )
    low_order_shellization_policy_requested =
        low_order_shellization_policy_contract.low_order_shellization_policy_requested
    low_order_shellization_policy_resolved =
        low_order_shellization_policy_contract.low_order_shellization_policy_resolved
    low_order_shellization_policy_source =
        low_order_shellization_policy_contract.low_order_shellization_policy_source
    low_order_shellization_policy_status =
        low_order_shellization_policy_contract.low_order_shellization_policy_status
    low_order_shellization_policy_blocker =
        low_order_shellization_policy_contract.low_order_shellization_policy_blocker
    route_configured_one_center_materializer_requested =
        probe_route_configured_one_center_materializer ||
        (
            materialize_route &&
            route_family == :white_lindsey_low_order &&
            route_configured_system_classification == :one_center
        )
    route_configured_one_center_materializer_probe =
        _pqs_source_box_route_driver_one_center_materializer_probe(
            route_configured_materializer_config;
            probe_route_configured_one_center_materializer =
                route_configured_one_center_materializer_requested,
            white_lindsey_expansion,
        )
    route_configured_one_center_materializer_probe_requested =
        route_configured_one_center_materializer_probe.requested
    route_configured_one_center_materializer_probe_status =
        route_configured_one_center_materializer_probe.status
    route_configured_one_center_materializer_probe_materialized =
        route_configured_one_center_materializer_probe.materialized
    route_configured_one_center_materializer_probe_consumed =
        route_configured_one_center_materializer_probe.route_configured_shellization_consumed
    route_configured_one_center_materializer_probe_blocker =
        route_configured_one_center_materializer_probe.blocker
    route_configured_diatomic_materializer_requested =
        materialize_route &&
        route_family == :white_lindsey_low_order &&
        route_configured_system_classification == :bond_aligned_diatomic &&
        low_order_shellization_policy_status == :available_low_order_shellization_policy &&
        low_order_shellization_policy_resolved == :legacy_diatomic_source
    route_materializer_payload =
        hasproperty(report, :route_materializer_payload) ?
        report.route_materializer_payload :
        nothing
    route_configured_diatomic_shared_shell_layer_policy =
        route_configured_diatomic_materializer_requested &&
        route_configured_materializer_config.materializer_backend_requested ==
        :pgdg_localized_experimental ?
        :endcap_panel_owned :
        nothing
    route_configured_diatomic_packet_kernel =
        route_configured_diatomic_materializer_requested ?
        :factorized_direct :
        nothing
    route_configured_diatomic_policy_source =
        isnothing(route_configured_diatomic_shared_shell_layer_policy) ?
        nothing :
        :existing_endcap_panel_owned_pgdg_route
    route_configured_diatomic_materializer_probe =
        _pqs_source_box_route_driver_diatomic_materializer_probe(
            route_configured_materializer_config;
            probe_route_configured_diatomic_materializer =
                route_configured_diatomic_materializer_requested,
            route_materializer_payload,
            white_lindsey_expansion,
            shared_shell_layer_policy =
                route_configured_diatomic_shared_shell_layer_policy,
            packet_kernel = route_configured_diatomic_packet_kernel,
        )
    route_configured_diatomic_materializer_probe_requested =
        route_configured_diatomic_materializer_probe.requested
    route_configured_diatomic_materializer_probe_status =
        route_configured_diatomic_materializer_probe.status
    route_configured_diatomic_materializer_probe_materialized =
        route_configured_diatomic_materializer_probe.materialized
    route_configured_diatomic_materializer_probe_consumed =
        route_configured_diatomic_materializer_probe.route_configured_shellization_consumed
    route_configured_diatomic_materializer_probe_blocker =
        route_configured_diatomic_materializer_probe.blocker
    route_configured_diatomic_materializer_missing_contract =
        route_configured_diatomic_materializer_probe.missing_contract
    route_configured_diatomic_materializer_payload_available =
        !isnothing(route_materializer_payload)
    route_configured_diatomic_parent_qw_basis_object_handoff_available =
        route_configured_diatomic_materializer_probe.parent_qw_basis_object_handoff_available
    route_configured_diatomic_parent_axis_bundle_object_handoff_available =
        route_configured_diatomic_materializer_probe.parent_axis_bundle_object_handoff_available
    route_configured_diatomic_axis_bundle_backend_handoff_available =
        route_configured_diatomic_materializer_probe.axis_bundle_backend_handoff_available
    route_configured_diatomic_axis_bundle_backend_handoff =
        isnothing(route_materializer_payload) ?
        nothing :
        route_materializer_payload.axis_bundle_backend
    route_configured_diatomic_seed_fallback =
        route_configured_diatomic_materializer_probe_requested &&
        !route_configured_diatomic_materializer_probe_consumed
    route_configured_diatomic_atom_growth_materializer_requested =
        materialize_route &&
        route_family == :white_lindsey_low_order &&
        route_configured_system_classification == :bond_aligned_diatomic &&
        low_order_shellization_policy_status == :available_low_order_shellization_policy &&
        low_order_shellization_policy_resolved == :atom_growth_complete_rectangular
    route_configured_diatomic_atom_growth_materializer_probe_requested_input =
        low_order_shellization_policy_status == :available_low_order_shellization_policy &&
        (
            probe_route_configured_diatomic_atom_growth_materializer ||
            route_configured_diatomic_atom_growth_materializer_requested
        )
    route_configured_diatomic_atom_growth_materializer_probe =
        _pqs_source_box_route_driver_diatomic_atom_growth_materializer_probe(
            route_configured_materializer_config;
            probe_route_configured_diatomic_atom_growth_materializer =
                route_configured_diatomic_atom_growth_materializer_probe_requested_input,
            route_materializer_payload,
            white_lindsey_expansion,
            packet_kernel = route_configured_diatomic_packet_kernel,
        )
    route_configured_diatomic_atom_growth_materializer_probe_requested =
        route_configured_diatomic_atom_growth_materializer_probe.requested
    route_configured_diatomic_atom_growth_materializer_probe_status =
        route_configured_diatomic_atom_growth_materializer_probe.status
    route_configured_diatomic_atom_growth_materializer_probe_materialized =
        route_configured_diatomic_atom_growth_materializer_probe.materialized
    route_configured_diatomic_atom_growth_materializer_probe_consumed =
        route_configured_diatomic_atom_growth_materializer_probe.atom_growth_shellification_consumed
    route_configured_diatomic_atom_growth_shellification_consumed =
        route_configured_diatomic_atom_growth_materializer_probe.route_configured_diatomic_atom_growth_shellification_consumed
    route_configured_diatomic_atom_growth_materializer_probe_blocker =
        route_configured_diatomic_atom_growth_materializer_probe.blocker
    route_configured_diatomic_atom_growth_materializer_missing_contract =
        route_configured_diatomic_atom_growth_materializer_probe.missing_contract
    route_configured_diatomic_atom_growth_sequence_available =
        route_configured_diatomic_atom_growth_materializer_probe.sequence_available
    route_configured_diatomic_atom_growth_retained_dimension =
        route_configured_diatomic_atom_growth_materializer_probe.retained_dimension
    route_configured_diatomic_atom_growth_support_count =
        route_configured_diatomic_atom_growth_materializer_probe.support_count
    route_configured_diatomic_atom_growth_coverage_complete =
        route_configured_diatomic_atom_growth_materializer_probe.coverage_complete
    route_configured_diatomic_atom_growth_active_source_authority =
        route_configured_diatomic_atom_growth_materializer_probe.active_source_authority
    route_configured_diatomic_atom_growth_plan_authority =
        route_configured_diatomic_atom_growth_materializer_probe.atom_growth_construction_plan_authority
    route_configured_diatomic_atom_growth_calls_shellification_plan_materializer =
        route_configured_diatomic_atom_growth_materializer_probe.calls_shellification_plan_materializer
    route_configured_diatomic_atom_growth_fixed_block_available =
        route_configured_diatomic_atom_growth_materializer_probe.fixed_block_available
    route_configured_diatomic_atom_growth_fixed_block_status =
        route_configured_diatomic_atom_growth_materializer_probe.fixed_block_status
    route_configured_diatomic_atom_growth_basis_adapter_status =
        route_configured_diatomic_atom_growth_materializer_probe.basis_adapter_status
    route_configured_diatomic_atom_growth_basis_adapter_blocker =
        route_configured_diatomic_atom_growth_materializer_probe.basis_adapter_blocker
    route_configured_diatomic_atom_growth_final_integral_weights_status =
        route_configured_diatomic_atom_growth_materializer_probe.final_integral_weights_status
    route_configured_diatomic_atom_growth_ham_adapter_status =
        route_configured_diatomic_atom_growth_materializer_probe.ham_adapter_status
    route_configured_diatomic_atom_growth_ham_adapter_blocker =
        route_configured_diatomic_atom_growth_materializer_probe.ham_adapter_blocker
    route_configured_materializer_backend_requested =
        route_configured_materializer_config.materializer_backend_requested
    route_configured_materializer_backend_source =
        route_configured_materializer_config.materializer_backend_source
    route_configured_materializer_backend_status =
        route_configured_materializer_config.materializer_backend_status
    route_configured_materializer_d_requested =
        route_configured_materializer_config.materializer_d_requested
    route_configured_materializer_d_source =
        route_configured_materializer_config.materializer_d_source
    route_configured_materializer_d_status =
        route_configured_materializer_config.materializer_d_status
    route_configured_materializer_nside_requested =
        route_configured_materializer_config.materializer_nside_requested
    route_configured_materializer_nside_source =
        route_configured_materializer_config.materializer_nside_source
    route_configured_materializer_nside_status =
        route_configured_materializer_config.materializer_nside_status
    route_configured_materializer_reference_spacing_requested =
        route_configured_materializer_config.materializer_reference_spacing_requested
    route_configured_materializer_reference_spacing_source =
        route_configured_materializer_config.materializer_reference_spacing_source
    route_configured_materializer_reference_spacing_status =
        route_configured_materializer_config.materializer_reference_spacing_status
    route_configured_materializer_tail_spacing_requested =
        route_configured_materializer_config.materializer_tail_spacing_requested
    route_configured_materializer_tail_spacing_source =
        route_configured_materializer_config.materializer_tail_spacing_source
    route_configured_materializer_tail_spacing_status =
        route_configured_materializer_config.materializer_tail_spacing_status
    route_configured_materializer_options_ready =
        route_configured_materializer_config.materializer_options_ready
    route_configured_materializer_missing_options =
        route_configured_materializer_config.missing_materializer_options
    route_configured_materializer_option_blocker =
        route_configured_materializer_config.materializer_option_blocker
    route_configured_materializer_consumed_options =
        route_configured_one_center_materializer_probe_materialized ?
        route_configured_one_center_materializer_probe.materialization.materializer_options :
        route_configured_diatomic_materializer_probe_materialized ?
        route_configured_diatomic_materializer_probe.materialization.materializer_options :
        route_configured_diatomic_atom_growth_materializer_probe_materialized ?
        route_configured_diatomic_atom_growth_materializer_probe.materializer_options :
        nothing
    route_configured_materializer_backend_consumed =
        isnothing(route_configured_materializer_consumed_options) ?
        nothing :
        route_configured_materializer_consumed_options.gausslet_backend
    route_configured_materializer_d_consumed =
        isnothing(route_configured_materializer_consumed_options) ?
        nothing :
        route_configured_materializer_consumed_options.d
    route_configured_materializer_nside_consumed =
        isnothing(route_configured_materializer_consumed_options) ?
        nothing :
        route_configured_materializer_consumed_options.nside
    route_configured_materializer_reference_spacing_consumed =
        isnothing(route_configured_materializer_consumed_options) ?
        nothing :
        route_configured_materializer_consumed_options.reference_spacing
    route_configured_materializer_tail_spacing_consumed =
        isnothing(route_configured_materializer_consumed_options) ?
        nothing :
        route_configured_materializer_consumed_options.tail_spacing
    route_configured_materializer_contract = (;
        route_configured_materializer_backend_requested,
        route_configured_materializer_backend_source,
        route_configured_materializer_backend_status,
        route_configured_materializer_backend_consumed,
        route_configured_materializer_d_requested,
        route_configured_materializer_d_source,
        route_configured_materializer_d_status,
        route_configured_materializer_d_consumed,
        route_configured_materializer_nside_requested,
        route_configured_materializer_nside_source,
        route_configured_materializer_nside_status,
        route_configured_materializer_nside_consumed,
        route_configured_materializer_reference_spacing_requested,
        route_configured_materializer_reference_spacing_source,
        route_configured_materializer_reference_spacing_status,
        route_configured_materializer_reference_spacing_consumed,
        route_configured_materializer_tail_spacing_requested,
        route_configured_materializer_tail_spacing_source,
        route_configured_materializer_tail_spacing_status,
        route_configured_materializer_tail_spacing_consumed,
        route_configured_materializer_options_ready,
        route_configured_materializer_missing_options,
        route_configured_materializer_option_blocker,
        low_order_shellization_policy_requested,
        low_order_shellization_policy_resolved,
        low_order_shellization_policy_source,
        low_order_shellization_policy_status,
        low_order_shellization_policy_blocker,
    )
    route_configured_diatomic_materializer_contract = (;
        route_configured_diatomic_materializer_probe,
        route_configured_diatomic_materializer_probe_requested,
        route_configured_diatomic_materializer_probe_status,
        route_configured_diatomic_materializer_probe_materialized,
        route_configured_diatomic_materializer_probe_consumed,
        route_configured_diatomic_materializer_probe_blocker,
        route_configured_diatomic_materializer_missing_contract,
        route_configured_diatomic_materializer_payload_available,
        route_configured_diatomic_parent_qw_basis_object_handoff_available,
        route_configured_diatomic_parent_axis_bundle_object_handoff_available,
        route_configured_diatomic_axis_bundle_backend_handoff_available,
        route_configured_diatomic_axis_bundle_backend_handoff,
        route_configured_diatomic_shared_shell_layer_policy,
        route_configured_diatomic_packet_kernel,
        route_configured_diatomic_policy_source,
        route_configured_diatomic_seed_fallback,
    )
    route_configured_diatomic_atom_growth_materializer_contract = (;
        route_configured_diatomic_atom_growth_materializer_probe_requested,
        route_configured_diatomic_atom_growth_materializer_probe_status,
        route_configured_diatomic_atom_growth_materializer_probe_materialized,
        route_configured_diatomic_atom_growth_materializer_probe_consumed,
        route_configured_diatomic_atom_growth_shellification_consumed,
        route_configured_diatomic_atom_growth_materializer_probe_blocker,
        route_configured_diatomic_atom_growth_materializer_missing_contract,
        route_configured_diatomic_atom_growth_sequence_available,
        route_configured_diatomic_atom_growth_retained_dimension,
        route_configured_diatomic_atom_growth_support_count,
        route_configured_diatomic_atom_growth_coverage_complete,
        route_configured_diatomic_atom_growth_active_source_authority,
        route_configured_diatomic_atom_growth_plan_authority,
        route_configured_diatomic_atom_growth_calls_shellification_plan_materializer,
        route_configured_diatomic_atom_growth_fixed_block_available,
        route_configured_diatomic_atom_growth_fixed_block_status,
        route_configured_diatomic_atom_growth_basis_adapter_status,
        route_configured_diatomic_atom_growth_basis_adapter_blocker,
        route_configured_diatomic_atom_growth_final_integral_weights_status,
        route_configured_diatomic_atom_growth_ham_adapter_status,
        route_configured_diatomic_atom_growth_ham_adapter_blocker,
    )

    if !materialize_route
        return (;
            object_kind = :cartesian_nesting_route_driver_materialization,
            route_family,
            private_development_only = true,
            materialize_route_requested = false,
            save_basis_artifact_requested = save_basis_artifact,
            save_ham_artifact_requested = save_ham_artifact,
            route_configured_diatomic_ham_interaction_treatment_requested = nothing,
            route_configured_diatomic_ham_interaction_treatment_consumed = nothing,
            route_configured_diatomic_ham_interaction_treatment_status =
                :not_applicable,
            status = :not_requested_metadata_only,
            materialized_report = nothing,
            materialized_report_kind = nothing,
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
            shellization_summary = nothing,
            shellization_summary_available = false,
            shellization_source =
                route_family == :white_lindsey_low_order ?
                :white_lindsey_one_center_seed_not_materialized :
                nothing,
            route_configured_shellization_consumed = false,
            route_configured_legacy_diatomic_source_consumed = false,
            materialized_shellization_stage = :not_checked_metadata_only,
            seed_materialization_status =
                route_family == :white_lindsey_low_order ?
                :not_requested_seed_materialization :
                :not_applicable,
            retained_dimension = report.retained_dimension,
            final_integral_weights_status = :not_checked_metadata_only,
            one_body_operator_status = :not_checked_metadata_only,
            basis_bundle_export_status = :not_requested,
            basis_artifact_status =
                save_basis_artifact ?
                :not_written_materialization_not_requested :
                :not_requested,
            basis_artifact_written = false,
            basisfile,
            basis_artifact_path = nothing,
            basis_export_blocker =
                save_basis_artifact ? :materialize_route_false : nothing,
            ham_preflight_status = :not_checked_metadata_only,
            ham_missing_builder = nothing,
            ham_operator_payload_status = :not_checked_metadata_only,
            ham_interaction_status = :not_checked_metadata_only,
            ham_bundle_export_status = :not_requested,
            ham_artifact_status =
                save_ham_artifact ?
                :not_written_materialization_not_requested :
                :not_requested,
            ham_artifact_written = false,
            hamfile,
            ham_export_blocker =
                save_ham_artifact ? :materialize_route_false : nothing,
            ham_preflight = nothing,
            pqs_materialization_status =
                route_family == :pqs_source_box ?
                :pending_source_box_retained_route :
                :not_applicable,
        )
    end

    if low_order_shellization_policy_status != :available_low_order_shellization_policy
        return (;
            object_kind = :cartesian_nesting_route_driver_materialization,
            route_family,
            private_development_only = true,
            materialize_route_requested = true,
            save_basis_artifact_requested = save_basis_artifact,
            save_ham_artifact_requested = save_ham_artifact,
            route_configured_diatomic_ham_interaction_treatment_requested = nothing,
            route_configured_diatomic_ham_interaction_treatment_consumed = nothing,
            route_configured_diatomic_ham_interaction_treatment_status =
                :not_applicable,
            status = low_order_shellization_policy_status,
            materialized_report = nothing,
            materialized_report_kind = nothing,
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
            shellization_summary = nothing,
            shellization_summary_available = false,
            shellization_source = :blocked_low_order_shellization_policy,
            route_configured_shellization_consumed = false,
            route_configured_legacy_diatomic_source_consumed = false,
            materialized_shellization_stage = :blocked_low_order_shellization_policy,
            seed_materialization_status = :not_applicable,
            retained_dimension = report.retained_dimension,
            final_integral_weights_status = :not_checked_low_order_policy_blocked,
            one_body_operator_status = :not_checked_low_order_policy_blocked,
            basis_bundle_export_status = :not_requested,
            basis_artifact_status =
                save_basis_artifact ?
                :not_written_invalid_low_order_shellization_policy :
                :not_requested,
            basis_artifact_written = false,
            basisfile,
            basis_artifact_path = nothing,
            basis_export_blocker =
                save_basis_artifact ? low_order_shellization_policy_blocker : nothing,
            ham_preflight_status = :not_checked_low_order_policy_blocked,
            ham_missing_builder =
                save_ham_artifact ? low_order_shellization_policy_blocker : nothing,
            ham_operator_payload_status = :not_checked_low_order_policy_blocked,
            ham_interaction_status = :not_checked_low_order_policy_blocked,
            ham_bundle_export_status = :not_requested,
            ham_artifact_status =
                save_ham_artifact ?
                :not_written_invalid_low_order_shellization_policy :
                :not_requested,
            ham_artifact_written = false,
            hamfile,
            ham_export_blocker =
                save_ham_artifact ? low_order_shellization_policy_blocker : nothing,
            ham_preflight = nothing,
            pqs_materialization_status =
                route_family == :pqs_source_box ?
                :pending_source_box_retained_route :
                :not_applicable,
        )
    end

    if route_family == :white_lindsey_low_order
        route_configured_one_center_report_required =
            materialize_route &&
            route_configured_system_classification == :one_center
        if route_configured_one_center_report_required &&
           !route_configured_one_center_materializer_probe_materialized
            throw(
                ArgumentError(
                    "route-configured one-center materializer failed: " *
                    string(route_configured_one_center_materializer_probe_blocker) *
                    " " *
                    string(route_configured_one_center_materializer_probe.error_message),
                ),
            )
        end
        use_route_configured_one_center_report =
            route_configured_one_center_report_required &&
            route_configured_one_center_materializer_probe_materialized
        use_route_configured_diatomic_shellization =
            route_configured_diatomic_materializer_requested &&
            route_configured_diatomic_materializer_probe_consumed
        route_configured_diatomic_atom_growth_materialization_requested =
            route_configured_diatomic_atom_growth_materializer_probe_requested
        if route_configured_diatomic_atom_growth_materialization_requested
            atom_growth_materialization_context = (;
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
            )
            return _pqs_source_box_route_driver_diatomic_atom_growth_materialization(
                atom_growth_materialization_context,
            )
        end
        if use_route_configured_diatomic_shellization
            diatomic_materialization =
                route_configured_diatomic_materializer_probe.materialization
            diatomic_basis_adapter =
                _pqs_source_box_route_driver_route_configured_diatomic_basis_adapter(
                    diatomic_materialization,
                )
            diatomic_basis_adapter_summary =
                _pqs_source_box_route_driver_route_configured_diatomic_basis_adapter_summary(
                    diatomic_basis_adapter,
                )
            diatomic_basis_adapter_available =
                diatomic_basis_adapter.status ==
                :available_route_configured_diatomic_basis_adapter
            diatomic_ham_adapter =
                save_ham_artifact ?
                _pqs_source_box_route_driver_route_configured_diatomic_ham_adapter(
                    diatomic_basis_adapter,
                    white_lindsey_expansion;
                    gausslet_backend =
                        route_configured_materializer_backend_requested,
                    interaction_treatment =
                        route_configured_diatomic_ham_interaction_treatment,
                ) :
                nothing
            diatomic_ham_adapter_summary =
                isnothing(diatomic_ham_adapter) ?
                nothing :
                _pqs_source_box_route_driver_route_configured_diatomic_ham_adapter_summary(
                    diatomic_ham_adapter,
                )
            diatomic_ham_adapter_available =
                !isnothing(diatomic_ham_adapter) &&
                diatomic_ham_adapter.status ==
                :available_route_configured_diatomic_ham_adapter
            diatomic_ham_interaction_treatment_consumed =
                diatomic_ham_adapter_available ?
                diatomic_ham_adapter.interaction_treatment :
                nothing
            diatomic_ham_interaction_treatment_status =
                diatomic_ham_adapter_available ?
                :available_route_configured_diatomic_ham_interaction_treatment :
                save_ham_artifact && !isnothing(diatomic_ham_adapter) ?
                diatomic_ham_adapter.interaction_status :
                :not_requested
            shellization_summary = diatomic_materialization.shellization_summary
            shellization_summary_available = !isnothing(shellization_summary)
            basis_artifact_status =
                save_basis_artifact ?
                (
                    diatomic_basis_adapter_available ?
                    :written_route_configured_diatomic_basis_only_bundle :
                    :not_written_route_configured_diatomic_basis_adapter_blocked
                ) :
                :not_requested
            ham_artifact_status =
                save_ham_artifact ?
                (
                    diatomic_ham_adapter_available ?
                    :written_route_configured_diatomic_ham_bundle :
                    :not_written_route_configured_diatomic_ham_adapter_blocked
                ) :
                :not_requested
            diatomic_ham_adapter_blocker =
                isnothing(diatomic_ham_adapter) ? nothing : diatomic_ham_adapter.blocker
            diatomic_ham_bundle_export_status =
                diatomic_ham_adapter_available ?
                :available_route_configured_diatomic_ham_bundle_payload :
                save_ham_artifact ?
                something(
                    diatomic_ham_adapter_blocker,
                    :pending_route_configured_diatomic_ham_export,
                ) :
                :not_requested
            ham_export_blocker =
                diatomic_ham_adapter_available || !save_ham_artifact ?
                nothing :
                something(
                    diatomic_ham_adapter_blocker,
                    :pending_route_configured_diatomic_ham_export,
                )
            basis_companion_ham_artifact_status =
                save_ham_artifact ?
                (
                    diatomic_ham_adapter_available ?
                    :companion_route_configured_diatomic_ham_artifact_ready :
                    something(
                        diatomic_ham_adapter_blocker,
                        :pending_route_configured_diatomic_ham_export,
                    )
                ) :
                :not_requested
            basis_companion_ham_export_status =
                save_ham_artifact ? diatomic_ham_bundle_export_status : :not_requested
            basis_companion_ham_export_blocker =
                save_ham_artifact ? ham_export_blocker : nothing
            basis_artifact_written = false
            if save_basis_artifact && diatomic_basis_adapter_available
                write_cartesian_basis_bundle_jld2(
                    basisfile,
                    diatomic_basis_adapter.fixed_block;
                    include_ham = false,
                    meta = (;
                        route_family,
                        route_kind = report.recipe_metadata.route_kind,
                        benchmark_role = report.recipe_metadata.benchmark_role,
                        materialized_report_kind = diatomic_materialization.object_kind,
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
                        route_configured_diatomic_basis_adapter_status =
                            diatomic_basis_adapter.status,
                        route_configured_diatomic_basis_adapter_retained_dimension =
                            diatomic_basis_adapter.retained_dimension,
                        route_configured_diatomic_basis_adapter_final_integral_weights_status =
                            diatomic_basis_adapter.final_integral_weights_status,
                        route_configured_diatomic_basis_adapter_label_status =
                            diatomic_basis_adapter.label_status,
                        route_configured_diatomic_basis_adapter_grouping_status =
                            diatomic_basis_adapter.grouping_status,
                        route_configured_diatomic_ham_interaction_treatment_requested =
                            route_configured_diatomic_ham_interaction_treatment,
                        route_configured_diatomic_ham_interaction_treatment_consumed =
                            diatomic_ham_interaction_treatment_consumed,
                        route_configured_diatomic_ham_interaction_treatment_status =
                            diatomic_ham_interaction_treatment_status,
                        shellization_summary_available,
                        shellization_source =
                            :route_configured_bond_aligned_diatomic_source,
                        route_configured_shellization_consumed = true,
                        route_configured_legacy_diatomic_source_consumed = true,
                        materialized_shellization_stage =
                            shellization_summary.shellization_stage,
                        seed_materialization_status =
                            :not_seed_route_configured_diatomic_shellization,
                        export_status = :basis_only,
                        basis_export_status =
                            :supported_route_configured_diatomic_basis_only_fixed_block,
                        ham_export_status =
                            :artifact_local_basis_only_no_ham_payload,
                        ham_export_blocker = nothing,
                        companion_ham_artifact_requested = save_ham_artifact,
                        companion_ham_artifact_status =
                            basis_companion_ham_artifact_status,
                        companion_ham_export_status =
                            basis_companion_ham_export_status,
                        companion_ham_export_blocker =
                            basis_companion_ham_export_blocker,
                        private_development_only = true,
                    ),
                )
                basis_artifact_written = true
            end
            ham_artifact_written = false
            if save_ham_artifact && diatomic_ham_adapter_available
                write_cartesian_basis_bundle_jld2(
                    hamfile,
                    diatomic_ham_adapter.operators;
                    include_ham = true,
                    meta = (;
                        route_family,
                        route_kind = report.recipe_metadata.route_kind,
                        benchmark_role = report.recipe_metadata.benchmark_role,
                        materialized_report_kind = diatomic_materialization.object_kind,
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
                        route_configured_diatomic_basis_adapter_status =
                            diatomic_basis_adapter.status,
                        route_configured_diatomic_basis_adapter_retained_dimension =
                            diatomic_basis_adapter.retained_dimension,
                        route_configured_diatomic_basis_adapter_final_integral_weights_status =
                            diatomic_basis_adapter.final_integral_weights_status,
                        route_configured_diatomic_ham_adapter_status =
                            diatomic_ham_adapter.status,
                        route_configured_diatomic_ham_adapter_operator_payload_status =
                            diatomic_ham_adapter.operator_payload_status,
                        route_configured_diatomic_ham_adapter_interaction_status =
                            diatomic_ham_adapter.interaction_status,
                        route_configured_diatomic_ham_adapter_nuclear_metadata_status =
                            diatomic_ham_adapter.nuclear_metadata_status,
                        route_configured_diatomic_ham_interaction_treatment_requested =
                            route_configured_diatomic_ham_interaction_treatment,
                        route_configured_diatomic_ham_interaction_treatment_consumed =
                            diatomic_ham_interaction_treatment_consumed,
                        route_configured_diatomic_ham_interaction_treatment_status =
                            diatomic_ham_interaction_treatment_status,
                        shellization_summary_available,
                        shellization_source =
                            :route_configured_bond_aligned_diatomic_source,
                        route_configured_shellization_consumed = true,
                        route_configured_legacy_diatomic_source_consumed = true,
                        materialized_shellization_stage =
                            shellization_summary.shellization_stage,
                        seed_materialization_status =
                            :not_seed_route_configured_diatomic_shellization,
                        export_status = :basis_and_ham,
                        basis_export_status =
                            :supported_route_configured_diatomic_basis_only_fixed_block,
                        ham_preflight_status =
                            :available_route_configured_diatomic_ham_adapter,
                        ham_operator_payload_status =
                            :available_route_configured_diatomic_operator_payload,
                        ham_interaction_status =
                            :available_route_configured_diatomic_density_density_interaction_matrix,
                        ham_export_status =
                            :available_route_configured_diatomic_ham_bundle_payload,
                        ham_export_blocker = nothing,
                        private_development_only = true,
                    ),
                )
                ham_artifact_written = true
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
                    diatomic_ham_interaction_treatment_consumed,
                route_configured_diatomic_ham_interaction_treatment_status =
                    diatomic_ham_interaction_treatment_status,
                status = :materialized_route_configured_diatomic_shellization_available,
                materialized_report = nothing,
                materialized_report_kind = diatomic_materialization.object_kind,
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
                    diatomic_basis_adapter_summary,
                route_configured_diatomic_ham_adapter_summary =
                    diatomic_ham_adapter_summary,
                shellization_summary,
                shellization_summary_available,
                shellization_source = :route_configured_bond_aligned_diatomic_source,
                route_configured_shellization_consumed = true,
                route_configured_legacy_diatomic_source_consumed = true,
                materialized_shellization_stage =
                    shellization_summary.shellization_stage,
                seed_materialization_status =
                    :not_seed_route_configured_diatomic_shellization,
                retained_dimension = diatomic_materialization.retained_dimension,
                final_integral_weights_status =
                    diatomic_basis_adapter.final_integral_weights_status,
                one_body_operator_status =
                    :pending_route_configured_diatomic_operator_inventory,
                basis_bundle_export_status =
                    diatomic_basis_adapter_available ?
                    :supported_route_configured_diatomic_basis_only_fixed_block :
                    :pending_route_configured_diatomic_basis_export,
                basis_artifact_status,
                basis_artifact_written,
                basisfile,
                basis_artifact_path = basis_artifact_written ? basisfile : nothing,
                basis_export_blocker =
                    diatomic_basis_adapter_available ?
                    nothing :
                    :pending_route_configured_diatomic_basis_adapter_contract,
                ham_preflight_status =
                    diatomic_ham_adapter_available ?
                    :available_route_configured_diatomic_ham_adapter :
                    save_ham_artifact ?
                    :blocked_route_configured_diatomic_ham_adapter :
                    :not_requested,
                ham_missing_builder =
                    diatomic_ham_adapter_available ?
                    nothing :
                    save_ham_artifact ?
                    something(
                        diatomic_ham_adapter_blocker,
                        :pending_route_configured_diatomic_ham_adapter,
                    ) :
                    nothing,
                ham_operator_payload_status =
                    isnothing(diatomic_ham_adapter) ?
                    :not_requested :
                    diatomic_ham_adapter.operator_payload_status,
                ham_interaction_status =
                    isnothing(diatomic_ham_adapter) ?
                    :not_requested :
                    diatomic_ham_adapter.interaction_status,
                ham_bundle_export_status =
                    diatomic_ham_bundle_export_status,
                ham_artifact_status,
                ham_artifact_written,
                hamfile,
                ham_export_blocker,
                ham_preflight = nothing,
                pqs_materialization_status = :not_applicable,
            )
        end
        materialized_report =
            use_route_configured_one_center_report ?
            _pqs_source_box_route_driver_route_configured_one_center_report(
                route_configured_one_center_materializer_probe.materialization,
            ) :
            _white_lindsey_low_order_materialized_seed_report()
        basis_export_status =
            use_route_configured_one_center_report ?
            :supported_route_configured_one_center_basis_only_fixed_block :
            :supported_basis_only_fixed_block
        shellization_summary = materialized_report.shellization_summary
        shellization_summary_available = materialized_report.shellization_summary_available
        shellization_source = materialized_report.shellization_source
        route_configured_shellization_consumed =
            materialized_report.route_configured_shellization_consumed
        materialized_shellization_stage = materialized_report.materialized_shellization_stage
        seed_materialization_status = materialized_report.seed_materialization_status
        ham_bundle_adapter = nothing
        if save_ham_artifact
            isnothing(white_lindsey_expansion) && throw(
                ArgumentError(
                    "White-Lindsey Ham artifact export requires explicit white_lindsey_expansion",
                ),
            )
            isnothing(white_lindsey_Z) && throw(
                ArgumentError("White-Lindsey Ham artifact export requires explicit white_lindsey_Z"),
            )
            ham_bundle_adapter =
                _white_lindsey_low_order_materialized_seed_ham_bundle_adapter(
                    materialized_report;
                    expansion = white_lindsey_expansion,
                    Z = white_lindsey_Z,
                )
        end
        ham_preflight =
            isnothing(ham_bundle_adapter) ?
            _pqs_source_box_route_driver_white_lindsey_ham_preflight(materialized_report) :
            _pqs_source_box_route_driver_white_lindsey_ham_preflight(
                materialized_report;
                ham_bundle_adapter = ham_bundle_adapter,
            )
        ham_operator_payload_status = ham_preflight.ham_operator_payload_status
        ham_interaction_status = ham_preflight.ham_interaction_status
        ham_bundle_export_status = ham_preflight.ham_bundle_export_status
        ham_export_blocker = ham_preflight.missing_builder
        basis_artifact_written = false
        basis_artifact_status =
            save_basis_artifact ?
            (
                use_route_configured_one_center_report ?
                :written_route_configured_one_center_basis_only_bundle :
                :written_basis_only_bundle
            ) :
            :not_requested
        if save_basis_artifact
            write_cartesian_basis_bundle_jld2(
                basisfile,
                materialized_report.fixture.fixed_block;
                include_ham = false,
                meta = (;
                    route_family,
                    route_kind = report.recipe_metadata.route_kind,
                    benchmark_role = report.recipe_metadata.benchmark_role,
                    materialized_report_kind = materialized_report.object_kind,
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
                    shellization_summary_available,
                    shellization_source,
                    route_configured_shellization_consumed,
                    route_configured_legacy_diatomic_source_consumed = false,
                    materialized_shellization_stage,
                    seed_materialization_status,
                    export_status = :basis_only,
                    basis_export_status,
                    ham_preflight_status = ham_preflight.status,
                    ham_missing_builder = ham_preflight.missing_builder,
                    ham_operator_payload_status,
                    ham_interaction_status,
                    ham_export_status = ham_bundle_export_status,
                    ham_export_blocker,
                    private_development_only = true,
                ),
            )
            basis_artifact_written = true
        end
        ham_artifact_written = false
        ham_artifact_status =
            save_ham_artifact ?
            :not_written_private_white_lindsey_ham_adapter_not_ready :
            :not_requested
        if save_ham_artifact && ham_preflight.full_ham_export_ready
            write_cartesian_basis_bundle_jld2(
                hamfile,
                ham_bundle_adapter;
                include_ham = true,
                meta = (;
                    route_family,
                    route_kind = report.recipe_metadata.route_kind,
                    benchmark_role = report.recipe_metadata.benchmark_role,
                    materialized_report_kind = materialized_report.object_kind,
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
                    shellization_summary_available,
                    shellization_source,
                    route_configured_shellization_consumed,
                    route_configured_legacy_diatomic_source_consumed = false,
                    materialized_shellization_stage,
                    seed_materialization_status,
                    export_status = :basis_and_ham,
                    basis_export_status,
                    ham_preflight_status = ham_preflight.status,
                    ham_missing_builder = ham_preflight.missing_builder,
                    ham_operator_payload_status,
                    ham_interaction_status,
                    ham_export_status = ham_bundle_export_status,
                    ham_export_blocker,
                    private_development_only = true,
                    private_writer_adapter =
                        :_WhiteLindseyLowOrderHamBundleAdapter,
                    ham_payload_candidate_status =
                        ham_bundle_adapter.candidate.status,
                ),
            )
            ham_artifact_written = true
            ham_artifact_status =
                use_route_configured_one_center_report ?
                :written_route_configured_one_center_ham_bundle :
                :written_white_lindsey_low_order_ham_bundle
        end
        return (;
            object_kind = :cartesian_nesting_route_driver_materialization,
            route_family,
            private_development_only = true,
            materialize_route_requested = true,
            save_basis_artifact_requested = save_basis_artifact,
            save_ham_artifact_requested = save_ham_artifact,
            route_configured_diatomic_ham_interaction_treatment_requested = nothing,
            route_configured_diatomic_ham_interaction_treatment_consumed = nothing,
            route_configured_diatomic_ham_interaction_treatment_status =
                :not_applicable,
            status =
                use_route_configured_one_center_report ?
                :materialized_route_configured_one_center_report_available :
                :materialized_seed_report_available,
            materialized_report,
            materialized_report_kind = materialized_report.object_kind,
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
            shellization_summary,
            shellization_summary_available,
            shellization_source,
            route_configured_shellization_consumed,
            route_configured_legacy_diatomic_source_consumed = false,
            materialized_shellization_stage,
            seed_materialization_status,
            retained_dimension = materialized_report.retained_dimension,
            final_integral_weights_status =
                materialized_report.inventory.retained_basis_integral_weights_ready ?
                :available_retained_basis_integral_weights :
                :not_ready,
            one_body_operator_status =
                materialized_report.operator_inventory.all_finite ?
                :materialized_finite_one_body_inventory :
                :not_ready,
            basis_bundle_export_status = basis_export_status,
            basis_artifact_status,
            basis_artifact_written,
            basisfile,
            basis_artifact_path = basis_artifact_written ? basisfile : nothing,
            basis_export_blocker = nothing,
            ham_preflight_status = ham_preflight.status,
            ham_missing_builder = ham_preflight.missing_builder,
            ham_operator_payload_status,
            ham_interaction_status,
            ham_bundle_export_status,
            ham_artifact_status,
            ham_artifact_written,
            hamfile,
            ham_export_blocker,
            ham_preflight,
            pqs_materialization_status = :not_applicable,
        )
    end

    ham_export_blocker = :pending_source_box_retained_route
    return (;
        object_kind = :cartesian_nesting_route_driver_materialization,
        route_family,
        private_development_only = true,
        materialize_route_requested = true,
        save_basis_artifact_requested = save_basis_artifact,
        save_ham_artifact_requested = save_ham_artifact,
        route_configured_diatomic_ham_interaction_treatment_requested = nothing,
        route_configured_diatomic_ham_interaction_treatment_consumed = nothing,
        route_configured_diatomic_ham_interaction_treatment_status = :not_applicable,
        status = :pending_source_box_retained_route,
        materialized_report = nothing,
        materialized_report_kind = nothing,
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
        shellization_summary = nothing,
        shellization_summary_available = false,
        shellization_source = :pending_source_box_route_shellization,
        route_configured_shellization_consumed = false,
        route_configured_legacy_diatomic_source_consumed = false,
        materialized_shellization_stage = :pending_source_box_retained_route,
        seed_materialization_status = :not_applicable,
        retained_dimension = report.retained_dimension,
        final_integral_weights_status = :pending_final_ida_weights,
        one_body_operator_status = :pending_source_box_retained_blocks,
        basis_bundle_export_status = :pending_final_retained_basis,
        basis_artifact_status =
            save_basis_artifact ?
            :not_written_pending_final_retained_basis :
            :not_requested,
        basis_artifact_written = false,
        basisfile,
        basis_artifact_path = nothing,
        basis_export_blocker =
            save_basis_artifact ? :pending_final_retained_basis : nothing,
        ham_preflight_status = :not_applicable_to_pqs_source_box_route,
        ham_missing_builder = :pending_source_box_retained_route,
        ham_operator_payload_status = :pending_source_box_retained_operator_payload,
        ham_interaction_status = :pending_source_box_retained_density_density_blocks,
        ham_bundle_export_status = :pending_source_box_retained_route,
        ham_artifact_status =
            save_ham_artifact ?
            :not_written_pending_source_box_retained_route :
            :not_requested,
        ham_artifact_written = false,
        hamfile,
        ham_export_blocker,
        ham_preflight = nothing,
        pqs_materialization_status = :pending_source_box_retained_route,
    )
end

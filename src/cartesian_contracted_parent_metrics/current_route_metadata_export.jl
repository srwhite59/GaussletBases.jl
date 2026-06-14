function _pqs_inventory_wrapped_staged_unit(
    unit::_CartesianNestedProductStagedByCenterUnit3D;
    category::Symbol,
    support_source_semantics::Symbol,
    safe_term_capability::Symbol,
    active_representation_stage::Symbol,
    source_helper::Symbol,
    raw_box_auxiliary_metadata = nothing,
)
    source_axis_center_vectors =
        hasproperty(unit.provenance, :source_axis_center_vectors) ?
        unit.provenance.source_axis_center_vectors : nothing
    source_center_convention =
        hasproperty(unit.provenance, :source_center_convention) ?
        unit.provenance.source_center_convention : :unavailable
    source_center_status =
        hasproperty(unit.provenance, :source_center_status) ?
        unit.provenance.source_center_status : :unavailable
    return (
        role = unit.role,
        category = category,
        kind = unit.kind,
        column_range = unit.column_range,
        retained_count = length(unit.column_range),
        support_count = length(unit.support_indices),
        support_indices = copy(unit.support_indices),
        support_states = copy(unit.support_states),
        coefficient_matrix = unit.coefficient_matrix,
        staged_unit = unit,
        support_source_semantics = support_source_semantics,
        safe_term_capability = safe_term_capability,
        active_representation_stage = active_representation_stage,
        raw_box_auxiliary_metadata = raw_box_auxiliary_metadata,
        source_axis_center_vectors = source_axis_center_vectors,
        source_center_convention = source_center_convention,
        source_center_status = source_center_status,
        raw_product_box_operator_contract =
            category == :product_doside && active_representation_stage == :product_doside_bridge,
        route_descriptor_emitted = false,
        construction_mutated = false,
        sidecar_installation = false,
        packet_adoption = false,
        provenance = (
            source = source_helper,
            original_provenance = unit.provenance,
        ),
        diagnostics = (
            source = :pqs_current_route_retained_unit_inventory_unit,
            original_diagnostics = unit.diagnostics,
            private_diagnostic_only = true,
            active_current_route_unit = true,
            retained_weight_semantics = :not_positive_quadrature_weights,
            ida_weight_division_allowed = false,
        ),
    )
end

function _pqs_shared_shell_realized_fixture(
    construction::_BondAlignedDiatomicHighOrderRecipeSourceConstruction3D,
    fact,
    build;
    shared_shell_index::Int,
    shared_shell_count::Int,
)
    fact.primitive_family == :projected_q_shell || throw(
        ArgumentError("PQS current-route inventory shared shell must be projected_q_shell"),
    )
    fact.classification == :out_of_scope || throw(
        ArgumentError("PQS current-route inventory shared shell must remain out-of-scope for body units"),
    )
    build.region_role == fact.role || throw(
        ArgumentError("PQS current-route inventory shared shell role mismatch"),
    )
    build.primitive_family == :projected_q_shell || throw(
        ArgumentError("PQS current-route inventory shared shell build must be projected_q_shell"),
    )
    layer = build.built_object
    descriptor = hasproperty(build.metadata, :pqs_staged_unit_descriptor) ?
        build.metadata.pqs_staged_unit_descriptor :
        _nested_projected_q_shell_staged_unit_descriptor(layer)
    descriptor isa _CartesianNestedProjectedQShellStagedUnitDescriptor3D || throw(
        ArgumentError("PQS current-route inventory shared shell requires a staged descriptor"),
    )
    dims = _nested_axis_lengths(construction.axis_bundles)
    parent_dimension = prod(dims)
    size(layer.coefficient_matrix, 1) == parent_dimension || throw(
        ArgumentError("PQS current-route inventory shared shell coefficient matrix has inconsistent parent dimension"),
    )
    size(layer.coefficient_matrix, 2) == build.retained_count == fact.retained_count ||
        throw(
            ArgumentError("PQS current-route inventory shared shell retained count mismatch"),
        )
    support_indices = copy(layer.support_indices)
    support_states = copy(layer.support_states)
    support_indices == descriptor.support_indices || throw(
        ArgumentError("PQS current-route inventory shared shell support indices do not match descriptor"),
    )
    support_states == descriptor.support_states || throw(
        ArgumentError("PQS current-route inventory shared shell support states do not match descriptor"),
    )
    support_count = length(support_indices)
    support_count == build.built_support_count == fact.support_count == descriptor.support_count ||
        throw(
            ArgumentError("PQS current-route inventory shared shell support count mismatch"),
        )
    build.column_range == fact.column_range || throw(
        ArgumentError("PQS current-route inventory shared shell column range mismatch"),
    )
    local_coefficients =
        Matrix{Float64}(layer.coefficient_matrix[support_indices, :])
    size(local_coefficients) == descriptor.support_local_coefficient_shape || throw(
        ArgumentError("PQS current-route inventory shared shell local coefficient shape mismatch"),
    )
    parent_coefficients = _pqs_support_local_parent_coefficient_matrix(
        support_indices,
        local_coefficients,
        parent_dimension,
    )
    parent_difference =
        Matrix{Float64}(parent_coefficients) -
        Matrix{Float64}(layer.coefficient_matrix)
    max_parent_coefficient_error =
        isempty(parent_difference) ? 0.0 : maximum(abs, parent_difference)
    raw_box_auxiliary_reference_available =
        !isempty(descriptor.axis_intervals) &&
        !isempty(descriptor.axis_local_coefficients) &&
        descriptor.mode_count == length(descriptor.boundary_mode_indices)
    shell_realization_transform_fact =
        _pqs_current_route_shell_realization_transform_fact(
            descriptor;
            role = shared_shell_count == 1 ?
                   fact.role :
                   Symbol(string(fact.role), "_", string(shared_shell_index)),
            original_role = fact.role,
            column_range = build.column_range,
            support_indices = support_indices,
            support_states = support_states,
            support_local_coefficient_matrix = local_coefficients,
        )
    unique_role =
        shared_shell_count == 1 ?
        fact.role :
        Symbol(string(fact.role), "_", string(shared_shell_index))
    return (
        role = unique_role,
        original_role = fact.role,
        original_column_range = fact.column_range,
        shared_shell_index = shared_shell_index,
        shared_shell_count = shared_shell_count,
        category = :shell_realized_pqs_fixture,
        kind = :projected_q_shell,
        column_range = build.column_range,
        retained_count = build.retained_count,
        support_count = support_count,
        support_indices = support_indices,
        support_states = support_states,
        support_local_coefficient_matrix = local_coefficients,
        descriptor = descriptor,
        shell_realization_transform_fact = shell_realization_transform_fact,
        support_source_semantics = :shell_realized_support_local_coefficients,
        safe_term_capability = :support_local_oracle_for_shell_realization,
        active_representation_stage = :shell_realized_pqs_fixture,
        raw_product_box_operator_contract = false,
        raw_box_auxiliary_metadata = (
            available = raw_box_auxiliary_reference_available,
            source_mode_dims = (descriptor.q, descriptor.q, descriptor.L),
            raw_q = descriptor.q,
            raw_L = descriptor.L,
            mode_count = descriptor.mode_count,
            boundary_mode_count = length(descriptor.boundary_mode_indices),
            reference_only = true,
            active_current_route_contract = false,
        ),
        equivalence = (
            support_indices_match = support_indices == descriptor.support_indices,
            support_states_match = support_states == descriptor.support_states,
            retained_count_match =
                build.retained_count == fact.retained_count == descriptor.retained_count,
            column_range_match = build.column_range == fact.column_range,
            original_role_match = build.region_role == fact.role,
            coefficient_matrix_matches_active_shell =
                max_parent_coefficient_error == 0.0,
            max_parent_coefficient_error = max_parent_coefficient_error,
            local_support_coefficient_shape = size(local_coefficients),
        ),
        route_descriptor_emitted = false,
        construction_mutated = false,
        sidecar_installation = false,
        packet_adoption = false,
        provenance = (
            source = :pqs_current_route_shared_shell_realized_fixture,
            original_role = fact.role,
            original_column_range = fact.column_range,
            shared_shell_index = shared_shell_index,
            shared_shell_count = shared_shell_count,
            primitive_family = build.primitive_family,
            mapped_primitive = build.mapped_primitive,
        ),
        diagnostics = (
            source = :pqs_current_route_shared_shell_realized_fixture,
            private_diagnostic_only = true,
            active_current_route_unit = true,
            unique_role = unique_role,
            original_role = fact.role,
            original_column_range = fact.column_range,
            shared_shell_index = shared_shell_index,
            shared_shell_count = shared_shell_count,
            representation_stage = :shell_realized_pqs_fixture,
            shell_projection_lowdin_realization = true,
            raw_product_box_operator_contract = false,
            source_box_algorithm_available = false,
            support_local_oracle_used = true,
            shell_row_oracle_only = true,
            shell_realization_transform_fact_available = true,
            source_box_operator_application_ready =
                shell_realization_transform_fact.source_box_operator_application_ready,
            compact_source_space_transform_available =
                shell_realization_transform_fact.compact_source_space_transform.available,
            raw_box_auxiliary_reference_available =
                raw_box_auxiliary_reference_available,
            retained_weight_semantics = :not_positive_quadrature_weights,
            ida_weight_division_allowed = false,
        ),
    )
end

function _pqs_current_route_shell_realization_transform_fact(
    descriptor::_CartesianNestedProjectedQShellStagedUnitDescriptor3D;
    role = descriptor.role,
    original_role = descriptor.role,
    column_range = nothing,
    support_indices = descriptor.support_indices,
    support_states = descriptor.support_states,
    support_local_coefficient_matrix = nothing,
    metrics = nothing,
    coefficient_atol::Real = 1.0e-12,
    isometry_atol::Real = 1.0e-8,
)
    descriptor.kind == :projected_q_shell || throw(
        ArgumentError("PQS shell-realization transform fact requires a projected q-shell descriptor"),
    )
    length(support_indices) == descriptor.support_count || throw(
        ArgumentError("PQS shell-realization transform fact support index count mismatch"),
    )
    length(support_states) == descriptor.support_count || throw(
        ArgumentError("PQS shell-realization transform fact support state count mismatch"),
    )
    source_mode_dims =
        ntuple(axis -> size(descriptor.axis_local_coefficients[axis], 2), 3)
    source_mode_count = prod(source_mode_dims)
    boundary_mode_count = length(descriptor.boundary_mode_indices)
    boundary_mode_count == descriptor.mode_count || throw(
        DimensionMismatch("PQS shell-realization boundary mode count must match descriptor mode_count"),
    )
    length(descriptor.boundary_column_indices) == descriptor.mode_count || throw(
        DimensionMismatch("PQS shell-realization boundary column count must match descriptor mode_count"),
    )
    cleanup_transform = Matrix{Float64}(descriptor.cleanup_transform)
    size(cleanup_transform) == (descriptor.mode_count, descriptor.retained_count) ||
        throw(
            DimensionMismatch("PQS shell-realization cleanup transform shape mismatch"),
        )
    shell_projection_matrix =
        _nested_projected_q_shell_descriptor_seed_coefficients(descriptor)
    size(shell_projection_matrix) == (descriptor.support_count, descriptor.mode_count) ||
        throw(
            DimensionMismatch("PQS shell-realization projection shape mismatch"),
        )

    stored_shape = isnothing(support_local_coefficient_matrix) ?
                   descriptor.support_local_coefficient_shape :
                   size(support_local_coefficient_matrix)
    stored_shape == descriptor.support_local_coefficient_shape || throw(
        DimensionMismatch("PQS shell-realization support-local coefficient shape mismatch"),
    )
    realized_support_coefficients = shell_projection_matrix * cleanup_transform
    size(realized_support_coefficients) == descriptor.support_local_coefficient_shape ||
        throw(
            DimensionMismatch("PQS shell-realization realized support coefficient shape mismatch"),
        )
    max_support_local_coefficient_error = nothing
    coefficient_matches_descriptor_realization = nothing
    if !isnothing(support_local_coefficient_matrix)
        stored_coefficients = Matrix{Float64}(support_local_coefficient_matrix)
        difference = realized_support_coefficients - stored_coefficients
        max_support_local_coefficient_error =
            isempty(difference) ? 0.0 : maximum(abs, difference)
        coefficient_matches_descriptor_realization =
            max_support_local_coefficient_error <= coefficient_atol
    end

    isometry_checked = false
    isometry_error = nothing
    isometric = nothing
    if !isnothing(metrics)
        shell_plan = _pqs_shell_realization_plan(
            descriptor,
            metrics;
            isometry_atol = isometry_atol,
        )
        isometry_checked = true
        isometry_error = shell_plan.isometry_error
        isometric = shell_plan.isometric
    end

    compact_missing_reason =
        :shell_projection_maps_selected_modes_to_shell_rows_not_compact_raw_mode_space
    return (
        object_kind = :pqs_current_route_shell_realization_transform_fact,
        status = :metadata_precursor,
        role = role,
        original_role = original_role,
        column_range = column_range,
        representation_stage = :shell_realized_pqs_fixture,
        transform_contract = :raw_box_boundary_selection_shell_projection_lowdin,
        source_box = (
            source_mode_dims = source_mode_dims,
            source_mode_count = source_mode_count,
            axis_intervals = descriptor.axis_intervals,
            axis_local_coefficient_shapes =
                ntuple(axis -> size(descriptor.axis_local_coefficients[axis]), 3),
            source_mode_dims_are_total_lengths = true,
        ),
        boundary_selection = (
            stage = :boundary_source_mode_selection,
            mode_indices = copy(descriptor.boundary_mode_indices),
            column_indices = copy(descriptor.boundary_column_indices),
            mode_count = descriptor.mode_count,
            retained_mode_count = boundary_mode_count,
            selector_matrix_shape = (source_mode_count, descriptor.mode_count),
            selector_matrix_materialized = false,
        ),
        shell_projection = (
            stage = :shell_projection_to_shell_rows,
            support_indices = copy(support_indices),
            support_states = copy(support_states),
            support_count = descriptor.support_count,
            matrix_shape = size(shell_projection_matrix),
            matrix_available_from_descriptor = true,
            matrix_materialized_in_fact = false,
        ),
        lowdin_cleanup = (
            stage = :full_rank_symmetric_lowdin_cleanup,
            method = descriptor.cleanup_method,
            transform_shape = size(cleanup_transform),
            rank_count = descriptor.cleanup_rank_count,
            rank_drop_count = descriptor.cleanup_rank_drop_count,
            cutoff = descriptor.cleanup_cutoff,
            eigenvalue_count = length(descriptor.cleanup_eigenvalues),
        ),
        retained_columns = (
            retained_count = descriptor.retained_count,
            support_local_coefficient_shape = descriptor.support_local_coefficient_shape,
            coefficient_matches_descriptor_realization =
                coefficient_matches_descriptor_realization,
            max_support_local_coefficient_error =
                max_support_local_coefficient_error,
        ),
        shell_realization = (
            shell_projection_used = true,
            lowdin_cleanup_used = true,
            shell_projection_lowdin_realization = true,
            isometry_checked = isometry_checked,
            isometry_error = isometry_error,
            isometric = isometric,
        ),
        compact_source_space_transform = (
            available = false,
            materialized = false,
            shape = nothing,
            missing_reason = compact_missing_reason,
            boundary_selector_cleanup_shape =
                (source_mode_count, descriptor.retained_count),
            boundary_selector_cleanup_is_not_shell_realization = true,
        ),
        source_box_operator_application_ready = false,
        missing_for_source_box_operator_application = (
            :exact_compact_source_space_shell_realization_transform,
        ),
        diagnostics = (
            source = :pqs_current_route_shell_realization_transform_fact,
            private_metadata_only = true,
            metadata_precursor = true,
            representation_stage = :shell_realized_pqs_fixture,
            support_local_oracle_used = true,
            shell_row_oracle_only = true,
            shell_projection_used = true,
            lowdin_cleanup_used = true,
            compact_source_space_transform_available = false,
            compact_source_space_transform_missing_reason =
                compact_missing_reason,
            source_box_operator_application_ready = false,
            raw_product_box_operator_contract = false,
            packet_adoption = false,
            fixed_block_construction_changed = false,
            qwhamiltonian_changed = false,
            retained_weight_semantics = :not_positive_quadrature_weights,
            ida_weight_division_allowed = false,
        ),
    )
end

function _pqs_current_route_shell_realization_transform_fact(
    unit;
    metrics = nothing,
    coefficient_atol::Real = 1.0e-12,
    isometry_atol::Real = 1.0e-8,
)
    hasproperty(unit, :category) && unit.category == :shell_realized_pqs_fixture ||
        throw(
            ArgumentError("PQS shell-realization transform fact requires a shell-realized PQS fixture unit"),
        )
    return _pqs_current_route_shell_realization_transform_fact(
        unit.descriptor;
        role = unit.role,
        original_role = hasproperty(unit, :original_role) ? unit.original_role : unit.role,
        column_range = unit.column_range,
        support_indices = unit.support_indices,
        support_states = unit.support_states,
        support_local_coefficient_matrix = unit.support_local_coefficient_matrix,
        metrics = metrics,
        coefficient_atol = coefficient_atol,
        isometry_atol = isometry_atol,
    )
end

function _pqs_shared_shell_realized_fixtures(
    construction::_BondAlignedDiatomicHighOrderRecipeSourceConstruction3D,
    audit,
)
    shared_facts = sort(
        [
            fact for fact in audit.unit_facts
            if fact.role == :regular_shared_molecular_shell
        ];
        by = fact -> first(fact.column_range),
    )
    !isempty(shared_facts) || throw(
        ArgumentError("PQS current-route inventory requires at least one shared PQS shell fact"),
    )
    builds = sort(
        [
            build for build in construction.region_builds
            if build.region_role == :regular_shared_molecular_shell
        ];
        by = build -> first(build.column_range),
    )
    length(builds) == length(shared_facts) || throw(
        ArgumentError("PQS current-route inventory shared PQS shell fact/build count mismatch"),
    )
    shared_shell_count = length(shared_facts)
    return Tuple(
        _pqs_shared_shell_realized_fixture(
            construction,
            fact,
            build;
            shared_shell_index = index,
            shared_shell_count = shared_shell_count,
        ) for (index, (fact, build)) in enumerate(zip(shared_facts, builds))
    )
end

function _pqs_current_route_inventory_pair_policies()
    return (
        (
            pair_type = :product_product,
            policy = :product_doside_source_box_path,
            active_current_route = true,
            active_algorithmic_policy = true,
            source_box_algorithm_available = true,
            support_local_oracle_used = false,
            shell_row_oracle_only = false,
        ),
        (
            pair_type = :support_support,
            policy = :support_local_fallback,
            active_current_route = true,
            active_algorithmic_policy = true,
            source_box_algorithm_available = false,
            support_local_oracle_used = false,
            shell_row_oracle_only = false,
        ),
        (
            pair_type = :support_product,
            policy = :support_local_fallback_unless_both_product_doside,
            active_current_route = true,
            active_algorithmic_policy = true,
            source_box_algorithm_available = false,
            support_local_oracle_used = false,
            shell_row_oracle_only = false,
        ),
        (
            pair_type = :shell_realized_pqs_product,
            policy = :support_local_oracle_for_shell_realization,
            active_current_route = false,
            active_algorithmic_policy = false,
            source_box_algorithm_available = false,
            support_local_oracle_used = true,
            shell_row_oracle_only = true,
            raw_box_pqs_helper_status = :reference_shadow_only,
        ),
        (
            pair_type = :shell_realized_pqs_support,
            policy = :support_local_oracle_for_shell_realization,
            active_current_route = false,
            active_algorithmic_policy = false,
            source_box_algorithm_available = false,
            support_local_oracle_used = true,
            shell_row_oracle_only = true,
        ),
        (
            pair_type = :shell_realized_pqs_pqs,
            policy = :support_local_oracle_for_shell_realization,
            active_current_route = false,
            active_algorithmic_policy = false,
            source_box_algorithm_available = false,
            support_local_oracle_used = true,
            shell_row_oracle_only = true,
        ),
        (
            pair_type = :raw_box_pqs_helpers,
            policy = :reference_shadow_only_not_active_current_route,
            active_current_route = false,
            active_algorithmic_policy = false,
            source_box_algorithm_available = true,
            support_local_oracle_used = false,
            shell_row_oracle_only = false,
        ),
    )
end

function _pqs_current_route_inventory_coverage(units)
    ordered_units = sort(collect(units); by = unit -> first(unit.column_range))
    non_overlapping = true
    contiguous = true
    expected_first = first(first(ordered_units).column_range)
    represented_count = 0
    for unit in ordered_units
        first(unit.column_range) == expected_first || (contiguous = false)
        represented_count += length(unit.column_range)
        expected_first = last(unit.column_range) + 1
    end
    for index in 1:(length(ordered_units) - 1)
        last(ordered_units[index].column_range) < first(ordered_units[index + 1].column_range) ||
            (non_overlapping = false)
    end
    return (
        first_column = first(first(ordered_units).column_range),
        last_column = last(last(ordered_units).column_range),
        represented_count = represented_count,
        unit_count = length(ordered_units),
        contiguous = contiguous,
        non_overlapping = non_overlapping,
        covers_every_column_once = contiguous && non_overlapping,
        ordered_roles = Tuple(unit.role for unit in ordered_units),
        ordered_column_ranges = Tuple(unit.column_range for unit in ordered_units),
    )
end

function _pqs_current_route_retained_unit_inventory(
    construction::_BondAlignedDiatomicHighOrderRecipeSourceConstruction3D;
    audit = _pqs_route_retained_unit_fact_audit(
        construction;
        include_support_indices = true,
    ),
)
    outer_fixture =
        _pqs_outer_mismatch_product_doside_units(construction; audit = audit)
    atom_fixture =
        _pqs_atom_box_support_dense_units(construction; audit = audit)
    contact_fixture =
        _pqs_contact_cap_product_doside_unit(construction; audit = audit)
    shared_fixtures =
        _pqs_shared_shell_realized_fixtures(construction, audit)

    units = Any[]
    append!(
        units,
        _pqs_inventory_wrapped_staged_unit(
            unit;
            category = :product_doside,
            support_source_semantics = :identity_selector_boundary_slab,
            safe_term_capability = :product_doside_source_box_safe_terms,
            active_representation_stage = :product_doside_bridge,
            source_helper = :_pqs_outer_mismatch_product_doside_units,
        ) for unit in outer_fixture.units
    )
    append!(
        units,
        _pqs_inventory_wrapped_staged_unit(
            unit;
            category = :support_dense,
            support_source_semantics = :support_local_direct_rows,
            safe_term_capability = :support_local_fallback_safe_terms,
            active_representation_stage = :support_dense_direct_support,
            source_helper = :_pqs_atom_box_support_dense_units,
        ) for unit in atom_fixture.units
    )
    push!(
        units,
        _pqs_inventory_wrapped_staged_unit(
            contact_fixture.unit;
            category = :product_doside,
            support_source_semantics = :identity_selector_contact_cap_slab,
            safe_term_capability = :product_doside_source_box_safe_terms,
            active_representation_stage = :product_doside_bridge,
            source_helper = :_pqs_contact_cap_product_doside_unit,
        ),
    )
    append!(units, shared_fixtures)
    ordered_units = Tuple(sort(units; by = unit -> first(unit.column_range)))
    by_role = Dict(unit.role => unit for unit in ordered_units)
    length(by_role) == length(ordered_units) || throw(
        ArgumentError("PQS current-route retained-unit inventory requires unique unit roles"),
    )
    coverage = _pqs_current_route_inventory_coverage(ordered_units)
    q4_single_shared_expected_roles = (
        :outer_mismatch_z_low_slab,
        :outer_mismatch_z_high_slab,
        :left_atom_box,
        :right_atom_box,
        :contact_cap_slab,
        :regular_shared_molecular_shell,
    )
    fixed_dimension = sum(build.retained_count for build in construction.region_builds)
    coverage.first_column == 1 || throw(
        ArgumentError("PQS current-route retained-unit inventory must start at column 1"),
    )
    coverage.last_column == fixed_dimension || throw(
        ArgumentError("PQS current-route retained-unit inventory must end at fixed dimension"),
    )
    coverage.represented_count == fixed_dimension || throw(
        ArgumentError("PQS current-route retained-unit inventory represented count mismatch"),
    )
    coverage.covers_every_column_once || throw(
        ArgumentError("PQS current-route retained-unit inventory must cover every column once"),
    )
    raw_box_available = all(
        fixture -> fixture.raw_box_auxiliary_metadata.available,
        shared_fixtures,
    )
    diagnostics = (
        source = :pqs_current_route_retained_unit_inventory,
        private_diagnostic_only = true,
        current_route_inventory = true,
        route_descriptor_emitted = false,
        construction_mutated = false,
        sidecar_installation = false,
        packet_adoption = false,
        fixed_block_construction_changed = false,
        qwhamiltonian_changed = false,
        ida_weight_division_allowed = false,
        retained_weight_semantics = :not_positive_quadrature_weights,
        local_ecp_gaussian_mwg_interaction_changed = false,
        shared_pqs_unit_count = length(shared_fixtures),
        shared_pqs_roles = Tuple(fixture.role for fixture in shared_fixtures),
        shared_pqs_original_roles =
            Tuple(fixture.original_role for fixture in shared_fixtures),
        shared_pqs_column_ranges =
            Tuple(fixture.column_range for fixture in shared_fixtures),
        shared_pqs_original_column_ranges =
            Tuple(fixture.original_column_range for fixture in shared_fixtures),
        shared_pqs_active_representation = :shell_realized_pqs_fixture,
        shared_pqs_raw_box_operator_contract = false,
        shared_pqs_shell_realization_transform_fact_count =
            count(
                fixture -> hasproperty(fixture, :shell_realization_transform_fact),
                shared_fixtures,
            ),
        shared_pqs_source_box_operator_application_ready_count =
            count(
                fixture ->
                    fixture.shell_realization_transform_fact.source_box_operator_application_ready,
                shared_fixtures,
            ),
        raw_box_pqs_auxiliary_reference_available = raw_box_available,
        raw_box_pqs_auxiliary_reference_unavailable_reason =
            raw_box_available ? nothing : :shared_pqs_descriptor_missing_raw_plan_facts,
        whole_route_safe_term_matrix_consumer = false,
        fixed_dimension = fixed_dimension,
        unit_count = length(ordered_units),
        coverage_complete = coverage.covers_every_column_once,
        q4_single_shared_role_order_preserved =
            coverage.ordered_roles == q4_single_shared_expected_roles,
    )
    return (
        object_kind = :pqs_current_route_retained_unit_inventory_fixture,
        status = :private_diagnostic_only,
        units = ordered_units,
        by_role = by_role,
        coverage = coverage,
        pair_policies = _pqs_current_route_inventory_pair_policies(),
        source_fixtures = (
            outer_mismatch = outer_fixture,
            atom_boxes = atom_fixture,
            contact_cap = contact_fixture,
            shared_pqs = shared_fixtures,
        ),
        diagnostics = diagnostics,
    )
end

function _pqs_fixed_side_unit_get(object, field::Symbol, default = nothing)
    return hasproperty(object, field) ? getproperty(object, field) : default
end

function _pqs_fixed_side_unit_provenance_get(unit, field::Symbol, default = nothing)
    provenance = _pqs_fixed_side_unit_get(unit, :provenance, (;))
    return hasproperty(provenance, field) ? getproperty(provenance, field) : default
end

function _pqs_fixed_side_unit_diagnostic_get(unit, field::Symbol, default = nothing)
    diagnostics = _pqs_fixed_side_unit_get(unit, :diagnostics, (;))
    return hasproperty(diagnostics, field) ? getproperty(diagnostics, field) : default
end

function _pqs_fixed_side_unit_source_axis_center_vectors(unit, local_dims::NTuple{3,Int})
    source_axis_center_vectors =
        _pqs_fixed_side_unit_get(unit, :source_axis_center_vectors, nothing)
    isnothing(source_axis_center_vectors) && return nothing
    length(source_axis_center_vectors) == 3 || throw(
        ArgumentError("product/doside source center metadata requires three axis center vectors"),
    )
    return ntuple(axis -> begin
        centers = Float64.(source_axis_center_vectors[axis])
        length(centers) == local_dims[axis] || throw(
            DimensionMismatch("product/doside source center vector length does not match local axis dimension"),
        )
        centers
    end, 3)
end

function _pqs_fixed_side_unit_source_center_convention(unit)
    return _pqs_fixed_side_unit_get(
        unit,
        :source_center_convention,
        :unavailable,
    )
end

function _pqs_fixed_side_unit_source_center_status(unit)
    return _pqs_fixed_side_unit_get(
        unit,
        :source_center_status,
        :unavailable,
    )
end

function _pqs_source_mode_tuple3(value, context::AbstractString)
    length(value) == 3 || throw(
        ArgumentError("$(context) must be a three-axis source-mode tuple"),
    )
    return (Int(value[1]), Int(value[2]), Int(value[3]))
end

function _pqs_source_mode_tuple_label(shell_id::Int, mode_tuple::NTuple{3,Int})
    return string(
        "source_mode:",
        shell_id,
        ":",
        mode_tuple[1],
        ",",
        mode_tuple[2],
        ",",
        mode_tuple[3],
    )
end

function _pqs_current_route_unit_support_states(unit)
    raw_states = _pqs_fixed_side_unit_get(unit, :support_states, nothing)
    isnothing(raw_states) && return nothing
    states = NTuple{3,Int}[
        _pqs_source_mode_tuple3(state, "support-dense source support state")
        for state in raw_states
    ]
    isempty(states) && return nothing
    return states
end

function _pqs_source_mode_axis_bounds(states::AbstractVector{<:NTuple{3,Int}})
    isempty(states) && throw(
        ArgumentError("source-mode axis bounds require at least one source state"),
    )
    starts = ntuple(axis -> minimum(state[axis] for state in states), 3)
    stops = ntuple(axis -> maximum(state[axis] for state in states), 3)
    dims = ntuple(axis -> stops[axis] - starts[axis] + 1, 3)
    return (starts = starts, stops = stops, dims = dims)
end

function _pqs_current_route_shell_realized_source_mode_metadata(unit)
    if _pqs_fixed_side_unit_get(unit, :descriptor, nothing) isa
       _CartesianNestedProjectedQShellStagedUnitDescriptor3D
        descriptor = unit.descriptor
        source_mode_dims =
            ntuple(axis -> size(descriptor.axis_local_coefficients[axis], 2), 3)
        mode_indices = NTuple{3,Int}[
            _pqs_source_mode_tuple3(
                mode,
                "shell-realized PQS boundary source mode",
            )
            for mode in descriptor.boundary_mode_indices
        ]
        axis_intervals = descriptor.axis_intervals
        selection_rule = descriptor.selection_rule
        source = :pqs_projected_q_shell_descriptor
    else
        fact =
            _pqs_fixed_side_unit_get(unit, :shell_realization_transform_fact, nothing)
        raw_metadata =
            _pqs_fixed_side_unit_get(unit, :raw_box_auxiliary_metadata, nothing)
        if !isnothing(fact) &&
           hasproperty(fact, :source_box) &&
           hasproperty(fact, :boundary_selection)
            source_mode_dims =
                _pqs_source_mode_tuple3(
                    fact.source_box.source_mode_dims,
                    "shell-realized PQS source-mode dimensions",
                )
            mode_indices = NTuple{3,Int}[
                _pqs_source_mode_tuple3(
                    mode,
                    "shell-realized PQS boundary source mode",
                )
                for mode in fact.boundary_selection.mode_indices
            ]
            axis_intervals = hasproperty(fact.source_box, :axis_intervals) ?
                fact.source_box.axis_intervals :
                ntuple(axis -> 1:source_mode_dims[axis], 3)
            selection_rule =
                hasproperty(fact.boundary_selection, :selection_rule) ?
                fact.boundary_selection.selection_rule :
                :boundary_source_mode_selection
            source = :pqs_shell_realization_transform_fact
        elseif !isnothing(raw_metadata) &&
               hasproperty(raw_metadata, :available) &&
               Bool(raw_metadata.available) &&
               hasproperty(raw_metadata, :source_mode_dims) &&
               hasproperty(raw_metadata, :boundary_mode_indices)
            source_mode_dims =
                _pqs_source_mode_tuple3(
                    raw_metadata.source_mode_dims,
                    "shell-realized PQS source-mode dimensions",
                )
            mode_indices = NTuple{3,Int}[
                _pqs_source_mode_tuple3(
                    mode,
                    "shell-realized PQS boundary source mode",
                )
                for mode in raw_metadata.boundary_mode_indices
            ]
            axis_intervals =
                hasproperty(raw_metadata, :axis_intervals) ?
                raw_metadata.axis_intervals :
                ntuple(axis -> 1:source_mode_dims[axis], 3)
            selection_rule =
                hasproperty(raw_metadata, :selection_rule) ?
                raw_metadata.selection_rule :
                :boundary_source_mode_selection
            source = :pqs_raw_box_auxiliary_metadata
        else
            return nothing
        end
    end

    isempty(mode_indices) && return nothing
    length(axis_intervals) == 3 || throw(
        ArgumentError("shell-realized PQS source-box axis intervals must be three-dimensional"),
    )
    axis_intervals = (axis_intervals[1], axis_intervals[2], axis_intervals[3])
    for axis in 1:3
        isempty(axis_intervals[axis]) && throw(
            ArgumentError("shell-realized PQS source-box axis intervals must be nonempty"),
        )
    end
    for mode in mode_indices
        all(axis -> 1 <= mode[axis] <= source_mode_dims[axis], 1:3) ||
            throw(
                ArgumentError("shell-realized PQS boundary source mode exceeds source-mode dimensions"),
            )
    end
    return (
        source = source,
        source_mode_dims = source_mode_dims,
        axis_intervals = axis_intervals,
        mode_indices = mode_indices,
        selection_rule = selection_rule,
    )
end

function _pqs_fixed_side_unit_raw_box_auxiliary_metadata(unit)
    raw_box = _pqs_fixed_side_unit_get(unit, :raw_box_auxiliary_metadata, nothing)
    return isnothing(raw_box) ? nothing : raw_box
end

function _pqs_fixed_side_unit_source_mode_dims(unit)
    raw_box = _pqs_fixed_side_unit_raw_box_auxiliary_metadata(unit)
    !isnothing(raw_box) && hasproperty(raw_box, :source_mode_dims) &&
        return raw_box.source_mode_dims
    transform_fact =
        _pqs_fixed_side_unit_get(unit, :shell_realization_transform_fact, nothing)
    if !isnothing(transform_fact) &&
       hasproperty(transform_fact, :source_box) &&
       hasproperty(transform_fact.source_box, :source_mode_dims)
        return transform_fact.source_box.source_mode_dims
    end
    return nothing
end

function _pqs_fixed_side_unit_source_dimension(unit)
    source_mode_dims = _pqs_fixed_side_unit_source_mode_dims(unit)
    isnothing(source_mode_dims) || return prod(source_mode_dims)
    return _pqs_fixed_side_unit_get(unit, :source_dimension, nothing)
end

function _pqs_fixed_side_unit_class(category::Symbol)
    category == :product_doside && return :product_doside
    category == :support_dense && return :support_dense_direct_support
    category == :shell_realized_pqs_fixture &&
        return :shell_realized_pqs_fixture
    return :unknown
end

function _pqs_current_route_product_axis_source_index(axis, local_index::Integer)
    if axis.kind == :fixed
        local_index == 1 || throw(
            ArgumentError("fixed product-axis relation requires local index 1"),
        )
        !isnothing(axis.fixed_index) || throw(
            ArgumentError("fixed product-axis relation requires a fixed index"),
        )
        return Int(axis.fixed_index)
    elseif axis.kind == :active
        interval = axis.interval
        !isnothing(interval) || throw(
            ArgumentError("active product-axis relation requires an interval"),
        )
        1 <= local_index <= length(interval) || throw(
            ArgumentError("active product-axis local index $(local_index) outside $(interval)"),
        )
        return first(interval) + Int(local_index) - 1
    end
    throw(ArgumentError("unsupported product-axis kind $(axis.kind)"))
end

function _pqs_current_route_product_axis_local_dims(staged_unit)
    hasproperty(staged_unit, :axes) && hasproperty(staged_unit, :axis_function_indices) ||
        throw(
            ArgumentError("product/doside source-shell modes require axes and axis_function_indices"),
        )
    length(staged_unit.axes) == 3 || throw(
        ArgumentError("product/doside source-shell modes require three staged axes"),
    )
    isempty(staged_unit.axis_function_indices) && throw(
        ArgumentError("product/doside source-shell modes require at least one axis tuple"),
    )
    return ntuple(
        axis -> maximum(tuple -> Int(tuple[axis]), staged_unit.axis_function_indices),
        3,
    )
end

function _pqs_current_route_source_shell_mode_inventory(
    inventory;
    provenance = (;),
    strict::Bool = true,
)
    inventory.object_kind == :pqs_current_route_retained_unit_inventory_fixture ||
        throw(
            ArgumentError("source-shell/source-mode inventory requires _pqs_current_route_retained_unit_inventory output"),
        )
    units = Tuple(inventory.units)
    isempty(units) && throw(
        ArgumentError("source-shell/source-mode inventory requires at least one retained unit"),
    )
    coverage = inventory.coverage
    coverage_complete = hasproperty(coverage, :covers_every_column_once) ?
        Bool(coverage.covers_every_column_once) : false
    strict && !coverage_complete && throw(
        ArgumentError("source-shell/source-mode inventory requires complete fixed-side column coverage"),
    )

    product_unit_entries = Tuple(
        (unit_index = unit_index, unit = unit)
        for (unit_index, unit) in enumerate(units)
        if unit.category == :product_doside && unit.kind == :product_doside
    )
    support_dense_unit_entries = NamedTuple[]
    shell_realized_unit_entries = NamedTuple[]
    product_shell_count = length(product_unit_entries)
    product_mode_count = 0
    support_dense_source_shell_count = 0
    support_dense_source_mode_count = 0
    shell_realized_source_shell_count = 0
    shell_realized_source_mode_count = 0
    support_dense_unavailable_unit_count = 0
    support_dense_unavailable_column_count = 0
    shell_realized_unavailable_unit_count = 0
    shell_realized_unavailable_column_count = 0
    other_unavailable_unit_count = 0
    other_unavailable_column_count = 0
    for (unit_index, unit) in enumerate(units)
        unit_count = length(unit.column_range)
        if unit.category == :product_doside && unit.kind == :product_doside
            hasproperty(unit, :staged_unit) && !isnothing(unit.staged_unit) ||
                throw(
                    ArgumentError("product/doside source-shell modes require staged-axis metadata"),
                )
            product_mode_count += length(unit.staged_unit.axis_function_indices)
        elseif unit.category == :support_dense
            support_states = _pqs_current_route_unit_support_states(unit)
            if isnothing(support_states)
                support_dense_unavailable_unit_count += 1
                support_dense_unavailable_column_count += unit_count
            else
                length(support_states) == Int(unit.support_count) || throw(
                    DimensionMismatch("support-dense source-mode labels require support_states to match support_count"),
                )
                push!(
                    support_dense_unit_entries,
                    (
                        unit_index = unit_index,
                        unit = unit,
                        support_states = support_states,
                    ),
                )
                support_dense_source_shell_count += 1
                support_dense_source_mode_count += length(support_states)
            end
        elseif unit.category == :shell_realized_pqs_fixture
            shell_metadata =
                _pqs_current_route_shell_realized_source_mode_metadata(unit)
            if isnothing(shell_metadata)
                shell_realized_unavailable_unit_count += 1
                shell_realized_unavailable_column_count += unit_count
            else
                push!(
                    shell_realized_unit_entries,
                    (
                        unit_index = unit_index,
                        unit = unit,
                        metadata = shell_metadata,
                    ),
                )
                shell_realized_source_shell_count += 1
                shell_realized_source_mode_count +=
                    length(shell_metadata.mode_indices)
            end
        else
            other_unavailable_unit_count += 1
            other_unavailable_column_count += unit_count
        end
    end

    source_shell_count =
        product_shell_count +
        support_dense_source_shell_count +
        shell_realized_source_shell_count
    source_mode_count =
        product_mode_count +
        support_dense_source_mode_count +
        shell_realized_source_mode_count

    source_shell_ids = collect(1:source_shell_count)
    source_shell_unit_indices = Vector{Int}(undef, source_shell_count)
    source_shell_unit_labels = Vector{Symbol}(undef, source_shell_count)
    source_shell_unit_categories = Vector{Symbol}(undef, source_shell_count)
    source_shell_unit_kinds = Vector{Symbol}(undef, source_shell_count)
    source_shell_retained_starts = Vector{Int}(undef, source_shell_count)
    source_shell_retained_stops = Vector{Int}(undef, source_shell_count)
    source_shell_labels = Vector{Symbol}(undef, source_shell_count)
    source_shell_statuses = Vector{Symbol}(undef, source_shell_count)
    source_shell_construction_kinds = Vector{Symbol}(undef, source_shell_count)
    source_shell_axis_kinds = Matrix{Symbol}(undef, source_shell_count, 3)
    source_shell_axis_starts = zeros(Int, source_shell_count, 3)
    source_shell_axis_stops = zeros(Int, source_shell_count, 3)
    source_shell_fixed_axis_indices = zeros(Int, source_shell_count, 3)
    source_shell_contracted_dims = zeros(Int, source_shell_count, 3)
    source_shell_mode_counts = zeros(Int, source_shell_count)
    source_shell_mode_orderings = Vector{Symbol}(undef, source_shell_count)
    source_shell_center_definitions = fill(:unavailable, source_shell_count)
    source_shell_center_statuses = fill(:unavailable, source_shell_count)
    source_shell_lowdin_correction_applied = falses(source_shell_count)
    source_shell_shell_label_statuses = fill(:unavailable, source_shell_count)
    source_shell_ray_label_statuses = fill(:unavailable, source_shell_count)
    source_shell_radial_order_statuses =
        fill(:unavailable, source_shell_count)

    mode_source_shell_ids = Vector{Int}(undef, source_mode_count)
    mode_indices = Vector{Int}(undef, source_mode_count)
    mode_unit_labels = Vector{Symbol}(undef, source_mode_count)
    native_source_id_labels = Vector{String}(undef, source_mode_count)
    local_axis_function_indices = zeros(Int, source_mode_count, 3)
    source_axis_indices = zeros(Int, source_mode_count, 3)
    parent_lattice_axis_indices = zeros(Int, source_mode_count, 3)
    source_mode_statuses = Vector{Symbol}(undef, source_mode_count)
    source_axis_tuple_statuses = Vector{Symbol}(undef, source_mode_count)
    parent_lattice_axis_statuses = fill(:unavailable, source_mode_count)
    mode_center_coordinates = fill(NaN, source_mode_count, 3)
    mode_center_definitions = fill(:unavailable, source_mode_count)
    mode_center_statuses = fill(:unavailable, source_mode_count)
    mode_lowdin_correction_applied = falses(source_mode_count)
    mode_shell_label_statuses = fill(:unavailable, source_mode_count)
    mode_ray_label_statuses = fill(:unavailable, source_mode_count)
    mode_radial_order_statuses = fill(:unavailable, source_mode_count)
    inferred_from_centers = falses(source_mode_count)
    inferred_from_nearest_grid = falses(source_mode_count)
    inferred_from_support_order = falses(source_mode_count)
    inferred_from_support_indices = falses(source_mode_count)
    inferred_from_raw_to_final_support = falses(source_mode_count)

    native_center_shell_count = 0
    native_center_mode_count = 0
    mode_row = 1
    for (shell_id, entry) in enumerate(product_unit_entries)
        unit = entry.unit
        staged_unit = unit.staged_unit
        local_dims = _pqs_current_route_product_axis_local_dims(staged_unit)
        length(staged_unit.axis_function_indices) == length(unit.column_range) ||
            throw(
                ArgumentError("product/doside source mode count must match unit retained range"),
            )
        source_axis_center_vectors =
            _pqs_fixed_side_unit_source_axis_center_vectors(unit, local_dims)
        source_center_convention =
            _pqs_fixed_side_unit_source_center_convention(unit)
        source_center_status =
            _pqs_fixed_side_unit_source_center_status(unit)
        source_centers_available = !isnothing(source_axis_center_vectors)
        if source_centers_available
            source_center_convention == :comx_construction || throw(
                ArgumentError("product/doside source center metadata must use :comx_construction"),
            )
            source_center_status == :native_representative || throw(
                ArgumentError("product/doside source center metadata must use :native_representative status"),
            )
            native_center_shell_count += 1
            native_center_mode_count += length(staged_unit.axis_function_indices)
            source_shell_center_definitions[shell_id] = source_center_convention
            source_shell_center_statuses[shell_id] = source_center_status
        end
        source_shell_unit_indices[shell_id] = entry.unit_index
        source_shell_unit_labels[shell_id] = unit.role
        source_shell_unit_categories[shell_id] = unit.category
        source_shell_unit_kinds[shell_id] = unit.kind
        source_shell_retained_starts[shell_id] = first(unit.column_range)
        source_shell_retained_stops[shell_id] = last(unit.column_range)
        source_shell_labels[shell_id] = unit.role
        source_shell_statuses[shell_id] = :native_product_doside_source_box
        source_shell_construction_kinds[shell_id] = :product_doside
        source_shell_mode_counts[shell_id] =
            length(staged_unit.axis_function_indices)
        source_shell_mode_orderings[shell_id] = :axis_function_indices_order
        for axis in 1:3
            staged_axis = staged_unit.axes[axis]
            interval = _staged_axis_interval(staged_axis)
            source_shell_axis_kinds[shell_id, axis] = staged_axis.kind
            source_shell_axis_starts[shell_id, axis] = first(interval)
            source_shell_axis_stops[shell_id, axis] = last(interval)
            source_shell_fixed_axis_indices[shell_id, axis] =
                staged_axis.kind == :fixed ? Int(staged_axis.fixed_index) : 0
            source_shell_contracted_dims[shell_id, axis] = local_dims[axis]
        end
        for (local_col, local_tuple) in enumerate(staged_unit.axis_function_indices)
            source_tuple = ntuple(
                axis -> _pqs_current_route_product_axis_source_index(
                    staged_unit.axes[axis],
                    local_tuple[axis],
                ),
                3,
            )
            mode_source_shell_ids[mode_row] = shell_id
            mode_indices[mode_row] = local_col
            mode_unit_labels[mode_row] = unit.role
            native_source_id_labels[mode_row] =
                _pqs_source_mode_tuple_label(shell_id, local_tuple)
            local_axis_function_indices[mode_row, :] .= collect(local_tuple)
            source_axis_indices[mode_row, :] .= collect(source_tuple)
            parent_lattice_axis_indices[mode_row, :] .= collect(source_tuple)
            source_mode_statuses[mode_row] =
                :native_product_doside_source_mode
            source_axis_tuple_statuses[mode_row] =
                :native_product_axis_tuple
            parent_lattice_axis_statuses[mode_row] =
                :native_product_parent_lattice_axis_tuple
            if source_centers_available
                for axis in 1:3
                    mode_center_coordinates[mode_row, axis] =
                        source_axis_center_vectors[axis][local_tuple[axis]]
                end
                mode_center_definitions[mode_row] = source_center_convention
                mode_center_statuses[mode_row] = source_center_status
            end
            mode_row += 1
        end
    end

    shell_row = product_shell_count + 1
    for entry in support_dense_unit_entries
        unit = entry.unit
        shell_id = shell_row
        support_states = entry.support_states
        bounds = _pqs_source_mode_axis_bounds(support_states)

        source_shell_unit_indices[shell_id] = entry.unit_index
        source_shell_unit_labels[shell_id] = unit.role
        source_shell_unit_categories[shell_id] = unit.category
        source_shell_unit_kinds[shell_id] = unit.kind
        source_shell_retained_starts[shell_id] = first(unit.column_range)
        source_shell_retained_stops[shell_id] = last(unit.column_range)
        source_shell_labels[shell_id] = unit.role
        source_shell_statuses[shell_id] =
            :native_support_dense_source_support_states
        source_shell_construction_kinds[shell_id] =
            :support_dense_direct_support
        source_shell_mode_counts[shell_id] = length(support_states)
        source_shell_mode_orderings[shell_id] =
            :construction_support_state_order
        for axis in 1:3
            source_shell_axis_kinds[shell_id, axis] =
                :parent_lattice_support_state
            source_shell_axis_starts[shell_id, axis] = bounds.starts[axis]
            source_shell_axis_stops[shell_id, axis] = bounds.stops[axis]
            source_shell_fixed_axis_indices[shell_id, axis] = 0
            source_shell_contracted_dims[shell_id, axis] = bounds.dims[axis]
        end
        for (local_col, support_state) in enumerate(support_states)
            local_tuple = ntuple(
                axis -> support_state[axis] - bounds.starts[axis] + 1,
                3,
            )
            mode_source_shell_ids[mode_row] = shell_id
            mode_indices[mode_row] = local_col
            mode_unit_labels[mode_row] = unit.role
            native_source_id_labels[mode_row] =
                _pqs_source_mode_tuple_label(shell_id, local_tuple)
            local_axis_function_indices[mode_row, :] .= collect(local_tuple)
            source_axis_indices[mode_row, :] .= collect(support_state)
            parent_lattice_axis_indices[mode_row, :] .= collect(support_state)
            source_mode_statuses[mode_row] =
                :native_support_dense_source_support_state
            source_axis_tuple_statuses[mode_row] =
                :native_parent_lattice_support_state
            parent_lattice_axis_statuses[mode_row] =
                :native_parent_lattice_support_state
            mode_row += 1
        end
        shell_row += 1
    end

    for entry in shell_realized_unit_entries
        unit = entry.unit
        metadata = entry.metadata
        shell_id = shell_row

        source_shell_unit_indices[shell_id] = entry.unit_index
        source_shell_unit_labels[shell_id] = unit.role
        source_shell_unit_categories[shell_id] = unit.category
        source_shell_unit_kinds[shell_id] = unit.kind
        source_shell_retained_starts[shell_id] = first(unit.column_range)
        source_shell_retained_stops[shell_id] = last(unit.column_range)
        source_shell_labels[shell_id] = unit.role
        source_shell_statuses[shell_id] =
            :native_shell_realized_boundary_source_box
        source_shell_construction_kinds[shell_id] =
            :shell_realized_pqs_fixture
        source_shell_mode_counts[shell_id] = length(metadata.mode_indices)
        source_shell_mode_orderings[shell_id] = :boundary_mode_indices_order
        for axis in 1:3
            interval = metadata.axis_intervals[axis]
            source_shell_axis_kinds[shell_id, axis] = :raw_box_axis
            source_shell_axis_starts[shell_id, axis] = first(interval)
            source_shell_axis_stops[shell_id, axis] = last(interval)
            source_shell_fixed_axis_indices[shell_id, axis] = 0
            source_shell_contracted_dims[shell_id, axis] =
                metadata.source_mode_dims[axis]
        end
        for (local_col, mode_tuple) in enumerate(metadata.mode_indices)
            mode_source_shell_ids[mode_row] = shell_id
            mode_indices[mode_row] = local_col
            mode_unit_labels[mode_row] = unit.role
            native_source_id_labels[mode_row] =
                _pqs_source_mode_tuple_label(shell_id, mode_tuple)
            local_axis_function_indices[mode_row, :] .= collect(mode_tuple)
            source_axis_indices[mode_row, :] .= collect(mode_tuple)
            source_mode_statuses[mode_row] =
                :native_shell_realized_boundary_source_mode
            source_axis_tuple_statuses[mode_row] =
                :native_boundary_source_mode_tuple
            mode_row += 1
        end
        shell_row += 1
    end

    shell_row == source_shell_count + 1 || throw(
        AssertionError("source-shell inventory did not fill every source shell"),
    )
    mode_row == source_mode_count + 1 || throw(
        AssertionError("source-shell inventory did not fill every source mode"),
    )
    for row in 1:source_mode_count
        shell_id = mode_source_shell_ids[row]
        all(
            axis -> begin
                local_axis = local_axis_function_indices[row, axis]
                1 <= local_axis <= source_shell_contracted_dims[shell_id, axis]
            end,
            1:3,
        ) || throw(
            ArgumentError("source-mode local labels must be within source-shell contracted dimensions"),
        )
    end

    fixed_dimension = hasproperty(inventory.diagnostics, :fixed_dimension) ?
        Int(inventory.diagnostics.fixed_dimension) : Int(coverage.last_column)
    total_unavailable_unit_count =
        support_dense_unavailable_unit_count +
        shell_realized_unavailable_unit_count +
        other_unavailable_unit_count
    total_unavailable_column_count =
        support_dense_unavailable_column_count +
        shell_realized_unavailable_column_count +
        other_unavailable_column_count
    center_status =
        native_center_mode_count == 0 ?
        :unavailable_missing_native_comx_center_facts :
        native_center_mode_count == source_mode_count ?
        :native_representative :
        :partial_native_representative_product_doside
    center_definition =
        native_center_mode_count == 0 ? :unavailable : :comx_construction
    non_product_source_shell_count =
        support_dense_source_shell_count + shell_realized_source_shell_count
    inventory_status =
        non_product_source_shell_count > 0 ?
        (
            total_unavailable_unit_count == 0 ?
            :native_current_route_source_shell_modes :
            :partial_current_route_source_shell_modes
        ) :
        :product_doside_source_shell_modes_only
    covered_categories = Symbol[]
    product_shell_count > 0 && push!(covered_categories, :product_doside)
    support_dense_source_shell_count > 0 &&
        push!(covered_categories, :support_dense)
    shell_realized_source_shell_count > 0 &&
        push!(covered_categories, :shell_realized_pqs_fixture)
    non_product_source_mode_status =
        non_product_source_shell_count == 0 ?
        :unavailable_missing_native_non_product_source_mode_producer :
        total_unavailable_unit_count == 0 ?
        :native_non_product_source_shell_mode_labels :
        :partial_native_non_product_source_shell_mode_labels
    source_mode_label_status =
        non_product_source_shell_count == 0 ?
        :native_product_doside_source_mode_indices_only :
        :native_source_mode_tuple_labels_shell_ray_radial_unavailable
    available_native_facts =
        native_center_mode_count == 0 ?
        "product/doside retained ranges, staged axes, fixed/active-axis intervals, axis_function_indices, support-dense support-state source tuples when present, and shell-realized PQS boundary source-mode tuples when present" :
        "product/doside retained ranges, staged axes, fixed/active-axis intervals, axis_function_indices, native COMX/source-transform representative center vectors, support-dense support-state source tuples when present, and shell-realized PQS boundary source-mode tuples when present"
    unavailable_native_facts =
        native_center_mode_count == 0 ?
        "native COMX representative centers, shell-realized PQS compact source relations, and relation weights or spans" :
        "native representative centers for non-product source modes, shell-realized PQS compact source relations, and relation weights or spans"
    return (
        object_kind = :pqs_current_route_source_shell_mode_inventory,
        status = inventory_status,
        schema_version = :pqs_source_shell_modes_private_v1,
        inventory_object_kind = inventory.object_kind,
        fixed_dimension = fixed_dimension,
        source_shell_count = source_shell_count,
        source_mode_count = source_mode_count,
        source_shells = (
            source_shell_ids = source_shell_ids,
            unit_indices = source_shell_unit_indices,
            unit_labels = source_shell_unit_labels,
            unit_categories = source_shell_unit_categories,
            unit_kinds = source_shell_unit_kinds,
            retained_starts = source_shell_retained_starts,
            retained_stops = source_shell_retained_stops,
            source_shell_labels = source_shell_labels,
            source_shell_statuses = source_shell_statuses,
            construction_kinds = source_shell_construction_kinds,
            axis_kinds = source_shell_axis_kinds,
            axis_starts = source_shell_axis_starts,
            axis_stops = source_shell_axis_stops,
            fixed_axis_indices = source_shell_fixed_axis_indices,
            contracted_dims = source_shell_contracted_dims,
            source_mode_counts = source_shell_mode_counts,
            source_mode_orderings = source_shell_mode_orderings,
            center_definitions = source_shell_center_definitions,
            center_statuses = source_shell_center_statuses,
            lowdin_correction_applied =
                source_shell_lowdin_correction_applied,
            shell_label_statuses = source_shell_shell_label_statuses,
            ray_label_statuses = source_shell_ray_label_statuses,
            radial_order_statuses = source_shell_radial_order_statuses,
        ),
        source_modes = (
            source_shell_ids = mode_source_shell_ids,
            mode_indices = mode_indices,
            unit_labels = mode_unit_labels,
            native_source_id_labels = native_source_id_labels,
            local_axis_function_indices = local_axis_function_indices,
            source_axis_indices = source_axis_indices,
            parent_lattice_axis_indices = parent_lattice_axis_indices,
            source_mode_statuses = source_mode_statuses,
            source_axis_tuple_statuses = source_axis_tuple_statuses,
            parent_lattice_axis_statuses = parent_lattice_axis_statuses,
            center_coordinates = mode_center_coordinates,
            center_definitions = mode_center_definitions,
            center_statuses = mode_center_statuses,
            lowdin_correction_applied = mode_lowdin_correction_applied,
            shell_label_statuses = mode_shell_label_statuses,
            ray_label_statuses = mode_ray_label_statuses,
            radial_order_statuses = mode_radial_order_statuses,
            inferred_from_centers = inferred_from_centers,
            inferred_from_nearest_grid = inferred_from_nearest_grid,
            inferred_from_support_order = inferred_from_support_order,
            inferred_from_support_indices = inferred_from_support_indices,
            inferred_from_raw_to_final_support =
                inferred_from_raw_to_final_support,
        ),
        center_status = center_status,
        covered_unit_categories = Tuple(covered_categories),
        non_product_source_mode_status = non_product_source_mode_status,
        source_mode_label_status = source_mode_label_status,
        available_native_facts = available_native_facts,
        unavailable_native_facts = unavailable_native_facts,
        absences_by_contract = (
            repo_ray_grouping_policy = true,
            product_axis_tuples_not_interpreted_as_ray_labels = true,
            representative_centers_as_identity_labels = true,
            native_comx_centers = center_status != :native_representative,
            support_dense_source_shell_modes =
                support_dense_unavailable_unit_count > 0,
            shell_realized_pqs_source_relations = true,
            lowdin_mixture_weights_or_spans = true,
            coordinate_or_nearest_grid_reconstruction = true,
            support_row_order_or_support_index_inference = true,
            raw_to_final_support_inference = true,
            retained_weight_ida_division = true,
            route_hamiltonian_adoption = true,
        ),
        provenance = merge(
            (
                source = :pqs_current_route_source_shell_mode_inventory,
                inventory_object_kind = inventory.object_kind,
            ),
            provenance,
        ),
        diagnostics = (
            source = :pqs_current_route_source_shell_mode_inventory,
            private_reporting_only = true,
            product_doside_source_shell_modes_only =
                inventory_status == :product_doside_source_shell_modes_only,
            repo_exports_native_facts_not_ray_policy = true,
            fixed_dimension = fixed_dimension,
            source_shell_count = source_shell_count,
            source_mode_count = source_mode_count,
            product_doside_source_shell_count = product_shell_count,
            product_doside_source_mode_count = product_mode_count,
            support_dense_source_shell_count = support_dense_source_shell_count,
            support_dense_source_mode_count = support_dense_source_mode_count,
            shell_realized_pqs_source_shell_count =
                shell_realized_source_shell_count,
            shell_realized_pqs_source_mode_count =
                shell_realized_source_mode_count,
            native_center_shell_count = native_center_shell_count,
            native_center_mode_count = native_center_mode_count,
            support_dense_unavailable_unit_count =
                support_dense_unavailable_unit_count,
            support_dense_unavailable_column_count =
                support_dense_unavailable_column_count,
            shell_realized_pqs_unavailable_unit_count =
                shell_realized_unavailable_unit_count,
            shell_realized_pqs_unavailable_column_count =
                shell_realized_unavailable_column_count,
            other_unavailable_unit_count = other_unavailable_unit_count,
            other_unavailable_column_count = other_unavailable_column_count,
            total_unavailable_unit_count = total_unavailable_unit_count,
            total_unavailable_column_count = total_unavailable_column_count,
            coverage_complete = coverage_complete,
            non_product_source_mode_status = non_product_source_mode_status,
            source_mode_label_status = source_mode_label_status,
            normalized_local_axis_labels = true,
            parent_lattice_axis_coordinates_explicit = true,
            center_status = center_status,
            center_definition = center_definition,
            lowdin_correction_applied = false,
            shell_label_status = :unavailable,
            ray_label_status = :unavailable,
            radial_order_status = :unavailable,
            product_axis_tuples_not_interpreted_as_ray_labels = true,
            product_axis_tuples_are_ray_labels = false,
            representative_centers_are_identity_labels = false,
            inferred_from_centers = false,
            inferred_from_nearest_grid = false,
            inferred_from_support_order = false,
            inferred_from_support_indices = false,
            inferred_from_raw_to_final_support = false,
            route_construction_changed = false,
            construction_mutated = false,
            sidecar_installation = false,
            packet_adoption = false,
            fixed_block_routing = false,
            qwhamiltonian_changed = false,
            hamiltonian_matrix_built = false,
            public_default_consumes = false,
            mwg_ida_semantics_changed = false,
            ecp_terms_implemented = false,
            scf_hf_validation_claim = false,
            cr2_science_status_changed = false,
            retained_weight_or_ida_division = false,
        ),
    )
end

function _pqs_current_route_source_shell_mode_inventory(
    construction::_BondAlignedDiatomicHighOrderRecipeSourceConstruction3D;
    inventory = _pqs_current_route_retained_unit_inventory(construction),
    provenance = (;),
    strict::Bool = true,
)
    return _pqs_current_route_source_shell_mode_inventory(
        inventory;
        provenance,
        strict,
    )
end

"""
    _pqs_source_metadata_export_contract()

Private source metadata sidecar contract for `source_shells` and
`source_modes`.

The tables export construction-native provenance only. `source_shells` rows
identify source boxes/shells and their source-axis metadata. `source_modes`
rows identify native source functions by `(source_shell_id, ix, iy, iz)`,
where `local_axis_*` is normalized to `1:nx`, `1:ny`, and `1:nz` for that
source shell's contracted dimensions. Parent-lattice coordinates are exported
only in the explicit `parent_lattice_axis_*` fields when native facts exist.
The older `source_axis_*` fields are retained as private native-tuple
compatibility columns; consumers should use `local_axis_*` for identity labels
and `parent_lattice_axis_*` for parent/support coordinates.
For shell-realized PQS those parent-lattice fields remain unavailable because
the local labels are source-mode coordinates, not parent lattice coordinates.
These labels are not fixed-column-to-source-mode decomposition relations, do
not define `ray_id`, and do not carry relation weights or spans.
"""
function _pqs_source_metadata_export_contract()
    return (
        schema_version = :pqs_source_shell_modes_private_v1,
        status = :private_source_metadata_export_contract,
        source_shells_header = (
            "source_shell_id",
            "unit_index",
            "unit_label",
            "unit_category",
            "unit_kind",
            "retained_start",
            "retained_stop",
            "source_shell_label",
            "source_shell_status",
            "construction_kind",
            "axis_kind_x",
            "axis_kind_y",
            "axis_kind_z",
            "axis_start_x",
            "axis_start_y",
            "axis_start_z",
            "axis_stop_x",
            "axis_stop_y",
            "axis_stop_z",
            "fixed_axis_index_x",
            "fixed_axis_index_y",
            "fixed_axis_index_z",
            "contracted_dim_x",
            "contracted_dim_y",
            "contracted_dim_z",
            "source_mode_count",
            "source_mode_ordering",
            "center_definition",
            "center_status",
            "lowdin_correction_applied",
            "shell_label_status",
            "ray_label_status",
            "radial_order_status",
        ),
        source_modes_header = (
            "source_shell_id",
            "mode_index",
            "unit_label",
            "native_source_id_label",
            "local_axis_x",
            "local_axis_y",
            "local_axis_z",
            "source_axis_x",
            "source_axis_y",
            "source_axis_z",
            "parent_lattice_axis_x",
            "parent_lattice_axis_y",
            "parent_lattice_axis_z",
            "source_mode_status",
            "source_axis_tuple_status",
            "parent_lattice_axis_status",
            "center_x",
            "center_y",
            "center_z",
            "center_definition",
            "center_status",
            "lowdin_correction_applied",
            "shell_label_status",
            "ray_label_status",
            "radial_order_status",
            "inferred_from_centers",
            "inferred_from_nearest_grid",
            "inferred_from_support_order",
            "inferred_from_support_indices",
            "inferred_from_raw_to_final_support",
        ),
        label_semantics = :construction_native_identifiers_not_relations,
        source_mode_local_axis_semantics =
            :normalized_source_shell_local_coordinates,
        parent_lattice_axis_coordinate_policy =
            :explicit_columns_when_native_available,
        shell_realized_pqs_source_axis_indices =
            :local_native_source_mode_coordinates,
        repo_ray_id_policy = :not_exported,
        relation_weight_span_policy = :not_in_source_metadata_tables,
        retained_weight_ida_division = :forbidden,
        route_operator_public_adoption = :forbidden,
    )
end

function _pqs_source_metadata_export_string(value)
    return isnothing(value) ? "" : string(value)
end

function _pqs_source_metadata_export_tsv_row(io, values)
    println(io, join(_pqs_source_metadata_export_string.(values), '\t'))
    return nothing
end

function _pqs_source_metadata_export_validate_inventory(source_shell_mode_inventory)
    source_shell_mode_inventory.object_kind ==
        :pqs_current_route_source_shell_mode_inventory || throw(
        ArgumentError("source metadata export requires _pqs_current_route_source_shell_mode_inventory output"),
    )
    diagnostics = source_shell_mode_inventory.diagnostics
    !diagnostics.inferred_from_centers || throw(
        ArgumentError("source metadata export cannot include center-inferred labels"),
    )
    !diagnostics.inferred_from_nearest_grid || throw(
        ArgumentError("source metadata export cannot include nearest-grid-inferred labels"),
    )
    !diagnostics.inferred_from_support_order || throw(
        ArgumentError("source metadata export cannot include support-order-inferred labels"),
    )
    !diagnostics.inferred_from_support_indices || throw(
        ArgumentError("source metadata export cannot include support-index-inferred labels"),
    )
    !diagnostics.inferred_from_raw_to_final_support || throw(
        ArgumentError("source metadata export cannot include raw_to_final-inferred labels"),
    )
    !diagnostics.retained_weight_or_ida_division || throw(
        ArgumentError("source metadata export cannot include retained-weight or IDA division"),
    )
    !diagnostics.route_construction_changed || throw(
        ArgumentError("source metadata export cannot change route construction"),
    )
    !diagnostics.packet_adoption || throw(
        ArgumentError("source metadata export cannot adopt packet/fixed-block behavior"),
    )
    !diagnostics.qwhamiltonian_changed || throw(
        ArgumentError("source metadata export cannot change QW/Hamiltonian behavior"),
    )
    !diagnostics.public_default_consumes || throw(
        ArgumentError("source metadata export cannot change public/default routing"),
    )
    return nothing
end

function _write_pqs_source_shells_table(io::IO, source_shell_mode_inventory)
    _pqs_source_metadata_export_validate_inventory(source_shell_mode_inventory)
    contract = _pqs_source_metadata_export_contract()
    source_shells = source_shell_mode_inventory.source_shells
    _pqs_source_metadata_export_tsv_row(io, contract.source_shells_header)
    for index in eachindex(source_shells.source_shell_ids)
        _pqs_source_metadata_export_tsv_row(
            io,
            (
                source_shells.source_shell_ids[index],
                source_shells.unit_indices[index],
                source_shells.unit_labels[index],
                source_shells.unit_categories[index],
                source_shells.unit_kinds[index],
                source_shells.retained_starts[index],
                source_shells.retained_stops[index],
                source_shells.source_shell_labels[index],
                source_shells.source_shell_statuses[index],
                source_shells.construction_kinds[index],
                source_shells.axis_kinds[index, 1],
                source_shells.axis_kinds[index, 2],
                source_shells.axis_kinds[index, 3],
                source_shells.axis_starts[index, 1],
                source_shells.axis_starts[index, 2],
                source_shells.axis_starts[index, 3],
                source_shells.axis_stops[index, 1],
                source_shells.axis_stops[index, 2],
                source_shells.axis_stops[index, 3],
                source_shells.fixed_axis_indices[index, 1],
                source_shells.fixed_axis_indices[index, 2],
                source_shells.fixed_axis_indices[index, 3],
                source_shells.contracted_dims[index, 1],
                source_shells.contracted_dims[index, 2],
                source_shells.contracted_dims[index, 3],
                source_shells.source_mode_counts[index],
                source_shells.source_mode_orderings[index],
                source_shells.center_definitions[index],
                source_shells.center_statuses[index],
                source_shells.lowdin_correction_applied[index],
                source_shells.shell_label_statuses[index],
                source_shells.ray_label_statuses[index],
                source_shells.radial_order_statuses[index],
            ),
        )
    end
    return nothing
end

function _write_pqs_source_shells_table(
    path::AbstractString,
    source_shell_mode_inventory,
)
    open(path, "w") do io
        _write_pqs_source_shells_table(io, source_shell_mode_inventory)
    end
    return nothing
end

function _write_pqs_source_modes_table(io::IO, source_shell_mode_inventory)
    _pqs_source_metadata_export_validate_inventory(source_shell_mode_inventory)
    contract = _pqs_source_metadata_export_contract()
    source_modes = source_shell_mode_inventory.source_modes
    _pqs_source_metadata_export_tsv_row(io, contract.source_modes_header)
    for index in eachindex(source_modes.source_shell_ids)
        _pqs_source_metadata_export_tsv_row(
            io,
            (
                source_modes.source_shell_ids[index],
                source_modes.mode_indices[index],
                source_modes.unit_labels[index],
                source_modes.native_source_id_labels[index],
                source_modes.local_axis_function_indices[index, 1],
                source_modes.local_axis_function_indices[index, 2],
                source_modes.local_axis_function_indices[index, 3],
                source_modes.source_axis_indices[index, 1],
                source_modes.source_axis_indices[index, 2],
                source_modes.source_axis_indices[index, 3],
                source_modes.parent_lattice_axis_indices[index, 1],
                source_modes.parent_lattice_axis_indices[index, 2],
                source_modes.parent_lattice_axis_indices[index, 3],
                source_modes.source_mode_statuses[index],
                source_modes.source_axis_tuple_statuses[index],
                source_modes.parent_lattice_axis_statuses[index],
                source_modes.center_coordinates[index, 1],
                source_modes.center_coordinates[index, 2],
                source_modes.center_coordinates[index, 3],
                source_modes.center_definitions[index],
                source_modes.center_statuses[index],
                source_modes.lowdin_correction_applied[index],
                source_modes.shell_label_statuses[index],
                source_modes.ray_label_statuses[index],
                source_modes.radial_order_statuses[index],
                source_modes.inferred_from_centers[index],
                source_modes.inferred_from_nearest_grid[index],
                source_modes.inferred_from_support_order[index],
                source_modes.inferred_from_support_indices[index],
                source_modes.inferred_from_raw_to_final_support[index],
            ),
        )
    end
    return nothing
end

function _write_pqs_source_modes_table(
    path::AbstractString,
    source_shell_mode_inventory,
)
    open(path, "w") do io
        _write_pqs_source_modes_table(io, source_shell_mode_inventory)
    end
    return nothing
end

function _pqs_current_route_fixed_column_source_relation_inventory(
    inventory;
    provenance = (;),
    strict::Bool = true,
)
    inventory.object_kind == :pqs_current_route_retained_unit_inventory_fixture ||
        throw(
            ArgumentError("fixed-column source relations require _pqs_current_route_retained_unit_inventory output"),
        )
    units = Tuple(inventory.units)
    isempty(units) && throw(
        ArgumentError("fixed-column source relations require at least one retained unit"),
    )
    coverage = inventory.coverage
    coverage_complete = hasproperty(coverage, :covers_every_column_once) ?
        Bool(coverage.covers_every_column_once) : false
    strict && !coverage_complete && throw(
        ArgumentError("fixed-column source relations require complete fixed-side column coverage"),
    )

    product_units = Tuple(
        unit for unit in units
        if unit.category == :product_doside && unit.kind == :product_doside
    )
    product_row_count = 0
    support_dense_unavailable_count = 0
    shell_realized_unavailable_count = 0
    other_unavailable_count = 0
    for unit in units
        unit_count = length(unit.column_range)
        if unit.category == :product_doside && unit.kind == :product_doside
            product_row_count += unit_count
        elseif unit.category == :support_dense
            support_dense_unavailable_count += unit_count
        elseif unit.category == :shell_realized_pqs_fixture
            shell_realized_unavailable_count += unit_count
        else
            other_unavailable_count += unit_count
        end
    end

    fixed_cols = Vector{Int}(undef, product_row_count)
    relation_indices = ones(Int, product_row_count)
    relation_kinds = fill(:product_axis_tuple, product_row_count)
    source_unit_labels = Vector{Symbol}(undef, product_row_count)
    source_mode_labels = Vector{String}(undef, product_row_count)
    source_axis_indices = zeros(Int, product_row_count, 3)
    local_axis_function_indices = zeros(Int, product_row_count, 3)
    relation_statuses = fill(:native_product_axis_tuple, product_row_count)
    shell_label_statuses = fill(:unavailable, product_row_count)
    ray_label_statuses = fill(:unavailable, product_row_count)
    radial_order_statuses = fill(:unavailable, product_row_count)
    coefficient_statuses = fill(:unavailable, product_row_count)
    weight_statuses = fill(:unavailable, product_row_count)
    span_statuses = fill(:unavailable, product_row_count)
    inferred_from_centers = falses(product_row_count)
    inferred_from_nearest_grid = falses(product_row_count)
    inferred_from_support_order = falses(product_row_count)
    inferred_from_support_indices = falses(product_row_count)
    inferred_from_raw_to_final_support = falses(product_row_count)

    row = 1
    for unit in product_units
        hasproperty(unit, :staged_unit) && !isnothing(unit.staged_unit) ||
            throw(
                ArgumentError("product/doside fixed-column source relations require staged-axis metadata"),
            )
        staged_unit = unit.staged_unit
        hasproperty(staged_unit, :axes) && hasproperty(staged_unit, :axis_function_indices) ||
            throw(
                ArgumentError("product/doside fixed-column source relations require axes and axis_function_indices"),
            )
        length(staged_unit.axis_function_indices) == length(unit.column_range) ||
            throw(
                ArgumentError("product/doside relation axis tuple count must match unit retained range"),
            )
        length(staged_unit.axes) == 3 || throw(
            ArgumentError("product/doside relation staged axes must be three-dimensional"),
        )
        for (local_col, fixed_col) in enumerate(unit.column_range)
            local_tuple = staged_unit.axis_function_indices[local_col]
            source_tuple = ntuple(
                axis -> _pqs_current_route_product_axis_source_index(
                    staged_unit.axes[axis],
                    local_tuple[axis],
                ),
                3,
            )
            fixed_cols[row] = fixed_col
            source_unit_labels[row] = unit.role
            source_mode_labels[row] = string(
                "product_axis_tuple:",
                source_tuple[1],
                ",",
                source_tuple[2],
                ",",
                source_tuple[3],
            )
            source_axis_indices[row, :] .= collect(source_tuple)
            local_axis_function_indices[row, :] .= collect(local_tuple)
            row += 1
        end
    end

    fixed_dimension = hasproperty(inventory.diagnostics, :fixed_dimension) ?
        Int(inventory.diagnostics.fixed_dimension) : Int(coverage.last_column)
    total_unavailable_count =
        support_dense_unavailable_count +
        shell_realized_unavailable_count +
        other_unavailable_count
    return (
        object_kind = :pqs_current_route_fixed_column_source_relation_inventory,
        status = :product_doside_axis_tuple_relations_only,
        schema_version = :pqs_fixed_column_source_relations_private_v1,
        inventory_object_kind = inventory.object_kind,
        fixed_dimension = fixed_dimension,
        row_count = product_row_count,
        relation_rows_available = product_row_count > 0,
        fixed_cols = fixed_cols,
        relation_indices = relation_indices,
        relation_kinds = relation_kinds,
        source_unit_labels = source_unit_labels,
        source_mode_labels = source_mode_labels,
        source_axis_indices = source_axis_indices,
        local_axis_function_indices = local_axis_function_indices,
        relation_statuses = relation_statuses,
        shell_label_statuses = shell_label_statuses,
        ray_label_statuses = ray_label_statuses,
        radial_order_statuses = radial_order_statuses,
        coefficient_statuses = coefficient_statuses,
        weight_statuses = weight_statuses,
        span_statuses = span_statuses,
        inferred_from_centers = inferred_from_centers,
        inferred_from_nearest_grid = inferred_from_nearest_grid,
        inferred_from_support_order = inferred_from_support_order,
        inferred_from_support_indices = inferred_from_support_indices,
        inferred_from_raw_to_final_support = inferred_from_raw_to_final_support,
        covered_unit_categories = (:product_doside,),
        non_product_relation_status =
            :unavailable_missing_native_non_product_relation_producer,
        relation_label_status =
            :native_product_axis_tuple_only_shell_ray_radial_unavailable,
        relation_weight_status = :unavailable,
        relation_span_status = :unavailable,
        missing_producer =
            :construction_native_non_product_shell_ray_relation_producer,
        available_native_facts =
            "product/doside retained ranges, staged axes, fixed/active-axis intervals, and axis_function_indices",
        unavailable_native_facts =
            "support-dense atom-box shell/ray/radial labels, shell-realized PQS compact source relations, and relation weights or spans",
        absences_by_contract = (
            product_axis_tuples_not_interpreted_as_ray_labels = true,
            shell_ray_radial_labels_for_product_rows = true,
            support_dense_relation_rows = true,
            shell_realized_pqs_relation_rows = true,
            relation_weights_or_spans = true,
            coordinate_or_nearest_grid_reconstruction = true,
            support_row_order_or_support_index_inference = true,
            raw_to_final_support_inference = true,
            retained_weight_ida_division = true,
            route_hamiltonian_adoption = true,
        ),
        provenance = merge(
            (
                source = :pqs_current_route_fixed_column_source_relation_inventory,
                inventory_object_kind = inventory.object_kind,
            ),
            provenance,
        ),
        diagnostics = (
            source = :pqs_current_route_fixed_column_source_relation_inventory,
            private_reporting_only = true,
            product_doside_relation_rows_only = true,
            product_axis_tuples_not_interpreted_as_ray_labels = true,
            product_axis_tuples_are_ray_labels = false,
            fixed_dimension = fixed_dimension,
            row_count = product_row_count,
            product_doside_row_count = product_row_count,
            support_dense_unavailable_column_count =
                support_dense_unavailable_count,
            shell_realized_pqs_unavailable_column_count =
                shell_realized_unavailable_count,
            other_unavailable_column_count = other_unavailable_count,
            total_unavailable_column_count = total_unavailable_count,
            coverage_complete = coverage_complete,
            shell_label_status = :unavailable,
            ray_label_status = :unavailable,
            radial_order_status = :unavailable,
            coefficient_status = :unavailable,
            weight_status = :unavailable,
            span_status = :unavailable,
            inferred_from_centers = false,
            inferred_from_nearest_grid = false,
            inferred_from_support_order = false,
            inferred_from_support_indices = false,
            inferred_from_raw_to_final_support = false,
            route_construction_changed = false,
            construction_mutated = false,
            sidecar_installation = false,
            packet_adoption = false,
            fixed_block_routing = false,
            qwhamiltonian_changed = false,
            hamiltonian_matrix_built = false,
            public_default_consumes = false,
            mwg_ida_semantics_changed = false,
            ecp_terms_implemented = false,
            scf_hf_validation_claim = false,
            cr2_science_status_changed = false,
            retained_weight_or_ida_division = false,
        ),
    )
end

function _pqs_current_route_fixed_column_source_relation_inventory(
    construction::_BondAlignedDiatomicHighOrderRecipeSourceConstruction3D;
    inventory = _pqs_current_route_retained_unit_inventory(construction),
    provenance = (;),
    strict::Bool = true,
)
    return _pqs_current_route_fixed_column_source_relation_inventory(
        inventory;
        provenance,
        strict,
    )
end

function _pqs_fixed_side_retained_unit_metadata_record(unit)
    category = _pqs_fixed_side_unit_get(unit, :category, :unknown)
    role = _pqs_fixed_side_unit_get(unit, :role, :unknown)
    column_range = _pqs_fixed_side_unit_get(unit, :column_range, 1:0)
    support_count = _pqs_fixed_side_unit_get(
        unit,
        :support_count,
        hasproperty(unit, :support_indices) ? length(unit.support_indices) : nothing,
    )
    source_mode_dims = _pqs_fixed_side_unit_source_mode_dims(unit)
    shell_transform =
        _pqs_fixed_side_unit_get(unit, :shell_realization_transform_fact, nothing)
    source_box_operator_application_ready =
        !isnothing(shell_transform) &&
        hasproperty(shell_transform, :source_box_operator_application_ready) ?
        Bool(shell_transform.source_box_operator_application_ready) : false
    compact_source_space_transform_available =
        !isnothing(shell_transform) &&
        hasproperty(shell_transform, :compact_source_space_transform) &&
        hasproperty(shell_transform.compact_source_space_transform, :available) ?
        Bool(shell_transform.compact_source_space_transform.available) : false
    is_shell_realized = category == :shell_realized_pqs_fixture

    return (
        unit_key = role,
        stable_unit_label = role,
        role = role,
        original_role = _pqs_fixed_side_unit_get(unit, :original_role, role),
        category = category,
        kind = _pqs_fixed_side_unit_get(unit, :kind, :unknown),
        unit_class = _pqs_fixed_side_unit_class(category),
        retained_range = column_range,
        retained_count = _pqs_fixed_side_unit_get(
            unit,
            :retained_count,
            length(column_range),
        ),
        support_count = support_count,
        source_dimensions = source_mode_dims,
        source_mode_dims = source_mode_dims,
        source_dimension = _pqs_fixed_side_unit_source_dimension(unit),
        primitive_family =
            _pqs_fixed_side_unit_provenance_get(unit, :primitive_family, nothing),
        representation_kind =
            _pqs_fixed_side_unit_get(unit, :active_representation_stage, nothing),
        support_source_semantics =
            _pqs_fixed_side_unit_get(unit, :support_source_semantics, nothing),
        safe_term_capability =
            _pqs_fixed_side_unit_get(unit, :safe_term_capability, nothing),
        is_product_doside = category == :product_doside,
        is_support_dense_direct_support = category == :support_dense,
        is_shell_realized_pqs_fixture = is_shell_realized,
        shell_realized_pqs_metadata_oracle_fixture = is_shell_realized,
        shell_realized_pqs_source_box_operator_ready =
            source_box_operator_application_ready,
        compact_source_space_transform_available =
            compact_source_space_transform_available,
        shell_row_oracle_only =
            _pqs_fixed_side_unit_diagnostic_get(
                unit,
                :shell_row_oracle_only,
                false,
            ),
        support_local_oracle_used =
            _pqs_fixed_side_unit_diagnostic_get(
                unit,
                :support_local_oracle_used,
                false,
            ),
        raw_product_box_operator_contract =
            _pqs_fixed_side_unit_get(
                unit,
                :raw_product_box_operator_contract,
                false,
            ),
        retained_weight_semantics =
            _pqs_fixed_side_unit_diagnostic_get(
                unit,
                :retained_weight_semantics,
                :not_positive_quadrature_weights,
            ),
        ida_weight_division_allowed =
            _pqs_fixed_side_unit_diagnostic_get(
                unit,
                :ida_weight_division_allowed,
                false,
            ),
        route_descriptor_emitted =
            _pqs_fixed_side_unit_get(unit, :route_descriptor_emitted, false),
        construction_mutated =
            _pqs_fixed_side_unit_get(unit, :construction_mutated, false),
        sidecar_installation =
            _pqs_fixed_side_unit_get(unit, :sidecar_installation, false),
        packet_adoption =
            _pqs_fixed_side_unit_get(unit, :packet_adoption, false),
    )
end

function _pqs_unavailable_symbol_labels(row_count::Integer)
    labels = Vector{Union{Nothing,Symbol}}(undef, Int(row_count))
    fill!(labels, nothing)
    return labels
end

function _pqs_current_route_fixed_column_label_inventory(
    inventory;
    provenance = (;),
    strict::Bool = true,
)
    inventory.object_kind == :pqs_current_route_retained_unit_inventory_fixture ||
        throw(
            ArgumentError("fixed-column label inventory requires _pqs_current_route_retained_unit_inventory output"),
        )
    units = Tuple(inventory.units)
    isempty(units) && throw(
        ArgumentError("fixed-column label inventory requires at least one retained unit"),
    )
    coverage = inventory.coverage
    coverage_complete = hasproperty(coverage, :covers_every_column_once) ?
        Bool(coverage.covers_every_column_once) : false
    strict && !coverage_complete && throw(
        ArgumentError("fixed-column label inventory requires complete fixed-side column coverage"),
    )

    records = Tuple(_pqs_fixed_side_retained_unit_metadata_record(unit) for unit in units)
    fixed_dimension = hasproperty(inventory.diagnostics, :fixed_dimension) ?
        Int(inventory.diagnostics.fixed_dimension) : Int(coverage.last_column)
    row_count = sum(record.retained_count for record in records)
    fixed_cols = Vector{Int}(undef, row_count)
    unit_indices = Vector{Int}(undef, row_count)
    unit_labels = Vector{Symbol}(undef, row_count)
    unit_categories = Vector{Symbol}(undef, row_count)
    unit_kinds = Vector{Symbol}(undef, row_count)
    unit_retained_starts = Vector{Int}(undef, row_count)
    unit_retained_stops = Vector{Int}(undef, row_count)
    source_region_labels = Vector{Symbol}(undef, row_count)
    source_region_label_statuses =
        fill(:retained_unit_region_label, row_count)
    source_box_labels = _pqs_unavailable_symbol_labels(row_count)
    source_box_label_statuses = fill(:unavailable, row_count)
    owner_labels = _pqs_unavailable_symbol_labels(row_count)
    owner_label_statuses = fill(:unavailable, row_count)
    shell_label_statuses = fill(:unavailable, row_count)
    shell_indices = zeros(Int, row_count)
    ray_label_statuses = fill(:unavailable, row_count)
    ray_ids = _pqs_unavailable_symbol_labels(row_count)
    ray_family_labels = _pqs_unavailable_symbol_labels(row_count)
    radial_order_statuses = fill(:unavailable, row_count)
    radial_orders = zeros(Int, row_count)
    inferred_from_centers = falses(row_count)
    inferred_from_nearest_grid = falses(row_count)
    inferred_from_support_order = falses(row_count)
    inferred_from_support_indices = falses(row_count)
    inferred_from_raw_to_final_support = falses(row_count)

    row = 1
    for (unit_index, record) in enumerate(records)
        for fixed_col in record.retained_range
            fixed_cols[row] = fixed_col
            unit_indices[row] = unit_index
            unit_labels[row] = record.unit_key
            unit_categories[row] = record.category
            unit_kinds[row] = record.kind
            unit_retained_starts[row] = first(record.retained_range)
            unit_retained_stops[row] = last(record.retained_range)
            source_region_labels[row] = record.unit_key
            row += 1
        end
    end

    fixed_cols_cover = fixed_cols == collect(1:fixed_dimension)
    strict && row_count != fixed_dimension && throw(
        ArgumentError("fixed-column label inventory row count does not match fixed dimension"),
    )
    strict && !fixed_cols_cover && throw(
        ArgumentError("fixed-column label inventory does not cover 1:fixed_dimension in order"),
    )
    unit_ranges_match_inventory = all(
        index -> begin
            record = records[unit_indices[index]]
            unit_labels[index] == record.unit_key &&
                unit_categories[index] == record.category &&
                unit_kinds[index] == record.kind &&
                unit_retained_starts[index] == first(record.retained_range) &&
                unit_retained_stops[index] == last(record.retained_range) &&
                unit_retained_starts[index] <= fixed_cols[index] <=
                unit_retained_stops[index]
        end,
        eachindex(fixed_cols),
    )
    shell_realized_rows = findall(
        ==(:shell_realized_pqs_fixture),
        unit_categories,
    )
    shell_realized_unavailable =
        all(row -> shell_label_statuses[row] == :unavailable, shell_realized_rows) &&
        all(row -> ray_label_statuses[row] == :unavailable, shell_realized_rows) &&
        all(
            row -> radial_order_statuses[row] == :unavailable,
            shell_realized_rows,
        )

    return (
        object_kind = :pqs_current_route_fixed_column_label_inventory,
        status = :private_fixed_column_label_inventory,
        schema_version = :pqs_fixed_column_labels_private_v1,
        inventory_object_kind = inventory.object_kind,
        fixed_dimension = fixed_dimension,
        row_count = row_count,
        fixed_cols = fixed_cols,
        unit_indices = unit_indices,
        unit_labels = unit_labels,
        unit_categories = unit_categories,
        unit_kinds = unit_kinds,
        unit_retained_starts = unit_retained_starts,
        unit_retained_stops = unit_retained_stops,
        source_region_labels = source_region_labels,
        source_region_label_statuses = source_region_label_statuses,
        source_box_labels = source_box_labels,
        source_box_label_statuses = source_box_label_statuses,
        owner_labels = owner_labels,
        owner_label_statuses = owner_label_statuses,
        shell_label_statuses = shell_label_statuses,
        shell_indices = shell_indices,
        ray_label_statuses = ray_label_statuses,
        ray_ids = ray_ids,
        ray_family_labels = ray_family_labels,
        radial_order_statuses = radial_order_statuses,
        radial_orders = radial_orders,
        inferred_from_centers = inferred_from_centers,
        inferred_from_nearest_grid = inferred_from_nearest_grid,
        inferred_from_support_order = inferred_from_support_order,
        inferred_from_support_indices = inferred_from_support_indices,
        inferred_from_raw_to_final_support = inferred_from_raw_to_final_support,
        label_status = (
            source_region_label_status = :retained_unit_region_label,
            source_box_label_status = :unavailable,
            owner_label_status = :unavailable,
            shell_label_status = :unavailable,
            ray_label_status = :unavailable,
            radial_order_status = :unavailable,
        ),
        absences_by_contract = (
            source_box_labels = true,
            owner_labels = true,
            shell_labels = true,
            ray_labels = true,
            radial_order_labels = true,
            coordinate_or_nearest_grid_reconstruction = true,
            support_row_order_or_support_index_inference = true,
            raw_to_final_support_inference = true,
            retained_weight_ida_division = true,
            route_hamiltonian_adoption = true,
        ),
        provenance = merge(
            (
                source = :pqs_current_route_fixed_column_label_inventory,
                inventory_object_kind = inventory.object_kind,
            ),
            provenance,
        ),
        diagnostics = (
            source = :pqs_current_route_fixed_column_label_inventory,
            private_reporting_only = true,
            fixed_dimension = fixed_dimension,
            row_count = row_count,
            coverage_complete = coverage_complete,
            fixed_cols_cover_1_to_fixed_dimension = fixed_cols_cover,
            unit_ranges_match_inventory = unit_ranges_match_inventory,
            source_region_labels_match_unit_labels =
                source_region_labels == unit_labels &&
                all(
                    ==(:retained_unit_region_label),
                    source_region_label_statuses,
                ),
            product_doside_row_count =
                count(==(:product_doside), unit_categories),
            support_dense_row_count =
                count(==(:support_dense), unit_categories),
            shell_realized_pqs_row_count = length(shell_realized_rows),
            shell_realized_pqs_shell_ray_radial_unavailable =
                shell_realized_unavailable,
            source_box_label_status = :unavailable,
            owner_label_status = :unavailable,
            shell_label_status = :unavailable,
            ray_label_status = :unavailable,
            radial_order_status = :unavailable,
            inferred_from_centers = false,
            inferred_from_nearest_grid = false,
            inferred_from_support_order = false,
            inferred_from_support_indices = false,
            inferred_from_raw_to_final_support = false,
            route_construction_changed = false,
            construction_mutated = false,
            sidecar_installation = false,
            packet_adoption = false,
            fixed_block_routing = false,
            qwhamiltonian_changed = false,
            hamiltonian_matrix_built = false,
            public_default_consumes = false,
            mwg_ida_semantics_changed = false,
            ecp_terms_implemented = false,
            scf_hf_validation_claim = false,
            cr2_science_status_changed = false,
            retained_weight_or_ida_division = false,
        ),
    )
end

function _pqs_current_route_fixed_column_label_inventory(
    construction::_BondAlignedDiatomicHighOrderRecipeSourceConstruction3D;
    inventory = _pqs_current_route_retained_unit_inventory(construction),
    provenance = (;),
    strict::Bool = true,
)
    return _pqs_current_route_fixed_column_label_inventory(
        inventory;
        provenance,
        strict,
    )
end

function _pqs_current_route_fixed_side_retained_unit_metadata(
    inventory;
    provenance = (;),
    strict::Bool = true,
)
    inventory.object_kind == :pqs_current_route_retained_unit_inventory_fixture ||
        throw(
            ArgumentError("fixed-side retained-unit metadata requires _pqs_current_route_retained_unit_inventory output"),
        )
    units = Tuple(inventory.units)
    isempty(units) && throw(
        ArgumentError("fixed-side retained-unit metadata requires at least one retained unit"),
    )
    records = Tuple(_pqs_fixed_side_retained_unit_metadata_record(unit) for unit in units)
    labels = Tuple(record.unit_key for record in records)
    coverage = inventory.coverage
    diagnostics = inventory.diagnostics
    fixed_dimension = hasproperty(diagnostics, :fixed_dimension) ?
        Int(diagnostics.fixed_dimension) : Int(coverage.last_column)
    coverage_complete = hasproperty(coverage, :covers_every_column_once) ?
        Bool(coverage.covers_every_column_once) : false
    shell_records = Tuple(
        record for record in records if record.is_shell_realized_pqs_fixture
    )
    source_box_ready_count = count(
        record -> record.shell_realized_pqs_source_box_operator_ready,
        shell_records,
    )

    if strict
        coverage_complete || throw(
            ArgumentError("fixed-side retained-unit metadata requires complete fixed-side column coverage"),
        )
        all(record -> !record.route_descriptor_emitted, records) || throw(
            ArgumentError("fixed-side retained-unit metadata cannot include route descriptor emission"),
        )
        all(record -> !record.construction_mutated, records) || throw(
            ArgumentError("fixed-side retained-unit metadata cannot include construction mutation"),
        )
        all(record -> !record.sidecar_installation, records) || throw(
            ArgumentError("fixed-side retained-unit metadata cannot include sidecar installation"),
        )
        all(record -> !record.packet_adoption, records) || throw(
            ArgumentError("fixed-side retained-unit metadata cannot include packet adoption"),
        )
        all(record -> !record.ida_weight_division_allowed, records) || throw(
            ArgumentError("fixed-side retained-unit metadata cannot include retained-weight/IDA division"),
        )
        source_box_ready_count == 0 || throw(
            ArgumentError("fixed-side retained-unit metadata cannot mark shell-realized PQS fixtures as source-box-operator-ready without an explicit framework update"),
        )
    end

    return (
        object_kind = :pqs_current_route_fixed_side_retained_unit_metadata,
        status = :private_fixed_side_retained_unit_metadata,
        schema_version = :pqs_fixed_side_retained_units_private_v1,
        inventory_object_kind = inventory.object_kind,
        fixed_dimension = fixed_dimension,
        unit_count = length(records),
        retained_units = records,
        labels = (
            source_unit_label_status = :explicit_inventory_unit_keys,
            source_unit_labels = labels,
            shell_label_status = :unavailable,
            shell_labels = (),
            label_reconstruction_from_centers = false,
            nearest_grid_or_center_label_heuristic = false,
        ),
        absences_by_contract = (
            route_construction = true,
            packet_fixed_block_qw_hamiltonian_adoption = true,
            mwg_ida_semantic_change = true,
            ecp_scf_hf_cr2_science_claim = true,
            retained_weight_ida_division = true,
            shell_label_reconstruction_from_centers = true,
            nearest_grid_or_center_label_heuristic = true,
        ),
        provenance = merge(
            (
                source = :pqs_current_route_fixed_side_retained_unit_metadata,
                inventory_object_kind = inventory.object_kind,
            ),
            provenance,
        ),
        diagnostics = (
            source = :pqs_current_route_fixed_side_retained_unit_metadata,
            private_reporting_only = true,
            fixed_side_records_are_explicit_metadata = true,
            unit_count = length(records),
            fixed_dimension = fixed_dimension,
            coverage_complete = coverage_complete,
            first_column = coverage.first_column,
            last_column = coverage.last_column,
            represented_count = coverage.represented_count,
            source_unit_label_status = :explicit_inventory_unit_keys,
            shell_label_status = :unavailable,
            label_reconstruction_from_centers = false,
            nearest_grid_or_center_label_heuristic = false,
            shell_realized_pqs_fixture_count = length(shell_records),
            shell_realized_pqs_source_box_operator_ready_count =
                source_box_ready_count,
            shell_realized_pqs_fixtures_are_metadata_oracle_only =
                all(
                    record -> record.shell_realized_pqs_metadata_oracle_fixture &&
                              !record.shell_realized_pqs_source_box_operator_ready,
                    shell_records,
                ),
            route_construction_changed = false,
            packet_adoption = false,
            fixed_block_routing = false,
            qwhamiltonian_changed = false,
            hamiltonian_matrix_built = false,
            public_default_consumes = false,
            mwg_ida_semantics_changed = false,
            ecp_terms_implemented = false,
            scf_hf_validation_claim = false,
            cr2_science_status_changed = false,
            retained_weight_or_ida_division = false,
        ),
    )
end

function _pqs_current_route_retained_pair_policy(left, right)
    categories = (left.category, right.category)
    if categories == (:product_doside, :product_doside)
        return (
            pair_group = :product_product,
            policy = :product_doside_source_box_path,
            active_current_route = true,
            active_algorithmic_policy = true,
            source_box_algorithm_available = true,
            support_local_oracle_used = false,
            shell_row_oracle_only = false,
            raw_box_pqs_active_pair_policy = false,
        )
    elseif categories == (:support_dense, :support_dense)
        return (
            pair_group = :support_support,
            policy = :support_local_fallback,
            active_current_route = true,
            active_algorithmic_policy = true,
            source_box_algorithm_available = false,
            support_local_oracle_used = false,
            shell_row_oracle_only = false,
            raw_box_pqs_active_pair_policy = false,
        )
    elseif (:support_dense in categories) && (:product_doside in categories)
        return (
            pair_group = :support_product,
            policy = :support_local_fallback,
            active_current_route = true,
            active_algorithmic_policy = true,
            source_box_algorithm_available = false,
            support_local_oracle_used = false,
            shell_row_oracle_only = false,
            raw_box_pqs_active_pair_policy = false,
        )
    elseif (:shell_realized_pqs_fixture in categories) &&
           (:product_doside in categories)
        return (
            pair_group = :shell_realized_pqs_product,
            policy = :support_local_oracle_for_shell_realization,
            active_current_route = false,
            active_algorithmic_policy = false,
            source_box_algorithm_available = false,
            support_local_oracle_used = true,
            shell_row_oracle_only = true,
            raw_box_pqs_active_pair_policy = false,
        )
    elseif (:shell_realized_pqs_fixture in categories) &&
           (:support_dense in categories)
        return (
            pair_group = :shell_realized_pqs_support,
            policy = :support_local_oracle_for_shell_realization,
            active_current_route = false,
            active_algorithmic_policy = false,
            source_box_algorithm_available = false,
            support_local_oracle_used = true,
            shell_row_oracle_only = true,
            raw_box_pqs_active_pair_policy = false,
        )
    elseif categories == (:shell_realized_pqs_fixture, :shell_realized_pqs_fixture)
        return (
            pair_group = :shell_realized_pqs_pqs,
            policy = :support_local_oracle_for_shell_realization,
            active_current_route = false,
            active_algorithmic_policy = false,
            source_box_algorithm_available = false,
            support_local_oracle_used = true,
            shell_row_oracle_only = true,
            raw_box_pqs_active_pair_policy = false,
        )
    end
    throw(
        ArgumentError(
            "unsupported PQS current-route retained pair categories $(categories)",
        ),
    )
end

function _pqs_current_route_retained_pair_counts(pairs)
    group_counts = Dict{Symbol,Int}()
    policy_counts = Dict{Symbol,Int}()
    raw_box_pqs_active_pair_policy_count = 0
    algorithmic_policy_count = 0
    source_box_algorithm_available_count = 0
    support_local_oracle_pair_count = 0
    shell_row_oracle_pair_count = 0
    for pair in pairs
        group_counts[pair.pair_group] = get(group_counts, pair.pair_group, 0) + 1
        policy_counts[pair.policy] = get(policy_counts, pair.policy, 0) + 1
        pair.raw_box_pqs_active_pair_policy &&
            (raw_box_pqs_active_pair_policy_count += 1)
        pair.active_algorithmic_policy && (algorithmic_policy_count += 1)
        pair.source_box_algorithm_available &&
            (source_box_algorithm_available_count += 1)
        pair.support_local_oracle_used && (support_local_oracle_pair_count += 1)
        pair.shell_row_oracle_only && (shell_row_oracle_pair_count += 1)
    end
    return (
        pair_count = length(pairs),
        product_product = get(group_counts, :product_product, 0),
        support_support = get(group_counts, :support_support, 0),
        support_product = get(group_counts, :support_product, 0),
        shell_realized_pqs_product =
            get(group_counts, :shell_realized_pqs_product, 0),
        shell_realized_pqs_support =
            get(group_counts, :shell_realized_pqs_support, 0),
        shell_realized_pqs_pqs = get(group_counts, :shell_realized_pqs_pqs, 0),
        raw_box_pqs_active = raw_box_pqs_active_pair_policy_count,
        active_algorithmic_policy = algorithmic_policy_count,
        source_box_algorithm_available = source_box_algorithm_available_count,
        support_local_oracle_for_shell_realization =
            get(policy_counts, :support_local_oracle_for_shell_realization, 0),
        support_local_oracle_pair_count = support_local_oracle_pair_count,
        shell_row_oracle_pair_count = shell_row_oracle_pair_count,
        product_doside_source_box_path =
            get(policy_counts, :product_doside_source_box_path, 0),
        support_local_fallback = get(policy_counts, :support_local_fallback, 0),
    )
end

function _pqs_current_route_retained_pair_inventory(inventory)
    inventory.object_kind == :pqs_current_route_retained_unit_inventory_fixture ||
        throw(
            ArgumentError("PQS current-route pair inventory requires unit inventory fixture"),
        )
    inventory.status == :private_diagnostic_only || throw(
        ArgumentError("PQS current-route pair inventory requires private diagnostic inventory"),
    )
    inventory.coverage.covers_every_column_once || throw(
        ArgumentError("PQS current-route pair inventory requires complete unit coverage"),
    )
    units = inventory.units
    !isempty(units) || throw(
        ArgumentError("PQS current-route pair inventory requires retained units"),
    )

    pairs = Any[]
    pair_index = 0
    for left_index in eachindex(units)
        left = units[left_index]
        for right_index in left_index:length(units)
            right = units[right_index]
            pair_index += 1
            policy = _pqs_current_route_retained_pair_policy(left, right)
            push!(
                pairs,
                (
                    pair_index = pair_index,
                    left_unit_index = left_index,
                    right_unit_index = right_index,
                    left_role = left.role,
                    right_role = right.role,
                    left_category = left.category,
                    right_category = right.category,
                    left_kind = left.kind,
                    right_kind = right.kind,
                    left_column_range = left.column_range,
                    right_column_range = right.column_range,
                    left_retained_count = left.retained_count,
                    right_retained_count = right.retained_count,
                    pair_shape = (left.retained_count, right.retained_count),
                    pair_group = policy.pair_group,
                    policy = policy.policy,
                    active_current_route = policy.active_current_route,
                    active_algorithmic_policy = policy.active_algorithmic_policy,
                    source_box_algorithm_available =
                        policy.source_box_algorithm_available,
                    support_local_oracle_used = policy.support_local_oracle_used,
                    shell_row_oracle_only = policy.shell_row_oracle_only,
                    raw_box_pqs_active_pair_policy =
                        policy.raw_box_pqs_active_pair_policy,
                ),
            )
        end
    end
    pair_tuple = Tuple(pairs)
    counts = _pqs_current_route_retained_pair_counts(pair_tuple)
    expected_pair_count = div(length(units) * (length(units) + 1), 2)
    counts.pair_count == expected_pair_count || throw(
        ArgumentError("PQS current-route pair inventory pair count mismatch"),
    )
    counts.raw_box_pqs_active == 0 || throw(
        ArgumentError("PQS current-route pair inventory must not activate raw-box PQS policy"),
    )
    diagnostics = (
        source = :pqs_current_route_retained_pair_inventory,
        private_diagnostic_only = true,
        current_route_pair_inventory = true,
        unit_inventory_complete = inventory.coverage.covers_every_column_once,
        upper_triangular_pairs = true,
        pair_count = counts.pair_count,
        expected_pair_count = expected_pair_count,
        unit_count = length(units),
        raw_box_pqs_active_pair_policy_count = counts.raw_box_pqs_active,
        active_algorithmic_policy_pair_count = counts.active_algorithmic_policy,
        source_box_algorithm_available_pair_count =
            counts.source_box_algorithm_available,
        support_local_oracle_for_shell_realization_pair_count =
            counts.support_local_oracle_for_shell_realization,
        shell_row_oracle_pair_count = counts.shell_row_oracle_pair_count,
        shell_realized_pqs_pairs_are_oracle_only =
            counts.shell_row_oracle_pair_count ==
            counts.shell_realized_pqs_product +
            counts.shell_realized_pqs_support +
            counts.shell_realized_pqs_pqs,
        route_descriptor_emitted = false,
        construction_mutated = false,
        sidecar_installation = false,
        packet_adoption = false,
        fixed_block_construction_changed = false,
        qwhamiltonian_changed = false,
        ida_weight_division_allowed = false,
        retained_weight_semantics = :not_positive_quadrature_weights,
        local_ecp_gaussian_mwg_interaction_changed = false,
        whole_route_safe_term_matrix_consumer = false,
    )
    return (
        object_kind = :pqs_current_route_retained_pair_inventory_fixture,
        status = :private_diagnostic_only,
        unit_inventory = inventory,
        pairs = pair_tuple,
        counts = counts,
        diagnostics = diagnostics,
    )
end

function _pqs_current_route_retained_pair_inventory(
    construction::_BondAlignedDiatomicHighOrderRecipeSourceConstruction3D;
    inventory = _pqs_current_route_retained_unit_inventory(construction),
)
    return _pqs_current_route_retained_pair_inventory(inventory)
end

function _pqs_current_route_safe_term_axis_factor_terms(
    metrics::NamedTuple{(:x,:y,:z)},
    term::Symbol,
)
    term in _PRODUCT_DOSIDE_SOURCE_BOX_REFERENCE_TERMS || throw(
        ArgumentError("PQS current-route safe-term matrix received unsupported term $(term)"),
    )
    return Tuple(
        ntuple(
            axis -> _product_doside_axis_metric_matrix(metrics, axis, factor_kinds[axis]),
            3,
        )
        for factor_kinds in _source_box_separable_term_factor_kinds(term)
    )
end

function _pqs_current_route_inventory_unit_entries(unit)
    coefficients =
        unit.category == :shell_realized_pqs_fixture ?
        unit.support_local_coefficient_matrix :
        unit.coefficient_matrix
    return _support_local_retained_entries(
        unit.column_range,
        unit.support_states,
        coefficients,
    )
end

function _pqs_current_route_cached_unit_entries!(entries_cache, unit_index::Int, unit)
    cached = entries_cache[unit_index]
    !isnothing(cached) && return cached
    entries = _pqs_current_route_inventory_unit_entries(unit)
    entries_cache[unit_index] = entries
    return entries
end

function _pqs_current_route_safe_term_matrices_payload(
    inventory,
    pair_inventory,
    metrics::NamedTuple{(:x,:y,:z)},
    selected_terms::Tuple;
    atol::Real,
)
    pair_inventory.unit_inventory === inventory || throw(
        ArgumentError("PQS current-route safe-term matrices require matching pair inventory"),
    )
    retained_dimension = inventory.coverage.last_column
    retained_dimension == inventory.coverage.represented_count || throw(
        ArgumentError("PQS current-route safe-term matrices require complete one-based coverage"),
    )
    expected_pair_count = div(length(inventory.units) * (length(inventory.units) + 1), 2)
    pair_inventory.counts.pair_count == expected_pair_count || throw(
        ArgumentError("PQS current-route safe-term matrices retained-unit pair count mismatch"),
    )
    pair_inventory.counts.raw_box_pqs_active == 0 || throw(
        ArgumentError("PQS current-route safe-term matrices cannot use active raw-box PQS pair policy"),
    )

    matrices = Dict{Symbol,Matrix{Float64}}()
    oracle_matrices = Dict{Symbol,Matrix{Float64}}()
    term_errors = Dict{Symbol,Float64}()
    entries_cache = Vector{Any}(undef, length(inventory.units))
    fill!(entries_cache, nothing)
    output_finite = true
    product_source_box_pair_count = 0
    support_local_fallback_pair_count = 0
    shell_realized_pqs_oracle_pair_count = 0

    for term in selected_terms
        axis_factor_terms =
            _pqs_current_route_safe_term_axis_factor_terms(metrics, term)
        matrix = zeros(Float64, retained_dimension, retained_dimension)
        oracle_matrix = zeros(Float64, retained_dimension, retained_dimension)
        for pair in pair_inventory.pairs
            left_unit = inventory.units[pair.left_unit_index]
            right_unit = inventory.units[pair.right_unit_index]
            if pair.pair_group == :product_product
                product_reference = _product_doside_source_box_reference_block(
                    left_unit.staged_unit,
                    right_unit.staged_unit,
                    metrics;
                    term,
                    atol,
                )
                block = product_reference.block
                product_source_box_pair_count += 1
            else
                left_entries = _pqs_current_route_cached_unit_entries!(
                    entries_cache,
                    pair.left_unit_index,
                    left_unit,
                )
                right_entries = _pqs_current_route_cached_unit_entries!(
                    entries_cache,
                    pair.right_unit_index,
                    right_unit,
                )
                block = _fallback_staged_separable_sum_block(
                    left_entries,
                    right_entries,
                    axis_factor_terms,
                )
                support_local_fallback_pair_count += 1
                (:shell_realized_pqs_fixture in (pair.left_category, pair.right_category)) &&
                    (shell_realized_pqs_oracle_pair_count += 1)
            end

            left_entries = _pqs_current_route_cached_unit_entries!(
                entries_cache,
                pair.left_unit_index,
                left_unit,
            )
            right_entries = _pqs_current_route_cached_unit_entries!(
                entries_cache,
                pair.right_unit_index,
                right_unit,
            )
            oracle_block = _fallback_staged_separable_sum_block(
                left_entries,
                right_entries,
                axis_factor_terms,
            )
            size(block) == pair.pair_shape || throw(
                DimensionMismatch("PQS current-route safe-term block shape mismatch"),
            )
            size(oracle_block) == pair.pair_shape || throw(
                DimensionMismatch("PQS current-route safe-term oracle block shape mismatch"),
            )
            all(isfinite, block) || (output_finite = false)
            all(isfinite, oracle_block) || (output_finite = false)
            matrix[pair.left_column_range, pair.right_column_range] .= block
            oracle_matrix[pair.left_column_range, pair.right_column_range] .=
                oracle_block
            if pair.left_unit_index != pair.right_unit_index
                matrix[pair.right_column_range, pair.left_column_range] .=
                    transpose(block)
                oracle_matrix[pair.right_column_range, pair.left_column_range] .=
                    transpose(oracle_block)
            end
        end
        all(isfinite, matrix) || (output_finite = false)
        all(isfinite, oracle_matrix) || (output_finite = false)
        term_error = LinearAlgebra.norm(matrix - oracle_matrix, Inf)
        term_errors[term] = term_error
        matrices[term] = matrix
        oracle_matrices[term] = oracle_matrix
    end

    output_finite || throw(
        ArgumentError("PQS current-route safe-term matrices produced non-finite entries"),
    )
    global_max_error = maximum(values(term_errors))
    global_max_error <= Float64(atol) || throw(
        ArgumentError("PQS current-route safe-term matrices disagree with support-local oracle"),
    )

    diagnostics = (
        source = :pqs_current_route_safe_term_matrices,
        private_diagnostic_only = true,
        whole_route_safe_term_matrix_consumer = true,
        retained_dimension = retained_dimension,
        terms_checked = selected_terms,
        pair_count = pair_inventory.counts.pair_count,
        expected_pair_count = expected_pair_count,
        unit_count = length(inventory.units),
        product_product_pair_count = pair_inventory.counts.product_product,
        support_support_pair_count = pair_inventory.counts.support_support,
        support_product_pair_count = pair_inventory.counts.support_product,
        shell_realized_pqs_product_pair_count =
            pair_inventory.counts.shell_realized_pqs_product,
        shell_realized_pqs_support_pair_count =
            pair_inventory.counts.shell_realized_pqs_support,
        shell_realized_pqs_pqs_pair_count =
            pair_inventory.counts.shell_realized_pqs_pqs,
        product_source_box_pair_count = div(
            product_source_box_pair_count,
            length(selected_terms),
        ),
        support_local_fallback_pair_count = div(
            support_local_fallback_pair_count,
            length(selected_terms),
        ),
        support_local_oracle_for_shell_realization_pair_count = div(
            shell_realized_pqs_oracle_pair_count,
            length(selected_terms),
        ),
        shell_row_oracle_only = true,
        source_box_algorithm_available_for_shell_realized_pqs = false,
        support_local_oracle_used = true,
        raw_box_pqs_active_pair_policy_count =
            pair_inventory.counts.raw_box_pqs_active,
        term_errors = term_errors,
        global_max_error = global_max_error,
        finite_output = output_finite,
        support_local_oracle_compared = true,
        support_local_oracle_is_debug_validation = true,
        shell_realized_pqs_pairs_use_oracle_not_algorithm = true,
        raw_box_pqs_active_policy_used = false,
        route_descriptor_emitted = false,
        construction_mutated = false,
        sidecar_installation = false,
        packet_adoption = false,
        fixed_block_construction_changed = false,
        qwhamiltonian_changed = false,
        ida_weight_division_allowed = false,
        retained_weight_semantics = :not_positive_quadrature_weights,
        local_ecp_gaussian_mwg_interaction_changed = false,
    )
    return (
        matrices = matrices,
        oracle_matrices = oracle_matrices,
        term_errors = term_errors,
        global_max_error = global_max_error,
        diagnostics = diagnostics,
    )
end

function _pqs_current_route_safe_term_matrices(
    construction::_BondAlignedDiatomicHighOrderRecipeSourceConstruction3D,
    metrics::NamedTuple{(:x,:y,:z)};
    inventory = _pqs_current_route_retained_unit_inventory(construction),
    pair_inventory = _pqs_current_route_retained_pair_inventory(inventory),
    terms = _PRODUCT_DOSIDE_SOURCE_BOX_REFERENCE_TERMS,
    atol::Real = 1.0e-12,
)
    selected_terms = Tuple(Symbol(term) for term in terms)
    !isempty(selected_terms) || throw(
        ArgumentError("PQS current-route safe-term matrices require at least one term"),
    )
    for term in selected_terms
        term in _PRODUCT_DOSIDE_SOURCE_BOX_REFERENCE_TERMS || throw(
            ArgumentError("PQS current-route safe-term matrix received unsupported term $(term)"),
        )
    end

    timed = @timed _pqs_current_route_safe_term_matrices_payload(
        inventory,
        pair_inventory,
        metrics,
        selected_terms;
        atol,
    )
    payload = timed.value
    diagnostics = merge(
        payload.diagnostics,
        (
            elapsed_seconds = Float64(timed.time),
            allocated_bytes = Int(timed.bytes),
            gc_time_seconds = Float64(timed.gctime),
        ),
    )
    return (
        object_kind = :pqs_current_route_safe_term_matrices_fixture,
        status = :private_diagnostic_only,
        terms = selected_terms,
        matrices = payload.matrices,
        oracle_matrices = payload.oracle_matrices,
        term_errors = payload.term_errors,
        global_max_error = payload.global_max_error,
        inventory = inventory,
        pair_inventory = pair_inventory,
        diagnostics = diagnostics,
    )
end

function _pqs_current_route_authority_matrix_candidate(
    source::Symbol,
    object,
    term::Symbol,
    retained_dimension::Int,
)
    field = term
    isnothing(object) && return (
        available = false,
        source = source,
        field = field,
        matrix = nothing,
        shape = nothing,
        failure = (
            source = source,
            field = field,
            reason = :source_unavailable,
        ),
    )
    !hasproperty(object, field) && return (
        available = false,
        source = source,
        field = field,
        matrix = nothing,
        shape = nothing,
        failure = (
            source = source,
            field = field,
            reason = :field_absent,
        ),
    )

    matrix = getproperty(object, field)
    !(matrix isa AbstractMatrix) && return (
        available = false,
        source = source,
        field = field,
        matrix = nothing,
        shape = nothing,
        failure = (
            source = source,
            field = field,
            reason = :field_is_not_matrix,
            value_type = typeof(matrix),
        ),
    )

    shape = size(matrix)
    expected_shape = (retained_dimension, retained_dimension)
    shape != expected_shape && return (
        available = false,
        source = source,
        field = field,
        matrix = nothing,
        shape = shape,
        failure = (
            source = source,
            field = field,
            reason = :wrong_retained_shape,
            shape = shape,
            expected_shape = expected_shape,
        ),
    )
    !all(isfinite, matrix) && return (
        available = false,
        source = source,
        field = field,
        matrix = nothing,
        shape = shape,
        failure = (
            source = source,
            field = field,
            reason = :nonfinite_authority_matrix,
            shape = shape,
        ),
    )
    return (
        available = true,
        source = source,
        field = field,
        matrix = matrix,
        shape = shape,
        failure = nothing,
    )
end

function _pqs_current_route_authority_matrix(
    candidates,
    term::Symbol,
    retained_dimension::Int,
)
    failures = Any[]
    for candidate in candidates
        result = _pqs_current_route_authority_matrix_candidate(
            candidate.source,
            candidate.object,
            term,
            retained_dimension,
        )
        result.available && return result
        push!(failures, result.failure)
    end
    return (
        available = false,
        source = nothing,
        field = term,
        matrix = nothing,
        shape = nothing,
        failure = (
            term = term,
            reason = :no_authoritative_retained_space_field,
            checked_sources = Tuple(failures),
        ),
    )
end

function _pqs_current_route_safe_term_authority_comparison_payload(
    safe_terms,
    candidates;
    atol::Real,
)
    retained_dimension = safe_terms.diagnostics.retained_dimension
    term_errors = Dict{Symbol,Float64}()
    authority_sources = Dict{Symbol,Symbol}()
    authority_fields = Dict{Symbol,Symbol}()
    authority_shapes = Dict{Symbol,Tuple{Int,Int}}()
    unavailable = Any[]
    compared_terms = Symbol[]

    for term in safe_terms.terms
        !haskey(safe_terms.matrices, term) && begin
            push!(
                unavailable,
                (
                    term = term,
                    reason = :safe_term_matrix_absent,
                    checked_sources = (),
                ),
            )
            continue
        end
        safe_matrix = safe_terms.matrices[term]
        authority = _pqs_current_route_authority_matrix(
            candidates,
            term,
            retained_dimension,
        )
        !authority.available && begin
            push!(
                unavailable,
                (
                    term = term,
                    reason = authority.failure.reason,
                    safe_shape = size(safe_matrix),
                    checked_sources = authority.failure.checked_sources,
                ),
            )
            continue
        end

        error = LinearAlgebra.norm(safe_matrix - authority.matrix, Inf)
        term_errors[term] = error
        authority_sources[term] = authority.source
        authority_fields[term] = authority.field
        authority_shapes[term] = authority.shape
        push!(compared_terms, term)
    end

    required_authority_terms = (:overlap, :kinetic)
    missing_required = Tuple(
        term for term in required_authority_terms if !(term in compared_terms)
    )
    isempty(missing_required) || throw(
        ArgumentError(
            "PQS current-route authority comparison could not find authoritative retained-space fields for $(missing_required)",
        ),
    )

    max_authority_error = maximum(values(term_errors))
    max_authority_error <= atol || throw(
        ArgumentError(
            "PQS current-route safe-term authority comparison exceeded tolerance: max error $(max_authority_error), tolerance $(atol)",
        ),
    )
    diagnostics = (
        private_diagnostic_only = true,
        current_route_safe_term_authority_comparison = true,
        retained_dimension = retained_dimension,
        terms_requested = safe_terms.terms,
        compared_terms = Tuple(compared_terms),
        compared_term_count = length(compared_terms),
        unavailable_terms = Tuple(entry.term for entry in unavailable),
        unavailable_term_count = length(unavailable),
        authoritative_sources = Tuple(unique(values(authority_sources))),
        authority_fixed_block_or_sequence_packet_only = true,
        support_local_oracle_secondary = true,
        support_local_oracle_global_max_error = safe_terms.global_max_error,
        max_authority_error = max_authority_error,
        finite_output = all(
            term -> all(isfinite, safe_terms.matrices[term]),
            safe_terms.terms,
        ),
        construction_mutated = false,
        sidecar_installation = false,
        packet_adoption = false,
        fixed_block_construction_changed = false,
        qwhamiltonian_changed = false,
        ida_weight_division_allowed = false,
        retained_weight_semantics = :not_positive_quadrature_weights,
        local_ecp_gaussian_mwg_interaction_changed = false,
    )
    return (
        terms = safe_terms.terms,
        compared_terms = Tuple(compared_terms),
        unavailable_terms = Tuple(unavailable),
        term_errors = term_errors,
        max_authority_error = max_authority_error,
        authority_sources = authority_sources,
        authority_fields = authority_fields,
        authority_shapes = authority_shapes,
        diagnostics = diagnostics,
    )
end

function _pqs_current_route_safe_term_authority_comparison(
    construction::_BondAlignedDiatomicHighOrderRecipeSourceConstruction3D,
    metrics::NamedTuple{(:x,:y,:z)};
    inventory = _pqs_current_route_retained_unit_inventory(construction),
    pair_inventory = _pqs_current_route_retained_pair_inventory(inventory),
    safe_terms = _pqs_current_route_safe_term_matrices(
        construction,
        metrics;
        inventory,
        pair_inventory,
    ),
    fixed_block = nothing,
    sequence_packet = construction.sequence.packet,
    atol::Real = 1.0e-8,
)
    authoritative_fixed_block = isnothing(fixed_block) ?
        _nested_bond_aligned_diatomic_high_order_recipe_source_fixed_block(construction) :
        fixed_block
    candidates = (
        (source = :fixed_block, object = authoritative_fixed_block),
        (source = :sequence_packet, object = sequence_packet),
    )
    timed = @timed _pqs_current_route_safe_term_authority_comparison_payload(
        safe_terms,
        candidates;
        atol,
    )
    payload = timed.value
    diagnostics = merge(
        payload.diagnostics,
        (
            elapsed_seconds = Float64(timed.time),
            allocated_bytes = Int(timed.bytes),
            gc_time_seconds = Float64(timed.gctime),
        ),
    )
    return (
        object_kind = :pqs_current_route_safe_term_authority_comparison_fixture,
        status = :private_diagnostic_only,
        terms = payload.terms,
        compared_terms = payload.compared_terms,
        unavailable_terms = payload.unavailable_terms,
        term_errors = payload.term_errors,
        max_authority_error = payload.max_authority_error,
        authority_sources = payload.authority_sources,
        authority_fields = payload.authority_fields,
        authority_shapes = payload.authority_shapes,
        safe_terms = safe_terms,
        diagnostics = diagnostics,
    )
end

function _pqs_pqs_product_route_descriptor_diagnostic(
    construction::_BondAlignedDiatomicHighOrderRecipeSourceConstruction3D,
    metrics = nothing;
    supported_terms = _PQS_PRODUCT_SOURCE_BOX_SHADOW_TERMS,
    orthogonality_atol::Real = 1.0e-8,
)
    selected_terms = _pqs_pqs_product_supported_safe_terms(supported_terms)
    descriptors = _pqs_route_staged_descriptors(construction)
    convertibility = _pqs_route_raw_plan_convertibility(
        descriptors,
        construction,
        metrics;
        orthogonality_atol,
    )
    direct_or_support_mismatches =
        _pqs_route_direct_or_support_body_mismatches(construction)
    product_doside_unit_count = 0

    missing = _pqs_route_descriptor_missing_symbols(
        min(length(descriptors), convertibility.converted_count),
        product_doside_unit_count,
    )
    !convertibility.checked && !isempty(descriptors) &&
        push!(missing, :axis_metrics_for_raw_plan_conversion)
    convertibility.checked &&
        convertibility.converted_count < length(descriptors) &&
        push!(missing, :raw_pqs_plan_conversion)
    unique!(missing)

    route_shape_mismatches = Any[
        (
            reason = :shared_pqs_descriptors_are_not_route_left_right_group,
            pqs_descriptor_count = length(descriptors),
        ),
    ]
    append!(route_shape_mismatches, direct_or_support_mismatches)
    append!(route_shape_mismatches, collect(convertibility.failures))

    diagnostics = _pqs_route_descriptor_diagnostic_common(
        source = :pqs_pqs_product_route_descriptor_diagnostic,
        route_kind = :bond_aligned_diatomic_high_order_recipe_source_construction,
        pqs_descriptor_count = length(descriptors),
        pqs_raw_plan_convertible_count = convertibility.converted_count,
        product_doside_unit_count = product_doside_unit_count,
        direct_or_support_body_piece_count =
            length(direct_or_support_mismatches),
        descriptor_emitted = false,
        supported_terms = selected_terms,
        extra = (
            status = :descriptor_unavailable,
            descriptor_available = false,
            descriptor_unavailable = true,
            metrics_supplied = !isnothing(metrics),
            raw_plan_convertibility_checked = convertibility.checked,
            raw_plan_conversion_failure_count =
                length(convertibility.failures),
            shared_shell_layer_count = length(construction.shared_shell_layers),
            region_build_count = length(construction.region_builds),
            high_order_recipe_source_construction = true,
            current_route_contains_pqs_descriptors =
                !isempty(descriptors),
            current_route_contains_explicit_product_doside_body_unit = false,
        ),
    )
    return (
        status = :descriptor_unavailable,
        descriptor = nothing,
        missing = Tuple(missing),
        mismatches = Tuple(route_shape_mismatches),
        diagnostics = diagnostics,
    )
end

function _pqs_standard_setup_charges(nuclear_charges)
    charges = nuclear_charges isa Real ?
        (Float64(nuclear_charges),) :
        Tuple(Float64(charge) for charge in nuclear_charges)
    !isempty(charges) || throw(
        ArgumentError("PQS standard source-box route setup requires at least one nuclear charge"),
    )
    all(isfinite, charges) || throw(
        ArgumentError("PQS standard source-box route setup nuclear charges must be finite"),
    )
    all(charge -> charge > 0.0, charges) || throw(
        ArgumentError("PQS standard source-box route setup nuclear charges must be positive"),
    )
    return charges
end

function _pqs_standard_setup_atom_locations(atom_locations)
    locations = Tuple(
        begin
            length(location) == 3 || throw(
                ArgumentError("PQS standard source-box route setup atom locations must be 3-vectors"),
            )
            (Float64(location[1]), Float64(location[2]), Float64(location[3]))
        end for location in atom_locations
    )
    !isempty(locations) || throw(
        ArgumentError("PQS standard source-box route setup requires at least one atom location"),
    )
    all(location -> all(isfinite, location), locations) || throw(
        ArgumentError("PQS standard source-box route setup atom locations must be finite"),
    )
    return locations
end

function _pqs_standard_setup_physical_parent_box(atom_locations, radius::Float64)
    return (
        x = (
            minimum(location -> location[1], atom_locations) - radius,
            maximum(location -> location[1], atom_locations) + radius,
        ),
        y = (
            minimum(location -> location[2], atom_locations) - radius,
            maximum(location -> location[2], atom_locations) + radius,
        ),
        z = (
            minimum(location -> location[3], atom_locations) - radius,
            maximum(location -> location[3], atom_locations) + radius,
        ),
    )
end

function _pqs_standard_n_s_core_spacing_default(n_s::Integer, q_to_core_spacing_rule::Symbol)
    q_to_core_spacing_rule == :explicit_core_spacing_only && return nothing
    _pqs_standard_validate_q_to_core_spacing_rule(q_to_core_spacing_rule)
    n_s_value = Int(n_s)
    n_s_value > 3 || throw(
        ArgumentError(
            "standard PQS n_s core-spacing default requires n_s > 3",
        ),
    )
    return (
        core_spacing = 1.2 / (4.0 * (n_s_value - 3)),
        q_to_core_spacing_rule_status = :standard_n_s_core_spacing_default,
        provenance =
            :white_lindsey_shared_shell_policy_core_spacing_1p2_over_4_ns_minus_3,
        formula = :core_spacing_equals_1p2_over_4_times_n_s_minus_3,
    )
end

function _pqs_standard_validate_q_to_core_spacing_rule(q_to_core_spacing_rule::Symbol)
    q_to_core_spacing_rule in (
        :explicit_core_spacing_only,
        :standard_pqs_ns_equals_q,
        :standard_pqs_n_s_default,
    ) || throw(
        ArgumentError(
            "unsupported PQS q-to-core-spacing rule $(q_to_core_spacing_rule)",
        ),
    )
    return nothing
end

function _pqs_standard_setup_spacing(
    nuclear_charges;
    core_spacing,
    n_s::Integer,
    reference_spacing::Float64,
    tail_spacing::Float64,
    q_to_core_spacing_rule::Symbol,
)
    _pqs_standard_validate_q_to_core_spacing_rule(q_to_core_spacing_rule)
    default_spacing = isnothing(core_spacing) ?
        _pqs_standard_n_s_core_spacing_default(n_s, q_to_core_spacing_rule) :
        nothing
    if isnothing(core_spacing) && isnothing(default_spacing)
        return (
            core_spacing = nothing,
            d = nothing,
            mapping_s = nothing,
            mapping_s_by_atom = nothing,
            core_range_by_atom = nothing,
            n_s = Int(n_s),
            reference_spacing = reference_spacing,
            tail_spacing = tail_spacing,
            q_to_core_spacing_rule = q_to_core_spacing_rule,
            q_to_core_spacing_rule_status =
                :explicit_core_spacing_required,
            provenance = :explicit_core_spacing_only_rule_selected,
            core_spacing_source = :unavailable,
            core_spacing_default_formula = nothing,
            white_lindsey_formula_available_when_d_is_explicit = true,
            non_optimality_claim = :not_claimed,
            replaceable = true,
        )
    end

    d_value = isnothing(core_spacing) ?
        Float64(default_spacing.core_spacing) :
        Float64(core_spacing)
    isfinite(d_value) && d_value > 0.0 || throw(
        ArgumentError("PQS standard source-box route setup core_spacing must be positive when provided"),
    )
    s_by_atom = Tuple(sqrt(d_value * charge) for charge in nuclear_charges)
    core_ranges = Tuple(sqrt(d_value / charge) for charge in nuclear_charges)
    same_charge = all(charge -> charge == first(nuclear_charges), nuclear_charges)
    explicit_override = !isnothing(core_spacing)
    return (
        core_spacing = d_value,
        d = d_value,
        mapping_s = same_charge ? first(s_by_atom) : nothing,
        mapping_s_by_atom = s_by_atom,
        core_range_by_atom = core_ranges,
        n_s = Int(n_s),
        reference_spacing = reference_spacing,
        tail_spacing = tail_spacing,
        q_to_core_spacing_rule = q_to_core_spacing_rule,
        q_to_core_spacing_rule_status = explicit_override ?
            :explicit_core_spacing_override :
            default_spacing.q_to_core_spacing_rule_status,
        provenance = explicit_override ?
            :explicit_core_spacing_with_white_lindsey_mapping_s_sqrt_dZ :
            default_spacing.provenance,
        core_spacing_source = explicit_override ?
            :explicit_core_spacing_override :
            :standard_n_s_default,
        core_spacing_default_formula = explicit_override ?
            nothing :
            default_spacing.formula,
        white_lindsey_formula_available_when_d_is_explicit = true,
        non_optimality_claim = :not_claimed,
        replaceable = true,
    )
end

function _pqs_standard_source_box_route_setup(;
    nuclear_charges,
    atom_locations,
    q::Integer,
    radius,
    reference_spacing::Real = 1.0,
    tail_spacing::Real = 10.0,
    q_to_core_spacing_rule::Symbol = :standard_pqs_ns_equals_q,
    core_spacing = nothing,
    n_s::Integer = q,
)
    q_value = Int(q)
    q_value >= 2 || throw(
        ArgumentError("PQS standard source-box route setup requires q >= 2"),
    )
    n_s_value = Int(n_s)
    n_s_value >= 2 || throw(
        ArgumentError("PQS standard source-box route setup requires n_s >= 2"),
    )
    radius_value = Float64(radius)
    isfinite(radius_value) && radius_value > 0.0 || throw(
        ArgumentError("PQS standard source-box route setup requires radius > 0"),
    )
    reference_spacing_value = Float64(reference_spacing)
    isfinite(reference_spacing_value) && reference_spacing_value > 0.0 || throw(
        ArgumentError("PQS standard source-box route setup requires reference_spacing > 0"),
    )
    tail_spacing_value = Float64(tail_spacing)
    isfinite(tail_spacing_value) && tail_spacing_value > 0.0 || throw(
        ArgumentError("PQS standard source-box route setup requires tail_spacing > 0"),
    )

    charges = _pqs_standard_setup_charges(nuclear_charges)
    locations = _pqs_standard_setup_atom_locations(atom_locations)
    length(charges) == length(locations) || throw(
        DimensionMismatch("PQS standard source-box route setup charge and atom-location counts must match"),
    )

    core_cube_side = isodd(q_value) ? q_value : q_value + 1
    parent_box = _pqs_standard_setup_physical_parent_box(
        locations,
        radius_value,
    )
    parent_box_lengths = (
        x = parent_box.x[2] - parent_box.x[1],
        y = parent_box.y[2] - parent_box.y[1],
        z = parent_box.z[2] - parent_box.z[1],
    )
    spacing = _pqs_standard_setup_spacing(
        charges;
        core_spacing,
        n_s = n_s_value,
        reference_spacing = reference_spacing_value,
        tail_spacing = tail_spacing_value,
        q_to_core_spacing_rule,
    )
    return (
        object_kind = :pqs_standard_source_box_route_setup,
        status = :private_development_setup,
        nuclear_charges = charges,
        atom_locations = locations,
        atom_count = length(locations),
        q = q_value,
        n_s = n_s_value,
        n_s_source = n_s_value == q_value ? :q_default : :explicit_override,
        radius = radius_value,
        core_cube_side = core_cube_side,
        core_cube_side_rule =
            :q_for_odd_q_q_plus_one_for_even_q,
        parent_box = parent_box,
        parent_box_lengths = parent_box_lengths,
        parent_box_rule =
            :minimal_axis_aligned_box_enclosing_radius_pads_around_atoms,
        spacing = spacing,
        core_spacing = spacing.core_spacing,
        d = spacing.d,
        mapping_s = spacing.mapping_s,
        mapping_s_by_atom = spacing.mapping_s_by_atom,
        reference_spacing = reference_spacing_value,
        tail_spacing = tail_spacing_value,
        q_to_core_spacing_rule = q_to_core_spacing_rule,
        diagnostics = (
            source = :pqs_standard_source_box_route_setup,
            private_development_only = true,
            production_route = false,
            n_s_equals_q = n_s_value == q_value,
            n_s_source = n_s_value == q_value ? :q_default : :explicit_override,
            q = q_value,
            n_s = n_s_value,
            core_cube_side = core_cube_side,
            core_cube_side_rule =
                :q_for_odd_q_q_plus_one_for_even_q,
            physical_parent_box_minimal_radius_pad = true,
            parent_box_rule =
                :minimal_axis_aligned_box_enclosing_radius_pads_around_atoms,
            radius = radius_value,
            atom_count = length(locations),
            charges_positive = true,
            atom_locations_finite = true,
            q_to_core_spacing_rule = q_to_core_spacing_rule,
            q_to_core_spacing_rule_status =
                spacing.q_to_core_spacing_rule_status,
            q_to_core_spacing_provenance = spacing.provenance,
            core_spacing_source = spacing.core_spacing_source,
            core_spacing_default_formula = spacing.core_spacing_default_formula,
            q_to_core_spacing_non_optimality_claim =
                spacing.non_optimality_claim,
            q_to_core_spacing_replaceable = spacing.replaceable,
            explicit_core_spacing_override_used = !isnothing(core_spacing),
            standard_n_s_default_core_spacing_used =
                spacing.core_spacing_source == :standard_n_s_default,
            public_default_consumes = false,
            packet_adoption = false,
            fixed_block_routing = false,
            qwhamiltonian_consumes = false,
            hamiltonian_matrix_built = false,
            shell_projection_used = false,
            lowdin_cleanup_used = false,
            support_local_shell_row_algorithm = false,
            support_coefficient_matrix_used = false,
            retained_pqs_weights_used = false,
            retained_weight_division_allowed = false,
            repo_side_ray_id = false,
            mwg_ida_semantics_changed = false,
            ecp_terms_implemented = false,
            cr2_science_status_changed = false,
        ),
    )
end

function _pqs_source_box_route_skeleton_axis_counts(parent_axis_counts)
    if parent_axis_counts isa NamedTuple
        all(axis -> hasproperty(parent_axis_counts, axis), (:x, :y, :z)) || throw(
            ArgumentError("PQS source-box route skeleton requires x/y/z parent axis counts"),
        )
        counts = (
            x = Int(parent_axis_counts.x),
            y = Int(parent_axis_counts.y),
            z = Int(parent_axis_counts.z),
        )
    elseif parent_axis_counts isa Tuple && length(parent_axis_counts) == 3
        counts = (
            x = Int(parent_axis_counts[1]),
            y = Int(parent_axis_counts[2]),
            z = Int(parent_axis_counts[3]),
        )
    else
        throw(
            ArgumentError(
                "PQS source-box route skeleton requires a NamedTuple x/y/z or a 3-tuple of parent axis counts",
            ),
        )
    end
    all(axis -> getproperty(counts, axis) > 0, (:x, :y, :z)) || throw(
        ArgumentError("PQS source-box route skeleton parent axis counts must be positive"),
    )
    return counts
end

function _pqs_standard_axis_aligned_diatomic_geometry(atom_locations; atol::Real = 1.0e-12)
    axis_names = (:x, :y, :z)
    length(atom_locations) == 2 || return (
        status = :not_diatomic,
        atom_count = length(atom_locations),
        bond_axis = nothing,
        bond_length = nothing,
        origin_centered_on_bond_axis = false,
        transverse_coordinates_zero = false,
        existing_bond_aligned_api_geometry_ready = false,
    )

    atol_value = Float64(atol)
    left = atom_locations[1]
    right = atom_locations[2]
    deltas = ntuple(index -> right[index] - left[index], 3)
    nonzero_indices = Tuple(index for index in 1:3 if abs(deltas[index]) > atol_value)
    length(nonzero_indices) == 1 || return (
        status = isempty(nonzero_indices) ?
            :coincident_or_degenerate_diatomic :
            :not_axis_aligned_diatomic,
        atom_count = 2,
        bond_axis = nothing,
        bond_length = nothing,
        origin_centered_on_bond_axis = false,
        transverse_coordinates_zero = false,
        existing_bond_aligned_api_geometry_ready = false,
    )

    bond_index = first(nonzero_indices)
    origin_centered_on_bond_axis = abs(left[bond_index] + right[bond_index]) <= atol_value
    transverse_coordinates_zero = all(
        index -> index == bond_index ||
                 (abs(left[index]) <= atol_value && abs(right[index]) <= atol_value),
        1:3,
    )
    return (
        status = :axis_aligned_diatomic,
        atom_count = 2,
        bond_axis = axis_names[bond_index],
        bond_length = abs(deltas[bond_index]),
        origin_centered_on_bond_axis = origin_centered_on_bond_axis,
        transverse_coordinates_zero = transverse_coordinates_zero,
        existing_bond_aligned_api_geometry_ready =
            origin_centered_on_bond_axis && transverse_coordinates_zero,
    )
end

function _pqs_standard_parent_axis_extent_candidates(setup, geometry)
    geometry.existing_bond_aligned_api_geometry_ready || return (
        available = false,
        xmax_parallel = nothing,
        xmax_transverse = nothing,
        derivation = :unavailable_without_origin_centered_axis_aligned_diatomic_geometry,
    )

    bond_axis = geometry.bond_axis
    parent_box = setup.parent_box
    parallel_interval = getproperty(parent_box, bond_axis)
    transverse_axes = Tuple(axis for axis in (:x, :y, :z) if axis != bond_axis)
    transverse_extent = maximum(
        begin
            interval = getproperty(parent_box, axis)
            max(abs(interval[1]), abs(interval[2]))
        end for axis in transverse_axes
    )
    return (
        available = true,
        xmax_parallel = max(abs(parallel_interval[1]), abs(parallel_interval[2])),
        xmax_transverse = transverse_extent,
        derivation = :from_radius_parent_box_symmetric_extent_candidates,
    )
end

function _pqs_standard_parent_axis_counts_readiness(parent_axis_counts)
    isnothing(parent_axis_counts) && return (
        parent_axis_counts = nothing,
        status = :pending_helper_or_documented_rule,
        manual_fixture = false,
        derived = false,
        derivation = :unavailable,
    )
    counts = _pqs_source_box_route_skeleton_axis_counts(parent_axis_counts)
    return (
        parent_axis_counts = counts,
        status = :manual_fixture,
        manual_fixture = true,
        derived = false,
        derivation = :manual_driver_fixture_not_standard_setup_derivation,
    )
end

function _pqs_standard_parent_axis_pending_facts(
    setup,
    geometry,
    axis_counts,
    homonuclear::Bool,
)
    pending = Symbol[]
    isnothing(setup.core_spacing) && push!(
        pending,
        :explicit_core_spacing_or_documented_q_to_core_spacing_rule,
    )
    axis_counts.derived || push!(
        pending,
        axis_counts.manual_fixture ?
        :standard_parent_axis_count_rule_replacing_manual_fixture :
        :parent_axis_counts_or_documented_axis_count_rule,
    )
    geometry.existing_bond_aligned_api_geometry_ready || push!(
        pending,
        :axis_aligned_origin_centered_diatomic_geometry_or_general_parent_api,
    )
    homonuclear || push!(pending, :atom_symbol_labels_for_heteronuclear_parent_api)
    push!(pending, :reviewed_parent_axis_constructor_call)
    return Tuple(pending)
end

function _pqs_standard_parent_axis_construction_readiness(
    setup;
    parent_axis_counts = nothing,
)
    charges = setup.nuclear_charges
    homonuclear = all(charge -> charge == first(charges), charges)
    charge_family = homonuclear ? :homonuclear : :heteronuclear
    core_spacing_available = !isnothing(setup.core_spacing)
    d_available = !isnothing(setup.d)
    white_lindsey_spacing_facts_available =
        core_spacing_available &&
        d_available &&
        !isnothing(setup.spacing.mapping_s_by_atom)
    geometry = _pqs_standard_axis_aligned_diatomic_geometry(setup.atom_locations)
    extent_candidates = _pqs_standard_parent_axis_extent_candidates(setup, geometry)
    axis_counts = _pqs_standard_parent_axis_counts_readiness(parent_axis_counts)
    homonuclear_api_appears_applicable =
        homonuclear &&
        core_spacing_available &&
        geometry.existing_bond_aligned_api_geometry_ready &&
        extent_candidates.available
    heteronuclear_api_appears_applicable =
        !homonuclear &&
        core_spacing_available &&
        geometry.existing_bond_aligned_api_geometry_ready &&
        extent_candidates.available
    existing_parent_api_appears_applicable =
        homonuclear_api_appears_applicable || heteronuclear_api_appears_applicable
    pending_facts = _pqs_standard_parent_axis_pending_facts(
        setup,
        geometry,
        axis_counts,
        homonuclear,
    )
    standard_parent_axis_rule_ready =
        core_spacing_available &&
        axis_counts.derived &&
        geometry.existing_bond_aligned_api_geometry_ready
    return (
        object_kind = :pqs_standard_parent_axis_construction_readiness,
        status = standard_parent_axis_rule_ready ?
            :standard_parent_axis_rule_ready :
            :not_ready_pending_facts,
        setup_object_kind = setup.object_kind,
        core_spacing_available = core_spacing_available,
        d_available = d_available,
        white_lindsey_spacing_facts_available =
            white_lindsey_spacing_facts_available,
        core_spacing = setup.core_spacing,
        d = setup.d,
        mapping_s = setup.mapping_s,
        mapping_s_by_atom = setup.mapping_s_by_atom,
        charge_family = charge_family,
        homonuclear = homonuclear,
        heteronuclear = !homonuclear,
        geometry = geometry,
        extent_candidates = extent_candidates,
        parent_axis_counts = axis_counts.parent_axis_counts,
        parent_axis_counts_status = axis_counts.status,
        parent_axis_counts_manual_fixture = axis_counts.manual_fixture,
        parent_axis_counts_derived = axis_counts.derived,
        parent_axis_counts_derivation = axis_counts.derivation,
        existing_parent_api_candidates = (
            bond_aligned_homonuclear_qw_basis = (
                appears_applicable = homonuclear_api_appears_applicable,
                requires_core_spacing = true,
                requires_extent_inputs = true,
                requires_homonuclear_two_center_axis_aligned_geometry = true,
                constructor_builds_axis_counts_from_mapped_extents = true,
            ),
            bond_aligned_heteronuclear_qw_basis = (
                appears_applicable = heteronuclear_api_appears_applicable,
                requires_core_spacings = true,
                requires_extent_inputs = true,
                requires_atom_symbol_labels = true,
                constructor_builds_axis_counts_from_mapped_extents = true,
            ),
            cartesian_parent_gausslet_basis = (
                appears_applicable_after_qw_basis = existing_parent_api_appears_applicable,
                wraps_existing_axis_bases = true,
            ),
        ),
        existing_parent_api_appears_applicable =
            existing_parent_api_appears_applicable,
        standard_parent_axis_rule_ready = standard_parent_axis_rule_ready,
        parent_axis_construction_ready = standard_parent_axis_rule_ready,
        parent_axis_metadata_constructed = false,
        construction_decision =
            :readiness_only_no_parent_axis_construction_added,
        pending_facts = pending_facts,
        diagnostics = (
            source = :pqs_standard_parent_axis_construction_readiness,
            private_development_only = true,
            production_route = false,
            core_spacing_available = core_spacing_available,
            white_lindsey_spacing_facts_available =
                white_lindsey_spacing_facts_available,
            existing_parent_api_appears_applicable =
                existing_parent_api_appears_applicable,
            parent_axis_counts_status = axis_counts.status,
            parent_axis_counts_manual_fixture = axis_counts.manual_fixture,
            parent_axis_counts_derived = axis_counts.derived,
            standard_parent_axis_rule_ready = standard_parent_axis_rule_ready,
            parent_axis_metadata_constructed = false,
            pending_facts = pending_facts,
            public_default_consumes = false,
            packet_adoption = false,
            fixed_block_routing = false,
            qwhamiltonian_consumes = false,
            hamiltonian_matrix_built = false,
            shell_projection_used = false,
            lowdin_cleanup_used = false,
            support_local_shell_row_algorithm = false,
            support_coefficient_matrix_used = false,
            retained_pqs_weights_used = false,
            retained_weight_division_allowed = false,
            repo_side_ray_id = false,
            mwg_ida_semantics_changed = false,
            ecp_terms_implemented = false,
            cr2_science_status_changed = false,
        ),
    )
end

function _pqs_explicit_core_spacing_probe_backend_role(gausslet_backend::Symbol)
    gausslet_backend == :pgdg_localized_experimental &&
        return :default_development_pgdg_localized
    gausslet_backend == :numerical_reference &&
        return :explicit_non_default_numerical_reference
    return :unsupported_backend
end

function _pqs_explicit_core_spacing_parent_axis_probe_pending_facts(
    readiness;
    construct_axis_bundles::Bool,
    gausslet_backend::Symbol,
)
    pending = Symbol[]
    readiness.core_spacing_available || push!(pending, :explicit_core_spacing)
    readiness.homonuclear || push!(pending, :homonuclear_setup)
    readiness.geometry.existing_bond_aligned_api_geometry_ready || push!(
        pending,
        :origin_centered_axis_aligned_diatomic_geometry,
    )
    readiness.extent_candidates.available || push!(
        pending,
        :physical_extent_inputs_for_bond_aligned_qw_basis,
    )
    construct_axis_bundles || push!(pending, :probe_parent_axis_construction_flag)
    gausslet_backend in (:pgdg_localized_experimental, :numerical_reference) || push!(
        pending,
        :pgdg_localized_or_explicit_numerical_reference_backend,
    )
    return Tuple(pending)
end

function _pqs_explicit_core_spacing_parent_axis_probe(
    setup;
    expansion = nothing,
    gausslet_backend::Symbol = :pgdg_localized_experimental,
    family = :G10,
    construct_axis_bundles::Bool = true,
    carry_objects::Bool = false,
)
    readiness = _pqs_standard_parent_axis_construction_readiness(setup)
    spacing_source = setup.spacing.core_spacing_source
    explicit_spacing_probe_only =
        spacing_source == :explicit_core_spacing_override
    default_standard_rule =
        spacing_source == :standard_n_s_default
    pending_facts = _pqs_explicit_core_spacing_parent_axis_probe_pending_facts(
        readiness;
        construct_axis_bundles,
        gausslet_backend,
    )
    construction_safe = isempty(pending_facts)
    if !construction_safe
        return (
            object_kind = :pqs_explicit_core_spacing_parent_axis_probe,
            status = :not_constructed_pending_facts,
            readiness = readiness,
            basis_metadata = nothing,
            axis_bundle_metadata = (
                object_kind = nothing,
                status = :not_constructed,
                axis_lengths = nothing,
            ),
            axis_lengths = nothing,
            physical_extent_inputs = readiness.extent_candidates,
            core_spacing = setup.core_spacing,
            reference_spacing = setup.reference_spacing,
            tail_spacing = setup.tail_spacing,
            gausslet_backend = gausslet_backend,
            gausslet_backend_role =
                _pqs_explicit_core_spacing_probe_backend_role(gausslet_backend),
            expansion_source = isnothing(expansion) ? :not_used : :explicit,
            explicit_spacing_probe_only = explicit_spacing_probe_only,
            default_standard_rule = default_standard_rule,
            core_spacing_source = spacing_source,
            parent_axis_metadata_constructed = false,
            carry_objects_requested = carry_objects,
            basis_object_available = false,
            axis_bundle_object_available = false,
            basis_object_type_label = "unavailable",
            axis_bundle_object_type_label = "unavailable",
            basis_object = nothing,
            axis_bundle_object = nothing,
            pending_facts = pending_facts,
            diagnostics = (
                source = :pqs_explicit_core_spacing_parent_axis_probe,
                private_development_only = true,
                production_route = false,
                explicit_spacing_probe_only = explicit_spacing_probe_only,
                default_standard_rule = default_standard_rule,
                core_spacing_source = spacing_source,
                parent_axis_metadata_constructed = false,
                carry_objects_requested = carry_objects,
                basis_object_available = false,
                axis_bundle_object_available = false,
                gausslet_backend = gausslet_backend,
                gausslet_backend_role =
                    _pqs_explicit_core_spacing_probe_backend_role(gausslet_backend),
                pending_facts = pending_facts,
                public_default_consumes = false,
                packet_adoption = false,
                fixed_block_routing = false,
                qwhamiltonian_consumes = false,
                hamiltonian_matrix_built = false,
                shell_projection_used = false,
                lowdin_cleanup_used = false,
                support_local_shell_row_algorithm = false,
                support_coefficient_matrix_used = false,
                retained_pqs_weights_used = false,
                retained_weight_division_allowed = false,
                repo_side_ray_id = false,
                mwg_ida_semantics_changed = false,
                ecp_terms_implemented = false,
                cr2_science_status_changed = false,
            ),
        )
    end

    expansion_value = if isnothing(expansion)
        coulomb_gaussian_expansion(doacc = false)
    elseif expansion isa CoulombGaussianExpansion
        expansion
    else
        throw(ArgumentError("explicit-core-spacing parent-axis probe expansion must be a CoulombGaussianExpansion"))
    end
    expansion_source = isnothing(expansion) ?
        :default_coulomb_gaussian_expansion_doacc_false :
        :explicit
    geometry = readiness.geometry
    extents = readiness.extent_candidates
    basis = bond_aligned_homonuclear_qw_basis(
        ;
        family,
        bond_length = geometry.bond_length,
        core_spacing = setup.core_spacing,
        xmax_parallel = extents.xmax_parallel,
        xmax_transverse = extents.xmax_transverse,
        bond_axis = geometry.bond_axis,
        nuclear_charge = first(setup.nuclear_charges),
        reference_spacing = setup.reference_spacing,
        tail_spacing = setup.tail_spacing,
    )
    axis_bundles = _qwrg_bond_aligned_axis_bundles(
        basis,
        expansion_value;
        gausslet_backend,
    )
    axis_lengths = _nested_axis_lengths(axis_bundles)
    basis_axis_lengths = (
        x = length(basis.basis_x),
        y = length(basis.basis_y),
        z = length(basis.basis_z),
    )
    basis_object = carry_objects ? basis : nothing
    axis_bundle_object = carry_objects ? axis_bundles : nothing
    return (
        object_kind = :pqs_explicit_core_spacing_parent_axis_probe,
        status = :constructed_explicit_core_spacing_parent_axis_metadata,
        readiness = readiness,
        basis_metadata = (
            object_kind = :BondAlignedDiatomicQWBasis3D,
            constructor = :bond_aligned_homonuclear_qw_basis,
            family = family,
            bond_axis = basis.bond_axis,
            bond_length = geometry.bond_length,
            nuclear_charge = first(setup.nuclear_charges),
            nuclei = Tuple(basis.nuclei),
            nuclear_charges = Tuple(basis.nuclear_charges),
            target_core_spacing = basis.target_core_spacing,
            axis_lengths = basis_axis_lengths,
        ),
        axis_bundle_metadata = (
            object_kind = :_CartesianNestedAxisBundles3D,
            constructor = :_qwrg_bond_aligned_axis_bundles,
            status = :constructed,
            axis_lengths = axis_lengths,
            gausslet_backend = gausslet_backend,
            expansion_source = expansion_source,
        ),
        axis_lengths = axis_lengths,
        physical_extent_inputs = (
            xmax_parallel = extents.xmax_parallel,
            xmax_transverse = extents.xmax_transverse,
            derivation = extents.derivation,
        ),
        core_spacing = setup.core_spacing,
        reference_spacing = setup.reference_spacing,
        tail_spacing = setup.tail_spacing,
        gausslet_backend = gausslet_backend,
        gausslet_backend_role =
            _pqs_explicit_core_spacing_probe_backend_role(gausslet_backend),
        expansion_source = expansion_source,
        explicit_spacing_probe_only = explicit_spacing_probe_only,
        default_standard_rule = default_standard_rule,
        core_spacing_source = spacing_source,
        parent_axis_metadata_constructed = true,
        carry_objects_requested = carry_objects,
        basis_object_available = !isnothing(basis_object),
        axis_bundle_object_available = !isnothing(axis_bundle_object),
        basis_object_type_label = string(nameof(typeof(basis))),
        axis_bundle_object_type_label = string(nameof(typeof(axis_bundles))),
        basis_object = basis_object,
        axis_bundle_object = axis_bundle_object,
        pending_facts = (),
        diagnostics = (
            source = :pqs_explicit_core_spacing_parent_axis_probe,
            private_development_only = true,
            production_route = false,
            explicit_spacing_probe_only = explicit_spacing_probe_only,
            default_standard_rule = default_standard_rule,
            core_spacing_source = spacing_source,
            parent_axis_metadata_constructed = true,
            carry_objects_requested = carry_objects,
            basis_object_available = !isnothing(basis_object),
            axis_bundle_object_available = !isnothing(axis_bundle_object),
            gausslet_backend = gausslet_backend,
            gausslet_backend_role =
                _pqs_explicit_core_spacing_probe_backend_role(gausslet_backend),
            axis_lengths = axis_lengths,
            public_default_consumes = false,
            packet_adoption = false,
            fixed_block_routing = false,
            qwhamiltonian_consumes = false,
            hamiltonian_matrix_built = false,
            shell_projection_used = false,
            lowdin_cleanup_used = false,
            support_local_shell_row_algorithm = false,
            support_coefficient_matrix_used = false,
            retained_pqs_weights_used = false,
            retained_weight_division_allowed = false,
            repo_side_ray_id = false,
            mwg_ida_semantics_changed = false,
            ecp_terms_implemented = false,
            cr2_science_status_changed = false,
        ),
    )
end

function _pqs_source_box_route_parent_axis_counts_for_skeleton(
    setup,
    parent_axis_readiness,
    parent_axis_probe;
    manual_parent_axis_counts = nothing,
)
    probe_constructed =
        !isnothing(parent_axis_probe) &&
        hasproperty(parent_axis_probe, :parent_axis_metadata_constructed) &&
        parent_axis_probe.parent_axis_metadata_constructed
    if probe_constructed
        counts = _pqs_source_box_route_skeleton_axis_counts(parent_axis_probe.axis_lengths)
        return (
            object_kind = :pqs_source_box_route_parent_axis_counts_for_skeleton,
            status = :available,
            parent_axis_counts = counts,
            parent_axis_counts_source = :constructed_parent_axis_probe,
            parent_axis_counts_derived = true,
            parent_axis_counts_manual_fixture = false,
            parent_axis_probe_status = parent_axis_probe.status,
            parent_axis_readiness_status = parent_axis_readiness.status,
            setup_object_kind = setup.object_kind,
            q = setup.q,
            q_minimum_satisfied =
                counts.x >= setup.q && counts.y >= setup.q && counts.z >= setup.q,
            pending_facts = (),
            diagnostics = (
                source = :pqs_source_box_route_parent_axis_counts_for_skeleton,
                private_development_only = true,
                production_route = false,
                parent_axis_counts_source = :constructed_parent_axis_probe,
                parent_axis_counts_derived = true,
                parent_axis_counts_manual_fixture = false,
                public_default_consumes = false,
                packet_adoption = false,
                fixed_block_routing = false,
                qwhamiltonian_consumes = false,
                hamiltonian_matrix_built = false,
                shell_projection_used = false,
                lowdin_cleanup_used = false,
                support_local_shell_row_algorithm = false,
                support_coefficient_matrix_used = false,
                retained_pqs_weights_used = false,
                retained_weight_division_allowed = false,
                repo_side_ray_id = false,
                mwg_ida_semantics_changed = false,
                ecp_terms_implemented = false,
                cr2_science_status_changed = false,
            ),
        )
    elseif !isnothing(manual_parent_axis_counts)
        counts = _pqs_source_box_route_skeleton_axis_counts(manual_parent_axis_counts)
        return (
            object_kind = :pqs_source_box_route_parent_axis_counts_for_skeleton,
            status = :available,
            parent_axis_counts = counts,
            parent_axis_counts_source = :manual_fixture,
            parent_axis_counts_derived = false,
            parent_axis_counts_manual_fixture = true,
            parent_axis_probe_status =
                isnothing(parent_axis_probe) ? :not_requested : parent_axis_probe.status,
            parent_axis_readiness_status = parent_axis_readiness.status,
            setup_object_kind = setup.object_kind,
            q = setup.q,
            q_minimum_satisfied =
                counts.x >= setup.q && counts.y >= setup.q && counts.z >= setup.q,
            pending_facts = (),
            diagnostics = (
                source = :pqs_source_box_route_parent_axis_counts_for_skeleton,
                private_development_only = true,
                production_route = false,
                parent_axis_counts_source = :manual_fixture,
                parent_axis_counts_derived = false,
                parent_axis_counts_manual_fixture = true,
                public_default_consumes = false,
                packet_adoption = false,
                fixed_block_routing = false,
                qwhamiltonian_consumes = false,
                hamiltonian_matrix_built = false,
                shell_projection_used = false,
                lowdin_cleanup_used = false,
                support_local_shell_row_algorithm = false,
                support_coefficient_matrix_used = false,
                retained_pqs_weights_used = false,
                retained_weight_division_allowed = false,
                repo_side_ray_id = false,
                mwg_ida_semantics_changed = false,
                ecp_terms_implemented = false,
                cr2_science_status_changed = false,
            ),
        )
    end

    probe_pending = isnothing(parent_axis_probe) ? () : parent_axis_probe.pending_facts
    pending_facts = (
        :manual_parent_axis_counts_or_constructed_parent_axis_probe,
        probe_pending...,
    )
    return (
        object_kind = :pqs_source_box_route_parent_axis_counts_for_skeleton,
        status = :not_available_pending_facts,
        parent_axis_counts = nothing,
        parent_axis_counts_source = :unavailable,
        parent_axis_counts_derived = false,
        parent_axis_counts_manual_fixture = false,
        parent_axis_probe_status =
            isnothing(parent_axis_probe) ? :not_requested : parent_axis_probe.status,
        parent_axis_readiness_status = parent_axis_readiness.status,
        setup_object_kind = setup.object_kind,
        q = setup.q,
        q_minimum_satisfied = false,
        pending_facts = pending_facts,
        diagnostics = (
            source = :pqs_source_box_route_parent_axis_counts_for_skeleton,
            private_development_only = true,
            production_route = false,
            parent_axis_counts_source = :unavailable,
            parent_axis_counts_derived = false,
            parent_axis_counts_manual_fixture = false,
            pending_facts = pending_facts,
            public_default_consumes = false,
            packet_adoption = false,
            fixed_block_routing = false,
            qwhamiltonian_consumes = false,
            hamiltonian_matrix_built = false,
            shell_projection_used = false,
            lowdin_cleanup_used = false,
            support_local_shell_row_algorithm = false,
            support_coefficient_matrix_used = false,
            retained_pqs_weights_used = false,
            retained_weight_division_allowed = false,
            repo_side_ray_id = false,
            mwg_ida_semantics_changed = false,
            ecp_terms_implemented = false,
            cr2_science_status_changed = false,
        ),
    )
end

function _pqs_source_box_route_named_box_tuple(box)
    return (box.x, box.y, box.z)
end

function _pqs_raw_product_box_plan_probe_metadata(unit_key::Symbol, source_box, plan)
    return (
        unit_key = unit_key,
        source_box = source_box,
        source_box_tuple = plan.source_box,
        source_dimensions = Tuple(length(interval) for interval in plan.source_box),
        source_mode_dims = plan.source_mode_dims,
        source_mode_count = plan.source_mode_count,
        integration_contract = plan.diagnostics.integration_contract,
        integration_contract_label = plan.diagnostics.integration_contract_label,
        numerical_reference_fallback = plan.diagnostics.numerical_reference_fallback,
        max_axis_overlap_error = plan.diagnostics.max_axis_overlap_error,
        source_product_modes_orthogonal =
            plan.diagnostics.source_product_modes_orthogonal,
        retained_rule_attached = plan.diagnostics.retained_rule_attached,
        packet_adoption = plan.diagnostics.packet_adoption,
    )
end

function _pqs_explicit_core_spacing_route_raw_product_box_plan_probe(
    setup,
    route_skeleton;
    expansion = nothing,
    gausslet_backend::Symbol = :pgdg_localized_experimental,
    family = :G10,
    construct_raw_product_box_plans::Bool = true,
)
    construct_raw_product_box_plans || return (
        object_kind = :pqs_explicit_core_spacing_route_raw_product_box_plan_probe,
        status = :not_constructed_pending_facts,
        parent_axis_probe = nothing,
        unit_plan_metadata = (),
        raw_product_box_plan_count = 0,
        all_pgdg_exact = false,
        any_numerical_reference_fallback = false,
        max_axis_overlap_error = nothing,
        gausslet_backend = gausslet_backend,
        gausslet_backend_role =
            _pqs_explicit_core_spacing_probe_backend_role(gausslet_backend),
        pending_facts = (:probe_raw_product_box_plans_flag,),
        diagnostics = (
            source = :pqs_explicit_core_spacing_route_raw_product_box_plan_probe,
            private_development_only = true,
            production_route = false,
            raw_product_box_plans_constructed = false,
            pending_facts = (:probe_raw_product_box_plans_flag,),
            gausslet_backend = gausslet_backend,
            gausslet_backend_role =
                _pqs_explicit_core_spacing_probe_backend_role(gausslet_backend),
            public_default_consumes = false,
            packet_adoption = false,
            fixed_block_routing = false,
            qwhamiltonian_consumes = false,
            hamiltonian_matrix_built = false,
            shell_projection_used = false,
            lowdin_cleanup_used = false,
            support_local_shell_row_algorithm = false,
            support_coefficient_matrix_used = false,
            retained_pqs_weights_used = false,
            retained_weight_division_allowed = false,
            repo_side_ray_id = false,
            mwg_ida_semantics_changed = false,
            ecp_terms_implemented = false,
            cr2_science_status_changed = false,
        ),
    )

    parent_axis_probe = _pqs_explicit_core_spacing_parent_axis_probe(
        setup;
        expansion,
        gausslet_backend,
        family,
        construct_axis_bundles = true,
    )
    pending = Symbol[]
    parent_axis_probe.parent_axis_metadata_constructed || append!(
        pending,
        parent_axis_probe.pending_facts,
    )
    if parent_axis_probe.parent_axis_metadata_constructed
        expected_counts =
            _pqs_source_box_route_skeleton_axis_counts(parent_axis_probe.axis_lengths)
        route_counts =
            _pqs_source_box_route_skeleton_axis_counts(route_skeleton.parent_axis_counts)
        route_counts == expected_counts || push!(
            pending,
            :route_skeleton_parent_axis_counts_from_constructed_probe,
        )
    end
    if !isempty(pending)
        pending_facts = Tuple(pending)
        return (
            object_kind = :pqs_explicit_core_spacing_route_raw_product_box_plan_probe,
            status = :not_constructed_pending_facts,
            parent_axis_probe = parent_axis_probe,
            unit_plan_metadata = (),
            raw_product_box_plan_count = 0,
            all_pgdg_exact = false,
            any_numerical_reference_fallback = false,
            max_axis_overlap_error = nothing,
            gausslet_backend = gausslet_backend,
            gausslet_backend_role =
                _pqs_explicit_core_spacing_probe_backend_role(gausslet_backend),
            pending_facts = pending_facts,
            diagnostics = (
                source = :pqs_explicit_core_spacing_route_raw_product_box_plan_probe,
                private_development_only = true,
                production_route = false,
                raw_product_box_plans_constructed = false,
                pending_facts = pending_facts,
                gausslet_backend = gausslet_backend,
                gausslet_backend_role =
                    _pqs_explicit_core_spacing_probe_backend_role(gausslet_backend),
                public_default_consumes = false,
                packet_adoption = false,
                fixed_block_routing = false,
                qwhamiltonian_consumes = false,
                hamiltonian_matrix_built = false,
                shell_projection_used = false,
                lowdin_cleanup_used = false,
                support_local_shell_row_algorithm = false,
                support_coefficient_matrix_used = false,
                retained_pqs_weights_used = false,
                retained_weight_division_allowed = false,
                repo_side_ray_id = false,
                mwg_ida_semantics_changed = false,
                ecp_terms_implemented = false,
                cr2_science_status_changed = false,
            ),
        )
    end

    expansion_value = if isnothing(expansion)
        coulomb_gaussian_expansion(doacc = false)
    elseif expansion isa CoulombGaussianExpansion
        expansion
    else
        throw(ArgumentError("raw product-box plan probe expansion must be a CoulombGaussianExpansion"))
    end
    readiness = parent_axis_probe.readiness
    geometry = readiness.geometry
    extents = readiness.extent_candidates
    basis = bond_aligned_homonuclear_qw_basis(
        ;
        family,
        bond_length = geometry.bond_length,
        core_spacing = setup.core_spacing,
        xmax_parallel = extents.xmax_parallel,
        xmax_transverse = extents.xmax_transverse,
        bond_axis = geometry.bond_axis,
        nuclear_charge = first(setup.nuclear_charges),
        reference_spacing = setup.reference_spacing,
        tail_spacing = setup.tail_spacing,
    )
    axis_bundles = _qwrg_bond_aligned_axis_bundles(
        basis,
        expansion_value;
        gausslet_backend,
    )
    unit_keys = (:pqs_left, :pqs_right, :product)
    unit_plan_metadata = map(unit_keys) do unit_key
        source_box = getproperty(route_skeleton.source_boxes, unit_key)
        source_dimensions = getproperty(route_skeleton.source_dimensions, unit_key)
        plan = _cartesian_raw_product_box_plan(
            axis_bundles,
            _pqs_source_box_route_named_box_tuple(source_box),
            source_dimensions;
            enforce_symmetric_odd = false,
        )
        _pqs_raw_product_box_plan_probe_metadata(unit_key, source_box, plan)
    end
    all_pgdg_exact =
        all(metadata -> metadata.integration_contract == :pgdg_exact, unit_plan_metadata)
    any_numerical_reference_fallback =
        any(metadata -> metadata.numerical_reference_fallback, unit_plan_metadata)
    max_axis_overlap_error =
        maximum(metadata.max_axis_overlap_error for metadata in unit_plan_metadata)
    return (
        object_kind = :pqs_explicit_core_spacing_route_raw_product_box_plan_probe,
        status = :constructed_raw_product_box_plan_metadata,
        parent_axis_probe = parent_axis_probe,
        unit_plan_metadata = unit_plan_metadata,
        raw_product_box_plan_count = length(unit_plan_metadata),
        all_pgdg_exact = all_pgdg_exact,
        any_numerical_reference_fallback = any_numerical_reference_fallback,
        max_axis_overlap_error = max_axis_overlap_error,
        gausslet_backend = gausslet_backend,
        gausslet_backend_role =
            _pqs_explicit_core_spacing_probe_backend_role(gausslet_backend),
        pending_facts = (),
        diagnostics = (
            source = :pqs_explicit_core_spacing_route_raw_product_box_plan_probe,
            private_development_only = true,
            production_route = false,
            raw_product_box_plans_constructed = true,
            raw_product_box_plan_count = length(unit_plan_metadata),
            all_pgdg_exact = all_pgdg_exact,
            any_numerical_reference_fallback = any_numerical_reference_fallback,
            max_axis_overlap_error = max_axis_overlap_error,
            gausslet_backend = gausslet_backend,
            gausslet_backend_role =
                _pqs_explicit_core_spacing_probe_backend_role(gausslet_backend),
            public_default_consumes = false,
            packet_adoption = false,
            fixed_block_routing = false,
            qwhamiltonian_consumes = false,
            hamiltonian_matrix_built = false,
            shell_projection_used = false,
            lowdin_cleanup_used = false,
            support_local_shell_row_algorithm = false,
            support_coefficient_matrix_used = false,
            retained_pqs_weights_used = false,
            retained_weight_division_allowed = false,
            repo_side_ray_id = false,
            mwg_ida_semantics_changed = false,
            ecp_terms_implemented = false,
            cr2_science_status_changed = false,
        ),
    )
end

function _pqs_source_box_route_skeleton_source_dimension(box)
    return prod(length(getproperty(box, axis)) for axis in (:x, :y, :z))
end

function _pqs_source_box_route_skeleton_range(offset::Int, count::Int)
    count > 0 || throw(ArgumentError("PQS source-box route skeleton counts must be positive"))
    return (offset + 1):(offset + count)
end

function _pqs_source_box_route_skeleton_product_slab(
    parent_axis_counts,
    q::Int,
    product_body_rule,
)
    length_value = nothing
    derivation = :pending_helper
    if product_body_rule == :centered_single_z_slab
        length_value = 1
        derivation = :derived_from_centered_single_z_slab_rule
    elseif product_body_rule isa NamedTuple &&
           get(product_body_rule, :kind, nothing) == :centered_z_slab
        rule_length = get(product_body_rule, :length, :single_plane)
        if rule_length == :single_plane
            length_value = 1
            derivation = :derived_from_named_centered_z_slab_single_plane_rule
        elseif rule_length isa Integer
            length_value = Int(rule_length)
            derivation = :manual_fixture_named_rule_length
        else
            throw(
                ArgumentError(
                    "unsupported centered_z_slab length rule $(rule_length)",
                ),
            )
        end
    else
        throw(ArgumentError("unsupported product/body slab rule $(product_body_rule)"))
    end

    length_value > 0 || throw(
        ArgumentError("product/body slab length must be positive"),
    )
    parent_axis_counts.z >= length_value || throw(
        ArgumentError("product/body slab length exceeds parent z axis count"),
    )
    parent_axis_counts.x >= q && parent_axis_counts.y >= q || throw(
        ArgumentError("product/body slab transverse source dimensions require parent x/y counts >= q"),
    )
    start_z = div(parent_axis_counts.z - length_value, 2) + 1
    z_range = start_z:(start_z + length_value - 1)
    return (
        source_box = (x = 1:q, y = 1:q, z = z_range),
        source_dimensions = (q, q, length_value),
        length = length_value,
        axis = :z,
        placement = :centered,
        derivation = derivation,
    )
end

function _pqs_source_box_route_skeleton_pair_family(
    left_kind::Symbol,
    right_kind::Symbol,
)
    if left_kind == :product_doside && right_kind == :product_doside
        return :product_product
    elseif left_kind == :pqs && right_kind == :pqs
        return :pqs_pqs
    elseif left_kind == :pqs && right_kind == :product_doside
        return :pqs_product
    elseif left_kind == :product_doside && right_kind == :pqs
        return :product_pqs
    end
    throw(ArgumentError("unsupported source-box route unit pair $(left_kind), $(right_kind)"))
end

function _pqs_source_box_route_skeleton_density_density_helper(
    pair_family::Symbol,
    pair_factor_normalization::Symbol,
)
    if pair_family == :pqs_pqs
        return pair_factor_normalization == :raw_weighted ?
            :_pqs_pqs_source_box_raw_weighted_density_density_interaction_block :
            :_pqs_pqs_source_box_density_density_interaction_block
    elseif pair_family in (:pqs_product, :product_pqs)
        return pair_factor_normalization == :raw_weighted ?
            :_pqs_product_source_box_raw_weighted_density_density_interaction_block :
            :_pqs_product_source_box_density_density_interaction_block
    elseif pair_family == :product_product
        return pair_factor_normalization == :raw_weighted ?
            :_product_doside_source_box_raw_weighted_density_density_interaction_block :
            :_product_doside_source_box_density_density_interaction_block
    end
    throw(ArgumentError("unsupported density-density pair family $(pair_family)"))
end

function _pqs_pqs_product_source_box_route_skeleton(;
    q::Integer,
    parent_axis_counts,
    route_shape = (:pqs_left, :product, :pqs_right),
    retained_unit_order = (:pqs_left, :pqs_right, :product),
    product_body_rule = :centered_single_z_slab,
    pqs_retained_rule = :boundary_comx_product_mode_selection,
    product_retained_rule = :product_doside_retained_unit,
    pair_factor_normalization::Symbol = :density_normalized,
)
    q_value = Int(q)
    q_value >= 2 || throw(
        ArgumentError("PQS source-box route skeleton requires q >= 2"),
    )
    pair_factor_normalization in (:density_normalized, :raw_weighted) || throw(
        ArgumentError("PQS source-box route skeleton requires density_normalized or raw_weighted pair-factor mode"),
    )
    route_shape == (:pqs_left, :product, :pqs_right) || throw(
        ArgumentError("PQS source-box route skeleton currently supports route shape (:pqs_left, :product, :pqs_right)"),
    )
    retained_unit_order == (:pqs_left, :pqs_right, :product) || throw(
        ArgumentError("PQS source-box route skeleton currently supports retained unit order (:pqs_left, :pqs_right, :product)"),
    )
    pqs_retained_rule == :boundary_comx_product_mode_selection || throw(
        ArgumentError("PQS source-box route skeleton currently derives only boundary COMX-product retained counts"),
    )
    product_retained_rule == :product_doside_retained_unit || throw(
        ArgumentError("PQS source-box route skeleton currently derives only product/doside retained counts"),
    )

    counts = _pqs_source_box_route_skeleton_axis_counts(parent_axis_counts)
    counts.x >= q_value && counts.y >= q_value && counts.z >= q_value || throw(
        ArgumentError("PQS source-box route skeleton parent axis counts must be >= q"),
    )
    pqs_source_mode_dims = (q_value, q_value, q_value)
    pqs_boundary_selector =
        _pqs_raw_product_box_boundary_selector(pqs_source_mode_dims)
    pqs_count = pqs_boundary_selector.selected_count
    product_slab = _pqs_source_box_route_skeleton_product_slab(
        counts,
        q_value,
        product_body_rule,
    )
    product_count = prod(product_slab.source_dimensions)

    source_boxes = (
        pqs_left = (x = 1:q_value, y = 1:q_value, z = 1:q_value),
        product = product_slab.source_box,
        pqs_right = (
            x = 1:q_value,
            y = 1:q_value,
            z = (counts.z - q_value + 1):counts.z,
        ),
    )
    source_dimensions = (
        pqs_left = pqs_source_mode_dims,
        pqs_right = pqs_source_mode_dims,
        product = product_slab.source_dimensions,
    )
    retained_counts = (
        pqs_left = pqs_count,
        pqs_right = pqs_count,
        product = product_count,
    )
    ranges = (
        pqs_left = _pqs_source_box_route_skeleton_range(0, retained_counts.pqs_left),
        pqs_right = _pqs_source_box_route_skeleton_range(
            retained_counts.pqs_left,
            retained_counts.pqs_right,
        ),
        product = _pqs_source_box_route_skeleton_range(
            retained_counts.pqs_left + retained_counts.pqs_right,
            retained_counts.product,
        ),
    )
    retained_dimension = last(ranges.product)

    retained_units = (
        (
            unit_key = :pqs_left,
            unit_role = :left_mode_selected_raw_box_pqs_unit,
            retained_unit_kind = :pqs,
            source_family = :mode_selected_raw_product_box,
            source_box = source_boxes.pqs_left,
            source_dimensions = source_dimensions.pqs_left,
            source_dimension =
                _pqs_source_box_route_skeleton_source_dimension(source_boxes.pqs_left),
            retained_rule_kind = pqs_retained_rule,
            retained_rule_derivation =
                :derived_from_boundary_comx_product_mode_selector,
            retained_range = ranges.pqs_left,
            retained_count = retained_counts.pqs_left,
            provenance_label = :pqs_left_source_modes,
            weight_semantics = :retained_columns_not_positive_quadrature_weights,
        ),
        (
            unit_key = :pqs_right,
            unit_role = :right_mode_selected_raw_box_pqs_unit,
            retained_unit_kind = :pqs,
            source_family = :mode_selected_raw_product_box,
            source_box = source_boxes.pqs_right,
            source_dimensions = source_dimensions.pqs_right,
            source_dimension =
                _pqs_source_box_route_skeleton_source_dimension(source_boxes.pqs_right),
            retained_rule_kind = pqs_retained_rule,
            retained_rule_derivation =
                :derived_from_boundary_comx_product_mode_selector,
            retained_range = ranges.pqs_right,
            retained_count = retained_counts.pqs_right,
            provenance_label = :pqs_right_source_modes,
            weight_semantics = :retained_columns_not_positive_quadrature_weights,
        ),
        (
            unit_key = :product,
            unit_role = :middle_product_doside_slab_unit,
            retained_unit_kind = :product_doside,
            source_family = :product_doside,
            source_box = source_boxes.product,
            source_dimensions = source_dimensions.product,
            source_dimension =
                _pqs_source_box_route_skeleton_source_dimension(source_boxes.product),
            retained_rule_kind = product_retained_rule,
            retained_rule_derivation =
                :derived_from_product_source_dimension,
            retained_range = ranges.product,
            retained_count = retained_counts.product,
            provenance_label = :product_doside_source_modes,
            weight_semantics = :product_source_weights_owned_by_source_box_helpers,
        ),
    )
    unit_by_key = Dict(unit.unit_key => unit for unit in retained_units)
    pair_entries = Any[]
    for i in eachindex(retained_unit_order), j in i:length(retained_unit_order)
        left_key = retained_unit_order[i]
        right_key = retained_unit_order[j]
        left_unit = unit_by_key[left_key]
        right_unit = unit_by_key[right_key]
        pair_family = _pqs_source_box_route_skeleton_pair_family(
            left_unit.retained_unit_kind,
            right_unit.retained_unit_kind,
        )
        push!(
            pair_entries,
            (
                pair_key = (left_key, right_key),
                pair_family = pair_family,
                pair_kind =
                    pair_family == :pqs_pqs ?
                    :pqs_pqs_source_box_density_density_pair :
                    pair_family == :pqs_product ?
                    :pqs_product_source_box_density_density_pair :
                    :product_doside_source_box_density_density_pair,
                density_density_helper =
                    _pqs_source_box_route_skeleton_density_density_helper(
                        pair_family,
                        pair_factor_normalization,
                    ),
                source_box_algorithmic_path = true,
                fallback_oracle_path = false,
                transpose_policy = left_key == right_key ? :none :
                    :lower_block_uses_transpose_when_pair_factors_are_symmetric,
                output_representation = :retained_two_index_density_density,
            ),
        )
    end
    pair_entries = Tuple(pair_entries)
    pair_family_counts = (
        pqs_pqs = count(entry -> entry.pair_family == :pqs_pqs, pair_entries),
        pqs_product = count(entry -> entry.pair_family == :pqs_product, pair_entries),
        product_pqs = count(entry -> entry.pair_family == :product_pqs, pair_entries),
        product_product =
            count(entry -> entry.pair_family == :product_product, pair_entries),
    )
    helper_by_pair_family = (
        pqs_pqs = _pqs_source_box_route_skeleton_density_density_helper(
            :pqs_pqs,
            pair_factor_normalization,
        ),
        pqs_product = _pqs_source_box_route_skeleton_density_density_helper(
            :pqs_product,
            pair_factor_normalization,
        ),
        product_pqs = :transpose_of_pqs_product_helper_for_lower_blocks_only,
        product_product = _pqs_source_box_route_skeleton_density_density_helper(
            :product_product,
            pair_factor_normalization,
        ),
    )
    pending_facts = (
        :parent_axis_transform_objects,
        :raw_product_box_plan_objects,
        :operator_pair_factor_data,
        :materialized_retained_operator_blocks,
    )
    return (
        object_kind = :pqs_pqs_product_source_box_route_skeleton,
        status = :private_development_skeleton,
        route_shape = route_shape,
        retained_unit_order = retained_unit_order,
        q = q_value,
        parent_axis_counts = counts,
        source_boxes = source_boxes,
        source_dimensions = source_dimensions,
        retained_units = retained_units,
        retained_counts = retained_counts,
        ranges = ranges,
        retained_dimension = retained_dimension,
        pair_entries = pair_entries,
        pair_family_counts = pair_family_counts,
        helper_by_pair_family = helper_by_pair_family,
        product_body = product_slab,
        pqs_boundary_selector = pqs_boundary_selector,
        pending_facts = pending_facts,
        diagnostics = (
            source = :pqs_pqs_product_source_box_route_skeleton,
            private_development_only = true,
            production_route = false,
            route_shape = route_shape,
            retained_unit_order = retained_unit_order,
            pqs_retained_count_derivation =
                :boundary_comx_product_mode_selection,
            pqs_boundary_selection_rule = pqs_boundary_selector.selection_rule,
            pqs_boundary_column_count = pqs_count,
            product_retained_count_derivation =
                :product_source_dimension_under_product_doside_rule,
            product_body_rule = product_body_rule,
            product_slab_length = product_slab.length,
            product_slab_length_derivation = product_slab.derivation,
            product_slab_placement = product_slab.placement,
            derived_retained_counts = true,
            derived_retained_ranges = true,
            derived_pair_inventory = true,
            pending_facts = pending_facts,
            pair_factor_normalization = pair_factor_normalization,
            raw_weight_division_owner =
                pair_factor_normalization == :raw_weighted ?
                :explicit_source_quadrature_weight_outer_products :
                :caller_supplied_density_normalized_pair_factors,
            source_box_first = true,
            source_box_algorithmic_path_true_for_every_pair =
                all(entry -> entry.source_box_algorithmic_path, pair_entries),
            pair_count = length(pair_entries),
            pair_family_counts = pair_family_counts,
            retained_dimension = retained_dimension,
            retained_unit_count = length(retained_units),
            output_representation = :retained_two_index_density_density,
            four_index_galerkin_tensor = false,
            product_pqs_explicit_helper_required = false,
            product_pqs_transpose_requires_symmetric_pair_factors = true,
            product_pqs_lower_block_count = 2,
            retained_pqs_weights_used = false,
            retained_weight_division_allowed = false,
            packet_adoption = false,
            fixed_block_routing = false,
            qwhamiltonian_consumes = false,
            public_default_consumes = false,
            mwg_ida_semantics_changed = false,
            repo_side_ray_id = false,
            ecp_terms_implemented = false,
            cr2_science_status_changed = false,
            shell_projection_used = false,
            lowdin_cleanup_used = false,
            support_local_shell_row_algorithm = false,
            support_coefficient_matrix_used = false,
        ),
    )
end

function _pqs_parent_coefficient_matrix_from_raw_plan(raw_plan, parent_dims::NTuple{3,Int})
    plan = _pqs_raw_product_box_plan_view(raw_plan)
    plan.representation == :orthogonal_raw_product_box || throw(
        ArgumentError("PQS parent coefficient projection requires a raw product-box plan"),
    )
    intervals = plan.axis_intervals
    axis_coefficients =
        ntuple(axis -> Matrix{Float64}(plan.axis_local_coefficients[axis]), 3)
    for axis in 1:3
        interval = intervals[axis]
        first(interval) >= 1 && last(interval) <= parent_dims[axis] || throw(
            ArgumentError("PQS parent coefficient projection interval exceeds parent dimensions"),
        )
        size(axis_coefficients[axis], 1) == length(interval) || throw(
            DimensionMismatch("PQS parent coefficient projection axis rows must match source interval length"),
        )
        size(axis_coefficients[axis], 2) == plan.source_mode_dims[axis] || throw(
            DimensionMismatch("PQS parent coefficient projection axis columns must match source modes"),
        )
    end
    modes = plan.boundary_selector.mode_indices
    retained_count = plan.boundary_selector.selected_count
    length(modes) == retained_count || throw(
        DimensionMismatch("PQS parent coefficient projection boundary mode count must match retained count"),
    )
    coefficients = zeros(Float64, prod(parent_dims), retained_count)
    x_interval, y_interval, z_interval = intervals
    cx, cy, cz = axis_coefficients
    @inbounds for (column, mode) in pairs(modes)
        mx, my, mz = mode
        for (local_x, parent_x) in enumerate(x_interval)
            x_value = cx[local_x, mx]
            iszero(x_value) && continue
            for (local_y, parent_y) in enumerate(y_interval)
                xy_value = x_value * cy[local_y, my]
                iszero(xy_value) && continue
                for (local_z, parent_z) in enumerate(z_interval)
                    value = xy_value * cz[local_z, mz]
                    iszero(value) && continue
                    row = _cartesian_flat_index(
                        parent_x,
                        parent_y,
                        parent_z,
                        parent_dims,
                    )
                    coefficients[row, column] += value
                end
            end
        end
    end
    return coefficients
end

function _source_box_nuclear_attraction_center_labels(center_count::Int, center_labels)
    center_count > 0 || throw(
        ArgumentError("source-box nuclear attraction center labels require a positive center count"),
    )
    if isnothing(center_labels)
        return ntuple(index -> Symbol("center_", index), center_count)
    end
    labels = Tuple(Symbol(label) for label in center_labels)
    length(labels) == center_count || throw(
        ArgumentError("source-box nuclear attraction center label count must match center count"),
    )
    length(unique(labels)) == center_count || throw(
        ArgumentError("source-box nuclear attraction center labels must be unique"),
    )
    return labels
end

function _pqs_pqs_product_nuclear_attraction_all_pairs_inventory(
    left_raw_plan,
    right_raw_plan,
    product_unit::_CartesianNestedProductStagedByCenterUnit3D,
    ranges,
)
    product_retained_unit_plan = _product_doside_retained_unit_plan(product_unit)
    retained_units = (
        (
            unit_key = :pqs_left,
            retained_unit_kind = :pqs,
            source_family = :mode_selected_raw_product_box,
            retained_rule_kind = left_raw_plan.retained_rule_kind,
            retained_range = ranges.pqs_left,
            source_dimensions = left_raw_plan.source_mode_dims,
            source_dimension = left_raw_plan.source_mode_count,
            retained_count = left_raw_plan.boundary_selector.selected_count,
            supported_one_body_terms = (:electron_nuclear_attraction,),
        ),
        (
            unit_key = :pqs_right,
            retained_unit_kind = :pqs,
            source_family = :mode_selected_raw_product_box,
            retained_rule_kind = right_raw_plan.retained_rule_kind,
            retained_range = ranges.pqs_right,
            source_dimensions = right_raw_plan.source_mode_dims,
            source_dimension = right_raw_plan.source_mode_count,
            retained_count = right_raw_plan.boundary_selector.selected_count,
            supported_one_body_terms = (:electron_nuclear_attraction,),
        ),
        (
            unit_key = :product,
            retained_unit_kind = :product_doside,
            source_family = :product_doside,
            retained_rule_kind = product_retained_unit_plan.retained_rule_kind,
            retained_range = ranges.product,
            source_dimensions = product_retained_unit_plan.source_axis_lengths,
            source_dimension = product_retained_unit_plan.source_dimension,
            retained_count = product_retained_unit_plan.retained_count,
            supported_one_body_terms = (:electron_nuclear_attraction,),
        ),
    )
    pqs_pqs_helper = :_pqs_pqs_source_box_nuclear_attraction_by_center
    pqs_product_helper = :_pqs_product_source_box_nuclear_attraction_by_center
    product_product_helper = :_product_doside_source_box_nuclear_attraction_by_center
    pair_entries = (
        (
            pair_key = (:pqs_left, :pqs_left),
            pair_family = :pqs_pqs,
            pair_kind = :pqs_pqs_source_box_nuclear_attraction_pair,
            block_helper = pqs_pqs_helper,
            upper_triangular = true,
            transpose_only_lower_block = false,
            pair_policy = :source_box_algorithm_available,
            source_box_algorithmic = true,
        ),
        (
            pair_key = (:pqs_left, :pqs_right),
            pair_family = :pqs_pqs,
            pair_kind = :pqs_pqs_source_box_nuclear_attraction_pair,
            block_helper = pqs_pqs_helper,
            upper_triangular = true,
            transpose_only_lower_block = true,
            pair_policy = :source_box_algorithm_available,
            source_box_algorithmic = true,
        ),
        (
            pair_key = (:pqs_left, :product),
            pair_family = :pqs_product,
            pair_kind = :pqs_product_source_box_nuclear_attraction_pair,
            block_helper = pqs_product_helper,
            upper_triangular = true,
            transpose_only_lower_block = true,
            pair_policy = :source_box_algorithm_available,
            source_box_algorithmic = true,
        ),
        (
            pair_key = (:pqs_right, :pqs_right),
            pair_family = :pqs_pqs,
            pair_kind = :pqs_pqs_source_box_nuclear_attraction_pair,
            block_helper = pqs_pqs_helper,
            upper_triangular = true,
            transpose_only_lower_block = false,
            pair_policy = :source_box_algorithm_available,
            source_box_algorithmic = true,
        ),
        (
            pair_key = (:pqs_right, :product),
            pair_family = :pqs_product,
            pair_kind = :pqs_product_source_box_nuclear_attraction_pair,
            block_helper = pqs_product_helper,
            upper_triangular = true,
            transpose_only_lower_block = true,
            pair_policy = :source_box_algorithm_available,
            source_box_algorithmic = true,
        ),
        (
            pair_key = (:product, :product),
            pair_family = :product_product,
            pair_kind = :product_doside_source_box_nuclear_attraction_pair,
            block_helper = product_product_helper,
            upper_triangular = true,
            transpose_only_lower_block = false,
            pair_policy = :source_box_algorithm_available,
            source_box_algorithmic = true,
        ),
    )
    pair_family_counts = (
        pqs_pqs = count(entry -> entry.pair_family == :pqs_pqs, pair_entries),
        pqs_product = count(entry -> entry.pair_family == :pqs_product, pair_entries),
        product_product =
            count(entry -> entry.pair_family == :product_product, pair_entries),
    )
    return (
        object_kind = :pqs_pqs_product_nuclear_attraction_all_pairs_inventory,
        retained_units = retained_units,
        pair_entries = pair_entries,
        pair_family_counts = pair_family_counts,
        supported_one_body_terms = (:electron_nuclear_attraction,),
        diagnostics = (
            source = :pqs_pqs_product_nuclear_attraction_all_pairs_inventory,
            all_pairs_inventory_private = true,
            pair_inventory_complete_for_units = (:pqs_left, :pqs_right, :product),
            retained_unit_count = length(retained_units),
            upper_triangular_pair_count = length(pair_entries),
            expected_upper_triangular_pair_count = 6,
            pair_family_counts = pair_family_counts,
            pair_keys = map(entry -> entry.pair_key, pair_entries),
            pair_families = map(entry -> entry.pair_family, pair_entries),
            block_helpers = map(entry -> entry.block_helper, pair_entries),
            block_helper_by_family = (
                pqs_pqs = pqs_pqs_helper,
                pqs_product = pqs_product_helper,
                product_product = product_product_helper,
            ),
            pair_policies = map(entry -> entry.pair_policy, pair_entries),
            every_pair_uses_source_box_algorithmic_policy = all(
                entry -> entry.source_box_algorithmic,
                pair_entries,
            ),
            source_box_algorithmic_pair_count =
                count(entry -> entry.source_box_algorithmic, pair_entries),
            private_shadow_only = true,
            physical_operator = :electron_nuclear_attraction,
            positive_gaussian_sum_component = true,
            center_contributions_preserved = true,
            counterpoise_center_identity_preserved = true,
            lower_triangular_cross_blocks_transpose_only = true,
            product_pqs_blocks_transpose_only = true,
            product_pqs_explicit_helper_used = false,
            packet_adoption = false,
            fixed_block_routing = false,
            qwhamiltonian_consumes = false,
            public_default_consumes = false,
            cr2_science_status_changed = false,
            ecp_terms_implemented = false,
            electron_electron_terms_implemented = false,
            mwg_interaction_implemented = false,
            ida_mwg_semantics_changed = false,
            mwg_ida_semantics_changed = false,
            shell_projection_used = false,
            lowdin_cleanup_used = false,
            support_local_oracle_used = false,
            support_local_pqs_oracle_used = false,
            support_coefficient_matrix_used = false,
            shell_row_algorithm = false,
            retained_pqs_weights_used = false,
            retained_pqs_weights_positive_checked = false,
            retained_weight_semantics = :not_positive_quadrature_weights,
            retained_weight_division_allowed = false,
            retained_pqs_weight_division_allowed = false,
            ida_weight_division_allowed = false,
            dense_raw_source_box_pair_matrix_materialized = false,
            dense_raw_pair_storage_avoided = true,
        ),
    )
end

function _pqs_pqs_product_nuclear_attraction_pair_block(
    entry,
    units;
    axis_layers::NamedTuple{(:x,:y,:z)},
    expansion::CoulombGaussianExpansion,
    centers,
    nuclear_charges,
)
    left_role, right_role = entry.pair_key
    if entry.pair_family == :pqs_pqs
        left_unit = getproperty(units, left_role)
        right_unit = getproperty(units, right_role)
        return _pqs_pqs_source_box_nuclear_attraction_by_center(
            left_unit,
            right_unit,
            axis_layers,
            expansion;
            centers,
            nuclear_charges,
        )
    elseif entry.pair_family == :pqs_product
        right_role == :product || throw(
            ArgumentError("nuclear-attraction route only uses product/PQS as transpose of PQS/product"),
        )
        pqs_unit = getproperty(units, left_role)
        product_unit = getproperty(units, right_role)
        return _pqs_product_source_box_nuclear_attraction_by_center(
            pqs_unit,
            product_unit,
            axis_layers,
            expansion;
            centers,
            nuclear_charges,
        )
    elseif entry.pair_family == :product_product
        left_role == :product && right_role == :product || throw(
            ArgumentError("nuclear-attraction route product/product entry must use the product unit on both sides"),
        )
        product_unit = getproperty(units, :product)
        return _product_doside_source_box_nuclear_attraction_by_center(
            product_unit,
            product_unit,
            axis_layers,
            expansion;
            centers,
            nuclear_charges,
        )
    end
    throw(ArgumentError("unsupported nuclear-attraction route pair family $(entry.pair_family)"))
end

function _pqs_pqs_product_route_shaped_nuclear_attraction_by_center(
    route_units,
    axis_layers::NamedTuple{(:x,:y,:z)},
    expansion::CoulombGaussianExpansion;
    centers,
    nuclear_charges,
    center_labels = nothing,
    symmetry_atol::Real = 1.0e-10,
)
    hasproperty(route_units, :route_kind) || throw(
        ArgumentError("route-shaped nuclear-attraction consumer requires route_kind"),
    )
    hasproperty(route_units, :units) || throw(
        ArgumentError("route-shaped nuclear-attraction consumer requires units"),
    )
    route_kind = route_units.route_kind
    route_kind in _PQS_PQS_PRODUCT_SAFE_TERM_ROUTE_KINDS || throw(
        ArgumentError("unsupported route-shaped nuclear-attraction consumer route_kind $(route_kind)"),
    )
    roles = hasproperty(route_units, :roles) ?
        Tuple(route_units.roles) :
        (:pqs_left, :pqs_right, :product)
    roles == (:pqs_left, :pqs_right, :product) || throw(
        ArgumentError("route-shaped nuclear-attraction consumer requires roles (:pqs_left, :pqs_right, :product)"),
    )
    units = route_units.units
    for role in roles
        hasproperty(units, role) || throw(
            ArgumentError("route-shaped nuclear-attraction consumer missing unit $(role)"),
        )
    end
    center_values = _source_box_nuclear_attraction_center_values(centers)
    charge_values = _source_box_nuclear_attraction_charge_values(nuclear_charges)
    length(center_values) == length(charge_values) || throw(
        ArgumentError("route-shaped nuclear-attraction center and charge counts must match"),
    )
    labels = _source_box_nuclear_attraction_center_labels(
        length(center_values),
        center_labels,
    )

    left_raw_plan = _pqs_raw_product_box_plan_view(units.pqs_left)
    right_raw_plan = _pqs_raw_product_box_plan_view(units.pqs_right)
    left_raw_plan.representation == :orthogonal_raw_product_box || throw(
        ArgumentError("route-shaped nuclear-attraction consumer requires a raw product-box left PQS plan"),
    )
    right_raw_plan.representation == :orthogonal_raw_product_box || throw(
        ArgumentError("route-shaped nuclear-attraction consumer requires a raw product-box right PQS plan"),
    )
    units.product.kind == :product_doside || throw(
        ArgumentError("route-shaped nuclear-attraction consumer requires a product_doside middle unit"),
    )

    product_retained_unit_plan = _product_doside_retained_unit_plan(units.product)
    left_count = left_raw_plan.boundary_selector.selected_count
    right_count = right_raw_plan.boundary_selector.selected_count
    product_count = product_retained_unit_plan.retained_count
    ranges = (
        pqs_left = 1:left_count,
        pqs_right = (left_count + 1):(left_count + right_count),
        product = (left_count + right_count + 1):
                  (left_count + right_count + product_count),
    )
    retained_dimension = left_count + right_count + product_count
    inventory = _pqs_pqs_product_nuclear_attraction_all_pairs_inventory(
        left_raw_plan,
        right_raw_plan,
        units.product,
        ranges,
    )
    inventory.diagnostics.every_pair_uses_source_box_algorithmic_policy ||
        throw(ArgumentError("nuclear-attraction route requires all pairs to use source-box algorithms"))

    descriptor_expected_ranges_checked = false
    descriptor_retained_dimension_checked = false
    descriptor_pair_count_checked = false
    if hasproperty(route_units, :expected_ranges)
        route_units.expected_ranges == ranges || throw(
            ArgumentError("route descriptor expected ranges disagree with nuclear-attraction route ranges"),
        )
        descriptor_expected_ranges_checked = true
    end
    if hasproperty(route_units, :retained_dimension)
        route_units.retained_dimension == retained_dimension || throw(
            DimensionMismatch("route descriptor retained dimension disagrees with nuclear-attraction route"),
        )
        descriptor_retained_dimension_checked = true
    end
    if hasproperty(route_units, :expected_pair_count)
        route_units.expected_pair_count == length(inventory.pair_entries) || throw(
            ArgumentError("route descriptor expected pair count disagrees with nuclear-attraction route"),
        )
        descriptor_pair_count_checked = true
    end

    timed = @timed begin
        center_blocks = [
            zeros(Float64, retained_dimension, retained_dimension) for
            _ in eachindex(center_values)
        ]
        pair_block_results = Dict{Tuple{Symbol,Symbol},Any}()
        for entry in inventory.pair_entries
            pair_result = _pqs_pqs_product_nuclear_attraction_pair_block(
                entry,
                units;
                axis_layers,
                expansion,
                centers = center_values,
                nuclear_charges = charge_values,
            )
            length(pair_result.blocks_by_center) == length(center_values) || throw(
                DimensionMismatch("nuclear-attraction route pair center count mismatch"),
            )
            left_role, right_role = entry.pair_key
            left_range = getproperty(ranges, left_role)
            right_range = getproperty(ranges, right_role)
            for center_index in eachindex(center_values)
                center_block =
                    pair_result.blocks_by_center[center_index].block
                center_blocks[center_index][left_range, right_range] .=
                    center_block
                if left_role != right_role
                    entry.transpose_only_lower_block || throw(
                        ArgumentError("nuclear-attraction route lower block requires transpose-only policy"),
                    )
                    center_blocks[center_index][right_range, left_range] .=
                        transpose(center_block)
                end
            end
            pair_block_results[entry.pair_key] = pair_result
        end
        (center_blocks = center_blocks, pair_block_results = pair_block_results)
    end
    assembled = timed.value
    center_blocks = assembled.center_blocks
    pair_block_results = assembled.pair_block_results
    total_block = zeros(Float64, retained_dimension, retained_dimension)
    for center_block in center_blocks
        all(isfinite, center_block) || throw(
            ArgumentError("route-shaped nuclear-attraction center block produced non-finite entries"),
        )
        total_block .+= center_block
    end
    all(isfinite, total_block) || throw(
        ArgumentError("route-shaped nuclear-attraction total block produced non-finite entries"),
    )
    center_symmetry_errors = ntuple(
        center_index -> LinearAlgebra.norm(
            center_blocks[center_index] - transpose(center_blocks[center_index]),
            Inf,
        ),
        length(center_values),
    )
    symmetry_error = LinearAlgebra.norm(total_block - transpose(total_block), Inf)
    symmetry_atol_value = Float64(symmetry_atol)
    all(error -> error <= symmetry_atol_value, center_symmetry_errors) || throw(
        ArgumentError("route-shaped nuclear-attraction center symmetry error exceeded tolerance"),
    )
    symmetry_error <= symmetry_atol_value || throw(
        ArgumentError("route-shaped nuclear-attraction total symmetry error exceeded tolerance"),
    )
    total_from_center_blocks = zeros(Float64, retained_dimension, retained_dimension)
    for center_block in center_blocks
        total_from_center_blocks .+= center_block
    end
    total_from_center_error =
        LinearAlgebra.norm(total_block - total_from_center_blocks, Inf)

    center_component_blocks = ntuple(center_index -> (
        pqs_left_pqs_left =
            pair_block_results[(:pqs_left, :pqs_left)].blocks_by_center[center_index].block,
        pqs_left_pqs_right =
            pair_block_results[(:pqs_left, :pqs_right)].blocks_by_center[center_index].block,
        pqs_right_pqs_left =
            transpose(pair_block_results[(:pqs_left, :pqs_right)].blocks_by_center[center_index].block),
        pqs_left_product =
            pair_block_results[(:pqs_left, :product)].blocks_by_center[center_index].block,
        product_pqs_left =
            transpose(pair_block_results[(:pqs_left, :product)].blocks_by_center[center_index].block),
        pqs_right_pqs_right =
            pair_block_results[(:pqs_right, :pqs_right)].blocks_by_center[center_index].block,
        pqs_right_product =
            pair_block_results[(:pqs_right, :product)].blocks_by_center[center_index].block,
        product_pqs_right =
            transpose(pair_block_results[(:pqs_right, :product)].blocks_by_center[center_index].block),
        product_product =
            pair_block_results[(:product, :product)].blocks_by_center[center_index].block,
    ), length(center_values))
    blocks_by_center = ntuple(center_index -> (
        center_index = center_index,
        center_label = labels[center_index],
        center = center_values[center_index],
        nuclear_charge = charge_values[center_index],
        sign_charge_scale = -charge_values[center_index],
        block = center_blocks[center_index],
        component_blocks = center_component_blocks[center_index],
        symmetry_error = center_symmetry_errors[center_index],
    ), length(center_values))
    component_block_provenance = (
        pqs_left_pqs_left =
            inventory.diagnostics.block_helper_by_family.pqs_pqs,
        pqs_left_pqs_right =
            inventory.diagnostics.block_helper_by_family.pqs_pqs,
        pqs_right_pqs_left = :transpose_of_pqs_left_pqs_right,
        pqs_left_product =
            inventory.diagnostics.block_helper_by_family.pqs_product,
        product_pqs_left = :transpose_of_pqs_left_product,
        pqs_right_pqs_right =
            inventory.diagnostics.block_helper_by_family.pqs_pqs,
        pqs_right_product =
            inventory.diagnostics.block_helper_by_family.pqs_product,
        product_pqs_right = :transpose_of_pqs_right_product,
        product_product =
            inventory.diagnostics.block_helper_by_family.product_product,
    )
    metadata = hasproperty(route_units, :metadata) ? route_units.metadata : (;)
    provenance = hasproperty(route_units, :provenance) ? route_units.provenance : (;)
    route_name = hasproperty(route_units, :route_name) ? route_units.route_name : route_kind
    route_descriptor_object_kind =
        hasproperty(route_units, :object_kind) ? route_units.object_kind : :legacy_route_units
    performance = (
        elapsed_seconds = Float64(timed.time),
        allocated_bytes = Int(timed.bytes),
        gc_time_seconds = Float64(timed.gctime),
        retained_dimension = retained_dimension,
        pair_count = length(inventory.pair_entries),
        center_count = length(center_values),
        dense_raw_source_box_pair_matrix_materialized = false,
        dense_raw_pair_storage_avoided = true,
    )
    return (
        path = :pqs_pqs_product_route_shaped_nuclear_attraction_by_center,
        route_kind = route_kind,
        route_name = route_name,
        route_units = route_units,
        physical_operator = :electron_nuclear_attraction,
        retained_units = inventory.retained_units,
        all_pairs_inventory = inventory,
        blocks_by_center = blocks_by_center,
        center_blocks = blocks_by_center,
        per_center_matrices = blocks_by_center,
        total_block = total_block,
        block = total_block,
        complete_retained_space_matrix = total_block,
        component_blocks_by_center = center_component_blocks,
        component_block_provenance = component_block_provenance,
        pair_block_results = pair_block_results,
        ranges = ranges,
        retained_dimension = retained_dimension,
        pair_count = length(inventory.pair_entries),
        pair_family_counts = inventory.pair_family_counts,
        centers = Tuple(center_values),
        center_labels = labels,
        nuclear_charges = Tuple(charge_values),
        output_finite = true,
        symmetry_error = symmetry_error,
        center_symmetry_errors = center_symmetry_errors,
        total_from_center_error = total_from_center_error,
        performance = performance,
        metadata = metadata,
        provenance = provenance,
        diagnostics = merge(
            inventory.diagnostics,
            (
                source = :pqs_pqs_product_route_shaped_nuclear_attraction_by_center,
                route_shaped_consumer = true,
                route_shape = (:pqs_left, :pqs_right, :product),
                route_kind = route_kind,
                route_name = route_name,
                route_descriptor_object_kind = route_descriptor_object_kind,
                route_roles = roles,
                retained_ranges = ranges,
                retained_dimension = retained_dimension,
                retained_unit_count = length(inventory.retained_units),
                pair_count = length(inventory.pair_entries),
                pair_family_counts = inventory.pair_family_counts,
                physical_operator = :electron_nuclear_attraction,
                one_body_operator = true,
                electron_electron_terms_implemented = false,
                helper_used_for_pair_families =
                    inventory.diagnostics.block_helper_by_family,
                block_helpers = inventory.diagnostics.block_helpers,
                source_box_first = true,
                source_box_algorithmic_path_true_for_every_pair = true,
                every_pair_uses_source_box_algorithmic_policy =
                    inventory.diagnostics.every_pair_uses_source_box_algorithmic_policy,
                source_box_algorithmic_pair_count =
                    inventory.diagnostics.source_box_algorithmic_pair_count,
                centers = Tuple(center_values),
                center_labels = labels,
                nuclear_charges = Tuple(charge_values),
                center_records = ntuple(
                    center_index -> (
                        label = labels[center_index],
                        center = center_values[center_index],
                        nuclear_charge = charge_values[center_index],
                    ),
                    length(center_values),
                ),
                center_count = length(center_values),
                center_contributions_preserved = true,
                counterpoise_center_identity_preserved = true,
                total_block_is_explicit_sum_of_center_pieces = true,
                total_from_center_error = total_from_center_error,
                lower_triangular_cross_blocks_transpose_only = true,
                product_pqs_blocks_transpose_only = true,
                output_finite = true,
                symmetry_error = symmetry_error,
                center_symmetry_errors = center_symmetry_errors,
                symmetry_atol = symmetry_atol_value,
                descriptor_expected_ranges_checked =
                    descriptor_expected_ranges_checked,
                descriptor_retained_dimension_checked =
                    descriptor_retained_dimension_checked,
                descriptor_pair_count_checked = descriptor_pair_count_checked,
                positive_gaussian_sum_component = true,
                positive_gaussian_sum_convention = true,
                nuclear_charge_applied = true,
                nuclear_attraction_sign_applied = true,
                nuclear_charge_sign_applied = true,
                local_gaussian_source_box_terms = true,
                ecp = false,
                ecp_terms_implemented = false,
                mwg_interaction_implemented = false,
                ida_mwg_semantics_changed = false,
                mwg_ida_semantics_changed = false,
                retained_pqs_weights_used = false,
                retained_pqs_weights_positive_checked = false,
                retained_weight_division_allowed = false,
                retained_pqs_weight_division_allowed = false,
                ida_weight_division_allowed = false,
                retained_weight_semantics = :not_positive_quadrature_weights,
                shell_projection_used = false,
                lowdin_cleanup_used = false,
                support_local_oracle_used = false,
                support_local_pqs_oracle_used = false,
                support_local_shell_row_algorithm = false,
                support_coefficient_matrix_used = false,
                shell_row_algorithm = false,
                packet_adoption = false,
                fixed_block_routing = false,
                qwhamiltonian_consumes = false,
                public_default_consumes = false,
                cr2_science_status_changed = false,
                dense_raw_source_box_pair_matrix_materialized = false,
                dense_raw_pair_storage_avoided = true,
                performance_recorded = true,
                elapsed_seconds = performance.elapsed_seconds,
                allocated_bytes = performance.allocated_bytes,
                gc_time_seconds = performance.gc_time_seconds,
            ),
        ),
    )
end

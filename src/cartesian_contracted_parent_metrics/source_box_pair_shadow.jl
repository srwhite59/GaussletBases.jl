function _pqs_raw_product_box_source_mode_indices(source_mode_dims::NTuple{3,Int})
    return _cartesian_raw_product_box_source_mode_indices(source_mode_dims)
end

function _pqs_raw_product_box_mode_matrix(
    mode_indices::AbstractVector{<:NTuple{3,Int}},
    axis_matrices::NTuple{3,AbstractMatrix{<:Real}},
)
    mode_count = length(mode_indices)
    block = Matrix{Float64}(undef, mode_count, mode_count)
    @inbounds for col in 1:mode_count
        jx, jy, jz = mode_indices[col]
        for row in 1:mode_count
            ix, iy, iz = mode_indices[row]
            block[row, col] =
                axis_matrices[1][ix, jx] *
                axis_matrices[2][iy, jy] *
                axis_matrices[3][iz, jz]
        end
    end
    return block
end

function _pqs_raw_product_box_pair_mode_matrix(
    left_mode_indices::AbstractVector{<:NTuple{3,Int}},
    right_mode_indices::AbstractVector{<:NTuple{3,Int}},
    axis_matrices::NTuple{3,AbstractMatrix{<:Real}},
)
    block = Matrix{Float64}(undef, length(left_mode_indices), length(right_mode_indices))
    @inbounds for col in eachindex(right_mode_indices)
        jx, jy, jz = right_mode_indices[col]
        for row in eachindex(left_mode_indices)
            ix, iy, iz = left_mode_indices[row]
            block[row, col] =
                axis_matrices[1][ix, jx] *
                axis_matrices[2][iy, jy] *
                axis_matrices[3][iz, jz]
        end
    end
    return block
end

function _pqs_raw_product_box_low_order_axis_kinds(term::Symbol)
    term == :overlap && return (:overlap, :overlap, :overlap)
    term == :position_x && return (:position, :overlap, :overlap)
    term == :position_y && return (:overlap, :position, :overlap)
    term == :position_z && return (:overlap, :overlap, :position)
    term == :x2_x && return (:x2, :overlap, :overlap)
    term == :x2_y && return (:overlap, :x2, :overlap)
    term == :x2_z && return (:overlap, :overlap, :x2)
    throw(
        ArgumentError(
            "PQS raw product-box reference block supports only :overlap, :position_x/y/z, :x2_x/y/z, and :kinetic",
        ),
    )
end

function _pqs_raw_product_box_selected_mode_matrix(
    raw_plan,
    axis_matrices::NTuple{3,AbstractMatrix{<:Real}},
)
    mode_count = length(raw_plan.boundary_selector.mode_indices)
    mode_count == raw_plan.boundary_selector.selected_count || throw(
        ArgumentError("PQS raw product-box boundary mode count must match selected count"),
    )
    return _pqs_raw_product_box_mode_matrix(
        raw_plan.boundary_selector.mode_indices,
        axis_matrices,
    )
end

function _pqs_raw_product_box_factor_tuple(raw_plan)
    raw_plan.operator_factors_available || throw(
        ArgumentError("PQS raw product-box plan does not carry 1D operator factors"),
    )
    factors = raw_plan.one_dimensional_operator_factors
    return (factors.x, factors.y, factors.z)
end

function _pqs_raw_product_box_low_order_mode_matrix(raw_plan; term::Symbol)
    axis_kinds = _pqs_raw_product_box_low_order_axis_kinds(term)
    factors = _pqs_raw_product_box_factor_tuple(raw_plan)
    axis_matrices = ntuple(axis -> getproperty(factors[axis], axis_kinds[axis]), 3)
    return _pqs_raw_product_box_selected_mode_matrix(raw_plan, axis_matrices)
end

function _pqs_raw_product_box_kinetic_mode_matrix(raw_plan)
    factors = _pqs_raw_product_box_factor_tuple(raw_plan)
    overlap_matrices = ntuple(axis -> factors[axis].overlap, 3)
    kinetic_matrices = ntuple(axis -> factors[axis].kinetic, 3)
    mode_matrix = zeros(
        Float64,
        raw_plan.boundary_selector.selected_count,
        raw_plan.boundary_selector.selected_count,
    )
    for active_axis in 1:3
        axis_matrices = ntuple(
            axis -> axis == active_axis ? kinetic_matrices[axis] : overlap_matrices[axis],
            3,
        )
        mode_matrix .+= _pqs_raw_product_box_selected_mode_matrix(
            raw_plan,
            axis_matrices,
        )
    end
    return mode_matrix
end

function _pqs_raw_product_box_overlap_diagnostics(
    raw_plan,
    metrics::NamedTuple{(:x,:y,:z)},
)
    axis_overlap_matrices = ntuple(
        axis -> _pqs_raw_product_box_axis_operator(
            raw_plan,
            metrics,
            axis,
            :overlap,
        ),
        3,
    )
    axis_overlap_errors = ntuple(
        axis -> LinearAlgebra.norm(
            axis_overlap_matrices[axis] -
            Matrix{Float64}(LinearAlgebra.I, size(axis_overlap_matrices[axis], 1), size(axis_overlap_matrices[axis], 2)),
            Inf,
        ),
        3,
    )
    source_mode_dims = ntuple(axis -> size(axis_overlap_matrices[axis], 1), 3)
    product_overlap = _pqs_raw_product_box_mode_matrix(
        raw_plan.source_mode_indices,
        axis_overlap_matrices,
    )
    source_mode_count = size(product_overlap, 1)
    max_product_overlap_error = LinearAlgebra.norm(
        product_overlap -
        Matrix{Float64}(LinearAlgebra.I, source_mode_count, source_mode_count),
        Inf,
    )
    boundary_overlap = _pqs_raw_product_box_selected_mode_matrix(
        raw_plan,
        axis_overlap_matrices,
    )
    retained_count = size(boundary_overlap, 1)
    selected_overlap_error = LinearAlgebra.norm(
        boundary_overlap -
        Matrix{Float64}(LinearAlgebra.I, retained_count, retained_count),
        Inf,
    )
    return (
        axis_overlap_errors = axis_overlap_errors,
        max_1d_source_overlap_error = maximum(axis_overlap_errors),
        max_product_overlap_error = max_product_overlap_error,
        selected_overlap_error = selected_overlap_error,
        overlap_identity_error = selected_overlap_error,
    )
end

function _pqs_product_box_one_dimensional_operator_factors(
    raw_plan,
    metrics::NamedTuple{(:x,:y,:z)},
)
    return ntuple(axis -> (
        overlap = _pqs_raw_product_box_axis_operator(
            raw_plan,
            metrics,
            axis,
            :overlap,
        ),
        position = _pqs_raw_product_box_axis_operator(
            raw_plan,
            metrics,
            axis,
            :position,
        ),
        x2 = _pqs_raw_product_box_axis_operator(
            raw_plan,
            metrics,
            axis,
            :x2,
        ),
        kinetic = _pqs_raw_product_box_axis_operator(
            raw_plan,
            metrics,
            axis,
            :kinetic,
        ),
    ), 3)
end

function _pqs_raw_product_box_axis_operator(
    raw_plan,
    metrics::NamedTuple{(:x,:y,:z)},
    axis::Int,
    kind::Symbol,
)
    return _pqs_raw_product_box_axis_operator(
        raw_plan,
        getproperty(metrics, (:x, :y, :z)[axis]),
        axis,
        kind,
    )
end

function _pqs_raw_product_box_axis_operator(
    raw_plan,
    axis_data,
    axis::Int,
    kind::Symbol,
)
    axis_name = (:x, :y, :z)[axis]
    hasproperty(axis_data, kind) || throw(
        ArgumentError("PQS raw product-box axis $(axis_name) metric data is missing $(kind)"),
    )
    axis_operator = Matrix{Float64}(getproperty(axis_data, kind))
    interval = raw_plan.axis_intervals[axis]
    first(interval) >= 1 && last(interval) <= size(axis_operator, 1) || throw(
        ArgumentError("PQS raw product-box axis $(axis_name) interval exceeds $(kind) matrix"),
    )
    coefficients = Matrix{Float64}(raw_plan.axis_local_coefficients[axis])
    size(coefficients, 1) == length(interval) || throw(
        DimensionMismatch("PQS raw product-box axis $(axis_name) coefficients must match interval length"),
    )
    size(coefficients, 2) >= 1 || throw(
        ArgumentError("PQS raw product-box axis $(axis_name) must retain at least one source mode"),
    )
    @views return Matrix{Float64}(
        transpose(coefficients) * axis_operator[interval, interval] * coefficients,
    )
end

function _pqs_product_box_support_overlap_matrix(
    support_states::AbstractVector{<:NTuple{3,Int}},
    metrics::NamedTuple{(:x,:y,:z)},
)
    overlap_x = metrics.x.overlap
    overlap_y = metrics.y.overlap
    overlap_z = metrics.z.overlap
    support_count = length(support_states)
    overlap = Matrix{Float64}(undef, support_count, support_count)
    @inbounds for col in 1:support_count
        jx, jy, jz = support_states[col]
        for row in 1:support_count
            ix, iy, iz = support_states[row]
            overlap[row, col] =
                overlap_x[ix, jx] *
                overlap_y[iy, jy] *
                overlap_z[iz, jz]
        end
    end
    return overlap
end

function _pqs_shared_raw_product_box_plan_status(
    shared_raw_product_box_plan,
)
    isnothing(shared_raw_product_box_plan) &&
        return :unavailable_descriptor_only_path_has_no_axis_bundles
    return :available
end

function _pqs_shared_raw_product_box_plan_unavailable_reason(
    shared_raw_product_box_plan,
)
    isnothing(shared_raw_product_box_plan) &&
        return :descriptor_only_path_has_no_axis_bundles
    return nothing
end

function _pqs_validated_shared_raw_product_box_plan(
    descriptor::_CartesianNestedProjectedQShellStagedUnitDescriptor3D,
    shared_raw_product_box_plan,
)
    isnothing(shared_raw_product_box_plan) && return nothing
    shared_raw_product_box_plan.object_kind == :cartesian_raw_product_box_plan_3d ||
        throw(ArgumentError("PQS raw plan requires a cartesian_raw_product_box_plan_3d"))
    source_mode_dims = ntuple(
        axis -> size(descriptor.axis_local_coefficients[axis], 2),
        3,
    )
    shared_raw_product_box_plan.axis_intervals == descriptor.axis_intervals ||
        throw(DimensionMismatch("shared raw product-box plan axis intervals must match PQS descriptor"))
    shared_raw_product_box_plan.source_mode_dims == source_mode_dims || throw(
        DimensionMismatch("shared raw product-box plan source dimensions must match PQS descriptor"),
    )
    shared_raw_product_box_plan.source_mode_count == prod(source_mode_dims) || throw(
        DimensionMismatch("shared raw product-box plan source mode count must match dimensions"),
    )
    length(shared_raw_product_box_plan.source_mode_indices) ==
        shared_raw_product_box_plan.source_mode_count || throw(
            DimensionMismatch("shared raw product-box plan source mode ordering length is inconsistent"),
        )
    all(axis -> size(shared_raw_product_box_plan.axis_local_coefficients[axis]) ==
                size(descriptor.axis_local_coefficients[axis]), 1:3) ||
        throw(DimensionMismatch("shared raw product-box plan coefficient shapes must match PQS descriptor"))
    all(axis -> isapprox(
        shared_raw_product_box_plan.axis_local_coefficients[axis],
        descriptor.axis_local_coefficients[axis];
        atol = 1.0e-10,
        rtol = 1.0e-10,
    ), 1:3) || throw(
        ArgumentError("shared raw product-box plan coefficients disagree with PQS descriptor"),
    )
    return shared_raw_product_box_plan
end

function _pqs_validated_raw_product_box_plan(raw_product_box_plan)
    raw_product_box_plan.object_kind == :cartesian_raw_product_box_plan_3d ||
        throw(ArgumentError("PQS raw-box retained rule requires a cartesian_raw_product_box_plan_3d"))
    raw_product_box_plan.source_mode_count == prod(raw_product_box_plan.source_mode_dims) ||
        throw(DimensionMismatch("raw product-box plan source mode count must match dimensions"))
    length(raw_product_box_plan.source_mode_indices) ==
        raw_product_box_plan.source_mode_count || throw(
            DimensionMismatch("raw product-box plan source mode ordering length is inconsistent"),
        )
    length(raw_product_box_plan.axis_local_coefficients) == 3 ||
        throw(DimensionMismatch("raw product-box plan must carry three axis transforms"))
    all(axis -> size(raw_product_box_plan.axis_local_coefficients[axis]) ==
                (length(raw_product_box_plan.axis_intervals[axis]),
                 raw_product_box_plan.source_mode_dims[axis]), 1:3) ||
        throw(DimensionMismatch("raw product-box plan coefficient shapes must match intervals and source dimensions"))
    return raw_product_box_plan
end

function _pqs_raw_product_box_boundary_selector(source_mode_dims::NTuple{3,Int})
    all(dim -> dim >= 2, source_mode_dims) || throw(
        ArgumentError("PQS raw-box boundary selector requires at least two source modes per axis"),
    )
    source_mode_indices = _pqs_raw_product_box_source_mode_indices(source_mode_dims)
    mode_indices = NTuple{3,Int}[]
    column_indices = Int[]
    for (column_index, mode_index) in pairs(source_mode_indices)
        if any(axis -> mode_index[axis] == 1 ||
                       mode_index[axis] == source_mode_dims[axis], 1:3)
            push!(mode_indices, mode_index)
            push!(column_indices, column_index)
        end
    end
    return (
        mode_indices = mode_indices,
        column_indices = column_indices,
        selection_rule = :any_axis_mode_index_first_or_last,
        retained_rule_contract = :RetainedRule,
        retained_rule_kind = :boundary_comx_product_mode_selection,
        selected_count = length(mode_indices),
        preserves_orthogonality = true,
    )
end

function _pqs_raw_product_box_structural_plan(
    descriptor::_CartesianNestedProjectedQShellStagedUnitDescriptor3D,
    shared_raw_product_box_plan = nothing,
)
    descriptor.kind == :projected_q_shell || throw(
        ArgumentError("PQS raw product-box plan requires a projected_q_shell descriptor"),
    )
    shared_plan =
        _pqs_validated_shared_raw_product_box_plan(descriptor, shared_raw_product_box_plan)
    source_mode_dims = isnothing(shared_plan) ?
                       ntuple(axis -> size(descriptor.axis_local_coefficients[axis], 2), 3) :
                       shared_plan.source_mode_dims
    for mode in descriptor.boundary_mode_indices
        all(axis -> 1 <= mode[axis] <= source_mode_dims[axis], 1:3) || throw(
            ArgumentError("PQS raw product-box boundary mode index exceeds source-mode dimensions"),
        )
    end
    boundary_selector = (
        mode_indices = descriptor.boundary_mode_indices,
        column_indices = descriptor.boundary_column_indices,
        selection_rule = descriptor.selection_rule,
        retained_rule_contract = :RetainedRule,
        retained_rule_kind = :boundary_comx_product_mode_selection,
        selected_count = descriptor.mode_count,
        preserves_orthogonality = true,
    )
    return (
        path = :pqs_raw_product_box_plan,
        representation = :orthogonal_raw_product_box,
        source_box_plan_contract = :RawProductBoxPlan,
        retained_rule_contract = :RetainedRule,
        retained_rule_kind = :boundary_comx_product_mode_selection,
        retained_rule_algorithmic = true,
        source_mode_dims = source_mode_dims,
        source_mode_count = isnothing(shared_plan) ? prod(source_mode_dims) : shared_plan.source_mode_count,
        axis_intervals = isnothing(shared_plan) ? descriptor.axis_intervals : shared_plan.axis_intervals,
        axis_local_coefficients =
            isnothing(shared_plan) ? descriptor.axis_local_coefficients : shared_plan.axis_local_coefficients,
        source_mode_indices = isnothing(shared_plan) ?
                              _pqs_raw_product_box_source_mode_indices(source_mode_dims) :
                              shared_plan.source_mode_indices,
        source_mode_ordering = isnothing(shared_plan) ?
                               :x_major_y_major_z_fast :
                               shared_plan.source_mode_ordering,
        shared_raw_product_box_plan = shared_plan,
        shared_raw_product_box_plan_available = !isnothing(shared_plan),
        shared_raw_product_box_plan_used = !isnothing(shared_plan),
        boundary_selector = boundary_selector,
        one_dimensional_operator_factors = nothing,
        operator_factors_available = false,
        axis_overlap_errors = nothing,
        max_1d_source_overlap_error = nothing,
        max_product_overlap_error = nothing,
        selected_overlap_error = nothing,
        overlap_identity_error = nothing,
        source_product_modes_orthogonal = nothing,
        row_projected_shell_support = false,
        lowdin_cleanup_used = false,
        diagnostics = (
            source = :pqs_raw_product_box_structural_plan,
            private_shadow_only = true,
            pqs_representation = :mode_selected_raw_product_box,
            source_box_plan_contract = :RawProductBoxPlan,
            retained_rule_contract = :RetainedRule,
            retained_rule_kind = :boundary_comx_product_mode_selection,
            retained_rule_algorithmic = true,
            source_mode_dims_are_total_lengths = true,
            source_mode_ordering = :x_major_y_major_z_fast,
            shared_raw_product_box_plan_available = !isnothing(shared_plan),
            shared_raw_product_box_plan_used = !isnothing(shared_plan),
            shared_raw_product_box_plan_status =
                _pqs_shared_raw_product_box_plan_status(shared_plan),
            shared_raw_product_box_plan_unavailable_reason =
                _pqs_shared_raw_product_box_plan_unavailable_reason(shared_plan),
            boundary_column_selection_only = true,
            raw_product_box_operators_use_1d_factors = false,
            operator_factors_available = false,
            shell_projection_used = false,
            lowdin_cleanup_used = false,
            support_coefficient_matrix_used = false,
            retained_weight_semantics = :not_positive_quadrature_weights,
            ida_weight_division_allowed = false,
            packet_adoption = false,
        ),
    )
end

function _pqs_raw_product_box_structural_plan_from_raw_product_box_plan(
    raw_product_box_plan,
)
    raw_plan = _pqs_validated_raw_product_box_plan(raw_product_box_plan)
    source_mode_dims = raw_plan.source_mode_dims
    boundary_selector = _pqs_raw_product_box_boundary_selector(source_mode_dims)
    return (
        path = :pqs_raw_product_box_plan,
        representation = :orthogonal_raw_product_box,
        source_box_plan_contract = :RawProductBoxPlan,
        retained_rule_contract = :RetainedRule,
        retained_rule_kind = :boundary_comx_product_mode_selection,
        retained_rule_algorithmic = true,
        source_mode_dims = source_mode_dims,
        source_mode_count = raw_plan.source_mode_count,
        axis_intervals = raw_plan.axis_intervals,
        axis_local_coefficients = raw_plan.axis_local_coefficients,
        source_mode_indices = raw_plan.source_mode_indices,
        source_mode_ordering = raw_plan.source_mode_ordering,
        shared_raw_product_box_plan = raw_plan,
        shared_raw_product_box_plan_available = true,
        shared_raw_product_box_plan_used = true,
        boundary_selector = boundary_selector,
        one_dimensional_operator_factors = nothing,
        operator_factors_available = false,
        axis_overlap_errors = nothing,
        max_1d_source_overlap_error = nothing,
        max_product_overlap_error = nothing,
        selected_overlap_error = nothing,
        overlap_identity_error = nothing,
        source_product_modes_orthogonal = nothing,
        row_projected_shell_support = false,
        lowdin_cleanup_used = false,
        diagnostics = (
            source = :pqs_raw_product_box_structural_plan_from_raw_product_box_plan,
            private_shadow_only = true,
            pqs_representation = :mode_selected_raw_product_box,
            source_box_plan_contract = :RawProductBoxPlan,
            retained_rule_contract = :RetainedRule,
            retained_rule_kind = :boundary_comx_product_mode_selection,
            retained_rule_algorithmic = true,
            source_mode_dims_are_total_lengths = true,
            source_mode_ordering = raw_plan.source_mode_ordering,
            shared_raw_product_box_plan_available = true,
            shared_raw_product_box_plan_used = true,
            shared_raw_product_box_plan_status = :available,
            shared_raw_product_box_plan_unavailable_reason = nothing,
            descriptor_required = false,
            boundary_column_selection_only = true,
            raw_product_box_operators_use_1d_factors = false,
            operator_factors_available = false,
            shell_projection_used = false,
            lowdin_cleanup_used = false,
            support_coefficient_matrix_used = false,
            retained_weight_semantics = :not_positive_quadrature_weights,
            ida_weight_division_allowed = false,
            packet_adoption = false,
            fixed_block_routing = false,
            qwhamiltonian_consumes = false,
            public_default_consumes = false,
            cr2_science_status_changed = false,
        ),
    )
end

function _pqs_raw_product_box_plan(
    descriptor::_CartesianNestedProjectedQShellStagedUnitDescriptor3D,
)
    return _pqs_raw_product_box_structural_plan(descriptor)
end

function _pqs_raw_product_box_plan(
    descriptor::_CartesianNestedProjectedQShellStagedUnitDescriptor3D,
    shared_raw_product_box_plan,
)
    return _pqs_raw_product_box_structural_plan(
        descriptor,
        shared_raw_product_box_plan,
    )
end

function _pqs_raw_product_box_operator_plan_from_structural_plan(
    structural_plan,
    metrics::NamedTuple{(:x,:y,:z)};
    orthogonality_atol::Real = 1.0e-8,
    source::Symbol = :pqs_raw_product_box_plan,
)
    structural_plan.representation == :orthogonal_raw_product_box || throw(
        ArgumentError("PQS raw product-box operator plan requires an orthogonal raw product-box structural plan"),
    )
    one_dimensional_operator_factors =
        _pqs_product_box_one_dimensional_operator_factors(structural_plan, metrics)
    overlap_diagnostics = _pqs_raw_product_box_overlap_diagnostics(
        structural_plan,
        metrics,
    )
    source_product_modes_orthogonal =
        overlap_diagnostics.max_1d_source_overlap_error <= orthogonality_atol &&
        overlap_diagnostics.max_product_overlap_error <= orthogonality_atol &&
        overlap_diagnostics.selected_overlap_error <= orthogonality_atol
    return merge(
        structural_plan,
        (
            one_dimensional_operator_factors = (
                x = one_dimensional_operator_factors[1],
                y = one_dimensional_operator_factors[2],
                z = one_dimensional_operator_factors[3],
            ),
            operator_factors_available = true,
            axis_overlap_errors = overlap_diagnostics.axis_overlap_errors,
            max_1d_source_overlap_error =
                overlap_diagnostics.max_1d_source_overlap_error,
            max_product_overlap_error =
                overlap_diagnostics.max_product_overlap_error,
            selected_overlap_error = overlap_diagnostics.selected_overlap_error,
            overlap_identity_error = overlap_diagnostics.overlap_identity_error,
            source_product_modes_orthogonal = source_product_modes_orthogonal,
            diagnostics = merge(
                structural_plan.diagnostics,
                (
                    source = source,
                    raw_product_box_operators_use_1d_factors = true,
                    operator_factors_available = true,
                    operator_factor_source = :explicit_axis_metric_data,
                    operator_metric_sources =
                        _cartesian_source_box_metric_sources(metrics),
                    raw_product_box_integration_contract =
                        isnothing(structural_plan.shared_raw_product_box_plan) ?
                        :unavailable_descriptor_only_path :
                        structural_plan.shared_raw_product_box_plan.diagnostics.integration_contract,
                    raw_product_box_numerical_reference_fallback =
                        isnothing(structural_plan.shared_raw_product_box_plan) ?
                        nothing :
                        structural_plan.shared_raw_product_box_plan.diagnostics.numerical_reference_fallback,
                    source_product_modes_orthogonal =
                        source_product_modes_orthogonal,
                    max_1d_source_overlap_error =
                        overlap_diagnostics.max_1d_source_overlap_error,
                    max_product_overlap_error =
                        overlap_diagnostics.max_product_overlap_error,
                    selected_overlap_error =
                        overlap_diagnostics.selected_overlap_error,
                ),
            ),
        ),
    )
end

function _pqs_raw_product_box_plan(
    descriptor::_CartesianNestedProjectedQShellStagedUnitDescriptor3D,
    metrics::NamedTuple{(:x,:y,:z)};
    orthogonality_atol::Real = 1.0e-8,
    shared_raw_product_box_plan = nothing,
)
    structural_plan = _pqs_raw_product_box_structural_plan(
        descriptor,
        shared_raw_product_box_plan,
    )
    return _pqs_raw_product_box_operator_plan_from_structural_plan(
        structural_plan,
        metrics;
        orthogonality_atol,
        source = :pqs_raw_product_box_plan,
    )
end

function _pqs_raw_product_box_plan(
    descriptor::_CartesianNestedProjectedQShellStagedUnitDescriptor3D,
    shared_raw_product_box_plan,
    metrics::NamedTuple{(:x,:y,:z)};
    orthogonality_atol::Real = 1.0e-8,
)
    return _pqs_raw_product_box_plan(
        descriptor,
        metrics;
        orthogonality_atol,
        shared_raw_product_box_plan,
    )
end

function _pqs_raw_product_box_plan_from_raw_product_box_plan(
    raw_product_box_plan,
    metrics::NamedTuple{(:x,:y,:z)};
    orthogonality_atol::Real = 1.0e-8,
)
    structural_plan =
        _pqs_raw_product_box_structural_plan_from_raw_product_box_plan(
            raw_product_box_plan,
        )
    return _pqs_raw_product_box_operator_plan_from_structural_plan(
        structural_plan,
        metrics;
        orthogonality_atol,
        source = :pqs_raw_product_box_plan_from_raw_product_box_plan,
    )
end

function _pqs_shell_realization_plan(
    descriptor::_CartesianNestedProjectedQShellStagedUnitDescriptor3D,
    metrics::NamedTuple{(:x,:y,:z)};
    isometry_atol::Real = 1.0e-8,
)
    descriptor.kind == :projected_q_shell || throw(
        ArgumentError("PQS shell-realization plan requires a projected_q_shell descriptor"),
    )
    shell_projection_matrix =
        _nested_projected_q_shell_descriptor_seed_coefficients(descriptor)
    size(shell_projection_matrix) == (descriptor.support_count, descriptor.mode_count) ||
        throw(
            DimensionMismatch("PQS shell projection matrix must map selected product-box modes to shell support rows"),
        )
    lowdin_cleanup = Matrix{Float64}(descriptor.cleanup_transform)
    size(lowdin_cleanup) == (descriptor.mode_count, descriptor.retained_count) ||
        throw(
            DimensionMismatch("PQS shell realization Lowdin cleanup dimensions must match selected/retained counts"),
        )
    shell_overlap_matrix = _pqs_product_box_support_overlap_matrix(
        descriptor.support_states,
        metrics,
    )
    shell_projection_gram =
        transpose(shell_projection_matrix) * shell_overlap_matrix * shell_projection_matrix
    shell_isometry_matrix = shell_projection_matrix * lowdin_cleanup
    realized_overlap =
        transpose(shell_isometry_matrix) * shell_overlap_matrix * shell_isometry_matrix
    isometry_error = LinearAlgebra.norm(
        realized_overlap -
        Matrix{Float64}(
            LinearAlgebra.I,
            descriptor.retained_count,
            descriptor.retained_count,
        ),
        Inf,
    )
    isfinite(isometry_error) || throw(
        ArgumentError("PQS shell-realization plan produced a non-finite isometry error"),
    )
    return (
        path = :pqs_shell_realization_plan,
        representation = :shell_projection_lowdin_isometry,
        shell_projection_used = true,
        lowdin_cleanup_used = true,
        shell_projection_matrix = shell_projection_matrix,
        shell_projection_gram = shell_projection_gram,
        lowdin_cleanup = lowdin_cleanup,
        cleanup_method = descriptor.cleanup_method,
        cleanup_matrix_size = descriptor.cleanup_matrix_size,
        shell_isometry_matrix = shell_isometry_matrix,
        realized_overlap = realized_overlap,
        isometry_error = isometry_error,
        isometric = isometry_error <= isometry_atol,
        diagnostics = (
            source = :pqs_shell_realization_plan,
            private_shadow_only = true,
            shell_projection_used = true,
            lowdin_cleanup_used = true,
            shell_projection_realization_requires_lowdin = true,
            shell_projection_realization_isometric =
                isometry_error <= isometry_atol,
            raw_product_box_operators_use_1d_factors = false,
            packet_adoption = false,
        ),
    )
end

function _pqs_product_box_realization_plan(
    descriptor::_CartesianNestedProjectedQShellStagedUnitDescriptor3D,
    metrics::NamedTuple{(:x,:y,:z)};
    orthogonality_atol::Real = 1.0e-8,
    isometry_atol::Real = 1.0e-8,
    shared_raw_product_box_plan = nothing,
)
    raw_plan = _pqs_raw_product_box_plan(
        descriptor,
        metrics;
        orthogonality_atol = orthogonality_atol,
        shared_raw_product_box_plan = shared_raw_product_box_plan,
    )
    shell_plan = _pqs_shell_realization_plan(
        descriptor,
        metrics;
        isometry_atol = isometry_atol,
    )
    return (
        raw_product_box_plan = raw_plan,
        source_box_plan = raw_plan,
        boundary_selector = raw_plan.boundary_selector,
        one_dimensional_operator_factors = raw_plan.one_dimensional_operator_factors,
        shell_realization_plan = shell_plan,
        shell_projection_matrix = shell_plan.shell_projection_matrix,
        lowdin_cleanup = shell_plan.lowdin_cleanup,
        isometry_error = shell_plan.isometry_error,
        diagnostics = (
            source = :pqs_product_box_realization_plan,
            private_shadow_only = true,
            raw_product_box_stage_lowdin_cleanup_used = false,
            raw_product_box_operators_use_1d_factors = true,
            shell_projection_used_for_raw_box_operators = false,
            shell_projection_realization_available = true,
            shell_projection_realization_requires_lowdin = true,
            shell_projection_realization_isometric =
                shell_plan.isometry_error <= isometry_atol,
            retained_functions_live_in_product_box_mode_span = true,
            shell_realized_functions_live_in_shell_row_support_subspace = true,
            retained_weight_semantics = :not_positive_quadrature_weights,
            ida_weight_division_allowed = false,
            packet_adoption = false,
            fixed_block_sidecar_installation = false,
            qwhamiltonian_consumes = false,
            public_default_consumes = false,
            cr2_science_status_changed = false,
            generic_retained_unit_framework = false,
        ),
    )
end

function _pqs_raw_product_box_plan_view(pqs_plan)
    if hasproperty(pqs_plan, :representation) &&
       pqs_plan.representation == :orthogonal_raw_product_box
        return pqs_plan
    elseif hasproperty(pqs_plan, :raw_product_box_plan)
        return pqs_plan.raw_product_box_plan
    elseif hasproperty(pqs_plan, :source_box_plan)
        return pqs_plan.source_box_plan
    end
    throw(ArgumentError("expected a PQS raw product-box plan"))
end

function _pqs_product_source_box_project_axis_terms(
    operator_terms::AbstractArray{<:Real,3},
    pqs_plan,
    product_axis,
    axis::Int;
    label::AbstractString,
)
    raw_plan = _pqs_raw_product_box_plan_view(pqs_plan)
    nterms = size(operator_terms, 1)
    nterms > 0 || throw(
        ArgumentError("$(label) requires at least one term"),
    )
    pqs_interval = raw_plan.axis_intervals[axis]
    product_interval = _staged_axis_interval(product_axis)
    first(pqs_interval) >= 1 && last(pqs_interval) <= size(operator_terms, 2) ||
        throw(
            ArgumentError("$(label) PQS interval exceeds axis $(axis) term table"),
        )
    first(product_interval) >= 1 && last(product_interval) <= size(operator_terms, 3) ||
        throw(
            ArgumentError("$(label) product interval exceeds axis $(axis) term table"),
        )
    pqs_coefficients = Matrix{Float64}(raw_plan.axis_local_coefficients[axis])
    product_coefficients = Matrix{Float64}(product_axis.coefficient_matrix)
    size(pqs_coefficients, 1) == length(pqs_interval) || throw(
        DimensionMismatch("$(label) PQS axis coefficients must match interval length"),
    )
    size(product_coefficients, 1) == length(product_interval) || throw(
        DimensionMismatch("$(label) product axis coefficients must match interval length"),
    )
    projected = Array{Float64,3}(
        undef,
        nterms,
        size(pqs_coefficients, 2),
        size(product_coefficients, 2),
    )
    @inbounds for term in 1:nterms
        term_matrix = @view operator_terms[term, pqs_interval, product_interval]
        projected[term, :, :] .=
            transpose(pqs_coefficients) * term_matrix * product_coefficients
    end
    return projected
end

function _pqs_product_source_box_project_local_gaussian_axis_terms(
    operator_terms::AbstractArray{<:Real,3},
    pqs_plan,
    product_axis,
    axis::Int,
)
    return _pqs_product_source_box_project_axis_terms(
        operator_terms,
        pqs_plan,
        product_axis,
        axis;
        label = "PQS/product source-box local-Gaussian terms",
    )
end

function _pqs_product_source_box_local_gaussian_axis_factors(
    pqs_plan,
    product_unit::_CartesianNestedProductStagedByCenterUnit3D,
    axis_gaussian_terms::NamedTuple{(:x,:y,:z)},
)
    return ntuple(axis -> begin
        terms = getproperty(axis_gaussian_terms, (:x, :y, :z)[axis])
        _pqs_product_source_box_project_local_gaussian_axis_terms(
            terms,
            pqs_plan,
            product_unit.axes[axis],
            axis,
        )
    end, 3)
end

function _pqs_product_source_box_density_density_axis_factors(
    pqs_plan,
    product_unit::_CartesianNestedProductStagedByCenterUnit3D,
    axis_pair_factor_terms::NamedTuple{(:x,:y,:z)},
)
    return ntuple(axis -> begin
        terms = getproperty(axis_pair_factor_terms, (:x, :y, :z)[axis])
        _pqs_product_source_box_project_axis_terms(
            terms,
            pqs_plan,
            product_unit.axes[axis],
            axis;
            label = "PQS/product source-box density-density pair-factor terms",
        )
    end, 3)
end

function _pqs_product_source_box_density_density_weight_views(
    raw_plan,
    product_retained_unit_plan,
    axis_weights::NamedTuple{(:x,:y,:z)},
)
    return (
        pqs = ntuple(axis -> _source_box_axis_positive_weights(
            getproperty(axis_weights, (:x, :y, :z)[axis]),
            raw_plan.axis_intervals[axis];
            axis_name = (:x, :y, :z)[axis],
            side = :pqs,
        ), 3),
        product = ntuple(axis -> _source_box_axis_positive_weights(
            getproperty(axis_weights, (:x, :y, :z)[axis]),
            product_retained_unit_plan.source_axis_intervals[axis];
            axis_name = (:x, :y, :z)[axis],
            side = :product,
        ), 3),
    )
end

function _pqs_product_source_box_density_density_interaction_block(
    pqs_plan,
    product_unit::_CartesianNestedProductStagedByCenterUnit3D;
    term_coefficients::AbstractVector{<:Real},
    axis_pair_factor_terms::NamedTuple{(:x,:y,:z)},
    axis_weights::NamedTuple{(:x,:y,:z)},
    pair_factor_normalization::Symbol = :density_normalized,
)
    pair_factor_normalization == :density_normalized || throw(
        ArgumentError(
            "PQS/product source-box density-density fixture currently requires density-normalized pair factors",
        ),
    )
    raw_plan = _pqs_raw_product_box_plan_view(pqs_plan)
    raw_plan.representation == :orthogonal_raw_product_box || throw(
        ArgumentError("PQS/product source-box density-density block requires a raw product-box PQS plan"),
    )
    product_unit.kind == :product_doside || throw(
        ArgumentError("PQS/product source-box density-density block requires a product_doside unit"),
    )
    length(product_unit.axis_function_indices) == length(product_unit.column_range) ||
        throw(
            ArgumentError("PQS/product source-box density-density product unit axis metadata does not match its column range"),
        )
    nterms = length(term_coefficients)
    nterms > 0 || throw(
        ArgumentError("PQS/product source-box density-density block requires at least one term coefficient"),
    )
    for axis_name in (:x, :y, :z)
        size(getproperty(axis_pair_factor_terms, axis_name), 1) == nterms || throw(
            ArgumentError("PQS/product source-box density-density pair-factor term count mismatch on $(axis_name)"),
        )
    end
    product_retained_unit_plan = _product_doside_retained_unit_plan(product_unit)
    source_weights = _pqs_product_source_box_density_density_weight_views(
        raw_plan,
        product_retained_unit_plan,
        axis_weights,
    )
    projected_terms = _pqs_product_source_box_density_density_axis_factors(
        raw_plan,
        product_unit,
        axis_pair_factor_terms,
    )
    coeffs = Float64[Float64(value) for value in term_coefficients]
    pqs_modes = raw_plan.boundary_selector.mode_indices
    product_modes = product_unit.axis_function_indices
    length(pqs_modes) == raw_plan.boundary_selector.selected_count || throw(
        ArgumentError("PQS/product source-box density-density boundary mode count must match selected count"),
    )
    block = zeros(Float64, length(pqs_modes), length(product_modes))
    projected_x, projected_y, projected_z = projected_terms
    @inbounds for col in eachindex(product_modes)
        px, py, pz = product_modes[col]
        for row in eachindex(pqs_modes)
            qx, qy, qz = pqs_modes[row]
            value = 0.0
            @simd for term in 1:nterms
                value +=
                    coeffs[term] *
                    projected_x[term, qx, px] *
                    projected_y[term, qy, py] *
                    projected_z[term, qz, pz]
            end
            block[row, col] = value
        end
    end
    all(isfinite, block) || throw(
        ArgumentError("PQS/product source-box density-density block produced non-finite entries"),
    )
    pair_plan = (
        pair_kind = :pqs_product_source_box_density_density_pair,
        left_source_family = :mode_selected_raw_product_box,
        right_source_family = :product_doside,
        left_retained_rule_kind = raw_plan.retained_rule_kind,
        right_retained_rule_kind = product_retained_unit_plan.retained_rule_kind,
        left_source_dimensions = raw_plan.source_mode_dims,
        right_source_dimensions = product_retained_unit_plan.source_axis_lengths,
        left_source_dimension = raw_plan.source_mode_count,
        right_source_dimension = product_retained_unit_plan.source_dimension,
        left_retained_count = raw_plan.boundary_selector.selected_count,
        right_retained_count = product_retained_unit_plan.retained_count,
        left_column_range = nothing,
        right_column_range = product_retained_unit_plan.column_range,
        axis_intervals = (
            pqs = raw_plan.axis_intervals,
            product = product_retained_unit_plan.source_axis_intervals,
        ),
        pqs_boundary_mode_selector = raw_plan.boundary_selector,
        product_retained_unit_plan = product_retained_unit_plan,
        product_retained_transform = product_retained_unit_plan,
        source_weights = source_weights,
        supported_terms = (:pair_sum,),
        diagnostics = (
            source = :pqs_product_source_box_density_density_pair_plan,
            private_shadow_only = true,
            source_box_pair_plan = true,
            interaction_operator = :electron_electron_density_density,
            output_representation = :two_index_density_density,
            four_index_galerkin_tensor = false,
            pqs_representation = :mode_selected_raw_product_box,
            product_representation = :product_doside,
            pqs_boundary_mode_selection_used = true,
            product_doside_retained_transform_used = true,
            product_doside_retained_unit_plan_used = true,
            shared_raw_product_box_plan_available =
                raw_plan.shared_raw_product_box_plan_available,
            shared_raw_product_box_plan_used =
                raw_plan.shared_raw_product_box_plan_used,
            source_mode_ordering = raw_plan.source_mode_ordering,
            operator_factor_source = :explicit_density_normalized_pair_factor_terms,
            input_pair_factor_data = :caller_supplied_explicit_data,
            input_pair_factor_data_pgdg_checked = false,
            raw_source_weights_available = true,
            density_normalized_pair_factors = true,
            raw_weighted_pair_factors = false,
            pair_factor_normalization = pair_factor_normalization,
            source_weight_division_owner = :caller_supplied_density_normalized_pair_factors,
            source_weight_division_applied_by_helper = false,
            raw_product_box_operators_use_1d_factors = true,
            shell_projection_used = false,
            lowdin_cleanup_used = false,
            support_coefficient_matrix_used = false,
            shell_row_algorithm = false,
            support_local_oracle_used = false,
            support_local_pqs_oracle_used = false,
            retained_pqs_weights_used = false,
            retained_pqs_weights_positive_checked = false,
            retained_weight_division_allowed = false,
            retained_pqs_weight_division_allowed = false,
            ida_weight_division_allowed = false,
            ida_mwg_semantics_changed = false,
            mwg_ida_semantics_changed = false,
            retained_weight_semantics_changed = false,
            packet_adoption = false,
            fixed_block_routing = false,
            qwhamiltonian_consumes = false,
            public_default_consumes = false,
            cr2_science_status_changed = false,
            numerical_reference_fallback = false,
            ecp_terms_implemented = false,
            local_gaussian_one_body_implemented = false,
            mwg_interaction_implemented = false,
            generic_retained_unit_framework = false,
        ),
    )
    return (
        path = :pqs_product_source_box_density_density_interaction,
        interaction_operator = :electron_electron_density_density,
        block = block,
        pair_plan = pair_plan,
        one_dimensional_pair_factors = (
            x = projected_x,
            y = projected_y,
            z = projected_z,
        ),
        diagnostics = merge(
            pair_plan.diagnostics,
            (
                source = :pqs_product_source_box_density_density_interaction_block,
                source_box_first = true,
                output_finite = true,
                electron_electron_terms_implemented = true,
                term_count = nterms,
                axis_term_dimensions = (
                    x = size(axis_pair_factor_terms.x),
                    y = size(axis_pair_factor_terms.y),
                    z = size(axis_pair_factor_terms.z),
                ),
                projected_axis_term_dimensions = (
                    x = size(projected_x),
                    y = size(projected_y),
                    z = size(projected_z),
                ),
            ),
        ),
    )
end

function _pqs_product_source_box_raw_weighted_density_density_interaction_block(
    pqs_plan,
    product_unit::_CartesianNestedProductStagedByCenterUnit3D;
    term_coefficients::AbstractVector{<:Real},
    raw_axis_pair_factor_terms::NamedTuple{(:x,:y,:z)},
    axis_weights::NamedTuple{(:x,:y,:z)},
)
    density_normalized_terms = _source_box_density_normalized_axis_pair_terms(
        raw_axis_pair_factor_terms,
        axis_weights,
    )
    density_normalized_core =
        _pqs_product_source_box_density_density_interaction_block(
            pqs_plan,
            product_unit;
            term_coefficients,
            axis_pair_factor_terms = density_normalized_terms,
            axis_weights,
            pair_factor_normalization = :density_normalized,
        )
    return (
        path = :pqs_product_source_box_raw_weighted_density_density_interaction,
        interaction_operator = :electron_electron_density_density,
        block = density_normalized_core.block,
        density_normalized_core = density_normalized_core,
        normalized_axis_pair_factor_terms = density_normalized_terms,
        raw_axis_pair_factor_terms = raw_axis_pair_factor_terms,
        pair_plan = density_normalized_core.pair_plan,
        diagnostics = merge(
            density_normalized_core.diagnostics,
            (
                source = :pqs_product_source_box_raw_weighted_density_density_interaction_block,
                path = :pqs_product_source_box_raw_weighted_density_density_interaction,
                pair_factor_normalization = :raw_weighted,
                raw_weighted_pair_factors = true,
                density_normalized_pair_factors = false,
                density_normalized_pair_factors_generated = true,
                source_weight_division_owner = :source_box_raw_weights,
                source_weight_division_applied_by_helper = true,
                source_weight_division_shape = :axis_pair_weight_outer,
                density_normalized_core_helper =
                    :_pqs_product_source_box_density_density_interaction_block,
                retained_pqs_weights_used = false,
                retained_weight_division_allowed = false,
                retained_pqs_weight_division_allowed = false,
                ida_weight_division_allowed = false,
                ida_mwg_semantics_changed = false,
                mwg_ida_semantics_changed = false,
                packet_adoption = false,
                qwhamiltonian_consumes = false,
                numerical_reference_fallback = false,
            ),
        ),
    )
end

function _pqs_product_source_box_local_gaussian_sum_block(
    pqs_plan,
    product_unit::_CartesianNestedProductStagedByCenterUnit3D;
    term_coefficients::AbstractVector{<:Real},
    axis_gaussian_terms::NamedTuple{(:x,:y,:z)},
)
    raw_plan = _pqs_raw_product_box_plan_view(pqs_plan)
    raw_plan.representation == :orthogonal_raw_product_box || throw(
        ArgumentError("PQS/product source-box local-Gaussian block requires a raw product-box PQS plan"),
    )
    product_unit.kind == :product_doside || throw(
        ArgumentError("PQS/product source-box local-Gaussian block requires a product_doside unit"),
    )
    length(product_unit.axis_function_indices) == length(product_unit.column_range) ||
        throw(
            ArgumentError("PQS/product source-box local-Gaussian product unit axis metadata does not match its column range"),
        )
    nterms = length(term_coefficients)
    nterms > 0 || throw(
        ArgumentError("PQS/product source-box local-Gaussian block requires at least one term coefficient"),
    )
    for axis_name in (:x, :y, :z)
        size(getproperty(axis_gaussian_terms, axis_name), 1) == nterms || throw(
            ArgumentError("PQS/product source-box local-Gaussian axis term count mismatch on $(axis_name)"),
        )
    end
    product_retained_unit_plan = _product_doside_retained_unit_plan(product_unit)
    projected_terms =
        _pqs_product_source_box_local_gaussian_axis_factors(
            raw_plan,
            product_unit,
            axis_gaussian_terms,
        )
    coeffs = Float64[Float64(value) for value in term_coefficients]
    pqs_modes = raw_plan.boundary_selector.mode_indices
    product_modes = product_unit.axis_function_indices
    length(pqs_modes) == raw_plan.boundary_selector.selected_count || throw(
        ArgumentError("PQS/product source-box local-Gaussian boundary mode count must match selected count"),
    )
    block = zeros(Float64, length(pqs_modes), length(product_modes))
    projected_x, projected_y, projected_z = projected_terms
    @inbounds for col in eachindex(product_modes)
        px, py, pz = product_modes[col]
        for row in eachindex(pqs_modes)
            qx, qy, qz = pqs_modes[row]
            value = 0.0
            @simd for term in 1:nterms
                value +=
                    coeffs[term] *
                    projected_x[term, qx, px] *
                    projected_y[term, qy, py] *
                    projected_z[term, qz, pz]
            end
            block[row, col] = value
        end
    end
    all(isfinite, block) || throw(
        ArgumentError("PQS/product source-box local-Gaussian block produced non-finite entries"),
    )
    pair_plan = (
        pair_kind = :pqs_product_source_box_local_gaussian_pair,
        left_source_family = :mode_selected_raw_product_box,
        right_source_family = :product_doside,
        left_retained_rule_kind = raw_plan.retained_rule_kind,
        right_retained_rule_kind = product_retained_unit_plan.retained_rule_kind,
        left_source_dimensions = raw_plan.source_mode_dims,
        right_source_dimensions = product_retained_unit_plan.source_axis_lengths,
        left_source_dimension = raw_plan.source_mode_count,
        right_source_dimension = product_retained_unit_plan.source_dimension,
        left_retained_count = raw_plan.boundary_selector.selected_count,
        right_retained_count = product_retained_unit_plan.retained_count,
        left_column_range = nothing,
        right_column_range = product_retained_unit_plan.column_range,
        axis_intervals = (
            pqs = raw_plan.axis_intervals,
            product = product_retained_unit_plan.source_axis_intervals,
        ),
        pqs_boundary_mode_selector = raw_plan.boundary_selector,
        product_retained_unit_plan = product_retained_unit_plan,
        product_retained_transform = product_retained_unit_plan,
        supported_terms = (:gaussian_sum,),
        diagnostics = (
            source = :pqs_product_source_box_local_gaussian_pair_plan,
            private_shadow_only = true,
            source_box_pair_plan = true,
            pqs_representation = :mode_selected_raw_product_box,
            pqs_boundary_mode_selection_used = true,
            product_doside_retained_transform_used = true,
            product_doside_retained_unit_plan_used = true,
            shared_raw_product_box_plan_available =
                raw_plan.shared_raw_product_box_plan_available,
            shared_raw_product_box_plan_used =
                raw_plan.shared_raw_product_box_plan_used,
            source_mode_ordering = raw_plan.source_mode_ordering,
            operator_factor_source = :explicit_local_gaussian_axis_terms,
            input_local_gaussian_data = :caller_supplied_explicit_data,
            input_local_gaussian_data_pgdg_checked = false,
            pgdg_analytic_operator_provenance_claimed = false,
            raw_product_box_numerical_reference_fallback =
                hasproperty(raw_plan.diagnostics, :raw_product_box_numerical_reference_fallback) ?
                raw_plan.diagnostics.raw_product_box_numerical_reference_fallback :
                nothing,
            numerical_reference_fallback = false,
            raw_product_box_operators_use_1d_factors = true,
            shell_projection_used = false,
            lowdin_cleanup_used = false,
            support_coefficient_matrix_used = false,
            shell_row_algorithm = false,
            support_local_oracle_used = false,
            support_local_pqs_oracle_used = false,
            retained_pqs_weights_used = false,
            retained_pqs_weights_positive_checked = false,
            retained_pqs_weight_division_allowed = false,
            ida_weight_division_allowed = false,
            ida_mwg_semantics_changed = false,
            retained_weight_semantics_changed = false,
            packet_adoption = false,
            fixed_block_routing = false,
            qwhamiltonian_consumes = false,
            public_default_consumes = false,
            cr2_science_status_changed = false,
            ecp_terms_implemented = false,
            electron_electron_terms_implemented = false,
            mwg_interaction_implemented = false,
            generic_retained_unit_framework = false,
        ),
    )
    return (
        path = :pqs_product_source_box_local_gaussian_sum,
        block = block,
        pair_plan = pair_plan,
        one_dimensional_gaussian_factors = (
            x = projected_x,
            y = projected_y,
            z = projected_z,
        ),
        diagnostics = merge(
            pair_plan.diagnostics,
            (
                source = :pqs_product_source_box_local_gaussian_sum_block,
                source_box_first = true,
                local_gaussian_source_box_terms = true,
                positive_gaussian_sum_convention = true,
                nuclear_charge_applied = false,
                nuclear_attraction_sign_applied = false,
                term_coefficients_source = :caller_supplied_explicit_data,
                axis_gaussian_terms_source = :caller_supplied_explicit_data,
                input_local_gaussian_data_pgdg_checked = false,
                pgdg_analytic_operator_provenance_claimed = false,
                numerical_reference_fallback = false,
                term_count = nterms,
                axis_term_dimensions = (
                    x = size(axis_gaussian_terms.x),
                    y = size(axis_gaussian_terms.y),
                    z = size(axis_gaussian_terms.z),
                ),
                projected_axis_term_dimensions = (
                    x = size(projected_x),
                    y = size(projected_y),
                    z = size(projected_z),
                ),
                operator_factor_source = :explicit_local_gaussian_axis_terms,
                dense_raw_source_box_pair_matrix_materialized = false,
                dense_raw_pair_storage_avoided = true,
                retained_block_assembled_directly_from_1d_factors = true,
                source_box_pair_storage_scaling =
                    :one_dimensional_factors_plus_retained_block,
                pqs_representation = :mode_selected_raw_product_box,
                shell_projection_used = false,
                lowdin_cleanup_used = false,
                support_coefficient_matrix_used = false,
                shell_row_algorithm = false,
                support_local_oracle_used = false,
                support_local_pqs_oracle_used = false,
                retained_pqs_weight_division_allowed = false,
                ida_weight_division_allowed = false,
                ida_mwg_semantics_changed = false,
                ecp_terms_implemented = false,
                electron_electron_terms_implemented = false,
                mwg_interaction_implemented = false,
                local_gaussian_one_body_implemented = true,
                product_staged_nuclear_execution_changed = false,
                packet_adoption = false,
                fixed_block_routing = false,
                qwhamiltonian_consumes = false,
                public_default_consumes = false,
                cr2_science_status_changed = false,
                output_finite = true,
            ),
        ),
    )
end

function _source_box_nuclear_attraction_center_values(centers)
    center_items = collect(centers)
    !isempty(center_items) || throw(
        ArgumentError("source-box nuclear attraction requires at least one center"),
    )
    return map(center_items) do center
        value = try
            ntuple(axis -> Float64(center[axis]), 3)
        catch err
            throw(
                ArgumentError("source-box nuclear attraction centers must be 3-coordinate tuples"),
            )
        end
        all(isfinite, value) || throw(
            ArgumentError("source-box nuclear attraction centers must be finite"),
        )
        value
    end
end

function _source_box_nuclear_attraction_charge_values(nuclear_charges)
    charges = Float64[Float64(charge) for charge in nuclear_charges]
    !isempty(charges) || throw(
        ArgumentError("source-box nuclear attraction requires at least one nuclear charge"),
    )
    all(isfinite, charges) || throw(
        ArgumentError("source-box nuclear attraction charges must be finite"),
    )
    all(>=(0.0), charges) || throw(
        ArgumentError("source-box nuclear attraction charges must be nonnegative"),
    )
    return charges
end

function _source_box_nuclear_attraction_by_center(
    centered_positive_block_builder;
    centers,
    nuclear_charges,
    path::Symbol,
    source::Symbol,
)
    center_values = _source_box_nuclear_attraction_center_values(centers)
    charge_values = _source_box_nuclear_attraction_charge_values(nuclear_charges)
    length(center_values) == length(charge_values) || throw(
        ArgumentError("source-box nuclear attraction center and charge counts must match"),
    )
    contributions = map(eachindex(center_values)) do center_index
        positive_block = centered_positive_block_builder(center_values[center_index])
        hasproperty(positive_block, :block) || throw(
            ArgumentError("source-box nuclear attraction positive helper must return a block"),
        )
        charge = charge_values[center_index]
        block = Matrix{Float64}((-charge) .* positive_block.block)
        all(isfinite, block) || throw(
            ArgumentError("source-box nuclear attraction block produced non-finite entries"),
        )
        (
            center_index = center_index,
            center = center_values[center_index],
            nuclear_charge = charge,
            sign_charge_scale = -charge,
            gaussian_sum_block = positive_block.block,
            block = block,
            positive_gaussian_sum = positive_block,
            diagnostics = merge(
                positive_block.diagnostics,
                (
                    source = source,
                    physical_operator = :electron_nuclear_attraction,
                    positive_gaussian_sum_component = true,
                    gaussian_sum_path = positive_block.path,
                    nuclear_charge_applied = true,
                    nuclear_attraction_sign_applied = true,
                    nuclear_charge_sign_applied = true,
                    nuclear_charge = charge,
                    sign_charge_scale = -charge,
                    center = center_values[center_index],
                    center_index = center_index,
                    center_contribution = true,
                    center_contributions_preserved = true,
                    counterpoise_center_identity_preserved = true,
                    ecp = false,
                    ecp_terms_implemented = false,
                    electron_electron_terms_implemented = false,
                    mwg_interaction_implemented = false,
                    packet_adoption = false,
                    fixed_block_routing = false,
                    qwhamiltonian_consumes = false,
                    public_default_consumes = false,
                    cr2_science_status_changed = false,
                    output_finite = true,
                ),
            ),
        )
    end
    total_block = zeros(Float64, size(first(contributions).block))
    for contribution in contributions
        size(contribution.block) == size(total_block) || throw(
            DimensionMismatch("source-box nuclear attraction center blocks must have matching sizes"),
        )
        total_block .+= contribution.block
    end
    all(isfinite, total_block) || throw(
        ArgumentError("source-box nuclear attraction total block produced non-finite entries"),
    )
    return (
        path = path,
        physical_operator = :electron_nuclear_attraction,
        blocks_by_center = contributions,
        center_blocks = contributions,
        total_block = total_block,
        block = total_block,
        centers = Tuple(center_values),
        nuclear_charges = Tuple(charge_values),
        diagnostics = (
            source = source,
            private_shadow_only = true,
            source_box_first = true,
            local_gaussian_source_box_terms = true,
            physical_operator = :electron_nuclear_attraction,
            positive_gaussian_sum_component = true,
            positive_gaussian_sum_paths =
                Tuple(contribution.positive_gaussian_sum.path for contribution in contributions),
            nuclear_charge_applied = true,
            nuclear_attraction_sign_applied = true,
            nuclear_charge_sign_applied = true,
            center_contributions_preserved = true,
            counterpoise_center_identity_preserved = true,
            primary_result = :blocks_by_center,
            total_block_is_derived_convenience = true,
            block_field_is_derived_total = true,
            center_count = length(center_values),
            block_shape = size(total_block),
            positive_gaussian_sum_convention = true,
            ecp = false,
            ecp_terms_implemented = false,
            electron_electron_terms_implemented = false,
            mwg_interaction_implemented = false,
            retained_pqs_weight_division_allowed = false,
            ida_weight_division_allowed = false,
            ida_mwg_semantics_changed = false,
            numerical_reference_fallback = false,
            shell_row_algorithm = false,
            packet_adoption = false,
            fixed_block_routing = false,
            qwhamiltonian_consumes = false,
            public_default_consumes = false,
            cr2_science_status_changed = false,
            output_finite = true,
        ),
    )
end

function _product_doside_source_box_nuclear_attraction_by_center(
    left_unit::_CartesianNestedProductStagedByCenterUnit3D,
    right_unit::_CartesianNestedProductStagedByCenterUnit3D,
    axis_layers::NamedTuple{(:x,:y,:z)},
    expansion::CoulombGaussianExpansion;
    centers,
    nuclear_charges,
)
    return _source_box_nuclear_attraction_by_center(
        center -> _product_doside_source_box_centered_local_gaussian_sum_block(
            left_unit,
            right_unit,
            axis_layers,
            expansion;
            center,
        );
        centers,
        nuclear_charges,
        path = :product_doside_source_box_nuclear_attraction_by_center,
        source = :product_doside_source_box_nuclear_attraction_by_center,
    )
end

function _pqs_product_source_box_centered_local_gaussian_sum_block(
    pqs_plan,
    product_unit::_CartesianNestedProductStagedByCenterUnit3D,
    axis_layers::NamedTuple{(:x,:y,:z)},
    expansion::CoulombGaussianExpansion;
    center::NTuple{3,<:Real},
)
    term_data = _product_doside_source_box_centered_local_gaussian_term_data(
        axis_layers,
        expansion;
        center,
    )
    pqs_term_data = merge(
        term_data,
        (
            diagnostics = merge(
                term_data.diagnostics,
                (
                    source =
                        :pqs_product_source_box_centered_local_gaussian_term_data,
                    centered_term_data_helper_reused =
                        :_product_doside_source_box_centered_local_gaussian_term_data,
                ),
            ),
        ),
    )
    explicit_block = _pqs_product_source_box_local_gaussian_sum_block(
        pqs_plan,
        product_unit;
        term_coefficients = pqs_term_data.term_coefficients,
        axis_gaussian_terms = pqs_term_data.axis_gaussian_terms,
    )
    return (
        path = :pqs_product_source_box_centered_local_gaussian_sum,
        block = explicit_block.block,
        explicit_block = explicit_block,
        term_data = pqs_term_data,
        pair_plan = explicit_block.pair_plan,
        one_dimensional_gaussian_factors =
            explicit_block.one_dimensional_gaussian_factors,
        diagnostics = merge(
            explicit_block.diagnostics,
            pqs_term_data.diagnostics,
            (
                source =
                    :pqs_product_source_box_centered_local_gaussian_sum_block,
                source_box_first = true,
                local_gaussian_source_box_terms = true,
                centered_local_gaussian_terms_generated = true,
                centered_term_data_helper_reused =
                    :_product_doside_source_box_centered_local_gaussian_term_data,
                explicit_table_helper_used =
                    :_pqs_product_source_box_local_gaussian_sum_block,
                axis_gaussian_terms_source = :analytic_gaussian_factor_matrices,
                analytic_primitive_backend_required = true,
                analytic_primitive_backend_checked = true,
                pqs_representation = :mode_selected_raw_product_box,
                positive_gaussian_sum_convention = true,
                nuclear_charge_applied = false,
                nuclear_attraction_sign_applied = false,
                numerical_reference_fallback = false,
                shell_projection_used = false,
                lowdin_cleanup_used = false,
                support_coefficient_matrix_used = false,
                shell_row_algorithm = false,
                support_local_oracle_used = false,
                support_local_pqs_oracle_used = false,
                retained_pqs_weight_division_allowed = false,
                ida_weight_division_allowed = false,
                ida_mwg_semantics_changed = false,
                ecp_terms_implemented = false,
                electron_electron_terms_implemented = false,
                mwg_interaction_implemented = false,
                local_gaussian_one_body_implemented = true,
                product_staged_nuclear_execution_changed = false,
                packet_adoption = false,
                fixed_block_routing = false,
                qwhamiltonian_consumes = false,
                public_default_consumes = false,
                cr2_science_status_changed = false,
                output_finite = true,
            ),
        ),
    )
end

function _pqs_product_source_box_nuclear_attraction_by_center(
    pqs_plan,
    product_unit::_CartesianNestedProductStagedByCenterUnit3D,
    axis_layers::NamedTuple{(:x,:y,:z)},
    expansion::CoulombGaussianExpansion;
    centers,
    nuclear_charges,
)
    return _source_box_nuclear_attraction_by_center(
        center -> _pqs_product_source_box_centered_local_gaussian_sum_block(
            pqs_plan,
            product_unit,
            axis_layers,
            expansion;
            center,
        );
        centers,
        nuclear_charges,
        path = :pqs_product_source_box_nuclear_attraction_by_center,
        source = :pqs_product_source_box_nuclear_attraction_by_center,
    )
end

function _same_pqs_raw_product_box_plan(
    left_raw_plan,
    right_raw_plan;
    atol::Real = 1.0e-12,
)
    left_raw_plan.representation == :orthogonal_raw_product_box &&
        right_raw_plan.representation == :orthogonal_raw_product_box || return false
    left_raw_plan.source_mode_dims == right_raw_plan.source_mode_dims || return false
    left_raw_plan.source_mode_count == right_raw_plan.source_mode_count || return false
    left_raw_plan.axis_intervals == right_raw_plan.axis_intervals || return false
    left_raw_plan.source_mode_indices == right_raw_plan.source_mode_indices || return false
    left_raw_plan.source_mode_ordering == right_raw_plan.source_mode_ordering || return false
    left_raw_plan.boundary_selector.mode_indices ==
        right_raw_plan.boundary_selector.mode_indices || return false
    left_raw_plan.boundary_selector.column_indices ==
        right_raw_plan.boundary_selector.column_indices || return false
    left_raw_plan.boundary_selector.selection_rule ==
        right_raw_plan.boundary_selector.selection_rule || return false
    left_raw_plan.boundary_selector.selected_count ==
        right_raw_plan.boundary_selector.selected_count || return false
    return all(
        axis -> isapprox(
            left_raw_plan.axis_local_coefficients[axis],
            right_raw_plan.axis_local_coefficients[axis];
            atol,
            rtol = atol,
        ),
        1:3,
    )
end

function _compatible_cross_pqs_raw_product_box_plans(
    left_raw_plan,
    right_raw_plan,
)
    left_raw_plan.representation == :orthogonal_raw_product_box &&
        right_raw_plan.representation == :orthogonal_raw_product_box || return false
    left_raw_plan.source_mode_dims == right_raw_plan.source_mode_dims || return false
    left_raw_plan.source_mode_count == right_raw_plan.source_mode_count || return false
    left_raw_plan.source_mode_indices == right_raw_plan.source_mode_indices || return false
    left_raw_plan.source_mode_ordering == right_raw_plan.source_mode_ordering || return false
    left_raw_plan.boundary_selector.mode_indices ==
        right_raw_plan.boundary_selector.mode_indices || return false
    left_raw_plan.boundary_selector.column_indices ==
        right_raw_plan.boundary_selector.column_indices || return false
    left_raw_plan.boundary_selector.selection_rule ==
        right_raw_plan.boundary_selector.selection_rule || return false
    left_raw_plan.boundary_selector.selected_count ==
        right_raw_plan.boundary_selector.selected_count || return false
    return true
end

function _pqs_pqs_source_box_project_axis_terms(
    operator_terms::AbstractArray{<:Real,3},
    left_pqs_plan,
    right_pqs_plan,
    axis::Int;
    label::AbstractString,
)
    left_raw_plan = _pqs_raw_product_box_plan_view(left_pqs_plan)
    right_raw_plan = _pqs_raw_product_box_plan_view(right_pqs_plan)
    nterms = size(operator_terms, 1)
    nterms > 0 || throw(
        ArgumentError("$(label) requires at least one term"),
    )
    left_interval = left_raw_plan.axis_intervals[axis]
    right_interval = right_raw_plan.axis_intervals[axis]
    first(left_interval) >= 1 && last(left_interval) <= size(operator_terms, 2) ||
        throw(
            ArgumentError("$(label) left interval exceeds axis $(axis) term table"),
        )
    first(right_interval) >= 1 && last(right_interval) <= size(operator_terms, 3) ||
        throw(
            ArgumentError("$(label) right interval exceeds axis $(axis) term table"),
        )
    left_coefficients =
        Matrix{Float64}(left_raw_plan.axis_local_coefficients[axis])
    right_coefficients =
        Matrix{Float64}(right_raw_plan.axis_local_coefficients[axis])
    size(left_coefficients, 1) == length(left_interval) || throw(
        DimensionMismatch("$(label) left axis coefficients must match interval length"),
    )
    size(right_coefficients, 1) == length(right_interval) || throw(
        DimensionMismatch("$(label) right axis coefficients must match interval length"),
    )
    projected = Array{Float64,3}(
        undef,
        nterms,
        size(left_coefficients, 2),
        size(right_coefficients, 2),
    )
    @inbounds for term in 1:nterms
        term_matrix = @view operator_terms[term, left_interval, right_interval]
        projected[term, :, :] .=
            transpose(left_coefficients) * term_matrix * right_coefficients
    end
    return projected
end

function _pqs_pqs_source_box_project_local_gaussian_axis_terms(
    operator_terms::AbstractArray{<:Real,3},
    left_pqs_plan,
    right_pqs_plan,
    axis::Int,
)
    return _pqs_pqs_source_box_project_axis_terms(
        operator_terms,
        left_pqs_plan,
        right_pqs_plan,
        axis;
        label = "PQS/PQS source-box local-Gaussian terms",
    )
end

function _pqs_pqs_source_box_local_gaussian_axis_factors(
    left_pqs_plan,
    right_pqs_plan,
    axis_gaussian_terms::NamedTuple{(:x,:y,:z)},
)
    return ntuple(axis -> begin
        terms = getproperty(axis_gaussian_terms, (:x, :y, :z)[axis])
        _pqs_pqs_source_box_project_local_gaussian_axis_terms(
            terms,
            left_pqs_plan,
            right_pqs_plan,
            axis,
        )
    end, 3)
end

function _pqs_pqs_source_box_density_density_axis_factors(
    left_pqs_plan,
    right_pqs_plan,
    axis_pair_factor_terms::NamedTuple{(:x,:y,:z)},
)
    return ntuple(axis -> begin
        terms = getproperty(axis_pair_factor_terms, (:x, :y, :z)[axis])
        _pqs_pqs_source_box_project_axis_terms(
            terms,
            left_pqs_plan,
            right_pqs_plan,
            axis;
            label = "PQS/PQS source-box density-density pair-factor terms",
        )
    end, 3)
end

function _pqs_pqs_source_box_density_density_weight_views(
    left_raw_plan,
    right_raw_plan,
    axis_weights::NamedTuple{(:x,:y,:z)},
)
    return (
        left = ntuple(axis -> _source_box_axis_positive_weights(
            getproperty(axis_weights, (:x, :y, :z)[axis]),
            left_raw_plan.axis_intervals[axis];
            axis_name = (:x, :y, :z)[axis],
            side = :left_pqs,
        ), 3),
        right = ntuple(axis -> _source_box_axis_positive_weights(
            getproperty(axis_weights, (:x, :y, :z)[axis]),
            right_raw_plan.axis_intervals[axis];
            axis_name = (:x, :y, :z)[axis],
            side = :right_pqs,
        ), 3),
    )
end

function _pqs_pqs_source_box_density_density_interaction_block(
    left_pqs_plan,
    right_pqs_plan;
    term_coefficients::AbstractVector{<:Real},
    axis_pair_factor_terms::NamedTuple{(:x,:y,:z)},
    axis_weights::NamedTuple{(:x,:y,:z)},
    pair_factor_normalization::Symbol = :density_normalized,
)
    pair_factor_normalization == :density_normalized || throw(
        ArgumentError(
            "PQS/PQS source-box density-density fixture currently requires density-normalized pair factors",
        ),
    )
    left_raw_plan = _pqs_raw_product_box_plan_view(left_pqs_plan)
    right_raw_plan = _pqs_raw_product_box_plan_view(right_pqs_plan)
    left_raw_plan.representation == :orthogonal_raw_product_box || throw(
        ArgumentError("PQS/PQS source-box density-density block requires a raw product-box left PQS plan"),
    )
    right_raw_plan.representation == :orthogonal_raw_product_box || throw(
        ArgumentError("PQS/PQS source-box density-density block requires a raw product-box right PQS plan"),
    )
    _compatible_cross_pqs_raw_product_box_plans(left_raw_plan, right_raw_plan) ||
        throw(
            ArgumentError("PQS/PQS source-box density-density block requires equal source-mode dimensions, ordering, and boundary selectors"),
        )
    nterms = length(term_coefficients)
    nterms > 0 || throw(
        ArgumentError("PQS/PQS source-box density-density block requires at least one term coefficient"),
    )
    for axis_name in (:x, :y, :z)
        size(getproperty(axis_pair_factor_terms, axis_name), 1) == nterms || throw(
            ArgumentError("PQS/PQS source-box density-density pair-factor term count mismatch on $(axis_name)"),
        )
    end
    source_weights = _pqs_pqs_source_box_density_density_weight_views(
        left_raw_plan,
        right_raw_plan,
        axis_weights,
    )
    projected_terms = _pqs_pqs_source_box_density_density_axis_factors(
        left_raw_plan,
        right_raw_plan,
        axis_pair_factor_terms,
    )
    coeffs = Float64[Float64(value) for value in term_coefficients]
    left_modes = left_raw_plan.boundary_selector.mode_indices
    right_modes = right_raw_plan.boundary_selector.mode_indices
    length(left_modes) == left_raw_plan.boundary_selector.selected_count || throw(
        ArgumentError("PQS/PQS source-box density-density left boundary mode count must match selected count"),
    )
    length(right_modes) == right_raw_plan.boundary_selector.selected_count || throw(
        ArgumentError("PQS/PQS source-box density-density right boundary mode count must match selected count"),
    )
    block = zeros(Float64, length(left_modes), length(right_modes))
    projected_x, projected_y, projected_z = projected_terms
    @inbounds for col in eachindex(right_modes)
        rx, ry, rz = right_modes[col]
        for row in eachindex(left_modes)
            lx, ly, lz = left_modes[row]
            value = 0.0
            @simd for term in 1:nterms
                value +=
                    coeffs[term] *
                    projected_x[term, lx, rx] *
                    projected_y[term, ly, ry] *
                    projected_z[term, lz, rz]
            end
            block[row, col] = value
        end
    end
    all(isfinite, block) || throw(
        ArgumentError("PQS/PQS source-box density-density block produced non-finite entries"),
    )
    same_raw_product_box_plan =
        _same_pqs_raw_product_box_plan(left_raw_plan, right_raw_plan)
    pair_plan = (
        pair_kind = :pqs_pqs_source_box_density_density_pair,
        pair_family = :pqs_pqs,
        left_source_family = :mode_selected_raw_product_box,
        right_source_family = :mode_selected_raw_product_box,
        left_retained_rule_kind = left_raw_plan.retained_rule_kind,
        right_retained_rule_kind = right_raw_plan.retained_rule_kind,
        left_source_dimensions = left_raw_plan.source_mode_dims,
        right_source_dimensions = right_raw_plan.source_mode_dims,
        left_source_dimension = left_raw_plan.source_mode_count,
        right_source_dimension = right_raw_plan.source_mode_count,
        left_retained_count = left_raw_plan.boundary_selector.selected_count,
        right_retained_count = right_raw_plan.boundary_selector.selected_count,
        axis_intervals = (
            left = left_raw_plan.axis_intervals,
            right = right_raw_plan.axis_intervals,
        ),
        left_boundary_mode_selector = left_raw_plan.boundary_selector,
        right_boundary_mode_selector = right_raw_plan.boundary_selector,
        source_weights = source_weights,
        supported_terms = (:pair_sum,),
        diagnostics = (
            source = :pqs_pqs_source_box_density_density_pair_plan,
            private_shadow_only = true,
            source_box_pair_plan = true,
            pair_family = :pqs_pqs,
            interaction_operator = :electron_electron_density_density,
            output_representation = :two_index_density_density,
            four_index_galerkin_tensor = false,
            pqs_representation = :mode_selected_raw_product_box,
            left_pqs_boundary_mode_selection_used = true,
            right_pqs_boundary_mode_selection_used = true,
            same_raw_product_box_plan = same_raw_product_box_plan,
            cross_pqs_inputs_supported = !same_raw_product_box_plan,
            shared_raw_product_box_plan_used =
                left_raw_plan.shared_raw_product_box_plan_used &&
                right_raw_plan.shared_raw_product_box_plan_used,
            source_mode_ordering = left_raw_plan.source_mode_ordering,
            operator_factor_source = :explicit_density_normalized_pair_factor_terms,
            input_pair_factor_data = :caller_supplied_explicit_data,
            input_pair_factor_data_pgdg_checked = false,
            raw_source_weights_available = true,
            density_normalized_pair_factors = true,
            raw_weighted_pair_factors = false,
            pair_factor_normalization = pair_factor_normalization,
            source_weight_division_owner = :caller_supplied_density_normalized_pair_factors,
            source_weight_division_applied_by_helper = false,
            raw_product_box_operators_use_1d_factors = true,
            shell_projection_used = false,
            lowdin_cleanup_used = false,
            support_coefficient_matrix_used = false,
            shell_row_algorithm = false,
            support_local_oracle_used = false,
            support_local_pqs_oracle_used = false,
            retained_pqs_weights_used = false,
            retained_pqs_weights_positive_checked = false,
            retained_weight_division_allowed = false,
            retained_pqs_weight_division_allowed = false,
            ida_weight_division_allowed = false,
            ida_mwg_semantics_changed = false,
            mwg_ida_semantics_changed = false,
            retained_weight_semantics_changed = false,
            packet_adoption = false,
            fixed_block_routing = false,
            qwhamiltonian_consumes = false,
            public_default_consumes = false,
            cr2_science_status_changed = false,
            numerical_reference_fallback = false,
            ecp_terms_implemented = false,
            local_gaussian_one_body_implemented = false,
            mwg_interaction_implemented = false,
            generic_retained_unit_framework = false,
        ),
    )
    return (
        path = :pqs_pqs_source_box_density_density_interaction,
        interaction_operator = :electron_electron_density_density,
        block = block,
        pair_plan = pair_plan,
        one_dimensional_pair_factors = (
            x = projected_x,
            y = projected_y,
            z = projected_z,
        ),
        diagnostics = merge(
            pair_plan.diagnostics,
            (
                source = :pqs_pqs_source_box_density_density_interaction_block,
                source_box_first = true,
                output_finite = true,
                electron_electron_terms_implemented = true,
                term_count = nterms,
                axis_term_dimensions = (
                    x = size(axis_pair_factor_terms.x),
                    y = size(axis_pair_factor_terms.y),
                    z = size(axis_pair_factor_terms.z),
                ),
                projected_axis_term_dimensions = (
                    x = size(projected_x),
                    y = size(projected_y),
                    z = size(projected_z),
                ),
            ),
        ),
    )
end

function _pqs_pqs_source_box_raw_weighted_density_density_interaction_block(
    left_pqs_plan,
    right_pqs_plan;
    term_coefficients::AbstractVector{<:Real},
    raw_axis_pair_factor_terms::NamedTuple{(:x,:y,:z)},
    axis_weights::NamedTuple{(:x,:y,:z)},
)
    density_normalized_terms = _source_box_density_normalized_axis_pair_terms(
        raw_axis_pair_factor_terms,
        axis_weights,
    )
    density_normalized_core =
        _pqs_pqs_source_box_density_density_interaction_block(
            left_pqs_plan,
            right_pqs_plan;
            term_coefficients,
            axis_pair_factor_terms = density_normalized_terms,
            axis_weights,
            pair_factor_normalization = :density_normalized,
        )
    return (
        path = :pqs_pqs_source_box_raw_weighted_density_density_interaction,
        interaction_operator = :electron_electron_density_density,
        block = density_normalized_core.block,
        density_normalized_core = density_normalized_core,
        normalized_axis_pair_factor_terms = density_normalized_terms,
        raw_axis_pair_factor_terms = raw_axis_pair_factor_terms,
        pair_plan = density_normalized_core.pair_plan,
        diagnostics = merge(
            density_normalized_core.diagnostics,
            (
                source = :pqs_pqs_source_box_raw_weighted_density_density_interaction_block,
                path = :pqs_pqs_source_box_raw_weighted_density_density_interaction,
                pair_factor_normalization = :raw_weighted,
                raw_weighted_pair_factors = true,
                density_normalized_pair_factors = false,
                density_normalized_pair_factors_generated = true,
                source_weight_division_owner = :source_box_raw_weights,
                source_weight_division_applied_by_helper = true,
                source_weight_division_shape = :axis_pair_weight_outer,
                density_normalized_core_helper =
                    :_pqs_pqs_source_box_density_density_interaction_block,
                retained_pqs_weights_used = false,
                retained_weight_division_allowed = false,
                retained_pqs_weight_division_allowed = false,
                ida_weight_division_allowed = false,
                ida_mwg_semantics_changed = false,
                mwg_ida_semantics_changed = false,
                packet_adoption = false,
                qwhamiltonian_consumes = false,
                numerical_reference_fallback = false,
            ),
        ),
    )
end

function _pqs_pqs_source_box_local_gaussian_sum_block(
    left_pqs_plan,
    right_pqs_plan;
    term_coefficients::AbstractVector{<:Real},
    axis_gaussian_terms::NamedTuple{(:x,:y,:z)},
)
    left_raw_plan = _pqs_raw_product_box_plan_view(left_pqs_plan)
    right_raw_plan = _pqs_raw_product_box_plan_view(right_pqs_plan)
    left_raw_plan.representation == :orthogonal_raw_product_box || throw(
        ArgumentError("PQS/PQS source-box local-Gaussian block requires a raw product-box left PQS plan"),
    )
    right_raw_plan.representation == :orthogonal_raw_product_box || throw(
        ArgumentError("PQS/PQS source-box local-Gaussian block requires a raw product-box right PQS plan"),
    )
    _compatible_cross_pqs_raw_product_box_plans(left_raw_plan, right_raw_plan) ||
        throw(
            ArgumentError("PQS/PQS source-box local-Gaussian block requires equal source-mode dimensions, ordering, and boundary selectors"),
        )
    nterms = length(term_coefficients)
    nterms > 0 || throw(
        ArgumentError("PQS/PQS source-box local-Gaussian block requires at least one term coefficient"),
    )
    for axis_name in (:x, :y, :z)
        size(getproperty(axis_gaussian_terms, axis_name), 1) == nterms || throw(
            ArgumentError("PQS/PQS source-box local-Gaussian axis term count mismatch on $(axis_name)"),
        )
    end
    projected_terms =
        _pqs_pqs_source_box_local_gaussian_axis_factors(
            left_raw_plan,
            right_raw_plan,
            axis_gaussian_terms,
        )
    coeffs = Float64[Float64(value) for value in term_coefficients]
    left_modes = left_raw_plan.boundary_selector.mode_indices
    right_modes = right_raw_plan.boundary_selector.mode_indices
    length(left_modes) == left_raw_plan.boundary_selector.selected_count || throw(
        ArgumentError("PQS/PQS source-box local-Gaussian left boundary mode count must match selected count"),
    )
    length(right_modes) == right_raw_plan.boundary_selector.selected_count || throw(
        ArgumentError("PQS/PQS source-box local-Gaussian right boundary mode count must match selected count"),
    )
    block = zeros(Float64, length(left_modes), length(right_modes))
    projected_x, projected_y, projected_z = projected_terms
    @inbounds for col in eachindex(right_modes)
        rx, ry, rz = right_modes[col]
        for row in eachindex(left_modes)
            lx, ly, lz = left_modes[row]
            value = 0.0
            @simd for term in 1:nterms
                value +=
                    coeffs[term] *
                    projected_x[term, lx, rx] *
                    projected_y[term, ly, ry] *
                    projected_z[term, lz, rz]
            end
            block[row, col] = value
        end
    end
    all(isfinite, block) || throw(
        ArgumentError("PQS/PQS source-box local-Gaussian block produced non-finite entries"),
    )
    same_raw_product_box_plan =
        _same_pqs_raw_product_box_plan(left_raw_plan, right_raw_plan)
    pair_plan = (
        pair_kind = :pqs_pqs_source_box_local_gaussian_pair,
        left_source_family = :mode_selected_raw_product_box,
        right_source_family = :mode_selected_raw_product_box,
        left_retained_rule_kind = left_raw_plan.retained_rule_kind,
        right_retained_rule_kind = right_raw_plan.retained_rule_kind,
        left_source_dimensions = left_raw_plan.source_mode_dims,
        right_source_dimensions = right_raw_plan.source_mode_dims,
        left_source_dimension = left_raw_plan.source_mode_count,
        right_source_dimension = right_raw_plan.source_mode_count,
        left_retained_count = left_raw_plan.boundary_selector.selected_count,
        right_retained_count = right_raw_plan.boundary_selector.selected_count,
        axis_intervals = (
            left = left_raw_plan.axis_intervals,
            right = right_raw_plan.axis_intervals,
        ),
        left_boundary_mode_selector = left_raw_plan.boundary_selector,
        right_boundary_mode_selector = right_raw_plan.boundary_selector,
        supported_terms = (:gaussian_sum,),
        diagnostics = (
            source = :pqs_pqs_source_box_local_gaussian_pair_plan,
            private_shadow_only = true,
            source_box_pair_plan = true,
            pqs_representation = :mode_selected_raw_product_box,
            left_pqs_boundary_mode_selection_used = true,
            right_pqs_boundary_mode_selection_used = true,
            same_raw_product_box_plan = same_raw_product_box_plan,
            cross_pqs_inputs_supported = !same_raw_product_box_plan,
            shared_raw_product_box_plan_used =
                left_raw_plan.shared_raw_product_box_plan_used &&
                right_raw_plan.shared_raw_product_box_plan_used,
            source_mode_ordering = left_raw_plan.source_mode_ordering,
            operator_factor_source = :explicit_local_gaussian_axis_terms,
            input_local_gaussian_data = :caller_supplied_explicit_data,
            input_local_gaussian_data_pgdg_checked = false,
            pgdg_analytic_operator_provenance_claimed = false,
            raw_product_box_numerical_reference_fallback =
                hasproperty(left_raw_plan.diagnostics, :raw_product_box_numerical_reference_fallback) ?
                left_raw_plan.diagnostics.raw_product_box_numerical_reference_fallback :
                nothing,
            numerical_reference_fallback = false,
            raw_product_box_operators_use_1d_factors = true,
            shell_projection_used = false,
            lowdin_cleanup_used = false,
            support_coefficient_matrix_used = false,
            shell_row_algorithm = false,
            support_local_oracle_used = false,
            support_local_pqs_oracle_used = false,
            retained_pqs_weights_used = false,
            retained_pqs_weights_positive_checked = false,
            retained_pqs_weight_division_allowed = false,
            ida_weight_division_allowed = false,
            ida_mwg_semantics_changed = false,
            retained_weight_semantics_changed = false,
            packet_adoption = false,
            fixed_block_routing = false,
            qwhamiltonian_consumes = false,
            public_default_consumes = false,
            cr2_science_status_changed = false,
            ecp_terms_implemented = false,
            electron_electron_terms_implemented = false,
            mwg_interaction_implemented = false,
            generic_retained_unit_framework = false,
        ),
    )
    return (
        path = :pqs_pqs_source_box_local_gaussian_sum,
        block = block,
        pair_plan = pair_plan,
        one_dimensional_gaussian_factors = (
            x = projected_x,
            y = projected_y,
            z = projected_z,
        ),
        diagnostics = merge(
            pair_plan.diagnostics,
            (
                source = :pqs_pqs_source_box_local_gaussian_sum_block,
                source_box_first = true,
                local_gaussian_source_box_terms = true,
                positive_gaussian_sum_convention = true,
                nuclear_charge_applied = false,
                nuclear_attraction_sign_applied = false,
                term_coefficients_source = :caller_supplied_explicit_data,
                axis_gaussian_terms_source = :caller_supplied_explicit_data,
                input_local_gaussian_data_pgdg_checked = false,
                pgdg_analytic_operator_provenance_claimed = false,
                numerical_reference_fallback = false,
                term_count = nterms,
                axis_term_dimensions = (
                    x = size(axis_gaussian_terms.x),
                    y = size(axis_gaussian_terms.y),
                    z = size(axis_gaussian_terms.z),
                ),
                projected_axis_term_dimensions = (
                    x = size(projected_x),
                    y = size(projected_y),
                    z = size(projected_z),
                ),
                operator_factor_source = :explicit_local_gaussian_axis_terms,
                dense_raw_source_box_pair_matrix_materialized = false,
                dense_raw_pair_storage_avoided = true,
                retained_block_assembled_directly_from_1d_factors = true,
                source_box_pair_storage_scaling =
                    :one_dimensional_factors_plus_retained_block,
                pqs_representation = :mode_selected_raw_product_box,
                shell_projection_used = false,
                lowdin_cleanup_used = false,
                support_coefficient_matrix_used = false,
                shell_row_algorithm = false,
                support_local_oracle_used = false,
                support_local_pqs_oracle_used = false,
                retained_pqs_weight_division_allowed = false,
                ida_weight_division_allowed = false,
                ida_mwg_semantics_changed = false,
                ecp_terms_implemented = false,
                electron_electron_terms_implemented = false,
                mwg_interaction_implemented = false,
                local_gaussian_one_body_implemented = true,
                packet_adoption = false,
                fixed_block_routing = false,
                qwhamiltonian_consumes = false,
                public_default_consumes = false,
                cr2_science_status_changed = false,
                output_finite = true,
            ),
        ),
    )
end

function _pqs_pqs_source_box_centered_local_gaussian_sum_block(
    left_pqs_plan,
    right_pqs_plan,
    axis_layers::NamedTuple{(:x,:y,:z)},
    expansion::CoulombGaussianExpansion;
    center::NTuple{3,<:Real},
)
    term_data = _product_doside_source_box_centered_local_gaussian_term_data(
        axis_layers,
        expansion;
        center,
    )
    pqs_term_data = merge(
        term_data,
        (
            diagnostics = merge(
                term_data.diagnostics,
                (
                    source =
                        :pqs_pqs_source_box_centered_local_gaussian_term_data,
                    centered_term_data_helper_reused =
                        :_product_doside_source_box_centered_local_gaussian_term_data,
                ),
            ),
        ),
    )
    explicit_block = _pqs_pqs_source_box_local_gaussian_sum_block(
        left_pqs_plan,
        right_pqs_plan;
        term_coefficients = pqs_term_data.term_coefficients,
        axis_gaussian_terms = pqs_term_data.axis_gaussian_terms,
    )
    return (
        path = :pqs_pqs_source_box_centered_local_gaussian_sum,
        block = explicit_block.block,
        explicit_block = explicit_block,
        term_data = pqs_term_data,
        pair_plan = explicit_block.pair_plan,
        one_dimensional_gaussian_factors =
            explicit_block.one_dimensional_gaussian_factors,
        diagnostics = merge(
            explicit_block.diagnostics,
            pqs_term_data.diagnostics,
            (
                source = :pqs_pqs_source_box_centered_local_gaussian_sum_block,
                source_box_first = true,
                local_gaussian_source_box_terms = true,
                centered_local_gaussian_terms_generated = true,
                centered_term_data_helper_reused =
                    :_product_doside_source_box_centered_local_gaussian_term_data,
                explicit_table_helper_used =
                    :_pqs_pqs_source_box_local_gaussian_sum_block,
                axis_gaussian_terms_source = :analytic_gaussian_factor_matrices,
                analytic_primitive_backend_required = true,
                analytic_primitive_backend_checked = true,
                pqs_representation = :mode_selected_raw_product_box,
                positive_gaussian_sum_convention = true,
                nuclear_charge_applied = false,
                nuclear_attraction_sign_applied = false,
                numerical_reference_fallback = false,
                shell_projection_used = false,
                lowdin_cleanup_used = false,
                support_coefficient_matrix_used = false,
                shell_row_algorithm = false,
                support_local_oracle_used = false,
                support_local_pqs_oracle_used = false,
                retained_pqs_weight_division_allowed = false,
                ida_weight_division_allowed = false,
                ida_mwg_semantics_changed = false,
                ecp_terms_implemented = false,
                electron_electron_terms_implemented = false,
                mwg_interaction_implemented = false,
                local_gaussian_one_body_implemented = true,
                packet_adoption = false,
                fixed_block_routing = false,
                qwhamiltonian_consumes = false,
                public_default_consumes = false,
                cr2_science_status_changed = false,
                output_finite = true,
            ),
        ),
    )
end

function _pqs_pqs_source_box_nuclear_attraction_by_center(
    left_pqs_plan,
    right_pqs_plan,
    axis_layers::NamedTuple{(:x,:y,:z)},
    expansion::CoulombGaussianExpansion;
    centers,
    nuclear_charges,
)
    return _source_box_nuclear_attraction_by_center(
        center -> _pqs_pqs_source_box_centered_local_gaussian_sum_block(
            left_pqs_plan,
            right_pqs_plan,
            axis_layers,
            expansion;
            center,
        );
        centers,
        nuclear_charges,
        path = :pqs_pqs_source_box_nuclear_attraction_by_center,
        source = :pqs_pqs_source_box_nuclear_attraction_by_center,
    )
end

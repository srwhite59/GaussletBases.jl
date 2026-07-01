function _staged_axis_interval(axis)
    axis.kind == :fixed && return axis.fixed_index:axis.fixed_index
    axis.kind == :active && return axis.interval
    throw(ArgumentError("unknown product-staged axis kind $(axis.kind)"))
end

function _staged_axis_count(axis)
    axis.kind == :fixed && return 1
    axis.kind == :active && return size(axis.coefficient_matrix, 2)
    throw(ArgumentError("unknown product-staged axis kind $(axis.kind)"))
end

function _validate_staged_axis(axis, parent_count::Int)
    if axis.kind == :fixed
        index = something(axis.fixed_index, 0)
        1 <= index <= parent_count || throw(
            ArgumentError("product-staged metric axis fixed index lies outside parent dimensions"),
        )
        size(axis.coefficient_matrix) == (1, 1) || throw(
            ArgumentError("product-staged metric fixed-axis coefficient map must be 1x1"),
        )
    elseif axis.kind == :active
        interval = axis.interval
        !isnothing(interval) || throw(
            ArgumentError("product-staged metric active axis is missing its interval"),
        )
        first(interval) >= 1 && last(interval) <= parent_count || throw(
            ArgumentError("product-staged metric active-axis interval exceeds parent dimensions"),
        )
        size(axis.coefficient_matrix, 1) == length(interval) || throw(
            ArgumentError("product-staged metric active-axis coefficient rows must match interval length"),
        )
        size(axis.coefficient_matrix, 2) >= 1 || throw(
            ArgumentError("product-staged metric active-axis coefficient map must retain at least one column"),
        )
    else
        throw(ArgumentError("unknown product-staged metric axis kind $(axis.kind)"))
    end
    return nothing
end

function _project_staged_axis_matrix(
    left_axis,
    right_axis,
    matrix::AbstractMatrix{<:Real},
)
    parent_count = size(matrix, 1)
    size(matrix, 2) == parent_count || throw(
        ArgumentError("product-staged metric axis matrix must be square"),
    )
    _validate_staged_axis(left_axis, parent_count)
    _validate_staged_axis(right_axis, parent_count)
    left_interval = _staged_axis_interval(left_axis)
    right_interval = _staged_axis_interval(right_axis)
    local_matrix = Matrix{Float64}(matrix[left_interval, right_interval])
    return Matrix{Float64}(
        transpose(left_axis.coefficient_matrix) * local_matrix * right_axis.coefficient_matrix,
    )
end

function _project_staged_axis_vector(axis, values::AbstractVector{<:Real})
    _validate_staged_axis(axis, length(values))
    interval = _staged_axis_interval(axis)
    local_values = Float64[Float64(value) for value in values[interval]]
    return vec(transpose(axis.coefficient_matrix) * local_values)
end

function _product_doside_low_order_axis_matrix_kind(term::Symbol, axis::Int)
    term == :overlap && return :overlap
    term == :position_x && return axis == 1 ? :position : :overlap
    term == :position_y && return axis == 2 ? :position : :overlap
    term == :position_z && return axis == 3 ? :position : :overlap
    term == :x2_x && return axis == 1 ? :x2 : :overlap
    term == :x2_y && return axis == 2 ? :x2 : :overlap
    term == :x2_z && return axis == 3 ? :x2 : :overlap
    throw(
        ArgumentError(
            "product/doside retained low-order block supports only :overlap, :position_x/y/z, and :x2_x/y/z",
        ),
    )
end

function _product_doside_axis_metric_matrix(metrics::NamedTuple{(:x,:y,:z)}, axis::Int, kind::Symbol)
    axis_name = (:x, :y, :z)[axis]
    return getproperty(getproperty(metrics, axis_name), kind)
end

function _cartesian_source_box_metric_sources(metrics::NamedTuple{(:x,:y,:z)})
    return ntuple(axis -> begin
        axis_data = getproperty(metrics, (:x, :y, :z)[axis])
        hasproperty(axis_data, :source) ? axis_data.source : :unspecified
    end, 3)
end

function _require_product_doside_retained_block_unit(
    unit::_CartesianNestedProductStagedByCenterUnit3D;
    side::Symbol,
)
    unit.kind == :product_doside || throw(
        ArgumentError("product/doside retained low-order block requires a product_doside $(side) unit"),
    )
    length(unit.axis_function_indices) == length(unit.column_range) || throw(
        ArgumentError("product/doside retained low-order block $(side) unit axis metadata does not match its column range"),
    )
    return nothing
end

function _product_doside_retained_unit_plan(
    product_unit::_CartesianNestedProductStagedByCenterUnit3D,
)
    _require_product_doside_retained_block_unit(product_unit; side = :product)
    source_axis_intervals =
        ntuple(axis -> _staged_axis_interval(product_unit.axes[axis]), 3)
    source_axis_lengths = ntuple(axis -> length(source_axis_intervals[axis]), 3)
    retained_axis_counts =
        ntuple(axis -> _staged_axis_count(product_unit.axes[axis]), 3)
    retained_count = length(product_unit.column_range)
    return (
        object_kind = :product_doside_retained_unit_plan,
        kind = product_unit.kind,
        retained_rule_kind = :product_doside,
        source_axis_intervals = source_axis_intervals,
        axis_intervals = source_axis_intervals,
        source_axis_lengths = source_axis_lengths,
        source_dimension = prod(source_axis_lengths),
        retained_axis_counts = retained_axis_counts,
        column_range = product_unit.column_range,
        retained_count = retained_count,
        axes = product_unit.axes,
        axis_coefficient_matrices =
            ntuple(axis -> product_unit.axes[axis].coefficient_matrix, 3),
        axis_function_indices = product_unit.axis_function_indices,
        support_indices = product_unit.support_indices,
        support_states = product_unit.support_states,
        coefficient_matrix = product_unit.coefficient_matrix,
        provenance = product_unit.provenance,
        diagnostics = (
            source = :product_doside_retained_unit_plan,
            private_adapter = true,
            metadata_only = true,
            coefficients_rebuilt = false,
            block_math_changed = false,
            retained_rule_kind = :product_doside,
            packet_adoption = false,
            fixed_block_routing = false,
            qwhamiltonian_consumes = false,
            public_default_consumes = false,
            ida_weight_semantics_changed = false,
            retained_weight_division_allowed = false,
            generic_retained_unit_framework = false,
        ),
    )
end

function _product_doside_source_axis_center_metadata(
    construction::_BondAlignedDiatomicHighOrderRecipeSourceConstruction3D,
    axes::Tuple,
)
    source_axis_intervals = ntuple(axis -> _staged_axis_interval(axes[axis]), 3)
    source_mode_dims = ntuple(axis -> _staged_axis_count(axes[axis]), 3)
    raw_product_box_plan = _cartesian_raw_product_box_plan(
        construction.axis_bundles,
        source_axis_intervals,
        source_mode_dims;
        enforce_symmetric_odd = false,
    )
    return (
        source_axis_center_vectors = ntuple(
            axis -> Float64.(raw_product_box_plan.axis_transform_plan.axes[axis].localized_centers),
            3,
        ),
        source_center_convention = :comx_construction,
        source_center_status = :native_representative,
        source_center_source = :cartesian_raw_product_box_plan_axis_transform_localized_centers,
        raw_product_box_plan_object_kind = raw_product_box_plan.object_kind,
        source_axis_intervals = source_axis_intervals,
        source_mode_dims = source_mode_dims,
    )
end

function _product_doside_retained_low_order_block(
    left_unit::_CartesianNestedProductStagedByCenterUnit3D,
    right_unit::_CartesianNestedProductStagedByCenterUnit3D,
    metrics::NamedTuple{(:x,:y,:z)};
    term::Symbol,
)
    term in (:overlap, :position_x, :position_y, :position_z, :x2_x, :x2_y, :x2_z) || throw(
        ArgumentError("product/doside retained low-order block supports only :overlap, :position_x/y/z, and :x2_x/y/z"),
    )
    _require_product_doside_retained_block_unit(left_unit; side = :left)
    _require_product_doside_retained_block_unit(right_unit; side = :right)
    projected_axis_matrices = ntuple(
        axis -> _project_staged_axis_matrix(
            left_unit.axes[axis],
            right_unit.axes[axis],
            _product_doside_axis_metric_matrix(
                metrics,
                axis,
                _product_doside_low_order_axis_matrix_kind(term, axis),
            ),
        ),
        3,
    )
    block = zeros(Float64, length(left_unit.column_range), length(right_unit.column_range))
    @inbounds for local_col in eachindex(right_unit.axis_function_indices)
        xj, yj, zj = right_unit.axis_function_indices[local_col]
        for local_row in eachindex(left_unit.axis_function_indices)
            xi, yi, zi = left_unit.axis_function_indices[local_row]
            block[local_row, local_col] =
                projected_axis_matrices[1][xi, xj] *
                projected_axis_matrices[2][yi, yj] *
                projected_axis_matrices[3][zi, zj]
        end
    end
    return block
end

function _product_doside_retained_separable_sum_block(
    left_unit::_CartesianNestedProductStagedByCenterUnit3D,
    right_unit::_CartesianNestedProductStagedByCenterUnit3D,
    axis_factor_terms,
)
    _require_product_doside_retained_block_unit(left_unit; side = :left)
    _require_product_doside_retained_block_unit(right_unit; side = :right)
    !isempty(axis_factor_terms) || throw(
        ArgumentError("product/doside retained separable-sum block requires at least one axis-factor triple"),
    )
    block = zeros(Float64, length(left_unit.column_range), length(right_unit.column_range))
    for term in axis_factor_terms
        length(term) == 3 || throw(
            ArgumentError("product/doside retained separable-sum block terms must be axis-factor triples"),
        )
        projected_axis_matrices = ntuple(
            axis -> _project_staged_axis_matrix(
                left_unit.axes[axis],
                right_unit.axes[axis],
                term[axis],
            ),
            3,
        )
        @inbounds for local_col in eachindex(right_unit.axis_function_indices)
            xj, yj, zj = right_unit.axis_function_indices[local_col]
            for local_row in eachindex(left_unit.axis_function_indices)
                xi, yi, zi = left_unit.axis_function_indices[local_row]
                block[local_row, local_col] +=
                    projected_axis_matrices[1][xi, xj] *
                    projected_axis_matrices[2][yi, yj] *
                    projected_axis_matrices[3][zi, zj]
            end
        end
    end
    return block
end

function _product_doside_axis_operator_matrix(axis_ops::NamedTuple{(:x,:y,:z)}, axis::Int, kind::Symbol)
    axis_name = (:x, :y, :z)[axis]
    axis_data = getproperty(axis_ops, axis_name)
    hasproperty(axis_data, kind) || throw(
        ArgumentError("product/doside retained kinetic block axis $(axis_name) is missing $(kind) matrix"),
    )
    return getproperty(axis_data, kind)
end

function _product_doside_kinetic_axis_factor_terms(
    axis_ops::NamedTuple{(:x,:y,:z)},
)
    return (
        (
            _product_doside_axis_operator_matrix(axis_ops, 1, :kinetic),
            _product_doside_axis_operator_matrix(axis_ops, 2, :overlap),
            _product_doside_axis_operator_matrix(axis_ops, 3, :overlap),
        ),
        (
            _product_doside_axis_operator_matrix(axis_ops, 1, :overlap),
            _product_doside_axis_operator_matrix(axis_ops, 2, :kinetic),
            _product_doside_axis_operator_matrix(axis_ops, 3, :overlap),
        ),
        (
            _product_doside_axis_operator_matrix(axis_ops, 1, :overlap),
            _product_doside_axis_operator_matrix(axis_ops, 2, :overlap),
            _product_doside_axis_operator_matrix(axis_ops, 3, :kinetic),
        ),
    )
end

function _product_doside_retained_kinetic_block(
    left_unit::_CartesianNestedProductStagedByCenterUnit3D,
    right_unit::_CartesianNestedProductStagedByCenterUnit3D,
    axis_ops::NamedTuple{(:x,:y,:z)},
)
    return _product_doside_retained_separable_sum_block(
        left_unit,
        right_unit,
        _product_doside_kinetic_axis_factor_terms(axis_ops),
    )
end

function _pqs_product_low_order_axis_matrices(metrics::NamedTuple{(:x,:y,:z)}, term::Symbol)
    term == :overlap && return (metrics.x.overlap, metrics.y.overlap, metrics.z.overlap)
    term == :position_x && return (metrics.x.position, metrics.y.overlap, metrics.z.overlap)
    term == :position_y && return (metrics.x.overlap, metrics.y.position, metrics.z.overlap)
    term == :position_z && return (metrics.x.overlap, metrics.y.overlap, metrics.z.position)
    throw(
        ArgumentError(
            "PQS/product low-order reference block supports only :overlap and :position_x/y/z",
        ),
    )
end

function _pqs_factored_support_entries(
    pqs_payload::_CartesianExecutableProjectedQShellPayload3D,
    factored_coefficients::AbstractMatrix{<:Real},
)
    return _support_local_retained_entries(
        pqs_payload.column_range,
        pqs_payload.support_states,
        factored_coefficients,
    )
end

function _pqs_factored_support_coefficients(
    pqs_payload::_CartesianExecutableProjectedQShellPayload3D;
    atol::Real,
    context::AbstractString,
)
    seed = _nested_projected_q_shell_descriptor_seed_coefficients(pqs_payload.descriptor)
    size(seed) == (pqs_payload.descriptor.support_count, pqs_payload.descriptor.mode_count) ||
        throw(
            DimensionMismatch("$(context) seed shape does not match descriptor metadata"),
        )
    size(pqs_payload.cleanup_transform) ==
        (pqs_payload.descriptor.mode_count, pqs_payload.descriptor.retained_count) ||
        throw(
            DimensionMismatch("$(context) cleanup transform shape is inconsistent"),
        )
    factored_coefficients = Matrix{Float64}(seed * pqs_payload.cleanup_transform)
    stored_coefficients = Matrix{Float64}(pqs_payload.support_coefficient_matrix)
    size(factored_coefficients) == size(stored_coefficients) || throw(
        DimensionMismatch("$(context) factored and stored support coefficient shapes differ"),
    )
    coefficient_error = LinearAlgebra.norm(factored_coefficients - stored_coefficients, Inf)
    coefficient_error <= atol || throw(
        ArgumentError("$(context) factored support coefficients disagree with stored payload coefficients"),
    )
    return (
        factored_coefficients = factored_coefficients,
        stored_coefficients = stored_coefficients,
        coefficient_error = coefficient_error,
    )
end

function _pqs_product_low_order_reference_block(
    pqs_payload,
    product_unit,
    metrics::NamedTuple{(:x,:y,:z)};
    term::Symbol,
    atol::Real = 1.0e-10,
)
    pqs_payload isa _CartesianExecutableProjectedQShellPayload3D || throw(
        ArgumentError("PQS/product low-order reference block requires an executable PQS payload"),
    )
    product_unit isa _CartesianNestedProductStagedByCenterUnit3D || throw(
        ArgumentError("PQS/product low-order reference block requires a staged product unit"),
    )
    product_unit.kind == :product_doside || throw(
        ArgumentError("PQS/product low-order reference block requires a product_doside product unit"),
    )
    axis_matrices = _pqs_product_low_order_axis_matrices(metrics, term)
    coefficient_check = _pqs_factored_support_coefficients(
        pqs_payload;
        atol,
        context = "PQS/product low-order reference",
    )
    factored_coefficients = coefficient_check.factored_coefficients
    coefficient_error = coefficient_check.coefficient_error

    factored_entries = _pqs_factored_support_entries(pqs_payload, factored_coefficients)
    product_entries = _staged_unit_entries(product_unit)
    oracle_entries = _staged_unit_entries(pqs_payload)
    factored_block = _contract_pair_block(
        factored_entries,
        product_entries,
        axis_matrices[1],
        axis_matrices[2],
        axis_matrices[3],
    )
    oracle_block = _contract_pair_block(
        oracle_entries,
        product_entries,
        axis_matrices[1],
        axis_matrices[2],
        axis_matrices[3],
    )
    block_error = LinearAlgebra.norm(factored_block - oracle_block, Inf)
    block_error <= atol || throw(
        ArgumentError("PQS/product factored low-order block disagrees with support-local oracle"),
    )
    return (
        path = :pqs_product_low_order_reference,
        term = term,
        block = factored_block,
        oracle_block = oracle_block,
        coefficient_error = coefficient_error,
        block_error = block_error,
        diagnostics = (
            source = :pqs_product_low_order_reference_block,
            fixture_reference_only = true,
            production_supported = false,
            packet_adoption = false,
            fixed_block_sidecar_installation = false,
            qwhamiltonian_consumes = false,
            public_default_consumes = false,
            cr2_science_status_changed = false,
            ida_weight_division_allowed = false,
            retained_pqs_weights_role = :debug_reference_only,
            retained_pqs_weights_used = false,
            retained_pqs_weights_positive_checked = false,
            pqs_self_overlap_identity_shortcut_used = false,
            factored_pqs_transform_used = true,
            seed_reconstructed_from_descriptor = true,
            cleanup_transform_stage_applied = true,
            support_coefficient_matrix_compared = true,
            support_local_oracle_used = true,
            optimized_pqs_product_path = false,
            supported_terms = (:overlap, :position_x, :position_y, :position_z),
            unsupported_terms = (
                :kinetic,
                :weights,
                :first_moments,
                :x2,
                :nuclear_one_body,
                :local_coulomb_one_body,
                :local_ecp_one_body,
                :gaussian_local_terms,
                :mwg_interaction,
                :interaction,
            ),
            coefficient_error = coefficient_error,
            block_error = block_error,
        ),
    )
end

function _pqs_product_low_order_reference_block(
    product_unit::_CartesianNestedProductStagedByCenterUnit3D,
    pqs_payload::_CartesianExecutableProjectedQShellPayload3D,
    metrics::NamedTuple{(:x,:y,:z)};
    term::Symbol,
    atol::Real = 1.0e-10,
)
    forward = _pqs_product_low_order_reference_block(
        pqs_payload,
        product_unit,
        metrics;
        term,
        atol,
    )
    return (
        path = :product_pqs_low_order_reference,
        term = term,
        block = Matrix{Float64}(transpose(forward.block)),
        oracle_block = Matrix{Float64}(transpose(forward.oracle_block)),
        coefficient_error = forward.coefficient_error,
        block_error = forward.block_error,
        diagnostics = merge(
            forward.diagnostics,
            (
                source = :product_pqs_low_order_reference_block,
                transposed_from_pqs_product_reference = true,
                pqs_self_overlap_identity_shortcut_used = false,
            ),
        ),
    )
end

function _pqs_product_kinetic_reference_block(
    pqs_payload,
    product_unit,
    axis_ops::NamedTuple{(:x,:y,:z)};
    atol::Real = 1.0e-10,
)
    pqs_payload isa _CartesianExecutableProjectedQShellPayload3D || throw(
        ArgumentError("PQS/product kinetic reference block requires an executable PQS payload"),
    )
    product_unit isa _CartesianNestedProductStagedByCenterUnit3D || throw(
        ArgumentError("PQS/product kinetic reference block requires a staged product unit"),
    )
    product_unit.kind == :product_doside || throw(
        ArgumentError("PQS/product kinetic reference block requires a product_doside product unit"),
    )
    coefficient_check = _pqs_factored_support_coefficients(
        pqs_payload;
        atol,
        context = "PQS/product kinetic reference",
    )
    factored_entries = _pqs_factored_support_entries(
        pqs_payload,
        coefficient_check.factored_coefficients,
    )
    product_entries = _staged_unit_entries(product_unit)
    oracle_entries = _staged_unit_entries(pqs_payload)
    axis_factor_terms = _product_doside_kinetic_axis_factor_terms(axis_ops)
    factored_block = _fallback_staged_separable_sum_block(
        factored_entries,
        product_entries,
        axis_factor_terms,
    )
    oracle_block = _fallback_staged_separable_sum_block(
        oracle_entries,
        product_entries,
        axis_factor_terms,
    )
    block_error = LinearAlgebra.norm(factored_block - oracle_block, Inf)
    block_error <= atol || throw(
        ArgumentError("PQS/product factored kinetic block disagrees with support-local oracle"),
    )
    return (
        path = :pqs_product_kinetic_reference,
        term = :kinetic,
        block = factored_block,
        oracle_block = oracle_block,
        coefficient_error = coefficient_check.coefficient_error,
        block_error = block_error,
        diagnostics = (
            source = :pqs_product_kinetic_reference_block,
            fixture_reference_only = true,
            production_supported = false,
            signed_operator_reference = true,
            retained_weight_semantics = :not_used,
            retained_pqs_weights_role = :debug_reference_only,
            retained_pqs_weights_used = false,
            retained_pqs_weights_positive_checked = false,
            ida_weight_division_allowed = false,
            quadrature_weight_semantics_claimed = false,
            packet_adoption = false,
            fixed_block_sidecar_installation = false,
            qwhamiltonian_consumes = false,
            public_default_consumes = false,
            cr2_science_status_changed = false,
            factored_pqs_transform_used = true,
            seed_reconstructed_from_descriptor = true,
            cleanup_transform_stage_applied = true,
            support_coefficient_matrix_compared = true,
            support_local_oracle_used = true,
            optimized_pqs_product_path = false,
            kinetic_factor_form = (
                (:kinetic, :overlap, :overlap),
                (:overlap, :kinetic, :overlap),
                (:overlap, :overlap, :kinetic),
            ),
            supported_terms = (:kinetic,),
            unsupported_terms = (
                :overlap,
                :position_x,
                :position_y,
                :position_z,
                :weights,
                :first_moments,
                :x2,
                :nuclear_one_body,
                :local_coulomb_one_body,
                :local_ecp_one_body,
                :gaussian_local_terms,
                :mwg_interaction,
                :interaction,
            ),
            coefficient_error = coefficient_check.coefficient_error,
            block_error = block_error,
        ),
    )
end

function _pqs_product_kinetic_reference_block(
    product_unit::_CartesianNestedProductStagedByCenterUnit3D,
    pqs_payload::_CartesianExecutableProjectedQShellPayload3D,
    axis_ops::NamedTuple{(:x,:y,:z)};
    atol::Real = 1.0e-10,
)
    forward = _pqs_product_kinetic_reference_block(
        pqs_payload,
        product_unit,
        axis_ops;
        atol,
    )
    return (
        path = :product_pqs_kinetic_reference,
        term = :kinetic,
        block = Matrix{Float64}(transpose(forward.block)),
        oracle_block = Matrix{Float64}(transpose(forward.oracle_block)),
        coefficient_error = forward.coefficient_error,
        block_error = forward.block_error,
        diagnostics = merge(
            forward.diagnostics,
            (
                source = :product_pqs_kinetic_reference_block,
                transposed_from_pqs_product_reference = true,
            ),
        ),
    )
end

function _pqs_raw_product_box_reference_block(
    raw_product_box_plan;
    term::Symbol,
)
    raw_plan = _pqs_raw_product_box_plan_view(raw_product_box_plan)
    term in (
        :overlap,
        :position_x,
        :position_y,
        :position_z,
        :x2_x,
        :x2_y,
        :x2_z,
        :kinetic,
    ) || throw(
        ArgumentError("PQS raw product-box reference block received unsupported term $(term)"),
    )
    raw_plan.operator_factors_available || throw(
        ArgumentError("PQS raw product-box reference block requires an operator-backed raw plan"),
    )
    raw_plan.source_product_modes_orthogonal === true || throw(
        ArgumentError("PQS raw product-box source modes are not orthogonal within tolerance"),
    )
    block = term == :kinetic ?
            _pqs_raw_product_box_kinetic_mode_matrix(raw_plan) :
            _pqs_raw_product_box_low_order_mode_matrix(raw_plan; term)
    all(isfinite, block) || throw(
        ArgumentError("PQS raw product-box reference block produced non-finite entries"),
    )
    return (
        path = :pqs_mode_selected_raw_product_box_reference,
        term = term,
        block = block,
        raw_product_box_plan = raw_plan,
        diagnostics = (
            source = :pqs_raw_product_box_reference_block,
            pqs_representation = :mode_selected_raw_product_box,
            private_shadow_only = true,
            production_supported = false,
            row_projected_shell_support = false,
            shell_row_projection_used = false,
            lowdin_cleanup_used = false,
            raw_product_box_stage_lowdin_cleanup_used = false,
            shell_projection_realization_stage = :postponed,
            shell_projection_realization_applied = false,
            shell_projection_realization_requires_lowdin = true,
            lowdin_cleanup_scope = :shell_projection_realization_stage_only,
            descriptor_cleanup_transform_ignored = true,
            source_product_modes_orthogonal = raw_plan.source_product_modes_orthogonal,
            operator_factor_source = :pqs_raw_product_box_plan,
            operator_metric_sources =
                hasproperty(raw_plan.diagnostics, :operator_metric_sources) ?
                raw_plan.diagnostics.operator_metric_sources :
                (:unspecified, :unspecified, :unspecified),
            raw_product_box_integration_contract =
                hasproperty(raw_plan.diagnostics, :raw_product_box_integration_contract) ?
                raw_plan.diagnostics.raw_product_box_integration_contract :
                :unknown,
            raw_product_box_numerical_reference_fallback =
                hasproperty(raw_plan.diagnostics, :raw_product_box_numerical_reference_fallback) ?
                raw_plan.diagnostics.raw_product_box_numerical_reference_fallback :
                nothing,
            axis_overlap_errors = raw_plan.axis_overlap_errors,
            max_1d_source_overlap_error = raw_plan.max_1d_source_overlap_error,
            max_product_overlap_error = raw_plan.max_product_overlap_error,
            selected_overlap_error = raw_plan.selected_overlap_error,
            overlap_identity_error = raw_plan.overlap_identity_error,
            boundary_selection_preserves_orthogonality = true,
            product_box_column_selection_reference = true,
            boundary_column_selection_only = true,
            boundary_mode_selection_rule =
                raw_plan.boundary_selector.selection_rule,
            source_mode_dims = raw_plan.source_mode_dims,
            boundary_mode_count = raw_plan.boundary_selector.selected_count,
            retained_count = raw_plan.boundary_selector.selected_count,
            retained_functions_live_in_product_box_mode_span = true,
            retained_functions_live_in_shell_row_support_subspace = false,
            retained_columns_have_full_product_box_support = true,
            support_coefficient_matrix_oracle_used = false,
            shell_row_support_oracle_used = false,
            retained_weight_semantics = :not_positive_quadrature_weights,
            retained_pqs_weights_used = false,
            retained_pqs_weights_positive_checked = false,
            ida_weight_division_allowed = false,
            packet_adoption = false,
            fixed_block_sidecar_installation = false,
            qwhamiltonian_consumes = false,
            public_default_consumes = false,
            cr2_science_status_changed = false,
            local_ecp_gaussian_mwg_implemented = false,
            generic_retained_unit_framework = false,
            supported_terms = (
                :overlap,
                :position_x,
                :position_y,
                :position_z,
                :x2_x,
                :x2_y,
                :x2_z,
                :kinetic,
            ),
            unsupported_terms = (
                :weights,
                :first_moments,
                :nuclear_one_body,
                :local_coulomb_one_body,
                :local_ecp_one_body,
                :gaussian_local_terms,
                :gaussian_sum,
                :pair_sum,
                :mwg_interaction,
                :interaction,
            ),
        ),
    )
end

function _pqs_raw_product_box_reference_block(
    descriptor::_CartesianNestedProjectedQShellStagedUnitDescriptor3D,
    metrics::NamedTuple{(:x,:y,:z)};
    term::Symbol,
    orthogonality_atol::Real = 1.0e-8,
)
    raw_plan = _pqs_raw_product_box_plan(
        descriptor,
        metrics;
        orthogonality_atol = orthogonality_atol,
    )
    return _pqs_raw_product_box_reference_block(raw_plan; term)
end

function _fallback_staged_separable_sum_block(
    left_entries::Vector{Vector{_ParentCoefficientEntry3D}},
    right_entries::Vector{Vector{_ParentCoefficientEntry3D}},
    axis_factor_terms,
)
    !isempty(axis_factor_terms) || throw(
        ArgumentError("fallback staged separable-sum block requires at least one axis-factor triple"),
    )
    block = zeros(Float64, length(left_entries), length(right_entries))
    for term in axis_factor_terms
        length(term) == 3 || throw(
            ArgumentError("fallback staged separable-sum block terms must be axis-factor triples"),
        )
        block .+= _contract_pair_block(
            left_entries,
            right_entries,
            term[1],
            term[2],
            term[3],
        )
    end
    return block
end

function _fallback_staged_kinetic_block(
    left_entries::Vector{Vector{_ParentCoefficientEntry3D}},
    right_entries::Vector{Vector{_ParentCoefficientEntry3D}},
    axis_ops::NamedTuple{(:x,:y,:z)},
)
    return _fallback_staged_separable_sum_block(
        left_entries,
        right_entries,
        _product_doside_kinetic_axis_factor_terms(axis_ops),
    )
end

function _fill_product_staged_metric_blocks!(
    overlap::Matrix{Float64},
    position_x::Matrix{Float64},
    position_y::Matrix{Float64},
    position_z::Matrix{Float64},
    left_unit,
    right_unit,
    metrics::NamedTuple{(:x,:y,:z)},
)
    left_unit.kind == :product_doside && right_unit.kind == :product_doside || throw(
        ArgumentError("product-staged metric block requires product_doside units"),
    )
    length(left_unit.axis_function_indices) == length(left_unit.column_range) || throw(
        ArgumentError("product-staged metric left unit axis metadata does not match its column range"),
    )
    length(right_unit.axis_function_indices) == length(right_unit.column_range) || throw(
        ArgumentError("product-staged metric right unit axis metadata does not match its column range"),
    )
    overlap_block = _product_doside_retained_low_order_block(
        left_unit,
        right_unit,
        metrics;
        term = :overlap,
    )
    position_x_block = _product_doside_retained_low_order_block(
        left_unit,
        right_unit,
        metrics;
        term = :position_x,
    )
    position_y_block = _product_doside_retained_low_order_block(
        left_unit,
        right_unit,
        metrics;
        term = :position_y,
    )
    position_z_block = _product_doside_retained_low_order_block(
        left_unit,
        right_unit,
        metrics;
        term = :position_z,
    )
    left_range = left_unit.column_range
    right_range = right_unit.column_range
    same_unit = left_unit.column_range == right_unit.column_range
    overlap[left_range, right_range] .= overlap_block
    position_x[left_range, right_range] .= position_x_block
    position_y[left_range, right_range] .= position_y_block
    position_z[left_range, right_range] .= position_z_block
    if !same_unit
        overlap[right_range, left_range] .= transpose(overlap_block)
        position_x[right_range, left_range] .= transpose(position_x_block)
        position_y[right_range, left_range] .= transpose(position_y_block)
        position_z[right_range, left_range] .= transpose(position_z_block)
    end
    return nothing
end

function _support_local_retained_entries(
    column_range::UnitRange{Int},
    support_states::AbstractVector{<:NTuple{3,Int}},
    coefficients,
)
    length(column_range) == size(coefficients, 2) || throw(
        DimensionMismatch("support-local retained entries require one coefficient column per retained column"),
    )
    length(support_states) == size(coefficients, 1) || throw(
        DimensionMismatch("support-local retained entries require one support state per coefficient row"),
    )
    entries = [Vector{_ParentCoefficientEntry3D}() for _ in column_range]
    for local_col in axes(coefficients, 2)
        column_entries = entries[local_col]
        for local_row in axes(coefficients, 1)
            value = Float64(coefficients[local_row, local_col])
            iszero(value) && continue
            ix, iy, iz = support_states[local_row]
            push!(column_entries, _ParentCoefficientEntry3D(ix, iy, iz, value))
        end
    end
    return entries
end

function _staged_unit_entries(unit)
    return _support_local_retained_entries(
        unit.column_range,
        unit.support_states,
        unit.coefficient_matrix,
    )
end

function _staged_unit_entries(unit::_CartesianExecutableProjectedQShellPayload3D)
    return _support_local_retained_entries(
        unit.column_range,
        unit.support_states,
        unit.support_coefficient_matrix,
    )
end

function _contract_pair_block(
    left_entries::Vector{Vector{_ParentCoefficientEntry3D}},
    right_entries::Vector{Vector{_ParentCoefficientEntry3D}},
    mx::AbstractMatrix{<:Real},
    my::AbstractMatrix{<:Real},
    mz::AbstractMatrix{<:Real},
)
    result = zeros(Float64, length(left_entries), length(right_entries))
    for col_b in eachindex(right_entries), col_a in eachindex(left_entries)
        value = 0.0
        for left in left_entries[col_a], right in right_entries[col_b]
            value += _metric_value(left, right, mx, my, mz)
        end
        result[col_a, col_b] = value
    end
    return result
end

function _fallback_staged_metric_block(
    left_entries::Vector{Vector{_ParentCoefficientEntry3D}},
    right_entries::Vector{Vector{_ParentCoefficientEntry3D}},
    axis_matrices::NTuple{3,AbstractMatrix{<:Real}},
)
    return _contract_pair_block(
        left_entries,
        right_entries,
        axis_matrices[1],
        axis_matrices[2],
        axis_matrices[3],
    )
end

function _cached_staged_unit_entries!(
    entries_cache::Vector{Union{Nothing,Vector{Vector{_ParentCoefficientEntry3D}}}},
    index::Int,
    unit,
)
    cached = entries_cache[index]
    cached === nothing || return cached
    entries = _staged_unit_entries(unit)
    entries_cache[index] = entries
    return entries
end

function _fallback_staged_metric_blocks(
    left_entries::Vector{Vector{_ParentCoefficientEntry3D}},
    right_entries::Vector{Vector{_ParentCoefficientEntry3D}},
    metrics::NamedTuple{(:x,:y,:z)},
)
    return (
        overlap = _fallback_staged_metric_block(
            left_entries,
            right_entries,
            (metrics.x.overlap, metrics.y.overlap, metrics.z.overlap),
        ),
        position_x = _fallback_staged_metric_block(
            left_entries,
            right_entries,
            (metrics.x.position, metrics.y.overlap, metrics.z.overlap),
        ),
        position_y = _fallback_staged_metric_block(
            left_entries,
            right_entries,
            (metrics.x.overlap, metrics.y.position, metrics.z.overlap),
        ),
        position_z = _fallback_staged_metric_block(
            left_entries,
            right_entries,
            (metrics.x.overlap, metrics.y.overlap, metrics.z.position),
        ),
    )
end

function _fallback_staged_x2_blocks(
    left_entries::Vector{Vector{_ParentCoefficientEntry3D}},
    right_entries::Vector{Vector{_ParentCoefficientEntry3D}},
    metrics::NamedTuple{(:x,:y,:z)},
)
    return (
        x2_x = _fallback_staged_metric_block(
            left_entries,
            right_entries,
            (metrics.x.x2, metrics.y.overlap, metrics.z.overlap),
        ),
        x2_y = _fallback_staged_metric_block(
            left_entries,
            right_entries,
            (metrics.x.overlap, metrics.y.x2, metrics.z.overlap),
        ),
        x2_z = _fallback_staged_metric_block(
            left_entries,
            right_entries,
            (metrics.x.overlap, metrics.y.overlap, metrics.z.x2),
        ),
    )
end

function _product_doside_retained_linear_vectors(
    unit::_CartesianNestedProductStagedByCenterUnit3D,
    metrics::NamedTuple{(:x,:y,:z)},
)
    _require_product_doside_retained_block_unit(unit; side = :linear)
    wx = _project_staged_axis_vector(unit.axes[1], metrics.x.weights)
    wy = _project_staged_axis_vector(unit.axes[2], metrics.y.weights)
    wz = _project_staged_axis_vector(unit.axes[3], metrics.z.weights)
    xw = _project_staged_axis_vector(unit.axes[1], metrics.x.centers .* metrics.x.weights)
    yw = _project_staged_axis_vector(unit.axes[2], metrics.y.centers .* metrics.y.weights)
    zw = _project_staged_axis_vector(unit.axes[3], metrics.z.centers .* metrics.z.weights)
    weights = zeros(Float64, length(unit.column_range))
    first_moments = zeros(Float64, length(unit.column_range), 3)
    @inbounds for column in eachindex(unit.axis_function_indices)
        xi, yi, zi = unit.axis_function_indices[column]
        weights[column] = wx[xi] * wy[yi] * wz[zi]
        first_moments[column, 1] = xw[xi] * wy[yi] * wz[zi]
        first_moments[column, 2] = wx[xi] * yw[yi] * wz[zi]
        first_moments[column, 3] = wx[xi] * wy[yi] * zw[zi]
    end
    return weights, first_moments
end

function _staged_unit_linear_vectors(unit, metrics::NamedTuple{(:x,:y,:z)})
    if unit.kind == :product_doside
        return _product_doside_retained_linear_vectors(unit, metrics)
    end
    entries = _staged_unit_entries(unit)
    weights = _contract_linear_vector(
        entries,
        metrics.x.weights,
        metrics.y.weights,
        metrics.z.weights,
    )
    first_moments = _contracted_first_moments(entries, metrics)
    return weights, first_moments
end

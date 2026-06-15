const _PRODUCT_DOSIDE_SOURCE_BOX_REFERENCE_TERMS = (
    :overlap,
    :position_x,
    :position_y,
    :position_z,
    :x2_x,
    :x2_y,
    :x2_z,
    :kinetic,
)

# CCPM retirement note: this product/product source-box family remains only as
# a private oracle and route-shadow bridge while old diagnostics are retired.
# New safe one-body product-box checks should use CartesianCPBBlockProviders
# axis-product or sum-of-axis-products blocks directly.
function _product_doside_source_box_axis_centers(
    retained_unit_plan,
    metrics::NamedTuple{(:x,:y,:z)},
)
    return ntuple(axis -> begin
        interval = retained_unit_plan.source_axis_intervals[axis]
        Float64.(getproperty(getproperty(metrics, (:x, :y, :z)[axis]), :centers)[interval])
    end, 3)
end

function _product_doside_source_box_cross_factors(
    left_unit::_CartesianNestedProductStagedByCenterUnit3D,
    right_unit::_CartesianNestedProductStagedByCenterUnit3D,
    metrics::NamedTuple{(:x,:y,:z)},
)
    return ntuple(axis -> (
        overlap = _project_staged_axis_matrix(
            left_unit.axes[axis],
            right_unit.axes[axis],
            _product_doside_axis_metric_matrix(metrics, axis, :overlap),
        ),
        position = _project_staged_axis_matrix(
            left_unit.axes[axis],
            right_unit.axes[axis],
            _product_doside_axis_metric_matrix(metrics, axis, :position),
        ),
        x2 = _project_staged_axis_matrix(
            left_unit.axes[axis],
            right_unit.axes[axis],
            _product_doside_axis_metric_matrix(metrics, axis, :x2),
        ),
        kinetic = _project_staged_axis_matrix(
            left_unit.axes[axis],
            right_unit.axes[axis],
            _product_doside_axis_metric_matrix(metrics, axis, :kinetic),
        ),
    ), 3)
end

function _product_doside_source_box_pair_plan(
    left_unit::_CartesianNestedProductStagedByCenterUnit3D,
    right_unit::_CartesianNestedProductStagedByCenterUnit3D,
    metrics::NamedTuple{(:x,:y,:z)},
)
    left_retained_unit_plan = _product_doside_retained_unit_plan(left_unit)
    right_retained_unit_plan = _product_doside_retained_unit_plan(right_unit)
    cross_factors =
        _product_doside_source_box_cross_factors(left_unit, right_unit, metrics)
    left_axis_centers =
        _product_doside_source_box_axis_centers(left_retained_unit_plan, metrics)
    right_axis_centers =
        _product_doside_source_box_axis_centers(right_retained_unit_plan, metrics)
    return (
        pair_kind = :product_doside_source_box_pair,
        left_source_family = :product_doside,
        right_source_family = :product_doside,
        left_retained_rule_kind = left_retained_unit_plan.retained_rule_kind,
        right_retained_rule_kind = right_retained_unit_plan.retained_rule_kind,
        left_source_dimensions = left_retained_unit_plan.source_axis_lengths,
        right_source_dimensions = right_retained_unit_plan.source_axis_lengths,
        left_source_dimension = left_retained_unit_plan.source_dimension,
        right_source_dimension = right_retained_unit_plan.source_dimension,
        left_retained_axis_counts = left_retained_unit_plan.retained_axis_counts,
        right_retained_axis_counts = right_retained_unit_plan.retained_axis_counts,
        left_column_range = left_retained_unit_plan.column_range,
        right_column_range = right_retained_unit_plan.column_range,
        left_retained_count = left_retained_unit_plan.retained_count,
        right_retained_count = right_retained_unit_plan.retained_count,
        axis_intervals = (
            left = left_retained_unit_plan.source_axis_intervals,
            right = right_retained_unit_plan.source_axis_intervals,
        ),
        axis_centers = (
            left = left_axis_centers,
            right = right_axis_centers,
        ),
        left_retained_unit_plan = left_retained_unit_plan,
        right_retained_unit_plan = right_retained_unit_plan,
        left_retained_transform = left_retained_unit_plan,
        right_retained_transform = right_retained_unit_plan,
        one_dimensional_cross_factors = (
            x = cross_factors[1],
            y = cross_factors[2],
            z = cross_factors[3],
        ),
        supported_terms = _PRODUCT_DOSIDE_SOURCE_BOX_REFERENCE_TERMS,
        diagnostics = (
            source = :product_doside_source_box_pair_plan,
            private_shadow_only = true,
            source_box_pair_plan = true,
            product_doside_retained_unit_plan_used = true,
            existing_product_staged_retained_helpers_authoritative = true,
            product_staged_metric_execution_changed = false,
            product_doside_retained_block_math_changed = false,
            operator_factor_source = :explicit_metric_operator_data,
            operator_metric_sources =
                _cartesian_source_box_metric_sources(metrics),
            input_metric_operator_data = :caller_supplied_explicit_data,
            input_metric_operator_data_pgdg_checked = false,
            pgdg_analytic_operator_provenance_claimed = false,
            numerical_reference_fallback = false,
            raw_product_box_operators_use_1d_factors = true,
            shell_projection_used = false,
            lowdin_cleanup_used = false,
            retained_pqs_weights_used = false,
            ida_weight_division_allowed = false,
            ida_mwg_semantics_changed = false,
            retained_weight_semantics_changed = false,
            packet_adoption = false,
            fixed_block_routing = false,
            qwhamiltonian_consumes = false,
            public_default_consumes = false,
            cr2_science_status_changed = false,
            local_ecp_gaussian_mwg_implemented = false,
            generic_retained_unit_framework = false,
        ),
    )
end
function _product_doside_source_box_factor(
    pair_plan,
    axis::Int,
    kind::Symbol,
)
    axis_factors = getproperty(pair_plan.one_dimensional_cross_factors, (:x, :y, :z)[axis])
    return getproperty(axis_factors, kind)
end

function _product_doside_source_box_block_from_factors(
    pair_plan,
    term::Symbol,
)
    term in pair_plan.supported_terms || throw(
        ArgumentError("product/doside source-box reference block received unsupported term $(term)"),
    )
    left_modes = pair_plan.left_retained_transform.axis_function_indices
    right_modes = pair_plan.right_retained_transform.axis_function_indices
    block = zeros(
        Float64,
        pair_plan.left_retained_count,
        pair_plan.right_retained_count,
    )
    @inbounds for factor_kinds in _source_box_separable_term_factor_kinds(term)
        fx = _product_doside_source_box_factor(pair_plan, 1, factor_kinds[1])
        fy = _product_doside_source_box_factor(pair_plan, 2, factor_kinds[2])
        fz = _product_doside_source_box_factor(pair_plan, 3, factor_kinds[3])
        for col in eachindex(right_modes)
            xj, yj, zj = right_modes[col]
            for row in eachindex(left_modes)
                xi, yi, zi = left_modes[row]
                block[row, col] += fx[xi, xj] * fy[yi, yj] * fz[zi, zj]
            end
        end
    end
    return block
end

function _product_doside_source_box_authoritative_block(
    left_unit::_CartesianNestedProductStagedByCenterUnit3D,
    right_unit::_CartesianNestedProductStagedByCenterUnit3D,
    metrics::NamedTuple{(:x,:y,:z)},
    term::Symbol,
)
    term == :kinetic && return _product_doside_retained_kinetic_block(
        left_unit,
        right_unit,
        metrics,
    )
    return _product_doside_retained_low_order_block(
        left_unit,
        right_unit,
        metrics;
        term,
    )
end

function _product_doside_source_box_reference_block(
    left_unit::_CartesianNestedProductStagedByCenterUnit3D,
    right_unit::_CartesianNestedProductStagedByCenterUnit3D,
    metrics::NamedTuple{(:x,:y,:z)};
    term::Symbol,
    atol::Real = 1.0e-12,
)
    pair_plan = _product_doside_source_box_pair_plan(left_unit, right_unit, metrics)
    term in pair_plan.supported_terms || throw(
        ArgumentError("product/doside source-box reference block received unsupported term $(term)"),
    )
    block = _product_doside_source_box_block_from_factors(pair_plan, term)
    authoritative_block =
        _product_doside_source_box_authoritative_block(
            left_unit,
            right_unit,
            metrics,
            term,
        )
    block_error = LinearAlgebra.norm(block - authoritative_block, Inf)
    block_error <= atol || throw(
        ArgumentError("product/doside source-box reference block disagrees with authoritative retained helper"),
    )
    return (
        path = :product_doside_source_box_reference,
        term = term,
        block = block,
        authoritative_block = authoritative_block,
        block_error = block_error,
        pair_plan = pair_plan,
        diagnostics = merge(
            pair_plan.diagnostics,
            (
                source = :product_doside_source_box_reference_block,
                supported_terms = pair_plan.supported_terms,
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
                kinetic_factor_form = (
                    (:kinetic, :overlap, :overlap),
                    (:overlap, :kinetic, :overlap),
                    (:overlap, :overlap, :kinetic),
                ),
                authoritative_helper =
                    term == :kinetic ?
                    :_product_doside_retained_kinetic_block :
                    :_product_doside_retained_low_order_block,
                authoritative_block_compared = true,
                block_error = block_error,
            ),
        ),
    )
end

function _product_doside_source_box_project_axis_terms(
    operator_terms::AbstractArray{<:Real,3},
    left_axis,
    right_axis;
    label::AbstractString,
)
    nterms = size(operator_terms, 1)
    nterms > 0 || throw(
        ArgumentError("$(label) requires at least one term"),
    )
    projected = Array{Float64,3}(
        undef,
        nterms,
        _staged_axis_count(left_axis),
        _staged_axis_count(right_axis),
    )
    @inbounds for term in 1:nterms
        term_matrix = @view operator_terms[term, :, :]
        projected[term, :, :] .= _project_staged_axis_matrix(
            left_axis,
            right_axis,
            term_matrix,
        )
    end
    return projected
end

function _product_doside_source_box_project_local_gaussian_axis_terms(
    operator_terms::AbstractArray{<:Real,3},
    left_axis,
    right_axis,
)
    return _product_doside_source_box_project_axis_terms(
        operator_terms,
        left_axis,
        right_axis;
        label = "product/doside source-box local-Gaussian terms",
    )
end

function _product_doside_source_box_local_gaussian_axis_factors(
    left_unit::_CartesianNestedProductStagedByCenterUnit3D,
    right_unit::_CartesianNestedProductStagedByCenterUnit3D,
    axis_gaussian_terms::NamedTuple{(:x,:y,:z)},
)
    return ntuple(axis -> begin
        terms = getproperty(axis_gaussian_terms, (:x, :y, :z)[axis])
        _product_doside_source_box_project_local_gaussian_axis_terms(
            terms,
            left_unit.axes[axis],
            right_unit.axes[axis],
        )
    end, 3)
end

function _product_doside_source_box_local_gaussian_sum_block(
    left_unit::_CartesianNestedProductStagedByCenterUnit3D,
    right_unit::_CartesianNestedProductStagedByCenterUnit3D;
    term_coefficients::AbstractVector{<:Real},
    axis_gaussian_terms::NamedTuple{(:x,:y,:z)},
)
    left_unit.kind == :product_doside && right_unit.kind == :product_doside || throw(
        ArgumentError("product/doside source-box local-Gaussian block requires product_doside units"),
    )
    length(left_unit.axis_function_indices) == length(left_unit.column_range) || throw(
        ArgumentError("product/doside source-box local-Gaussian left unit axis metadata does not match its column range"),
    )
    length(right_unit.axis_function_indices) == length(right_unit.column_range) || throw(
        ArgumentError("product/doside source-box local-Gaussian right unit axis metadata does not match its column range"),
    )
    nterms = length(term_coefficients)
    nterms > 0 || throw(
        ArgumentError("product/doside source-box local-Gaussian block requires at least one term coefficient"),
    )
    for axis_name in (:x, :y, :z)
        size(getproperty(axis_gaussian_terms, axis_name), 1) == nterms || throw(
            ArgumentError("product/doside source-box local-Gaussian axis term count mismatch on $(axis_name)"),
        )
    end
    left_retained_unit_plan = _product_doside_retained_unit_plan(left_unit)
    right_retained_unit_plan = _product_doside_retained_unit_plan(right_unit)
    pair_plan = (
        pair_kind = :product_doside_source_box_local_gaussian_pair,
        left_source_family = :product_doside,
        right_source_family = :product_doside,
        left_retained_rule_kind = left_retained_unit_plan.retained_rule_kind,
        right_retained_rule_kind = right_retained_unit_plan.retained_rule_kind,
        left_source_dimensions = left_retained_unit_plan.source_axis_lengths,
        right_source_dimensions = right_retained_unit_plan.source_axis_lengths,
        left_source_dimension = left_retained_unit_plan.source_dimension,
        right_source_dimension = right_retained_unit_plan.source_dimension,
        left_retained_axis_counts = left_retained_unit_plan.retained_axis_counts,
        right_retained_axis_counts = right_retained_unit_plan.retained_axis_counts,
        left_column_range = left_retained_unit_plan.column_range,
        right_column_range = right_retained_unit_plan.column_range,
        left_retained_count = left_retained_unit_plan.retained_count,
        right_retained_count = right_retained_unit_plan.retained_count,
        axis_intervals = (
            left = left_retained_unit_plan.source_axis_intervals,
            right = right_retained_unit_plan.source_axis_intervals,
        ),
        left_retained_unit_plan = left_retained_unit_plan,
        right_retained_unit_plan = right_retained_unit_plan,
        left_retained_transform = left_retained_unit_plan,
        right_retained_transform = right_retained_unit_plan,
        supported_terms = (:gaussian_sum,),
        diagnostics = (
            source = :product_doside_source_box_local_gaussian_pair_plan,
            private_shadow_only = true,
            source_box_pair_plan = true,
            product_doside_retained_unit_plan_used = true,
            operator_factor_source = :explicit_local_gaussian_axis_terms,
            input_local_gaussian_data = :caller_supplied_explicit_data,
            input_local_gaussian_data_pgdg_checked = false,
            pgdg_analytic_operator_provenance_claimed = false,
            numerical_reference_fallback = false,
            raw_product_box_operators_use_1d_factors = true,
            shell_projection_used = false,
            lowdin_cleanup_used = false,
            retained_pqs_weights_used = false,
            retained_pqs_weight_division_allowed = false,
            ida_weight_division_allowed = false,
            ida_mwg_semantics_changed = false,
            retained_weight_semantics_changed = false,
            packet_adoption = false,
            fixed_block_routing = false,
            qwhamiltonian_consumes = false,
            public_default_consumes = false,
            cr2_science_status_changed = false,
            generic_retained_unit_framework = false,
        ),
    )
    projected_terms =
        _product_doside_source_box_local_gaussian_axis_factors(
            left_unit,
            right_unit,
            axis_gaussian_terms,
        )
    coeffs = Float64[Float64(value) for value in term_coefficients]
    left_modes = left_unit.axis_function_indices
    right_modes = right_unit.axis_function_indices
    block = zeros(Float64, length(left_modes), length(right_modes))
    projected_x, projected_y, projected_z = projected_terms
    @inbounds for col in eachindex(right_modes)
        xj, yj, zj = right_modes[col]
        for row in eachindex(left_modes)
            xi, yi, zi = left_modes[row]
            value = 0.0
            @simd for term in 1:nterms
                value +=
                    coeffs[term] *
                    projected_x[term, xi, xj] *
                    projected_y[term, yi, yj] *
                    projected_z[term, zi, zj]
            end
            block[row, col] = value
        end
    end
    all(isfinite, block) || throw(
        ArgumentError("product/doside source-box local-Gaussian block produced non-finite entries"),
    )
    return (
        path = :product_doside_source_box_local_gaussian_sum,
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
                source = :product_doside_source_box_local_gaussian_sum_block,
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
                shell_row_algorithm = false,
                support_local_oracle_used = false,
                support_local_oracle_only = false,
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

function _source_box_axis_positive_weights(
    axis_weights::AbstractVector{<:Real},
    interval::UnitRange{Int};
    axis_name::Symbol,
    side::Symbol,
)
    length(axis_weights) >= last(interval) || throw(
        DimensionMismatch("$(side) $(axis_name) source weights do not cover the source interval"),
    )
    values = Float64[Float64(value) for value in axis_weights[interval]]
    all(isfinite, values) || throw(
        ArgumentError("$(side) $(axis_name) source weights must be finite"),
    )
    all(>(0.0), values) || throw(
        ArgumentError("$(side) $(axis_name) source weights must be positive"),
    )
    return values
end

function _product_doside_source_box_pair_factor_axis_terms(
    left_unit::_CartesianNestedProductStagedByCenterUnit3D,
    right_unit::_CartesianNestedProductStagedByCenterUnit3D,
    axis_pair_factor_terms::NamedTuple{(:x,:y,:z)},
)
    return ntuple(axis -> begin
        terms = getproperty(axis_pair_factor_terms, (:x, :y, :z)[axis])
        _product_doside_source_box_project_axis_terms(
            terms,
            left_unit.axes[axis],
            right_unit.axes[axis];
            label = "product/doside source-box density-density pair-factor terms",
        )
    end, 3)
end

function _product_doside_source_box_density_density_weight_views(
    left_retained_unit_plan,
    right_retained_unit_plan,
    axis_weights::NamedTuple{(:x,:y,:z)},
)
    return (
        left = ntuple(axis -> _source_box_axis_positive_weights(
            getproperty(axis_weights, (:x, :y, :z)[axis]),
            left_retained_unit_plan.source_axis_intervals[axis];
            axis_name = (:x, :y, :z)[axis],
            side = :left,
        ), 3),
        right = ntuple(axis -> _source_box_axis_positive_weights(
            getproperty(axis_weights, (:x, :y, :z)[axis]),
            right_retained_unit_plan.source_axis_intervals[axis];
            axis_name = (:x, :y, :z)[axis],
            side = :right,
        ), 3),
    )
end

function _product_doside_source_box_density_density_interaction_block(
    left_unit::_CartesianNestedProductStagedByCenterUnit3D,
    right_unit::_CartesianNestedProductStagedByCenterUnit3D;
    term_coefficients::AbstractVector{<:Real},
    axis_pair_factor_terms::NamedTuple{(:x,:y,:z)},
    axis_weights::NamedTuple{(:x,:y,:z)},
    pair_factor_normalization::Symbol = :density_normalized,
)
    pair_factor_normalization == :density_normalized || throw(
        ArgumentError(
            "product/doside source-box density-density fixture currently requires density-normalized pair factors",
        ),
    )
    left_unit.kind == :product_doside && right_unit.kind == :product_doside || throw(
        ArgumentError("product/doside source-box density-density block requires product_doside units"),
    )
    length(left_unit.axis_function_indices) == length(left_unit.column_range) || throw(
        ArgumentError("product/doside source-box density-density left unit axis metadata does not match its column range"),
    )
    length(right_unit.axis_function_indices) == length(right_unit.column_range) || throw(
        ArgumentError("product/doside source-box density-density right unit axis metadata does not match its column range"),
    )
    nterms = length(term_coefficients)
    nterms > 0 || throw(
        ArgumentError("product/doside source-box density-density block requires at least one term coefficient"),
    )
    for axis_name in (:x, :y, :z)
        size(getproperty(axis_pair_factor_terms, axis_name), 1) == nterms || throw(
            ArgumentError("product/doside source-box density-density pair-factor term count mismatch on $(axis_name)"),
        )
    end
    left_retained_unit_plan = _product_doside_retained_unit_plan(left_unit)
    right_retained_unit_plan = _product_doside_retained_unit_plan(right_unit)
    source_weights = _product_doside_source_box_density_density_weight_views(
        left_retained_unit_plan,
        right_retained_unit_plan,
        axis_weights,
    )
    projected_terms =
        _product_doside_source_box_pair_factor_axis_terms(
            left_unit,
            right_unit,
            axis_pair_factor_terms,
        )
    coeffs = Float64[Float64(value) for value in term_coefficients]
    left_modes = left_unit.axis_function_indices
    right_modes = right_unit.axis_function_indices
    block = zeros(Float64, length(left_modes), length(right_modes))
    projected_x, projected_y, projected_z = projected_terms
    @inbounds for col in eachindex(right_modes)
        xj, yj, zj = right_modes[col]
        for row in eachindex(left_modes)
            xi, yi, zi = left_modes[row]
            value = 0.0
            @simd for term in 1:nterms
                value +=
                    coeffs[term] *
                    projected_x[term, xi, xj] *
                    projected_y[term, yi, yj] *
                    projected_z[term, zi, zj]
            end
            block[row, col] = value
        end
    end
    all(isfinite, block) || throw(
        ArgumentError("product/doside source-box density-density block produced non-finite entries"),
    )
    pair_plan = (
        pair_kind = :product_doside_source_box_density_density_pair,
        left_source_family = :product_doside,
        right_source_family = :product_doside,
        left_retained_rule_kind = left_retained_unit_plan.retained_rule_kind,
        right_retained_rule_kind = right_retained_unit_plan.retained_rule_kind,
        left_source_dimensions = left_retained_unit_plan.source_axis_lengths,
        right_source_dimensions = right_retained_unit_plan.source_axis_lengths,
        left_source_dimension = left_retained_unit_plan.source_dimension,
        right_source_dimension = right_retained_unit_plan.source_dimension,
        left_retained_axis_counts = left_retained_unit_plan.retained_axis_counts,
        right_retained_axis_counts = right_retained_unit_plan.retained_axis_counts,
        left_column_range = left_retained_unit_plan.column_range,
        right_column_range = right_retained_unit_plan.column_range,
        left_retained_count = left_retained_unit_plan.retained_count,
        right_retained_count = right_retained_unit_plan.retained_count,
        axis_intervals = (
            left = left_retained_unit_plan.source_axis_intervals,
            right = right_retained_unit_plan.source_axis_intervals,
        ),
        left_retained_unit_plan = left_retained_unit_plan,
        right_retained_unit_plan = right_retained_unit_plan,
        left_retained_transform = left_retained_unit_plan,
        right_retained_transform = right_retained_unit_plan,
        source_weights = source_weights,
        supported_terms = (:pair_sum,),
        diagnostics = (
            source = :product_doside_source_box_density_density_pair_plan,
            private_shadow_only = true,
            source_box_pair_plan = true,
            product_doside_retained_unit_plan_used = true,
            interaction_operator = :electron_electron_density_density,
            output_representation = :two_index_density_density,
            four_index_galerkin_tensor = false,
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
            retained_pqs_weights_used = false,
            retained_weight_division_allowed = false,
            retained_pqs_weight_division_allowed = false,
            ida_weight_division_allowed = false,
            ida_mwg_semantics_changed = false,
            mwg_ida_semantics_changed = false,
            packet_adoption = false,
            fixed_block_routing = false,
            qwhamiltonian_consumes = false,
            public_default_consumes = false,
            cr2_science_status_changed = false,
            numerical_reference_fallback = false,
            generic_retained_unit_framework = false,
        ),
    )
    return (
        path = :product_doside_source_box_density_density_interaction,
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
                source = :product_doside_source_box_density_density_interaction_block,
                source_box_first = true,
                output_finite = true,
                electron_electron_terms_implemented = true,
                local_gaussian_one_body_implemented = false,
                ecp_terms_implemented = false,
                mwg_interaction_implemented = false,
                support_local_oracle_used = false,
                shell_row_algorithm = false,
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

function _source_box_density_normalized_axis_pair_terms(
    raw_axis_pair_terms::AbstractArray{<:Real,3},
    axis_weights::AbstractVector{<:Real};
    axis_name::Symbol,
)
    nterms = size(raw_axis_pair_terms, 1)
    nterms > 0 || throw(
        ArgumentError("$(axis_name) raw pair-factor terms require at least one term"),
    )
    parent_count = size(raw_axis_pair_terms, 2)
    size(raw_axis_pair_terms, 3) == parent_count || throw(
        DimensionMismatch("$(axis_name) raw pair-factor terms must be square per term"),
    )
    length(axis_weights) >= parent_count || throw(
        DimensionMismatch("$(axis_name) source weights do not cover raw pair-factor dimensions"),
    )
    weights = Float64[Float64(value) for value in axis_weights[1:parent_count]]
    all(isfinite, weights) || throw(
        ArgumentError("$(axis_name) source weights must be finite"),
    )
    all(>(0.0), weights) || throw(
        ArgumentError("$(axis_name) source weights must be positive"),
    )
    weight_outer = weights * transpose(weights)
    normalized = Array{Float64,3}(undef, size(raw_axis_pair_terms))
    @inbounds for term in 1:nterms
        normalized[term, :, :] .= raw_axis_pair_terms[term, :, :] ./ weight_outer
    end
    return normalized
end

function _source_box_density_normalized_axis_pair_terms(
    raw_axis_pair_terms::NamedTuple{(:x,:y,:z)},
    axis_weights::NamedTuple{(:x,:y,:z)},
)
    return (
        x = _source_box_density_normalized_axis_pair_terms(
            raw_axis_pair_terms.x,
            axis_weights.x;
            axis_name = :x,
        ),
        y = _source_box_density_normalized_axis_pair_terms(
            raw_axis_pair_terms.y,
            axis_weights.y;
            axis_name = :y,
        ),
        z = _source_box_density_normalized_axis_pair_terms(
            raw_axis_pair_terms.z,
            axis_weights.z;
            axis_name = :z,
        ),
    )
end

function _pqs_source_box_ida_expected_term_count(expected_term_count)
    isnothing(expected_term_count) && return nothing
    term_count = Int(expected_term_count)
    term_count == expected_term_count || throw(
        ArgumentError("expected IDA pair-factor term count must be an integer"),
    )
    term_count > 0 || throw(
        ArgumentError("expected IDA pair-factor term count must be positive"),
    )
    return term_count
end

function _pqs_source_box_ida_axis_factor_provenance(
    pgdg::_MappedOrdinaryPGDGIntermediate1D;
    axis::Symbol = :x,
    expected_term_count = nothing,
    bundle_backend = nothing,
    provenance_source::Symbol = :mapped_ordinary_pgdg_intermediate,
)
    axis in (:x, :y, :z) || throw(
        ArgumentError("PQS source-box IDA axis provenance requires axis :x, :y, or :z"),
    )
    density_terms = pgdg.pair_factor_terms
    raw_terms = pgdg.pair_factor_terms_raw
    size(density_terms) == size(raw_terms) || throw(
        DimensionMismatch("PQS source-box IDA raw and density-normalized pair-factor terms must have the same shape"),
    )
    term_count = size(density_terms, 1)
    term_count > 0 || throw(
        ArgumentError("PQS source-box IDA provenance requires at least one pair-factor term"),
    )
    nrow = size(density_terms, 2)
    ncol = size(density_terms, 3)
    nrow == ncol || throw(
        DimensionMismatch("PQS source-box IDA pair-factor terms must be square per axis"),
    )
    nrow > 0 || throw(
        ArgumentError("PQS source-box IDA pair-factor dimension must be positive"),
    )
    expected_value = _pqs_source_box_ida_expected_term_count(expected_term_count)
    if !isnothing(expected_value)
        term_count == expected_value || throw(
            DimensionMismatch("PQS source-box IDA pair-factor term count does not match expected term count"),
        )
    end
    all(isfinite, density_terms) || throw(
        ArgumentError("PQS source-box IDA density-normalized pair-factor terms must be finite"),
    )
    all(isfinite, raw_terms) || throw(
        ArgumentError("PQS source-box IDA raw pair-factor terms must be finite"),
    )
    weights = Float64[Float64(weight) for weight in pgdg.weights]
    length(weights) == nrow || throw(
        DimensionMismatch("PQS source-box IDA source weights must match pair-factor dimension"),
    )
    all(isfinite, weights) || throw(
        ArgumentError("PQS source-box IDA source weights must be finite"),
    )
    all(weight -> abs(weight) > 1.0e-12, weights) || throw(
        ArgumentError("PQS source-box IDA source weights must be nonzero"),
    )
    center_values = Float64[Float64(point) for point in pgdg.centers]
    length(center_values) == nrow || throw(
        DimensionMismatch("PQS source-box IDA centers must match pair-factor dimension"),
    )
    all(isfinite, center_values) || throw(
        ArgumentError("PQS source-box IDA centers must be finite"),
    )

    weight_outer = weights * transpose(weights)
    reconstruction_error = 0.0
    @inbounds for term in 1:term_count
        for column in 1:nrow, row in 1:nrow
            reconstructed = Float64(raw_terms[term, row, column]) / weight_outer[row, column]
            reconstruction_error = max(
                reconstruction_error,
                abs(Float64(density_terms[term, row, column]) - reconstructed),
            )
        end
    end

    return (
        object_kind = :pqs_source_box_ida_axis_factor_provenance,
        axis = axis,
        interaction_path = :ida_gausslet_source_box,
        density_normalized_pair_factor_terms = density_terms,
        raw_pair_factor_terms = raw_terms,
        source_raw_quadrature_weights = weights,
        source_weights = weights,
        centers = center_values,
        term_count = term_count,
        factor_dimensions = (nrow, ncol),
        factor_shape = size(density_terms),
        backend = pgdg.backend,
        bundle_backend = bundle_backend,
        provenance_metadata = (
            source = provenance_source,
            pgdg_backend = pgdg.backend,
            bundle_backend = bundle_backend,
            raw_pair_factor_terms_field = :pair_factor_terms_raw,
            density_normalized_pair_factor_terms_field = :pair_factor_terms,
            source_weights_field = :weights,
            centers_field = :centers,
        ),
        diagnostics = (
            source = :pqs_source_box_ida_axis_factor_provenance,
            private_diagnostic_only = true,
            interaction_path = :ida_gausslet_source_box,
            mwg_supplement_residual_path = false,
            density_normalized_pair_factors = true,
            raw_weighted_pair_factors_available = true,
            density_normalized_from_raw_weights = true,
            source_weight_division_owner = :pgdg_auxiliary_source_weights,
            source_weight_division_shape = :axis_pair_weight_outer,
            raw_terms_field = :pair_factor_terms_raw,
            density_normalized_terms_field = :pair_factor_terms,
            source_weights_field = :weights,
            source_centers_field = :centers,
            source_weights_finite = true,
            source_weights_nonzero = true,
            source_weights_positive = all(>(0.0), weights),
            term_count = term_count,
            factor_dimensions = (nrow, ncol),
            density_normalized_reconstruction_error = reconstruction_error,
            retained_pqs_weights_used = false,
            retained_pqs_weights_positive_checked = false,
            retained_weight_division_allowed = false,
            retained_pqs_weight_division_allowed = false,
            ida_weight_division_allowed = false,
            packet_adoption = false,
            fixed_block_routing = false,
            qwhamiltonian_consumes = false,
            public_default_consumes = false,
            ecp_terms_implemented = false,
            cr2_science_status_changed = false,
            mwg_ida_semantics_changed = false,
            ida_mwg_semantics_changed = false,
            route_adapter_connected = false,
        ),
    )
end

function _pqs_source_box_ida_axis_factor_provenance(
    bundle::_MappedOrdinaryGausslet1DBundle;
    axis::Symbol = :x,
    expected_term_count = nothing,
)
    return _pqs_source_box_ida_axis_factor_provenance(
        bundle.pgdg_intermediate;
        axis,
        expected_term_count,
        bundle_backend = bundle.backend,
        provenance_source = :mapped_ordinary_gausslet_1d_bundle,
    )
end

function _pqs_source_box_ida_factor_provenance(
    bundles::_CartesianNestedAxisBundles3D;
    expected_term_count = nothing,
)
    axis_data = (
        x = _pqs_source_box_ida_axis_factor_provenance(
            _nested_axis_bundle(bundles, :x);
            axis = :x,
            expected_term_count,
        ),
        y = _pqs_source_box_ida_axis_factor_provenance(
            _nested_axis_bundle(bundles, :y);
            axis = :y,
            expected_term_count,
        ),
        z = _pqs_source_box_ida_axis_factor_provenance(
            _nested_axis_bundle(bundles, :z);
            axis = :z,
            expected_term_count,
        ),
    )
    term_count = axis_data.x.term_count
    term_count == axis_data.y.term_count == axis_data.z.term_count || throw(
        DimensionMismatch("PQS source-box IDA axis provenance requires matching x/y/z term counts"),
    )
    return (
        object_kind = :pqs_source_box_ida_factor_provenance,
        interaction_path = :ida_gausslet_source_box,
        axes = axis_data,
        axis_pair_factor_terms = (
            x = axis_data.x.density_normalized_pair_factor_terms,
            y = axis_data.y.density_normalized_pair_factor_terms,
            z = axis_data.z.density_normalized_pair_factor_terms,
        ),
        raw_axis_pair_factor_terms = (
            x = axis_data.x.raw_pair_factor_terms,
            y = axis_data.y.raw_pair_factor_terms,
            z = axis_data.z.raw_pair_factor_terms,
        ),
        axis_weights = (
            x = axis_data.x.source_weights,
            y = axis_data.y.source_weights,
            z = axis_data.z.source_weights,
        ),
        axis_centers = (
            x = axis_data.x.centers,
            y = axis_data.y.centers,
            z = axis_data.z.centers,
        ),
        term_count = term_count,
        factor_dimensions = (
            x = axis_data.x.factor_dimensions,
            y = axis_data.y.factor_dimensions,
            z = axis_data.z.factor_dimensions,
        ),
        diagnostics = (
            source = :pqs_source_box_ida_factor_provenance,
            private_diagnostic_only = true,
            interaction_path = :ida_gausslet_source_box,
            mwg_supplement_residual_path = false,
            axis_count = 3,
            term_count = term_count,
            factor_dimensions = (
                x = axis_data.x.factor_dimensions,
                y = axis_data.y.factor_dimensions,
                z = axis_data.z.factor_dimensions,
            ),
            source_weight_division_owner = :pgdg_auxiliary_source_weights,
            source_weight_division_shape = :axis_pair_weight_outer,
            density_normalized_pair_factors = true,
            raw_weighted_pair_factors_available = true,
            retained_pqs_weights_used = false,
            retained_pqs_weights_positive_checked = false,
            retained_weight_division_allowed = false,
            retained_pqs_weight_division_allowed = false,
            ida_weight_division_allowed = false,
            packet_adoption = false,
            fixed_block_routing = false,
            qwhamiltonian_consumes = false,
            public_default_consumes = false,
            ecp_terms_implemented = false,
            cr2_science_status_changed = false,
            mwg_ida_semantics_changed = false,
            ida_mwg_semantics_changed = false,
            route_adapter_connected = false,
        ),
    )
end

function _pqs_source_box_ida_factor_provenance(
    bundle_x::_MappedOrdinaryGausslet1DBundle,
    bundle_y::_MappedOrdinaryGausslet1DBundle,
    bundle_z::_MappedOrdinaryGausslet1DBundle;
    expected_term_count = nothing,
)
    return _pqs_source_box_ida_factor_provenance(
        _CartesianNestedAxisBundles3D(bundle_x, bundle_y, bundle_z);
        expected_term_count,
    )
end

function _product_doside_source_box_raw_weighted_density_density_interaction_block(
    left_unit::_CartesianNestedProductStagedByCenterUnit3D,
    right_unit::_CartesianNestedProductStagedByCenterUnit3D;
    term_coefficients::AbstractVector{<:Real},
    raw_axis_pair_factor_terms::NamedTuple{(:x,:y,:z)},
    axis_weights::NamedTuple{(:x,:y,:z)},
)
    density_normalized_terms = _source_box_density_normalized_axis_pair_terms(
        raw_axis_pair_factor_terms,
        axis_weights,
    )
    density_normalized_core =
        _product_doside_source_box_density_density_interaction_block(
            left_unit,
            right_unit;
            term_coefficients,
            axis_pair_factor_terms = density_normalized_terms,
            axis_weights,
            pair_factor_normalization = :density_normalized,
        )
    return (
        path = :product_doside_source_box_raw_weighted_density_density_interaction,
        interaction_operator = :electron_electron_density_density,
        block = density_normalized_core.block,
        density_normalized_core = density_normalized_core,
        normalized_axis_pair_factor_terms = density_normalized_terms,
        raw_axis_pair_factor_terms = raw_axis_pair_factor_terms,
        pair_plan = density_normalized_core.pair_plan,
        diagnostics = merge(
            density_normalized_core.diagnostics,
            (
                source = :product_doside_source_box_raw_weighted_density_density_interaction_block,
                path = :product_doside_source_box_raw_weighted_density_density_interaction,
                pair_factor_normalization = :raw_weighted,
                raw_weighted_pair_factors = true,
                density_normalized_pair_factors = false,
                density_normalized_pair_factors_generated = true,
                source_weight_division_owner = :source_box_raw_weights,
                source_weight_division_applied_by_helper = true,
                source_weight_division_shape = :axis_pair_weight_outer,
                density_normalized_core_helper =
                    :_product_doside_source_box_density_density_interaction_block,
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

function _cartesian_source_box_term_first_axis_array(
    matrices::AbstractVector,
    nterms::Int,
    axis_name::Symbol,
)
    length(matrices) == nterms || throw(
        ArgumentError("centered local-Gaussian $(axis_name) axis term matrix count does not match coefficients"),
    )
    nterms > 0 || throw(
        ArgumentError("centered local-Gaussian term arrays require at least one term"),
    )
    first_matrix = Matrix{Float64}(matrices[1])
    size(first_matrix, 1) == size(first_matrix, 2) || throw(
        ArgumentError("centered local-Gaussian $(axis_name) axis matrix must be square"),
    )
    parent_count = size(first_matrix, 1)
    terms = Array{Float64,3}(undef, nterms, parent_count, parent_count)
    terms[1, :, :] .= first_matrix
    @inbounds for term in 2:nterms
        matrix = Matrix{Float64}(matrices[term])
        size(matrix) == (parent_count, parent_count) || throw(
            ArgumentError("centered local-Gaussian $(axis_name) axis term matrices must have a common square size"),
        )
        terms[term, :, :] .= matrix
    end
    return terms
end

function _product_doside_source_box_centered_local_gaussian_term_data(
    axis_layers::NamedTuple{(:x,:y,:z)},
    expansion::CoulombGaussianExpansion;
    center::NTuple{3,<:Real},
)
    nterms = length(expansion.coefficients)
    nterms == length(expansion.exponents) || throw(
        ArgumentError("centered local-Gaussian expansion has mismatched coefficient/exponent counts"),
    )
    nterms > 0 || throw(
        ArgumentError("centered local-Gaussian expansion requires at least one term"),
    )
    center_values = ntuple(axis -> Float64(center[axis]), 3)
    term_coefficients = Float64[Float64(value) for value in expansion.coefficients]
    axis_terms = ntuple(axis -> begin
        axis_name = (:x, :y, :z)[axis]
        layer = getproperty(axis_layers, axis_name)
        _require_analytic_primitive_backend(
            primitive_set(layer),
            "product/doside source-box centered local-Gaussian $(axis_name) axis terms",
        )
        matrices = gaussian_factor_matrices(
            layer;
            exponents = expansion.exponents,
            center = center_values[axis],
        )
        _cartesian_source_box_term_first_axis_array(matrices, nterms, axis_name)
    end, 3)
    return (
        term_coefficients = term_coefficients,
        axis_gaussian_terms = (
            x = axis_terms[1],
            y = axis_terms[2],
            z = axis_terms[3],
        ),
        diagnostics = (
            source = :product_doside_source_box_centered_local_gaussian_term_data,
            axis_gaussian_terms_source = :analytic_gaussian_factor_matrices,
            analytic_primitive_backend_required = true,
            analytic_primitive_backend_checked = true,
            center = center_values,
            term_count = nterms,
            axis_term_dimensions = (
                x = size(axis_terms[1]),
                y = size(axis_terms[2]),
                z = size(axis_terms[3]),
            ),
            positive_gaussian_sum_convention = true,
            nuclear_charge_applied = false,
            nuclear_attraction_sign_applied = false,
            input_local_gaussian_data_pgdg_checked = false,
            pgdg_analytic_operator_provenance_claimed = false,
            numerical_reference_fallback = false,
            shell_row_algorithm = false,
            packet_adoption = false,
            fixed_block_routing = false,
            qwhamiltonian_consumes = false,
            public_default_consumes = false,
            cr2_science_status_changed = false,
        ),
    )
end

function _product_doside_source_box_centered_local_gaussian_sum_block(
    left_unit::_CartesianNestedProductStagedByCenterUnit3D,
    right_unit::_CartesianNestedProductStagedByCenterUnit3D,
    axis_layers::NamedTuple{(:x,:y,:z)},
    expansion::CoulombGaussianExpansion;
    center::NTuple{3,<:Real},
)
    term_data = _product_doside_source_box_centered_local_gaussian_term_data(
        axis_layers,
        expansion;
        center,
    )
    explicit_block = _product_doside_source_box_local_gaussian_sum_block(
        left_unit,
        right_unit;
        term_coefficients = term_data.term_coefficients,
        axis_gaussian_terms = term_data.axis_gaussian_terms,
    )
    return (
        path = :product_doside_source_box_centered_local_gaussian_sum,
        block = explicit_block.block,
        explicit_block = explicit_block,
        term_data = term_data,
        pair_plan = explicit_block.pair_plan,
        one_dimensional_gaussian_factors =
            explicit_block.one_dimensional_gaussian_factors,
        diagnostics = merge(
            explicit_block.diagnostics,
            term_data.diagnostics,
            (
                source = :product_doside_source_box_centered_local_gaussian_sum_block,
                source_box_first = true,
                local_gaussian_source_box_terms = true,
                centered_local_gaussian_terms_generated = true,
                axis_gaussian_terms_source = :analytic_gaussian_factor_matrices,
                explicit_table_helper_used =
                    :_product_doside_source_box_local_gaussian_sum_block,
                positive_gaussian_sum_convention = true,
                nuclear_charge_applied = false,
                nuclear_attraction_sign_applied = false,
                numerical_reference_fallback = false,
                shell_row_algorithm = false,
                support_local_oracle_used = false,
                retained_pqs_weight_division_allowed = false,
                ida_weight_division_allowed = false,
                ida_mwg_semantics_changed = false,
                ecp_terms_implemented = false,
                electron_electron_terms_implemented = false,
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

# Three-unit PQS/PQS/product shadow inventory remains route-diagnostic only.
# Do not extend it as a production placement or provider-layer contract.
function _pqs_pqs_product_source_box_all_pairs_inventory(
    left_raw_plan,
    right_raw_plan,
    product_unit::_CartesianNestedProductStagedByCenterUnit3D,
    ranges,
    supported_terms::Tuple{Vararg{Symbol}},
)
    product_retained_unit_plan = _product_doside_retained_unit_plan(product_unit)
    retained_units = (
        (
            unit_key = :pqs_left,
            retained_unit_kind = :pqs,
            source_family = :mode_selected_raw_product_box,
            retained_range = ranges.pqs_left,
            source_dimensions = left_raw_plan.source_mode_dims,
            source_dimension = left_raw_plan.source_mode_count,
            retained_count = left_raw_plan.boundary_selector.selected_count,
            supported_safe_terms = supported_terms,
        ),
        (
            unit_key = :pqs_right,
            retained_unit_kind = :pqs,
            source_family = :mode_selected_raw_product_box,
            retained_range = ranges.pqs_right,
            source_dimensions = right_raw_plan.source_mode_dims,
            source_dimension = right_raw_plan.source_mode_count,
            retained_count = right_raw_plan.boundary_selector.selected_count,
            supported_safe_terms = supported_terms,
        ),
        (
            unit_key = :product,
            retained_unit_kind = :product_doside,
            source_family = :product_doside,
            retained_range = ranges.product,
            source_dimensions = product_retained_unit_plan.source_axis_lengths,
            source_dimension = product_retained_unit_plan.source_dimension,
            retained_count = product_retained_unit_plan.retained_count,
            supported_safe_terms = supported_terms,
        ),
    )
    pair_entries = (
        (
            pair_key = (:pqs_left, :pqs_left),
            pair_kind = :pqs_pqs_source_box,
            block_helper = :_pqs_pqs_source_box_reference_blocks,
            upper_triangular = true,
            transpose_only_lower_block = false,
            pair_policy = :source_box_algorithm_available,
            source_box_algorithmic = true,
        ),
        (
            pair_key = (:pqs_left, :pqs_right),
            pair_kind = :pqs_pqs_source_box,
            block_helper = :_pqs_pqs_source_box_reference_blocks,
            upper_triangular = true,
            transpose_only_lower_block = true,
            pair_policy = :source_box_algorithm_available,
            source_box_algorithmic = true,
        ),
        (
            pair_key = (:pqs_left, :product),
            pair_kind = :pqs_product_source_box,
            block_helper = :_pqs_product_source_box_reference_blocks,
            upper_triangular = true,
            transpose_only_lower_block = true,
            pair_policy = :source_box_algorithm_available,
            source_box_algorithmic = true,
        ),
        (
            pair_key = (:pqs_right, :pqs_right),
            pair_kind = :pqs_pqs_source_box,
            block_helper = :_pqs_pqs_source_box_reference_blocks,
            upper_triangular = true,
            transpose_only_lower_block = false,
            pair_policy = :source_box_algorithm_available,
            source_box_algorithmic = true,
        ),
        (
            pair_key = (:pqs_right, :product),
            pair_kind = :pqs_product_source_box,
            block_helper = :_pqs_product_source_box_reference_blocks,
            upper_triangular = true,
            transpose_only_lower_block = true,
            pair_policy = :source_box_algorithm_available,
            source_box_algorithmic = true,
        ),
        (
            pair_key = (:product, :product),
            pair_kind = :product_doside_source_box_pair,
            block_helper = :_product_doside_source_box_reference_block,
            upper_triangular = true,
            transpose_only_lower_block = false,
            pair_policy = :source_box_algorithm_available,
            source_box_algorithmic = true,
        ),
    )
    return (
        object_kind = :pqs_pqs_product_source_box_all_pairs_inventory,
        retained_units = retained_units,
        pair_entries = pair_entries,
        supported_terms = supported_terms,
        diagnostics = (
            source = :pqs_pqs_product_source_box_all_pairs_inventory,
            all_pairs_inventory_private = true,
            pair_inventory_complete_for_units = (:pqs_left, :pqs_right, :product),
            retained_unit_count = length(retained_units),
            upper_triangular_pair_count = length(pair_entries),
            expected_upper_triangular_pair_count = 6,
            pair_keys = map(entry -> entry.pair_key, pair_entries),
            block_helpers = map(entry -> entry.block_helper, pair_entries),
            pair_policies = map(entry -> entry.pair_policy, pair_entries),
            every_pair_uses_source_box_algorithmic_policy = all(
                entry -> entry.source_box_algorithmic,
                pair_entries,
            ),
            source_box_algorithmic_pair_count =
                count(entry -> entry.source_box_algorithmic, pair_entries),
            private_shadow_only = true,
            packet_adoption = false,
            fixed_block_routing = false,
            qwhamiltonian_consumes = false,
            public_default_consumes = false,
            cr2_science_status_changed = false,
            shell_projection_used = false,
            lowdin_cleanup_used = false,
            support_local_pqs_oracle_used = false,
            retained_pqs_weights_used = false,
            retained_pqs_weights_positive_checked = false,
            retained_weight_semantics = :not_positive_quadrature_weights,
            ida_weight_division_allowed = false,
            dense_raw_source_box_pair_matrix_materialized = false,
            dense_raw_pair_storage_avoided = true,
            generic_retained_unit_framework = false,
            lower_triangular_cross_blocks_transpose_only = true,
        ),
    )
end

const _PQS_PQS_PRODUCT_SAFE_TERM_ROUTE_KINDS = (
    :pqs_pqs_product_source_box_safe_term_route,
    :homonuclear_pqs_product_source_box_safe_term_fixture,
)

function _pqs_pqs_product_safe_term_route_descriptor(
    left_pqs_plan,
    right_pqs_plan,
    product_unit::_CartesianNestedProductStagedByCenterUnit3D;
    route_name::Symbol = :pqs_pqs_product_source_box_safe_term_route,
    parent_dims = nothing,
    bond_axis = nothing,
    metadata = (;),
    provenance = (;),
    supported_terms = _PQS_PRODUCT_SOURCE_BOX_SHADOW_TERMS,
)
    left_raw_plan = _pqs_raw_product_box_plan_view(left_pqs_plan)
    right_raw_plan = _pqs_raw_product_box_plan_view(right_pqs_plan)
    left_raw_plan.representation == :orthogonal_raw_product_box ||
        throw(ArgumentError("PQS/PQS/product route descriptor requires a raw product-box left PQS plan"))
    right_raw_plan.representation == :orthogonal_raw_product_box ||
        throw(ArgumentError("PQS/PQS/product route descriptor requires a raw product-box right PQS plan"))
    _require_product_doside_retained_block_unit(product_unit; side = :product)
    selected_terms = _pqs_pqs_product_supported_safe_terms(supported_terms)
    product_retained_unit_plan = _product_doside_retained_unit_plan(product_unit)

    left_count = left_raw_plan.boundary_selector.selected_count
    right_count = right_raw_plan.boundary_selector.selected_count
    product_count = product_retained_unit_plan.retained_count
    ranges = (
        pqs_left = 1:left_count,
        pqs_right = (left_count + 1):(left_count + right_count),
        product = (left_count + right_count + 1):(left_count + right_count + product_count),
    )
    retained_dimension = left_count + right_count + product_count
    unit_summaries = (
        (
            unit_key = :pqs_left,
            retained_unit_kind = :pqs,
            source_family = :mode_selected_raw_product_box,
            retained_rule_kind = :boundary_comx_product_mode_selection,
            retained_range = ranges.pqs_left,
            source_dimensions = left_raw_plan.source_mode_dims,
            source_dimension = left_raw_plan.source_mode_count,
            axis_intervals = left_raw_plan.axis_intervals,
            source_mode_ordering = left_raw_plan.source_mode_ordering,
            boundary_mode_count = left_count,
            retained_count = left_count,
            supported_safe_terms = selected_terms,
        ),
        (
            unit_key = :pqs_right,
            retained_unit_kind = :pqs,
            source_family = :mode_selected_raw_product_box,
            retained_rule_kind = :boundary_comx_product_mode_selection,
            retained_range = ranges.pqs_right,
            source_dimensions = right_raw_plan.source_mode_dims,
            source_dimension = right_raw_plan.source_mode_count,
            axis_intervals = right_raw_plan.axis_intervals,
            source_mode_ordering = right_raw_plan.source_mode_ordering,
            boundary_mode_count = right_count,
            retained_count = right_count,
            supported_safe_terms = selected_terms,
        ),
        (
            unit_key = :product,
            retained_unit_kind = :product_doside,
            source_family = :product_doside,
            retained_rule_kind = :product_doside,
            retained_range = ranges.product,
            source_dimensions = product_retained_unit_plan.source_axis_lengths,
            source_dimension = product_retained_unit_plan.source_dimension,
            axis_intervals = product_retained_unit_plan.source_axis_intervals,
            retained_axis_counts = product_retained_unit_plan.retained_axis_counts,
            retained_count = product_count,
            supported_safe_terms = selected_terms,
        ),
    )
    route_metadata = merge((parent_dims = parent_dims, bond_axis = bond_axis), metadata)
    return (
        object_kind = :pqs_pqs_product_safe_term_route_descriptor,
        route_kind = :pqs_pqs_product_source_box_safe_term_route,
        route_name = route_name,
        roles = (:pqs_left, :pqs_right, :product),
        units = (
            pqs_left = left_pqs_plan,
            pqs_right = right_pqs_plan,
            product = product_unit,
        ),
        unit_summaries = unit_summaries,
        expected_ranges = ranges,
        ranges = ranges,
        retained_dimension = retained_dimension,
        retained_unit_count = length(unit_summaries),
        expected_pair_count = 6,
        pair_count = 6,
        supported_terms = selected_terms,
        metadata = route_metadata,
        provenance = provenance,
        diagnostics = (
            source = :pqs_pqs_product_safe_term_route_descriptor,
            route_descriptor_only = true,
            private_shadow_only = true,
            packet_adoption = false,
            fixed_block_routing = false,
            qwhamiltonian_consumes = false,
            public_default_consumes = false,
            cr2_science_status_changed = false,
            shell_projection_used = false,
            lowdin_cleanup_used = false,
            support_local_pqs_oracle_used = false,
            retained_weight_semantics = :not_positive_quadrature_weights,
            ida_weight_division_allowed = false,
            dense_raw_source_box_pair_matrix_materialized = false,
            dense_raw_pair_storage_avoided = true,
            generic_retained_unit_framework = false,
            operator_algebra_changed = false,
            geometry_construction = false,
            retained_unit_count = length(unit_summaries),
            expected_pair_count = 6,
            term_count = length(selected_terms),
            supported_terms = selected_terms,
        ),
    )
end

function _pqs_validate_source_box_inside_parent_dims(
    source_box::NTuple{3,UnitRange{Int}},
    parent_dims::NTuple{3,Int};
    role::Symbol,
)
    all(dim -> dim > 0, parent_dims) || throw(
        ArgumentError("$(role) parent dimensions must be positive"),
    )
    for axis in 1:3
        interval = source_box[axis]
        !isempty(interval) || throw(
            ArgumentError("$(role) source-box axis $(axis) must be nonempty"),
        )
        first(interval) >= 1 && last(interval) <= parent_dims[axis] || throw(
            ArgumentError("$(role) source-box axis $(axis) lies outside parent dimensions"),
        )
    end
    return source_box
end

function _pqs_product_doside_identity_slab_unit(
    product_source_box::NTuple{3,UnitRange{Int}},
    parent_dims::NTuple{3,Int};
    role::Symbol = :middle_body_product_slab,
    provenance = (;),
)
    _pqs_validate_source_box_inside_parent_dims(
        product_source_box,
        parent_dims;
        role,
    )
    fixed_axes = findall(axis -> length(product_source_box[axis]) == 1, 1:3)
    length(fixed_axes) == 1 || throw(
        ArgumentError("identity product/doside slab requires exactly one fixed source-box axis"),
    )
    fixed_axis = only(fixed_axes)
    active_axes = Tuple(axis for axis in 1:3 if axis != fixed_axis)
    first_axis, second_axis = active_axes
    first_interval = product_source_box[first_axis]
    second_interval = product_source_box[second_axis]
    fixed_index = only(product_source_box[fixed_axis])
    first_count = length(first_interval)
    second_count = length(second_interval)
    retained_count = first_count * second_count

    support_states = NTuple{3,Int}[]
    for first_state in first_interval, second_state in second_interval
        push!(
            support_states,
            ntuple(axis -> axis == fixed_axis ? fixed_index :
                            axis == first_axis ? first_state : second_state, 3),
        )
    end
    support_indices = [
        _cartesian_flat_index(state[1], state[2], state[3], parent_dims) for
        state in support_states
    ]
    axes = ntuple(axis -> begin
        if axis == fixed_axis
            _nested_product_staged_fixed_axis(fixed_index)
        elseif axis == first_axis
            _nested_product_staged_active_axis(
                first_interval,
                Matrix{Float64}(LinearAlgebra.I, first_count, first_count),
            )
        elseif axis == second_axis
            _nested_product_staged_active_axis(
                second_interval,
                Matrix{Float64}(LinearAlgebra.I, second_count, second_count),
            )
        else
            throw(ArgumentError("identity product/doside slab has inconsistent axes"))
        end
    end, 3)
    axis_function_indices = _nested_product_axis_function_indices(
        fixed_axis,
        first_axis,
        first_count,
        second_axis,
        second_count,
    )
    return _CartesianNestedProductStagedByCenterUnit3D(
        role,
        :product_doside,
        1:retained_count,
        support_indices,
        support_states,
        Matrix{Float64}(LinearAlgebra.I, retained_count, retained_count),
        axes,
        axis_function_indices,
        merge(
            (
                source = :pqs_product_doside_identity_slab_unit,
                product_source_box = product_source_box,
                product_retained_rule = :identity_product_doside_slab,
            ),
            provenance,
        ),
        (
            support_count = length(support_indices),
            retained_count = retained_count,
            source_axis_lengths = ntuple(axis -> length(product_source_box[axis]), 3),
            fixed_axis = fixed_axis,
            fixed_index = fixed_index,
            active_axes = active_axes,
            raw_product_box_plan_required = false,
            shell_projection_used = false,
            lowdin_cleanup_used = false,
            support_local_pqs_oracle_used = false,
            retained_weight_semantics = :not_positive_quadrature_weights,
            ida_weight_division_allowed = false,
        ),
    )
end

function _pqs_route_axis_index(axis::Symbol)
    axis == :x && return 1
    axis == :y && return 2
    axis == :z && return 3
    throw(ArgumentError("PQS raw-box geometry producer requires bond_axis = :x, :y, or :z"))
end

function _pqs_geometry_int3(value, label::Symbol)
    value_tuple = Tuple(value)
    length(value_tuple) == 3 || throw(
        ArgumentError("$(label) must contain exactly three entries"),
    )
    return ntuple(axis -> Int(value_tuple[axis]), 3)
end

function _pqs_geometry_source_mode_dims(;
    q = nothing,
    L = nothing,
    source_mode_dims = nothing,
    bond_axis_index::Int,
)
    if isnothing(source_mode_dims)
        !isnothing(q) || throw(
            ArgumentError("raw-box geometry facts require q when source_mode_dims is not supplied"),
        )
        q_int = Int(q)
        l_int = isnothing(L) ? q_int : Int(L)
        return ntuple(axis -> axis == bond_axis_index ? l_int : q_int, 3)
    end
    return _pqs_geometry_int3(source_mode_dims, :source_mode_dims)
end

function _pqs_pqs_product_raw_box_homonuclear_geometry_facts(;
    parent_dims::NTuple{3,Int},
    bond_axis::Symbol = :z,
    q = nothing,
    L = nothing,
    source_mode_dims = nothing,
    left_start = (1, 1, 1),
    right_shift,
    product_slab_fixed_index,
    route_name::Symbol = :homonuclear_pqs_product_source_box_geometry_fixture,
    metadata = (;),
    provenance = (;),
)
    all(dim -> dim > 0, parent_dims) || throw(
        ArgumentError("raw-box geometry facts parent dimensions must be positive"),
    )
    bond_axis_index = _pqs_route_axis_index(bond_axis)
    resolved_source_mode_dims = _pqs_geometry_source_mode_dims(
        q = q,
        L = L,
        source_mode_dims = source_mode_dims,
        bond_axis_index = bond_axis_index,
    )
    all(dim -> dim >= 2, resolved_source_mode_dims) || throw(
        ArgumentError("raw-box geometry facts source-mode dimensions must be total lengths of at least two per axis"),
    )
    left_start_tuple = _pqs_geometry_int3(left_start, :left_start)
    right_shift_tuple = _pqs_geometry_int3(right_shift, :right_shift)
    product_fixed_index = Int(product_slab_fixed_index)
    left_source_box = ntuple(
        axis -> left_start_tuple[axis]:(left_start_tuple[axis] + resolved_source_mode_dims[axis] - 1),
        3,
    )
    right_source_box = ntuple(
        axis -> (left_start_tuple[axis] + right_shift_tuple[axis]):
                (left_start_tuple[axis] + right_shift_tuple[axis] +
                 resolved_source_mode_dims[axis] - 1),
        3,
    )
    product_source_box = ntuple(
        axis -> axis == bond_axis_index ?
                (product_fixed_index:product_fixed_index) :
                left_source_box[axis],
        3,
    )
    _pqs_validate_source_box_inside_parent_dims(
        left_source_box,
        parent_dims;
        role = :pqs_left_geometry,
    )
    _pqs_validate_source_box_inside_parent_dims(
        right_source_box,
        parent_dims;
        role = :pqs_right_geometry,
    )
    _pqs_validate_source_box_inside_parent_dims(
        product_source_box,
        parent_dims;
        role = :product_geometry,
    )
    return (
        object_kind = :pqs_pqs_product_raw_box_homonuclear_geometry_facts,
        status = :private_fixture_geometry_producer,
        parent_dims = parent_dims,
        bond_axis = bond_axis,
        bond_axis_index = bond_axis_index,
        q = isnothing(q) ? nothing : Int(q),
        L = isnothing(L) ? nothing : Int(L),
        source_mode_dims = resolved_source_mode_dims,
        left_start = left_start_tuple,
        right_shift = right_shift_tuple,
        product_slab_fixed_index = product_fixed_index,
        left_source_box = left_source_box,
        right_source_box = right_source_box,
        product_source_box = product_source_box,
        route_name = route_name,
        metadata = merge(
            (
                parent_dims = parent_dims,
                bond_axis = bond_axis,
                pqs_left_box = left_source_box,
                pqs_right_box = right_source_box,
                product_source_box = product_source_box,
                pqs_source_mode_dims = resolved_source_mode_dims,
                geometry_producer =
                    :pqs_pqs_product_raw_box_homonuclear_geometry_facts,
            ),
            metadata,
        ),
        provenance = merge(
            (source = :pqs_pqs_product_raw_box_homonuclear_geometry_facts,),
            provenance,
        ),
        diagnostics = (
            source = :pqs_pqs_product_raw_box_homonuclear_geometry_facts,
            private_fixture_geometry_producer = true,
            public_route_builder = false,
            emits_explicit_source_box_inputs = true,
            geometry_policy = :small_homonuclear_fixture_offsets,
            source_mode_dims_are_total_lengths = true,
            source_boxes_inside_parent_dims = true,
            product_slab_fixed_axis = bond_axis_index,
            product_slab_has_exactly_one_fixed_axis = true,
            raw_box_route_producer_consumes = true,
            operator_algebra_defined = false,
            shell_projection_used = false,
            lowdin_cleanup_used = false,
            support_local_pqs_oracle_used = false,
            support_coefficient_matrix_used = false,
            retained_pqs_weights_used = false,
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

function _pqs_pqs_product_raw_box_route_from_geometry_facts(
    bundles,
    geometry_facts,
    metrics::NamedTuple{(:x,:y,:z)};
    supported_terms = _PQS_PRODUCT_SOURCE_BOX_SHADOW_TERMS,
    orthogonality_atol::Real = 1.0e-8,
)
    hasproperty(geometry_facts, :object_kind) &&
        geometry_facts.object_kind ==
            :pqs_pqs_product_raw_box_homonuclear_geometry_facts || throw(
                ArgumentError("raw-box geometry route producer requires homonuclear geometry facts"),
            )
    _nested_axis_lengths(bundles) == geometry_facts.parent_dims || throw(
        ArgumentError("raw-box geometry facts parent dimensions must match axis bundles"),
    )
    produced_route = _pqs_pqs_product_raw_box_route_producer(
        bundles,
        geometry_facts.left_source_box,
        geometry_facts.right_source_box,
        geometry_facts.product_source_box,
        metrics;
        source_mode_dims = geometry_facts.source_mode_dims,
        route_name = geometry_facts.route_name,
        parent_dims = geometry_facts.parent_dims,
        bond_axis = geometry_facts.bond_axis,
        metadata = geometry_facts.metadata,
        provenance = geometry_facts.provenance,
        supported_terms = supported_terms,
        orthogonality_atol = orthogonality_atol,
    )
    return merge(
        produced_route,
        (
            object_kind = :pqs_pqs_product_raw_box_geometry_route_producer,
            geometry_facts = geometry_facts,
            produced_route = produced_route,
            diagnostics = merge(
                produced_route.diagnostics,
                (
                    source = :pqs_pqs_product_raw_box_route_from_geometry_facts,
                    private_fixture_geometry_producer = true,
                    geometry_facts_consumed = true,
                    emits_explicit_source_box_inputs = true,
                    public_route_builder = false,
                    operator_algebra_defined = false,
                    shell_projection_used = false,
                    lowdin_cleanup_used = false,
                    support_local_pqs_oracle_used = false,
                    support_coefficient_matrix_used = false,
                    retained_pqs_weights_used = false,
                    retained_weight_semantics = :not_positive_quadrature_weights,
                    ida_weight_division_allowed = false,
                    packet_adoption = false,
                    fixed_block_routing = false,
                    qwhamiltonian_consumes = false,
                    public_default_consumes = false,
                    cr2_science_status_changed = false,
                ),
            ),
        ),
    )
end

function _pqs_pqs_product_raw_box_route_producer(
    bundles,
    left_source_box::NTuple{3,UnitRange{Int}},
    right_source_box::NTuple{3,UnitRange{Int}},
    product_source_box::NTuple{3,UnitRange{Int}},
    metrics::NamedTuple{(:x,:y,:z)};
    source_mode_dims::NTuple{3,Int},
    route_name::Symbol = :homonuclear_pqs_product_source_box_safe_term_fixture,
    parent_dims = _nested_axis_lengths(bundles),
    bond_axis = nothing,
    metadata = (;),
    provenance = (;),
    supported_terms = _PQS_PRODUCT_SOURCE_BOX_SHADOW_TERMS,
    orthogonality_atol::Real = 1.0e-8,
)
    selected_terms = _pqs_pqs_product_supported_safe_terms(supported_terms)
    parent_dims isa NTuple{3,Int} || throw(
        ArgumentError("raw-box route producer parent_dims must be a 3-tuple of Int"),
    )
    all(dim -> dim >= 2, source_mode_dims) || throw(
        ArgumentError("raw-box route producer PQS source-mode dimensions must be total lengths of at least two per axis"),
    )
    _pqs_validate_source_box_inside_parent_dims(
        left_source_box,
        parent_dims;
        role = :pqs_left,
    )
    _pqs_validate_source_box_inside_parent_dims(
        right_source_box,
        parent_dims;
        role = :pqs_right,
    )
    _pqs_validate_source_box_inside_parent_dims(
        product_source_box,
        parent_dims;
        role = :product,
    )
    left_raw_product_box_plan = _cartesian_raw_product_box_plan(
        bundles,
        left_source_box,
        source_mode_dims;
        enforce_symmetric_odd = false,
        orthogonality_atol,
    )
    right_raw_product_box_plan = _cartesian_raw_product_box_plan(
        bundles,
        right_source_box,
        source_mode_dims;
        enforce_symmetric_odd = false,
        orthogonality_atol,
    )
    left_pqs_plan = _pqs_raw_product_box_plan_from_raw_product_box_plan(
        left_raw_product_box_plan,
        metrics;
        orthogonality_atol,
    )
    right_pqs_plan = _pqs_raw_product_box_plan_from_raw_product_box_plan(
        right_raw_product_box_plan,
        metrics;
        orthogonality_atol,
    )
    product_unit = _pqs_product_doside_identity_slab_unit(
        product_source_box,
        parent_dims;
        role = :middle_body_product_slab,
        provenance = merge(
            (source = :pqs_pqs_product_raw_box_route_producer,),
            provenance,
        ),
    )
    route_metadata = merge(
        (
            parent_dims = parent_dims,
            bond_axis = bond_axis,
            pqs_left_box = left_source_box,
            pqs_right_box = right_source_box,
            product_source_box = product_source_box,
            pqs_source_mode_dims = source_mode_dims,
            route_producer = :pqs_pqs_product_raw_box_route_producer,
        ),
        metadata,
    )
    route_provenance = merge(
        (source = :pqs_pqs_product_raw_box_route_producer,),
        provenance,
    )
    descriptor = _pqs_pqs_product_safe_term_route_descriptor(
        left_pqs_plan,
        right_pqs_plan,
        product_unit;
        route_name,
        parent_dims,
        bond_axis,
        metadata = route_metadata,
        provenance = route_provenance,
        supported_terms = selected_terms,
    )
    all_pairs_inventory = _pqs_pqs_product_source_box_all_pairs_inventory(
        _pqs_raw_product_box_plan_view(left_pqs_plan),
        _pqs_raw_product_box_plan_view(right_pqs_plan),
        product_unit,
        descriptor.expected_ranges,
        selected_terms,
    )
    all_pairs_inventory.diagnostics.every_pair_uses_source_box_algorithmic_policy ||
        throw(ArgumentError("raw-box route producer requires all pairs to use source-box algorithms"))
    return (
        object_kind = :pqs_pqs_product_raw_box_route_producer,
        status = :private_shadow_only,
        descriptor = descriptor,
        route_descriptor = descriptor,
        raw_product_box_plans = (
            pqs_left = left_raw_product_box_plan,
            pqs_right = right_raw_product_box_plan,
        ),
        raw_pqs_plans = (
            pqs_left = left_pqs_plan,
            pqs_right = right_pqs_plan,
        ),
        product_unit = product_unit,
        retained_rules = (
            pqs_left = :boundary_comx_product_mode_selection,
            pqs_right = :boundary_comx_product_mode_selection,
            product = :identity_product_doside_slab,
        ),
        all_pairs_inventory = all_pairs_inventory,
        supported_terms = selected_terms,
        diagnostics = (
            source = :pqs_pqs_product_raw_box_route_producer,
            private_shadow_only = true,
            raw_product_box_plan_built = true,
            retained_rule_built = true,
            route_descriptor_emitted = true,
            route_descriptor_object_kind = descriptor.object_kind,
            source_box_plan_contract = :RawProductBoxPlan,
            retained_rule_contract = :RetainedRule,
            pqs_retained_rule_kind = :boundary_comx_product_mode_selection,
            product_retained_rule_kind = :identity_product_doside_slab,
            source_mode_dims_are_total_lengths = true,
            every_pair_uses_source_box_algorithmic_policy =
                all_pairs_inventory.diagnostics.every_pair_uses_source_box_algorithmic_policy,
            source_box_algorithmic_pair_count =
                all_pairs_inventory.diagnostics.source_box_algorithmic_pair_count,
            pair_policies = all_pairs_inventory.diagnostics.pair_policies,
            dense_raw_source_box_pair_matrix_materialized = false,
            dense_raw_source_box_pair_matrix_materialized_by_producer = false,
            dense_raw_source_box_pair_matrices_validation_only = true,
            shell_projection_used = false,
            lowdin_cleanup_used = false,
            support_local_pqs_oracle_used = false,
            support_coefficient_matrix_used = false,
            retained_pqs_weights_used = false,
            retained_pqs_weights_positive_checked = false,
            retained_weight_semantics = :not_positive_quadrature_weights,
            ida_weight_division_allowed = false,
            packet_adoption = false,
            fixed_block_routing = false,
            qwhamiltonian_consumes = false,
            public_default_consumes = false,
            cr2_science_status_changed = false,
            local_ecp_gaussian_mwg_implemented = false,
            generic_retained_unit_framework = false,
        ),
    )
end

function _pqs_route_descriptor_diagnostic_common(;
    source::Symbol,
    route_kind = nothing,
    pqs_descriptor_count::Int,
    pqs_raw_plan_convertible_count::Int,
    product_doside_unit_count::Int,
    direct_or_support_body_piece_count::Int,
    descriptor_emitted::Bool,
    supported_terms = _PQS_PRODUCT_SOURCE_BOX_SHADOW_TERMS,
    extra = (;),
)
    return merge(
        (
            source = source,
            private_diagnostic_only = true,
            route_fact_adapter = true,
            route_kind = route_kind,
            pqs_descriptor_count = pqs_descriptor_count,
            pqs_raw_plan_convertible_count = pqs_raw_plan_convertible_count,
            product_doside_unit_count = product_doside_unit_count,
            direct_or_support_body_piece_count =
                direct_or_support_body_piece_count,
            descriptor_emitted = descriptor_emitted,
            packet_adoption = false,
            fixed_block_construction_changed = false,
            qwhamiltonian_changed = false,
            shell_projection_used = false,
            lowdin_cleanup_used = false,
            support_local_pqs_oracle_used = false,
            retained_weight_semantics = :not_positive_quadrature_weights,
            ida_weight_division_allowed = false,
            sidecar_mutation = false,
            sidecar_installation = false,
            direct_support_reinterpreted_as_product_doside = false,
            local_ecp_gaussian_mwg_interaction_changed = false,
            public_default_behavior_changed = false,
            supported_terms =
                _pqs_pqs_product_supported_safe_terms(supported_terms),
        ),
        extra,
    )
end

function _pqs_route_descriptor_missing_symbols(
    pqs_plan_count::Int,
    product_doside_unit_count::Int,
)
    missing = Symbol[]
    if pqs_plan_count == 0
        push!(missing, :pqs_raw_plans)
    elseif pqs_plan_count == 1
        push!(missing, :second_pqs_raw_plan)
    end
    product_doside_unit_count == 0 &&
        push!(missing, :middle_product_doside_unit)
    return missing
end

function _pqs_route_descriptor_unit_mismatch(role::Symbol, unit, reason::Symbol)
    return (
        role = role,
        actual_type = nameof(typeof(unit)),
        actual_kind = hasproperty(unit, :kind) ? unit.kind : nothing,
        reason = reason,
    )
end

function _pqs_route_descriptor_available_diagnostic(
    descriptor,
    metrics = nothing;
    supported_terms = _PQS_PRODUCT_SOURCE_BOX_SHADOW_TERMS,
)
    route_kind = hasproperty(descriptor, :route_kind) ?
                 descriptor.route_kind : :legacy_route_units
    diagnostics = _pqs_route_descriptor_diagnostic_common(
        source = :pqs_pqs_product_route_descriptor_diagnostic,
        route_kind = route_kind,
        pqs_descriptor_count = 2,
        pqs_raw_plan_convertible_count = 2,
        product_doside_unit_count = 1,
        direct_or_support_body_piece_count = 0,
        descriptor_emitted = true,
        supported_terms = supported_terms,
        extra = (
            status = :descriptor_available,
            descriptor_available = true,
            descriptor_unavailable = false,
            metrics_supplied = !isnothing(metrics),
            descriptor_object_kind =
                hasproperty(descriptor, :object_kind) ?
                descriptor.object_kind : :legacy_route_units,
        ),
    )
    return (
        status = :descriptor_available,
        descriptor = descriptor,
        diagnostics = diagnostics,
    )
end

function _pqs_pqs_product_route_descriptor_diagnostic(
    route_like,
    metrics = nothing;
    route_name::Symbol = :pqs_pqs_product_source_box_safe_term_route,
    parent_dims = nothing,
    bond_axis = nothing,
    metadata = (;),
    provenance = (;),
    supported_terms = _PQS_PRODUCT_SOURCE_BOX_SHADOW_TERMS,
)
    selected_terms = _pqs_pqs_product_supported_safe_terms(supported_terms)
    missing = Symbol[]
    mismatches = Any[]
    pqs_plan_count = 0
    product_doside_unit_count = 0
    descriptor = nothing
    route_kind = hasproperty(route_like, :route_kind) ?
                 route_like.route_kind : nothing

    if !hasproperty(route_like, :units)
        push!(missing, :route_units)
    else
        units = route_like.units
        if hasproperty(units, :pqs_left)
            try
                left_raw_plan = _pqs_raw_product_box_plan_view(units.pqs_left)
                left_raw_plan.representation == :orthogonal_raw_product_box ||
                    push!(
                        mismatches,
                        _pqs_route_descriptor_unit_mismatch(
                            :pqs_left,
                            units.pqs_left,
                            :left_pqs_unit_is_not_orthogonal_raw_product_box,
                        ),
                    )
                pqs_plan_count += 1
            catch err
                push!(
                    mismatches,
                    _pqs_route_descriptor_unit_mismatch(
                        :pqs_left,
                        units.pqs_left,
                        :left_pqs_raw_plan_unavailable,
                    ),
                )
            end
        else
            push!(missing, :pqs_left)
        end

        if hasproperty(units, :pqs_right)
            try
                right_raw_plan = _pqs_raw_product_box_plan_view(units.pqs_right)
                right_raw_plan.representation == :orthogonal_raw_product_box ||
                    push!(
                        mismatches,
                        _pqs_route_descriptor_unit_mismatch(
                            :pqs_right,
                            units.pqs_right,
                            :right_pqs_unit_is_not_orthogonal_raw_product_box,
                        ),
                    )
                pqs_plan_count += 1
            catch err
                push!(
                    mismatches,
                    _pqs_route_descriptor_unit_mismatch(
                        :pqs_right,
                        units.pqs_right,
                        :right_pqs_raw_plan_unavailable,
                    ),
                )
            end
        else
            push!(missing, :pqs_right)
        end

        if hasproperty(units, :product)
            if units.product isa _CartesianNestedProductStagedByCenterUnit3D &&
               units.product.kind == :product_doside
                try
                    _require_product_doside_retained_block_unit(
                        units.product;
                        side = :product,
                    )
                    product_doside_unit_count = 1
                catch err
                    push!(
                        mismatches,
                        _pqs_route_descriptor_unit_mismatch(
                            :product,
                            units.product,
                            :product_doside_unit_metadata_invalid,
                        ),
                    )
                end
            else
                push!(
                    mismatches,
                    _pqs_route_descriptor_unit_mismatch(
                        :product,
                        units.product,
                        :product_unit_is_not_explicit_product_doside,
                    ),
                )
            end
        else
            push!(missing, :product)
        end
    end

    append!(
        missing,
        _pqs_route_descriptor_missing_symbols(
            pqs_plan_count,
            product_doside_unit_count,
        ),
    )
    unique!(missing)

    if isempty(missing) && isempty(mismatches)
        if hasproperty(route_like, :object_kind) &&
           route_like.object_kind == :pqs_pqs_product_safe_term_route_descriptor
            descriptor = route_like
        else
            units = route_like.units
            descriptor = _pqs_pqs_product_safe_term_route_descriptor(
                units.pqs_left,
                units.pqs_right,
                units.product;
                route_name,
                parent_dims,
                bond_axis,
                metadata,
                provenance,
                supported_terms = selected_terms,
            )
        end
        return _pqs_route_descriptor_available_diagnostic(
            descriptor,
            metrics;
            supported_terms = selected_terms,
        )
    end

    diagnostics = _pqs_route_descriptor_diagnostic_common(
        source = :pqs_pqs_product_route_descriptor_diagnostic,
        route_kind = route_kind,
        pqs_descriptor_count = pqs_plan_count,
        pqs_raw_plan_convertible_count = pqs_plan_count,
        product_doside_unit_count = product_doside_unit_count,
        direct_or_support_body_piece_count = 0,
        descriptor_emitted = false,
        supported_terms = selected_terms,
        extra = (
            status = :descriptor_unavailable,
            descriptor_available = false,
            descriptor_unavailable = true,
            metrics_supplied = !isnothing(metrics),
            route_has_units = hasproperty(route_like, :units),
            mismatch_count = length(mismatches),
            missing_count = length(missing),
        ),
    )
    return (
        status = :descriptor_unavailable,
        descriptor = nothing,
        missing = Tuple(missing),
        mismatches = Tuple(mismatches),
        diagnostics = diagnostics,
    )
end

function _pqs_route_staged_descriptor_signature(
    descriptor::_CartesianNestedProjectedQShellStagedUnitDescriptor3D,
)
    return (
        descriptor.current_box,
        descriptor.inner_box,
        descriptor.q,
        descriptor.L,
        descriptor.support_count,
        descriptor.mode_count,
        descriptor.retained_count,
    )
end

function _pqs_route_staged_descriptors(
    construction::_BondAlignedDiatomicHighOrderRecipeSourceConstruction3D,
)
    descriptors = _CartesianNestedProjectedQShellStagedUnitDescriptor3D[]
    seen = Set{Any}()
    function maybe_push_descriptor(descriptor)
        descriptor isa _CartesianNestedProjectedQShellStagedUnitDescriptor3D ||
            return nothing
        descriptor.kind == :projected_q_shell || return nothing
        signature = _pqs_route_staged_descriptor_signature(descriptor)
        signature in seen && return nothing
        push!(seen, signature)
        push!(descriptors, descriptor)
        return nothing
    end

    for layer in construction.shared_shell_layers
        try
            maybe_push_descriptor(_nested_projected_q_shell_staged_unit_descriptor(layer))
        catch err
            nothing
        end
    end
    for build in construction.region_builds
        metadata = build.metadata
        hasproperty(metadata, :pqs_staged_unit_descriptor) &&
            maybe_push_descriptor(metadata.pqs_staged_unit_descriptor)
    end
    return descriptors
end

function _pqs_route_raw_plan_convertibility(
    descriptors,
    construction::_BondAlignedDiatomicHighOrderRecipeSourceConstruction3D,
    metrics;
    orthogonality_atol::Real,
)
    checked = !isnothing(metrics)
    converted = 0
    failures = Any[]
    !checked && return (
        checked = false,
        converted_count = converted,
        failures = Tuple(failures),
    )

    for (index, descriptor) in pairs(descriptors)
        try
            source_mode_dims =
                ntuple(axis -> size(descriptor.axis_local_coefficients[axis], 2), 3)
            shared_plan = _cartesian_raw_product_box_plan(
                construction.axis_bundles,
                descriptor.axis_intervals,
                source_mode_dims;
                enforce_symmetric_odd = false,
            )
            raw_plan = _pqs_raw_product_box_plan(
                descriptor,
                shared_plan,
                metrics;
                orthogonality_atol,
            )
            _pqs_raw_product_box_plan_view(raw_plan)
            converted += 1
        catch err
            push!(
                failures,
                (
                    descriptor_index = index,
                    reason = :raw_pqs_plan_conversion_failed,
                    error_type = nameof(typeof(err)),
                    message = sprint(showerror, err),
                ),
            )
        end
    end
    return (
        checked = checked,
        converted_count = converted,
        failures = Tuple(failures),
    )
end

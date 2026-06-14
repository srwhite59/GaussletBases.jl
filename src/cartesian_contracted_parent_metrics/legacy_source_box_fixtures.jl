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

function _pqs_route_direct_or_support_body_mismatches(
    construction::_BondAlignedDiatomicHighOrderRecipeSourceConstruction3D,
)
    product_or_pqs_families = (
        :projected_q_shell,
        :shared_endcap_panel_shell_layer,
    )
    mismatches = Any[]
    for build in construction.region_builds
        build.primitive_family in product_or_pqs_families && continue
        push!(
            mismatches,
            (
                role = build.region_role,
                primitive_family = build.primitive_family,
                mapped_primitive = build.mapped_primitive,
                column_range = build.column_range,
                retained_count = build.retained_count,
                reason = :direct_or_support_piece_not_product_doside,
            ),
        )
    end
    return Tuple(mismatches)
end

_pqs_route_axis_symbol(axis::Int) = (:x, :y, :z)[axis]

function _pqs_route_policy_region_for_build(
    construction::_BondAlignedDiatomicHighOrderRecipeSourceConstruction3D,
    build,
)
    for region in construction.policy.construction_plan.regions
        region.order_index == build.region_order_index && return region
    end
    return nothing
end

function _pqs_route_realization_descriptor_for_build(
    construction::_BondAlignedDiatomicHighOrderRecipeSourceConstruction3D,
    build,
)
    for realization in construction.realization_audit.region_realizations
        realization.region_order_index == build.region_order_index &&
            return realization.realization_descriptor
    end
    return nothing
end

function _pqs_route_retained_unit_classification(build)
    build.primitive_family == :shared_endcap_panel_shell_layer &&
        return :already_product_doside
    build.primitive_family in (
        :contact_cap_owned_slab,
        :outer_mismatch_boundary_slab_set,
    ) && return :product_box_constructible
    build.primitive_family == :projected_q_shell && return :out_of_scope
    return :needs_direct_support_retained_unit_kind
end

function _pqs_route_coefficient_scope(classification::Symbol)
    classification == :already_product_doside &&
        return :product_doside_retained_unit
    classification == :product_box_constructible &&
        return :support_local_direct_rows
    classification == :needs_direct_support_retained_unit_kind &&
        return :support_local_direct_rows
    return :projected_q_shell_descriptor
end

function _pqs_route_safe_term_capability(classification::Symbol)
    classification == :already_product_doside &&
        return :product_doside_source_box_safe_terms
    classification == :product_box_constructible &&
        return :product_box_rule_available_not_instantiated
    classification == :needs_direct_support_retained_unit_kind &&
        return :support_local_reference_only
    return :not_body_retained_unit
end

function _pqs_route_slab_piece_rule(piece; slab_kind::Symbol, rule_kind::Symbol)
    fixed_axis_indices = findall(axis -> length(piece.box[axis]) == 1, 1:3)
    fixed_axis_index = length(fixed_axis_indices) == 1 ? only(fixed_axis_indices) : nothing
    active_axis_indices = Tuple(
        axis for axis in 1:3 if isnothing(fixed_axis_index) || axis != fixed_axis_index
    )
    return (
        rule_kind = rule_kind,
        slab_kind = slab_kind,
        piece_role = piece.role,
        piece_index = piece.piece_index,
        primitive_family = piece.primitive_family,
        source_box = piece.box,
        fixed_axis = isnothing(fixed_axis_index) ?
            nothing : _pqs_route_axis_symbol(fixed_axis_index),
        fixed_axis_index = fixed_axis_index,
        fixed_index = isnothing(fixed_axis_index) ?
            nothing : first(piece.box[fixed_axis_index]),
        active_axes =
            Tuple(_pqs_route_axis_symbol(axis) for axis in active_axis_indices),
        active_axis_indices = active_axis_indices,
        active_intervals = Tuple(piece.box[axis] for axis in active_axis_indices),
        active_transform_rule = :identity_selectors_preserve_direct_support,
        support_count = length(piece.support_indices),
        retained_count = length(piece.support_indices),
        product_doside_unit_created = false,
    )
end

function _pqs_route_product_box_construction_rule(build, descriptor)
    isnothing(descriptor) && return nothing
    piece_rules = Tuple(
        _pqs_route_slab_piece_rule(
            piece;
            slab_kind = build.region_role,
            rule_kind = build.primitive_family == :contact_cap_owned_slab ?
                :identity_selector_product_doside_slab :
                :identity_selector_product_doside_boundary_slab,
        )
        for piece in descriptor.pieces
    )
    build.primitive_family == :contact_cap_owned_slab && return merge(
        only(piece_rules),
        (
            slab_piece_count = 1,
            boundary_slab_set = false,
            product_doside_units_created = false,
            single_middle_product_body_unit = true,
        ),
    )
    build.primitive_family == :outer_mismatch_boundary_slab_set && return (
        rule_kind = :identity_selector_product_doside_boundary_slab_set,
        slab_kind = :outer_mismatch_shared_molecular_shell,
        boundary_slab_set = true,
        slab_piece_count = length(piece_rules),
        slab_piece_rules = piece_rules,
        product_doside_unit_created = false,
        product_doside_units_created = false,
        single_middle_product_body_unit = false,
        active_transform_rule = :identity_selectors_preserve_direct_support,
    )
    return nothing
end

function _pqs_route_retained_unit_fact(
    construction::_BondAlignedDiatomicHighOrderRecipeSourceConstruction3D,
    build;
    include_support_indices::Bool,
)
    region = _pqs_route_policy_region_for_build(construction, build)
    descriptor = _pqs_route_realization_descriptor_for_build(construction, build)
    classification = _pqs_route_retained_unit_classification(build)
    construction_rule = classification == :product_box_constructible ?
        _pqs_route_product_box_construction_rule(build, descriptor) :
        nothing
    return (
        object_kind = :pqs_route_retained_unit_fact,
        role = build.region_role,
        primitive_family = build.primitive_family,
        mapped_primitive = build.mapped_primitive,
        box = isnothing(region) ? nothing : region.box,
        inner_exclusion_box = isnothing(region) ? nothing : region.inner_exclusion_box,
        support_count = build.built_support_count,
        region_support_count = build.region_support_count,
        retained_count = build.retained_count,
        column_range = build.column_range,
        classification = classification,
        product_doside_unit = classification == :already_product_doside,
        raw_product_box_operator_contract =
            classification == :already_product_doside,
        product_box_construction_rule_available = !isnothing(construction_rule),
        coefficient_scope = _pqs_route_coefficient_scope(classification),
        safe_term_capability = _pqs_route_safe_term_capability(classification),
        construction_rule = construction_rule,
        support_indices = include_support_indices && !isnothing(region) ?
            copy(region.support_indices) : nothing,
        notes = classification == :out_of_scope ? (
            current_single_pqs_descriptor = true,
            body_retained_unit_vocabulary = false,
        ) : classification == :needs_direct_support_retained_unit_kind ? (
            direct_support_retained_unit_kind_needed = true,
            fixed_axis = nothing,
            product_doside_owned_unit_metadata = false,
        ) : (
            diagnostic_product_box_rule_only = true,
            product_doside_unit_created = false,
        ),
    )
end

function _pqs_route_product_box_slab_rule_count(unit_facts)
    total = 0
    for fact in unit_facts
        fact.classification == :product_box_constructible || continue
        rule = fact.construction_rule
        isnothing(rule) && continue
        total += hasproperty(rule, :slab_piece_count) ? rule.slab_piece_count : 1
    end
    return total
end

function _pqs_route_retained_unit_fact_audit(
    construction::_BondAlignedDiatomicHighOrderRecipeSourceConstruction3D;
    include_support_indices::Bool = false,
)
    unit_facts = Tuple(
        _pqs_route_retained_unit_fact(
            construction,
            build;
            include_support_indices,
        )
        for build in construction.region_builds
    )
    classification_counts = (
        already_product_doside = count(
            fact -> fact.classification == :already_product_doside,
            unit_facts,
        ),
        product_box_constructible = count(
            fact -> fact.classification == :product_box_constructible,
            unit_facts,
        ),
        needs_direct_support_retained_unit_kind = count(
            fact -> fact.classification == :needs_direct_support_retained_unit_kind,
            unit_facts,
        ),
        out_of_scope = count(
            fact -> fact.classification == :out_of_scope,
            unit_facts,
        ),
    )
    summary = (
        unit_count = length(unit_facts),
        classification_counts = classification_counts,
        product_doside_unit_count = classification_counts.already_product_doside,
        product_box_constructible_region_count =
            classification_counts.product_box_constructible,
        product_box_constructible_slab_rule_count =
            _pqs_route_product_box_slab_rule_count(unit_facts),
        direct_support_retained_unit_kind_count =
            classification_counts.needs_direct_support_retained_unit_kind,
        pqs_descriptor_region_count = classification_counts.out_of_scope,
    )
    diagnostics = (
        source = :pqs_route_retained_unit_fact_audit,
        private_diagnostic_only = true,
        descriptor_emitted = false,
        packet_adoption = false,
        fixed_block_construction_changed = false,
        qwhamiltonian_changed = false,
        sidecar_mutation = false,
        sidecar_installation = false,
        direct_support_reinterpreted_as_product_doside = false,
        retained_weight_semantics = :not_positive_quadrature_weights,
        ida_weight_division_allowed = false,
        local_ecp_gaussian_mwg_interaction_changed = false,
        include_support_indices = include_support_indices,
    )
    return (
        object_kind = :pqs_route_retained_unit_fact_audit,
        status = :audit_only,
        unit_facts = unit_facts,
        summary = summary,
        diagnostics = diagnostics,
    )
end

function _pqs_support_local_parent_coefficient_matrix(
    support_indices::AbstractVector{<:Integer},
    support_coefficients::AbstractMatrix{<:Real},
    parent_dimension::Int,
)
    nrows, ncols = size(support_coefficients)
    length(support_indices) == nrows || throw(
        DimensionMismatch("PQS support-local product/doside support rows must match support indices"),
    )
    row_indices = Int[]
    col_indices = Int[]
    values = Float64[]
    for col in 1:ncols, row in 1:nrows
        value = Float64(support_coefficients[row, col])
        iszero(value) && continue
        push!(row_indices, Int(support_indices[row]))
        push!(col_indices, col)
        push!(values, value)
    end
    return SparseArrays.sparse(
        row_indices,
        col_indices,
        values,
        parent_dimension,
        ncols,
    )
end

function _pqs_contact_cap_parent_coefficient_matrix(
    support_indices::AbstractVector{<:Integer},
    support_coefficients::AbstractMatrix{<:Real},
    parent_dimension::Int,
)
    return _pqs_support_local_parent_coefficient_matrix(
        support_indices,
        support_coefficients,
        parent_dimension,
    )
end

function _pqs_contact_cap_product_doside_unit(
    construction::_BondAlignedDiatomicHighOrderRecipeSourceConstruction3D;
    audit = _pqs_route_retained_unit_fact_audit(
        construction;
        include_support_indices = true,
    ),
)
    contact_facts = [
        fact for fact in audit.unit_facts if fact.role == :contact_cap
    ]
    length(contact_facts) == 1 || throw(
        ArgumentError("PQS contact-cap product/doside helper requires exactly one contact-cap fact"),
    )
    contact_fact = only(contact_facts)
    contact_fact.classification == :product_box_constructible || throw(
        ArgumentError("PQS contact-cap product/doside helper requires a product-box-constructible contact cap"),
    )
    rule = contact_fact.construction_rule
    isnothing(rule) && throw(
        ArgumentError("PQS contact-cap product/doside helper requires an explicit construction rule"),
    )
    rule.rule_kind == :identity_selector_product_doside_slab || throw(
        ArgumentError("PQS contact-cap product/doside helper requires an identity-selector slab rule"),
    )
    rule.slab_piece_count == 1 || throw(
        ArgumentError("PQS contact-cap product/doside helper requires a single slab rule"),
    )
    !isnothing(rule.fixed_axis_index) || throw(
        ArgumentError("PQS contact-cap product/doside helper requires one fixed axis"),
    )
    length(rule.active_axis_indices) == 2 || throw(
        ArgumentError("PQS contact-cap product/doside helper requires two active axes"),
    )
    active_lengths = Tuple(length(interval) for interval in rule.active_intervals)
    active_lengths[1] == active_lengths[2] || throw(
        ArgumentError("PQS contact-cap product/doside helper requires a q x q x 1 slab"),
    )

    contact_builds = [
        build for build in construction.region_builds if build.region_role == :contact_cap
    ]
    length(contact_builds) == 1 || throw(
        ArgumentError("PQS contact-cap product/doside helper requires exactly one contact-cap build"),
    )
    contact_build = only(contact_builds)
    descriptor =
        _pqs_route_realization_descriptor_for_build(construction, contact_build)
    isnothing(descriptor) && throw(
        ArgumentError("PQS contact-cap product/doside helper requires a realization descriptor"),
    )
    length(descriptor.pieces) == 1 || throw(
        ArgumentError("PQS contact-cap product/doside helper requires a single descriptor piece"),
    )
    piece = only(descriptor.pieces)
    support_indices = isnothing(contact_fact.support_indices) ?
        copy(piece.support_indices) : copy(contact_fact.support_indices)
    support_indices == piece.support_indices || throw(
        ArgumentError("PQS contact-cap audited support indices must match descriptor piece support"),
    )

    dims = _nested_axis_lengths(construction.axis_bundles)
    parent_dimension = prod(dims)
    support_states = [_cartesian_unflat_index(index, dims) for index in support_indices]
    expected_support_indices = Int[
        _cartesian_flat_index(state[1], state[2], state[3], dims)
        for state in support_states
    ]
    support_indices == expected_support_indices || throw(
        ArgumentError("PQS contact-cap support states do not round-trip through parent indexing"),
    )

    retained_count = contact_fact.retained_count
    length(contact_fact.column_range) == retained_count || throw(
        ArgumentError("PQS contact-cap column range must match retained count"),
    )
    retained_count == prod(active_lengths) || throw(
        ArgumentError("PQS contact-cap retained count must match active identity selector size"),
    )
    support_count = length(support_indices)
    support_count == retained_count || throw(
        ArgumentError("PQS contact-cap identity selector requires support count equal retained count"),
    )

    active_coefficients = ntuple(
        axis -> Matrix{Float64}(
            LinearAlgebra.I,
            active_lengths[axis],
            active_lengths[axis],
        ),
        2,
    )
    axes = ntuple(axis -> begin
        axis == rule.fixed_axis_index &&
            return _nested_product_staged_fixed_axis(rule.fixed_index)
        axis == rule.active_axis_indices[1] &&
            return _nested_product_staged_active_axis(
                rule.active_intervals[1],
                active_coefficients[1],
            )
        axis == rule.active_axis_indices[2] &&
            return _nested_product_staged_active_axis(
                rule.active_intervals[2],
                active_coefficients[2],
            )
        throw(ArgumentError("PQS contact-cap rule has inconsistent active/fixed axes"))
    end, 3)
    source_center_metadata =
        _product_doside_source_axis_center_metadata(construction, axes)
    local_coefficients = Matrix{Float64}(
        LinearAlgebra.I,
        support_count,
        retained_count,
    )
    unit = _CartesianNestedProductStagedByCenterUnit3D(
        rule.piece_role,
        :product_doside,
        contact_fact.column_range,
        support_indices,
        support_states,
        local_coefficients,
        axes,
        _nested_product_axis_function_indices(
            rule.fixed_axis_index,
            rule.active_axis_indices[1],
            active_lengths[1],
            rule.active_axis_indices[2],
            active_lengths[2],
        ),
        (
            source = :pqs_contact_cap_product_doside_unit,
            fact_role = contact_fact.role,
            primitive_family = contact_fact.primitive_family,
            mapped_primitive = contact_fact.mapped_primitive,
            source_axis_center_vectors =
                source_center_metadata.source_axis_center_vectors,
            source_center_convention =
                source_center_metadata.source_center_convention,
            source_center_status =
                source_center_metadata.source_center_status,
            source_center_source =
                source_center_metadata.source_center_source,
        ),
        (
            support_count = support_count,
            retained_count = retained_count,
            fixed_axis = rule.fixed_axis_index,
            fixed_index = rule.fixed_index,
            active_axes = rule.active_axis_indices,
            active_retained_counts = active_lengths,
            contact_cap_only = true,
            private_diagnostic_only = true,
            source_center_convention =
                source_center_metadata.source_center_convention,
            source_center_status =
                source_center_metadata.source_center_status,
        ),
    )

    parent_coefficients = _pqs_support_local_parent_coefficient_matrix(
        support_indices,
        local_coefficients,
        parent_dimension,
    )
    direct_coefficients = contact_build.built_object.coefficient_matrix
    parent_difference =
        Matrix{Float64}(parent_coefficients) - Matrix{Float64}(direct_coefficients)
    max_parent_coefficient_error =
        isempty(parent_difference) ? 0.0 : maximum(abs, parent_difference)
    equivalence = (
        support_indices_match =
            support_indices == contact_build.built_object.support_indices ==
            piece.support_indices,
        support_states_match =
            support_states == [_cartesian_unflat_index(index, dims) for index in support_indices],
        retained_count_match =
            retained_count == contact_build.retained_count == length(contact_fact.column_range),
        column_range_match =
            contact_fact.column_range == contact_build.column_range,
        coefficient_matrix_matches_direct_selector =
            max_parent_coefficient_error == 0.0,
        max_parent_coefficient_error = max_parent_coefficient_error,
    )
    diagnostics = (
        contact_cap_only = true,
        product_doside_unit_created = true,
        route_descriptor_emitted = false,
        construction_mutated = false,
        sidecar_installation = false,
        packet_adoption = false,
        fixed_block_construction_changed = false,
        qwhamiltonian_changed = false,
        ida_weight_division_allowed = false,
        retained_weight_semantics = :not_positive_quadrature_weights,
        product_box_construction_rule_available =
            contact_fact.product_box_construction_rule_available,
        input_fact_raw_product_box_operator_contract =
            contact_fact.raw_product_box_operator_contract,
        created_unit_raw_product_box_operator_contract = true,
    )
    return (
        object_kind = :pqs_contact_cap_product_doside_unit_fixture,
        status = :private_diagnostic_only,
        fact = contact_fact,
        unit = unit,
        equivalence = equivalence,
        diagnostics = diagnostics,
    )
end

function _pqs_outer_mismatch_product_doside_units(
    construction::_BondAlignedDiatomicHighOrderRecipeSourceConstruction3D;
    audit = _pqs_route_retained_unit_fact_audit(
        construction;
        include_support_indices = true,
    ),
)
    outer_facts = [
        fact for fact in audit.unit_facts if fact.role == :outer_mismatch_shared_molecular_shell
    ]
    length(outer_facts) == 1 || throw(
        ArgumentError("PQS outer-mismatch product/doside helper requires exactly one outer-mismatch fact"),
    )
    outer_fact = only(outer_facts)
    outer_fact.classification == :product_box_constructible || throw(
        ArgumentError("PQS outer-mismatch product/doside helper requires a product-box-constructible outer mismatch"),
    )
    rule = outer_fact.construction_rule
    isnothing(rule) && throw(
        ArgumentError("PQS outer-mismatch product/doside helper requires an explicit construction rule"),
    )
    rule.rule_kind == :identity_selector_product_doside_boundary_slab_set || throw(
        ArgumentError("PQS outer-mismatch product/doside helper requires a boundary slab-set rule"),
    )
    rule.boundary_slab_set || throw(
        ArgumentError("PQS outer-mismatch product/doside helper requires boundary_slab_set = true"),
    )
    rule.slab_piece_count == length(rule.slab_piece_rules) || throw(
        ArgumentError("PQS outer-mismatch slab-piece count must match its piece rules"),
    )

    outer_builds = [
        build for build in construction.region_builds if build.region_role == :outer_mismatch_shared_molecular_shell
    ]
    length(outer_builds) == 1 || throw(
        ArgumentError("PQS outer-mismatch product/doside helper requires exactly one outer-mismatch build"),
    )
    outer_build = only(outer_builds)
    descriptor =
        _pqs_route_realization_descriptor_for_build(construction, outer_build)
    isnothing(descriptor) && throw(
        ArgumentError("PQS outer-mismatch product/doside helper requires a realization descriptor"),
    )
    length(descriptor.pieces) == rule.slab_piece_count || throw(
        ArgumentError("PQS outer-mismatch descriptor piece count must match slab-piece count"),
    )

    dims = _nested_axis_lengths(construction.axis_bundles)
    parent_dimension = prod(dims)
    direct = outer_build.built_object
    size(direct.coefficient_matrix, 1) == parent_dimension || throw(
        ArgumentError("PQS outer-mismatch direct/support coefficient matrix has inconsistent parent dimension"),
    )
    size(direct.coefficient_matrix, 2) == outer_fact.retained_count || throw(
        ArgumentError("PQS outer-mismatch direct/support coefficient matrix has inconsistent retained dimension"),
    )

    units = _CartesianNestedProductStagedByCenterUnit3D[]
    piece_equivalences = NamedTuple[]
    parent_coefficient_pieces = SparseArrays.SparseMatrixCSC{Float64,Int}[]
    concatenated_support_indices = Int[]
    next_column = first(outer_fact.column_range)
    local_column_start = 1
    for (piece_rule, piece) in zip(rule.slab_piece_rules, descriptor.pieces)
        piece_rule.piece_index == piece.piece_index || throw(
            ArgumentError("PQS outer-mismatch slab rule order must match descriptor piece order"),
        )
        piece_rule.piece_role == piece.role || throw(
            ArgumentError("PQS outer-mismatch slab rule role must match descriptor piece role"),
        )
        piece_rule.primitive_family == piece.primitive_family || throw(
            ArgumentError("PQS outer-mismatch slab rule primitive must match descriptor piece primitive"),
        )
        !isnothing(piece_rule.fixed_axis_index) || throw(
            ArgumentError("PQS outer-mismatch product/doside slab requires one fixed axis"),
        )
        !isnothing(piece_rule.fixed_index) || throw(
            ArgumentError("PQS outer-mismatch product/doside slab requires one fixed index"),
        )
        length(piece_rule.active_axis_indices) == 2 || throw(
            ArgumentError("PQS outer-mismatch product/doside slab requires two active axes"),
        )
        length(piece_rule.active_intervals) == 2 || throw(
            ArgumentError("PQS outer-mismatch product/doside slab requires two active intervals"),
        )
        active_lengths = Tuple(length(interval) for interval in piece_rule.active_intervals)
        support_indices = copy(piece.support_indices)
        support_count = length(support_indices)
        support_count == piece_rule.support_count || throw(
            ArgumentError("PQS outer-mismatch slab support count must match its rule"),
        )
        support_count == piece_rule.retained_count || throw(
            ArgumentError("PQS outer-mismatch slab retained count must match support count"),
        )
        support_count == prod(active_lengths) || throw(
            ArgumentError("PQS outer-mismatch slab support count must equal active-axis product"),
        )
        piece_column_range = next_column:(next_column + support_count - 1)
        last(piece_column_range) <= last(outer_fact.column_range) || throw(
            ArgumentError("PQS outer-mismatch piece column range exceeds outer range"),
        )
        next_column = last(piece_column_range) + 1
        local_column_range = local_column_start:(local_column_start + support_count - 1)
        local_column_start = last(local_column_range) + 1

        support_states = [_cartesian_unflat_index(index, dims) for index in support_indices]
        expected_support_indices = Int[
            _cartesian_flat_index(state[1], state[2], state[3], dims)
            for state in support_states
        ]
        support_indices == expected_support_indices || throw(
            ArgumentError("PQS outer-mismatch support states do not round-trip through parent indexing"),
        )
        active_coefficients = ntuple(
            axis -> Matrix{Float64}(
                LinearAlgebra.I,
                active_lengths[axis],
                active_lengths[axis],
            ),
            2,
        )
        axes = ntuple(axis -> begin
            axis == piece_rule.fixed_axis_index &&
                return _nested_product_staged_fixed_axis(piece_rule.fixed_index)
            axis == piece_rule.active_axis_indices[1] &&
                return _nested_product_staged_active_axis(
                    piece_rule.active_intervals[1],
                    active_coefficients[1],
                )
            axis == piece_rule.active_axis_indices[2] &&
                return _nested_product_staged_active_axis(
                    piece_rule.active_intervals[2],
                    active_coefficients[2],
                )
            throw(ArgumentError("PQS outer-mismatch slab rule has inconsistent active/fixed axes"))
        end, 3)
        source_center_metadata =
            _product_doside_source_axis_center_metadata(construction, axes)
        local_coefficients = Matrix{Float64}(
            LinearAlgebra.I,
            support_count,
            support_count,
        )
        unit = _CartesianNestedProductStagedByCenterUnit3D(
            piece_rule.piece_role,
            :product_doside,
            piece_column_range,
            support_indices,
            support_states,
            local_coefficients,
            axes,
            _nested_product_axis_function_indices(
                piece_rule.fixed_axis_index,
                piece_rule.active_axis_indices[1],
                active_lengths[1],
                piece_rule.active_axis_indices[2],
                active_lengths[2],
            ),
            (
                source = :pqs_outer_mismatch_product_doside_units,
                fact_role = outer_fact.role,
                primitive_family = outer_fact.primitive_family,
                piece_role = piece_rule.piece_role,
                source_axis_center_vectors =
                    source_center_metadata.source_axis_center_vectors,
                source_center_convention =
                    source_center_metadata.source_center_convention,
                source_center_status =
                    source_center_metadata.source_center_status,
                source_center_source =
                    source_center_metadata.source_center_source,
            ),
            (
                support_count = support_count,
                retained_count = support_count,
                fixed_axis = piece_rule.fixed_axis_index,
                fixed_index = piece_rule.fixed_index,
                active_axes = piece_rule.active_axis_indices,
                active_retained_counts = active_lengths,
                outer_mismatch_piece = true,
                private_diagnostic_only = true,
                source_center_convention =
                    source_center_metadata.source_center_convention,
                source_center_status =
                    source_center_metadata.source_center_status,
            ),
        )
        parent_coefficients = _pqs_support_local_parent_coefficient_matrix(
            support_indices,
            local_coefficients,
            parent_dimension,
        )
        direct_piece_coefficients =
            Matrix{Float64}(direct.coefficient_matrix[:, local_column_range])
        parent_difference =
            Matrix{Float64}(parent_coefficients) - direct_piece_coefficients
        max_parent_coefficient_error =
            isempty(parent_difference) ? 0.0 : maximum(abs, parent_difference)
        push!(units, unit)
        push!(parent_coefficient_pieces, parent_coefficients)
        append!(concatenated_support_indices, support_indices)
        push!(
            piece_equivalences,
            (
                piece_role = piece_rule.piece_role,
                support_indices_match =
                    support_indices == direct.support_indices[local_column_range],
                retained_count_match = length(piece_column_range) == support_count,
                column_range = piece_column_range,
                local_column_range = local_column_range,
                coefficient_matrix_matches_direct_selector =
                    max_parent_coefficient_error == 0.0,
                max_parent_coefficient_error = max_parent_coefficient_error,
            ),
        )
    end
    next_column == last(outer_fact.column_range) + 1 || throw(
        ArgumentError("PQS outer-mismatch piece column ranges do not partition outer range"),
    )
    local_column_start == outer_fact.retained_count + 1 || throw(
        ArgumentError("PQS outer-mismatch local column ranges do not partition direct coefficients"),
    )
    concatenated_support_indices == direct.support_indices || throw(
        ArgumentError("PQS outer-mismatch concatenated support indices must match direct/support build"),
    )
    audited_support_set_match =
        isnothing(outer_fact.support_indices) ||
        sort(concatenated_support_indices) == sort(outer_fact.support_indices)
    audited_support_set_match || throw(
        ArgumentError("PQS outer-mismatch concatenated support indices must cover audited support"),
    )
    parent_coefficients = hcat(parent_coefficient_pieces...)
    parent_difference =
        Matrix{Float64}(parent_coefficients) -
        Matrix{Float64}(direct.coefficient_matrix)
    max_parent_coefficient_error =
        isempty(parent_difference) ? 0.0 : maximum(abs, parent_difference)
    aggregate_equivalence = (
        support_indices_match = concatenated_support_indices == direct.support_indices,
        audited_support_set_match = audited_support_set_match,
        retained_count_match =
            sum(length(unit.column_range) for unit in units) ==
            outer_build.retained_count == outer_fact.retained_count,
        column_range_partition =
            first(first(units).column_range) == first(outer_fact.column_range) &&
            last(last(units).column_range) == last(outer_fact.column_range) &&
            sum(length(unit.column_range) for unit in units) == length(outer_fact.column_range),
        coefficient_matrix_matches_direct_selector =
            max_parent_coefficient_error == 0.0,
        max_parent_coefficient_error = max_parent_coefficient_error,
    )
    diagnostics = (
        outer_mismatch_only = true,
        boundary_slab_set = true,
        product_doside_units_created = true,
        unit_count = length(units),
        slab_piece_count = rule.slab_piece_count,
        route_descriptor_emitted = false,
        construction_mutated = false,
        sidecar_installation = false,
        packet_adoption = false,
        fixed_block_construction_changed = false,
        qwhamiltonian_changed = false,
        ida_weight_division_allowed = false,
        retained_weight_semantics = :not_positive_quadrature_weights,
        input_fact_raw_product_box_operator_contract =
            outer_fact.raw_product_box_operator_contract,
        created_units_raw_product_box_operator_contract = true,
        descriptor_piece_order_defines_columns = true,
        audited_support_checked_as_set = true,
        product_box_construction_rule_available =
            outer_fact.product_box_construction_rule_available,
        local_ecp_gaussian_mwg_interaction_changed = false,
    )
    return (
        object_kind = :pqs_outer_mismatch_product_doside_units_fixture,
        status = :private_diagnostic_only,
        fact = outer_fact,
        units = Tuple(units),
        piece_equivalences = Tuple(piece_equivalences),
        aggregate_equivalence = aggregate_equivalence,
        diagnostics = diagnostics,
    )
end

const _PQS_CONTACT_CAP_SAFE_TERM_OPERATOR_TERMS =
    _PRODUCT_DOSIDE_SOURCE_BOX_REFERENCE_TERMS

const _PQS_CONTACT_CAP_UNSUPPORTED_OPERATOR_TERMS = (
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
)

function _pqs_contact_cap_axis_factor_terms(
    metrics::NamedTuple{(:x,:y,:z)},
    term::Symbol,
)
    term in _PQS_CONTACT_CAP_SAFE_TERM_OPERATOR_TERMS || throw(
        ArgumentError("PQS contact-cap safe-term comparison received unsupported term $(term)"),
    )
    return Tuple(
        ntuple(
            axis -> _product_doside_axis_metric_matrix(metrics, axis, factor_kinds[axis]),
            3,
        )
        for factor_kinds in _source_box_separable_term_factor_kinds(term)
    )
end

function _pqs_contact_cap_direct_support_oracle_entries(
    construction::_BondAlignedDiatomicHighOrderRecipeSourceConstruction3D,
    fixture,
)
    contact_builds = [
        build for build in construction.region_builds if build.region_role == :contact_cap
    ]
    length(contact_builds) == 1 || throw(
        ArgumentError("PQS contact-cap safe-term comparison requires exactly one contact-cap build"),
    )
    contact_build = only(contact_builds)
    direct = contact_build.built_object
    dims = _nested_axis_lengths(construction.axis_bundles)
    parent_dimension = prod(dims)
    support_indices = copy(direct.support_indices)
    support_indices == fixture.unit.support_indices || throw(
        ArgumentError("PQS contact-cap direct/support oracle support does not match product fixture"),
    )
    size(direct.coefficient_matrix, 1) == parent_dimension || throw(
        ArgumentError("PQS contact-cap direct/support coefficient matrix has inconsistent parent dimension"),
    )
    size(direct.coefficient_matrix, 2) == length(contact_build.column_range) || throw(
        ArgumentError("PQS contact-cap direct/support coefficient matrix has inconsistent retained dimension"),
    )
    contact_build.column_range == fixture.unit.column_range || throw(
        ArgumentError("PQS contact-cap direct/support column range does not match product fixture"),
    )
    support_states = [_cartesian_unflat_index(index, dims) for index in support_indices]
    local_direct_coefficients =
        Matrix{Float64}(direct.coefficient_matrix[support_indices, :])
    entries = _support_local_retained_entries(
        contact_build.column_range,
        support_states,
        local_direct_coefficients,
    )
    return (
        contact_build = contact_build,
        direct = direct,
        support_indices = support_indices,
        support_states = support_states,
        local_direct_coefficients = local_direct_coefficients,
        entries = entries,
    )
end

function _pqs_contact_cap_direct_support_oracle_block(
    direct_oracle,
    metrics::NamedTuple{(:x,:y,:z)},
    term::Symbol,
)
    return _fallback_staged_separable_sum_block(
        direct_oracle.entries,
        direct_oracle.entries,
        _pqs_contact_cap_axis_factor_terms(metrics, term),
    )
end

function _pqs_contact_cap_safe_term_operator_comparison(
    construction::_BondAlignedDiatomicHighOrderRecipeSourceConstruction3D,
    metrics::NamedTuple{(:x,:y,:z)};
    terms = _PQS_CONTACT_CAP_SAFE_TERM_OPERATOR_TERMS,
    atol::Real = 1.0e-12,
)
    selected_terms = Tuple(Symbol(term) for term in terms)
    !isempty(selected_terms) || throw(
        ArgumentError("PQS contact-cap safe-term comparison requires at least one term"),
    )
    for term in selected_terms
        term in _PQS_CONTACT_CAP_SAFE_TERM_OPERATOR_TERMS || throw(
            ArgumentError("PQS contact-cap safe-term comparison received unsupported term $(term)"),
        )
    end
    fixture = _pqs_contact_cap_product_doside_unit(construction)
    product_unit = fixture.unit
    direct_oracle =
        _pqs_contact_cap_direct_support_oracle_entries(construction, fixture)

    product_blocks = Dict{Symbol,Matrix{Float64}}()
    direct_oracle_blocks = Dict{Symbol,Matrix{Float64}}()
    term_errors = Dict{Symbol,Float64}()
    product_references = Dict{Symbol,Any}()
    output_finite = true
    shape_matches = true
    for term in selected_terms
        product_reference = _product_doside_source_box_reference_block(
            product_unit,
            product_unit,
            metrics;
            term,
            atol,
        )
        product_block = product_reference.block
        direct_block =
            _pqs_contact_cap_direct_support_oracle_block(direct_oracle, metrics, term)
        size(product_block) == size(direct_block) || (shape_matches = false)
        all(isfinite, product_block) || (output_finite = false)
        all(isfinite, direct_block) || (output_finite = false)
        block_error = LinearAlgebra.norm(product_block - direct_block, Inf)
        product_blocks[term] = product_block
        direct_oracle_blocks[term] = direct_block
        term_errors[term] = block_error
        product_references[term] = product_reference
    end
    shape_matches || throw(
        ArgumentError("PQS contact-cap safe-term comparison produced mismatched block shapes"),
    )
    output_finite || throw(
        ArgumentError("PQS contact-cap safe-term comparison produced non-finite entries"),
    )
    max_block_error = maximum(values(term_errors))
    max_block_error <= Float64(atol) || throw(
        ArgumentError("PQS contact-cap safe-term product and direct/support blocks disagree"),
    )

    diagnostics = (
        source = :pqs_contact_cap_safe_term_operator_comparison,
        contact_cap_only = true,
        private_diagnostic_only = true,
        terms_checked = selected_terms,
        supported_terms = _PQS_CONTACT_CAP_SAFE_TERM_OPERATOR_TERMS,
        unsupported_terms = _PQS_CONTACT_CAP_UNSUPPORTED_OPERATOR_TERMS,
        product_path = :_product_doside_source_box_reference_block,
        direct_oracle_path = :support_local_direct_selector_contract_pair_block,
        current_direct_support_selector_compared = true,
        product_doside_unit_created =
            fixture.diagnostics.product_doside_unit_created,
        input_fact_raw_product_box_operator_contract =
            fixture.diagnostics.input_fact_raw_product_box_operator_contract,
        created_unit_raw_product_box_operator_contract =
            fixture.diagnostics.created_unit_raw_product_box_operator_contract,
        product_box_construction_rule_available =
            fixture.diagnostics.product_box_construction_rule_available,
        route_descriptor_emitted = false,
        construction_mutated = false,
        sidecar_installation = false,
        packet_adoption = false,
        fixed_block_construction_changed = false,
        qwhamiltonian_changed = false,
        ida_weight_division_allowed = false,
        retained_weight_semantics = :not_positive_quadrature_weights,
        local_ecp_gaussian_mwg_interaction_changed = false,
        operator_factor_source = :explicit_metric_operator_data,
        operator_metric_sources = _cartesian_source_box_metric_sources(metrics),
        input_metric_operator_data = :caller_supplied_explicit_data,
        input_metric_operator_data_pgdg_checked = false,
        pgdg_analytic_operator_provenance_claimed = false,
        numerical_reference_fallback = false,
        product_source_box_reference_compared = true,
        direct_support_oracle_entries_built = true,
        direct_support_local_coefficient_shape =
            size(direct_oracle.local_direct_coefficients),
        retained_count = length(product_unit.column_range),
        support_count = length(product_unit.support_indices),
        column_range = product_unit.column_range,
        atol = Float64(atol),
        max_block_error = max_block_error,
        term_errors = term_errors,
        output_finite = output_finite,
    )
    return (
        object_kind = :pqs_contact_cap_safe_term_operator_comparison,
        status = :private_diagnostic_only,
        terms = selected_terms,
        fixture = fixture,
        product_blocks = product_blocks,
        direct_oracle_blocks = direct_oracle_blocks,
        product_references = product_references,
        term_errors = term_errors,
        max_block_error = max_block_error,
        diagnostics = diagnostics,
    )
end

const _PQS_OUTER_MISMATCH_SAFE_TERM_OPERATOR_TERMS =
    _PRODUCT_DOSIDE_SOURCE_BOX_REFERENCE_TERMS

const _PQS_OUTER_MISMATCH_UNSUPPORTED_OPERATOR_TERMS =
    _PQS_CONTACT_CAP_UNSUPPORTED_OPERATOR_TERMS

function _pqs_outer_mismatch_axis_factor_terms(
    metrics::NamedTuple{(:x,:y,:z)},
    term::Symbol,
)
    term in _PQS_OUTER_MISMATCH_SAFE_TERM_OPERATOR_TERMS || throw(
        ArgumentError("PQS outer-mismatch safe-term comparison received unsupported term $(term)"),
    )
    return Tuple(
        ntuple(
            axis -> _product_doside_axis_metric_matrix(metrics, axis, factor_kinds[axis]),
            3,
        )
        for factor_kinds in _source_box_separable_term_factor_kinds(term)
    )
end

function _pqs_outer_mismatch_direct_support_oracle_entries(
    construction::_BondAlignedDiatomicHighOrderRecipeSourceConstruction3D,
    fixture,
)
    outer_builds = [
        build for build in construction.region_builds if build.region_role == :outer_mismatch_shared_molecular_shell
    ]
    length(outer_builds) == 1 || throw(
        ArgumentError("PQS outer-mismatch safe-term comparison requires exactly one outer-mismatch build"),
    )
    outer_build = only(outer_builds)
    direct = outer_build.built_object
    dims = _nested_axis_lengths(construction.axis_bundles)
    parent_dimension = prod(dims)
    support_indices = copy(direct.support_indices)
    concatenated_support_indices =
        reduce(vcat, (unit.support_indices for unit in fixture.units); init = Int[])
    support_indices == concatenated_support_indices || throw(
        ArgumentError("PQS outer-mismatch direct/support oracle support does not match product fixtures"),
    )
    size(direct.coefficient_matrix, 1) == parent_dimension || throw(
        ArgumentError("PQS outer-mismatch direct/support coefficient matrix has inconsistent parent dimension"),
    )
    size(direct.coefficient_matrix, 2) == length(outer_build.column_range) || throw(
        ArgumentError("PQS outer-mismatch direct/support coefficient matrix has inconsistent retained dimension"),
    )
    outer_build.column_range == fixture.fact.column_range || throw(
        ArgumentError("PQS outer-mismatch direct/support column range does not match product fixture"),
    )
    support_states = [_cartesian_unflat_index(index, dims) for index in support_indices]
    local_direct_coefficients =
        Matrix{Float64}(direct.coefficient_matrix[support_indices, :])
    entries = _support_local_retained_entries(
        outer_build.column_range,
        support_states,
        local_direct_coefficients,
    )
    return (
        outer_build = outer_build,
        direct = direct,
        support_indices = support_indices,
        support_states = support_states,
        local_direct_coefficients = local_direct_coefficients,
        entries = entries,
    )
end

function _pqs_outer_mismatch_direct_support_oracle_block(
    direct_oracle,
    metrics::NamedTuple{(:x,:y,:z)},
    term::Symbol,
)
    return _fallback_staged_separable_sum_block(
        direct_oracle.entries,
        direct_oracle.entries,
        _pqs_outer_mismatch_axis_factor_terms(metrics, term),
    )
end

function _pqs_outer_mismatch_local_column_range(
    column_range::UnitRange{Int},
    full_range::UnitRange{Int},
)
    first(column_range) >= first(full_range) &&
        last(column_range) <= last(full_range) || throw(
            ArgumentError("PQS outer-mismatch unit column range must lie inside full range"),
        )
    return (first(column_range) - first(full_range) + 1):(
        last(column_range) - first(full_range) + 1
    )
end

function _pqs_outer_mismatch_product_block(
    fixture,
    metrics::NamedTuple{(:x,:y,:z)};
    term::Symbol,
    atol::Real,
)
    retained_dimension = length(fixture.fact.column_range)
    product_block = zeros(Float64, retained_dimension, retained_dimension)
    pair_references = NamedTuple[]
    cross_slab_pair_count = 0
    for left_unit in fixture.units
        left_range = _pqs_outer_mismatch_local_column_range(
            left_unit.column_range,
            fixture.fact.column_range,
        )
        for right_unit in fixture.units
            right_range = _pqs_outer_mismatch_local_column_range(
                right_unit.column_range,
                fixture.fact.column_range,
            )
            product_reference = _product_doside_source_box_reference_block(
                left_unit,
                right_unit,
                metrics;
                term,
                atol,
            )
            product_block[left_range, right_range] .= product_reference.block
            left_unit.role != right_unit.role && (cross_slab_pair_count += 1)
            push!(
                pair_references,
                (
                    left_role = left_unit.role,
                    right_role = right_unit.role,
                    left_column_range = left_unit.column_range,
                    right_column_range = right_unit.column_range,
                    local_left_range = left_range,
                    local_right_range = right_range,
                    block_shape = size(product_reference.block),
                    reference = product_reference,
                ),
            )
        end
    end
    return (
        block = product_block,
        pair_references = Tuple(pair_references),
        pair_block_count = length(pair_references),
        cross_slab_pair_count = cross_slab_pair_count,
    )
end

function _pqs_outer_mismatch_safe_term_operator_comparison(
    construction::_BondAlignedDiatomicHighOrderRecipeSourceConstruction3D,
    metrics::NamedTuple{(:x,:y,:z)};
    terms = _PQS_OUTER_MISMATCH_SAFE_TERM_OPERATOR_TERMS,
    atol::Real = 1.0e-12,
)
    selected_terms = Tuple(Symbol(term) for term in terms)
    !isempty(selected_terms) || throw(
        ArgumentError("PQS outer-mismatch safe-term comparison requires at least one term"),
    )
    for term in selected_terms
        term in _PQS_OUTER_MISMATCH_SAFE_TERM_OPERATOR_TERMS || throw(
            ArgumentError("PQS outer-mismatch safe-term comparison received unsupported term $(term)"),
        )
    end
    fixture = _pqs_outer_mismatch_product_doside_units(construction)
    direct_oracle =
        _pqs_outer_mismatch_direct_support_oracle_entries(construction, fixture)

    product_blocks = Dict{Symbol,Matrix{Float64}}()
    direct_oracle_blocks = Dict{Symbol,Matrix{Float64}}()
    product_references = Dict{Symbol,Any}()
    term_errors = Dict{Symbol,Float64}()
    pair_block_counts = Dict{Symbol,Int}()
    cross_slab_pair_counts = Dict{Symbol,Int}()
    output_finite = true
    shape_matches = true
    retained_dimension = length(fixture.fact.column_range)
    for term in selected_terms
        product_result = _pqs_outer_mismatch_product_block(
            fixture,
            metrics;
            term,
            atol,
        )
        product_block = product_result.block
        direct_block =
            _pqs_outer_mismatch_direct_support_oracle_block(direct_oracle, metrics, term)
        size(product_block) == size(direct_block) || (shape_matches = false)
        size(product_block) == (retained_dimension, retained_dimension) ||
            (shape_matches = false)
        all(isfinite, product_block) || (output_finite = false)
        all(isfinite, direct_block) || (output_finite = false)
        block_error = LinearAlgebra.norm(product_block - direct_block, Inf)
        product_blocks[term] = product_block
        direct_oracle_blocks[term] = direct_block
        product_references[term] = product_result
        term_errors[term] = block_error
        pair_block_counts[term] = product_result.pair_block_count
        cross_slab_pair_counts[term] = product_result.cross_slab_pair_count
    end
    shape_matches || throw(
        ArgumentError("PQS outer-mismatch safe-term comparison produced mismatched block shapes"),
    )
    output_finite || throw(
        ArgumentError("PQS outer-mismatch safe-term comparison produced non-finite entries"),
    )
    max_block_error = maximum(values(term_errors))
    max_block_error <= Float64(atol) || throw(
        ArgumentError("PQS outer-mismatch safe-term product and direct/support blocks disagree"),
    )
    cross_slab_blocks_included =
        all(count -> count > 0, values(cross_slab_pair_counts))

    diagnostics = (
        source = :pqs_outer_mismatch_safe_term_operator_comparison,
        outer_mismatch_only = true,
        boundary_slab_set = true,
        private_diagnostic_only = true,
        terms_checked = selected_terms,
        supported_terms = _PQS_OUTER_MISMATCH_SAFE_TERM_OPERATOR_TERMS,
        unsupported_terms = _PQS_OUTER_MISMATCH_UNSUPPORTED_OPERATOR_TERMS,
        product_path = :_product_doside_source_box_reference_block,
        direct_oracle_path = :support_local_direct_selector_contract_pair_block,
        product_doside_units_created =
            fixture.diagnostics.product_doside_units_created,
        complete_slab_set_block_assembled = true,
        cross_slab_blocks_included = cross_slab_blocks_included,
        direct_support_oracle_compared = true,
        route_descriptor_emitted = false,
        construction_mutated = false,
        sidecar_installation = false,
        packet_adoption = false,
        fixed_block_construction_changed = false,
        qwhamiltonian_changed = false,
        ida_weight_division_allowed = false,
        retained_weight_semantics = :not_positive_quadrature_weights,
        local_ecp_gaussian_mwg_interaction_changed = false,
        operator_factor_source = :explicit_metric_operator_data,
        operator_metric_sources = _cartesian_source_box_metric_sources(metrics),
        input_metric_operator_data = :caller_supplied_explicit_data,
        input_metric_operator_data_pgdg_checked = false,
        pgdg_analytic_operator_provenance_claimed = false,
        numerical_reference_fallback = false,
        product_source_box_reference_compared = true,
        direct_support_oracle_entries_built = true,
        direct_support_local_coefficient_shape =
            size(direct_oracle.local_direct_coefficients),
        retained_count = retained_dimension,
        support_count = length(direct_oracle.support_indices),
        unit_count = length(fixture.units),
        pair_block_counts = pair_block_counts,
        cross_slab_pair_counts = cross_slab_pair_counts,
        atol = Float64(atol),
        max_block_error = max_block_error,
        term_errors = term_errors,
        output_finite = output_finite,
    )
    return (
        object_kind = :pqs_outer_mismatch_safe_term_operator_comparison,
        status = :private_diagnostic_only,
        terms = selected_terms,
        fixture = fixture,
        product_blocks = product_blocks,
        direct_oracle_blocks = direct_oracle_blocks,
        product_references = product_references,
        term_errors = term_errors,
        max_block_error = max_block_error,
        diagnostics = diagnostics,
    )
end

function _pqs_atom_box_support_dense_units(
    construction::_BondAlignedDiatomicHighOrderRecipeSourceConstruction3D;
    audit = _pqs_route_retained_unit_fact_audit(
        construction;
        include_support_indices = true,
    ),
)
    atom_roles = (:left_atom_box, :right_atom_box)
    facts_by_role = Dict(fact.role => fact for fact in audit.unit_facts)
    builds_by_role = Dict(build.region_role => build for build in construction.region_builds)
    dims = _nested_axis_lengths(construction.axis_bundles)
    parent_dimension = prod(dims)

    units = _CartesianNestedProductStagedByCenterUnit3D[]
    equivalences = NamedTuple[]
    local_identity_errors = Dict{Symbol,Float64}()
    max_parent_coefficient_error = 0.0
    for role in atom_roles
        haskey(facts_by_role, role) || throw(
            ArgumentError("PQS atom-box support-dense helper requires $(role) fact"),
        )
        haskey(builds_by_role, role) || throw(
            ArgumentError("PQS atom-box support-dense helper requires $(role) build"),
        )
        fact = facts_by_role[role]
        build = builds_by_role[role]
        fact.primitive_family == :atom_local_complete_shell_sequence || throw(
            ArgumentError("PQS atom-box support-dense helper requires atom-local complete-shell sequence facts"),
        )
        build.primitive_family == :atom_local_complete_shell_sequence || throw(
            ArgumentError("PQS atom-box support-dense helper requires atom-local complete-shell sequence builds"),
        )
        fact.classification == :needs_direct_support_retained_unit_kind || throw(
            ArgumentError("PQS atom-box support-dense helper requires direct/support atom-box classification"),
        )
        fact.safe_term_capability == :support_local_reference_only || throw(
            ArgumentError("PQS atom-box support-dense helper requires support-local reference capability"),
        )
        direct = build.built_object
        size(direct.coefficient_matrix, 1) == parent_dimension || throw(
            ArgumentError("PQS atom-box direct/support coefficient matrix has inconsistent parent dimension"),
        )
        size(direct.coefficient_matrix, 2) == build.retained_count || throw(
            ArgumentError("PQS atom-box direct/support coefficient matrix has inconsistent retained dimension"),
        )
        support_indices = copy(direct.support_indices)
        support_count = length(support_indices)
        support_count == build.built_support_count == fact.support_count || throw(
            ArgumentError("PQS atom-box support count mismatch"),
        )
        build.retained_count == fact.retained_count || throw(
            ArgumentError("PQS atom-box retained count mismatch"),
        )
        length(build.column_range) == build.retained_count || throw(
            ArgumentError("PQS atom-box column range must match retained count"),
        )
        build.column_range == fact.column_range || throw(
            ArgumentError("PQS atom-box column range mismatch"),
        )
        audited_support_set_match =
            isnothing(fact.support_indices) ||
            sort(support_indices) == sort(fact.support_indices)
        audited_support_set_match || throw(
            ArgumentError("PQS atom-box support indices must cover audited support"),
        )
        support_states = [_cartesian_unflat_index(index, dims) for index in support_indices]
        local_coefficients =
            Matrix{Float64}(direct.coefficient_matrix[support_indices, :])
        identity_matrix = Matrix{Float64}(
            LinearAlgebra.I,
            size(local_coefficients, 1),
            size(local_coefficients, 2),
        )
        local_identity_error =
            LinearAlgebra.norm(local_coefficients - identity_matrix, Inf)
        support_count_equals_retained_count = support_count == build.retained_count
        direct_support_unit =
            support_count_equals_retained_count && local_identity_error <= 1.0e-12
        local_identity_errors[role] = local_identity_error
        axes = ntuple(_axis -> _nested_product_staged_fixed_axis(1), 3)
        axis_function_indices = fill((1, 1, 1), length(build.column_range))
        unit = _CartesianNestedProductStagedByCenterUnit3D(
            role,
            :support_dense,
            build.column_range,
            support_indices,
            support_states,
            local_coefficients,
            axes,
            axis_function_indices,
            (
                source = :pqs_atom_box_support_dense_units,
                fact_role = fact.role,
                primitive_family = fact.primitive_family,
                mapped_primitive = fact.mapped_primitive,
            ),
            (
                atom_box_only = true,
                support_count = support_count,
                retained_count = build.retained_count,
                coefficient_shape = size(local_coefficients),
                support_count_equals_retained_count =
                    support_count_equals_retained_count,
                support_dense_direct_support_unit = direct_support_unit,
                support_local_retained_unit = true,
                support_local_reference_only = true,
                product_doside_unit = false,
                raw_product_box_operator_contract = false,
                local_identity_error = local_identity_error,
                private_diagnostic_only = true,
            ),
        )
        parent_coefficients = _pqs_support_local_parent_coefficient_matrix(
            support_indices,
            local_coefficients,
            parent_dimension,
        )
        parent_difference =
            Matrix{Float64}(parent_coefficients) -
            Matrix{Float64}(direct.coefficient_matrix)
        parent_coefficient_error =
            isempty(parent_difference) ? 0.0 : maximum(abs, parent_difference)
        max_parent_coefficient_error =
            max(max_parent_coefficient_error, parent_coefficient_error)
        push!(units, unit)
        push!(
            equivalences,
            (
                role = role,
                primitive_family = fact.primitive_family,
                support_indices_match = support_indices == direct.support_indices,
                audited_support_set_match = audited_support_set_match,
                support_states_match =
                    support_states == [_cartesian_unflat_index(index, dims) for index in support_indices],
                retained_count_match =
                    length(build.column_range) == build.retained_count ==
                    fact.retained_count == size(local_coefficients, 2),
                support_count_equals_retained_count =
                    support_count_equals_retained_count,
                column_range = build.column_range,
                coefficient_matrix_matches_direct_support =
                    parent_coefficient_error == 0.0,
                coefficient_matrix_matches_direct_parent_coefficients =
                    parent_coefficient_error == 0.0,
                max_parent_coefficient_error = parent_coefficient_error,
                local_identity_error = local_identity_error,
                local_support_coefficient_shape = size(local_coefficients),
                support_local_retained_unit = true,
                direct_support_unit = direct_support_unit,
            ),
        )
    end
    direct_support_units_created =
        all(unit -> unit.diagnostics.support_dense_direct_support_unit, units)
    support_count_equals_retained_count_by_role = Dict(
        unit.role => unit.diagnostics.support_count_equals_retained_count
        for unit in units
    )
    local_coefficient_shapes = Dict(
        unit.role => unit.diagnostics.coefficient_shape for unit in units
    )

    diagnostics = (
        source = :pqs_atom_box_support_dense_units,
        atom_box_only = true,
        support_dense_direct_support_units_created = direct_support_units_created,
        support_local_retained_units_created = true,
        product_doside_units_created = false,
        raw_product_box_operator_contract = false,
        support_local_reference_only = true,
        product_box_construction_rule_available = false,
        route_descriptor_emitted = false,
        construction_mutated = false,
        sidecar_installation = false,
        packet_adoption = false,
        fixed_block_construction_changed = false,
        qwhamiltonian_changed = false,
        ida_weight_division_allowed = false,
        retained_weight_semantics = :not_positive_quadrature_weights,
        local_ecp_gaussian_mwg_interaction_changed = false,
        unit_count = length(units),
        roles = atom_roles,
        max_parent_coefficient_error = max_parent_coefficient_error,
        local_identity_errors = local_identity_errors,
        local_coefficient_shapes = local_coefficient_shapes,
        support_count_equals_retained_count_by_role =
            support_count_equals_retained_count_by_role,
        local_identity_is_product_box_claim = false,
        safe_term_operator_comparison_added = false,
    )
    return (
        object_kind = :pqs_atom_box_support_dense_units_fixture,
        status = :private_diagnostic_only,
        units = Tuple(units),
        equivalences = Tuple(equivalences),
        diagnostics = diagnostics,
    )
end

const _PQS_ATOM_BOX_SAFE_TERM_OPERATOR_TERMS =
    _PRODUCT_DOSIDE_SOURCE_BOX_REFERENCE_TERMS

const _PQS_ATOM_BOX_UNSUPPORTED_OPERATOR_TERMS =
    _PQS_CONTACT_CAP_UNSUPPORTED_OPERATOR_TERMS

function _pqs_atom_box_axis_factor_terms(
    metrics::NamedTuple{(:x,:y,:z)},
    term::Symbol,
)
    term in _PQS_ATOM_BOX_SAFE_TERM_OPERATOR_TERMS || throw(
        ArgumentError("PQS atom-box safe-term comparison received unsupported term $(term)"),
    )
    return Tuple(
        ntuple(
            axis -> _product_doside_axis_metric_matrix(metrics, axis, factor_kinds[axis]),
            3,
        )
        for factor_kinds in _source_box_separable_term_factor_kinds(term)
    )
end

function _pqs_atom_box_local_column_range(
    column_range::UnitRange{Int},
    full_range::UnitRange{Int},
)
    first(column_range) >= first(full_range) &&
        last(column_range) <= last(full_range) || throw(
            ArgumentError("PQS atom-box unit column range must lie inside full range"),
        )
    return (first(column_range) - first(full_range) + 1):(
        last(column_range) - first(full_range) + 1
    )
end

function _pqs_atom_box_support_dense_unit_entries(fixture)
    return Dict(
        unit.role => _support_local_retained_entries(
            unit.column_range,
            unit.support_states,
            unit.coefficient_matrix,
        )
        for unit in fixture.units
    )
end

function _pqs_atom_box_direct_support_oracle_entries(
    construction::_BondAlignedDiatomicHighOrderRecipeSourceConstruction3D,
    fixture,
)
    atom_roles = Tuple(unit.role for unit in fixture.units)
    builds_by_role = Dict(build.region_role => build for build in construction.region_builds)
    dims = _nested_axis_lengths(construction.axis_bundles)
    role_entries = Dict{Symbol,Vector{Vector{_ParentCoefficientEntry3D}}}()
    direct_support_shapes = Dict{Symbol,Tuple{Int,Int}}()
    for unit in fixture.units
        role = unit.role
        haskey(builds_by_role, role) || throw(
            ArgumentError("PQS atom-box direct/support oracle requires $(role) build"),
        )
        build = builds_by_role[role]
        build.primitive_family == :atom_local_complete_shell_sequence || throw(
            ArgumentError("PQS atom-box direct/support oracle requires atom-local complete-shell builds"),
        )
        direct = build.built_object
        direct.support_indices == unit.support_indices || throw(
            ArgumentError("PQS atom-box direct/support oracle support does not match support-dense fixture"),
        )
        build.column_range == unit.column_range || throw(
            ArgumentError("PQS atom-box direct/support oracle column range does not match support-dense fixture"),
        )
        support_states = [_cartesian_unflat_index(index, dims) for index in direct.support_indices]
        support_states == unit.support_states || throw(
            ArgumentError("PQS atom-box direct/support oracle support states do not match support-dense fixture"),
        )
        local_direct_coefficients =
            Matrix{Float64}(direct.coefficient_matrix[direct.support_indices, :])
        role_entries[role] = _support_local_retained_entries(
            build.column_range,
            support_states,
            local_direct_coefficients,
        )
        direct_support_shapes[role] = size(local_direct_coefficients)
    end
    return (
        atom_roles = atom_roles,
        entries_by_role = role_entries,
        direct_support_shapes = direct_support_shapes,
    )
end

function _pqs_atom_box_support_local_block(
    units,
    entries_by_role,
    full_range::UnitRange{Int},
    metrics::NamedTuple{(:x,:y,:z)};
    term::Symbol,
)
    retained_dimension = length(full_range)
    block = zeros(Float64, retained_dimension, retained_dimension)
    pair_blocks = NamedTuple[]
    cross_atom_pair_count = 0
    axis_factor_terms = _pqs_atom_box_axis_factor_terms(metrics, term)
    for left_unit in units
        left_range =
            _pqs_atom_box_local_column_range(left_unit.column_range, full_range)
        for right_unit in units
            right_range =
                _pqs_atom_box_local_column_range(right_unit.column_range, full_range)
            pair_block = _fallback_staged_separable_sum_block(
                entries_by_role[left_unit.role],
                entries_by_role[right_unit.role],
                axis_factor_terms,
            )
            block[left_range, right_range] .= pair_block
            left_unit.role != right_unit.role && (cross_atom_pair_count += 1)
            push!(
                pair_blocks,
                (
                    left_role = left_unit.role,
                    right_role = right_unit.role,
                    left_column_range = left_unit.column_range,
                    right_column_range = right_unit.column_range,
                    local_left_range = left_range,
                    local_right_range = right_range,
                    block_shape = size(pair_block),
                ),
            )
        end
    end
    return (
        block = block,
        pair_blocks = Tuple(pair_blocks),
        pair_block_count = length(pair_blocks),
        cross_atom_pair_count = cross_atom_pair_count,
    )
end

function _pqs_atom_box_safe_term_operator_comparison(
    construction::_BondAlignedDiatomicHighOrderRecipeSourceConstruction3D,
    metrics::NamedTuple{(:x,:y,:z)};
    terms = _PQS_ATOM_BOX_SAFE_TERM_OPERATOR_TERMS,
    atol::Real = 1.0e-12,
)
    selected_terms = Tuple(Symbol(term) for term in terms)
    !isempty(selected_terms) || throw(
        ArgumentError("PQS atom-box safe-term comparison requires at least one term"),
    )
    for term in selected_terms
        term in _PQS_ATOM_BOX_SAFE_TERM_OPERATOR_TERMS || throw(
            ArgumentError("PQS atom-box safe-term comparison received unsupported term $(term)"),
        )
    end
    fixture = _pqs_atom_box_support_dense_units(construction)
    units = fixture.units
    first(first(units).column_range) <= last(last(units).column_range) || throw(
        ArgumentError("PQS atom-box safe-term comparison requires ordered atom-box unit ranges"),
    )
    full_range = first(first(units).column_range):last(last(units).column_range)
    sum(length(unit.column_range) for unit in units) == length(full_range) || throw(
        ArgumentError("PQS atom-box safe-term comparison requires contiguous atom-box ranges"),
    )
    unit_entries = _pqs_atom_box_support_dense_unit_entries(fixture)
    direct_oracle =
        _pqs_atom_box_direct_support_oracle_entries(construction, fixture)

    support_dense_blocks = Dict{Symbol,Matrix{Float64}}()
    direct_oracle_blocks = Dict{Symbol,Matrix{Float64}}()
    support_dense_references = Dict{Symbol,Any}()
    direct_oracle_references = Dict{Symbol,Any}()
    term_errors = Dict{Symbol,Float64}()
    pair_block_counts = Dict{Symbol,Int}()
    cross_atom_pair_counts = Dict{Symbol,Int}()
    output_finite = true
    shape_matches = true
    retained_dimension = length(full_range)
    for term in selected_terms
        support_dense_result = _pqs_atom_box_support_local_block(
            units,
            unit_entries,
            full_range,
            metrics;
            term,
        )
        direct_result = _pqs_atom_box_support_local_block(
            units,
            direct_oracle.entries_by_role,
            full_range,
            metrics;
            term,
        )
        support_dense_block = support_dense_result.block
        direct_block = direct_result.block
        size(support_dense_block) == size(direct_block) || (shape_matches = false)
        size(support_dense_block) == (retained_dimension, retained_dimension) ||
            (shape_matches = false)
        all(isfinite, support_dense_block) || (output_finite = false)
        all(isfinite, direct_block) || (output_finite = false)
        block_error = LinearAlgebra.norm(support_dense_block - direct_block, Inf)
        support_dense_blocks[term] = support_dense_block
        direct_oracle_blocks[term] = direct_block
        support_dense_references[term] = support_dense_result
        direct_oracle_references[term] = direct_result
        term_errors[term] = block_error
        pair_block_counts[term] = support_dense_result.pair_block_count
        cross_atom_pair_counts[term] = support_dense_result.cross_atom_pair_count
    end
    shape_matches || throw(
        ArgumentError("PQS atom-box safe-term comparison produced mismatched block shapes"),
    )
    output_finite || throw(
        ArgumentError("PQS atom-box safe-term comparison produced non-finite entries"),
    )
    max_block_error = maximum(values(term_errors))
    max_block_error <= Float64(atol) || throw(
        ArgumentError("PQS atom-box support-dense and direct/support blocks disagree"),
    )
    cross_atom_blocks_included =
        all(count -> count > 0, values(cross_atom_pair_counts))

    diagnostics = (
        source = :pqs_atom_box_safe_term_operator_comparison,
        atom_box_only = true,
        support_dense_direct_support_units_created =
            fixture.diagnostics.support_dense_direct_support_units_created,
        support_local_fallback_operator_comparison = true,
        product_doside_units_created = false,
        raw_product_box_operator_contract = false,
        product_box_construction_rule_available = false,
        complete_atom_box_block_assembled = true,
        cross_atom_blocks_included = cross_atom_blocks_included,
        direct_support_oracle_compared = true,
        route_descriptor_emitted = false,
        construction_mutated = false,
        sidecar_installation = false,
        packet_adoption = false,
        fixed_block_construction_changed = false,
        qwhamiltonian_changed = false,
        ida_weight_division_allowed = false,
        retained_weight_semantics = :not_positive_quadrature_weights,
        local_ecp_gaussian_mwg_interaction_changed = false,
        operator_factor_source = :explicit_metric_operator_data,
        operator_metric_sources = _cartesian_source_box_metric_sources(metrics),
        input_metric_operator_data = :caller_supplied_explicit_data,
        input_metric_operator_data_pgdg_checked = false,
        pgdg_analytic_operator_provenance_claimed = false,
        numerical_reference_fallback = false,
        supported_terms = _PQS_ATOM_BOX_SAFE_TERM_OPERATOR_TERMS,
        unsupported_terms = _PQS_ATOM_BOX_UNSUPPORTED_OPERATOR_TERMS,
        terms_checked = selected_terms,
        retained_count = retained_dimension,
        support_count = sum(length(unit.support_indices) for unit in units),
        unit_count = length(units),
        pair_block_counts = pair_block_counts,
        cross_atom_pair_counts = cross_atom_pair_counts,
        direct_support_shapes = direct_oracle.direct_support_shapes,
        atol = Float64(atol),
        max_block_error = max_block_error,
        term_errors = term_errors,
        output_finite = output_finite,
    )
    return (
        object_kind = :pqs_atom_box_safe_term_operator_comparison,
        status = :private_diagnostic_only,
        terms = selected_terms,
        fixture = fixture,
        support_dense_blocks = support_dense_blocks,
        direct_oracle_blocks = direct_oracle_blocks,
        support_dense_references = support_dense_references,
        direct_oracle_references = direct_oracle_references,
        term_errors = term_errors,
        max_block_error = max_block_error,
        diagnostics = diagnostics,
    )
end

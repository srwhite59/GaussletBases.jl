module CartesianContractedParentMetrics

import LinearAlgebra
import SparseArrays

import ..GaussletBases: _NestedFixedBlock3D,
                         _require_analytic_primitive_backend,
                         centers,
                         contract_primitive_matrix,
                         integral_weights,
                         overlap_matrix,
                         position_matrix,
                         primitive_set
import ..GaussletBases.CartesianContractedParents:
    CartesianContractedParent3D,
    _CartesianExecutableProjectedQShellPayload3D,
    _CartesianProjectedQShellSidecarFixture3D,
    _cartesian_resolved_contraction_payload,
    _cartesian_resolved_contraction_payloads,
    cartesian_contracted_parent,
    contracted_parent_contraction_rules,
    contracted_parent_basis,
    contracted_parent_coefficients,
    contracted_parent_dimension,
    contracted_parent_parent_dimension,
    contracted_parent_units
import ..GaussletBases.CartesianParentGaussletBases:
    CartesianParentGaussletBasis3D,
    axis_basis,
    parent_axis_counts,
    parent_dimension

export CartesianContractedParentMetricPacket3D,
       cartesian_contracted_parent_metric_packet,
       cartesian_contracted_parent_metric_packet_dense_reference,
       contracted_parent_metric_packet_parent,
       contracted_parent_metric_packet_overlap,
       contracted_parent_metric_packet_weights,
       contracted_parent_metric_packet_centers,
       contracted_parent_metric_packet_diagnostics

"""
    CartesianContractedParentMetricPacket3D

Internal retained-space metric packet for a `CartesianContractedParent3D`.

The packet stores one-body metric diagnostics only: retained overlap,
contracted integral weights, retained position matrices/centers, and
construction diagnostics. It deliberately does not store backend state,
Gaussian supplements, residual columns, Coulomb packets, or any dense
parent-dimension-by-parent-dimension 3D matrix.
"""
struct CartesianContractedParentMetricPacket3D{P<:CartesianContractedParent3D,D}
    contracted_parent::P
    overlap::Matrix{Float64}
    weights::Vector{Float64}
    position_x::Matrix{Float64}
    position_y::Matrix{Float64}
    position_z::Matrix{Float64}
    centers::Matrix{Float64}
    first_moments::Matrix{Float64}
    diagnostics::D
end

struct _AxisMetricData1D
    overlap::Matrix{Float64}
    position::Matrix{Float64}
    weights::Vector{Float64}
    centers::Vector{Float64}
    source::Symbol
end

struct _ParentCoefficientEntry3D
    ix::Int
    iy::Int
    iz::Int
    value::Float64
end

contracted_parent_metric_packet_parent(packet::CartesianContractedParentMetricPacket3D) =
    packet.contracted_parent
contracted_parent_metric_packet_overlap(packet::CartesianContractedParentMetricPacket3D) =
    packet.overlap
contracted_parent_metric_packet_weights(packet::CartesianContractedParentMetricPacket3D) =
    packet.weights
contracted_parent_metric_packet_centers(packet::CartesianContractedParentMetricPacket3D) =
    packet.centers
contracted_parent_metric_packet_diagnostics(packet::CartesianContractedParentMetricPacket3D) =
    packet.diagnostics

function _axis_metric_data(axis)
    primitives = primitive_set(axis)
    _require_analytic_primitive_backend(
        primitives,
        "contracted Cartesian parent default metric packet axis data",
    )
    return _AxisMetricData1D(
        Matrix{Float64}(contract_primitive_matrix(axis, overlap_matrix(primitives))),
        Matrix{Float64}(contract_primitive_matrix(axis, position_matrix(primitives))),
        Float64[Float64(value) for value in integral_weights(axis)],
        Float64[Float64(value) for value in centers(axis)],
        :basis_primitive_contraction,
    )
end

function _axis_metric_data(axis, data)
    data === nothing && return _axis_metric_data(axis)
    overlap = Matrix{Float64}(getproperty(data, :overlap))
    position = Matrix{Float64}(getproperty(data, :position))
    weights = Float64[Float64(value) for value in getproperty(data, :weights)]
    center_values = Float64[Float64(value) for value in getproperty(data, :centers)]
    source = hasproperty(data, :source) ? Symbol(getproperty(data, :source)) : :explicit_axis_metrics
    return _AxisMetricData1D(overlap, position, weights, center_values, source)
end

function _parent_axis_metric_data(
    parent::CartesianParentGaussletBasis3D;
    axis_metrics = nothing,
)
    axes = (axis_basis(parent, :x), axis_basis(parent, :y), axis_basis(parent, :z))
    if axis_metrics === nothing
        return (
            x = _axis_metric_data(axes[1]),
            y = _axis_metric_data(axes[2]),
            z = _axis_metric_data(axes[3]),
        )
    end
    return (
        x = _axis_metric_data(axes[1], getproperty(axis_metrics, :x)),
        y = _axis_metric_data(axes[2], getproperty(axis_metrics, :y)),
        z = _axis_metric_data(axes[3], getproperty(axis_metrics, :z)),
    )
end

function _validate_axis_metric_data(
    metrics::NamedTuple{(:x,:y,:z)},
    dims::NTuple{3,Int},
)
    for (axis_index, axis_name) in enumerate((:x, :y, :z))
        data = getproperty(metrics, axis_name)
        n = dims[axis_index]
        size(data.overlap) == (n, n) || throw(
            DimensionMismatch("$(axis_name)-axis overlap metric size must match parent axis count"),
        )
        size(data.position) == (n, n) || throw(
            DimensionMismatch("$(axis_name)-axis position metric size must match parent axis count"),
        )
        length(data.weights) == n || throw(
            DimensionMismatch("$(axis_name)-axis integral weights must match parent axis count"),
        )
        length(data.centers) == n || throw(
            DimensionMismatch("$(axis_name)-axis centers must match parent axis count"),
        )
    end
    return metrics
end

function _unflat_parent_index(index::Int, dims::NTuple{3,Int})
    nx, ny, nz = dims
    1 <= index <= nx * ny * nz || throw(
        ArgumentError("coefficient row index lies outside the Cartesian parent dimension"),
    )
    shifted = index - 1
    plane = ny * nz
    ix = shifted ÷ plane + 1
    remainder = shifted % plane
    iy = remainder ÷ nz + 1
    iz = remainder % nz + 1
    return ix, iy, iz
end

function _coefficient_column_entries(
    coefficients::SparseArrays.SparseMatrixCSC{Float64,Int},
    dims::NTuple{3,Int},
)
    entries = [Vector{_ParentCoefficientEntry3D}() for _ in axes(coefficients, 2)]
    rowvals = SparseArrays.rowvals(coefficients)
    nzvals = SparseArrays.nonzeros(coefficients)
    colptr = SparseArrays.getcolptr(coefficients)
    for col in axes(coefficients, 2)
        column_entries = entries[col]
        for ptr in colptr[col]:(colptr[col + 1] - 1)
            value = nzvals[ptr]
            iszero(value) && continue
            ix, iy, iz = _unflat_parent_index(rowvals[ptr], dims)
            push!(column_entries, _ParentCoefficientEntry3D(ix, iy, iz, value))
        end
    end
    return entries
end

function _coefficient_column_entries(
    coefficients::AbstractMatrix{<:Real},
    dims::NTuple{3,Int},
)
    entries = [Vector{_ParentCoefficientEntry3D}() for _ in axes(coefficients, 2)]
    for col in axes(coefficients, 2)
        column_entries = entries[col]
        for row in axes(coefficients, 1)
            value = Float64(coefficients[row, col])
            iszero(value) && continue
            ix, iy, iz = _unflat_parent_index(Int(row), dims)
            push!(column_entries, _ParentCoefficientEntry3D(ix, iy, iz, value))
        end
    end
    return entries
end

function _coefficient_nnz(coefficients::SparseArrays.SparseMatrixCSC{Float64,Int})
    return SparseArrays.nnz(coefficients)
end

function _coefficient_nnz(coefficients::AbstractMatrix{<:Real})
    count = 0
    for value in coefficients
        iszero(value) || (count += 1)
    end
    return count
end

function _metric_value(
    left::_ParentCoefficientEntry3D,
    right::_ParentCoefficientEntry3D,
    mx::AbstractMatrix{<:Real},
    my::AbstractMatrix{<:Real},
    mz::AbstractMatrix{<:Real},
)
    return left.value * right.value *
           Float64(mx[left.ix, right.ix]) *
           Float64(my[left.iy, right.iy]) *
           Float64(mz[left.iz, right.iz])
end

function _contract_pair_matrix(
    entries::Vector{Vector{_ParentCoefficientEntry3D}},
    mx::AbstractMatrix{<:Real},
    my::AbstractMatrix{<:Real},
    mz::AbstractMatrix{<:Real},
)
    ncols = length(entries)
    result = zeros(Float64, ncols, ncols)
    for col_b in 1:ncols
        right_entries = entries[col_b]
        for col_a in 1:col_b
            left_entries = entries[col_a]
            value = 0.0
            for left in left_entries, right in right_entries
                value += _metric_value(left, right, mx, my, mz)
            end
            result[col_a, col_b] = value
            result[col_b, col_a] = value
        end
    end
    return result
end

function _contract_linear_vector(
    entries::Vector{Vector{_ParentCoefficientEntry3D}},
    vx::AbstractVector{<:Real},
    vy::AbstractVector{<:Real},
    vz::AbstractVector{<:Real},
)
    result = zeros(Float64, length(entries))
    for col in eachindex(entries)
        value = 0.0
        for entry in entries[col]
            value += entry.value *
                     Float64(vx[entry.ix]) *
                     Float64(vy[entry.iy]) *
                     Float64(vz[entry.iz])
        end
        result[col] = value
    end
    return result
end

function _contracted_first_moments(
    entries::Vector{Vector{_ParentCoefficientEntry3D}},
    metrics::NamedTuple{(:x,:y,:z)},
)
    wx, wy, wz = metrics.x.weights, metrics.y.weights, metrics.z.weights
    xw = metrics.x.centers .* metrics.x.weights
    yw = metrics.y.centers .* metrics.y.weights
    zw = metrics.z.centers .* metrics.z.weights
    return hcat(
        _contract_linear_vector(entries, xw, wy, wz),
        _contract_linear_vector(entries, wx, yw, wz),
        _contract_linear_vector(entries, wx, wy, zw),
    )
end

function _staged_unit(contraction_unit)
    metadata = contraction_unit.metadata
    hasproperty(metadata, :staged_by_center_unit) &&
        return getproperty(metadata, :staged_by_center_unit)
    return nothing
end

function _metric_dispatch_unit_path_from_payload(staged_unit)
    if staged_unit === nothing
        return (
            family = :unsupported,
            kind = :missing_staged_payload,
            metric_capability = :none,
            linear_vector_path = :unsupported,
            block_role = :unsupported,
            unsupported = true,
            prototype = false,
        )
    elseif staged_unit.kind == :product_doside
        return (
            family = :product_owned_unit,
            kind = staged_unit.kind,
            metric_capability = :product_staged_metric_contraction,
            linear_vector_path = :product_staged_axis_projection,
            block_role = :product,
            unsupported = false,
            prototype = false,
        )
    elseif staged_unit.kind == :support_dense
        return (
            family = :support_dense_fallback,
            kind = staged_unit.kind,
            metric_capability = :support_local_product,
            linear_vector_path = :support_local_fallback,
            block_role = :fallback,
            unsupported = false,
            prototype = false,
        )
    end
    return (
        family = :unsupported,
        kind = staged_unit.kind,
        metric_capability = :none,
        linear_vector_path = :unsupported,
        block_role = :unsupported,
        unsupported = true,
        prototype = false,
    )
end

function _metric_dispatch_unit_path_from_rule(rule)
    if rule.rule_family == :product_owned_unit &&
       rule.metric_capability == :product_staged_metric_contraction
        return (
            family = rule.rule_family,
            kind = rule.kind,
            metric_capability = rule.metric_capability,
            linear_vector_path = :product_staged_axis_projection,
            block_role = :product,
            unsupported = false,
            prototype = false,
        )
    elseif rule.rule_family == :support_dense_fallback &&
           rule.metric_capability == :support_local_product
        return (
            family = rule.rule_family,
            kind = rule.kind,
            metric_capability = rule.metric_capability,
            linear_vector_path = :support_local_fallback,
            block_role = :fallback,
            unsupported = false,
            prototype = false,
        )
    end
    prototype =
        rule.rule_family == :projected_q_shell_boundary_modes ||
        rule.metric_capability == :pqs_low_order_product_metric_prototype
    return (
        family = rule.rule_family,
        kind = rule.kind,
        metric_capability = rule.metric_capability,
        linear_vector_path = :unsupported,
        block_role = prototype ? :unsupported_prototype : :unsupported,
        unsupported = true,
        prototype = prototype,
    )
end

function _metric_dispatch_unit_path_from_resolved_payload(payload)
    if payload.ready_for_metric_execution &&
       payload.payload_kind == :product_doside &&
       payload.metric_path == :product_staged_metric_contraction
        return (
            family = :product_owned_unit,
            kind = payload.payload_kind,
            metric_capability = payload.diagnostics.metric_capability,
            linear_vector_path = payload.diagnostics.linear_vector_path,
            block_role = payload.diagnostics.block_role,
            unsupported = false,
            prototype = false,
        )
    elseif payload.ready_for_metric_execution &&
           payload.payload_kind == :support_dense &&
           payload.metric_path == :support_local_product
        return (
            family = :support_dense_fallback,
            kind = payload.payload_kind,
            metric_capability = payload.diagnostics.metric_capability,
            linear_vector_path = payload.diagnostics.linear_vector_path,
            block_role = payload.diagnostics.block_role,
            unsupported = false,
            prototype = false,
        )
    elseif payload.ready_for_metric_execution &&
           payload.payload_kind == :projected_q_shell &&
           payload.metric_path == :pqs_low_order_support_local_reference
        return (
            family = :projected_q_shell_boundary_modes,
            kind = payload.payload_kind,
            metric_capability = payload.diagnostics.metric_capability,
            linear_vector_path = payload.diagnostics.linear_vector_path,
            block_role = payload.diagnostics.block_role,
            unsupported = false,
            prototype = false,
        )
    end
    prototype =
        payload.metric_path == :unsupported_prototype ||
        getproperty(payload.diagnostics, :prototype)
    return (
        family = getproperty(payload.diagnostics, :rule_family),
        kind = payload.payload_kind,
        metric_capability = getproperty(payload.diagnostics, :metric_capability),
        linear_vector_path = :unsupported,
        block_role = prototype ? :unsupported_prototype : :unsupported,
        unsupported = true,
        prototype = prototype,
    )
end

function _metric_dispatch_block_path(left_path, right_path)
    (left_path.unsupported || right_path.unsupported) && return :unsupported
    left_path.block_role == :product && right_path.block_role == :product &&
        return :product_product
    left_path.block_role == :pqs && right_path.block_role == :pqs &&
        return :pqs_pqs_low_order_reference
    if (left_path.block_role == :pqs && right_path.block_role == :product) ||
       (left_path.block_role == :product && right_path.block_role == :pqs)
        return :unsupported_pqs_product_optimized
    end
    if (left_path.block_role == :pqs && right_path.block_role == :fallback) ||
       (left_path.block_role == :fallback && right_path.block_role == :pqs)
        return :pqs_support_local_reference
    end
    return :support_local_fallback
end

function _pqs_product_mixed_block_policy(; explicit_reference_requested::Bool = false)
    return (
        pair_kind = :pqs_product_mixed,
        optimized_metric_path = :unsupported_pqs_product_optimized,
        optimized_supported = false,
        support_local_reference_allowed = explicit_reference_requested,
        support_local_reference_path = explicit_reference_requested ?
                                       :support_local_reference_explicit_only :
                                       :not_available_without_explicit_request,
        fixture_only = true,
        production_supported = false,
        reason = :pqs_product_optimized_metric_not_implemented,
        required_next_step = :separate_reference_or_optimized_helper_design,
    )
end

function _metric_dispatch_plan_from_unit_paths(unit_paths; source::Symbol)
    product_unit_count = count(path -> path.block_role == :product, unit_paths)
    fallback_unit_count = count(path -> path.block_role == :fallback, unit_paths)
    pqs_unit_count = count(path -> path.block_role == :pqs, unit_paths)
    unsupported_unit_count = count(path -> path.unsupported, unit_paths)
    prototype_rule_count = count(path -> path.prototype, unit_paths)
    block_paths = NamedTuple[]
    product_block_count = 0
    fallback_block_count = 0
    pqs_pqs_block_count = 0
    pqs_support_block_count = 0
    pqs_product_unsupported_block_count = 0
    unsupported_block_count = 0
    for right_index in eachindex(unit_paths)
        for left_index in 1:right_index
            path = _metric_dispatch_block_path(unit_paths[left_index], unit_paths[right_index])
            path == :product_product && (product_block_count += 1)
            path == :support_local_fallback && (fallback_block_count += 1)
            path == :pqs_pqs_low_order_reference && (pqs_pqs_block_count += 1)
            path == :pqs_support_local_reference && begin
                pqs_support_block_count += 1
                fallback_block_count += 1
            end
            path == :unsupported_pqs_product_optimized && begin
                pqs_product_unsupported_block_count += 1
                unsupported_block_count += 1
            end
            path == :unsupported && (unsupported_block_count += 1)
            push!(
                block_paths,
                (left_index = left_index, right_index = right_index, path = path),
            )
        end
    end
    return (
        source = source,
        unit_count = length(unit_paths),
        unit_paths = Tuple(unit_paths),
        block_paths = Tuple(block_paths),
        product_unit_count = product_unit_count,
        support_fallback_unit_count = fallback_unit_count,
        pqs_unit_count = pqs_unit_count,
        unsupported_unit_count = unsupported_unit_count,
        prototype_rule_count = prototype_rule_count,
        product_product_block_count = product_block_count,
        fallback_block_count = fallback_block_count,
        pqs_pqs_block_count = pqs_pqs_block_count,
        pqs_support_block_count = pqs_support_block_count,
        pqs_product_unsupported_block_count = pqs_product_unsupported_block_count,
        unsupported_block_count = unsupported_block_count,
        plan_supported = unsupported_unit_count == 0 && unsupported_block_count == 0,
    )
end

function _metric_dispatch_plan_from_resolved_payloads(
    payloads::AbstractVector;
    source::Symbol = :resolved_payload_fixture,
)
    return _metric_dispatch_plan_from_unit_paths(
        [_metric_dispatch_unit_path_from_resolved_payload(payload) for payload in payloads];
        source,
    )
end

function _contracted_parent_resolved_payloads(
    contracted_parent::CartesianContractedParent3D,
)
    parent_dim = contracted_parent_parent_dimension(contracted_parent)
    return [
        _cartesian_resolved_contraction_payload(unit; parent_dimension = parent_dim)
        for unit in contracted_parent_units(contracted_parent)
    ]
end

function _contracted_parent_metric_dispatch_plan_from_payload(
    contracted_parent::CartesianContractedParent3D,
)
    unit_paths = [
        _metric_dispatch_unit_path_from_payload(_staged_unit(unit)) for
        unit in contracted_parent_units(contracted_parent)
    ]
    return _metric_dispatch_plan_from_unit_paths(
        unit_paths;
        source = :staged_payload,
    )
end

function _contracted_parent_metric_dispatch_plan_from_resolved_payloads(
    contracted_parent::CartesianContractedParent3D,
)
    unit_paths = [
        _metric_dispatch_unit_path_from_resolved_payload(payload)
        for payload in _contracted_parent_resolved_payloads(contracted_parent)
    ]
    return _metric_dispatch_plan_from_unit_paths(
        unit_paths;
        source = :resolved_payload,
    )
end

function _contracted_parent_metric_dispatch_plan_from_rules(rules::AbstractVector)
    unit_paths = [_metric_dispatch_unit_path_from_rule(rule) for rule in rules]
    return _metric_dispatch_plan_from_unit_paths(
        unit_paths;
        source = :contraction_rules,
    )
end

function _contracted_parent_metric_dispatch_plan_from_rules(
    contracted_parent::CartesianContractedParent3D,
)
    return _contracted_parent_metric_dispatch_plan_from_rules(
        contracted_parent_contraction_rules(contracted_parent),
    )
end

function _metric_dispatch_plan_agreement(payload_plan, rule_plan)
    mismatch_fields = Symbol[]
    for field in (
        :unit_count,
        :product_unit_count,
        :support_fallback_unit_count,
        :unsupported_unit_count,
        :prototype_rule_count,
        :product_product_block_count,
        :fallback_block_count,
        :unsupported_block_count,
        :plan_supported,
    )
        getproperty(payload_plan, field) == getproperty(rule_plan, field) ||
            push!(mismatch_fields, field)
    end
    payload_units = [
        (path.family, path.kind, path.metric_capability, path.linear_vector_path, path.block_role) for
        path in payload_plan.unit_paths
    ]
    rule_units = [
        (path.family, path.kind, path.metric_capability, path.linear_vector_path, path.block_role) for
        path in rule_plan.unit_paths
    ]
    payload_units == rule_units || push!(mismatch_fields, :unit_paths)
    payload_blocks = [path.path for path in payload_plan.block_paths]
    rule_blocks = [path.path for path in rule_plan.block_paths]
    payload_blocks == rule_blocks || push!(mismatch_fields, :block_paths)
    return (
        agree = isempty(mismatch_fields),
        mismatch_fields = Tuple(mismatch_fields),
    )
end

function _contracted_parent_metric_dispatch_shadow_plan(
    contracted_parent::CartesianContractedParent3D,
)
    payload_plan = _contracted_parent_metric_dispatch_plan_from_payload(contracted_parent)
    rule_plan = _contracted_parent_metric_dispatch_plan_from_rules(contracted_parent)
    resolved_plan =
        _contracted_parent_metric_dispatch_plan_from_resolved_payloads(contracted_parent)
    return (
        payload_plan = payload_plan,
        rule_plan = rule_plan,
        resolved_plan = resolved_plan,
        comparison = _metric_dispatch_plan_agreement(payload_plan, rule_plan),
        resolved_comparison = _metric_dispatch_plan_agreement(
            payload_plan,
            resolved_plan,
        ),
    )
end

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
    projected_x_overlap = _project_staged_axis_matrix(
        left_unit.axes[1],
        right_unit.axes[1],
        metrics.x.overlap,
    )
    projected_y_overlap = _project_staged_axis_matrix(
        left_unit.axes[2],
        right_unit.axes[2],
        metrics.y.overlap,
    )
    projected_z_overlap = _project_staged_axis_matrix(
        left_unit.axes[3],
        right_unit.axes[3],
        metrics.z.overlap,
    )
    projected_x_position = _project_staged_axis_matrix(
        left_unit.axes[1],
        right_unit.axes[1],
        metrics.x.position,
    )
    projected_y_position = _project_staged_axis_matrix(
        left_unit.axes[2],
        right_unit.axes[2],
        metrics.y.position,
    )
    projected_z_position = _project_staged_axis_matrix(
        left_unit.axes[3],
        right_unit.axes[3],
        metrics.z.position,
    )
    same_unit = left_unit.column_range == right_unit.column_range
    @inbounds for local_col in eachindex(right_unit.axis_function_indices)
        xj, yj, zj = right_unit.axis_function_indices[local_col]
        global_col = right_unit.column_range[local_col]
        for local_row in eachindex(left_unit.axis_function_indices)
            xi, yi, zi = left_unit.axis_function_indices[local_row]
            global_row = left_unit.column_range[local_row]
            y_overlap = projected_y_overlap[yi, yj]
            z_overlap = projected_z_overlap[zi, zj]
            overlap_value = projected_x_overlap[xi, xj] * y_overlap * z_overlap
            position_x_value = projected_x_position[xi, xj] * y_overlap * z_overlap
            position_y_value =
                projected_x_overlap[xi, xj] * projected_y_position[yi, yj] * z_overlap
            position_z_value =
                projected_x_overlap[xi, xj] * y_overlap * projected_z_position[zi, zj]
            overlap[global_row, global_col] = overlap_value
            position_x[global_row, global_col] = position_x_value
            position_y[global_row, global_col] = position_y_value
            position_z[global_row, global_col] = position_z_value
            if !same_unit
                overlap[global_col, global_row] = overlap_value
                position_x[global_col, global_row] = position_x_value
                position_y[global_col, global_row] = position_y_value
                position_z[global_col, global_row] = position_z_value
            end
        end
    end
    return nothing
end

function _staged_unit_entries(unit)
    entries = [Vector{_ParentCoefficientEntry3D}() for _ in unit.column_range]
    coefficients = unit.coefficient_matrix
    for local_col in axes(coefficients, 2)
        column_entries = entries[local_col]
        for local_row in axes(coefficients, 1)
            value = Float64(coefficients[local_row, local_col])
            iszero(value) && continue
            ix, iy, iz = unit.support_states[local_row]
            push!(column_entries, _ParentCoefficientEntry3D(ix, iy, iz, value))
        end
    end
    return entries
end

function _staged_unit_entries(unit::_CartesianExecutableProjectedQShellPayload3D)
    entries = [Vector{_ParentCoefficientEntry3D}() for _ in unit.column_range]
    coefficients = unit.support_coefficient_matrix
    for local_col in axes(coefficients, 2)
        column_entries = entries[local_col]
        for local_row in axes(coefficients, 1)
            value = Float64(coefficients[local_row, local_col])
            iszero(value) && continue
            ix, iy, iz = unit.support_states[local_row]
            push!(column_entries, _ParentCoefficientEntry3D(ix, iy, iz, value))
        end
    end
    return entries
end

function _entries_from_resolved_payload(payload)
    payload.ready_for_metric_execution || throw(
        ArgumentError("resolved payload $(payload.payload_kind) is not ready for metric execution"),
    )
    payload.payload_kind in (:projected_q_shell, :support_dense, :product_doside) || throw(
        ArgumentError("unsupported resolved low-order metric payload kind $(payload.payload_kind)"),
    )
    return _staged_unit_entries(payload.payload)
end

function _resolved_payload_low_order_metric_block(left, right, metrics::NamedTuple{(:x,:y,:z)})
    left_path = _metric_dispatch_unit_path_from_resolved_payload(left)
    right_path = _metric_dispatch_unit_path_from_resolved_payload(right)
    block_path = _metric_dispatch_block_path(left_path, right_path)
    block_path == :unsupported && throw(
        ArgumentError("resolved payload low-order metric block is unsupported"),
    )
    block_path == :unsupported_pqs_product_optimized && throw(
        ArgumentError(
            "PQS/product optimized metric blocks are explicitly unsupported; " *
            "support-local reference requires a separately named explicit reference helper",
        ),
    )

    left_entries = _entries_from_resolved_payload(left)
    right_entries = _entries_from_resolved_payload(right)
    support_overlap = _contract_pair_block(
        left_entries,
        right_entries,
        metrics.x.overlap,
        metrics.y.overlap,
        metrics.z.overlap,
    )
    same_pqs_payload =
        left.payload_kind == :projected_q_shell &&
        right.payload_kind == :projected_q_shell &&
        left.payload === right.payload
    overlap =
        same_pqs_payload ?
        Matrix{Float64}(LinearAlgebra.I, length(left_entries), length(right_entries)) :
        support_overlap
    position_x = _contract_pair_block(
        left_entries,
        right_entries,
        metrics.x.position,
        metrics.y.overlap,
        metrics.z.overlap,
    )
    position_y = _contract_pair_block(
        left_entries,
        right_entries,
        metrics.x.overlap,
        metrics.y.position,
        metrics.z.overlap,
    )
    position_z = _contract_pair_block(
        left_entries,
        right_entries,
        metrics.x.overlap,
        metrics.y.overlap,
        metrics.z.position,
    )
    weights = same_pqs_payload ?
              _contract_linear_vector(
                  left_entries,
                  metrics.x.weights,
                  metrics.y.weights,
                  metrics.z.weights,
              ) : nothing
    first_moments = same_pqs_payload ?
                    _contracted_first_moments(left_entries, metrics) : nothing
    return (
        path = block_path,
        overlap = overlap,
        position_x = position_x,
        position_y = position_y,
        position_z = position_z,
        weights = weights,
        first_moments = first_moments,
        diagnostics = (
            source = :resolved_payload_low_order_metric_block,
            block_path = block_path,
            left_kind = left.payload_kind,
            right_kind = right.payload_kind,
            support_local_reference_used =
                block_path in (:pqs_support_local_reference, :support_local_fallback),
            pqs_self_overlap_invariant_applied = same_pqs_payload,
            pqs_cross_overlap_identity_shortcut_used = false,
            support_overlap_debug_error =
                same_pqs_payload ?
                LinearAlgebra.norm(support_overlap - LinearAlgebra.I, Inf) : nothing,
            production_optimized_pqs_product_path = false,
            supported_terms = (:overlap, :weights, :position_x, :position_y, :position_z),
            unsupported_terms = (
                :kinetic,
                :x2,
                :nuclear_one_body,
                :gaussian_sum,
                :pair_sum,
                :interaction,
            ),
        ),
    )
end

function _projected_q_shell_sidecar_low_order_metric_reference(
    sidecar::_CartesianProjectedQShellSidecarFixture3D,
    metrics::NamedTuple{(:x,:y,:z)};
    mixed_payloads = (),
)
    resolved_payloads = _cartesian_resolved_contraction_payloads(sidecar)
    isempty(resolved_payloads) && throw(
        ArgumentError("PQS sidecar metric reference requires at least one payload"),
    )
    all(payload -> payload.payload_kind == :projected_q_shell, resolved_payloads) || throw(
        ArgumentError("PQS sidecar metric reference only accepts projected_q_shell payloads"),
    )
    self_blocks = [
        _resolved_payload_low_order_metric_block(payload, payload, metrics)
        for payload in resolved_payloads
    ]
    mixed_blocks = [
        _resolved_payload_low_order_metric_block(payload, mixed, metrics)
        for payload in resolved_payloads
        for mixed in mixed_payloads
    ]
    return (
        resolved_payloads = resolved_payloads,
        self_blocks = self_blocks,
        mixed_blocks = mixed_blocks,
        diagnostics = (
            source = :projected_q_shell_sidecar_low_order_metric_reference,
            fixture_only = true,
            reference_scoped = true,
            production_supported = false,
            support_local_reference_only = true,
            fixed_block_sidecar_installed = sidecar.diagnostics.fixed_block_sidecar_installed,
            default_builder_consumes = sidecar.diagnostics.default_builder_consumes,
            qw_consumes = sidecar.diagnostics.qw_consumes,
            hamiltonian_consumes = sidecar.diagnostics.hamiltonian_consumes,
            payload_count = length(resolved_payloads),
            self_block_count = length(self_blocks),
            mixed_block_count = length(mixed_blocks),
            metric_capability = :pqs_low_order_support_local_reference,
            supported_terms = (:overlap, :weights, :position_x, :position_y, :position_z),
            unsupported_terms = (
                :kinetic,
                :x2,
                :nuclear_one_body,
                :gaussian_sum,
                :pair_sum,
                :interaction,
            ),
            pqs_product_optimized_path_ready = false,
        ),
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

function _staged_unit_linear_vectors(unit, metrics::NamedTuple{(:x,:y,:z)})
    if unit.kind == :product_doside
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

function _packet_diagnostics(;
    construction_path::Symbol,
    dense_parent_matrix_used::Bool,
    dense_reference_parent_dimension_limit = nothing,
    contracted_parent::CartesianContractedParent3D,
    entries::Vector{Vector{_ParentCoefficientEntry3D}},
    metrics::NamedTuple{(:x,:y,:z)},
    overlap::AbstractMatrix{<:Real},
)
    column_nnz = Int[length(column) for column in entries]
    return (
        construction_path = construction_path,
        dense_parent_matrix_used = dense_parent_matrix_used,
        dense_reference_parent_dimension_limit = dense_reference_parent_dimension_limit,
        parent_dimension = contracted_parent_parent_dimension(contracted_parent),
        contracted_dimension = contracted_parent_dimension(contracted_parent),
        coefficient_storage = nameof(typeof(contracted_parent_coefficients(contracted_parent))),
        coefficient_nnz = sum(column_nnz),
        max_column_nnz = isempty(column_nnz) ? 0 : maximum(column_nnz),
        axis_metric_sources = (
            x = metrics.x.source,
            y = metrics.y.source,
            z = metrics.z.source,
        ),
        overlap_symmetry_error = LinearAlgebra.norm(overlap - transpose(overlap), Inf),
        overlap_identity_error = LinearAlgebra.norm(overlap - LinearAlgebra.I, Inf),
    )
end

function _packet_diagnostics_from_counts(;
    construction_path::Symbol,
    dense_parent_matrix_used::Bool,
    contracted_parent::CartesianContractedParent3D,
    metrics::NamedTuple{(:x,:y,:z)},
    overlap::AbstractMatrix{<:Real},
    coefficient_nnz::Int,
    max_column_nnz::Int,
    product_unit_count::Int = 0,
    generic_unit_count::Int = 0,
    product_block_count::Int = 0,
    fallback_block_count::Int = 0,
    max_unit_support_count::Int = 0,
    max_unit_retained_count::Int = 0,
)
    return (
        construction_path = construction_path,
        dense_parent_matrix_used = dense_parent_matrix_used,
        dense_reference_parent_dimension_limit = nothing,
        parent_dimension = contracted_parent_parent_dimension(contracted_parent),
        contracted_dimension = contracted_parent_dimension(contracted_parent),
        coefficient_storage = nameof(typeof(contracted_parent_coefficients(contracted_parent))),
        coefficient_nnz = coefficient_nnz,
        max_column_nnz = max_column_nnz,
        product_unit_count = product_unit_count,
        generic_unit_count = generic_unit_count,
        product_block_count = product_block_count,
        fallback_block_count = fallback_block_count,
        max_unit_support_count = max_unit_support_count,
        max_unit_retained_count = max_unit_retained_count,
        axis_metric_sources = (
            x = metrics.x.source,
            y = metrics.y.source,
            z = metrics.z.source,
        ),
        overlap_symmetry_error = LinearAlgebra.norm(overlap - transpose(overlap), Inf),
        overlap_identity_error = LinearAlgebra.norm(overlap - LinearAlgebra.I, Inf),
    )
end

function _metric_packet_from_matrices(
    contracted_parent::CartesianContractedParent3D,
    metrics::NamedTuple{(:x,:y,:z)},
    entries::Vector{Vector{_ParentCoefficientEntry3D}},
    overlap::Matrix{Float64},
    position_x::Matrix{Float64},
    position_y::Matrix{Float64},
    position_z::Matrix{Float64},
    weights::Vector{Float64},
    first_moments::Matrix{Float64};
    diagnostics,
)
    centers = hcat(
        LinearAlgebra.diag(position_x),
        LinearAlgebra.diag(position_y),
        LinearAlgebra.diag(position_z),
    )
    return CartesianContractedParentMetricPacket3D(
        contracted_parent,
        Matrix{Float64}(LinearAlgebra.Symmetric(overlap)),
        weights,
        Matrix{Float64}(LinearAlgebra.Symmetric(position_x)),
        Matrix{Float64}(LinearAlgebra.Symmetric(position_y)),
        Matrix{Float64}(LinearAlgebra.Symmetric(position_z)),
        Matrix{Float64}(centers),
        first_moments,
        diagnostics,
    )
end

function _support_local_metric_packet(
    contracted_parent::CartesianContractedParent3D,
    metrics::NamedTuple{(:x,:y,:z)},
    dims::NTuple{3,Int},
)
    coefficients = contracted_parent_coefficients(contracted_parent)
    entries = _coefficient_column_entries(coefficients, dims)

    overlap = _contract_pair_matrix(
        entries,
        metrics.x.overlap,
        metrics.y.overlap,
        metrics.z.overlap,
    )
    position_x = _contract_pair_matrix(
        entries,
        metrics.x.position,
        metrics.y.overlap,
        metrics.z.overlap,
    )
    position_y = _contract_pair_matrix(
        entries,
        metrics.x.overlap,
        metrics.y.position,
        metrics.z.overlap,
    )
    position_z = _contract_pair_matrix(
        entries,
        metrics.x.overlap,
        metrics.y.overlap,
        metrics.z.position,
    )
    weights = _contract_linear_vector(
        entries,
        metrics.x.weights,
        metrics.y.weights,
        metrics.z.weights,
    )
    first_moments = _contracted_first_moments(entries, metrics)
    diagnostics = _packet_diagnostics(
        construction_path = :support_local_product,
        dense_parent_matrix_used = false,
        contracted_parent = contracted_parent,
        entries = entries,
        metrics = metrics,
        overlap = overlap,
    )
    return _metric_packet_from_matrices(
        contracted_parent,
        metrics,
        entries,
        overlap,
        position_x,
        position_y,
        position_z,
        weights,
        first_moments;
        diagnostics,
    )
end

function _product_unit_metric_packet(
    contracted_parent::CartesianContractedParent3D,
    metrics::NamedTuple{(:x,:y,:z)},
    dims::NTuple{3,Int},
)
    return _resolved_payload_product_staged_metric_packet(
        contracted_parent,
        metrics;
        construction_path = :product_staged_metric_contraction,
        extra_diagnostics = (;),
        include_resolved_payload_count = false,
    )
end

function _resolved_payload_product_staged_metric_packet(
    contracted_parent::CartesianContractedParent3D,
    metrics::NamedTuple{(:x,:y,:z)},
    ;
    construction_path::Symbol = :resolved_payload_product_staged_metric_contraction,
    extra_diagnostics = (
        source = :resolved_payload_metric_shadow,
        default_metric_execution_changed = false,
    ),
    include_resolved_payload_count::Bool = true,
)
    resolved_payloads = _contracted_parent_resolved_payloads(contracted_parent)
    isempty(resolved_payloads) && throw(
        ArgumentError("resolved-payload metric contraction requires contracted-parent units"),
    )
    staged_units = Any[]
    product_unit_count = 0
    generic_unit_count = 0
    for resolved in resolved_payloads
        resolved.ready_for_metric_execution || throw(
            ArgumentError(
                "resolved-payload metric contraction cannot consume $(resolved.payload_kind): missing $(resolved.missing_fields)",
            ),
        )
        resolved.payload_kind == :product_doside && (product_unit_count += 1)
        resolved.payload_kind == :support_dense && (generic_unit_count += 1)
        resolved.payload_kind in (:product_doside, :support_dense) || throw(
            ArgumentError("unsupported resolved metric payload kind $(resolved.payload_kind)"),
        )
        push!(staged_units, resolved.payload)
    end
    product_unit_count == 0 && throw(
        ArgumentError("resolved-payload metric contraction found no product_doside units"),
    )
    contracted_dimension = contracted_parent_dimension(contracted_parent)
    overlap = zeros(Float64, contracted_dimension, contracted_dimension)
    position_x = zeros(Float64, contracted_dimension, contracted_dimension)
    position_y = zeros(Float64, contracted_dimension, contracted_dimension)
    position_z = zeros(Float64, contracted_dimension, contracted_dimension)
    product_block_count = 0
    fallback_block_count = 0
    entries_cache = Vector{Union{Nothing,Vector{Vector{_ParentCoefficientEntry3D}}}}(
        undef,
        length(staged_units),
    )
    fill!(entries_cache, nothing)
    for right_index in eachindex(staged_units)
        right = staged_units[right_index]
        right_range = right.column_range
        for left_index in 1:right_index
            left = staged_units[left_index]
            left_range = left.column_range
            if left.kind == :product_doside && right.kind == :product_doside
                _fill_product_staged_metric_blocks!(
                    overlap,
                    position_x,
                    position_y,
                    position_z,
                    left,
                    right,
                    metrics,
                )
                product_block_count += 1
            else
                left_entries = _cached_staged_unit_entries!(entries_cache, left_index, left)
                right_entries = _cached_staged_unit_entries!(entries_cache, right_index, right)
                blocks = _fallback_staged_metric_blocks(left_entries, right_entries, metrics)
                fallback_block_count += 1
                overlap[left_range, right_range] .= blocks.overlap
                position_x[left_range, right_range] .= blocks.position_x
                position_y[left_range, right_range] .= blocks.position_y
                position_z[left_range, right_range] .= blocks.position_z
                if left_index != right_index
                    overlap[right_range, left_range] .= transpose(blocks.overlap)
                    position_x[right_range, left_range] .= transpose(blocks.position_x)
                    position_y[right_range, left_range] .= transpose(blocks.position_y)
                    position_z[right_range, left_range] .= transpose(blocks.position_z)
                end
            end
        end
    end
    weights = zeros(Float64, contracted_dimension)
    first_moments = zeros(Float64, contracted_dimension, 3)
    for unit in staged_units
        unit_weights, unit_first_moments = _staged_unit_linear_vectors(unit, metrics)
        weights[unit.column_range] .= unit_weights
        first_moments[unit.column_range, :] .= unit_first_moments
    end
    retained_counts = Int[length(unit.column_range) for unit in staged_units]
    support_counts = Int[length(unit.support_indices) for unit in staged_units]
    diagnostics = _packet_diagnostics_from_counts(
        construction_path = construction_path,
        dense_parent_matrix_used = false,
        contracted_parent = contracted_parent,
        metrics = metrics,
        overlap = overlap,
        coefficient_nnz = _coefficient_nnz(contracted_parent_coefficients(contracted_parent)),
        max_column_nnz = 0,
        product_unit_count = product_unit_count,
        generic_unit_count = generic_unit_count,
        product_block_count = product_block_count,
        fallback_block_count = fallback_block_count,
        max_unit_support_count = isempty(support_counts) ? 0 : maximum(support_counts),
        max_unit_retained_count = isempty(retained_counts) ? 0 : maximum(retained_counts),
    )
    return _metric_packet_from_matrices(
        contracted_parent,
        metrics,
        Vector{_ParentCoefficientEntry3D}[],
        overlap,
        position_x,
        position_y,
        position_z,
        weights,
        first_moments;
        diagnostics = merge(
            diagnostics,
            include_resolved_payload_count ?
            (resolved_payload_count = length(resolved_payloads),) : (;),
            extra_diagnostics,
        ),
    )
end

function _resolved_payload_product_staged_metric_packet(
    contracted_parent::CartesianContractedParent3D;
    axis_metrics = nothing,
)
    parent = contracted_parent_basis(contracted_parent)
    dims = parent_axis_counts(parent)
    metrics = _validate_axis_metric_data(
        _parent_axis_metric_data(parent; axis_metrics),
        dims,
    )
    return _resolved_payload_product_staged_metric_packet(contracted_parent, metrics)
end

"""
    cartesian_contracted_parent_metric_packet(contracted_parent; axis_metrics = nothing, construction_path = :support_local_product)

Build a retained metric packet from a contracted Cartesian parent without
forming dense full-parent 3D matrices. The default `axis_metrics` are derived
from the parent axes; callers with staged/product route data may pass explicit
axis metric named tuples with `overlap`, `position`, `weights`, and `centers`.

`construction_path = :support_local_product` is the generic sparse/local
fallback. `:product_staged_metric_contraction` consumes existing nested
product-staged sidecar units, using separable product contraction for
product/product blocks and support-local contraction for generic or mixed
blocks.
"""
function cartesian_contracted_parent_metric_packet(
    contracted_parent::CartesianContractedParent3D;
    axis_metrics = nothing,
    construction_path::Symbol = :support_local_product,
)
    parent = contracted_parent_basis(contracted_parent)
    dims = parent_axis_counts(parent)
    metrics = _validate_axis_metric_data(
        _parent_axis_metric_data(parent; axis_metrics),
        dims,
    )
    construction_path == :support_local_product &&
        return _support_local_metric_packet(contracted_parent, metrics, dims)
    construction_path == :product_staged_metric_contraction &&
        return _product_unit_metric_packet(contracted_parent, metrics, dims)
    throw(
        ArgumentError(
            "unsupported contracted-parent metric packet construction_path $(repr(construction_path)); " *
            "supported values are :support_local_product and :product_staged_metric_contraction",
        ),
    )
end

function cartesian_contracted_parent_metric_packet(
    fixed_block::_NestedFixedBlock3D;
    kwargs...,
)
    return cartesian_contracted_parent_metric_packet(
        cartesian_contracted_parent(fixed_block);
        kwargs...,
    )
end

function _dense_parent_product_matrix(
    mx::AbstractMatrix{<:Real},
    my::AbstractMatrix{<:Real},
    mz::AbstractMatrix{<:Real},
)
    return LinearAlgebra.kron(
        LinearAlgebra.kron(Matrix{Float64}(mx), Matrix{Float64}(my)),
        Matrix{Float64}(mz),
    )
end

function _dense_parent_product_vector(
    vx::AbstractVector{<:Real},
    vy::AbstractVector{<:Real},
    vz::AbstractVector{<:Real},
)
    result = Vector{Float64}(undef, length(vx) * length(vy) * length(vz))
    cursor = 1
    for ix in eachindex(vx), iy in eachindex(vy), iz in eachindex(vz)
        result[cursor] = Float64(vx[ix]) * Float64(vy[iy]) * Float64(vz[iz])
        cursor += 1
    end
    return result
end

"""
    cartesian_contracted_parent_metric_packet_dense_reference(contracted_parent; max_parent_dimension = 125)

Tiny-oracle constructor for tests. It deliberately forms dense full-parent 3D
metric matrices and therefore refuses parent dimensions above
`max_parent_dimension`.
"""
function cartesian_contracted_parent_metric_packet_dense_reference(
    contracted_parent::CartesianContractedParent3D;
    axis_metrics = nothing,
    max_parent_dimension::Integer = 125,
)
    parent = contracted_parent_basis(contracted_parent)
    dims = parent_axis_counts(parent)
    parent_dim = parent_dimension(parent)
    parent_dim <= Int(max_parent_dimension) || throw(
        ArgumentError(
            "dense contracted-parent metric oracle refuses parent dimension $parent_dim; " *
            "increase max_parent_dimension only for tiny reference tests",
        ),
    )
    metrics = _validate_axis_metric_data(
        _parent_axis_metric_data(parent; axis_metrics),
        dims,
    )
    coefficients = Matrix{Float64}(contracted_parent_coefficients(contracted_parent))
    entries = _coefficient_column_entries(contracted_parent_coefficients(contracted_parent), dims)

    overlap_parent = _dense_parent_product_matrix(
        metrics.x.overlap,
        metrics.y.overlap,
        metrics.z.overlap,
    )
    position_x_parent = _dense_parent_product_matrix(
        metrics.x.position,
        metrics.y.overlap,
        metrics.z.overlap,
    )
    position_y_parent = _dense_parent_product_matrix(
        metrics.x.overlap,
        metrics.y.position,
        metrics.z.overlap,
    )
    position_z_parent = _dense_parent_product_matrix(
        metrics.x.overlap,
        metrics.y.overlap,
        metrics.z.position,
    )
    weights_parent = _dense_parent_product_vector(
        metrics.x.weights,
        metrics.y.weights,
        metrics.z.weights,
    )
    first_moment_x_parent = _dense_parent_product_vector(
        metrics.x.centers .* metrics.x.weights,
        metrics.y.weights,
        metrics.z.weights,
    )
    first_moment_y_parent = _dense_parent_product_vector(
        metrics.x.weights,
        metrics.y.centers .* metrics.y.weights,
        metrics.z.weights,
    )
    first_moment_z_parent = _dense_parent_product_vector(
        metrics.x.weights,
        metrics.y.weights,
        metrics.z.centers .* metrics.z.weights,
    )

    overlap = Matrix{Float64}(transpose(coefficients) * overlap_parent * coefficients)
    position_x = Matrix{Float64}(transpose(coefficients) * position_x_parent * coefficients)
    position_y = Matrix{Float64}(transpose(coefficients) * position_y_parent * coefficients)
    position_z = Matrix{Float64}(transpose(coefficients) * position_z_parent * coefficients)
    weights = vec(transpose(coefficients) * weights_parent)
    first_moments = hcat(
        vec(transpose(coefficients) * first_moment_x_parent),
        vec(transpose(coefficients) * first_moment_y_parent),
        vec(transpose(coefficients) * first_moment_z_parent),
    )
    diagnostics = _packet_diagnostics(
        construction_path = :dense_reference_oracle,
        dense_parent_matrix_used = true,
        dense_reference_parent_dimension_limit = Int(max_parent_dimension),
        contracted_parent = contracted_parent,
        entries = entries,
        metrics = metrics,
        overlap = overlap,
    )
    return _metric_packet_from_matrices(
        contracted_parent,
        metrics,
        entries,
        overlap,
        position_x,
        position_y,
        position_z,
        weights,
        first_moments;
        diagnostics,
    )
end

function cartesian_contracted_parent_metric_packet_dense_reference(
    fixed_block::_NestedFixedBlock3D;
    kwargs...,
)
    return cartesian_contracted_parent_metric_packet_dense_reference(
        cartesian_contracted_parent(fixed_block);
        kwargs...,
    )
end

end

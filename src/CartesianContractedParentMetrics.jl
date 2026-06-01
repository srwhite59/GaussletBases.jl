module CartesianContractedParentMetrics

import LinearAlgebra
import SparseArrays

import ..GaussletBases: _NestedFixedBlock3D,
                         _BondAlignedDiatomicHighOrderRecipeSourceConstruction3D,
                         _CartesianNestedProjectedQShellStagedUnitDescriptor3D,
                         _CartesianNestedProductStagedByCenterUnit3D,
                         _cartesian_flat_index,
                         _cartesian_raw_product_box_plan,
                         _cartesian_raw_product_box_source_mode_indices,
                         _cartesian_unflat_index,
                         _nested_axis_lengths,
                         _nested_product_axis_function_indices,
                         _nested_product_staged_active_axis,
                         _nested_product_staged_fixed_axis,
                         _nested_projected_q_shell_staged_unit_descriptor,
                         _nested_projected_q_shell_descriptor_seed_coefficients,
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

const _PRODUCT_DOSIDE_SOURCE_BOX_SHADOW_TERMS =
    _PRODUCT_DOSIDE_SOURCE_BOX_REFERENCE_TERMS

function _product_doside_source_box_shadow_blocks(
    left_unit::_CartesianNestedProductStagedByCenterUnit3D,
    right_unit::_CartesianNestedProductStagedByCenterUnit3D,
    metrics::NamedTuple{(:x,:y,:z)};
    terms = _PRODUCT_DOSIDE_SOURCE_BOX_SHADOW_TERMS,
)
    selected_terms = Tuple(Symbol(term) for term in terms)
    !isempty(selected_terms) || throw(
        ArgumentError("product/doside source-box shadow requires at least one term"),
    )
    for term in selected_terms
        term in _PRODUCT_DOSIDE_SOURCE_BOX_SHADOW_TERMS || throw(
            ArgumentError("product/doside source-box shadow received unsupported term $(term)"),
        )
    end
    left_count = length(left_unit.column_range)
    right_count = length(right_unit.column_range)
    left_range = 1:left_count
    right_range = (left_count + 1):(left_count + right_count)
    retained_dimension = left_count + right_count
    blocks = Dict{Symbol,Matrix{Float64}}()
    component_blocks = Dict{Symbol,NamedTuple}()
    transpose_errors = Dict{Symbol,Float64}()
    for term in selected_terms
        left_left =
            _product_doside_source_box_reference_block(
                left_unit,
                left_unit,
                metrics;
                term,
            )
        left_right =
            _product_doside_source_box_reference_block(
                left_unit,
                right_unit,
                metrics;
                term,
            )
        right_left =
            _product_doside_source_box_reference_block(
                right_unit,
                left_unit,
                metrics;
                term,
            )
        right_right =
            _product_doside_source_box_reference_block(
                right_unit,
                right_unit,
                metrics;
                term,
            )
        transpose_error =
            LinearAlgebra.norm(right_left.block - transpose(left_right.block), Inf)
        block = zeros(Float64, retained_dimension, retained_dimension)
        block[left_range, left_range] .= left_left.block
        block[left_range, right_range] .= left_right.block
        block[right_range, left_range] .= right_left.block
        block[right_range, right_range] .= right_right.block
        all(isfinite, block) || throw(
            ArgumentError("product/doside source-box shadow produced non-finite entries"),
        )
        blocks[term] = block
        component_blocks[term] = (
            left_left = left_left.block,
            left_right = left_right.block,
            right_left = right_left.block,
            right_right = right_right.block,
            right_left_transpose_error = transpose_error,
        )
        transpose_errors[term] = transpose_error
    end
    max_transpose_error =
        isempty(transpose_errors) ? 0.0 : maximum(values(transpose_errors))
    return (
        path = :product_doside_source_box_shadow_blocks,
        blocks = blocks,
        component_blocks = component_blocks,
        terms = selected_terms,
        ranges = (left = left_range, right = right_range),
        retained_dimension = retained_dimension,
        transpose_errors = transpose_errors,
        diagnostics = (
            source = :product_doside_source_box_shadow_blocks,
            source_box_shadow_only = true,
            private_shadow_only = true,
            product_doside_retained_unit_plan_used = true,
            existing_product_staged_retained_helpers_authoritative = true,
            operator_factor_source = :explicit_metric_operator_data,
            operator_metric_sources =
                _cartesian_source_box_metric_sources(metrics),
            input_metric_operator_data = :caller_supplied_explicit_data,
            input_metric_operator_data_pgdg_checked = false,
            pgdg_analytic_operator_provenance_claimed = false,
            numerical_reference_fallback = false,
            component_reference_helper =
                :_product_doside_source_box_reference_block,
            right_left_block_source =
                :_product_doside_source_box_reference_block,
            symmetric_real_transpose_checked = true,
            max_right_left_transpose_error = max_transpose_error,
            product_staged_metric_execution_changed = false,
            product_doside_retained_block_math_changed = false,
            packet_adoption = false,
            fixed_block_routing = false,
            qwhamiltonian_consumes = false,
            public_default_consumes = false,
            ida_mwg_semantics_changed = false,
            retained_weight_semantics_changed = false,
            ida_weight_division_allowed = false,
            shell_projection_used = false,
            lowdin_cleanup_used = false,
            local_ecp_gaussian_mwg_implemented = false,
            generic_retained_unit_framework = false,
            supported_terms = _PRODUCT_DOSIDE_SOURCE_BOX_SHADOW_TERMS,
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
            output_finite = true,
        ),
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
        selected_count = descriptor.mode_count,
        preserves_orthogonality = true,
    )
    return (
        path = :pqs_raw_product_box_plan,
        representation = :orthogonal_raw_product_box,
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
                    source = :pqs_raw_product_box_plan,
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

function _pqs_product_source_box_axis_cross_factor(
    pqs_plan,
    product_axis,
    metrics::NamedTuple{(:x,:y,:z)},
    axis::Int,
    kind::Symbol,
)
    raw_plan = _pqs_raw_product_box_plan_view(pqs_plan)
    pqs_interval = raw_plan.axis_intervals[axis]
    product_interval = _staged_axis_interval(product_axis)
    axis_data = getproperty(metrics, (:x, :y, :z)[axis])
    hasproperty(axis_data, kind) || throw(
        ArgumentError("PQS/product source-box axis $(axis) is missing $(kind)"),
    )
    axis_operator = Matrix{Float64}(getproperty(axis_data, kind))
    size(axis_operator, 1) == size(axis_operator, 2) || throw(
        ArgumentError("PQS/product source-box axis $(axis) $(kind) matrix must be square"),
    )
    first(pqs_interval) >= 1 && last(pqs_interval) <= size(axis_operator, 1) ||
        throw(
            ArgumentError("PQS/product source-box PQS interval exceeds axis $(axis) $(kind) matrix"),
        )
    first(product_interval) >= 1 &&
        last(product_interval) <= size(axis_operator, 2) ||
        throw(
            ArgumentError("PQS/product source-box product interval exceeds axis $(axis) $(kind) matrix"),
        )
    pqs_coefficients =
        Matrix{Float64}(raw_plan.axis_local_coefficients[axis])
    product_coefficients = Matrix{Float64}(product_axis.coefficient_matrix)
    size(pqs_coefficients, 1) == length(pqs_interval) || throw(
        DimensionMismatch("PQS/product source-box PQS axis coefficients must match interval length"),
    )
    size(product_coefficients, 1) == length(product_interval) || throw(
        DimensionMismatch("PQS/product source-box product axis coefficients must match interval length"),
    )
    @views return Matrix{Float64}(
        transpose(pqs_coefficients) *
        axis_operator[pqs_interval, product_interval] *
        product_coefficients,
    )
end

function _pqs_product_source_box_cross_factors(
    pqs_plan,
    product_unit::_CartesianNestedProductStagedByCenterUnit3D,
    metrics::NamedTuple{(:x,:y,:z)},
)
    return ntuple(axis -> (
        overlap = _pqs_product_source_box_axis_cross_factor(
            pqs_plan,
            product_unit.axes[axis],
            metrics,
            axis,
            :overlap,
        ),
        position = _pqs_product_source_box_axis_cross_factor(
            pqs_plan,
            product_unit.axes[axis],
            metrics,
            axis,
            :position,
        ),
        x2 = _pqs_product_source_box_axis_cross_factor(
            pqs_plan,
            product_unit.axes[axis],
            metrics,
            axis,
            :x2,
        ),
        kinetic = _pqs_product_source_box_axis_cross_factor(
            pqs_plan,
            product_unit.axes[axis],
            metrics,
            axis,
            :kinetic,
        ),
    ), 3)
end

const _PQS_PRODUCT_SOURCE_BOX_REFERENCE_TERMS = (
    :overlap,
    :position_x,
    :position_y,
    :position_z,
    :x2_x,
    :x2_y,
    :x2_z,
    :kinetic,
)

function _pqs_product_source_box_pair_plan(
    pqs_plan,
    product_unit::_CartesianNestedProductStagedByCenterUnit3D,
    metrics::NamedTuple{(:x,:y,:z)},
)
    raw_plan = _pqs_raw_product_box_plan_view(pqs_plan)
    raw_plan.representation == :orthogonal_raw_product_box ||
        throw(
            ArgumentError("PQS/product source-box pair plan requires a raw product-box PQS plan"),
        )
    product_unit.kind == :product_doside || throw(
        ArgumentError("PQS/product source-box pair plan requires a product_doside unit"),
    )
    product_retained_unit_plan =
        _product_doside_retained_unit_plan(product_unit)
    cross_factors =
        _pqs_product_source_box_cross_factors(raw_plan, product_unit, metrics)
    product_axis_intervals = product_retained_unit_plan.source_axis_intervals
    product_axis_source_lengths = product_retained_unit_plan.source_axis_lengths
    pqs_axis_centers = ntuple(axis -> begin
        interval = raw_plan.axis_intervals[axis]
        Float64.(getproperty(getproperty(metrics, (:x, :y, :z)[axis]), :centers)[interval])
    end, 3)
    product_axis_centers = ntuple(axis -> begin
        interval = product_axis_intervals[axis]
        Float64.(getproperty(getproperty(metrics, (:x, :y, :z)[axis]), :centers)[interval])
    end, 3)
    return (
        pair_kind = :pqs_product_source_box,
        left_source_family = :mode_selected_raw_product_box,
        right_source_family = :product_doside,
        left_source_dimensions = raw_plan.source_mode_dims,
        right_source_dimensions = product_axis_source_lengths,
        left_source_dimension = raw_plan.source_mode_count,
        right_source_dimension = prod(product_axis_source_lengths),
        left_retained_count = raw_plan.boundary_selector.selected_count,
        right_retained_count = product_retained_unit_plan.retained_count,
        axis_intervals = (
            pqs = raw_plan.axis_intervals,
            product = product_axis_intervals,
        ),
        axis_centers = (
            pqs = pqs_axis_centers,
            product = product_axis_centers,
        ),
        pqs_boundary_mode_selector = raw_plan.boundary_selector,
        product_retained_unit_plan = product_retained_unit_plan,
        product_retained_transform = product_retained_unit_plan,
        one_dimensional_cross_factors = (
            x = cross_factors[1],
            y = cross_factors[2],
            z = cross_factors[3],
        ),
        supported_terms = _PQS_PRODUCT_SOURCE_BOX_REFERENCE_TERMS,
        diagnostics = (
            source = :pqs_product_source_box_pair_plan,
            private_shadow_only = true,
            pqs_representation = :mode_selected_raw_product_box,
            raw_product_box_plan_used = true,
            pqs_raw_product_box_plan_used = true,
            shared_raw_product_box_plan_available =
                raw_plan.shared_raw_product_box_plan_available,
            shared_raw_product_box_plan_used =
                raw_plan.shared_raw_product_box_plan_used,
            source_mode_ordering = raw_plan.source_mode_ordering,
            pqs_boundary_mode_selection_used = true,
            product_doside_retained_transform_used = true,
            product_doside_retained_unit_plan_used = true,
            raw_product_box_operators_use_1d_factors = true,
            operator_factor_source = :explicit_axis_metric_data,
            operator_metric_sources =
                _cartesian_source_box_metric_sources(metrics),
            raw_product_box_numerical_reference_fallback =
                hasproperty(raw_plan.diagnostics, :raw_product_box_numerical_reference_fallback) ?
                raw_plan.diagnostics.raw_product_box_numerical_reference_fallback :
                nothing,
            shell_projection_used = false,
            lowdin_cleanup_used = false,
            support_coefficient_matrix_used = false,
            support_local_pqs_oracle_used = false,
            retained_pqs_weights_used = false,
            retained_pqs_weights_positive_checked = false,
            ida_weight_division_allowed = false,
            packet_adoption = false,
            fixed_block_routing = false,
            qwhamiltonian_consumes = false,
            public_default_consumes = false,
            cr2_science_status_changed = false,
            local_ecp_gaussian_mwg_implemented = false,
            generic_retained_unit_framework = false,
            dense_raw_source_box_pair_matrix_materialized = false,
            dense_raw_pair_storage_avoided = true,
            retained_block_assembled_directly_from_1d_factors = true,
            source_box_pair_storage_scaling =
                :one_dimensional_factors_plus_retained_block,
        ),
    )
end

function _source_box_separable_term_factor_kinds(term::Symbol)
    term == :overlap && return ((:overlap, :overlap, :overlap),)
    term == :position_x && return ((:position, :overlap, :overlap),)
    term == :position_y && return ((:overlap, :position, :overlap),)
    term == :position_z && return ((:overlap, :overlap, :position),)
    term == :x2_x && return ((:x2, :overlap, :overlap),)
    term == :x2_y && return ((:overlap, :x2, :overlap),)
    term == :x2_z && return ((:overlap, :overlap, :x2),)
    term == :kinetic && return (
        (:kinetic, :overlap, :overlap),
        (:overlap, :kinetic, :overlap),
        (:overlap, :overlap, :kinetic),
    )
    throw(
        ArgumentError("source-box separable reference block received unsupported term $(term)"),
    )
end

function _pqs_product_source_box_factor(
    pair_plan,
    axis::Int,
    kind::Symbol,
)
    axis_factors = getproperty(pair_plan.one_dimensional_cross_factors, (:x, :y, :z)[axis])
    return getproperty(axis_factors, kind)
end

function _pqs_product_source_box_block_from_factors(pair_plan, term::Symbol)
    factor_terms = _source_box_separable_term_factor_kinds(term)
    block = zeros(
        Float64,
        pair_plan.left_retained_count,
        pair_plan.right_retained_count,
    )
    pqs_modes = pair_plan.pqs_boundary_mode_selector.mode_indices
    product_modes = pair_plan.product_retained_transform.axis_function_indices
    @inbounds for factor_kinds in factor_terms
        fx = _pqs_product_source_box_factor(pair_plan, 1, factor_kinds[1])
        fy = _pqs_product_source_box_factor(pair_plan, 2, factor_kinds[2])
        fz = _pqs_product_source_box_factor(pair_plan, 3, factor_kinds[3])
        for col in eachindex(product_modes)
            px, py, pz = product_modes[col]
            for row in eachindex(pqs_modes)
                qx, qy, qz = pqs_modes[row]
                block[row, col] += fx[qx, px] * fy[qy, py] * fz[qz, pz]
            end
        end
    end
    return block
end

function _pqs_product_source_box_reference_blocks_from_pair_plan(
    pair_plan;
    terms = pair_plan.supported_terms,
)
    pair_plan.pair_kind == :pqs_product_source_box || throw(
        ArgumentError("PQS/product source-box reference blocks require a PQS/product pair plan"),
    )
    selected_terms = Tuple(Symbol(term) for term in terms)
    !isempty(selected_terms) || throw(
        ArgumentError("PQS/product source-box reference blocks require at least one term"),
    )
    for term in selected_terms
        term in pair_plan.supported_terms || throw(
            ArgumentError("PQS/product source-box reference blocks received unsupported term $(term)"),
        )
    end
    blocks = Dict{Symbol,Matrix{Float64}}()
    for term in selected_terms
        block = _pqs_product_source_box_block_from_factors(pair_plan, term)
        all(isfinite, block) || throw(
            ArgumentError("PQS/product source-box reference blocks produced non-finite entries for $(term)"),
        )
        blocks[term] = block
    end
    return (
        path = :pqs_product_source_box_reference_blocks,
        terms = selected_terms,
        blocks = blocks,
        pair_plan = pair_plan,
        diagnostics = merge(
            pair_plan.diagnostics,
            (
                source = :pqs_product_source_box_reference_blocks,
                pair_plan_reused_for_terms = true,
                pair_plan_reuse_term_count = length(selected_terms),
                dense_raw_source_box_pair_matrix_materialized = false,
                dense_raw_pair_storage_avoided = true,
                retained_block_assembled_directly_from_1d_factors = true,
                source_box_pair_storage_scaling =
                    :one_dimensional_factors_plus_retained_block,
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
            ),
        ),
    )
end

function _pqs_product_source_box_reference_blocks(
    pqs_plan,
    product_unit::_CartesianNestedProductStagedByCenterUnit3D,
    metrics::NamedTuple{(:x,:y,:z)};
    terms = _PQS_PRODUCT_SOURCE_BOX_REFERENCE_TERMS,
)
    pair_plan = _pqs_product_source_box_pair_plan(pqs_plan, product_unit, metrics)
    return _pqs_product_source_box_reference_blocks_from_pair_plan(
        pair_plan;
        terms,
    )
end

function _pqs_product_source_box_reference_block(
    pqs_plan,
    product_unit::_CartesianNestedProductStagedByCenterUnit3D,
    metrics::NamedTuple{(:x,:y,:z)};
    term::Symbol,
)
    pair_plan = _pqs_product_source_box_pair_plan(pqs_plan, product_unit, metrics)
    term in pair_plan.supported_terms || throw(
        ArgumentError("PQS/product source-box reference block received unsupported term $(term)"),
    )
    block = _pqs_product_source_box_block_from_factors(pair_plan, term)
    all(isfinite, block) || throw(
        ArgumentError("PQS/product source-box reference block produced non-finite entries"),
    )
    return (
        path = :pqs_product_source_box_reference,
        term = term,
        block = block,
        pair_plan = pair_plan,
        diagnostics = merge(
            pair_plan.diagnostics,
            (
                source = :pqs_product_source_box_reference_block,
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
                pair_plan_reused_for_terms = false,
                dense_raw_source_box_pair_matrix_materialized = false,
                dense_raw_pair_storage_avoided = true,
                retained_block_assembled_directly_from_1d_factors = true,
                source_box_pair_storage_scaling =
                    :one_dimensional_factors_plus_retained_block,
            ),
        ),
    )
end

const _PQS_PQS_SOURCE_BOX_REFERENCE_TERMS =
    _PQS_PRODUCT_SOURCE_BOX_REFERENCE_TERMS

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

function _pqs_pqs_source_box_axis_cross_factor(
    left_raw_plan,
    right_raw_plan,
    metrics::NamedTuple{(:x,:y,:z)},
    axis::Int,
    kind::Symbol,
)
    axis_name = (:x, :y, :z)[axis]
    axis_data = getproperty(metrics, axis_name)
    hasproperty(axis_data, kind) || throw(
        ArgumentError("PQS/PQS source-box axis $(axis_name) is missing $(kind)"),
    )
    axis_operator = Matrix{Float64}(getproperty(axis_data, kind))
    size(axis_operator, 1) == size(axis_operator, 2) || throw(
        ArgumentError("PQS/PQS source-box axis $(axis_name) $(kind) matrix must be square"),
    )
    left_interval = left_raw_plan.axis_intervals[axis]
    right_interval = right_raw_plan.axis_intervals[axis]
    first(left_interval) >= 1 && last(left_interval) <= size(axis_operator, 1) ||
        throw(
            ArgumentError("PQS/PQS source-box left interval exceeds axis $(axis_name) $(kind) matrix"),
        )
    first(right_interval) >= 1 && last(right_interval) <= size(axis_operator, 2) ||
        throw(
            ArgumentError("PQS/PQS source-box right interval exceeds axis $(axis_name) $(kind) matrix"),
        )
    left_coefficients = Matrix{Float64}(left_raw_plan.axis_local_coefficients[axis])
    right_coefficients = Matrix{Float64}(right_raw_plan.axis_local_coefficients[axis])
    size(left_coefficients, 1) == length(left_interval) || throw(
        DimensionMismatch("PQS/PQS source-box left axis coefficients must match interval length"),
    )
    size(right_coefficients, 1) == length(right_interval) || throw(
        DimensionMismatch("PQS/PQS source-box right axis coefficients must match interval length"),
    )
    @views return Matrix{Float64}(
        transpose(left_coefficients) *
        axis_operator[left_interval, right_interval] *
        right_coefficients,
    )
end

function _pqs_pqs_source_box_cross_factors(
    left_raw_plan,
    right_raw_plan,
    metrics::NamedTuple{(:x,:y,:z)},
)
    return ntuple(axis -> (
        overlap = _pqs_pqs_source_box_axis_cross_factor(
            left_raw_plan,
            right_raw_plan,
            metrics,
            axis,
            :overlap,
        ),
        position = _pqs_pqs_source_box_axis_cross_factor(
            left_raw_plan,
            right_raw_plan,
            metrics,
            axis,
            :position,
        ),
        x2 = _pqs_pqs_source_box_axis_cross_factor(
            left_raw_plan,
            right_raw_plan,
            metrics,
            axis,
            :x2,
        ),
        kinetic = _pqs_pqs_source_box_axis_cross_factor(
            left_raw_plan,
            right_raw_plan,
            metrics,
            axis,
            :kinetic,
        ),
    ), 3)
end

function _pqs_pqs_source_box_pair_plan(
    left_pqs_plan,
    right_pqs_plan,
    metrics::NamedTuple{(:x,:y,:z)};
    same_plan_atol::Real = 1.0e-12,
)
    left_raw_plan = _pqs_raw_product_box_plan_view(left_pqs_plan)
    right_raw_plan = _pqs_raw_product_box_plan_view(right_pqs_plan)
    left_raw_plan.representation == :orthogonal_raw_product_box ||
        throw(ArgumentError("PQS/PQS source-box pair plan requires a raw product-box left PQS plan"))
    right_raw_plan.representation == :orthogonal_raw_product_box ||
        throw(ArgumentError("PQS/PQS source-box pair plan requires a raw product-box right PQS plan"))
    left_raw_plan.operator_factors_available && right_raw_plan.operator_factors_available ||
        throw(ArgumentError("PQS/PQS source-box pair plan requires operator-backed raw plans"))
    same_raw_product_box_plan = _same_pqs_raw_product_box_plan(
        left_raw_plan,
        right_raw_plan;
        atol = same_plan_atol,
    )
    _compatible_cross_pqs_raw_product_box_plans(left_raw_plan, right_raw_plan) ||
        throw(
            ArgumentError("PQS/PQS source-box pair plan requires equal source-mode dimensions, ordering, and boundary selectors"),
        )
    source_product_modes_orthogonal =
        left_raw_plan.source_product_modes_orthogonal === true &&
        right_raw_plan.source_product_modes_orthogonal === true
    source_product_modes_orthogonal || throw(
        ArgumentError("PQS/PQS source-box pair plan requires orthogonal source product modes"),
    )
    left_axis_centers = ntuple(axis -> begin
        interval = left_raw_plan.axis_intervals[axis]
        Float64.(getproperty(getproperty(metrics, (:x, :y, :z)[axis]), :centers)[interval])
    end, 3)
    right_axis_centers = ntuple(axis -> begin
        interval = right_raw_plan.axis_intervals[axis]
        Float64.(getproperty(getproperty(metrics, (:x, :y, :z)[axis]), :centers)[interval])
    end, 3)
    cross_factors =
        _pqs_pqs_source_box_cross_factors(left_raw_plan, right_raw_plan, metrics)
    return (
        pair_kind = :pqs_pqs_source_box,
        left_source_family = :mode_selected_raw_product_box,
        right_source_family = :mode_selected_raw_product_box,
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
        axis_centers = (
            left = left_axis_centers,
            right = right_axis_centers,
        ),
        left_raw_product_box_plan = left_raw_plan,
        right_raw_product_box_plan = right_raw_plan,
        left_boundary_mode_selector = left_raw_plan.boundary_selector,
        right_boundary_mode_selector = right_raw_plan.boundary_selector,
        one_dimensional_cross_factors = (
            x = cross_factors[1],
            y = cross_factors[2],
            z = cross_factors[3],
        ),
        supported_terms = _PQS_PQS_SOURCE_BOX_REFERENCE_TERMS,
        diagnostics = (
            source = :pqs_pqs_source_box_pair_plan,
            pair_kind = :pqs_pqs_source_box,
            private_shadow_only = true,
            self_same_plan_only = same_raw_product_box_plan,
            cross_pqs_inputs_supported = !same_raw_product_box_plan,
            same_raw_product_box_plan = same_raw_product_box_plan,
            equal_source_mode_dims_required = true,
            pqs_representation = :mode_selected_raw_product_box,
            raw_product_box_plan_used = true,
            source_mode_ordering = left_raw_plan.source_mode_ordering,
            left_source_mode_dims = left_raw_plan.source_mode_dims,
            right_source_mode_dims = right_raw_plan.source_mode_dims,
            left_retained_count =
                left_raw_plan.boundary_selector.selected_count,
            right_retained_count =
                right_raw_plan.boundary_selector.selected_count,
            left_boundary_mode_selection_used = true,
            right_boundary_mode_selection_used = true,
            raw_product_box_operators_use_1d_factors = true,
            operator_factor_source = :explicit_axis_metric_data,
            operator_metric_sources =
                _cartesian_source_box_metric_sources(metrics),
            raw_product_box_numerical_reference_fallback =
                hasproperty(left_raw_plan.diagnostics, :raw_product_box_numerical_reference_fallback) ?
                left_raw_plan.diagnostics.raw_product_box_numerical_reference_fallback :
                nothing,
            shell_projection_used = false,
            lowdin_cleanup_used = false,
            support_coefficient_matrix_used = false,
            support_local_pqs_oracle_used = false,
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
            dense_raw_source_box_pair_matrix_materialized = false,
            dense_raw_pair_storage_avoided = true,
            retained_block_assembled_directly_from_1d_factors = true,
            source_box_pair_storage_scaling =
                :one_dimensional_factors_plus_retained_block,
        ),
    )
end

function _pqs_pqs_source_box_factor(
    pair_plan,
    axis::Int,
    kind::Symbol,
)
    axis_factors = getproperty(pair_plan.one_dimensional_cross_factors, (:x, :y, :z)[axis])
    return getproperty(axis_factors, kind)
end

function _pqs_pqs_source_box_block_from_factors(pair_plan, term::Symbol)
    term in pair_plan.supported_terms || throw(
        ArgumentError("PQS/PQS source-box reference blocks received unsupported term $(term)"),
    )
    block = zeros(
        Float64,
        pair_plan.left_retained_count,
        pair_plan.right_retained_count,
    )
    left_modes = pair_plan.left_boundary_mode_selector.mode_indices
    right_modes = pair_plan.right_boundary_mode_selector.mode_indices
    @inbounds for factor_kinds in _source_box_separable_term_factor_kinds(term)
        fx = _pqs_pqs_source_box_factor(pair_plan, 1, factor_kinds[1])
        fy = _pqs_pqs_source_box_factor(pair_plan, 2, factor_kinds[2])
        fz = _pqs_pqs_source_box_factor(pair_plan, 3, factor_kinds[3])
        for col in eachindex(right_modes)
            rx, ry, rz = right_modes[col]
            for row in eachindex(left_modes)
                lx, ly, lz = left_modes[row]
                block[row, col] += fx[lx, rx] * fy[ly, ry] * fz[lz, rz]
            end
        end
    end
    return block
end

function _pqs_pqs_source_box_reference_blocks_from_pair_plan(
    pair_plan;
    terms = pair_plan.supported_terms,
    atol::Real = 1.0e-10,
)
    pair_plan.pair_kind == :pqs_pqs_source_box || throw(
        ArgumentError("PQS/PQS source-box reference blocks require a PQS/PQS pair plan"),
    )
    selected_terms = Tuple(Symbol(term) for term in terms)
    !isempty(selected_terms) || throw(
        ArgumentError("PQS/PQS source-box reference blocks require at least one term"),
    )
    for term in selected_terms
        term in pair_plan.supported_terms || throw(
            ArgumentError("PQS/PQS source-box reference blocks received unsupported term $(term)"),
        )
    end
    blocks = Dict{Symbol,Matrix{Float64}}()
    raw_box_reference_blocks = Dict{Symbol,Matrix{Float64}}()
    block_errors = Dict{Symbol,Float64}()
    same_raw_product_box_plan = Bool(pair_plan.diagnostics.same_raw_product_box_plan)
    for term in selected_terms
        block = _pqs_pqs_source_box_block_from_factors(pair_plan, term)
        if same_raw_product_box_plan
            raw_reference =
                _pqs_raw_product_box_reference_block(
                    pair_plan.left_raw_product_box_plan;
                    term,
                ).block
            block_error = LinearAlgebra.norm(block - raw_reference, Inf)
            block_error <= atol || throw(
                ArgumentError("PQS/PQS source-box reference block disagrees with raw-box self reference"),
            )
            raw_box_reference_blocks[term] = raw_reference
            block_errors[term] = block_error
        end
        all(isfinite, block) || throw(
            ArgumentError("PQS/PQS source-box reference blocks produced non-finite entries for $(term)"),
        )
        blocks[term] = block
    end
    max_block_error =
        isempty(block_errors) ? nothing : maximum(values(block_errors))
    return (
        path = :pqs_pqs_source_box_reference_blocks,
        terms = selected_terms,
        blocks = blocks,
        raw_box_reference_blocks = raw_box_reference_blocks,
        block_errors = block_errors,
        max_block_error = max_block_error,
        pair_plan = pair_plan,
        diagnostics = merge(
            pair_plan.diagnostics,
            (
                source = :pqs_pqs_source_box_reference_blocks,
                pair_plan_reused_for_terms = true,
                pair_plan_reuse_term_count = length(selected_terms),
                raw_box_self_reference_compared = same_raw_product_box_plan,
                raw_box_self_reference_helper =
                    same_raw_product_box_plan ?
                    :_pqs_raw_product_box_reference_block : nothing,
                explicit_source_box_oracle_tested = false,
                max_block_error = max_block_error,
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
            ),
        ),
    )
end

function _pqs_pqs_source_box_reference_blocks(
    left_pqs_plan,
    right_pqs_plan,
    metrics::NamedTuple{(:x,:y,:z)};
    terms = _PQS_PQS_SOURCE_BOX_REFERENCE_TERMS,
)
    pair_plan = _pqs_pqs_source_box_pair_plan(
        left_pqs_plan,
        right_pqs_plan,
        metrics,
    )
    return _pqs_pqs_source_box_reference_blocks_from_pair_plan(
        pair_plan;
        terms,
    )
end

function _pqs_pqs_source_box_reference_block(
    left_pqs_plan,
    right_pqs_plan,
    metrics::NamedTuple{(:x,:y,:z)};
    term::Symbol,
)
    reference = _pqs_pqs_source_box_reference_blocks(
        left_pqs_plan,
        right_pqs_plan,
        metrics;
        terms = (term,),
    )
    return (
        path = :pqs_pqs_source_box_reference,
        term = term,
        block = reference.blocks[term],
        raw_box_reference_block = get(reference.raw_box_reference_blocks, term, nothing),
        block_error = get(reference.block_errors, term, nothing),
        pair_plan = reference.pair_plan,
        diagnostics = merge(
            reference.diagnostics,
            (
                source = :pqs_pqs_source_box_reference_block,
                pair_plan_reused_for_terms = false,
            ),
        ),
    )
end

const _PQS_PRODUCT_SOURCE_BOX_SHADOW_TERMS =
    _PQS_PRODUCT_SOURCE_BOX_REFERENCE_TERMS

function _pqs_pqs_product_supported_safe_terms(terms)
    selected_terms = Tuple(Symbol(term) for term in terms)
    !isempty(selected_terms) || throw(
        ArgumentError("PQS/PQS/product source-box route requires at least one term"),
    )
    for term in selected_terms
        term in _PQS_PRODUCT_SOURCE_BOX_SHADOW_TERMS || throw(
            ArgumentError("PQS/PQS/product source-box route received unsupported term $(term)"),
        )
    end
    return selected_terms
end

function _pqs_product_source_box_product_block(
    product_unit::_CartesianNestedProductStagedByCenterUnit3D,
    metrics::NamedTuple{(:x,:y,:z)},
    term::Symbol,
)
    term == :kinetic && return _product_doside_retained_kinetic_block(
        product_unit,
        product_unit,
        metrics,
    )
    return _product_doside_retained_low_order_block(
        product_unit,
        product_unit,
        metrics;
        term,
    )
end

function _pqs_product_source_box_all_pairs_inventory(
    raw_plan,
    product_unit::_CartesianNestedProductStagedByCenterUnit3D,
    pqs_product_pair_plan,
    pqs_range::UnitRange{Int},
    product_range::UnitRange{Int},
    supported_terms::Tuple{Vararg{Symbol}},
)
    product_retained_unit_plan = pqs_product_pair_plan.product_retained_unit_plan
    retained_units = (
        (
            unit_key = :pqs,
            retained_unit_kind = :pqs,
            source_family = :mode_selected_raw_product_box,
            retained_range = pqs_range,
            source_dimensions = raw_plan.source_mode_dims,
            source_dimension = raw_plan.source_mode_count,
            retained_count = raw_plan.boundary_selector.selected_count,
            supported_safe_terms = supported_terms,
        ),
        (
            unit_key = :product,
            retained_unit_kind = :product_doside,
            source_family = :product_doside,
            retained_range = product_range,
            source_dimensions = product_retained_unit_plan.source_axis_lengths,
            source_dimension = product_retained_unit_plan.source_dimension,
            retained_count = product_retained_unit_plan.retained_count,
            supported_safe_terms = supported_terms,
        ),
    )
    pair_entries = (
        (
            pair_key = (:pqs, :pqs),
            pair_kind = :pqs_pqs_source_box,
            block_helper = :_pqs_pqs_source_box_reference_blocks,
            upper_triangular = true,
            transpose_only_lower_block = false,
        ),
        (
            pair_key = (:pqs, :product),
            pair_kind = :pqs_product_source_box,
            block_helper = :_pqs_product_source_box_reference_blocks,
            upper_triangular = true,
            transpose_only_lower_block = true,
        ),
        (
            pair_key = (:product, :product),
            pair_kind = :product_doside_source_box_pair,
            block_helper = :product_doside_retained_block_helpers,
            upper_triangular = true,
            transpose_only_lower_block = false,
        ),
    )
    return (
        object_kind = :pqs_product_source_box_all_pairs_inventory,
        retained_units = retained_units,
        pair_entries = pair_entries,
        supported_terms = supported_terms,
        diagnostics = (
            source = :pqs_product_source_box_all_pairs_inventory,
            all_pairs_inventory_private = true,
            pair_inventory_complete_for_units = (:pqs, :product),
            retained_unit_count = length(retained_units),
            upper_triangular_pair_count = length(pair_entries),
            expected_upper_triangular_pair_count = 3,
            pair_keys = map(entry -> entry.pair_key, pair_entries),
            block_helpers = map(entry -> entry.block_helper, pair_entries),
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
            product_lower_triangle_transpose_only = true,
        ),
    )
end

function _pqs_product_source_box_shadow_blocks(
    pqs_plan,
    product_unit::_CartesianNestedProductStagedByCenterUnit3D,
    metrics::NamedTuple{(:x,:y,:z)};
    terms = _PQS_PRODUCT_SOURCE_BOX_SHADOW_TERMS,
)
    raw_plan = _pqs_raw_product_box_plan_view(pqs_plan)
    raw_plan.representation == :orthogonal_raw_product_box ||
        throw(
            ArgumentError("PQS/product source-box shadow requires a raw product-box PQS plan"),
        )
    selected_terms = Tuple(Symbol(term) for term in terms)
    !isempty(selected_terms) || throw(
        ArgumentError("PQS/product source-box shadow requires at least one term"),
    )
    for term in selected_terms
        term in _PQS_PRODUCT_SOURCE_BOX_SHADOW_TERMS || throw(
            ArgumentError("PQS/product source-box shadow received unsupported term $(term)"),
        )
    end
    product_unit.kind == :product_doside || throw(
        ArgumentError("PQS/product source-box shadow requires a product_doside unit"),
    )
    _require_product_doside_retained_block_unit(product_unit; side = :right)

    pqs_count = raw_plan.boundary_selector.selected_count
    product_count = length(product_unit.column_range)
    pqs_range = 1:pqs_count
    product_range = (pqs_count + 1):(pqs_count + product_count)
    retained_dimension = pqs_count + product_count
    blocks = Dict{Symbol,Matrix{Float64}}()
    component_blocks = Dict{Symbol,NamedTuple}()
    pqs_pqs_references = _pqs_pqs_source_box_reference_blocks(
        raw_plan,
        raw_plan,
        metrics;
        terms = selected_terms,
    )
    pqs_product_references = _pqs_product_source_box_reference_blocks(
        raw_plan,
        product_unit,
        metrics;
        terms = selected_terms,
    )
    all_pairs_inventory = _pqs_product_source_box_all_pairs_inventory(
        raw_plan,
        product_unit,
        pqs_product_references.pair_plan,
        pqs_range,
        product_range,
        selected_terms,
    )
    for term in selected_terms
        pqs_self = pqs_pqs_references.blocks[term]
        pqs_product = pqs_product_references.blocks[term]
        product_self =
            _pqs_product_source_box_product_block(product_unit, metrics, term)
        block = zeros(Float64, retained_dimension, retained_dimension)
        block[pqs_range, pqs_range] .= pqs_self
        block[pqs_range, product_range] .= pqs_product
        block[product_range, pqs_range] .= transpose(pqs_product)
        block[product_range, product_range] .= product_self
        all(isfinite, block) || throw(
            ArgumentError("PQS/product source-box shadow produced non-finite entries"),
        )
        blocks[term] = block
        component_blocks[term] = (
            pqs_pqs = pqs_self,
            pqs_product = pqs_product,
            product_pqs = transpose(pqs_product),
            product_product = product_self,
        )
    end
    return (
        path = :pqs_product_source_box_shadow_blocks,
        blocks = blocks,
        component_blocks = component_blocks,
        terms = selected_terms,
        ranges = (pqs = pqs_range, product = product_range),
        retained_dimension = retained_dimension,
        all_pairs_inventory = all_pairs_inventory,
        pqs_pqs_reference_blocks = pqs_pqs_references,
        pqs_product_reference_blocks = pqs_product_references,
        diagnostics = (
            source = :pqs_product_source_box_shadow_blocks,
            raw_plan_first_path = true,
            descriptor_wrapper = false,
            source_box_shadow_only = true,
            private_shadow_only = true,
            all_pairs_inventory_private = true,
            pair_inventory_complete_for_units = (:pqs, :product),
            all_pairs_inventory_pair_count =
                length(all_pairs_inventory.pair_entries),
            packet_adoption = false,
            fixed_block_routing = false,
            pqs_representation = :mode_selected_raw_product_box,
            shared_raw_product_box_plan_available =
                raw_plan.shared_raw_product_box_plan_available,
            shared_raw_product_box_plan_used =
                raw_plan.shared_raw_product_box_plan_used,
            source_mode_ordering = raw_plan.source_mode_ordering,
            product_doside_retained_transform_used = true,
            pqs_pqs_block_source = :pqs_pqs_source_box_reference_blocks,
            pqs_pqs_raw_box_self_reference_compared =
                pqs_pqs_references.diagnostics.raw_box_self_reference_compared,
            pqs_product_block_source = :pqs_product_source_box_reference_blocks,
            pqs_product_pair_plan_reused_for_terms = true,
            pair_plan_reused_for_terms = true,
            product_product_block_source = :product_doside_retained_block_helpers,
            reverse_pqs_product_transpose_only = true,
            shell_projection_used = false,
            lowdin_cleanup_used = false,
            support_coefficient_matrix_used = false,
            support_local_pqs_oracle_used = false,
            retained_weight_semantics = :not_positive_quadrature_weights,
            retained_pqs_weights_used = false,
            retained_pqs_weights_positive_checked = false,
            ida_weight_division_allowed = false,
            qwhamiltonian_consumes = false,
            public_default_consumes = false,
            cr2_science_status_changed = false,
            local_ecp_gaussian_mwg_implemented = false,
            generic_retained_unit_framework = false,
            supported_terms = _PQS_PRODUCT_SOURCE_BOX_SHADOW_TERMS,
            dense_raw_source_box_pair_matrix_materialized = false,
            dense_raw_pair_storage_avoided = true,
            retained_block_assembled_directly_from_1d_factors = true,
            source_box_pair_storage_scaling =
                :one_dimensional_factors_plus_retained_block,
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

function _pqs_product_source_box_shadow_blocks(
    descriptor::_CartesianNestedProjectedQShellStagedUnitDescriptor3D,
    pqs_plan,
    product_unit::_CartesianNestedProductStagedByCenterUnit3D,
    metrics::NamedTuple{(:x,:y,:z)};
    terms = _PQS_PRODUCT_SOURCE_BOX_SHADOW_TERMS,
)
    descriptor.kind == :projected_q_shell || throw(
        ArgumentError("PQS/product source-box shadow requires a projected_q_shell descriptor"),
    )
    raw_plan = _pqs_raw_product_box_plan_view(pqs_plan)
    raw_plan.boundary_selector.selected_count == descriptor.mode_count || throw(
        DimensionMismatch("PQS/product source-box shadow PQS plan count must match descriptor mode count"),
    )
    raw_plan.boundary_selector.mode_indices == descriptor.boundary_mode_indices || throw(
        ArgumentError("PQS/product source-box shadow raw plan boundary modes disagree with descriptor"),
    )
    raw_plan.boundary_selector.column_indices == descriptor.boundary_column_indices || throw(
        ArgumentError("PQS/product source-box shadow raw plan boundary columns disagree with descriptor"),
    )
    raw_plan.source_mode_dims ==
        ntuple(axis -> size(descriptor.axis_local_coefficients[axis], 2), 3) ||
        throw(
            DimensionMismatch("PQS/product source-box shadow raw plan source-mode dimensions disagree with descriptor"),
        )
    return _pqs_product_source_box_shadow_blocks(
        pqs_plan,
        product_unit,
        metrics;
        terms = terms,
    )
end

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
        ),
        (
            pair_key = (:pqs_left, :pqs_right),
            pair_kind = :pqs_pqs_source_box,
            block_helper = :_pqs_pqs_source_box_reference_blocks,
            upper_triangular = true,
            transpose_only_lower_block = true,
        ),
        (
            pair_key = (:pqs_left, :product),
            pair_kind = :pqs_product_source_box,
            block_helper = :_pqs_product_source_box_reference_blocks,
            upper_triangular = true,
            transpose_only_lower_block = true,
        ),
        (
            pair_key = (:pqs_right, :pqs_right),
            pair_kind = :pqs_pqs_source_box,
            block_helper = :_pqs_pqs_source_box_reference_blocks,
            upper_triangular = true,
            transpose_only_lower_block = false,
        ),
        (
            pair_key = (:pqs_right, :product),
            pair_kind = :pqs_product_source_box,
            block_helper = :_pqs_product_source_box_reference_blocks,
            upper_triangular = true,
            transpose_only_lower_block = true,
        ),
        (
            pair_key = (:product, :product),
            pair_kind = :product_doside_source_box_pair,
            block_helper = :product_doside_retained_block_helpers,
            upper_triangular = true,
            transpose_only_lower_block = false,
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
        support_count == build.retained_count == fact.retained_count || throw(
            ArgumentError("PQS atom-box retained count mismatch"),
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
                support_dense_direct_support_unit = true,
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
                retained_count_match = length(build.column_range) == support_count,
                column_range = build.column_range,
                coefficient_matrix_matches_direct_support =
                    parent_coefficient_error == 0.0,
                max_parent_coefficient_error = parent_coefficient_error,
                local_identity_error = local_identity_error,
                local_support_coefficient_shape = size(local_coefficients),
            ),
        )
    end

    diagnostics = (
        source = :pqs_atom_box_support_dense_units,
        atom_box_only = true,
        support_dense_direct_support_units_created = true,
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

function _pqs_pqs_product_source_box_shadow_blocks(
    left_pqs_plan,
    right_pqs_plan,
    product_unit::_CartesianNestedProductStagedByCenterUnit3D,
    metrics::NamedTuple{(:x,:y,:z)};
    terms = _PQS_PRODUCT_SOURCE_BOX_SHADOW_TERMS,
)
    left_raw_plan = _pqs_raw_product_box_plan_view(left_pqs_plan)
    right_raw_plan = _pqs_raw_product_box_plan_view(right_pqs_plan)
    left_raw_plan.representation == :orthogonal_raw_product_box ||
        throw(ArgumentError("three-unit source-box shadow requires a raw product-box left PQS plan"))
    right_raw_plan.representation == :orthogonal_raw_product_box ||
        throw(ArgumentError("three-unit source-box shadow requires a raw product-box right PQS plan"))
    selected_terms = _pqs_pqs_product_supported_safe_terms(terms)
    product_unit.kind == :product_doside || throw(
        ArgumentError("three-unit source-box shadow requires a product_doside unit"),
    )
    _require_product_doside_retained_block_unit(product_unit; side = :right)

    left_count = left_raw_plan.boundary_selector.selected_count
    right_count = right_raw_plan.boundary_selector.selected_count
    product_count = length(product_unit.column_range)
    left_range = 1:left_count
    right_range = (last(left_range) + 1):(last(left_range) + right_count)
    product_range = (last(right_range) + 1):(last(right_range) + product_count)
    retained_dimension = last(product_range)
    ranges = (
        pqs_left = left_range,
        pqs_right = right_range,
        product = product_range,
    )

    pqs_left_left_references = _pqs_pqs_source_box_reference_blocks(
        left_raw_plan,
        left_raw_plan,
        metrics;
        terms = selected_terms,
    )
    pqs_left_right_references = _pqs_pqs_source_box_reference_blocks(
        left_raw_plan,
        right_raw_plan,
        metrics;
        terms = selected_terms,
    )
    pqs_right_right_references = _pqs_pqs_source_box_reference_blocks(
        right_raw_plan,
        right_raw_plan,
        metrics;
        terms = selected_terms,
    )
    pqs_left_product_references = _pqs_product_source_box_reference_blocks(
        left_raw_plan,
        product_unit,
        metrics;
        terms = selected_terms,
    )
    pqs_right_product_references = _pqs_product_source_box_reference_blocks(
        right_raw_plan,
        product_unit,
        metrics;
        terms = selected_terms,
    )
    all_pairs_inventory = _pqs_pqs_product_source_box_all_pairs_inventory(
        left_raw_plan,
        right_raw_plan,
        product_unit,
        ranges,
        selected_terms,
    )

    blocks = Dict{Symbol,Matrix{Float64}}()
    component_blocks = Dict{Symbol,NamedTuple}()
    for term in selected_terms
        pqs_left_left = pqs_left_left_references.blocks[term]
        pqs_left_right = pqs_left_right_references.blocks[term]
        pqs_right_right = pqs_right_right_references.blocks[term]
        pqs_left_product = pqs_left_product_references.blocks[term]
        pqs_right_product = pqs_right_product_references.blocks[term]
        product_product =
            _pqs_product_source_box_product_block(product_unit, metrics, term)

        block = zeros(Float64, retained_dimension, retained_dimension)
        block[left_range, left_range] .= pqs_left_left
        block[left_range, right_range] .= pqs_left_right
        block[right_range, left_range] .= transpose(pqs_left_right)
        block[left_range, product_range] .= pqs_left_product
        block[product_range, left_range] .= transpose(pqs_left_product)
        block[right_range, right_range] .= pqs_right_right
        block[right_range, product_range] .= pqs_right_product
        block[product_range, right_range] .= transpose(pqs_right_product)
        block[product_range, product_range] .= product_product
        all(isfinite, block) || throw(
            ArgumentError("three-unit source-box shadow produced non-finite entries"),
        )
        blocks[term] = block
        component_blocks[term] = (
            pqs_left_pqs_left = pqs_left_left,
            pqs_left_pqs_right = pqs_left_right,
            pqs_right_pqs_left = transpose(pqs_left_right),
            pqs_left_product = pqs_left_product,
            product_pqs_left = transpose(pqs_left_product),
            pqs_right_pqs_right = pqs_right_right,
            pqs_right_product = pqs_right_product,
            product_pqs_right = transpose(pqs_right_product),
            product_product = product_product,
        )
    end

    return (
        path = :pqs_pqs_product_source_box_shadow_blocks,
        blocks = blocks,
        component_blocks = component_blocks,
        terms = selected_terms,
        ranges = ranges,
        retained_dimension = retained_dimension,
        all_pairs_inventory = all_pairs_inventory,
        pqs_left_left_reference_blocks = pqs_left_left_references,
        pqs_left_right_reference_blocks = pqs_left_right_references,
        pqs_right_right_reference_blocks = pqs_right_right_references,
        pqs_left_product_reference_blocks = pqs_left_product_references,
        pqs_right_product_reference_blocks = pqs_right_product_references,
        diagnostics = (
            source = :pqs_pqs_product_source_box_shadow_blocks,
            source_box_shadow_only = true,
            private_shadow_only = true,
            all_pairs_inventory_private = true,
            pair_inventory_complete_for_units = (:pqs_left, :pqs_right, :product),
            all_pairs_inventory_pair_count =
                length(all_pairs_inventory.pair_entries),
            retained_unit_count = length(all_pairs_inventory.retained_units),
            packet_adoption = false,
            fixed_block_routing = false,
            pqs_representation = :mode_selected_raw_product_box,
            pqs_left_source_mode_ordering = left_raw_plan.source_mode_ordering,
            pqs_right_source_mode_ordering = right_raw_plan.source_mode_ordering,
            product_doside_retained_transform_used = true,
            pqs_self_block_source = :pqs_pqs_source_box_reference_blocks,
            cross_pqs_block_source = :pqs_pqs_source_box_reference_blocks,
            pqs_product_block_source = :pqs_product_source_box_reference_blocks,
            product_product_block_source = :product_doside_retained_block_helpers,
            pqs_left_right_raw_box_self_reference_compared =
                pqs_left_right_references.diagnostics.raw_box_self_reference_compared,
            pqs_left_left_raw_box_self_reference_compared =
                pqs_left_left_references.diagnostics.raw_box_self_reference_compared,
            pqs_right_right_raw_box_self_reference_compared =
                pqs_right_right_references.diagnostics.raw_box_self_reference_compared,
            pair_plan_reused_for_terms = true,
            lower_triangular_cross_blocks_transpose_only = true,
            explicit_source_box_oracle_tested = false,
            shell_projection_used = false,
            lowdin_cleanup_used = false,
            support_coefficient_matrix_used = false,
            support_local_pqs_oracle_used = false,
            retained_weight_semantics = :not_positive_quadrature_weights,
            retained_pqs_weights_used = false,
            retained_pqs_weights_positive_checked = false,
            ida_weight_division_allowed = false,
            qwhamiltonian_consumes = false,
            public_default_consumes = false,
            cr2_science_status_changed = false,
            local_ecp_gaussian_mwg_implemented = false,
            generic_retained_unit_framework = false,
            supported_terms = _PQS_PRODUCT_SOURCE_BOX_SHADOW_TERMS,
            dense_raw_source_box_pair_matrix_materialized = false,
            dense_raw_pair_storage_avoided = true,
            retained_block_assembled_directly_from_1d_factors = true,
            source_box_pair_storage_scaling =
                :one_dimensional_factors_plus_retained_block,
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

function _pqs_pqs_product_route_shaped_safe_term_consumer(
    route_units,
    metrics::NamedTuple{(:x,:y,:z)};
    terms = _PQS_PRODUCT_SOURCE_BOX_SHADOW_TERMS,
)
    hasproperty(route_units, :route_kind) || throw(
        ArgumentError("route-shaped safe-term consumer requires route_kind"),
    )
    hasproperty(route_units, :units) || throw(
        ArgumentError("route-shaped safe-term consumer requires units"),
    )
    route_kind = route_units.route_kind
    route_kind in _PQS_PQS_PRODUCT_SAFE_TERM_ROUTE_KINDS || throw(
        ArgumentError("unsupported route-shaped safe-term consumer route_kind $(route_kind)"),
    )
    roles = hasproperty(route_units, :roles) ?
        Tuple(route_units.roles) :
        (:pqs_left, :pqs_right, :product)
    roles == (:pqs_left, :pqs_right, :product) || throw(
        ArgumentError("route-shaped safe-term consumer requires roles (:pqs_left, :pqs_right, :product)"),
    )
    units = route_units.units
    for role in roles
        hasproperty(units, role) || throw(
            ArgumentError("route-shaped safe-term consumer missing unit $(role)"),
        )
    end

    timed = @timed begin
        _pqs_pqs_product_source_box_shadow_blocks(
            units.pqs_left,
            units.pqs_right,
            units.product,
            metrics;
            terms,
        )
    end
    shadow = timed.value
    pair_count = length(shadow.all_pairs_inventory.pair_entries)
    term_count = length(shadow.terms)
    performance = (
        elapsed_seconds = Float64(timed.time),
        allocated_bytes = Int(timed.bytes),
        gc_time_seconds = Float64(timed.gctime),
        retained_dimension = shadow.retained_dimension,
        pair_count = pair_count,
        term_count = term_count,
        dense_raw_source_box_pair_matrix_materialized = false,
        dense_raw_pair_storage_avoided = true,
    )
    component_block_provenance = (
        pqs_left_pqs_left = :pqs_pqs_source_box_reference_blocks,
        pqs_left_pqs_right = :pqs_pqs_source_box_reference_blocks,
        pqs_right_pqs_left = :transpose_of_pqs_left_pqs_right,
        pqs_left_product = :pqs_product_source_box_reference_blocks,
        product_pqs_left = :transpose_of_pqs_left_product,
        pqs_right_pqs_right = :pqs_pqs_source_box_reference_blocks,
        pqs_right_product = :pqs_product_source_box_reference_blocks,
        product_pqs_right = :transpose_of_pqs_right_product,
        product_product = :product_doside_retained_block_helpers,
    )
    metadata = hasproperty(route_units, :metadata) ? route_units.metadata : (;)
    provenance = hasproperty(route_units, :provenance) ? route_units.provenance : (;)
    route_name = hasproperty(route_units, :route_name) ? route_units.route_name : route_kind
    route_descriptor_object_kind =
        hasproperty(route_units, :object_kind) ? route_units.object_kind : :legacy_route_units
    descriptor_expected_ranges_checked = false
    descriptor_retained_dimension_checked = false
    descriptor_pair_count_checked = false
    descriptor_supported_terms_checked = false
    if hasproperty(route_units, :expected_ranges)
        route_units.expected_ranges == shadow.ranges || throw(
            ArgumentError("route descriptor expected ranges disagree with source-box shadow ranges"),
        )
        descriptor_expected_ranges_checked = true
    end
    if hasproperty(route_units, :retained_dimension)
        route_units.retained_dimension == shadow.retained_dimension || throw(
            DimensionMismatch("route descriptor retained dimension disagrees with source-box shadow"),
        )
        descriptor_retained_dimension_checked = true
    end
    if hasproperty(route_units, :expected_pair_count)
        route_units.expected_pair_count == pair_count || throw(
            ArgumentError("route descriptor expected pair count disagrees with source-box shadow"),
        )
        descriptor_pair_count_checked = true
    end
    if hasproperty(route_units, :supported_terms)
        advertised_terms = Tuple(Symbol(term) for term in route_units.supported_terms)
        for term in shadow.terms
            term in advertised_terms || throw(
                ArgumentError("route descriptor did not advertise requested term $(term)"),
            )
        end
        descriptor_supported_terms_checked = true
    end

    return (
        path = :pqs_pqs_product_route_shaped_safe_term_consumer,
        route_kind = route_kind,
        route_name = route_name,
        route_units = route_units,
        retained_units = shadow.all_pairs_inventory.retained_units,
        all_pairs_inventory = shadow.all_pairs_inventory,
        blocks = shadow.blocks,
        safe_term_matrices = shadow.blocks,
        complete_retained_space_matrices = shadow.blocks,
        component_blocks = shadow.component_blocks,
        component_block_provenance = component_block_provenance,
        pair_references = (
            pqs_left_left = shadow.pqs_left_left_reference_blocks,
            pqs_left_right = shadow.pqs_left_right_reference_blocks,
            pqs_right_right = shadow.pqs_right_right_reference_blocks,
            pqs_left_product = shadow.pqs_left_product_reference_blocks,
            pqs_right_product = shadow.pqs_right_product_reference_blocks,
        ),
        ranges = shadow.ranges,
        terms = shadow.terms,
        retained_dimension = shadow.retained_dimension,
        pair_count = pair_count,
        term_count = term_count,
        performance = performance,
        metadata = metadata,
        provenance = provenance,
        shadow = shadow,
        diagnostics = merge(
            shadow.diagnostics,
            (
                source = :pqs_pqs_product_route_shaped_safe_term_consumer,
                route_shaped_consumer = true,
                route_kind = route_kind,
                route_name = route_name,
                route_descriptor_object_kind = route_descriptor_object_kind,
                route_roles = roles,
                descriptor_expected_ranges_checked =
                    descriptor_expected_ranges_checked,
                descriptor_retained_dimension_checked =
                    descriptor_retained_dimension_checked,
                descriptor_pair_count_checked = descriptor_pair_count_checked,
                descriptor_supported_terms_checked =
                    descriptor_supported_terms_checked,
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
                performance_recorded = true,
                elapsed_seconds = performance.elapsed_seconds,
                allocated_bytes = performance.allocated_bytes,
                gc_time_seconds = performance.gc_time_seconds,
                retained_dimension = shadow.retained_dimension,
                retained_unit_count = length(shadow.all_pairs_inventory.retained_units),
                pair_count = pair_count,
                term_count = term_count,
                complete_retained_space_matrices_built = true,
                source_box_shadow_helper =
                    :_pqs_pqs_product_source_box_shadow_blocks,
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

function _product_doside_retained_kinetic_shadow_matrix(
    units,
    axis_ops::NamedTuple{(:x,:y,:z)};
    final_dimension = nothing,
)
    product_units = [unit for unit in units if unit.kind == :product_doside]
    !isempty(product_units) || throw(
        ArgumentError("product-only kinetic shadow matrix requires at least one product_doside unit"),
    )
    resolved_final_dimension = if isnothing(final_dimension)
        maximum(last(unit.column_range) for unit in units)
    else
        Int(final_dimension)
    end
    max_product_column = maximum(last(unit.column_range) for unit in product_units)
    max_product_column <= resolved_final_dimension || throw(
        ArgumentError("product-only kinetic shadow matrix final dimension does not cover product columns"),
    )
    kinetic = zeros(Float64, resolved_final_dimension, resolved_final_dimension)
    product_block_count = 0
    for right_index in eachindex(product_units)
        right_unit = product_units[right_index]
        for left_index in 1:right_index
            left_unit = product_units[left_index]
            block = _product_doside_retained_kinetic_block(
                left_unit,
                right_unit,
                axis_ops,
            )
            kinetic[left_unit.column_range, right_unit.column_range] .= block
            if left_index != right_index
                kinetic[right_unit.column_range, left_unit.column_range] .= transpose(block)
            end
            product_block_count += 1
        end
    end
    return (
        kinetic = kinetic,
        product_units = product_units,
        diagnostics = (
            source = :product_doside_retained_kinetic_shadow_matrix,
            product_only_shadow = true,
            full_kinetic_packet = false,
            support_dense_blocks_absent = true,
            mixed_blocks_absent = true,
            default_execution_changed = false,
            qwhamiltonian_consumes = false,
            backend_policy_changed = false,
            quadrature_policy_changed = false,
            cr2_science_status_changed = false,
            product_unit_count = length(product_units),
            product_block_count = product_block_count,
            final_dimension = resolved_final_dimension,
        ),
    )
end

function _staged_retained_kinetic_shadow_matrix(
    units,
    axis_ops::NamedTuple{(:x,:y,:z)};
    final_dimension = nothing,
)
    !isempty(units) || throw(
        ArgumentError("staged retained kinetic shadow matrix requires staged units"),
    )
    resolved_final_dimension = if isnothing(final_dimension)
        maximum(last(unit.column_range) for unit in units)
    else
        Int(final_dimension)
    end
    for unit in units
        first(unit.column_range) >= 1 &&
            last(unit.column_range) <= resolved_final_dimension ||
            throw(
                ArgumentError("staged retained kinetic shadow matrix final dimension does not cover unit columns"),
            )
    end
    kinetic = zeros(Float64, resolved_final_dimension, resolved_final_dimension)
    product_unit_count = count(unit -> unit.kind == :product_doside, units)
    generic_unit_count = length(units) - product_unit_count
    product_block_count = 0
    fallback_block_count = 0
    entries_cache = Vector{Union{Nothing,Vector{Vector{_ParentCoefficientEntry3D}}}}(
        undef,
        length(units),
    )
    fill!(entries_cache, nothing)
    for right_index in eachindex(units)
        right_unit = units[right_index]
        for left_index in 1:right_index
            left_unit = units[left_index]
            if left_unit.kind == :product_doside && right_unit.kind == :product_doside
                block = _product_doside_retained_kinetic_block(
                    left_unit,
                    right_unit,
                    axis_ops,
                )
                product_block_count += 1
            else
                left_entries = _cached_staged_unit_entries!(
                    entries_cache,
                    left_index,
                    left_unit,
                )
                right_entries = _cached_staged_unit_entries!(
                    entries_cache,
                    right_index,
                    right_unit,
                )
                block = _fallback_staged_kinetic_block(
                    left_entries,
                    right_entries,
                    axis_ops,
                )
                fallback_block_count += 1
            end
            kinetic[left_unit.column_range, right_unit.column_range] .= block
            if left_index != right_index
                kinetic[right_unit.column_range, left_unit.column_range] .= transpose(block)
            end
        end
    end
    return (
        kinetic = kinetic,
        diagnostics = (
            source = :staged_retained_kinetic_shadow_matrix,
            full_private_shadow_matrix = true,
            product_only_shadow = false,
            production_adoption = false,
            default_execution_changed = false,
            metric_packet_execution_changed = false,
            fixed_block_construction_changed = false,
            qwhamiltonian_consumes = false,
            backend_policy_changed = false,
            quadrature_policy_changed = false,
            ida_positive_weight_semantics_changed = false,
            cr2_science_status_changed = false,
            product_unit_count = product_unit_count,
            generic_unit_count = generic_unit_count,
            product_block_count = product_block_count,
            fallback_block_count = fallback_block_count,
            total_block_count = product_block_count + fallback_block_count,
            final_dimension = resolved_final_dimension,
        ),
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

const _PACKET_BUILD_SOURCE_SAFE_FIELDS = (
    :overlap,
    :position_x,
    :position_y,
    :position_z,
    :x2_x,
    :x2_y,
    :x2_z,
    :weights,
    :first_moments,
    :kinetic,
)

function _packet_build_source_safe_fields(fields)
    selected = Tuple(Symbol(field) for field in fields)
    isempty(selected) && throw(
        ArgumentError("packet-build source safe-field shadow requires at least one field"),
    )
    unsupported = setdiff(selected, _PACKET_BUILD_SOURCE_SAFE_FIELDS)
    isempty(unsupported) || throw(
        ArgumentError("packet-build source safe-field shadow does not support fields $(unsupported)"),
    )
    length(unique(selected)) == length(selected) || throw(
        ArgumentError("packet-build source safe-field shadow fields must be unique"),
    )
    return selected
end

function _packet_build_axis_property(axis_data, axis_name::Symbol, property_name::Symbol)
    hasproperty(axis_data, property_name) || throw(
        ArgumentError("packet-build source safe-field shadow axis $(axis_name) is missing $(property_name)"),
    )
    return getproperty(axis_data, property_name)
end

function _packet_build_axis_data(axis_data, axis_name::Symbol)
    overlap = Matrix{Float64}(
        _packet_build_axis_property(axis_data, axis_name, :overlap),
    )
    position = Matrix{Float64}(
        _packet_build_axis_property(axis_data, axis_name, :position),
    )
    x2 = Matrix{Float64}(
        _packet_build_axis_property(axis_data, axis_name, :x2),
    )
    kinetic = Matrix{Float64}(
        _packet_build_axis_property(axis_data, axis_name, :kinetic),
    )
    weights = Float64[
        Float64(value) for value in
        _packet_build_axis_property(axis_data, axis_name, :weights)
    ]
    centers = Float64[
        Float64(value) for value in
        _packet_build_axis_property(axis_data, axis_name, :centers)
    ]
    size(overlap, 1) == size(overlap, 2) || throw(
        ArgumentError("packet-build source safe-field shadow $(axis_name) overlap must be square"),
    )
    size(position) == size(overlap) || throw(
        ArgumentError("packet-build source safe-field shadow $(axis_name) position size must match overlap"),
    )
    size(x2) == size(overlap) || throw(
        ArgumentError("packet-build source safe-field shadow $(axis_name) x2 size must match overlap"),
    )
    size(kinetic) == size(overlap) || throw(
        ArgumentError("packet-build source safe-field shadow $(axis_name) kinetic size must match overlap"),
    )
    length(weights) == size(overlap, 1) || throw(
        ArgumentError("packet-build source safe-field shadow $(axis_name) weights length must match overlap"),
    )
    length(centers) == size(overlap, 1) || throw(
        ArgumentError("packet-build source safe-field shadow $(axis_name) centers length must match overlap"),
    )
    all(isfinite, overlap) &&
        all(isfinite, position) &&
        all(isfinite, x2) &&
        all(isfinite, kinetic) &&
        all(isfinite, weights) &&
        all(isfinite, centers) || throw(
            ArgumentError("packet-build source safe-field shadow $(axis_name) axis data must be finite"),
        )
    source = hasproperty(axis_data, :source) ?
        Symbol(getproperty(axis_data, :source)) :
        :explicit_packet_build_axis_data
    return (
        overlap = overlap,
        position = position,
        x2 = x2,
        weights = weights,
        centers = centers,
        kinetic = kinetic,
        source = source,
    )
end

function _packet_build_safe_axis_data(axis_data)
    for axis_name in (:x, :y, :z)
        hasproperty(axis_data, axis_name) || throw(
            ArgumentError("packet-build source safe-field shadow axis data is missing $(axis_name)"),
        )
    end
    return (
        x = _packet_build_axis_data(getproperty(axis_data, :x), :x),
        y = _packet_build_axis_data(getproperty(axis_data, :y), :y),
        z = _packet_build_axis_data(getproperty(axis_data, :z), :z),
    )
end

function _packet_build_safe_axis_ops(axis_data::NamedTuple{(:x,:y,:z)})
    return (
        x = (overlap = axis_data.x.overlap, kinetic = axis_data.x.kinetic),
        y = (overlap = axis_data.y.overlap, kinetic = axis_data.y.kinetic),
        z = (overlap = axis_data.z.overlap, kinetic = axis_data.z.kinetic),
    )
end

function _packet_build_assign_pair_block!(
    matrix::Matrix{Float64},
    left_range::UnitRange{Int},
    right_range::UnitRange{Int},
    block::AbstractMatrix{<:Real},
)
    matrix[left_range, right_range] .= block
    left_range == right_range || (matrix[right_range, left_range] .= transpose(block))
    return nothing
end

function _packet_build_source_units(source)
    payloads = collect(source.resolved_payloads)
    isempty(payloads) && throw(
        ArgumentError("packet-build source safe-field shadow requires at least one payload"),
    )
    for payload in payloads
        payload.ready_for_metric_execution || throw(
            ArgumentError("packet-build source safe-field shadow cannot consume non-executable payload $(payload.payload_kind)"),
        )
        payload.payload_kind in (:product_doside, :support_dense) || throw(
            ArgumentError("packet-build source safe-field shadow supports only product_doside/support_dense payloads; got $(payload.payload_kind)"),
        )
        isnothing(payload.column_range) && throw(
            ArgumentError("packet-build source safe-field shadow requires every payload to have a column range"),
        )
        isnothing(payload.support_indices) && throw(
            ArgumentError("packet-build source safe-field shadow requires every payload to have support indices"),
        )
    end
    return [payload.payload for payload in payloads]
end

function _cartesian_packet_build_source_safe_field_shadow(
    source,
    axis_data;
    fields = _PACKET_BUILD_SOURCE_SAFE_FIELDS,
)
    selected_fields = _packet_build_source_safe_fields(fields)
    for field in selected_fields
        field in source.candidate_packet_fields || throw(
            ArgumentError("packet-build source safe-field shadow field $(field) is not a candidate packet field"),
        )
        field in source.missing_packet_fields && throw(
            ArgumentError("packet-build source safe-field shadow field $(field) is marked missing/not implemented"),
        )
    end
    source.diagnostics.column_layout_ready_for_candidate_fields || throw(
        ArgumentError("packet-build source safe-field shadow requires complete non-overlapping column layout"),
    )

    units = _packet_build_source_units(source)
    axis_data = _packet_build_safe_axis_data(axis_data)
    axis_ops = _packet_build_safe_axis_ops(axis_data)
    contracted_dimension = source.contracted_dimension

    build_low_order = any(
        field -> field in selected_fields,
        (:overlap, :position_x, :position_y, :position_z),
    )
    build_x2 = any(field -> field in selected_fields, (:x2_x, :x2_y, :x2_z))
    build_linear = any(field -> field in selected_fields, (:weights, :first_moments))
    build_kinetic = :kinetic in selected_fields

    overlap = build_low_order ? zeros(Float64, contracted_dimension, contracted_dimension) : nothing
    position_x = build_low_order ? zeros(Float64, contracted_dimension, contracted_dimension) : nothing
    position_y = build_low_order ? zeros(Float64, contracted_dimension, contracted_dimension) : nothing
    position_z = build_low_order ? zeros(Float64, contracted_dimension, contracted_dimension) : nothing
    x2_x = build_x2 ? zeros(Float64, contracted_dimension, contracted_dimension) : nothing
    x2_y = build_x2 ? zeros(Float64, contracted_dimension, contracted_dimension) : nothing
    x2_z = build_x2 ? zeros(Float64, contracted_dimension, contracted_dimension) : nothing
    weights = build_linear ? zeros(Float64, contracted_dimension) : nothing
    first_moments = build_linear ? zeros(Float64, contracted_dimension, 3) : nothing
    kinetic = build_kinetic ? zeros(Float64, contracted_dimension, contracted_dimension) : nothing

    entries_cache = Vector{Union{Nothing,Vector{Vector{_ParentCoefficientEntry3D}}}}(
        undef,
        length(units),
    )
    fill!(entries_cache, nothing)

    low_order_product_block_count = 0
    low_order_fallback_block_count = 0
    x2_product_block_count = 0
    x2_fallback_block_count = 0
    kinetic_product_block_count = 0
    kinetic_fallback_block_count = 0
    total_block_count = length(units) * (length(units) + 1) ÷ 2
    for right_index in eachindex(units)
        right = units[right_index]
        for left_index in 1:right_index
            left = units[left_index]
            left_range = left.column_range
            right_range = right.column_range
            if build_low_order
                if left.kind == :product_doside && right.kind == :product_doside
                    _fill_product_staged_metric_blocks!(
                        overlap,
                        position_x,
                        position_y,
                        position_z,
                        left,
                        right,
                        axis_data,
                    )
                    low_order_product_block_count += 1
                else
                    left_entries = _cached_staged_unit_entries!(entries_cache, left_index, left)
                    right_entries = _cached_staged_unit_entries!(entries_cache, right_index, right)
                    blocks = _fallback_staged_metric_blocks(left_entries, right_entries, axis_data)
                    _packet_build_assign_pair_block!(overlap, left_range, right_range, blocks.overlap)
                    _packet_build_assign_pair_block!(
                        position_x,
                        left_range,
                        right_range,
                        blocks.position_x,
                    )
                    _packet_build_assign_pair_block!(
                        position_y,
                        left_range,
                        right_range,
                        blocks.position_y,
                    )
                    _packet_build_assign_pair_block!(
                        position_z,
                        left_range,
                        right_range,
                        blocks.position_z,
                    )
                    low_order_fallback_block_count += 1
                end
            end
            if build_x2
                if left.kind == :product_doside && right.kind == :product_doside
                    x2_x_block = _product_doside_retained_low_order_block(
                        left,
                        right,
                        axis_data;
                        term = :x2_x,
                    )
                    x2_y_block = _product_doside_retained_low_order_block(
                        left,
                        right,
                        axis_data;
                        term = :x2_y,
                    )
                    x2_z_block = _product_doside_retained_low_order_block(
                        left,
                        right,
                        axis_data;
                        term = :x2_z,
                    )
                    x2_product_block_count += 1
                else
                    left_entries = _cached_staged_unit_entries!(entries_cache, left_index, left)
                    right_entries = _cached_staged_unit_entries!(entries_cache, right_index, right)
                    blocks = _fallback_staged_x2_blocks(left_entries, right_entries, axis_data)
                    x2_x_block = blocks.x2_x
                    x2_y_block = blocks.x2_y
                    x2_z_block = blocks.x2_z
                    x2_fallback_block_count += 1
                end
                _packet_build_assign_pair_block!(x2_x, left_range, right_range, x2_x_block)
                _packet_build_assign_pair_block!(x2_y, left_range, right_range, x2_y_block)
                _packet_build_assign_pair_block!(x2_z, left_range, right_range, x2_z_block)
            end
            if build_kinetic
                if left.kind == :product_doside && right.kind == :product_doside
                    block = _product_doside_retained_kinetic_block(left, right, axis_ops)
                    kinetic_product_block_count += 1
                else
                    left_entries = _cached_staged_unit_entries!(entries_cache, left_index, left)
                    right_entries = _cached_staged_unit_entries!(entries_cache, right_index, right)
                    block = _fallback_staged_kinetic_block(left_entries, right_entries, axis_ops)
                    kinetic_fallback_block_count += 1
                end
                _packet_build_assign_pair_block!(kinetic, left_range, right_range, block)
            end
        end
    end

    if build_linear
        for unit in units
            unit_weights, unit_first_moments = _staged_unit_linear_vectors(unit, axis_data)
            weights[unit.column_range] .= unit_weights
            first_moments[unit.column_range, :] .= unit_first_moments
        end
    end

    product_unit_count = count(unit -> unit.kind == :product_doside, units)
    support_dense_unit_count = count(unit -> unit.kind == :support_dense, units)
    return (
        overlap = :overlap in selected_fields ? overlap : nothing,
        position_x = :position_x in selected_fields ? position_x : nothing,
        position_y = :position_y in selected_fields ? position_y : nothing,
        position_z = :position_z in selected_fields ? position_z : nothing,
        x2_x = :x2_x in selected_fields ? x2_x : nothing,
        x2_y = :x2_y in selected_fields ? x2_y : nothing,
        x2_z = :x2_z in selected_fields ? x2_z : nothing,
        weights = :weights in selected_fields ? weights : nothing,
        first_moments = :first_moments in selected_fields ? first_moments : nothing,
        kinetic = :kinetic in selected_fields ? kinetic : nothing,
        diagnostics = (
            source = :cartesian_packet_build_source_safe_field_shadow,
            source_driven_shadow_only = true,
            construction_adoption = false,
            current_builder_authority = :nested_shell_packet,
            nested_shell_packet_remains_authoritative = true,
            descriptor_drives_builder = false,
            numerical_packet_authority_changed = false,
            numerical_packet_matrices_built = true,
            operator_data_available_on_source = source.diagnostics.operator_data_available,
            packet_operator_data_checked_on_source = source.diagnostics.packet_operator_data_checked,
            fixed_block_construction_changed = false,
            metric_packet_execution_changed = false,
            qwhamiltonian_changed = false,
            backend_policy_changed = false,
            quadrature_policy_changed = false,
            ida_positive_weight_semantics_changed = false,
            cr2_science_status_changed = false,
            requested_fields = selected_fields,
            candidate_packet_fields = source.candidate_packet_fields,
            missing_packet_fields = source.missing_packet_fields,
            out_of_scope_fields = source.missing_packet_fields,
            x2_built = build_x2,
            x2_product_block_count = x2_product_block_count,
            x2_fallback_block_count = x2_fallback_block_count,
            gaussian_terms_built = false,
            nuclear_one_body_built = false,
            local_coulomb_one_body_built = false,
            local_ecp_one_body_built = false,
            pair_mwg_interaction_built = false,
            parent_dimension = source.parent_dimension,
            contracted_dimension = contracted_dimension,
            unit_count = length(units),
            product_unit_count = product_unit_count,
            support_dense_unit_count = support_dense_unit_count,
            total_upper_triangular_block_count = total_block_count,
            low_order_product_block_count = low_order_product_block_count,
            low_order_fallback_block_count = low_order_fallback_block_count,
            kinetic_product_block_count = kinetic_product_block_count,
            kinetic_fallback_block_count = kinetic_fallback_block_count,
            first_moment_reference_required = :contracted_parent_metric_packet_or_support_local_reference,
        ),
    )
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

module CartesianContractedParentMetrics

import LinearAlgebra
import SparseArrays

import ..GaussletBases: CoulombGaussianExpansion,
                         _NestedFixedBlock3D,
                         _BondAlignedDiatomicHighOrderRecipeSourceConstruction3D,
                         _CartesianNestedAxisBundles3D,
                         _CartesianNestedProjectedQShellStagedUnitDescriptor3D,
                         _CartesianNestedProductStagedByCenterUnit3D,
                         _MappedOrdinaryGausslet1DBundle,
                         _MappedOrdinaryPGDGIntermediate1D,
                         _cartesian_flat_index,
                         _cartesian_raw_product_box_plan,
                         _cartesian_raw_product_box_source_mode_indices,
                         _cartesian_unflat_index,
                         _nested_axis_bundle,
                         _nested_axis_lengths,
                         _nested_product_axis_function_indices,
                         _nested_product_staged_active_axis,
                         _nested_product_staged_fixed_axis,
                         _nested_bond_aligned_diatomic_high_order_recipe_source_fixed_block,
                         _nested_projected_q_shell_staged_unit_descriptor,
                         _nested_projected_q_shell_descriptor_seed_coefficients,
                         _require_analytic_primitive_backend,
                         centers,
                         contract_primitive_matrix,
                         gaussian_factor_matrices,
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
        object_contract = :SourceBoxPairOperatorPlan,
        pair_policy = :source_box_algorithm_available,
        left_source_family = :mode_selected_raw_product_box,
        right_source_family = :mode_selected_raw_product_box,
        left_raw_product_box_plan_contract = :RawProductBoxPlan,
        right_raw_product_box_plan_contract = :RawProductBoxPlan,
        left_retained_rule_contract = :RetainedRule,
        right_retained_rule_contract = :RetainedRule,
        left_retained_rule_kind = :boundary_comx_product_mode_selection,
        right_retained_rule_kind = :boundary_comx_product_mode_selection,
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
            source_box_pair_operator_plan_contract =
                :SourceBoxPairOperatorPlan,
            pair_policy = :source_box_algorithm_available,
            algorithmic_pair_policy = :source_box_algorithm_available,
            left_raw_product_box_plan_contract = :RawProductBoxPlan,
            right_raw_product_box_plan_contract = :RawProductBoxPlan,
            left_retained_rule_contract = :RetainedRule,
            right_retained_rule_contract = :RetainedRule,
            left_retained_rule_kind =
                :boundary_comx_product_mode_selection,
            right_retained_rule_kind =
                :boundary_comx_product_mode_selection,
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
            shell_realization_adapter_used = false,
            support_row_adapter_used = false,
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

function _pqs_pqs_source_box_explicit_boundary_selection_reference(
    pair_plan,
    term::Symbol,
)
    term in pair_plan.supported_terms || throw(
        ArgumentError("PQS/PQS source-box explicit reference received unsupported term $(term)"),
    )
    result = zeros(
        Float64,
        pair_plan.left_retained_count,
        pair_plan.right_retained_count,
    )
    left_columns = pair_plan.left_boundary_mode_selector.column_indices
    right_columns = pair_plan.right_boundary_mode_selector.column_indices
    @inbounds for factor_kinds in _source_box_separable_term_factor_kinds(term)
        axis_matrices = ntuple(
            axis -> _pqs_pqs_source_box_factor(pair_plan, axis, factor_kinds[axis]),
            3,
        )
        raw_pair_matrix = _pqs_raw_product_box_pair_mode_matrix(
            pair_plan.left_raw_product_box_plan.source_mode_indices,
            pair_plan.right_raw_product_box_plan.source_mode_indices,
            axis_matrices,
        )
        result .+= raw_pair_matrix[left_columns, right_columns]
    end
    return result
end

function _pqs_pqs_source_box_reference_blocks_from_pair_plan(
    pair_plan;
    terms = pair_plan.supported_terms,
    atol::Real = 1.0e-10,
    validate_explicit_raw_box_oracle::Bool = true,
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
        if validate_explicit_raw_box_oracle
            raw_reference = _pqs_pqs_source_box_explicit_boundary_selection_reference(
                pair_plan,
                term,
            )
            block_error = LinearAlgebra.norm(block - raw_reference, Inf)
            block_error <= atol || throw(
                ArgumentError("PQS/PQS source-box reference block disagrees with explicit raw product-box boundary-selection reference"),
            )
            raw_box_reference_blocks[term] = raw_reference
            block_errors[term] = block_error
        elseif same_raw_product_box_plan
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
                source_box_algorithm_formula_available = true,
                validation_reference_contract =
                    validate_explicit_raw_box_oracle ?
                    :explicit_raw_product_box_boundary_column_selection :
                    (
                        same_raw_product_box_plan ?
                        :raw_box_self_reference_helper :
                        :external_raw_product_box_boundary_column_selection_required
                    ),
                internal_validation_reference_compared =
                    validate_explicit_raw_box_oracle || same_raw_product_box_plan,
                explicit_raw_product_box_boundary_column_selection_reference_compared =
                    validate_explicit_raw_box_oracle,
                explicit_raw_product_box_boundary_column_selection_reference_helper =
                    validate_explicit_raw_box_oracle ?
                    :_pqs_pqs_source_box_explicit_boundary_selection_reference :
                    nothing,
                explicit_source_box_oracle_tested =
                    validate_explicit_raw_box_oracle,
                cross_box_external_raw_product_oracle_required =
                    !same_raw_product_box_plan && !validate_explicit_raw_box_oracle,
                cross_box_external_raw_product_oracle_compared_by_helper =
                    !same_raw_product_box_plan && validate_explicit_raw_box_oracle,
                dense_raw_source_box_pair_matrix_materialized_for_validation =
                    validate_explicit_raw_box_oracle,
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
    return _product_doside_source_box_reference_block(
        product_unit,
        product_unit,
        metrics;
        term,
    ).block
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
            pair_policy = :source_box_algorithm_available,
            source_box_algorithmic = true,
        ),
        (
            pair_key = (:pqs, :product),
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
            product_product_block_source =
                :_product_doside_source_box_reference_block,
            every_pair_uses_source_box_algorithmic_policy =
                all_pairs_inventory.diagnostics.every_pair_uses_source_box_algorithmic_policy,
            source_box_algorithmic_pair_count =
                all_pairs_inventory.diagnostics.source_box_algorithmic_pair_count,
            pair_policies = all_pairs_inventory.diagnostics.pair_policies,
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

function _pqs_inventory_wrapped_staged_unit(
    unit::_CartesianNestedProductStagedByCenterUnit3D;
    category::Symbol,
    support_source_semantics::Symbol,
    safe_term_capability::Symbol,
    active_representation_stage::Symbol,
    source_helper::Symbol,
    raw_box_auxiliary_metadata = nothing,
)
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
            product_product_block_source =
                :_product_doside_source_box_reference_block,
            pqs_left_right_raw_box_self_reference_compared =
                pqs_left_right_references.diagnostics.raw_box_self_reference_compared,
            pqs_left_left_raw_box_self_reference_compared =
                pqs_left_left_references.diagnostics.raw_box_self_reference_compared,
            pqs_right_right_raw_box_self_reference_compared =
                pqs_right_right_references.diagnostics.raw_box_self_reference_compared,
            pqs_left_left_explicit_source_box_oracle_tested =
                pqs_left_left_references.diagnostics.explicit_source_box_oracle_tested,
            pqs_left_right_explicit_source_box_oracle_tested =
                pqs_left_right_references.diagnostics.explicit_source_box_oracle_tested,
            pqs_right_right_explicit_source_box_oracle_tested =
                pqs_right_right_references.diagnostics.explicit_source_box_oracle_tested,
            pair_plan_reused_for_terms = true,
            lower_triangular_cross_blocks_transpose_only = true,
            explicit_source_box_oracle_tested =
                pqs_left_left_references.diagnostics.explicit_source_box_oracle_tested &&
                pqs_left_right_references.diagnostics.explicit_source_box_oracle_tested &&
                pqs_right_right_references.diagnostics.explicit_source_box_oracle_tested,
            pqs_cross_box_external_raw_product_oracle_required =
                pqs_left_right_references.diagnostics.cross_box_external_raw_product_oracle_required,
            pqs_cross_box_internal_raw_product_oracle_compared =
                pqs_left_right_references.diagnostics.cross_box_external_raw_product_oracle_compared_by_helper,
            every_pair_uses_source_box_algorithmic_policy =
                all_pairs_inventory.diagnostics.every_pair_uses_source_box_algorithmic_policy,
            source_box_algorithmic_pair_count =
                all_pairs_inventory.diagnostics.source_box_algorithmic_pair_count,
            pair_policies = all_pairs_inventory.diagnostics.pair_policies,
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
            dense_raw_source_box_pair_matrix_materialized_for_validation =
                pqs_left_left_references.diagnostics.dense_raw_source_box_pair_matrix_materialized_for_validation ||
                pqs_left_right_references.diagnostics.dense_raw_source_box_pair_matrix_materialized_for_validation ||
                pqs_right_right_references.diagnostics.dense_raw_source_box_pair_matrix_materialized_for_validation,
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
        dense_raw_source_box_pair_matrix_materialized_for_validation =
            shadow.diagnostics.dense_raw_source_box_pair_matrix_materialized_for_validation,
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
        product_product = :_product_doside_source_box_reference_block,
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
                dense_raw_source_box_pair_matrix_materialized_for_validation =
                    shadow.diagnostics.dense_raw_source_box_pair_matrix_materialized_for_validation,
                dense_raw_pair_storage_avoided = true,
                every_pair_uses_source_box_algorithmic_policy =
                    shadow.diagnostics.every_pair_uses_source_box_algorithmic_policy,
                source_box_algorithmic_pair_count =
                    shadow.diagnostics.source_box_algorithmic_pair_count,
                pair_policies = shadow.diagnostics.pair_policies,
                pqs_cross_box_external_raw_product_oracle_required =
                    shadow.diagnostics.pqs_cross_box_external_raw_product_oracle_required,
                pqs_cross_box_internal_raw_product_oracle_compared =
                    shadow.diagnostics.pqs_cross_box_internal_raw_product_oracle_compared,
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

const _PQS_PQS_PRODUCT_DENSITY_DENSITY_ROUTE_KINDS =
    _PQS_PQS_PRODUCT_SAFE_TERM_ROUTE_KINDS

function _source_box_axis_pair_terms_symmetric(
    axis_pair_factor_terms::AbstractArray{<:Real,3};
    axis_name::Symbol,
    atol::Real = 0.0,
)
    nterms = size(axis_pair_factor_terms, 1)
    nterms > 0 || throw(
        ArgumentError("$(axis_name) pair-factor terms require at least one term"),
    )
    size(axis_pair_factor_terms, 2) == size(axis_pair_factor_terms, 3) || throw(
        DimensionMismatch("$(axis_name) pair-factor terms must be square per term"),
    )
    atol_value = Float64(atol)
    for term in axes(axis_pair_factor_terms, 1)
        term_matrix = @view axis_pair_factor_terms[term, :, :]
        all(isfinite, term_matrix) || throw(
            ArgumentError("$(axis_name) pair-factor terms must be finite"),
        )
        isapprox(
            term_matrix,
            transpose(term_matrix);
            atol = atol_value,
            rtol = atol_value,
        ) || return false
    end
    return true
end

function _source_box_axis_pair_terms_symmetric(
    axis_pair_factor_terms::NamedTuple{(:x,:y,:z)};
    atol::Real = 0.0,
)
    return all(axis -> _source_box_axis_pair_terms_symmetric(
        getproperty(axis_pair_factor_terms, (:x, :y, :z)[axis]);
        axis_name = (:x, :y, :z)[axis],
        atol,
    ), 1:3)
end

function _pqs_pqs_product_density_density_route_ranges(
    left_raw_plan,
    right_raw_plan,
    product_unit::_CartesianNestedProductStagedByCenterUnit3D,
)
    product_retained_unit_plan = _product_doside_retained_unit_plan(product_unit)
    left_count = left_raw_plan.boundary_selector.selected_count
    right_count = right_raw_plan.boundary_selector.selected_count
    product_count = product_retained_unit_plan.retained_count
    ranges = (
        pqs_left = 1:left_count,
        pqs_right = (left_count + 1):(left_count + right_count),
        product = (left_count + right_count + 1):
                  (left_count + right_count + product_count),
    )
    return (
        ranges = ranges,
        retained_dimension = left_count + right_count + product_count,
    )
end

function _pqs_pqs_product_density_density_all_pairs_inventory(
    left_raw_plan,
    right_raw_plan,
    product_unit::_CartesianNestedProductStagedByCenterUnit3D,
    ranges,
    pair_factor_normalization::Symbol,
)
    pair_factor_normalization in (:density_normalized, :raw_weighted) || throw(
        ArgumentError("density-density route requires density_normalized or raw_weighted pair factors"),
    )
    product_retained_unit_plan = _product_doside_retained_unit_plan(product_unit)
    pqs_pqs_helper =
        pair_factor_normalization == :raw_weighted ?
        :_pqs_pqs_source_box_raw_weighted_density_density_interaction_block :
        :_pqs_pqs_source_box_density_density_interaction_block
    pqs_product_helper =
        pair_factor_normalization == :raw_weighted ?
        :_pqs_product_source_box_raw_weighted_density_density_interaction_block :
        :_pqs_product_source_box_density_density_interaction_block
    product_product_helper =
        pair_factor_normalization == :raw_weighted ?
        :_product_doside_source_box_raw_weighted_density_density_interaction_block :
        :_product_doside_source_box_density_density_interaction_block
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
            supported_interaction_terms = (:pair_sum,),
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
            supported_interaction_terms = (:pair_sum,),
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
            supported_interaction_terms = (:pair_sum,),
        ),
    )
    pair_entries = (
        (
            pair_key = (:pqs_left, :pqs_left),
            pair_family = :pqs_pqs,
            pair_kind = :pqs_pqs_source_box_density_density_pair,
            block_helper = pqs_pqs_helper,
            upper_triangular = true,
            transpose_only_lower_block = false,
            pair_policy = :source_box_algorithm_available,
            source_box_algorithmic = true,
        ),
        (
            pair_key = (:pqs_left, :pqs_right),
            pair_family = :pqs_pqs,
            pair_kind = :pqs_pqs_source_box_density_density_pair,
            block_helper = pqs_pqs_helper,
            upper_triangular = true,
            transpose_only_lower_block = true,
            pair_policy = :source_box_algorithm_available,
            source_box_algorithmic = true,
        ),
        (
            pair_key = (:pqs_left, :product),
            pair_family = :pqs_product,
            pair_kind = :pqs_product_source_box_density_density_pair,
            block_helper = pqs_product_helper,
            upper_triangular = true,
            transpose_only_lower_block = true,
            pair_policy = :source_box_algorithm_available,
            source_box_algorithmic = true,
        ),
        (
            pair_key = (:pqs_right, :pqs_right),
            pair_family = :pqs_pqs,
            pair_kind = :pqs_pqs_source_box_density_density_pair,
            block_helper = pqs_pqs_helper,
            upper_triangular = true,
            transpose_only_lower_block = false,
            pair_policy = :source_box_algorithm_available,
            source_box_algorithmic = true,
        ),
        (
            pair_key = (:pqs_right, :product),
            pair_family = :pqs_product,
            pair_kind = :pqs_product_source_box_density_density_pair,
            block_helper = pqs_product_helper,
            upper_triangular = true,
            transpose_only_lower_block = true,
            pair_policy = :source_box_algorithm_available,
            source_box_algorithmic = true,
        ),
        (
            pair_key = (:product, :product),
            pair_family = :product_product,
            pair_kind = :product_doside_source_box_density_density_pair,
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
        object_kind = :pqs_pqs_product_density_density_all_pairs_inventory,
        retained_units = retained_units,
        pair_entries = pair_entries,
        pair_family_counts = pair_family_counts,
        supported_interaction_terms = (:pair_sum,),
        diagnostics = (
            source = :pqs_pqs_product_density_density_all_pairs_inventory,
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
            pair_factor_normalization = pair_factor_normalization,
            private_shadow_only = true,
            output_representation = :two_index_density_density,
            four_index_galerkin_tensor = false,
            packet_adoption = false,
            fixed_block_routing = false,
            qwhamiltonian_consumes = false,
            public_default_consumes = false,
            cr2_science_status_changed = false,
            ecp_terms_implemented = false,
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
            lower_triangular_cross_blocks_transpose_only = true,
            product_pqs_blocks_transpose_only = true,
            product_pqs_explicit_helper_used = false,
            real_mwg_ida_pair_factor_provenance_adapted = false,
        ),
    )
end

function _pqs_pqs_product_density_density_pair_block(
    entry,
    units;
    term_coefficients::AbstractVector{<:Real},
    axis_pair_factor_terms,
    raw_axis_pair_factor_terms,
    axis_weights::NamedTuple{(:x,:y,:z)},
    pair_factor_normalization::Symbol,
)
    left_role, right_role = entry.pair_key
    if entry.pair_family == :pqs_pqs
        left_unit = getproperty(units, left_role)
        right_unit = getproperty(units, right_role)
        if pair_factor_normalization == :raw_weighted
            return _pqs_pqs_source_box_raw_weighted_density_density_interaction_block(
                left_unit,
                right_unit;
                term_coefficients,
                raw_axis_pair_factor_terms,
                axis_weights,
            )
        end
        return _pqs_pqs_source_box_density_density_interaction_block(
            left_unit,
            right_unit;
            term_coefficients,
            axis_pair_factor_terms,
            axis_weights,
        )
    elseif entry.pair_family == :pqs_product
        right_role == :product || throw(
            ArgumentError("density-density route only uses product/PQS as transpose of PQS/product"),
        )
        pqs_unit = getproperty(units, left_role)
        product_unit = getproperty(units, right_role)
        if pair_factor_normalization == :raw_weighted
            return _pqs_product_source_box_raw_weighted_density_density_interaction_block(
                pqs_unit,
                product_unit;
                term_coefficients,
                raw_axis_pair_factor_terms,
                axis_weights,
            )
        end
        return _pqs_product_source_box_density_density_interaction_block(
            pqs_unit,
            product_unit;
            term_coefficients,
            axis_pair_factor_terms,
            axis_weights,
        )
    elseif entry.pair_family == :product_product
        left_role == :product && right_role == :product || throw(
            ArgumentError("density-density route product/product entry must use the product unit on both sides"),
        )
        product_unit = getproperty(units, :product)
        if pair_factor_normalization == :raw_weighted
            return _product_doside_source_box_raw_weighted_density_density_interaction_block(
                product_unit,
                product_unit;
                term_coefficients,
                raw_axis_pair_factor_terms,
                axis_weights,
            )
        end
        return _product_doside_source_box_density_density_interaction_block(
            product_unit,
            product_unit;
            term_coefficients,
            axis_pair_factor_terms,
            axis_weights,
        )
    end
    throw(ArgumentError("unsupported density-density route pair family $(entry.pair_family)"))
end

function _pqs_pqs_product_route_shaped_density_density_consumer(
    route_units;
    term_coefficients::AbstractVector{<:Real},
    axis_pair_factor_terms = nothing,
    raw_axis_pair_factor_terms = nothing,
    axis_weights::NamedTuple{(:x,:y,:z)},
    pair_factor_normalization::Symbol = :density_normalized,
    pair_factor_symmetry_atol::Real = 1.0e-12,
    symmetry_atol::Real = 1.0e-10,
)
    hasproperty(route_units, :route_kind) || throw(
        ArgumentError("route-shaped density-density consumer requires route_kind"),
    )
    hasproperty(route_units, :units) || throw(
        ArgumentError("route-shaped density-density consumer requires units"),
    )
    route_kind = route_units.route_kind
    route_kind in _PQS_PQS_PRODUCT_DENSITY_DENSITY_ROUTE_KINDS || throw(
        ArgumentError("unsupported route-shaped density-density consumer route_kind $(route_kind)"),
    )
    roles = hasproperty(route_units, :roles) ?
        Tuple(route_units.roles) :
        (:pqs_left, :pqs_right, :product)
    roles == (:pqs_left, :pqs_right, :product) || throw(
        ArgumentError("route-shaped density-density consumer requires roles (:pqs_left, :pqs_right, :product)"),
    )
    units = route_units.units
    for role in roles
        hasproperty(units, role) || throw(
            ArgumentError("route-shaped density-density consumer missing unit $(role)"),
        )
    end
    pair_factor_normalization in (:density_normalized, :raw_weighted) || throw(
        ArgumentError("route-shaped density-density consumer requires density_normalized or raw_weighted pair factors"),
    )
    !isempty(term_coefficients) || throw(
        ArgumentError("route-shaped density-density consumer requires at least one pair-factor coefficient"),
    )
    coeffs = Float64[Float64(value) for value in term_coefficients]
    all(isfinite, coeffs) || throw(
        ArgumentError("route-shaped density-density term coefficients must be finite"),
    )
    selected_axis_pair_factor_terms =
        pair_factor_normalization == :raw_weighted ?
        raw_axis_pair_factor_terms : axis_pair_factor_terms
    selected_axis_pair_factor_terms isa NamedTuple{(:x,:y,:z)} || throw(
        ArgumentError("route-shaped density-density consumer requires caller-supplied $(pair_factor_normalization) pair factors"),
    )
    pair_factor_terms_symmetric = _source_box_axis_pair_terms_symmetric(
        selected_axis_pair_factor_terms;
        atol = pair_factor_symmetry_atol,
    )
    pair_factor_terms_symmetric || throw(
        ArgumentError("route-shaped density-density consumer requires symmetric pair-factor terms before transpose-only lower blocks are valid"),
    )

    left_raw_plan = _pqs_raw_product_box_plan_view(units.pqs_left)
    right_raw_plan = _pqs_raw_product_box_plan_view(units.pqs_right)
    left_raw_plan.representation == :orthogonal_raw_product_box || throw(
        ArgumentError("route-shaped density-density consumer requires a raw product-box left PQS plan"),
    )
    right_raw_plan.representation == :orthogonal_raw_product_box || throw(
        ArgumentError("route-shaped density-density consumer requires a raw product-box right PQS plan"),
    )
    units.product.kind == :product_doside || throw(
        ArgumentError("route-shaped density-density consumer requires a product_doside middle unit"),
    )

    range_info = _pqs_pqs_product_density_density_route_ranges(
        left_raw_plan,
        right_raw_plan,
        units.product,
    )
    ranges = range_info.ranges
    retained_dimension = range_info.retained_dimension
    inventory = _pqs_pqs_product_density_density_all_pairs_inventory(
        left_raw_plan,
        right_raw_plan,
        units.product,
        ranges,
        pair_factor_normalization,
    )
    inventory.diagnostics.every_pair_uses_source_box_algorithmic_policy ||
        throw(ArgumentError("density-density route requires all pairs to use source-box algorithms"))

    descriptor_expected_ranges_checked = false
    descriptor_retained_dimension_checked = false
    descriptor_pair_count_checked = false
    if hasproperty(route_units, :expected_ranges)
        route_units.expected_ranges == ranges || throw(
            ArgumentError("route descriptor expected ranges disagree with density-density route ranges"),
        )
        descriptor_expected_ranges_checked = true
    end
    if hasproperty(route_units, :retained_dimension)
        route_units.retained_dimension == retained_dimension || throw(
            DimensionMismatch("route descriptor retained dimension disagrees with density-density route"),
        )
        descriptor_retained_dimension_checked = true
    end
    if hasproperty(route_units, :expected_pair_count)
        route_units.expected_pair_count == length(inventory.pair_entries) || throw(
            ArgumentError("route descriptor expected pair count disagrees with density-density route"),
        )
        descriptor_pair_count_checked = true
    end

    timed = @timed begin
        block = zeros(Float64, retained_dimension, retained_dimension)
        pair_block_results = Dict{Tuple{Symbol,Symbol},Any}()
        for entry in inventory.pair_entries
            pair_result = _pqs_pqs_product_density_density_pair_block(
                entry,
                units;
                term_coefficients = coeffs,
                axis_pair_factor_terms,
                raw_axis_pair_factor_terms,
                axis_weights,
                pair_factor_normalization,
            )
            left_role, right_role = entry.pair_key
            left_range = getproperty(ranges, left_role)
            right_range = getproperty(ranges, right_role)
            block[left_range, right_range] .= pair_result.block
            if left_role != right_role
                entry.transpose_only_lower_block || throw(
                    ArgumentError("density-density route lower block requires transpose-only policy"),
                )
                block[right_range, left_range] .= transpose(pair_result.block)
            end
            pair_block_results[entry.pair_key] = pair_result
        end
        (block = block, pair_block_results = pair_block_results)
    end
    assembled = timed.value
    block = assembled.block
    pair_block_results = assembled.pair_block_results
    all(isfinite, block) || throw(
        ArgumentError("route-shaped density-density consumer produced non-finite entries"),
    )
    symmetry_error = LinearAlgebra.norm(block - transpose(block), Inf)
    symmetry_error <= symmetry_atol || throw(
        ArgumentError("route-shaped density-density consumer symmetry error exceeded tolerance"),
    )

    component_blocks = (
        pqs_left_pqs_left = pair_block_results[(:pqs_left, :pqs_left)].block,
        pqs_left_pqs_right = pair_block_results[(:pqs_left, :pqs_right)].block,
        pqs_right_pqs_left =
            transpose(pair_block_results[(:pqs_left, :pqs_right)].block),
        pqs_left_product = pair_block_results[(:pqs_left, :product)].block,
        product_pqs_left =
            transpose(pair_block_results[(:pqs_left, :product)].block),
        pqs_right_pqs_right = pair_block_results[(:pqs_right, :pqs_right)].block,
        pqs_right_product = pair_block_results[(:pqs_right, :product)].block,
        product_pqs_right =
            transpose(pair_block_results[(:pqs_right, :product)].block),
        product_product = pair_block_results[(:product, :product)].block,
    )
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
        term_count = length(coeffs),
        dense_raw_source_box_pair_matrix_materialized = false,
        dense_raw_pair_storage_avoided = true,
    )
    return (
        path = :pqs_pqs_product_route_shaped_density_density_consumer,
        route_kind = route_kind,
        route_name = route_name,
        route_units = route_units,
        retained_units = inventory.retained_units,
        all_pairs_inventory = inventory,
        block = block,
        density_density_matrix = block,
        complete_retained_space_matrix = block,
        component_blocks = component_blocks,
        component_block_provenance = component_block_provenance,
        pair_block_results = pair_block_results,
        ranges = ranges,
        retained_dimension = retained_dimension,
        pair_count = length(inventory.pair_entries),
        pair_family_counts = inventory.pair_family_counts,
        term_coefficients = coeffs,
        term_count = length(coeffs),
        pair_factor_normalization = pair_factor_normalization,
        output_finite = true,
        symmetry_error = symmetry_error,
        performance = performance,
        metadata = metadata,
        provenance = provenance,
        diagnostics = merge(
            inventory.diagnostics,
            (
                source = :pqs_pqs_product_route_shaped_density_density_consumer,
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
                term_count = length(coeffs),
                pair_factor_normalization = pair_factor_normalization,
                density_normalized_pair_factors =
                    pair_factor_normalization == :density_normalized,
                raw_weighted_pair_factors =
                    pair_factor_normalization == :raw_weighted,
                density_normalized_pair_factors_generated =
                    pair_factor_normalization == :raw_weighted,
                source_weight_division_owner =
                    pair_factor_normalization == :raw_weighted ?
                    :source_box_raw_weights :
                    :caller_supplied_density_normalized_pair_factors,
                source_weight_division_applied_by_helper =
                    pair_factor_normalization == :raw_weighted,
                source_weight_division_shape =
                    pair_factor_normalization == :raw_weighted ?
                    :axis_pair_weight_outer : nothing,
                output_representation = :two_index_density_density,
                four_index_galerkin_tensor = false,
                interaction_operator = :electron_electron_density_density,
                source_box_first = true,
                source_box_algorithmic_path_true_for_every_pair = true,
                every_pair_uses_source_box_algorithmic_policy =
                    inventory.diagnostics.every_pair_uses_source_box_algorithmic_policy,
                source_box_algorithmic_pair_count =
                    inventory.diagnostics.source_box_algorithmic_pair_count,
                helper_used_for_pair_families =
                    inventory.diagnostics.block_helper_by_family,
                block_helpers = inventory.diagnostics.block_helpers,
                product_pqs_blocks_transpose_only = true,
                lower_triangular_cross_blocks_transpose_only = true,
                pair_factor_terms_symmetric = pair_factor_terms_symmetric,
                pair_factor_symmetry_atol = Float64(pair_factor_symmetry_atol),
                symmetric_same_route_input = true,
                output_finite = true,
                symmetry_error = symmetry_error,
                symmetry_atol = Float64(symmetry_atol),
                descriptor_expected_ranges_checked =
                    descriptor_expected_ranges_checked,
                descriptor_retained_dimension_checked =
                    descriptor_retained_dimension_checked,
                descriptor_pair_count_checked = descriptor_pair_count_checked,
                private_shadow_only = true,
                input_pair_factor_data = :caller_supplied_explicit_data,
                input_pair_factor_data_pgdg_checked = false,
                real_mwg_ida_pair_factor_provenance_adapted = false,
                source_weights_are_raw_source_weights = true,
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
                ecp_terms_implemented = false,
                cr2_science_status_changed = false,
                ida_mwg_semantics_changed = false,
                mwg_ida_semantics_changed = false,
                mwg_interaction_implemented = false,
                complete_retained_space_matrix_built = true,
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

function _pqs_pqs_product_raw_box_density_density_route_producer(
    bundles,
    left_source_box::NTuple{3,UnitRange{Int}},
    right_source_box::NTuple{3,UnitRange{Int}},
    product_source_box::NTuple{3,UnitRange{Int}},
    metrics::NamedTuple{(:x,:y,:z)};
    source_mode_dims::NTuple{3,Int},
    term_coefficients::AbstractVector{<:Real},
    axis_weights::NamedTuple{(:x,:y,:z)},
    axis_pair_factor_terms = nothing,
    raw_axis_pair_factor_terms = nothing,
    pair_factor_normalization::Symbol = :density_normalized,
    pair_factor_symmetry_atol::Real = 1.0e-12,
    symmetry_atol::Real = 1.0e-10,
    route_name::Symbol =
        :homonuclear_pqs_product_source_box_density_density_fixture,
    parent_dims = _nested_axis_lengths(bundles),
    bond_axis = nothing,
    metadata = (;),
    provenance = (;),
    route_supported_terms = _PQS_PRODUCT_SOURCE_BOX_SHADOW_TERMS,
    orthogonality_atol::Real = 1.0e-8,
)
    pair_factor_normalization in (:density_normalized, :raw_weighted) || throw(
        ArgumentError("density-density route producer requires density_normalized or raw_weighted pair factors"),
    )
    producer_metadata = merge(
        (
            parent_dims = parent_dims,
            bond_axis = bond_axis,
            pqs_left_box = left_source_box,
            pqs_right_box = right_source_box,
            product_source_box = product_source_box,
            pqs_source_mode_dims = source_mode_dims,
            route_producer =
                :pqs_pqs_product_raw_box_density_density_route_producer,
            density_density_route_producer =
                :pqs_pqs_product_raw_box_density_density_route_producer,
            pair_factor_normalization = pair_factor_normalization,
            input_pair_factor_data = :caller_supplied_explicit_data,
            real_mwg_ida_pair_factor_provenance_adapted = false,
        ),
        metadata,
    )
    producer_provenance = merge(
        (source = :pqs_pqs_product_raw_box_density_density_route_producer,),
        provenance,
    )
    route = _pqs_pqs_product_raw_box_route_producer(
        bundles,
        left_source_box,
        right_source_box,
        product_source_box,
        metrics;
        source_mode_dims,
        route_name,
        parent_dims,
        bond_axis,
        metadata = producer_metadata,
        provenance = producer_provenance,
        supported_terms = route_supported_terms,
        orthogonality_atol,
    )
    consumer = _pqs_pqs_product_route_shaped_density_density_consumer(
        route.descriptor;
        term_coefficients,
        axis_pair_factor_terms,
        raw_axis_pair_factor_terms,
        axis_weights,
        pair_factor_normalization,
        pair_factor_symmetry_atol,
        symmetry_atol,
    )
    return (
        object_kind =
            :pqs_pqs_product_raw_box_density_density_route_producer,
        status = :private_density_density_reference_only,
        descriptor = route.descriptor,
        route_descriptor = route.descriptor,
        raw_box_route_producer = route,
        raw_product_box_plans = route.raw_product_box_plans,
        raw_pqs_plans = route.raw_pqs_plans,
        product_unit = route.product_unit,
        retained_rules = route.retained_rules,
        consumer = consumer,
        density_density_consumer = consumer,
        block = consumer.block,
        density_density_matrix = consumer.density_density_matrix,
        complete_retained_space_matrix = consumer.complete_retained_space_matrix,
        all_pairs_inventory = consumer.all_pairs_inventory,
        pair_block_results = consumer.pair_block_results,
        ranges = consumer.ranges,
        retained_dimension = consumer.retained_dimension,
        pair_count = consumer.pair_count,
        pair_family_counts = consumer.pair_family_counts,
        term_coefficients = consumer.term_coefficients,
        term_count = consumer.term_count,
        pair_factor_normalization = consumer.pair_factor_normalization,
        output_finite = consumer.output_finite,
        symmetry_error = consumer.symmetry_error,
        performance = consumer.performance,
        metadata = producer_metadata,
        provenance = producer_provenance,
        diagnostics = merge(
            route.diagnostics,
            consumer.diagnostics,
            (
                source =
                    :pqs_pqs_product_raw_box_density_density_route_producer,
                private_density_density_reference_only = true,
                private_shadow_only = true,
                raw_box_route_producer_called = true,
                route_descriptor_emitted = true,
                route_descriptor_source = route.diagnostics.source,
                route_descriptor_object_kind = route.descriptor.object_kind,
                route_descriptor_provenance = route.descriptor.provenance,
                route_descriptor_built_from_explicit_fixture_facts = true,
                density_density_consumer_called = true,
                density_density_consumer_path = consumer.path,
                returns_descriptor_and_density_density_consumer_result = true,
                retained_dimension = consumer.retained_dimension,
                pair_count = consumer.pair_count,
                pair_family_counts = consumer.pair_family_counts,
                pair_factor_normalization = consumer.pair_factor_normalization,
                input_pair_factor_data = :caller_supplied_explicit_data,
                synthetic_or_caller_supplied_pair_factors = true,
                real_mwg_ida_pair_factor_provenance_adapted = false,
                source_box_first = true,
                source_box_algorithmic_path_true_for_every_pair = true,
                every_pair_uses_source_box_algorithmic_policy =
                    consumer.diagnostics.every_pair_uses_source_box_algorithmic_policy,
                source_box_algorithmic_pair_count =
                    consumer.diagnostics.source_box_algorithmic_pair_count,
                shell_projection_used = false,
                lowdin_cleanup_used = false,
                support_local_oracle_used = false,
                support_local_pqs_oracle_used = false,
                support_local_shell_row_algorithm = false,
                support_coefficient_matrix_used = false,
                shell_row_algorithm = false,
                retained_pqs_weights_used = false,
                retained_pqs_weights_positive_checked = false,
                retained_weight_division_allowed = false,
                retained_pqs_weight_division_allowed = false,
                ida_weight_division_allowed = false,
                retained_weight_semantics = :not_positive_quadrature_weights,
                packet_adoption = false,
                fixed_block_routing = false,
                qwhamiltonian_consumes = false,
                public_default_consumes = false,
                ecp_terms_implemented = false,
                cr2_science_status_changed = false,
                ida_mwg_semantics_changed = false,
                mwg_ida_semantics_changed = false,
                mwg_interaction_implemented = false,
            ),
        ),
    )
end

function _pqs_pqs_product_raw_box_density_density_route_producer_from_ida_provenance(
    bundles,
    left_source_box::NTuple{3,UnitRange{Int}},
    right_source_box::NTuple{3,UnitRange{Int}},
    product_source_box::NTuple{3,UnitRange{Int}},
    metrics::NamedTuple{(:x,:y,:z)},
    ida_provenance;
    source_mode_dims::NTuple{3,Int},
    term_coefficients::AbstractVector{<:Real},
    pair_factor_normalization::Symbol = :density_normalized,
    pair_factor_symmetry_atol::Real = 1.0e-12,
    symmetry_atol::Real = 1.0e-10,
    route_name::Symbol =
        :homonuclear_pqs_product_source_box_density_density_ida_provenance_fixture,
    parent_dims = _nested_axis_lengths(bundles),
    bond_axis = nothing,
    metadata = (;),
    provenance = (;),
    route_supported_terms = _PQS_PRODUCT_SOURCE_BOX_SHADOW_TERMS,
    orthogonality_atol::Real = 1.0e-8,
)
    pair_factor_normalization in (:density_normalized, :raw_weighted) || throw(
        ArgumentError("IDA provenance density-density route producer requires density_normalized or raw_weighted pair factors"),
    )
    hasproperty(ida_provenance, :object_kind) &&
        ida_provenance.object_kind == :pqs_source_box_ida_factor_provenance || throw(
            ArgumentError("density-density IDA provenance route producer requires _pqs_source_box_ida_factor_provenance output"),
        )
    hasproperty(ida_provenance, :axis_pair_factor_terms) &&
        ida_provenance.axis_pair_factor_terms isa NamedTuple{(:x,:y,:z)} || throw(
            ArgumentError("density-density IDA provenance route producer requires density-normalized x/y/z pair-factor terms"),
        )
    hasproperty(ida_provenance, :raw_axis_pair_factor_terms) &&
        ida_provenance.raw_axis_pair_factor_terms isa NamedTuple{(:x,:y,:z)} || throw(
            ArgumentError("density-density IDA provenance route producer requires raw x/y/z pair-factor terms"),
        )
    hasproperty(ida_provenance, :axis_weights) &&
        ida_provenance.axis_weights isa NamedTuple{(:x,:y,:z)} || throw(
            ArgumentError("density-density IDA provenance route producer requires x/y/z source weights"),
        )
    hasproperty(ida_provenance, :diagnostics) || throw(
        ArgumentError("density-density IDA provenance route producer requires provenance diagnostics"),
    )
    ida_provenance.diagnostics.interaction_path == :ida_gausslet_source_box || throw(
        ArgumentError("density-density IDA provenance route producer requires IDA gausslet/source-box provenance"),
    )
    !ida_provenance.diagnostics.mwg_supplement_residual_path || throw(
        ArgumentError("density-density IDA provenance route producer cannot consume MWG supplement/residual provenance"),
    )
    coeffs = Float64[Float64(value) for value in term_coefficients]
    length(coeffs) == ida_provenance.term_count || throw(
        DimensionMismatch("density-density IDA provenance route producer term coefficients must match provenance term count"),
    )
    all(isfinite, coeffs) || throw(
        ArgumentError("density-density IDA provenance route producer term coefficients must be finite"),
    )

    axis_names = (:x, :y, :z)
    parent_dims_value = Int.(parent_dims)
    for axis in 1:3
        axis_name = axis_names[axis]
        density_terms = getproperty(ida_provenance.axis_pair_factor_terms, axis_name)
        raw_terms = getproperty(ida_provenance.raw_axis_pair_factor_terms, axis_name)
        weights = getproperty(ida_provenance.axis_weights, axis_name)
        size(density_terms) == size(raw_terms) || throw(
            DimensionMismatch("density-density IDA provenance raw and density-normalized $(axis_name) terms must have matching shapes"),
        )
        size(density_terms, 1) == length(coeffs) || throw(
            DimensionMismatch("density-density IDA provenance $(axis_name) term count must match coefficients"),
        )
        size(density_terms, 2) == size(density_terms, 3) || throw(
            DimensionMismatch("density-density IDA provenance $(axis_name) pair-factor terms must be square"),
        )
        size(density_terms, 2) == parent_dims_value[axis] || throw(
            DimensionMismatch("density-density IDA provenance $(axis_name) factor dimension must match parent_dims"),
        )
        length(weights) == parent_dims_value[axis] || throw(
            DimensionMismatch("density-density IDA provenance $(axis_name) weights must match parent_dims"),
        )
        all(isfinite, weights) || throw(
            ArgumentError("density-density IDA provenance $(axis_name) weights must be finite"),
        )
        all(weight -> abs(weight) > 1.0e-12, weights) || throw(
            ArgumentError("density-density IDA provenance $(axis_name) weights must be nonzero"),
        )
        for source_box in (left_source_box, right_source_box, product_source_box)
            first(source_box[axis]) >= 1 &&
                last(source_box[axis]) <= parent_dims_value[axis] || throw(
                    ArgumentError("density-density IDA provenance route source boxes must fit pair-factor dimensions"),
                )
        end
    end

    selected_axis_pair_factor_terms =
        pair_factor_normalization == :density_normalized ?
        ida_provenance.axis_pair_factor_terms : nothing
    selected_raw_axis_pair_factor_terms =
        pair_factor_normalization == :raw_weighted ?
        ida_provenance.raw_axis_pair_factor_terms : nothing
    adapter_metadata = merge(
        (
            input_pair_factor_data = :ida_gausslet_source_box_provenance,
            interaction_path = :ida_gausslet_source_box,
            mwg_supplement_residual_path = false,
            ida_provenance_object_kind = ida_provenance.object_kind,
            ida_provenance_term_count = ida_provenance.term_count,
            ida_provenance_factor_dimensions = ida_provenance.factor_dimensions,
            density_density_ida_provenance_route_adapter = true,
        ),
        metadata,
    )
    adapter_provenance = merge(
        (
            source =
                :pqs_pqs_product_raw_box_density_density_route_producer_from_ida_provenance,
            ida_provenance_object_kind = ida_provenance.object_kind,
            ida_provenance_diagnostics = ida_provenance.diagnostics,
        ),
        provenance,
    )
    explicit_route =
        _pqs_pqs_product_raw_box_density_density_route_producer(
            bundles,
            left_source_box,
            right_source_box,
            product_source_box,
            metrics;
            source_mode_dims,
            term_coefficients = coeffs,
            axis_weights = ida_provenance.axis_weights,
            axis_pair_factor_terms = selected_axis_pair_factor_terms,
            raw_axis_pair_factor_terms = selected_raw_axis_pair_factor_terms,
            pair_factor_normalization,
            pair_factor_symmetry_atol,
            symmetry_atol,
            route_name,
            parent_dims = parent_dims_value,
            bond_axis,
            metadata = adapter_metadata,
            provenance = adapter_provenance,
            route_supported_terms,
            orthogonality_atol,
        )
    adapter_diagnostics = merge(
        explicit_route.diagnostics,
        (
            source =
                :pqs_pqs_product_raw_box_density_density_route_producer_from_ida_provenance,
            private_density_density_reference_only = true,
            private_shadow_only = true,
            density_density_ida_provenance_route_adapter = true,
            explicit_input_route_producer_delegate =
                :_pqs_pqs_product_raw_box_density_density_route_producer,
            input_pair_factor_data = :ida_gausslet_source_box_provenance,
            input_pair_factor_data_pgdg_checked = true,
            ida_provenance_object_kind = ida_provenance.object_kind,
            ida_provenance_term_count = ida_provenance.term_count,
            ida_provenance_factor_dimensions = ida_provenance.factor_dimensions,
            interaction_path = :ida_gausslet_source_box,
            mwg_supplement_residual_path = false,
            mwg_supplement_residual_provenance_adapted = false,
            real_ida_gausslet_source_box_provenance_adapted = true,
            real_mwg_ida_pair_factor_provenance_adapted = false,
            synthetic_or_caller_supplied_pair_factors = false,
            source_weight_division_owner =
                pair_factor_normalization == :raw_weighted ?
                :pgdg_auxiliary_source_weights :
                :caller_supplied_density_normalized_pair_factors,
            source_weights_are_raw_source_weights = true,
            retained_pqs_weights_used = false,
            retained_pqs_weights_positive_checked = false,
            retained_weight_division_allowed = false,
            retained_pqs_weight_division_allowed = false,
            ida_weight_division_allowed = false,
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
            ecp_terms_implemented = false,
            cr2_science_status_changed = false,
            ida_mwg_semantics_changed = false,
            mwg_ida_semantics_changed = false,
            mwg_interaction_implemented = false,
        ),
    )
    return merge(
        explicit_route,
        (
            object_kind =
                :pqs_pqs_product_raw_box_density_density_route_producer_from_ida_provenance,
            path =
                :pqs_pqs_product_raw_box_density_density_route_producer_from_ida_provenance,
            status = :private_ida_provenance_density_density_reference_only,
            explicit_input_route_producer = explicit_route,
            ida_provenance = ida_provenance,
            metadata = adapter_metadata,
            provenance = adapter_provenance,
            diagnostics = adapter_diagnostics,
        ),
    )
end

function _pqs_route_parent_dims(route_units, parent_dims)
    !isnothing(parent_dims) && return _pqs_geometry_int3(parent_dims, :parent_dims)
    hasproperty(route_units, :metadata) &&
        hasproperty(route_units.metadata, :parent_dims) ||
        throw(ArgumentError("route parent coefficient projection requires parent_dims"))
    metadata_parent_dims = route_units.metadata.parent_dims
    !isnothing(metadata_parent_dims) ||
        throw(ArgumentError("route parent coefficient projection requires non-nothing parent_dims"))
    return _pqs_geometry_int3(metadata_parent_dims, :parent_dims)
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

function _product_doside_parent_coefficient_matrix(
    product_unit::_CartesianNestedProductStagedByCenterUnit3D,
    parent_dims::NTuple{3,Int},
)
    _require_product_doside_retained_block_unit(product_unit; side = :product)
    parent_dimension = prod(parent_dims)
    local_coefficients = Matrix{Float64}(product_unit.coefficient_matrix)
    size(local_coefficients, 1) == length(product_unit.support_indices) || throw(
        DimensionMismatch("product/doside parent coefficient projection support rows must match coefficient rows"),
    )
    size(local_coefficients, 2) == length(product_unit.column_range) || throw(
        DimensionMismatch("product/doside parent coefficient projection retained columns must match column range"),
    )
    coefficients = zeros(Float64, parent_dimension, size(local_coefficients, 2))
    @inbounds for (local_row, parent_row) in pairs(product_unit.support_indices)
        1 <= parent_row <= parent_dimension || throw(
            ArgumentError("product/doside parent coefficient projection support index exceeds parent dimension"),
        )
        for column in axes(local_coefficients, 2)
            value = local_coefficients[local_row, column]
            iszero(value) && continue
            coefficients[parent_row, column] += value
        end
    end
    return coefficients
end

function _pqs_pqs_product_route_parent_coefficient_matrix(route_units; parent_dims = nothing)
    hasproperty(route_units, :object_kind) &&
        route_units.object_kind == :pqs_pqs_product_safe_term_route_descriptor ||
        throw(ArgumentError("parent coefficient projection requires a PQS/PQS/product route descriptor"))
    parent_dims_value = _pqs_route_parent_dims(route_units, parent_dims)
    roles = hasproperty(route_units, :roles) ?
        Tuple(route_units.roles) :
        (:pqs_left, :pqs_right, :product)
    roles == (:pqs_left, :pqs_right, :product) || throw(
        ArgumentError("parent coefficient projection requires roles (:pqs_left, :pqs_right, :product)"),
    )
    hasproperty(route_units, :units) || throw(
        ArgumentError("parent coefficient projection requires route units"),
    )
    units = route_units.units
    ranges = hasproperty(route_units, :expected_ranges) ?
        route_units.expected_ranges : route_units.ranges
    retained_dimension = route_units.retained_dimension
    parent_dimension = prod(parent_dims_value)
    coefficients = zeros(Float64, parent_dimension, retained_dimension)
    unit_coefficients = (
        pqs_left = _pqs_parent_coefficient_matrix_from_raw_plan(
            units.pqs_left,
            parent_dims_value,
        ),
        pqs_right = _pqs_parent_coefficient_matrix_from_raw_plan(
            units.pqs_right,
            parent_dims_value,
        ),
        product = _product_doside_parent_coefficient_matrix(
            units.product,
            parent_dims_value,
        ),
    )
    size(unit_coefficients.pqs_left, 2) == length(ranges.pqs_left) ||
        throw(DimensionMismatch("left PQS coefficient columns must match route range"))
    size(unit_coefficients.pqs_right, 2) == length(ranges.pqs_right) ||
        throw(DimensionMismatch("right PQS coefficient columns must match route range"))
    size(unit_coefficients.product, 2) == length(ranges.product) ||
        throw(DimensionMismatch("product coefficient columns must match route range"))
    coefficients[:, ranges.pqs_left] .= unit_coefficients.pqs_left
    coefficients[:, ranges.pqs_right] .= unit_coefficients.pqs_right
    coefficients[:, ranges.product] .= unit_coefficients.product
    return (
        object_kind = :pqs_pqs_product_route_parent_coefficient_matrix,
        coefficient_matrix = coefficients,
        parent_dims = parent_dims_value,
        parent_dimension = parent_dimension,
        retained_dimension = retained_dimension,
        ranges = ranges,
        unit_coefficient_shapes = (
            pqs_left = size(unit_coefficients.pqs_left),
            pqs_right = size(unit_coefficients.pqs_right),
            product = size(unit_coefficients.product),
        ),
        diagnostics = (
            source = :pqs_pqs_product_route_parent_coefficient_matrix,
            route_parent_projection_private = true,
            source_box_algorithm_changed = false,
            support_local_algorithm_used = false,
            support_coefficient_matrix_used = false,
            shell_projection_used = false,
            lowdin_cleanup_used = false,
            retained_pqs_weights_used = false,
            ida_weight_division_allowed = false,
            packet_adoption = false,
            fixed_block_routing = false,
            qwhamiltonian_consumes = false,
            public_default_consumes = false,
        ),
    )
end

function _pqs_pqs_product_dense_parent_ida_authority_comparison(
    route_result,
    dense_parent_ida_matrix::AbstractMatrix{<:Real};
    parent_dims = nothing,
    dense_parent_matrix_source::Symbol = :caller_supplied_dense_parent_ida_matrix,
    comparison_atol::Real = 1.0e-10,
)
    hasproperty(route_result, :route_descriptor) || throw(
        ArgumentError("dense-parent IDA authority comparison requires a route_result with route_descriptor"),
    )
    hasproperty(route_result, :block) || throw(
        ArgumentError("dense-parent IDA authority comparison requires a route_result block"),
    )
    hasproperty(route_result, :diagnostics) || throw(
        ArgumentError("dense-parent IDA authority comparison requires route diagnostics"),
    )
    route_result.diagnostics.input_pair_factor_data ==
        :ida_gausslet_source_box_provenance || throw(
            ArgumentError("dense-parent IDA authority comparison requires IDA source-box provenance route output"),
        )
    route_result.diagnostics.interaction_path == :ida_gausslet_source_box || throw(
        ArgumentError("dense-parent IDA authority comparison requires IDA gausslet/source-box interaction path"),
    )
    !route_result.diagnostics.mwg_supplement_residual_path || throw(
        ArgumentError("dense-parent IDA authority comparison cannot consume MWG supplement/residual provenance"),
    )
    route_result.pair_factor_normalization == :density_normalized || throw(
        ArgumentError("dense-parent IDA authority comparison currently covers density-normalized route output"),
    )
    projection = _pqs_pqs_product_route_parent_coefficient_matrix(
        route_result.route_descriptor;
        parent_dims,
    )
    parent_matrix = Matrix{Float64}(dense_parent_ida_matrix)
    size(parent_matrix) == (projection.parent_dimension, projection.parent_dimension) ||
        throw(DimensionMismatch("dense-parent IDA authority matrix must match route parent dimension"))
    route_block = Matrix{Float64}(route_result.block)
    size(route_block) ==
        (projection.retained_dimension, projection.retained_dimension) ||
        throw(DimensionMismatch("route density-density block must match retained projection dimension"))
    coefficients = projection.coefficient_matrix
    projected_block = Matrix{Float64}(
        transpose(coefficients) * parent_matrix * coefficients,
    )
    max_error = LinearAlgebra.norm(projected_block - route_block, Inf)
    symmetry_error =
        LinearAlgebra.norm(projected_block - transpose(projected_block), Inf)
    output_finite =
        all(isfinite, projected_block) && all(isfinite, route_block) &&
        isfinite(max_error) && isfinite(symmetry_error)
    comparison_atol_value = Float64(comparison_atol)
    return (
        object_kind = :pqs_pqs_product_dense_parent_ida_authority_comparison,
        status = :private_validation_only,
        authority_kind = :dense_parent_ida_projection,
        route_result = route_result,
        parent_projection = projection,
        coefficient_matrix = coefficients,
        dense_parent_ida_matrix = parent_matrix,
        projected_block = projected_block,
        route_block = route_block,
        max_error = max_error,
        symmetry_error = symmetry_error,
        output_finite = output_finite,
        within_tolerance = output_finite && max_error <= comparison_atol_value,
        comparison_atol = comparison_atol_value,
        retained_dimension = projection.retained_dimension,
        parent_dimension = projection.parent_dimension,
        diagnostics = (
            source = :pqs_pqs_product_dense_parent_ida_authority_comparison,
            authority_kind = :dense_parent_ida_projection,
            dense_parent_matrix_source = dense_parent_matrix_source,
            dense_parent_matrix_used_for_validation = true,
            dense_parent_matrix_algorithmic = false,
            source_box_algorithm_changed = false,
            support_local_algorithm_used = false,
            support_local_oracle_used = false,
            support_coefficient_matrix_used = false,
            shell_projection_used = false,
            lowdin_cleanup_used = false,
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
            mwg_supplement_residual_path = false,
            mwg_supplement_residual_provenance_adapted = false,
            output_representation = :two_index_density_density,
            four_index_galerkin_tensor = false,
            parent_projection_shape = size(coefficients),
            route_block_shape = size(route_block),
            max_error = max_error,
            symmetry_error = symmetry_error,
            output_finite = output_finite,
            within_tolerance = output_finite && max_error <= comparison_atol_value,
        ),
    )
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

function _pqs_pqs_product_source_box_component_route_smoke(
    bundles,
    left_source_box::NTuple{3,UnitRange{Int}},
    right_source_box::NTuple{3,UnitRange{Int}},
    product_source_box::NTuple{3,UnitRange{Int}},
    metrics::NamedTuple{(:x,:y,:z)},
    ida_provenance,
    nuclear_axis_layers::NamedTuple{(:x,:y,:z)},
    nuclear_expansion::CoulombGaussianExpansion;
    source_mode_dims::NTuple{3,Int},
    centers,
    nuclear_charges,
    center_labels = nothing,
    term_coefficients::AbstractVector{<:Real},
    pair_factor_normalization::Symbol = :density_normalized,
    pair_factor_symmetry_atol::Real = 1.0e-12,
    symmetry_atol::Real = 1.0e-10,
    nuclear_symmetry_atol::Real = 1.0e-10,
    route_name::Symbol = :pqs_pqs_product_source_box_component_route_smoke,
    parent_dims = _nested_axis_lengths(bundles),
    bond_axis = nothing,
    metadata = (;),
    provenance = (;),
    dense_parent_ida_matrix = nothing,
    dense_parent_matrix_source::Symbol = :caller_supplied_dense_parent_ida_matrix,
    dense_parent_comparison_atol::Real = 1.0e-10,
    route_supported_terms = _PQS_PRODUCT_SOURCE_BOX_SHADOW_TERMS,
    orthogonality_atol::Real = 1.0e-8,
)
    pair_factor_normalization in (:density_normalized, :raw_weighted) || throw(
        ArgumentError("component route smoke requires density_normalized or raw_weighted IDA route output"),
    )
    center_values = _source_box_nuclear_attraction_center_values(centers)
    charge_values = _source_box_nuclear_attraction_charge_values(nuclear_charges)
    length(center_values) == length(charge_values) || throw(
        ArgumentError("component route smoke center and charge counts must match"),
    )
    labels = _source_box_nuclear_attraction_center_labels(
        length(center_values),
        center_labels,
    )
    component_metadata = merge(
        (
            component_route_smoke = true,
            nuclear_center_labels = labels,
            nuclear_center_count = length(center_values),
            electron_electron_source = :ida_gausslet_source_box,
        ),
        metadata,
    )
    component_provenance = merge(
        (source = :pqs_pqs_product_source_box_component_route_smoke,),
        provenance,
    )
    ida_route =
        _pqs_pqs_product_raw_box_density_density_route_producer_from_ida_provenance(
            bundles,
            left_source_box,
            right_source_box,
            product_source_box,
            metrics,
            ida_provenance;
            source_mode_dims,
            term_coefficients,
            pair_factor_normalization,
            pair_factor_symmetry_atol,
            symmetry_atol,
            route_name,
            parent_dims,
            bond_axis,
            metadata = component_metadata,
            provenance = component_provenance,
            route_supported_terms,
            orthogonality_atol,
        )
    nuclear_route = _pqs_pqs_product_route_shaped_nuclear_attraction_by_center(
        ida_route.route_descriptor,
        nuclear_axis_layers,
        nuclear_expansion;
        centers = center_values,
        nuclear_charges = charge_values,
        center_labels = labels,
        symmetry_atol = nuclear_symmetry_atol,
    )
    dense_parent_authority_supported =
        pair_factor_normalization == :density_normalized
    dense_parent_authority_skipped_reason =
        isnothing(dense_parent_ida_matrix) || dense_parent_authority_supported ?
        nothing : :density_normalized_authority_only
    dense_parent_authority =
        isnothing(dense_parent_ida_matrix) || !dense_parent_authority_supported ? nothing :
        _pqs_pqs_product_dense_parent_ida_authority_comparison(
            ida_route,
            dense_parent_ida_matrix;
            parent_dims,
            dense_parent_matrix_source,
            comparison_atol = dense_parent_comparison_atol,
        )
    output_finite =
        nuclear_route.output_finite &&
        ida_route.output_finite &&
        (isnothing(dense_parent_authority) || dense_parent_authority.output_finite)
    authority_max_errors = (
        electron_electron_dense_parent =
            isnothing(dense_parent_authority) ? nothing :
            dense_parent_authority.max_error,
        nuclear_total_from_center_sum = nuclear_route.total_from_center_error,
    )
    return (
        object_kind = :pqs_pqs_product_source_box_component_route_smoke,
        path = :pqs_pqs_product_source_box_component_route_smoke,
        status = :private_component_route_smoke,
        descriptor = ida_route.route_descriptor,
        route_descriptor = ida_route.route_descriptor,
        raw_box_route_producer = ida_route.raw_box_route_producer,
        route_shape = (:pqs_left, :pqs_right, :product),
        retained_units = ida_route.route_descriptor.unit_summaries,
        ranges = ida_route.ranges,
        retained_dimension = ida_route.retained_dimension,
        nuclear_attraction_by_center = nuclear_route,
        per_center_nuclear_matrices = nuclear_route.blocks_by_center,
        nuclear_total_matrix = nuclear_route.total_block,
        electron_electron_density_density = ida_route,
        electron_electron_matrix = ida_route.block,
        dense_parent_ida_authority = dense_parent_authority,
        center_labels = labels,
        centers = Tuple(center_values),
        nuclear_charges = Tuple(charge_values),
        ida_term_count = ida_route.term_count,
        pair_factor_normalization = pair_factor_normalization,
        output_finite = output_finite,
        nuclear_symmetry_error = nuclear_route.symmetry_error,
        electron_electron_symmetry_error = ida_route.symmetry_error,
        authority_max_errors = authority_max_errors,
        metadata = component_metadata,
        provenance = component_provenance,
        diagnostics = (
            source = :pqs_pqs_product_source_box_component_route_smoke,
            private_component_route_smoke = true,
            private_shadow_only = true,
            route_shape = (:pqs_left, :pqs_right, :product),
            route_name = ida_route.route_descriptor.route_name,
            route_descriptor_object_kind = ida_route.route_descriptor.object_kind,
            route_roles = ida_route.route_descriptor.roles,
            retained_ranges = ida_route.ranges,
            retained_dimension = ida_route.retained_dimension,
            retained_unit_count = ida_route.route_descriptor.retained_unit_count,
            unit_roles = map(unit -> unit.unit_key, ida_route.route_descriptor.unit_summaries),
            unit_retained_ranges =
                map(unit -> unit.retained_range, ida_route.route_descriptor.unit_summaries),
            nuclear_pair_count = nuclear_route.pair_count,
            nuclear_pair_family_counts = nuclear_route.pair_family_counts,
            electron_electron_pair_count = ida_route.pair_count,
            electron_electron_pair_family_counts = ida_route.pair_family_counts,
            helper_used_for_nuclear_pair_families =
                nuclear_route.diagnostics.helper_used_for_pair_families,
            helper_used_for_electron_electron_pair_families =
                ida_route.diagnostics.helper_used_for_pair_families,
            center_labels = labels,
            centers = Tuple(center_values),
            nuclear_charges = Tuple(charge_values),
            center_records = nuclear_route.diagnostics.center_records,
            center_count = length(center_values),
            by_center_nuclear_terms_preserved = true,
            nuclear_total_is_explicit_sum_of_center_pieces = true,
            nuclear_total_from_center_error = nuclear_route.total_from_center_error,
            ida_term_count = ida_route.term_count,
            ida_provenance_term_count = ida_provenance.term_count,
            pair_factor_normalization = pair_factor_normalization,
            density_normalized_pair_factors =
                pair_factor_normalization == :density_normalized,
            raw_weighted_pair_factors =
                pair_factor_normalization == :raw_weighted,
            density_normalized_pair_factors_generated =
                pair_factor_normalization == :raw_weighted,
            source_weight_division_owner =
                ida_route.diagnostics.source_weight_division_owner,
            source_weight_division_applied_by_helper =
                ida_route.diagnostics.source_weight_division_applied_by_helper,
            source_weight_division_shape =
                ida_route.diagnostics.source_weight_division_shape,
            source_weights_are_raw_source_weights =
                ida_route.diagnostics.source_weights_are_raw_source_weights,
            input_pair_factor_data = :ida_gausslet_source_box_provenance,
            interaction_path = :ida_gausslet_source_box,
            mwg_supplement_residual_path = false,
            real_ida_gausslet_source_box_provenance_adapted = true,
            real_mwg_ida_pair_factor_provenance_adapted = false,
            source_box_first = true,
            source_box_algorithmic_path_true_for_every_pair =
                nuclear_route.diagnostics.source_box_algorithmic_path_true_for_every_pair &&
                ida_route.diagnostics.source_box_algorithmic_path_true_for_every_pair,
            nuclear_source_box_algorithmic_pair_count =
                nuclear_route.diagnostics.source_box_algorithmic_pair_count,
            electron_electron_source_box_algorithmic_pair_count =
                ida_route.diagnostics.source_box_algorithmic_pair_count,
            output_finite = output_finite,
            nuclear_symmetry_error = nuclear_route.symmetry_error,
            electron_electron_symmetry_error = ida_route.symmetry_error,
            dense_parent_ida_authority_available = !isnothing(dense_parent_authority),
            dense_parent_ida_authority_supported_for_pair_factor_normalization =
                dense_parent_authority_supported,
            dense_parent_ida_authority_skipped_reason =
                dense_parent_authority_skipped_reason,
            dense_parent_ida_authority_max_error =
                isnothing(dense_parent_authority) ? nothing :
                dense_parent_authority.max_error,
            dense_parent_ida_authority_within_tolerance =
                isnothing(dense_parent_authority) ? nothing :
                dense_parent_authority.within_tolerance,
            dense_parent_projection_validation_only =
                !isnothing(dense_parent_authority),
            dense_parent_projection_algorithmic = false,
            output_representation = :component_matrices,
            electron_electron_output_representation = :two_index_density_density,
            electron_electron_four_index_galerkin_tensor = false,
            hamiltonian_matrix_built = false,
            packet_adoption = false,
            fixed_block_routing = false,
            qwhamiltonian_consumes = false,
            public_default_consumes = false,
            ecp_terms_implemented = false,
            cr2_science_status_changed = false,
            mwg_interaction_implemented = false,
            mwg_supplement_residual_provenance_adapted = false,
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
        ),
    )
end

function _pqs_pqs_product_component_route_smoke_summary(component)
    hasproperty(component, :object_kind) &&
        component.object_kind == :pqs_pqs_product_source_box_component_route_smoke ||
        throw(
            ArgumentError("component route smoke summary requires _pqs_pqs_product_source_box_component_route_smoke output"),
        )
    hasproperty(component, :nuclear_attraction_by_center) || throw(
        ArgumentError("component route smoke summary requires nuclear_attraction_by_center"),
    )
    hasproperty(component, :electron_electron_density_density) || throw(
        ArgumentError("component route smoke summary requires electron_electron_density_density"),
    )
    hasproperty(component, :diagnostics) || throw(
        ArgumentError("component route smoke summary requires diagnostics"),
    )
    nuclear = component.nuclear_attraction_by_center
    electron_electron = component.electron_electron_density_density
    diagnostics = component.diagnostics
    no_go_diagnostics = (
        source_box_first = diagnostics.source_box_first,
        source_box_algorithmic_path_true_for_every_pair =
            diagnostics.source_box_algorithmic_path_true_for_every_pair,
        shell_projection_used = diagnostics.shell_projection_used,
        lowdin_cleanup_used = diagnostics.lowdin_cleanup_used,
        support_local_oracle_used = diagnostics.support_local_oracle_used,
        support_local_pqs_oracle_used =
            diagnostics.support_local_pqs_oracle_used,
        support_local_shell_row_algorithm =
            diagnostics.support_local_shell_row_algorithm,
        support_coefficient_matrix_used =
            diagnostics.support_coefficient_matrix_used,
        shell_row_algorithm = diagnostics.shell_row_algorithm,
        retained_pqs_weights_used = diagnostics.retained_pqs_weights_used,
        retained_pqs_weights_positive_checked =
            diagnostics.retained_pqs_weights_positive_checked,
        retained_weight_division_allowed =
            diagnostics.retained_weight_division_allowed,
        retained_pqs_weight_division_allowed =
            diagnostics.retained_pqs_weight_division_allowed,
        ida_weight_division_allowed = diagnostics.ida_weight_division_allowed,
        packet_adoption = diagnostics.packet_adoption,
        fixed_block_routing = diagnostics.fixed_block_routing,
        qwhamiltonian_consumes = diagnostics.qwhamiltonian_consumes,
        hamiltonian_matrix_built = diagnostics.hamiltonian_matrix_built,
        public_default_consumes = diagnostics.public_default_consumes,
        mwg_supplement_residual_path =
            diagnostics.mwg_supplement_residual_path,
        mwg_supplement_residual_provenance_adapted =
            diagnostics.mwg_supplement_residual_provenance_adapted,
        ecp_terms_implemented = diagnostics.ecp_terms_implemented,
        cr2_science_status_changed = diagnostics.cr2_science_status_changed,
        dense_parent_projection_algorithmic =
            diagnostics.dense_parent_projection_algorithmic,
    )
    return (
        object_kind = :pqs_pqs_product_component_route_smoke_summary,
        status = :private_component_route_smoke_summary,
        component_object_kind = component.object_kind,
        component_status = component.status,
        route_shape = component.route_shape,
        retained_dimension = component.retained_dimension,
        center_labels = component.center_labels,
        nuclear_charges = component.nuclear_charges,
        center_count = length(component.center_labels),
        nuclear_pair_count = nuclear.pair_count,
        nuclear_pair_family_counts = nuclear.pair_family_counts,
        electron_electron_pair_count = electron_electron.pair_count,
        electron_electron_pair_family_counts =
            electron_electron.pair_family_counts,
        pair_factor_normalization = component.pair_factor_normalization,
        source_weight_division_owner =
            diagnostics.source_weight_division_owner,
        source_weight_division_applied_by_helper =
            diagnostics.source_weight_division_applied_by_helper,
        source_weight_division_shape =
            diagnostics.source_weight_division_shape,
        ida_term_count = component.ida_term_count,
        finite_checks = (
            output_finite = component.output_finite,
            nuclear_output_finite = nuclear.output_finite,
            electron_electron_output_finite = electron_electron.output_finite,
        ),
        symmetry_errors = (
            nuclear = component.nuclear_symmetry_error,
            electron_electron = component.electron_electron_symmetry_error,
        ),
        nuclear_total_from_center_error =
            nuclear.total_from_center_error,
        dense_parent_ida_authority = (
            available =
                diagnostics.dense_parent_ida_authority_available,
            max_error =
                diagnostics.dense_parent_ida_authority_max_error,
            within_tolerance =
                diagnostics.dense_parent_ida_authority_within_tolerance,
            skip_reason =
                diagnostics.dense_parent_ida_authority_skipped_reason,
            validation_only =
                diagnostics.dense_parent_projection_validation_only,
        ),
        no_go_diagnostics = no_go_diagnostics,
        performance = (
            nuclear = (
                elapsed_seconds = nuclear.performance.elapsed_seconds,
                allocated_bytes = nuclear.performance.allocated_bytes,
                gc_time_seconds = nuclear.performance.gc_time_seconds,
            ),
            electron_electron = (
                elapsed_seconds =
                    electron_electron.performance.elapsed_seconds,
                allocated_bytes =
                    electron_electron.performance.allocated_bytes,
                gc_time_seconds =
                    electron_electron.performance.gc_time_seconds,
            ),
        ),
        diagnostics = (
            source = :pqs_pqs_product_component_route_smoke_summary,
            private_component_route_smoke_summary = true,
            component_object_kind = component.object_kind,
            component_status = component.status,
            route_shape = component.route_shape,
            retained_dimension = component.retained_dimension,
            center_count = length(component.center_labels),
            pair_factor_normalization = component.pair_factor_normalization,
            output_finite = component.output_finite,
            source_box_first = no_go_diagnostics.source_box_first,
            source_box_algorithmic_path_true_for_every_pair =
                no_go_diagnostics.source_box_algorithmic_path_true_for_every_pair,
            retained_pqs_weights_used =
                no_go_diagnostics.retained_pqs_weights_used,
            retained_weight_division_allowed =
                no_go_diagnostics.retained_weight_division_allowed,
            ida_weight_division_allowed =
                no_go_diagnostics.ida_weight_division_allowed,
            packet_adoption = no_go_diagnostics.packet_adoption,
            fixed_block_routing = no_go_diagnostics.fixed_block_routing,
            qwhamiltonian_consumes = no_go_diagnostics.qwhamiltonian_consumes,
            public_default_consumes =
                no_go_diagnostics.public_default_consumes,
            mwg_supplement_residual_provenance_adapted =
                no_go_diagnostics.mwg_supplement_residual_provenance_adapted,
            ecp_terms_implemented = no_go_diagnostics.ecp_terms_implemented,
            cr2_science_status_changed =
                no_go_diagnostics.cr2_science_status_changed,
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

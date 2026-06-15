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

function _pqs_source_box_expected_ida_term_count(expected_term_count)
    isnothing(expected_term_count) && return nothing
    value = Int(expected_term_count)
    value > 0 ||
        throw(ArgumentError("PQS source-box IDA term count must be positive"))
    return value
end

function _pqs_source_box_ida_axis_factors(
    pgdg::_MappedOrdinaryPGDGIntermediate1D;
    axis::Symbol,
    expected_term_count = nothing,
)
    axis in (:x, :y, :z) ||
        throw(ArgumentError("PQS source-box IDA axis must be :x, :y, or :z"))
    density_terms = pgdg.pair_factor_terms
    raw_terms = pgdg.pair_factor_terms_raw
    size(density_terms) == size(raw_terms) || throw(
        DimensionMismatch("PQS source-box IDA raw and normalized terms differ"),
    )
    term_count, nrow, ncol = size(raw_terms)
    nrow == ncol ||
        throw(DimensionMismatch("PQS source-box IDA axis factors must be square"))
    expected = _pqs_source_box_expected_ida_term_count(expected_term_count)
    isnothing(expected) || term_count == expected || throw(
        DimensionMismatch("PQS source-box IDA term count mismatch"),
    )
    weights = Float64[Float64(weight) for weight in pgdg.weights]
    centers = Float64[Float64(center) for center in pgdg.centers]
    length(weights) == nrow ||
        throw(DimensionMismatch("PQS source-box IDA weights size mismatch"))
    length(centers) == nrow ||
        throw(DimensionMismatch("PQS source-box IDA centers size mismatch"))
    all(isfinite, raw_terms) ||
        throw(ArgumentError("PQS source-box IDA raw terms must be finite"))
    all(isfinite, weights) ||
        throw(ArgumentError("PQS source-box IDA weights must be finite"))
    all(isfinite, centers) ||
        throw(ArgumentError("PQS source-box IDA centers must be finite"))
    return (;
        axis,
        raw_pair_factor_terms = raw_terms,
        density_normalized_pair_factor_terms = density_terms,
        source_weights = weights,
        centers,
        term_count,
        factor_dimensions = (nrow, ncol),
    )
end

function _pqs_source_box_ida_axis_factors(
    bundle::_MappedOrdinaryGausslet1DBundle;
    axis::Symbol,
    expected_term_count = nothing,
)
    return _pqs_source_box_ida_axis_factors(
        bundle.pgdg_intermediate;
        axis,
        expected_term_count,
    )
end

function _pqs_source_box_ida_factor_provenance(
    bundles::_CartesianNestedAxisBundles3D;
    expected_term_count = nothing,
)
    axes = (;
        x = _pqs_source_box_ida_axis_factors(
            _nested_axis_bundle(bundles, :x);
            axis = :x,
            expected_term_count,
        ),
        y = _pqs_source_box_ida_axis_factors(
            _nested_axis_bundle(bundles, :y);
            axis = :y,
            expected_term_count,
        ),
        z = _pqs_source_box_ida_axis_factors(
            _nested_axis_bundle(bundles, :z);
            axis = :z,
            expected_term_count,
        ),
    )
    term_count = axes.x.term_count
    term_count == axes.y.term_count == axes.z.term_count || throw(
        DimensionMismatch("PQS source-box IDA axes have different term counts"),
    )
    return (;
        axes,
        raw_axis_pair_factor_terms = (;
            x = axes.x.raw_pair_factor_terms,
            y = axes.y.raw_pair_factor_terms,
            z = axes.z.raw_pair_factor_terms,
        ),
        axis_pair_factor_terms = (;
            x = axes.x.density_normalized_pair_factor_terms,
            y = axes.y.density_normalized_pair_factor_terms,
            z = axes.z.density_normalized_pair_factor_terms,
        ),
        axis_weights = (;
            x = axes.x.source_weights,
            y = axes.y.source_weights,
            z = axes.z.source_weights,
        ),
        axis_centers = (;
            x = axes.x.centers,
            y = axes.y.centers,
            z = axes.z.centers,
        ),
        term_count,
        factor_dimensions = (;
            x = axes.x.factor_dimensions,
            y = axes.y.factor_dimensions,
            z = axes.z.factor_dimensions,
        ),
    )
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

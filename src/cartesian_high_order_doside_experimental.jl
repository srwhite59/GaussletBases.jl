const _ExperimentalHighOrderCoefficientMap =
    Union{Matrix{Float64},SparseArrays.SparseMatrixCSC{Float64,Int}}

struct _ExperimentalHighOrderAxisData1D
    basis::MappedUniformBasis
    pgdg_intermediate
    one_body_layer
    backend::Symbol
    centers::Vector{Float64}
    weights::Vector{Float64}
    overlap::Matrix{Float64}
    position::Matrix{Float64}
    kinetic::Matrix{Float64}
    gaussian_factors::Vector{Matrix{Float64}}
    pair_factors_1d::Vector{Matrix{Float64}}
    one_body_cache::Dict{Tuple{Tuple{Vararg{Float64}},Float64},MappedOrdinaryOneBody1D}
end

struct _ExperimentalHighOrderBlock1D
    side::Int
    interval::UnitRange{Int}
    coefficients::_ExperimentalHighOrderCoefficientMap
    local_coefficients::Matrix{Float64}
    local_overlap::Matrix{Float64}
    local_position::Matrix{Float64}
    local_weights::Vector{Float64}
    local_centers::Vector{Float64}
    localized_centers::Vector{Float64}
    inner_axis_indices::Vector{Int}
    shell_axis_indices::Vector{Int}
end

struct _ExperimentalHighOrderTensorShell3D
    side::Int
    interval::NTuple{3,UnitRange{Int}}
    full_block_coefficients::_ExperimentalHighOrderCoefficientMap
    shell_coefficients::_ExperimentalHighOrderCoefficientMap
    shell_labels::Vector{NTuple{3,Int}}
    shell_kind_counts::NamedTuple{(:faces, :edges, :corners),NTuple{3,Int}}
    column_range::UnitRange{Int}
end

struct ExperimentalHighOrderDosideStack3D
    parent_basis::MappedUniformBasis
    backend::Symbol
    doside::Int
    sides::Vector{Int}
    parent_side::Int
    coefficient_matrix::_ExperimentalHighOrderCoefficientMap
    block_column_ranges::Vector{UnitRange{Int}}
    block_labels::Vector{Symbol}
    shell_layers::Vector{_ExperimentalHighOrderTensorShell3D}
    contracted_weights::Vector{Float64}
    diagnostics::NamedTuple
end

function _experimental_high_order_one_body_cache_key(
    exponents::AbstractVector{<:Real},
    center::Real,
)
    exponent_values = Float64[Float64(exponent) for exponent in exponents]
    center_value = Float64(center)
    return exponent_values, center_value, (Tuple(exponent_values), center_value)
end

function _experimental_high_order_one_body_from_pgdg_intermediate(
    pgdg_intermediate::_MappedOrdinaryPGDGIntermediate1D,
)
    return MappedOrdinaryOneBody1D(
        pgdg_intermediate.basis,
        pgdg_intermediate.backend,
        Matrix{Float64}(pgdg_intermediate.overlap),
        Matrix{Float64}(pgdg_intermediate.kinetic),
        Matrix{Float64}[Matrix{Float64}(factor) for factor in pgdg_intermediate.gaussian_factors],
        copy(pgdg_intermediate.exponents),
        Float64(pgdg_intermediate.center),
    )
end

function _experimental_high_order_centered_interval(
    n1d::Int,
    side::Int,
)
    side >= 1 || throw(ArgumentError("experimental high-order doside side lengths must be positive"))
    side <= n1d || throw(ArgumentError("experimental high-order doside side length must lie inside the parent basis"))
    start = (n1d - side) ÷ 2 + 1
    return start:(start + side - 1)
end

function _experimental_high_order_mapping_family(mapping_value::IdentityMapping)
    return :identity
end

function _experimental_high_order_mapping_family(mapping_value::AsinhMapping)
    c_value = mapping_value.a * mapping_value.s
    if isapprox(c_value, 0.2; atol = 1.0e-12, rtol = 1.0e-12) &&
       isapprox(mapping_value.s, sqrt(0.4); atol = 1.0e-12, rtol = 1.0e-12) &&
       isapprox(mapping_value.tail_spacing, 10.0; atol = 1.0e-12, rtol = 1.0e-12)
        return :white_lindsey_atomic_he_d0p2
    end
    throw(
        ArgumentError(
            "experimental high-order doside stack currently supports only IdentityMapping or the explicit distorted White-Lindsey He family AsinhMapping(c = 0.2, s = sqrt(0.4), tail_spacing = 10.0)",
        ),
    )
end

function _experimental_high_order_mapping_family(basis::MappedUniformBasis)
    return _experimental_high_order_mapping_family(mapping(basis))
end

function _experimental_high_order_centered_working_box(
    parent_side::Int,
    working_box_side::Int,
)
    interval = _experimental_high_order_centered_interval(parent_side, working_box_side)
    return (interval, interval, interval)
end

function _experimental_high_order_validate_request(
    basis::MappedUniformBasis,
    sides::AbstractVector{<:Integer},
    doside::Int,
)
    mapping_family = _experimental_high_order_mapping_family(basis)
    doside in (4, 5, 6) || throw(
        ArgumentError("experimental high-order doside stack currently supports only doside = 4, 5, or 6"),
    )
    n1d = length(basis)
    isodd(n1d) || throw(ArgumentError("experimental high-order doside stack currently requires an odd parent side length"))
    isempty(sides) && throw(ArgumentError("experimental high-order doside stack requires at least one side length"))
    side_values = Int[Int(side) for side in sides]
    first(side_values) == doside || throw(ArgumentError("experimental high-order doside side ladders must start at doside"))
    all(side -> isodd(side) == isodd(doside), side_values) || throw(
        ArgumentError("experimental high-order doside side ladders must keep the same parity as doside"),
    )
    all(side_values[index + 1] == side_values[index] + 2 for index in 1:(length(side_values) - 1)) || throw(
        ArgumentError("experimental high-order doside side ladders must increase by 2"),
    )
    maximum(side_values) <= n1d || throw(
        ArgumentError("experimental high-order doside stack currently requires maximum(sides) <= length(parent_basis)"),
    )
    return side_values, mapping_family
end

function _experimental_high_order_axis_data_1d(
    basis::MappedUniformBasis;
    backend::Symbol = :numerical_reference,
    prepared_bundle::Union{Nothing,_MappedOrdinaryGausslet1DBundle} = nothing,
    one_body_exponents::AbstractVector{<:Real} = Float64[],
    one_body_center::Real = 0.0,
)
    if backend == :numerical_reference
        representation = basis_representation(basis; operators = (:overlap, :position, :kinetic))
        return _ExperimentalHighOrderAxisData1D(
            basis,
            nothing,
            basis,
            backend,
            Float64[Float64(value) for value in centers(basis)],
            Float64[Float64(value) for value in integral_weights(basis)],
            Matrix{Float64}(representation.basis_matrices.overlap),
            Matrix{Float64}(representation.basis_matrices.position),
            Matrix{Float64}(representation.basis_matrices.kinetic),
            Matrix{Float64}[],
            Matrix{Float64}[],
            Dict{Tuple{Tuple{Vararg{Float64}},Float64},MappedOrdinaryOneBody1D}(),
        )
    end

    cache = Dict{Tuple{Tuple{Vararg{Float64}},Float64},MappedOrdinaryOneBody1D}()
    exponent_values, center_value, key = _experimental_high_order_one_body_cache_key(one_body_exponents, one_body_center)

    pgdg = if isnothing(prepared_bundle)
        pgdg_exponents = isempty(exponent_values) ? [1.0] : exponent_values
        pgdg_intermediate = _mapped_ordinary_pgdg_intermediate_1d(
            basis;
            exponents = pgdg_exponents,
            center = center_value,
            backend = backend,
        )
        if !isempty(exponent_values)
            cache[key] = _experimental_high_order_one_body_from_pgdg_intermediate(pgdg_intermediate)
        end
        pgdg_intermediate
    else
        prepared_bundle.backend == backend || throw(
            ArgumentError(
                "prepared high-order axis bundle backend $(prepared_bundle.backend) does not match requested backend $backend",
            ),
        )
        prepared_bundle.basis === basis || throw(
            ArgumentError("prepared high-order axis bundle must come from the same parent basis object"),
        )
        prepared = prepared_bundle.pgdg_intermediate
        if prepared_bundle.exponents == exponent_values &&
           isapprox(prepared_bundle.center, center_value; atol = 0.0, rtol = 0.0)
            cache[key] = _mapped_ordinary_one_body_from_bundle(prepared_bundle)
        end
        prepared
    end
    return _ExperimentalHighOrderAxisData1D(
        basis,
        pgdg,
        pgdg.auxiliary_layer,
        backend,
        Float64[Float64(value) for value in pgdg.centers],
        Float64[Float64(value) for value in pgdg.weights],
        Matrix{Float64}(pgdg.overlap),
        Matrix{Float64}(pgdg.position),
        Matrix{Float64}(pgdg.kinetic),
        Matrix{Float64}[Matrix{Float64}(factor) for factor in pgdg.gaussian_factors],
        Matrix{Float64}[Matrix{Float64}(factor) for factor in pgdg.pair_factors],
        cache,
    )
end

function _experimental_high_order_axis_one_body_1d(
    axis_data::_ExperimentalHighOrderAxisData1D;
    exponents::AbstractVector{<:Real} = Float64[],
    center::Real = 0.0,
)
    exponent_values, center_value, key = _experimental_high_order_one_body_cache_key(exponents, center)
    if haskey(axis_data.one_body_cache, key)
        return axis_data.one_body_cache[key]
    end

    gaussian_factors = Matrix{Float64}[
        Matrix{Float64}(factor) for factor in gaussian_factor_matrices(
            axis_data.one_body_layer;
            exponents = exponent_values,
            center = center_value,
        )
    ]

    one_body = MappedOrdinaryOneBody1D(
        axis_data.basis,
        axis_data.backend,
        axis_data.overlap,
        axis_data.kinetic,
        gaussian_factors,
        exponent_values,
        center_value,
    )
    axis_data.one_body_cache[key] = one_body
    return one_body
end

function _experimental_high_order_interval_data(
    axis_data::_ExperimentalHighOrderAxisData1D,
    interval::UnitRange{Int},
)
    if !isnothing(axis_data.pgdg_intermediate)
        return _nested_interval_data(axis_data.pgdg_intermediate, interval)
    end
    return (
        overlap = Matrix{Float64}(axis_data.overlap[interval, interval]),
        position = Matrix{Float64}(axis_data.position[interval, interval]),
        weights = Float64[axis_data.weights[index] for index in interval],
        centers = Float64[axis_data.centers[index] for index in interval],
        n1d = length(axis_data.centers),
    )
end

function _experimental_high_order_embed_local_coefficients(
    local_coefficients::AbstractMatrix{<:Real},
    interval::UnitRange{Int},
    n1d::Int,
)
    row_indices = Int[]
    col_indices = Int[]
    values = Float64[]
    for column in axes(local_coefficients, 2), (local_row, full_row) in enumerate(interval)
        value = Float64(local_coefficients[local_row, column])
        iszero(value) && continue
        push!(row_indices, full_row)
        push!(col_indices, column)
        push!(values, value)
    end
    return _nested_sparse_coefficient_map(
        row_indices,
        col_indices,
        values,
        n1d,
        size(local_coefficients, 2),
    )
end

function _experimental_high_order_parent_polynomial_images(
    local_weights::AbstractVector{<:Real},
    local_position::AbstractMatrix{<:Real},
    local_overlap::AbstractMatrix{<:Real},
    retained_count::Int,
)
    nlocal = length(local_weights)
    retained_count >= 1 || throw(
        ArgumentError("parent polynomial images require retained_count >= 1"),
    )
    retained_count <= nlocal || throw(
        ArgumentError("parent polynomial images require retained_count <= local interval size"),
    )
    overlap_factor = cholesky(Symmetric(Matrix{Float64}(local_overlap)))
    position_value = Matrix{Float64}(local_position)
    images = zeros(Float64, nlocal, retained_count)
    images[:, 1] .= overlap_factor \ Float64.(local_weights)
    for column in 2:retained_count
        previous = view(images, :, column - 1)
        images[:, column] .= overlap_factor \ (position_value * previous)
    end
    return images
end

function _experimental_high_order_metric_projection_summary(
    targets::AbstractMatrix{<:Real},
    basis::AbstractMatrix{<:Real},
    overlap::AbstractMatrix{<:Real},
)
    target_value = Matrix{Float64}(targets)
    basis_value = Matrix{Float64}(basis)
    couplings = Matrix{Float64}(transpose(basis_value) * overlap * target_value)
    projected = Matrix{Float64}(basis_value * couplings)
    residual = Matrix{Float64}(target_value - projected)
    errors = Float64[]
    capture_fractions = Float64[]
    for column in axes(target_value, 2)
        target_norm = _nested_metric_norm(view(target_value, :, column), overlap)
        projected_norm = _nested_metric_norm(view(projected, :, column), overlap)
        residual_norm = _nested_metric_norm(view(residual, :, column), overlap)
        error_value = residual_norm / max(target_norm, eps(Float64))
        capture_value = (projected_norm / max(target_norm, eps(Float64)))^2
        push!(errors, error_value)
        push!(capture_fractions, min(1.0, max(0.0, capture_value)))
    end
    return (
        errors = errors,
        capture_fractions = capture_fractions,
        max_error = isempty(errors) ? 0.0 : maximum(errors),
        min_capture_fraction = isempty(capture_fractions) ? 1.0 : minimum(capture_fractions),
    )
end

function _experimental_high_order_subspace_residual_error(
    target::AbstractMatrix{<:Real},
    basis::AbstractMatrix{<:Real},
    overlap::AbstractMatrix{<:Real},
)
    target_value = Matrix{Float64}(target)
    basis_value = Matrix{Float64}(basis)
    projected = Matrix{Float64}(basis_value * (transpose(basis_value) * overlap * target_value))
    residual = Matrix{Float64}(target_value - projected)
    return norm(transpose(residual) * overlap * residual, Inf)
end

function _experimental_high_order_axis_shell_indices(axis_count::Int)
    axis_count >= 1 || throw(ArgumentError("axis shell indices require axis_count >= 1"))
    if axis_count == 1
        return Int[]
    elseif axis_count == 2
        return Int[1, 2]
    end
    return collect(2:(axis_count - 1))
end

function _experimental_high_order_axis_boundary_indices(axis_count::Int)
    axis_count >= 1 || throw(ArgumentError("axis boundary indices require axis_count >= 1"))
    indices = Int[1]
    axis_count > 1 && push!(indices, axis_count)
    return unique(indices)
end

function _experimental_high_order_physical_block_1d(
    axis_data::_ExperimentalHighOrderAxisData1D,
    side::Int;
    doside::Int = 5,
)
    interval = _experimental_high_order_centered_interval(length(axis_data.centers), side)
    interval_data = _experimental_high_order_interval_data(axis_data, interval)
    parent_polynomial_images = _experimental_high_order_parent_polynomial_images(
        interval_data.weights,
        interval_data.position,
        interval_data.overlap,
        doside,
    )
    raw_basis = _nested_retained_span(
        interval_data.weights,
        interval_data.centers,
        interval_data.position,
        interval_data.overlap,
        doside,
    )
    overlap_seed = Matrix{Float64}(transpose(raw_basis) * interval_data.overlap * raw_basis)
    position_seed = Matrix{Float64}(transpose(raw_basis) * interval_data.position * raw_basis)
    sign_vector = vec(transpose(interval_data.weights) * raw_basis)
    transform, localized_centers = _cleanup_comx_transform(overlap_seed, position_seed, sign_vector)
    local_coefficients = Matrix{Float64}(raw_basis * transform)
    local_weights = vec(transpose(interval_data.weights) * local_coefficients)
    coefficients = _experimental_high_order_embed_local_coefficients(
        local_coefficients,
        interval,
        interval_data.n1d,
    )
    axis_count = size(local_coefficients, 2)
    block = _ExperimentalHighOrderBlock1D(
        side,
        interval,
        coefficients,
        local_coefficients,
        interval_data.overlap,
        interval_data.position,
        local_weights,
        interval_data.centers,
        localized_centers,
        _experimental_high_order_axis_shell_indices(axis_count),
        _experimental_high_order_axis_boundary_indices(axis_count),
    )
    polynomial_capture = _experimental_high_order_metric_projection_summary(
        parent_polynomial_images,
        local_coefficients,
        interval_data.overlap,
    )
    diagnostics = (
        polynomial_degrees = collect(0:(size(parent_polynomial_images, 2) - 1)),
        precleanup_overlap_spectrum = _experimental_high_order_positive_spectrum(overlap_seed; tol = 1.0e-12),
        overlap_error = _experimental_high_order_overlap_error(local_coefficients, interval_data.overlap),
        parent_polynomial_projection_errors = polynomial_capture.errors,
        parent_polynomial_capture_fractions = polynomial_capture.capture_fractions,
        max_parent_polynomial_projection_error = polynomial_capture.max_error,
        min_parent_polynomial_capture_fraction = polynomial_capture.min_capture_fraction,
    )
    return (
        block = block,
        raw_basis = raw_basis,
        parent_polynomial_images = parent_polynomial_images,
        diagnostics = diagnostics,
    )
end

function _experimental_high_order_block_1d(
    axis_data::_ExperimentalHighOrderAxisData1D,
    side::Int;
    doside::Int = 5,
)
    # Compatibility/debug helper retaining the earlier center-based local block
    # construction. The physical-coordinate block helper is the intended route
    # contract for high-order production-like paths.
    interval = _experimental_high_order_centered_interval(length(axis_data.centers), side)
    interval_data = _experimental_high_order_interval_data(axis_data, interval)
    raw_basis = _nested_retained_span(
        interval_data.weights,
        interval_data.centers,
        interval_data.position,
        interval_data.overlap,
        doside,
    )
    overlap_seed = Matrix{Float64}(transpose(raw_basis) * interval_data.overlap * raw_basis)
    position_seed = Matrix{Float64}(transpose(raw_basis) * interval_data.position * raw_basis)
    sign_vector = vec(transpose(raw_basis) * interval_data.weights)
    transform, localized_centers = _cleanup_comx_transform(overlap_seed, position_seed, sign_vector)
    local_coefficients = Matrix{Float64}(raw_basis * transform)
    local_weights = vec(transpose(interval_data.weights) * local_coefficients)
    coefficients = _experimental_high_order_embed_local_coefficients(
        local_coefficients,
        interval,
        interval_data.n1d,
    )
    return _ExperimentalHighOrderBlock1D(
        side,
        interval,
        coefficients,
        local_coefficients,
        interval_data.overlap,
        interval_data.position,
        local_weights,
        interval_data.centers,
        localized_centers,
        collect(2:(doside - 1)),
        Int[1, doside],
    )
end

function _experimental_high_order_product_coefficients(
    x_block::_ExperimentalHighOrderBlock1D,
    y_block::_ExperimentalHighOrderBlock1D,
    z_block::_ExperimentalHighOrderBlock1D,
    dims::NTuple{3,Int},
)
    nparent = prod(dims)
    row_indices = Int[]
    col_indices = Int[]
    values = Float64[]
    labels = NTuple{3,Int}[]
    column = 0
    for ixcol in axes(x_block.local_coefficients, 2),
        iycol in axes(y_block.local_coefficients, 2),
        izcol in axes(z_block.local_coefficients, 2)
        column += 1
        push!(labels, (ixcol, iycol, izcol))
        for (ix_local, ix) in enumerate(x_block.interval)
            vx = Float64(x_block.local_coefficients[ix_local, ixcol])
            iszero(vx) && continue
            for (iy_local, iy) in enumerate(y_block.interval)
                vy = Float64(y_block.local_coefficients[iy_local, iycol])
                iszero(vy) && continue
                for (iz_local, iz) in enumerate(z_block.interval)
                    vz = Float64(z_block.local_coefficients[iz_local, izcol])
                    iszero(vz) && continue
                    push!(row_indices, _cartesian_flat_index(ix, iy, iz, dims))
                    push!(col_indices, column)
                    push!(values, vx * vy * vz)
                end
            end
        end
    end
    coefficients = _nested_sparse_coefficient_map(
        row_indices,
        col_indices,
        values,
        nparent,
        length(labels),
    )
    return coefficients, labels
end

function _experimental_high_order_shell_kind_counts(
    labels::AbstractVector{<:NTuple{3,Int}};
    doside::Int = 5,
)
    faces = 0
    edges = 0
    corners = 0
    for label in labels
        boundary_axes = count(index -> index == 1 || index == doside, label)
        if boundary_axes == 1
            faces += 1
        elseif boundary_axes == 2
            edges += 1
        elseif boundary_axes == 3
            corners += 1
        end
    end
    return (faces = faces, edges = edges, corners = corners)
end

function _experimental_high_order_select_shell_indices_and_labels(
    labels::AbstractVector{<:NTuple{3,Int}};
    doside::Int = 5,
)
    shell_indices = Int[]
    shell_labels = NTuple{3,Int}[]
    for (column, label) in enumerate(labels)
        if any(index -> index == 1 || index == doside, label)
            push!(shell_indices, column)
            push!(shell_labels, label)
        end
    end
    return shell_indices, shell_labels
end

function _experimental_high_order_tensor_shell_indices_and_labels(
    doside::Int,
)
    labels = NTuple{3,Int}[]
    for ix in 1:doside, iy in 1:doside, iz in 1:doside
        push!(labels, (ix, iy, iz))
    end
    return _experimental_high_order_select_shell_indices_and_labels(labels; doside = doside)
end

function _experimental_high_order_tensor_shell_3d(
    block::_ExperimentalHighOrderBlock1D,
    parent_side::Int;
    doside::Int = 5,
    column_range::UnitRange{Int} = 1:0,
)
    full_block_coefficients, labels = _experimental_high_order_product_coefficients(
        block,
        block,
        block,
        (parent_side, parent_side, parent_side),
    )

    shell_indices, shell_labels = _experimental_high_order_select_shell_indices_and_labels(
        labels;
        doside = doside,
    )

    shell_coefficients = full_block_coefficients[:, shell_indices]
    return _ExperimentalHighOrderTensorShell3D(
        block.side,
        (block.interval, block.interval, block.interval),
        full_block_coefficients,
        shell_coefficients,
        shell_labels,
        _experimental_high_order_shell_kind_counts(shell_labels; doside = doside),
        column_range,
    )
end

function _experimental_high_order_tensor_shell_3d(
    axis_data::_ExperimentalHighOrderAxisData1D,
    side::Int;
    doside::Int = 5,
    column_range::UnitRange{Int} = 1:0,
)
    block = _experimental_high_order_block_1d(axis_data, side; doside = doside)
    return _experimental_high_order_tensor_shell_3d(
        block,
        length(axis_data.centers);
        doside = doside,
        column_range = column_range,
    )
end

function _experimental_high_order_physical_full_block_3d(
    axis_data::_ExperimentalHighOrderAxisData1D,
    side::Int;
    doside::Int = 5,
)
    parent_side = length(axis_data.centers)
    parent_overlap = _experimental_high_order_parent_overlap_3d(axis_data)
    physical_1d = _experimental_high_order_physical_block_1d(axis_data, side; doside = doside)
    physical_shell = _experimental_high_order_tensor_shell_3d(
        physical_1d.block,
        parent_side;
        doside = doside,
    )
    current_shell = _experimental_high_order_tensor_shell_3d(
        axis_data,
        side;
        doside = doside,
    )
    physical_full = Matrix{Float64}(physical_shell.full_block_coefficients)
    current_full = Matrix{Float64}(current_shell.full_block_coefficients)
    diagnostics = (
        full_block_dimension = size(physical_full, 2),
        full_block_overlap_error = _experimental_high_order_overlap_error(physical_full, parent_overlap),
        full_block_overlap_spectrum = _experimental_high_order_metric_spectrum(
            physical_full,
            parent_overlap;
            tol = 1.0e-12,
        ),
        current_in_physical_residual = _experimental_high_order_subspace_residual_error(
            current_full,
            physical_full,
            parent_overlap,
        ),
        physical_in_current_residual = _experimental_high_order_subspace_residual_error(
            physical_full,
            current_full,
            parent_overlap,
        ),
    )
    return (
        block_1d = physical_1d,
        shell = physical_shell,
        current_shell = current_shell,
        diagnostics = diagnostics,
    )
end

function _experimental_high_order_expected_shell_dimension(doside::Int)
    doside >= 1 || throw(ArgumentError("expected shell dimension requires doside >= 1"))
    inner = max(doside - 2, 0)
    return doside^3 - inner^3
end

function _experimental_high_order_physical_shell_3d(
    axis_data::_ExperimentalHighOrderAxisData1D,
    side::Int;
    doside::Int = 5,
)
    physical_full = _experimental_high_order_physical_full_block_3d(
        axis_data,
        side;
        doside = doside,
    )
    parent_overlap = _experimental_high_order_parent_overlap_3d(axis_data)
    physical_shell_coefficients = Matrix{Float64}(physical_full.shell.shell_coefficients)
    current_shell_coefficients = Matrix{Float64}(physical_full.current_shell.shell_coefficients)
    diagnostics = (
        shell_dimension = size(physical_shell_coefficients, 2),
        expected_shell_dimension = _experimental_high_order_expected_shell_dimension(doside),
        shell_overlap_error = _experimental_high_order_overlap_error(
            physical_shell_coefficients,
            parent_overlap,
        ),
        shell_overlap_spectrum = _experimental_high_order_metric_spectrum(
            physical_shell_coefficients,
            parent_overlap;
            tol = 1.0e-12,
        ),
        shell_kind_counts = physical_full.shell.shell_kind_counts,
        current_shell_in_physical_residual = _experimental_high_order_subspace_residual_error(
            current_shell_coefficients,
            physical_shell_coefficients,
            parent_overlap,
        ),
        physical_shell_in_current_residual = _experimental_high_order_subspace_residual_error(
            physical_shell_coefficients,
            current_shell_coefficients,
            parent_overlap,
        ),
    )
    return (
        full_block = physical_full,
        shell = physical_full.shell,
        current_shell = physical_full.current_shell,
        diagnostics = diagnostics,
    )
end

function _experimental_high_order_parent_overlap_3d(
    axis_data::_ExperimentalHighOrderAxisData1D,
)
    overlap = axis_data.overlap
    return Matrix{Float64}(kron(overlap, kron(overlap, overlap)))
end

function _experimental_high_order_parent_weights_3d(
    axis_data::_ExperimentalHighOrderAxisData1D,
)
    return _mapped_cartesian_weights(axis_data.weights)
end

function _experimental_high_order_parent_coordinates_3d(
    axis_data::_ExperimentalHighOrderAxisData1D,
)
    orbitals = _mapped_cartesian_orbitals(axis_data.centers)
    return (
        x = Float64[orbital.x for orbital in orbitals],
        y = Float64[orbital.y for orbital in orbitals],
        z = Float64[orbital.z for orbital in orbitals],
    )
end

function _experimental_high_order_sign_fix_columns!(
    coefficients::Matrix{Float64},
    sign_vector::AbstractVector{<:Real},
)
    for column in axes(coefficients, 2)
        marker = dot(sign_vector, view(coefficients, :, column))
        if abs(marker) <= 1.0e-12
            local_view = @view coefficients[:, column]
            pivot = argmax(abs.(local_view))
            marker = local_view[pivot]
        end
        if marker < 0.0
            coefficients[:, column] .*= -1.0
        end
    end
    return coefficients
end

function _experimental_high_order_metric_project_out(
    candidates::AbstractMatrix{<:Real},
    basis::AbstractMatrix{<:Real},
    overlap::AbstractMatrix{<:Real},
)
    basis_columns = size(basis, 2)
    basis_columns == 0 && return Matrix{Float64}(candidates)
    basis_value = Matrix{Float64}(basis)
    candidate_value = Matrix{Float64}(candidates)
    couplings = Matrix{Float64}(transpose(basis_value) * overlap * candidate_value)
    return Matrix{Float64}(candidate_value - basis_value * couplings)
end

function _experimental_high_order_lowdin_cleanup(
    coefficients::AbstractMatrix{<:Real},
    overlap::AbstractMatrix{<:Real};
    sign_vector::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    overlap_seed = Matrix{Float64}(transpose(coefficients) * overlap * coefficients)
    vectors, invhalf = _s_invsqrt_reduced(overlap_seed)
    cleaned = Matrix{Float64}(coefficients * (vectors * invhalf))
    if !isnothing(sign_vector)
        _experimental_high_order_sign_fix_columns!(cleaned, sign_vector)
    end
    return cleaned
end

function _experimental_high_order_overlap_error(
    coefficients::AbstractMatrix{<:Real},
    overlap::AbstractMatrix{<:Real},
)
    gram = Matrix{Float64}(transpose(coefficients) * overlap * coefficients)
    return norm(gram - I, Inf)
end

function _experimental_high_order_positive_spectrum(
    matrix::AbstractMatrix{<:Real};
    tol::Real = 1.0e-10,
)
    size(matrix, 1) == size(matrix, 2) || throw(ArgumentError("positive-spectrum helper requires a square matrix"))
    tol_value = Float64(tol)
    tol_value > 0.0 || throw(ArgumentError("positive-spectrum helper requires tol > 0"))
    values = Float64[
        Float64(real(value)) for value in eigvals(Symmetric(_symmetrize_ida_matrix(Matrix{Float64}(matrix))))
    ]
    kept = Float64[value for value in values if value > tol_value]
    if isempty(kept)
        return (
            minimum_eigenvalue = 0.0,
            maximum_eigenvalue = 0.0,
            condition_number = Inf,
            kept_rank = 0,
            dropped_rank = length(values),
        )
    end
    return (
        minimum_eigenvalue = minimum(kept),
        maximum_eigenvalue = maximum(kept),
        condition_number = maximum(kept) / minimum(kept),
        kept_rank = length(kept),
        dropped_rank = length(values) - length(kept),
    )
end

function _experimental_high_order_metric_spectrum(
    coefficients::AbstractMatrix{<:Real},
    overlap::AbstractMatrix{<:Real};
    tol::Real = 1.0e-10,
)
    gram = Matrix{Float64}(transpose(coefficients) * overlap * coefficients)
    return _experimental_high_order_positive_spectrum(gram; tol = tol)
end

function _experimental_high_order_axis_moment_summary(
    center::Float64,
    coordinates::AbstractVector{<:Real},
    signed_weights::AbstractVector{<:Real},
    absolute_weights::AbstractVector{<:Real},
    absolute_mass::Real,
)
    m1_signed = 0.0
    m2_signed = 0.0
    m3_signed = 0.0
    m4_signed = 0.0
    m1_absolute = 0.0
    m2_absolute = 0.0
    m3_absolute = 0.0
    m4_absolute = 0.0
    for index in eachindex(coordinates)
        dx = Float64(coordinates[index]) - center
        dx2 = dx * dx
        dx3 = dx2 * dx
        dx4 = dx2 * dx2
        signed_weight = Float64(signed_weights[index])
        absolute_weight = Float64(absolute_weights[index])
        adx = abs(dx)
        m1_signed += signed_weight * dx
        m2_signed += signed_weight * dx2
        m3_signed += signed_weight * dx3
        m4_signed += signed_weight * dx4
        m1_absolute += absolute_weight * adx
        m2_absolute += absolute_weight * dx2
        m3_absolute += absolute_weight * adx * dx2
        m4_absolute += absolute_weight * dx4
    end
    norm_value = Float64(absolute_mass)
    sigma = sqrt(max(m2_absolute / norm_value, eps(Float64)))
    normalized_mu3 = abs(m3_signed / norm_value) / sigma^3
    normalized_mu4 = abs(m4_absolute / norm_value) / sigma^4
    return (
        sigma = sigma,
        signed = (
            degree1 = m1_signed / norm_value,
            degree2 = m2_signed / norm_value,
            degree3 = m3_signed / norm_value,
            degree4 = m4_signed / norm_value,
        ),
        absolute = (
            degree1 = m1_absolute / norm_value,
            degree2 = m2_absolute / norm_value,
            degree3 = m3_absolute / norm_value,
            degree4 = m4_absolute / norm_value,
        ),
        normalized = (
            mu3 = normalized_mu3,
            mu4 = normalized_mu4,
        ),
    )
end

function _experimental_high_order_column_moment_diagnostic(
    column::AbstractVector{<:Real},
    parent_weights::AbstractVector{<:Real},
    parent_coordinates::NamedTuple;
    column_index::Int,
    block_index::Int,
    block_label::Symbol,
)
    length(column) == length(parent_weights) || throw(
        DimensionMismatch("column moment diagnostic requires parent weights that match the column length"),
    )
    length(parent_coordinates.x) == length(parent_weights) ||
        throw(DimensionMismatch("parent x coordinates must match the parent weight length"))
    length(parent_coordinates.y) == length(parent_weights) ||
        throw(DimensionMismatch("parent y coordinates must match the parent weight length"))
    length(parent_coordinates.z) == length(parent_weights) ||
        throw(DimensionMismatch("parent z coordinates must match the parent weight length"))

    signed_weights = Float64[Float64(parent_weights[index]) * Float64(column[index]) for index in eachindex(column)]
    absolute_weights = abs.(signed_weights)
    absolute_mass = sum(absolute_weights)
    absolute_mass > 0.0 || throw(ArgumentError("column moment diagnostic requires nonzero absolute parent weight"))

    center_x = dot(absolute_weights, parent_coordinates.x) / absolute_mass
    center_y = dot(absolute_weights, parent_coordinates.y) / absolute_mass
    center_z = dot(absolute_weights, parent_coordinates.z) / absolute_mass

    x_moments = _experimental_high_order_axis_moment_summary(
        center_x,
        parent_coordinates.x,
        signed_weights,
        absolute_weights,
        absolute_mass,
    )
    y_moments = _experimental_high_order_axis_moment_summary(
        center_y,
        parent_coordinates.y,
        signed_weights,
        absolute_weights,
        absolute_mass,
    )
    z_moments = _experimental_high_order_axis_moment_summary(
        center_z,
        parent_coordinates.z,
        signed_weights,
        absolute_weights,
        absolute_mass,
    )

    center_drift = sqrt(center_x^2 + center_y^2 + center_z^2)
    max_normalized_mu3 = maximum(
        (x_moments.normalized.mu3, y_moments.normalized.mu3, z_moments.normalized.mu3),
    )
    max_normalized_mu4 = maximum(
        (x_moments.normalized.mu4, y_moments.normalized.mu4, z_moments.normalized.mu4),
    )
    return (
        column_index = column_index,
        block_index = block_index,
        block_label = block_label,
        is_outer_shell = block_index > 1,
        center = (x = center_x, y = center_y, z = center_z),
        center_drift = center_drift,
        signed_centered_moments = (
            x = x_moments.signed,
            y = y_moments.signed,
            z = z_moments.signed,
        ),
        absolute_centered_moments = (
            x = x_moments.absolute,
            y = y_moments.absolute,
            z = z_moments.absolute,
        ),
        sigma = (
            x = x_moments.sigma,
            y = y_moments.sigma,
            z = z_moments.sigma,
        ),
        normalized_moment_ratios = (
            mu3 = (
                x = x_moments.normalized.mu3,
                y = y_moments.normalized.mu3,
                z = z_moments.normalized.mu3,
            ),
            mu4 = (
                x = x_moments.normalized.mu4,
                y = y_moments.normalized.mu4,
                z = z_moments.normalized.mu4,
            ),
        ),
        max_normalized_mu3 = max_normalized_mu3,
        max_normalized_mu4 = max_normalized_mu4,
        risk_score = max(max_normalized_mu3, max_normalized_mu4),
    )
end

function _experimental_high_order_stack_moment_risk(
    coefficients::AbstractMatrix{<:Real},
    parent_weights::AbstractVector{<:Real},
    parent_coordinates::NamedTuple,
    block_labels::AbstractVector{<:Symbol},
    block_column_ranges::AbstractVector{<:UnitRange{Int}},
)
    block_index_by_column = zeros(Int, size(coefficients, 2))
    block_label_by_column = fill(Symbol(""), size(coefficients, 2))
    for (block_index, column_range) in enumerate(block_column_ranges)
        for column in column_range
            block_index_by_column[column] = block_index
            block_label_by_column[column] = block_labels[block_index]
        end
    end

    column_diagnostics = NamedTuple[]
    for column in axes(coefficients, 2)
        push!(
            column_diagnostics,
            _experimental_high_order_column_moment_diagnostic(
                view(coefficients, :, column),
                parent_weights,
                parent_coordinates;
                column_index = column,
                block_index = block_index_by_column[column],
                block_label = block_label_by_column[column],
            ),
        )
    end

    outer_diagnostics = NamedTuple[
        diagnostic for diagnostic in column_diagnostics if diagnostic.is_outer_shell
    ]
    sorted_outer = sort(
        outer_diagnostics;
        by = diagnostic -> (-diagnostic.risk_score, -diagnostic.max_normalized_mu4, -diagnostic.max_normalized_mu3, diagnostic.column_index),
    )
    top_outer = first(sorted_outer, min(10, length(sorted_outer)))
    worst_outer = isempty(sorted_outer) ? nothing : first(sorted_outer)
    return (
        column_diagnostics = column_diagnostics,
        outer_column_count = length(sorted_outer),
        worst_outer = worst_outer,
        top_outer = top_outer,
    )
end

function _experimental_high_order_full_block_union_coefficients(
    axis_data::_ExperimentalHighOrderAxisData1D,
    sides::AbstractVector{<:Integer};
    doside::Int = 5,
)
    blocks = _ExperimentalHighOrderCoefficientMap[]
    for side in sides
        shell = _experimental_high_order_tensor_shell_3d(axis_data, Int(side); doside = doside)
        push!(blocks, shell.full_block_coefficients)
    end
    return _nested_hcat_coefficient_maps(blocks)
end

function _experimental_high_order_doside_stack_3d(
    basis::MappedUniformBasis;
    axis_data::Union{Nothing,_ExperimentalHighOrderAxisData1D} = nothing,
    backend::Symbol = :numerical_reference,
    doside::Int = 5,
    sides::AbstractVector{<:Integer} = [5, 7, 9, 11],
)
    side_values, mapping_family = _experimental_high_order_validate_request(basis, sides, doside)
    axis_data_value = isnothing(axis_data) ? _experimental_high_order_axis_data_1d(basis; backend = backend) : axis_data
    parent_overlap = _experimental_high_order_parent_overlap_3d(axis_data_value)
    parent_weights = _experimental_high_order_parent_weights_3d(axis_data_value)
    parent_coordinates = _experimental_high_order_parent_coordinates_3d(axis_data_value)

    first_block = _experimental_high_order_tensor_shell_3d(axis_data_value, first(side_values); doside = doside)
    accumulated = Matrix{Float64}(first_block.full_block_coefficients)
    _experimental_high_order_sign_fix_columns!(accumulated, parent_weights)

    block_column_ranges = UnitRange{Int}[1:size(accumulated, 2)]
    block_labels = Symbol[Symbol("side$(first(side_values))_full")]
    shell_layers = _ExperimentalHighOrderTensorShell3D[]
    shell_cleanup_spectra = NamedTuple[]
    next_column = size(accumulated, 2) + 1

    for side in side_values[2:end]
        shell = _experimental_high_order_tensor_shell_3d(axis_data_value, side; doside = doside)
        shell_residual = _experimental_high_order_metric_project_out(
            shell.shell_coefficients,
            accumulated,
            parent_overlap,
        )
        shell_spectrum = _experimental_high_order_metric_spectrum(
            shell_residual,
            parent_overlap;
            tol = 1.0e-10,
        )
        shell_clean = _experimental_high_order_lowdin_cleanup(
            shell_residual,
            parent_overlap;
            sign_vector = parent_weights,
        )
        shell_range = next_column:(next_column + size(shell_clean, 2) - 1)
        next_column = last(shell_range) + 1

        push!(
            shell_layers,
            _ExperimentalHighOrderTensorShell3D(
                shell.side,
                shell.interval,
                shell.full_block_coefficients,
                shell.shell_coefficients,
                shell.shell_labels,
                shell.shell_kind_counts,
                shell_range,
            ),
        )
        push!(
            shell_cleanup_spectra,
            (
                side = side,
                raw_dimension = size(shell_residual, 2),
                kept_rank = shell_spectrum.kept_rank,
                dropped_rank = shell_spectrum.dropped_rank,
                minimum_eigenvalue = shell_spectrum.minimum_eigenvalue,
                maximum_eigenvalue = shell_spectrum.maximum_eigenvalue,
                condition_number = shell_spectrum.condition_number,
            ),
        )
        accumulated = Matrix{Float64}(hcat(accumulated, shell_clean))
        push!(block_column_ranges, shell_range)
        push!(block_labels, Symbol("side$(side)_shell"))
    end

    contracted_weights = vec(transpose(parent_weights) * accumulated)
    overlap_spectrum = _experimental_high_order_metric_spectrum(accumulated, parent_overlap; tol = 1.0e-10)
    moment_risk = _experimental_high_order_stack_moment_risk(
        accumulated,
        parent_weights,
        parent_coordinates,
        block_labels,
        block_column_ranges,
    )
    diagnostics = (
        parent_mapping_family = mapping_family,
        parent_dimension = size(parent_overlap, 1),
        stack_dimension = size(accumulated, 2),
        parent_padding = length(basis) - maximum(side_values),
        shell_dimensions = Int[length(shell.shell_labels) for shell in shell_layers],
        shell_kind_counts = [shell.shell_kind_counts for shell in shell_layers],
        overlap_error = _experimental_high_order_overlap_error(accumulated, parent_overlap),
        overlap_spectrum = overlap_spectrum,
        shell_cleanup_spectra = shell_cleanup_spectra,
        moment_risk = moment_risk,
        contracted_weights_finite = all(isfinite, contracted_weights),
    )

    return ExperimentalHighOrderDosideStack3D(
        basis,
        backend,
        doside,
        side_values,
        length(basis),
        accumulated,
        block_column_ranges,
        block_labels,
        shell_layers,
        contracted_weights,
        diagnostics,
    )
end
